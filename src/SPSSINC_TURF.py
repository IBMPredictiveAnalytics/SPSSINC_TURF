from __future__ import with_statement
#/***********************************************************************
# * Licensed Materials - Property of IBM 
# *
# * IBM SPSS Products: Statistics Common
# *
# * (C) Copyright IBM Corp. 1989, 2012
# *
# * US Government Users Restricted Rights - Use, duplication or disclosure
# * restricted by GSA ADP Schedule Contract with IBM Corp. 
# ************************************************************************/

# TURF Analysis

__author__ = "SPSS, JKP"
__version__ = "2.5.0"

# history
# 05-Feb-2009  original version
# 07-Feb-2009  implement optional remove cycles
# 11-Feb-2009  add case weight support
# 23-Feb-2009  do percentage calculations in floating point for extra decimals
# 09-Mar-2009  evaluate candidates as accumulated to reduce memory and post-processing time
# 11-Apr-2009   don't use named tuples as they don't work with extended characters
# 15-Dec-2009  enable localization
# 05-Apr-2010  Add minimum parameter for maximum group size sequence
# 16-Apr-2010  Add logging capability
# 23-jun-2010  translatable procedure name
# 04-jun-2011  Add summary table and chart
# 04-feb-2012 Add importance weights

import spss, spssaux, spssdata
from extension import Template, Syntax, processcmd
from itertools import imap
from operator import mul

import copy, heapq, time, os, random

def turf(variables, bestn, number, threshold=0, criteria=None, removalcycles = 0, 
         maxsets=None, mintable=1, logfile=None, logfreq=500000, iweights=None, strict=True):
    """Calculate best reaching combinations eliminating duplicates
    
    variables is a sequence of variable names.  Each variable must be coded so that
    value=1 is the response to tabulate
    bestn is the depth, i.e., number of questions that can be combined.  bestn = 1 gives
    only single frequencies, bestn = 2 gives best sets of one or two questions.
    number is the maximum nfumber of groups to display
    threshold is minimum percentage required to consider a variable.  If the percentage
    of positives is lower, the variable is discarded.  Variables with 0 positive responses
    are always discarded.
    criteria is a sequence of values that define positive responses.
    number is the number of items to display
    removalcycless is the number of cycles to do after initial, removing best variables one at a time
    If maxsets is not None, the analysis is not run if the number of sets would exceed it
    iweights is a list of importance weights
    strict indicates the rules for validating iweights

    
    caseweights, if any, are honored.
"""
    
    info = NonProcPivotTable("INFORMATION", tabletitle=_("Information"),
        columnlabels=["Count"])
    if not (0. <= threshold <= 100.):
        raise ValueError(_("Threshold value must be between 0 and 100"))
    if mintable >bestn:
        raise ValueError(_("Invalid value for the minimum for the Maximum Group Size tables"))
    mintable -= 1
    
    if criteria is None:
        criteria = [1.]
    removelistnames = []  # accumulates variables removed on remove cycles
    variables = spssaux._buildvarlist(variables)
    iwt = Iweights(variables, iweights, strict)
    cycle = 0
    log = Writelog(logfile)
    log.write("Set operation count frequency: %d" % logfreq)
    try:
        StartProcedure(_("TURF analysis"), "SPSSINC TURF")
        # when doing remove cycles, there are two possible approaches.  One is to do a data pass
        # for each cycle and filter the excluded cases as sets are recreated.  The other is to subtract
        # the cases in the best set variable from each other question set without extra data passes.
        # Either could be faster under certain circumstances, but the set approach is cleaner and
        # probably faster in most typical cases, so it is used here.
        

        iwt.displayWeights()
        allsummaryresults = list()   # saving for charting
        while removalcycles >= 0:
            log.write("Beginning cycle.  removalcycles=%d" % removalcycles)
            removalcycles -= 1
            if removelistnames:  # in a remove cycle
                cycle += 1
                removeinfo = "\nSuperior Variables Removed: %s" % ", ".join(removelistnames)
                for qn in range(nquest):
                    if qn != removeindex:
                        # remove cases where current best var has positive response
                        qq[qn] = qq[qn] - qq[removeindex]
                totalcases -= len(qq[removeindex])
                qq.pop(removeindex)
                variables.pop(removeindex)
                iwt.pop(removeindex)
                nquest -= 1
            else:
                removeinfo = ""
                nquest = len(variables)
                qq = []
                for i in range(nquest):
                    qq.append(set())  # one set for each question.  Sets will hold case numbers of positive responses
                totalcases = float(0)
                
                caseweights = {}   # weighted counts are used for length if data are weighted
                caseweightvar = spss.GetWeightVar()
                if caseweightvar:
                    vv = variables + [caseweightvar]
                else:
                    vv = variables
                log.write("Beginning data pass")
                curs = spssdata.Spssdata(vv, names=False)  # to avoid problems with extended characters
                for casenum, case in enumerate(curs):
                    # For weighted data, cases with system-missing or nonpositive weights are automatically skipped.
                    if caseweightvar:
                        caseweights[casenum] = case[-1]
                        totalcases += case[-1]
                    else:
                        totalcases += 1
                    for i in range(nquest):
                        if case[i] in criteria:
                            qq[i].add(casenum)   # add casenum to set if value is 1
                curs.CClose()                
            
            wlen = Wlen(caseweights)
            # count total responses by question and prune questions below threshold or all zero
            log.write("Pruning questions")
            totalthresh = totalcases * threshold/100.
            totalresponses = float(0)
            qcounts = []
            for qn in range(nquest-1, -1, -1):
                numresp = wlen.wtlen(qq[qn])
                if numresp >= totalthresh and numresp > 0:
                    qcounts.insert(0, numresp)
                    totalresponses += numresp
                else:
                    qq.pop(qn)   # discard the question
                    variables.pop(qn)
                    iwt.pop(qn)
            nquest = len(qq)
            if nquest == 0:
                raise ValueError(_("Analysis stopped.  All variables have a positive count below the threshold of %s positive responses.") % totalthresh)
            bestn = min(bestn, nquest)
            mintable = min(mintable, bestn)   #always do at least one round
            
            unionsrequired = calcsets(nquest, bestn)
            info.addrow(rowlabel=_("Set unions to be calculated"),
                        cvalues=[unionsrequired])
            log.write("Unions required: %d" % unionsrequired)
            if maxsets is not None and unionsrequired > maxsets:
                raise ValueError(_("Analysis stopped because the required set union operations would exceed the specified limit"))
            
            visitor = Counter(qq, qcounts, bestn, number, wlen, log, logfreq, iwt, caseweights, casenum)
            for i in range(nquest):
                log.write("Visiting variable %d" % i)
                visitor.visit([], i)
            
            # Display results in pivot table, one for each number of items up to size bestn
    
            summarycats = list()   # for accumulating best results across group sizes
            summarystats = list()
            for b in range(bestn):
                log.write("Calculating best %i" % b)
                reach, freq, wtmean, wttotal = visitor.best(b+1, number)
                if b ==0:
                    removeindex = reach[0][1][0]
                    removelistnames.append(variables[removeindex])  # variable name of current best variable
                if b < mintable:
                    continue
    
                pt = spss.BasePivotTable(_("Maximum Group Size: %s.  Reach and Frequency.") % \
                    (b+1),
                    "REACHFREQ", _("Cycle: %s") % cycle, False, 
                    caption = _("Variables: ") + ", ".join(variables) + removeinfo)
                rowdim = pt.Append(spss.Dimension.Place.row, _("Variables"))
                coldim = pt.Append(spss.Dimension.Place.column, _("Statistics"))
                categories = [spss.CellText.String(_("Reach")), 
                    spss.CellText.String(_("Pct of Cases")),
                    spss.CellText.String(_("Frequency")),
                    spss.CellText.String(_("Pct of Responses"))]
                if iwt.iweights is not None:
                    categories.extend([spss.CellText.String(_("""Mean Importance of Positive Responses""")),
                        spss.CellText.String(_("""Total Importance of Positive Responses"""))])
                pt.SetCategories(coldim, categories)
                rowcats = []
                for ii, item in enumerate(reach):
                    names = ", ".join(sorted([variables[i] for i in item[1]]))
                    rowcats.append(spss.CellText.String(names))
                    if ii == 0:
                        summarycats.append(sorted([variables[i] for i in item[1]]))
                pt.SetCategories(rowdim, rowcats)
                ###summarycats.append(rowcats[0])
                
                for i, row in enumerate(rowcats):
                    values = [spss.CellText.Number(reach[i][0], spss.FormatSpec.Count), 
                        spss.CellText.Number(100 * float(reach[i][0]) / totalcases, spss.FormatSpec.PercentNoSign),
                        spss.CellText.Number(freq[i], spss.FormatSpec.Count),
                        spss.CellText.Number(100 * float(freq[i])/totalresponses, spss.FormatSpec.PercentNoSign)]
                    if iwt.iweights is not None:
                        values.extend([spss.CellText.Number(wtmean[i], spss.FormatSpec.GeneralStat),
                            spss.CellText.Number(wttotal[i], spss.FormatSpec.GeneralStat)])
                    pt.SetCellsByRow(row, values)
                summarystats.append([b + 1, reach[0][0], 100 * float(reach[0][0] / totalcases), 
                    freq[0], 100 * float(freq[0] / totalresponses)])
                
            info.addrow(rowlabel=_("Sets of variables analyzed (unions plus one-variable sets)"),
                        cvalues=[visitor.setops])
            
            # produce summary table for this cycle
            if bestn > 1:
                allsummaryresults.append(summarystats)   # collect for charting
                spt = spss.BasePivotTable(_("Best Reach and Frequency by Group Size"),
                    "SUMMARYREACHFREQ", _("Cycle: %s") % cycle, False)
                rowdim = spt.Append(spss.Dimension.Place.row, _("Variables"))
                coldim = spt.Append(spss.Dimension.Place.column, _("Statistics"))
                spt.SetCategories(coldim, [spss.CellText.String(_("Group Size")), spss.CellText.String(_("Reach")), 
                    spss.CellText.String(_("Pct of Cases")),
                    spss.CellText.String(_("Frequency")),
                    spss.CellText.String(_("Pct of Responses"))])
                summarycats = flagnew(summarycats)
                spt.SetCategories(rowdim, summarycats)

                for i, row in enumerate(summarycats):
                    spt.SetCellsByRow(row, [spss.CellText.Number(summarystats[i][0], spss.FormatSpec.Count),
                        spss.CellText.Number(summarystats[i][1], spss.FormatSpec.Count), 
                        spss.CellText.Number(summarystats[i][2], spss.FormatSpec.PercentNoSign),
                        spss.CellText.Number(summarystats[i][3], spss.FormatSpec.Count),
                        spss.CellText.Number(summarystats[i][4], spss.FormatSpec.PercentNoSign)])
    finally:
        info.generate()
        spss.EndProcedure()
        if bestn > 1:
            summaryCharts(allsummaryresults)
        log.close()
        
class Iweights(object):
    """Manage importance weights"""
    
    def __init__(self, variables, iweights, strict):
        """variables is the list of variables being analyzed
        iweights is the list of importance weights as pairs of name, value
        strict is True if all variables must appear in iweights"""
        
        if iweights is None:
            self.iweights = None
            return
        
        self.variables = variables
        # iweights is a dictionary keyed by cased variable name containing the importance weight
        if not strict:
            self.iweights = dict([(v, 1.) for v in variables]) # default weight is 1.
        else:
            self.iweights = {}
        lvariables = dict([(item.lower(), item) for item in variables])

        # validation
        for index in range(0, len(iweights), 2):
            name= iweights[index]
            try:
                wt = float(iweights[index+1])
            except:
                raise ValueError(_("""An importance weight is missing or not numeric.  Variable: %s""") % name)                
            if not name.lower() in lvariables:
                raise ValueError(_("""An importance weight was specified for a variable not in the analysis: %s""") % name)
            self.iweights[lvariables[name.lower()]] = wt
        if strict and len(self.iweights) != len(lvariables):
            raise ValueError(_("""Some variables were not assigned an importance weight: %s""")\
            % ", ".join(set(lvariables.values()) - set(self.iweights)))
        
        self.index = dict([(i, name) for i, name in enumerate(variables)])
        # eliminate indirect addressing
        for k, v in self.index.items():
            self.iweights[k] = self.iweights[v]
        for k in variables:
            del(self.iweights[k])
        
    def getwt(self, index):
        "Return weight for question"
        
        if self.iweights is None:
            return 1
        else:
            #return self.iweights[self.index[index]]  # TODO: direct index
            return self.iweights[index]  # TODO: direct index
        
    def pop(self, index):
        """Remove variable and adjust indexes"""
        
        # Not efficient, but these dicts are generally small and
        # removals are infrequent
        
        if self.iweights is None:
            return
        #vn = self.index[index]
        #del(self.iweights[vn])
        #del(self.index[index])
        #for key in sorted(self.index.keys()):
            #if key > index:
                #self.index[key-1] = self.index[key]
        del self.iweights[index]
        for key in sorted(self.iweights.keys()):
            if key > index:
                self.iweights[key-1] = self.iweights[key]

    def displayWeights(self):
        """Display the importance weights, if any
        
        This function is only valid if variable list has not be altered since the weights were created"""
        
        if self.iweights is None:  # no weights
            return
        wtbl = spss.BasePivotTable(_("""Importance Weights"""), "TURFIMPORTANCEWTS")
        cells = []
        rowlabels = []
        #for item in sorted(self.iweights):
        for i, item in enumerate(self.variables):
            rowlabels.append(item)
            #cells.append(self.iweights[item])
            cells.append(self.iweights[i])
        wtbl.SimplePivotTable(rowdim=_("""Variable"""), rowlabels = rowlabels,
            collabels=[_("""Importance Weights""")],
            cells = cells)

def flagnew(summarycats):
    """Categorize entries in list of names and return as a string
    
    summarycats is a list.  Each element in the list is the list of variable names for
    that group size."""
    
    # for each item in summarycats, mark the entries that were not in the previous entry

    previousset = set()
    for i, item in enumerate(summarycats):  # item is a list variable names
        currentset = set(item)
        added = currentset - previousset
        removed = previousset - currentset
        retained = currentset.intersection(previousset)
        item = "ADDED: " + ", ".join(sorted(added))
        if retained:
            item = item + "\nKEPT: " + ", ".join(sorted(retained))
        if removed:
            item = item + "\nDROPPED: " + ", ".join(sorted(removed))
        summarycats[i] = spss.CellText.String(item)
        previousset = currentset
    return summarycats
   
def summaryCharts(results):
    """Produce a summary line chart for each removal cycle
    
    results is a list of best statistics for each group size
    each item contains group size, reach, reach %, freq, freq %"""

    GroupSize = _("Group Size")
    Percentage = _("Percentage")
    title = _("Reach and Frequency by Group Size")
    PctofCases = _("Pct of Cases")
    PctofResponses = _("Pct of Responses")
    
    chartsyntax = r"""GGRAPH
  /GRAPHDATASET NAME="graphdataset" 
  VARIABLES=GroupSize MEAN(Reach) 
  MEAN(Frequency)
    TRANSFORM=VARSTOCASES(SUMMARY="#SUMMARY" INDEX="#INDEX")
  /GRAPHSPEC SOURCE=INLINE.
BEGIN GPL
  SOURCE: s=userSource(id("graphdataset"))
  DATA: GroupSize=col(source(s), name("GroupSize"), unit.category())
  DATA: SUMMARY=col(source(s), name("#SUMMARY"))
  DATA: INDEX=col(source(s), name("#INDEX"), unit.category())
  GUIDE: axis(dim(1), label("%(GroupSize)s"))
  GUIDE: axis(dim(2), label("%(Percentage)s"))
  GUIDE: legend(aesthetic(aesthetic.color.interior), label(""))
  GUIDE: text.title(label("%(title)s"))
  %(subtitle)s
  SCALE: linear(dim(2), include(0, 100))
  SCALE: cat(aesthetic(aesthetic.color.interior), include("0", "1"))
  ELEMENT: line(position(GroupSize*SUMMARY), color.interior(INDEX))
  ELEMENT: point(position(GroupSize*SUMMARY), color.interior(INDEX))
END GPL."""

    # get active dataset name - assign one if unnamed
    
    main = spss.ActiveDataset()
    if main == "*":
        main = "T" + str(random.uniform(0,1))
        spss.Submit("DATASET NAME %s" % main)
    
    # create a reach and frequency chart for each cycle    
    for cycle, result in enumerate(results):
        spss.StartDataStep()
        ds = spss.Dataset(name=None, hidden=False)
        dsname = ds.name
        ds.varlist.append("GroupSize")
        ds.varlist["GroupSize"].format = (5,5,0)
        ds.varlist.append("Reach")
        ds.varlist.append("Frequency")
        for case in result:
            ds.cases.append([case[0], case[2], case[4]])
        spss.EndDataStep()
        spss.Submit("DATASET ACTIVATE %s" % dsname)

        if len(results) > 1:
            subtitle = """GUIDE: text.subtitle(label("%s"))""" % (_("Cycle: %s") % cycle)
        else:
            subtitle = ""
        spss.Submit(chartsyntax % locals())
        spss.Submit("DATASET CLOSE %s" % dsname)
        
    spss.Submit("DATASET ACTIVATE %s" % main)
    
    
class Wlen(object):
    """Calculate weighted sum of set"""
    
    def __init__(self, caseweights):
        self.caseweights = caseweights
        self.weightsize = len(caseweights)
        
    def wtlen(self, qset):
        """Return weighted length of set
        
        qset is a set of integer case numbers.
        
        If data are unweighted, the case count (len) is returned.
        Otherwise, the sum of the weights of those cases is returned"""
        
        if self.weightsize == 0:
            return len(qset)
        else:
            return sum((self.caseweights[w] for w in qset))
        
        
class Counter(object):
    """Calculate reach counts for questions"""
    
    def __init__(self, qs, qcounts, bestn, number, wlen,log, logfreq, iwt, caseweights, ncases):
        """qs is a list of question sets
        qcounts is a list of the total positive response for each question
        bestn is the depth (number of questions) to calculate for
        wlen is the length-counting object
        iwt is the importance weight object
        caseweights is a dictionary of case weights - empty if unweighted
        ncases is the unweighted case count
        """
        
        # counts is a list with each element a duple.  First duple element is the count.
        # Second is a tuple of question numbers for that count
        # The question-number tuple will have up to bestn integers
        
        self.qs = qs
        self.nquest = len(qs)
        self.qcounts = qcounts   # response count for questions overall
        self.bestn = bestn
        self.counts = []
        self.hh = Countheaps(bestn, number)
        self.qsets = []
        self.wlen = wlen
        self.setops = 0
        self.log = log
        self.logfreq = logfreq
        self.iwt = iwt
        self.memomean = {}  # remembering importance-weighted means
        self.caseweights = caseweights
        self.casecount = ncases + 1
        
    def visit(self, nodelist, qnumber):
        """tally the current node and visit right children.
        
        nodelist is the path to the node
        qnumber is the question to add to the base"""
        
        # Each self.qs item is a list of positive cases for a question
        # The last set in self.qsets is the base set unless it is empty
        try:
            s = self.qsets[-1]
        except:
            s = set()
        #slen = self.wlen.wtlen(s)  # len before union
        self.qsets.append(s.union(self.qs[qnumber]))
        #slenafter = self.wlen.wtlen(self.qsets[-1])
        #delta = (slenafter - slen) * self.iwt.getwt(qnumber)
        #curlen = slen + delta
        newnodelist = nodelist + [qnumber]
        # len of set counts number of cases responding to any of the q's (reach)
        # sum of lengths of included questions will be the frequency

        #curlen = self.wlen.wtlen(self.qsets[-1])
        #curlen = curlen * self.iwt.getwt(qnumber)
        self.hh.add(self.wlen.wtlen(self.qsets[-1]), newnodelist)
        ###self.hh.add(curlen, newnodelist)
        self.setops += 1

        #if self.setops % 1000 == 1:
        #    self.debugf.write("setops: %d, qsets: %d\n" % (self.setops, len(self.qsets)))
        #    self.debugf.flush()
        if self.setops % self.logfreq == 1:
            self.log.write("Set union operations completed: %d" % self.setops)
        if len(newnodelist) < self.bestn:
            for i in range(qnumber+1, self.nquest):
                self.visit(newnodelist, i)
        self.qsets.pop()
        
    def best(self, depth, number):
        """Return best reaches and frequencies of specified depth
        Results are first sorted by reach and then by frequency
        
        depth is the maximum length of the counts vector
        number is the number of question combinations to return"""
        
        # Typically, number is << length of eligible items
        

        # combine all the count heaps where the number of variables is within the depth
        # limit screening out any zero counts that remain
        
        allatdepth = []        
        for i in range(depth):
            allatdepth.extend(filter(lambda x: x.count > 0, self.hh.heaplist[i]))
        
        
        allatdepth.sort(key=lambda x: x.count, reverse=True)
        # find all that could be eligible based on number and trim list to that number
        if number < len(allatdepth):
            lastelt = allatdepth[number].count
            for i in range(number, len(allatdepth)):
                if allatdepth[i].count < lastelt:    # can't be made eligible by a high freq
                    break
            allatdepth = allatdepth[:i]
        
        # construct a list of the total number of responses to the questions in each item
        # does not currently consider weights
        freqs = []
        for item in allatdepth:
            freqs.append(sum([self.qcounts[i] for i in item.varlist]))
        
        combined = [(item1, item2) for (item1, item2) in zip(allatdepth, freqs)]
        
        # sort by reach and then by freq
        combined.sort(key=sortkey, reverse=True)
        
        allatdepth = [(item[0].count, item[0].varlist) for item in combined]
        freqs = [item[-1] for item in combined]
        if self.iwt.iweights is not None:
            atmost = min(number, len(allatdepth))
            wtstats = [self.wtMean(allatdepth[i][-1]) for i in range(atmost)]
            return (allatdepth[:number], freqs[:number], [item[0] for item in wtstats], [item[1] for item in wtstats])
        else:
            return (allatdepth[:number], freqs[:number], None, None)

    def wtMean(self, qlist):
        """Return mean weighted reach and total for list of questions in qlist
        
        qlist is a list of question numbers
        mean is casewise mean for positive responses weighted by importance
        weights"""
        
        # memoize results, since a combination can recur
        # each q is represented by a list of cases with a positive response
        qlist = tuple(qlist)
        if qlist in self.memomean:
            return self.memomean[qlist]
        ssum = 0.
        scount = 0
        wts = [self.iwt.getwt(v) for v in qlist]
        for case in range(self.casecount):
            isin = [case in self.qs[q] for q in qlist]
            qkt = sum(isin)
            if qkt > 0:
                casewt = self.caseweights.get(case, 1.)
                wtd = sum(imap(mul, isin, wts))
                casemean = wtd / qkt
                ssum += casemean * casewt
                scount += casewt
        result = (ssum/scount, ssum)
        self.memomean[qlist] = result
        return result

def sortkey(x):
    return (x[0].count, x[1])


fac = lambda n: n and fac(n-1)*long(n) or 1

def binom(n, k):
    """binomial coefficient
    n and k should be longs"""
    
    # not worrying about overflow
    
    if k < 0 or k > n:
        return 0
    return fac(n)/(fac(k) * fac(n-k))
    

def calcsets(nvars, depth):
    """Return number of sets required for computation
    
    nvars is the number of variables to analyze
    depth is the maximum number of variables to be combined"""

    nvars = long(nvars)
    nfact = fac(nvars)
    count = 0L
    
    for i in range(1, depth):
        count += binom(nvars, i+1)
    return count
        


class C(object):
    """Class for making count objects work with heaps"""
    __slots__ = ["count", "varlist"]        # to reduce overhead as there will be many of these
    def __init__(self, count, varlist):
        self.count = count
        self.varlist = copy.copy(varlist)
    def __eq__(self, other):
        return self.count == other.count
    def __ne__(self, other):
        return self.count != other.count
    def __lt__(self, other):
        return self.count < other.count
    def __le__(self, other):
        return self.count <= other.count
    def __gt__(self, other):
        return self.count > other.count
    def __ge__(self, other):
        return self.count >= other.count
    
class Countheaps(object):
    """maintain a list of heaps for counts"""
    def __init__(self, depth, number):
        """depth is the number of heaps to be maintained: one for each number of variables
        number is the minimum number of items to be kept.
        
        Each heap object has a count and varlist."""
        
        self.heaplist = []
        self.depth = depth
        self.number = number
        for i in range(depth):
            self.heaplist.append([C(0,[])])
            heapq.heapify(self.heaplist[-1])
            
    def add(self, count, varlist):
        """Add an item to the appropriate heap if it qualifies
        count is the count data
        varlist is the list of variables it is based on"""
        
        # Selected heap is determined by the size of varlist - 1 is the first heap
        # If the heap size is less than self.number, always add the item
        # If the size is >= number and count is > min element, remove inferior elements
        # subject to min heap size, and add the new element

        h = self.heaplist[len(varlist)-1]  # pick the right heap

        if len(h) < self.number:
            heapq.heappush(h, C(count, varlist))
            return
        if count < h[0].count:
            return
        while count > h[0].count and len(h) >= self.number:
            heapq.heappop(h)
        heapq.heappush(h, C(count, varlist))

helptext="""
Calculate TURF values for specified variables.
TURF is "total unduplicated reach and frequency".

SPSSINC TURF VARIABLES=variable-list
/OPTIONS BESTN=n [MINTABLE=m] NUMBERTODISPLAY=m
[THRESHOLD = value]
[CRITERIA = {1*|list-of-values}]
[MAXSETS=number]
[REMOVALCYCLES={0*|n}] [LOGFILE=filespec] [LOGFREQ=value]
[/IMPORTANCE IWEIGHTS=varname wt varname wt ... [STRICT={YES*|NO}]
[/HELP]

Example:
SPSSINC TURF VARIABLES=q1 q2 q3 q4 q5
/OPTIONS BESTN=3 NUMBERTODISPLAY=10.

VARIABLES lists the variables to analyze.  All the variables must be coded
the same way, at least for the definition of positive response.

BESTN specifies the maximum group size, i.e., number of variables, that
will be analyzed.

MINTABLE can specify the minimum group size table to display.  Tables
will be displayed from MINTABLE to BESTN.  By default, tables start with 1.

NUMBERTODISPLAY indicates the maximum number of combinations to report
in each table.

THRESHOLD can optionally specify a minimum percentage of positive responses
required in order to analyze a variable.  If a variable has fewer than that, 
it is discarded from the analysis.  
Variables with no positive responses are always discarded.

CRITERIA optionally specifies one or more variable values that are considered a
positive response.  By default, values of 1 are considered positive.

Use REMOVALCYCLES to iterate the entire analysis by successively removing each
best variable.  The top of the group size 2 list identifies which variable combines best 
with the top of the size 1 list, but iterating the analysis with removal cycles would show,
for example, a table with the best single variables for cases that did not have a positive 
response to the top of the size 1 list and so on.  The precentages would then be 
calculated on the reduced samples.  By default, no removal cycles are carried out.

If MAXSETS is specified, the analysis is not run if the number of set union 
operations required exceeds it.

if a file specification is given in LOGFILE, an operations log will be written to that file.
The contents, which are timestamped in seconds, report various stages in the calculation.
For large problems, the time is dominated by the number of set operations required.
A log entry is made every LOGFREQ set union operations, so it should be possible to 
extrapolate to the approximate completion time.  Entries are flushed to disk i
mmediately, so you can watch the log as the job progresses.

LOGFREQ is the frequency of logging completed set operations.  The default is 500000.

The optional IMPORTANCE subcommand provides importance weights for variables.
By default, all variables have an importance weight of 1.
The IWEIGHTS keyword is a list of variable names and weights.
If STRICT is YES, every variable used must be listed with a weight.
If STRICT is NO, a weight of 1 is assumed for any variable not listed.

Computational notes:
The computational burden of this calculation increases very rapidly with the
number of variables and the value of BESTN.

All data on positive responses for these calculations must be held in memory.

If the dataset is weighted, only cases with positive weights are counted, and the 
counts are sums of weights with no rounding or truncation.

/HELP displays this help and does nothing else.
"""

def Run(args):
    """Execute the SPSSINC TURF command"""

    args = args[args.keys()[0]]
    ###print args   #debug
    
    ##debugging
    try:
        import wingdbstub
        if wingdbstub.debugger != None:
            import time
            wingdbstub.debugger.StopDebug()
            time.sleep(2)
            wingdbstub.debugger.StartDebug()
    except:
        pass

    oobj = Syntax([
        Template("VARIABLES", subc="",  ktype="existingvarlist", var="variables", islist=True),
        Template("BESTN", subc="OPTIONS",  ktype="int", var="bestn", vallist=[1]),
        Template("MINTABLE", subc="OPTIONS",  ktype="int", var="mintable", vallist=[1]),
        Template("NUMBERTODISPLAY", subc="OPTIONS", ktype="int", var="number", vallist=[1]),
        Template("THRESHOLD", subc="OPTIONS", ktype="float", var="threshold"),
        Template("CRITERIA", subc="OPTIONS", ktype="int", var="criteria", islist=True),
        Template("REMOVALCYCLES", subc="OPTIONS", ktype="int", var="removalcycles", vallist=[0]),
        Template("MAXSETS", subc="OPTIONS", ktype="int", var="maxsets"),
        Template("LOGFILE", subc="OPTIONS", ktype="literal", var="logfile"),
        Template("LOGFREQ", subc="OPTIONS", ktype="int", var="logfreq",vallist=[1]),
        Template("IWEIGHTS", subc="IMPORTANCE", ktype="str", var="iweights", islist=True),
        Template("STRICT", subc="IMPORTANCE", ktype="bool", var="strict"),
        Template("HELP", subc="", ktype="bool")])
    
        # ensure localization function is defined
    global _
    try:
        _("---")
    except:
        def _(msg):
            return msg

        # A HELP subcommand overrides all else
    if args.has_key("HELP"):
            print helptext
    else:
            processcmd(oobj, args, turf, vardict=spssaux.VariableDict())
            
class NonProcPivotTable(object):
    """Accumulate an object that can be turned into a basic pivot table once a procedure state can be established"""
    
    def __init__(self, omssubtype, outlinetitle="", tabletitle="", caption="", rowdim="", coldim="", columnlabels=[],
                 procname="Messages"):
        """omssubtype is the OMS table subtype.
        caption is the table caption.
        tabletitle is the table title.
        columnlabels is a sequence of column labels.
        If columnlabels is empty, this is treated as a one-column table, and the rowlabels are used as the values with
        the label column hidden
        
        procname is the procedure name.  It must not be translated."""
        
        attributesFromDict(locals())
        self.rowlabels = []
        self.columnvalues = []
        self.rowcount = 0

    def addrow(self, rowlabel=None, cvalues=None):
        """Append a row labelled rowlabel to the table and set value(s) from cvalues.
        
        rowlabel is a label for the stub.
        cvalues is a sequence of values with the same number of values are there are columns in the table."""
        
        if cvalues is None:
            cvalues = []
        self.rowcount += 1
        if rowlabel is None:
            self.rowlabels.append(str(self.rowcount))
        else:
            self.rowlabels.append(rowlabel)
        if not spssaux._isseq(cvalues):
            cvalues = [cvalues]
        self.columnvalues.extend(cvalues)
        
    def generate(self):
        """Produce the table assuming that a procedure state is now in effect if it has any rows."""
        
        privateproc = False
        if self.rowcount > 0:
            try:
                table = spss.BasePivotTable(self.tabletitle, self.omssubtype)
            except:
                StartProcedure(_("Messages"), self.procname)
                privateproc = True
                table = spss.BasePivotTable(self.tabletitle, self.omssubtype)
            if self.caption:
                table.Caption(self.caption)
            # Note: Unicode strings do not work as cell values in 18.0.1 and probably back to 16
            if self.columnlabels != []:
                table.SimplePivotTable(self.rowdim, self.rowlabels, self.coldim, self.columnlabels, self.columnvalues)
            else:
                table.Append(spss.Dimension.Place.row,"rowdim",hideName=True,hideLabels=True)
                table.Append(spss.Dimension.Place.column,"coldim",hideName=True,hideLabels=True)
                colcat = spss.CellText.String("Message")
                for r in self.rowlabels:
                    cellr = spss.CellText.String(r)
                    table[(cellr, colcat)] = cellr
            if privateproc:
                spss.EndProcedure()
                
def attributesFromDict(d):
    """build self attributes from a dictionary d."""
    self = d.pop('self')
    for name, value in d.iteritems():
        setattr(self, name, value)
        
class Writelog(object):
    """Manage a log file"""
    
    def __init__(self, logfile):
        """logfile is the file name or None"""

        self.logfile = logfile
        if self. logfile:
            self.file = open(logfile, "w")
            self.starttime = time.time()
            self.file.write("%.2f %s Starting log\n" % (time.time() - self.starttime, time.asctime()))
            
    def __enter__(self):
        return self
    
    def write(self, text):
        if self.logfile:
            self.file.write("%.2f: %s\n" % (time.time() - self.starttime,  text))
            self.file.flush()
            
    def close(self):
        if self.logfile:
            self.write("Closing log")
            self.file.close()

def StartProcedure(procname, omsid):
    """Start a procedure
    
    procname is the name that will appear in the Viewer outline.  It may be translated
    omsid is the OMS procedure identifier and should not be translated.
    
    Statistics versions prior to 19 support only a single term used for both purposes.
    For those versions, the omsid will be use for the procedure name.
    
    While the spss.StartProcedure function accepts the one argument, this function
    requires both."""
    
    try:
        spss.StartProcedure(procname, omsid)
    except TypeError:  #older version
        spss.StartProcedure(omsid)