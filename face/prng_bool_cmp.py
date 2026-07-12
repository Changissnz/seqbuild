"""
the measures in this file are built from <desi.multi_metric>. They output boolean values in 
comparing two PRNGs for qualities of randomness. The two primary statistical tests <seqbuild> 
is equipped with, QUALTEST and CHAINTEST, are the input source for these boolean values. To 
repeat on the values of these two tests: 

- QUALTEST 
[0] coverage : 0 <= x <= 1
    unidirectional weighted point distance: 0 <= x <= 1
    categorical entropy: 0 <= x <= 1
    K-complexity of 1st order ternary difference 
    K-complexity of 2nd order ternary difference
    K-complexity of (sequence - most common subsequence) 
[1] factor modular k-complexity: 
        factor -> (1st order k-complexity, 2nd order k-complexity)  
[2] element -> frequency map 

- CHAINTEST: 4 x 2 matrix 

    min(mean)               max(mean) 
    var(min)                var(max) 
    mean(F_i - F_{i+1})     mean(F_{i+1} - F_i) 
    var(F_i - F_{i+1})      var(F_{i+1} - F_i) 

"""
from desi.multi_metric import * 
from collections import deque 
from intigers.mod_prng import prg__iterable 

SCORE_TERMS = {"cov","uwpd","catent","tkc1","tkc2","mcstkc","ftkc","efreq",\
    "vmean","vvar","fmean","fvar"} 

"""
container holding scores from a (QUALTEST,CHAINTEST) pair, conducted on an 
arbitrary number sequence. Used for easy retrieval of values based on string 
terms. 
"""
class QCInfo: 

    def __init__(self,qc): 
        assert len(qc) == 2 

        assert len(qc[0]) == 3 
        #assert is_vector(qc[0][0]) and len(qc[0][0]) == 6 
        assert type(qc[0][1]) in {dict,defaultdict}
        assert type(qc[0][2]) in {dict,defaultdict}

        assert is_2dmatrix(qc[1])
        assert qc[1].shape == (4,2) 

        self.q = qc[0]
        x = list(self.q[0][0]) + list(self.q[0][1]) + [self.q[0][2]]
        self.q = (x,self.q[1],self.q[2])

        self.c = qc[1] 

    def retrieve(self,term):
        assert term in SCORE_TERMS 

        if term == "cov": 
            return self.q[0][0] 

        if term == "uwpd":
            return self.q[0][1]             
            
        if term == "catent":
            return self.q[0][2] 

        if term == "tkc1":
            return self.q[0][3] 
            
        if term == "tkc2": 
            return self.q[0][4] 
            
        if term == "mcstkc":
            return self.q[0][5] 

        if term == "ftkc":
            X = np.array(list(self.q[1].values()))
            return zero_div(np.sum(X),len(X),float('inf')) 
            
        if term == "efreq": 
            return len([k for k,v in self.q[2].items() if v > 0])

        if term == "vmean":
            return (self.c[0,1] - self.c[0,0])

        if term == "vvar":
            return (self.c[1,1] - self.c[1,0]) / 2.0

        if term == "fmean":
            return (self.c[2,1] - self.c[2,0])

        return (self.c[3,1] - self.c[3,0]) / 2.0

def cmp_PRNG_test_scores(s0,s1): 

    if s0 > s1: 
        return 1 
    elif s1 > s0: 
        return -1 
    return 0 

def PRNG_process_metrics(qc_info1,qc_info2,variables): 
    qc1 = QCInfo(qc_info1) 
    qc2 = QCInfo(qc_info2)

    l = [] 
    for v in variables: 
        s0 = qc1.retrieve(v)
        s1 = qc2.retrieve(v) 
        ##print("V: {}\t{},{}".format(v,s0,s1))
        x = cmp_PRNG_test_scores(s0,s1)
        l.append(x) 

    return to_trinary_relation_v2(sum(l),None,True,False)

"""
all measures, from both QUALTEST and CHAINTEST 
"""
def PRNG_QC_metric_1(qc_info1,qc_info2): 
    relevant_vars = sorted(SCORE_TERMS) 
    return PRNG_process_metrics(qc_info1,qc_info2,relevant_vars)

"""
only measures from CHAINTEST 
"""
def PRNG_QC_metric_2(qc_info1,qc_info2): 
    relevant_vars = ["vmean","vvar","fmean","fvar"] 
    return PRNG_process_metrics(qc_info1,qc_info2,relevant_vars)

"""
1st,2nd and MCS K-complexity measures from QUALTEST, all measures from 
CHAINTEST 
"""
def PRNG_QC_metric_3(qc_info1,qc_info2): 
    relevant_vars = ["tkc1","tkc2","mcstkc","vmean","vvar","fmean","fvar"]
    return PRNG_process_metrics(qc_info1,qc_info2,relevant_vars)

"""
all K-complexities and the number of unique elements from QUALTEST, all measures 
from CHAINTEST 
"""
def PRNG_QC_metric_4(qc_info1,qc_info2): 
    relevant_vars = ["tkc1","tkc2","mcstkc","ftkc","efreq",\
        "vmean","vvar","fmean","fvar"]
    return PRNG_process_metrics(qc_info1,qc_info2,relevant_vars)

"""
all measures except F-measures (modular factor complexities). 
"""
def PRNG_QC_metric_5(qc_info1,qc_info2):
    relevant_vars = ["cov","uwpd","catent","tkc1","tkc2","mcstkc","efreq",\
        "vmean","vvar"]  
    return PRNG_process_metrics(qc_info1,qc_info2,relevant_vars)

"""
[0] measures from QUALTEST
"""
def PRNG_QC_metric_6(qc_info1,qc_info2): 

    relevant_vars = ["cov","uwpd","catent","tkc1","tkc2","mcstkc"] 
    return PRNG_process_metrics(qc_info1,qc_info2,relevant_vars)

"""
all measures from QUALTEST
"""
def PRNG_QC_metric_7(qc_info1,qc_info2): 
    relevant_vars = ["cov","uwpd","catent","tkc1","tkc2","mcstkc","ftkc","efreq"] 
    return PRNG_process_metrics(qc_info1,qc_info2,relevant_vars)

PRNG_RANDOMNESS_METRIC_MAP = {1:PRNG_QC_metric_1,\
    2:PRNG_QC_metric_2,3:PRNG_QC_metric_3,\
    4:PRNG_QC_metric_4,5:PRNG_QC_metric_5,\
    6:PRNG_QC_metric_6,7:PRNG_QC_metric_7} 

"""
primary method in this file. 
"""
def is_PRNG_more_QC_random(qc_info1,qc_info2,metric_numbers:set): 
    l = []
    for m in metric_numbers: 
        assert m in PRNG_RANDOMNESS_METRIC_MAP

        F = PRNG_RANDOMNESS_METRIC_MAP[m] 
        x = F(qc_info1,qc_info2)
        l.append(x) 

    return to_trinary_relation_v2(sum(l),None,True,False)