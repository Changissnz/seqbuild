from morebs2.numerical_extras import * 
from morebs2.graph_basics import flatten_setseq

from .seq_struct import * 

DEFAULT_INT_MAX_THRESHOLD = 10 ** 6 

def factors_of_seq(seq): 
    assert type(seq) == IntSeq
    S = [] 
    for n in seq.l:      
        S.append(all_multiples(n)) 
    return S 

def intersection_disjunction__seq_factors(ms): 
    """
    :return: set of multiples shared by all integers, set of multiples not 
        shared by any integers
    :rtype: (set,set) 
    """

    if len(ms) == 0: return None 

    all_mult = flatten_setseq(ms) 
    q = ms[0] 
    for i in range(1,len(ms)): 
        q = q.intersection(ms[i]) 

    for i in range(0,len(ms) - 1): 
        for j in range(i+1,len(ms)): 
            x = ms[i].intersection(ms[j]) 
            all_mult -= x 
    return q,all_mult

"""
set operations of multiples for an integer sequence
"""
class ISFactorSetOps: 

    def __init__(self,l,int_limit=DEFAULT_INT_MAX_THRESHOLD):  
        assert len(l) >= 2
        assert abs(min(l)) <= DEFAULT_INT_MAX_THRESHOLD and \
            abs(max(l)) <= DEFAULT_INT_MAX_THRESHOLD
        self.iseq = IntSeq(l) 
        self.factors = factors_of_seq(self.iseq) 
        self.factor_count = None 
        self.max_cofactor_degree = None 

    def factor_count_(self): 
        self.factor_count = defaultdict(int) 
        self.max_cofactor_degree = 0 
        for mx in self.factors: 
            for q in mx: 
                self.factor_count[q] += 1 
                self.max_cofactor_degree = self.max_cofactor_degree if \
                    self.factor_count[q] < self.max_cofactor_degree else \
                    self.factor_count[q] 

    def cofactor_degree_set(self,d):
        assert d >= 1
        if d > self.max_cofactor_degree: return None 

        s = set() 
        for k,v in self.factor_count.items():
            if v == d: 
                s |= {k} 
        return s 

    def is_prime(self,x): 
        if x not in self.iseq.l: 
            return None 

        i = np.where(self.iseq.l == x)[0][0]
        return len(self.factors[i]) <= 2 

    def primes(self): 
        return set([x for x in self.iseq.l if self.is_prime(x)])