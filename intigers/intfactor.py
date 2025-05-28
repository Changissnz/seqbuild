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

# TODO: relocate 
def std_invert_map(m):
    assert type(m) in {dict,defaultdict}

    q = {}
    for k,v in m.items():
        if v in q:
            q[v] |= {k}
        else:
            q[v] = {k} 
    return q

def numberdict_subtraction(d1,d2):
    K = set(d1.keys()) | set(d2.keys()) 
    
    d3 = defaultdict(int) 
    for k in K: 
        x1 = d1[k]
        x2 = d2[k]
        d3[k] = x1 - x2
    return d3 


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
        # 
        self.factor_count = None 
        self.cfd_map = None 
        self.max_cofactor_degree = None 

    """
    main method 
    """
    def factor_count_(self): 
        self.factor_count = defaultdict(int) 
        self.max_cofactor_degree = 0 
        for mx in self.factors: 
            for q in mx: 
                self.factor_count[q] += 1 
                self.max_cofactor_degree = self.max_cofactor_degree if \
                    self.factor_count[q] < self.max_cofactor_degree else \
                    self.factor_count[q] 
        self.cfd_map_() 

    #--------------------------- updated methods related to removing elements 
    #--------------------------- from sequence  

    # TODO: test 
    def remove_seq_elements(self,elements): 
        element_indices = self.iseq.element_indices(elements) 
        self.iseq.remove_element_indices(element_indices)

        affected_factors = set() 
        factor_count_delta = defaultdict(int) # negative 
        for x in element_indices:
            for x_ in self.factors[x]: 
                factor_count_delta[x_] += 1
        self.update_factor_count(factor_count_delta)

        # NOTE: ineff 
        self.cfd_map_() 

        self.factors = [f for (i,f) in enumerate(self.factors) if \
            i not in element_indices] 

    def update_factor_count(self,factor_count_delta):  
        self.factor_count = numberdict_subtraction(self.factor_count,\
            factor_count_delta) 
        return

    #--------------------------- methods related to cofactor degree 

    def cfd_map_(self):
        assert type(self.factor_count) != type(None) 
        self.cfd_map = std_invert_map(self.factor_count)

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

# TODO: test 
class FactorClassifier: 

    def __init__(self,factor2class=dict()):
        assert type(factor2class) == dict 
        self.class_ctr = 0 
        if len(factor2class) > 0: 
            V = set(factor2class.values()) 
            assert min(V) == 0 and max(V) == len(V) - 1
            self.class_ctr = len(V) 
        self.f2c = factor2class
        self.f2c_inverted = dict([(v,k) for (k,v) in self.f2c.items()])

    def add_class(self,k): 
        assert k not in self.f2c
        self.f2c[k] = self.class_ctr 
        self.f2c_inverted[self.class_ctr] = k 
        self.class_ctr = self.class_ctr + 1 

    def classify(self,q):
        for i in range(self.class_ctr): 
            q2 = self.f2c_inverted[i]
            stat = q / q2 == q // q2 
            if stat: return i 
        return -1 