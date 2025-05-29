from morebs2.numerical_extras import * 
from morebs2.graph_basics import flatten_setseq
from morebs2.poly_factor import median_sort,select_median_in_sequence 
from morebs2.globalls import std_invert_map,numberdict_op,equal_intdicts
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

###

"""
set operations of multiples for an integer sequence
"""
class ISFactorSetOps: 

    def __init__(self,l,int_limit=DEFAULT_INT_MAX_THRESHOLD,str_mode_full:bool=True):  
        assert len(l) >= 2
        assert abs(min(l)) <= DEFAULT_INT_MAX_THRESHOLD and \
            abs(max(l)) <= DEFAULT_INT_MAX_THRESHOLD
        self.iseq = IntSeq(l) 
        # element in iseq -> factors 
        self.factors = factors_of_seq(self.iseq) 
        # factor -> degree 
        self.factor_count = None 
        # degree -> factor set  
        self.cfd_map = None 
        self.max_cofactor_degree = None 
        self.str_mode_full=str_mode_full

    """
    For the keys (elements of `iseq`) provided or all elements 
    of `iseq` if `pkeys` is None, sorts the keys by their co-factor 
    degree with respect to `pkeys`. 

    Outputs a sequence of pairs 
    (factor, co-factor degree)
    """
    def dsort(self,pkeys=None):
        dx = None 
        if type(pkeys) == type(None): 
            dx = self.factor_count
        else: 
            s = set() 
            element_indices = self.iseq.element_indices(pkeys) 
            for ei in element_indices: 
                s |= self.factors[ei] 
            dx = defaultdict(int) 
            for ei in element_indices:
                q = self.factors[ei]
                dx2 = dict([[q_,1] for q_ in q]) 
                dx = numberdict_op(dx,dx2,add)
        dx = sorted([[k,v] for (k,v) in dx.items()],key=lambda x:x[1])
        return dx 
        
    def median(self,pkeys=None,r=0.5): 
        qr = self.dsort(pkeys) 
        qr = [x[0] for x in qr]
        return select_median_in_sequence(qr,m=r)

    def median_sort(self,pkeys=None,r=0.5,fullpair_sequence:bool=True): 

        # get the degree to size map 
        qr = self.dsort(pkeys) 
        if not fullpair_sequence: 
            qr = [x[0] for x in qr]
        x = median_sort(qr,m=r)
        return x 

    def factor_to_keys(self,f): 
        ks = set()
        for (i,k) in enumerate(self.iseq.l): 
            if f in self.factors[i]: 
                ks |= {k} 
        return ks 

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

    def remove_seq_elements(self,elements): 
        element_indices = self.iseq.element_indices(elements) 
        self.iseq.remove_element_indices(element_indices)

        dx = self.factorcount_for_elementindices(element_indices) 
        self.update_factor_count(dx) 

        # NOTE: ineff 
        self.cfd_map_() 

        self.factors = [f for (i,f) in enumerate(self.factors) if \
            i not in element_indices] 

    def factorcount_for_elementindices(self,element_indices): 
        d = defaultdict(int)
        for ei in element_indices:
            for x_ in self.factors[ei]: 
                d[x_] += 1
        return d  

    # TODO: test 
    def factorset_for_elements(self,elements):
        s = set() 
        ei = self.iseq.element_indices(elements) 
        for ei_ in ei: s |= self.factors[ei_] 
        return s 

    def update_factor_count(self,factor_count_delta):  
        self.factor_count = numberdict_op(self.factor_count,\
            factor_count_delta,f=sub) 
        ks = set() 
        for k,v in self.factor_count.items(): 
            if v <= 0: ks |= {k} 
        for k in ks: del self.factor_count[k] 
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

    def coprimes_of(self,i1):
        coprimes = set()
        for l_ in self.iseq.l: 
            if l_ == i1: continue 
            stat = self.are_coprimes(i1,l_)
            if stat: coprimes |= {l_}
        return coprimes 

    def are_coprimes(self,i1,i2):
        assert i1 in self.iseq.l and i2 in self.iseq.l 
        index1 = np.where(self.iseq.l == i1)[0][0] 
        index2 = np.where(self.iseq.l == i2)[0][0]

        q = self.factors[index1].intersection(self.factors[index2])
        return q == {1} or q == {} 

    # TODO: test 
    def multiples_of(self,f): 
        s = set() 
        for (i,m) in enumerate(self.factor_count): 
            if f in m: 
                s |= {self.iseq[i]} 
        return s 

    def __str__(self): 
        s = ""
        for (i,f) in enumerate(self.factors): 
            s += "element {}\n".format(self.iseq[i])
            s += "factors\n"
            s += str(f) + "\n\n"

        if self.str_mode_full: 
            s += "\t\tdegree-to-factor map" + "\n"
            for (k,v) in self.cfd_map.items():
                s += "degree: {}\n".format(k)
                s += "factors\n"
                s += str(v) + "\n\n"
        return s 
