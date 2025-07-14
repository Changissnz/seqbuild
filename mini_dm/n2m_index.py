from morebs2.poly_struct import CEPoly 
from morebs2.matrix_methods import is_valid_range,vector_to_string,string_to_vector
from morebs2.numerical_generator import modulo_in_range,prg_choose_n
from math import floor 
from intigers.extraneous import zero_div0
from collections import defaultdict
from .csrch import * 
import numpy as np 

#--------------------------------------------------------------------------

def indexvalue_map_to_vector(m,sz): 
    l = np.zeros((sz,))
    for k,v in m.items(): 
        l[k] = v 
    return l 

def vector_to_indexvalue_map(v): 
    d = {}
    for (i,v_) in enumerate(v): 
        d[i] = v_
    return d 

def assert_nm(nm): 
    assert type(nm) == tuple and len(nm) == 2 
    assert type(nm[0]) in {int,np.int32,np.int64}
    assert type(nm[1]) in {int,np.int32,np.int64}
    assert min(nm) >= 0 

# NOTE: function is similar to <morebs2.HopPattern> 
def proportional_index__n2m(n,m,ni):
    assert type(n) in {int,np.int32,np.int64}
    assert type(m) in {int,np.int32,np.int64}
    assert n > 0
    assert m > 0 
    assert ni < n
    r = zero_div0(ni,(n -1)) 
    return int(round((m-1) * r))

# NOTE: function is similar to <morebs2.HopPattern> 
def paired_n2m_partition__orderedproportional(n,m):
    assert type(n) in {int,np.int32,np.int64}
    assert type(m) in {int,np.int32,np.int64}
    assert n > 0
    assert m > 0 

    q = n / m 
    n_ = [i for i in range(n)] 
    s = []
    ni = 0
    for i in range(m):
        ni2 = int(round(ni + q))
        
        ni_ = floor(ni)
        if ni2 != ni_:  
            qx = (set(n_[ni_:ni2]),{i})
        else: 
            qx = (set([n_[ni_]]),{i})
        s.append(qx) 
        ni = ni + q
    return s 

#---------------------------------------------------------------------

class N2MIndexMap:

    def __init__(self,nm,n2m_map=set()):
        assert_nm(nm) 
        self.nm = nm 
        self.n2m_map = n2m_map 
        self.check_map() 

    def check_map(self): 
        for m in self.n2m_map:
            m0,m1 = string_to_vector(m[0]),string_to_vector(m[1])
            self.check_element(m0,m1)
        return
    
    def check_element(self,m0,m1):
        assert min(m0) >= 0 and max(m0) < self.nm[0] 
        assert min(m1) >= 0 and max(m1) < self.nm[1] 
    
    def degree_of_mindex(self,j):
        c = 0
        for m in self.n2m_map: 
            m1 = string_to_vector(m[1]) 
            if j in m1: c += 1 
        return c 

    """
    main method 
    """
    def add(self,p): 
        assert type(p) == tuple and len(p) == 2 
        self.check_element(p[0],p[1])
        p0 = vector_to_string(sorted(p[0]))
        p1 = vector_to_string(sorted(p[1]))

        q = (p0,p1) 
        if q in self.n2m_map: return False 
        self.n2m_map |= {q}
        return True


    def mindex_degree_map(self):
        d = {} 
        for j in range(self.nm[1]): 
            d[j] = self.degree_of_mindex(j) 
        return d 

class N2MIndexMapGen: 

    def __init__(self,nm,index_degree_range,\
        nset_size_range,mset_size_range,prg,\
        index_degree_is_geq:bool=False): 
        assert_nm(nm) 
        assert is_valid_range(index_degree_range,is_int=True,inclusive=False)
        assert is_valid_range(nset_size_range,is_int=True,inclusive=False)
        assert is_valid_range(mset_size_range,is_int=True,inclusive=False)
        assert index_degree_range[0] >= 0 
        assert nset_size_range[0] > 0
        assert mset_size_range[0] > 0 
        assert type(prg) in {MethodType,FunctionType}
        assert type(index_degree_is_geq) == bool 

        self.nm = nm 
        self.id_range = index_degree_range
        self.ns_range = nset_size_range
        self.ms_range = mset_size_range
        self.prg = prg 
        self.is_geq = index_degree_is_geq

        self.n2m_imap = N2MIndexMap(self.nm) 
        self.target_index_degree_map = None
        self.dmap = defaultdict(int)
        self.preproc() 

    def preproc(self):
         self.target_index_degree_map = dict() 
         for i in range(self.nm[1]):
            self.target_index_degree_map[i] = \
                modulo_in_range(self.prg(),self.id_range)
            self.dmap[i] = 0 
         return
    
    """
    main function 
    """
    def make(self,attempts_per_relation:int = 10 ** 5): 
        stat = True 
        while stat: 
            r = self.one_new_relation(attempts_per_relation) 
            if type(r) == type(None):
                stat = False 
        return
    
    def one_new_relation(self,num_attempts = 10** 5):
        # mset is from available
        ai = self.available_indices()

        if len(ai) == 0: 
            return None 
        
        msize = int(modulo_in_range(self.prg(),self.ms_range)) 
        msize = min([msize,len(ai)])
        mset_index = prg_choose_n([i for i in range(len(ai))],\
            msize,self.prg,is_unique_picker=True) 
        mset = sorted([ai[mi] for mi in mset_index])
        mset = sorted(mset)

        # nset is from static even distribution 
        nsize = int(modulo_in_range(self.prg(),self.ns_range)) 
        sx = [i for i in range(self.nm[0])]
        nset = prg_choose_n(sx,nsize,self.prg,is_unique_picker=True)
        nset = sorted(nset) 

        # try adding the generated (nset,mset)
        stat = self.add_relation(nset,mset)
        if stat: 
            return (nset,mset) 

        def cfunc(ms):
            v1 = vector_to_string(nset) 
            v2 = vector_to_string(ms) 
            if (v1,v2) in self.n2m_imap.n2m_map: 
                return False 
            return True 

        # iteration through combinations in search 
        # of new (nset,mset) 
        ai2 = []
        q = set([i for i in range(msize)]) 
        for mi in mset_index:
            ai2.append(ai[mi]) 
            q -= {mi}
        q = sorted(q)
        ai2.extend(q) 
        ai = ai2 

        ccs = ClosestCombinationSearch(msize,ai,cfunc,num_attempts = num_attempts)
        mset,x = ccs.find() 

        if x: 
            stat = self.add_relation(nset,mset)
            return (nset,mset) 
        return None 
    
    def add_relation(self,nset,mset):
        stat = self.n2m_imap.add((nset,mset)) 
        if stat:
            for j in mset: 
                self.dmap[j] += 1  
        return stat 
    
    def available_indices(self):
        diffvec = []

        for k,v in self.dmap.items():
            v2 = self.target_index_degree_map[k]
            q = v2 - v 
            if not self.is_geq:
                if q <= 0:
                    continue 

            diffvec.append((k,q))
        diffvec = sorted(diffvec,key=lambda x:x[1],reverse=True) 
        return [d[0] for d in diffvec]

"""
def f(vec,index):
    return 
"""
class N2MVectorFunction:

    def __init__(self):
        return 