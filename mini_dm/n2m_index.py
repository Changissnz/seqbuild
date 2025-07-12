from morebs2.poly_struct import CEPoly 
from morebs2.matrix_methods import is_valid_range,vector_to_string
from morebs2.numerical_generator import modulo_in_range
from math import floor 
from intigers.extraneous import zero_div0
from types import MethodType,FunctionType
from collections import defaultdict
import numpy as np 

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
        return -1 
    
    def degree_of_mindex(self,j):
#def prg_choose_n(Q,n,prg,is_unique_picker:bool=False):

        return -1 

    def add(self,p): 
        assert type(p) == tuple and len(p) == 2 
        p0 = vector_to_string(sorted(p[0]))
        p1 = vector_to_string(sorted(p[1]))

        q = (p0,p1) 
        if q in self.n2m_map: return False 
        self.n2m_map |= {q}
        return True


    def mindex_degree_map(self):
        return -1 

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
    
    def one_new_relation(self):
        return -1 

    def add_relation(self,nset,mset):
        stat = self.n2m_imap.add((nset,mset)) 
        if stat:
            for j in mset: 
                self.dmap[j] += 1  
            return  
        
        return
    
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