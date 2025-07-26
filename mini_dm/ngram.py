import numpy as np 
from morebs2.matrix_methods import is_vector

#------------------------- subvector operations 

def subvec(l,start_index,length):
    assert is_vector(l) or type(l) == list 
    assert 0 <= start_index < len(l)
    assert 0 < length <= len(l) 

    q = list(l[start_index:start_index+length])

    l1 = len(q) 
    excess = length - l1
    q2 = [] 
    if excess > 0: 
        q2 = list(l[:excess]) 
    return q + q2 

#------------------------- arg. check 

def assert_nm(nm): 
    assert type(nm) == tuple and len(nm) == 2 
    assert type(nm[0]) in {int,np.int32,np.int64}
    assert type(nm[1]) in {int,np.int32,np.int64}
    assert min(nm) >= 0 

"""
class is an n-grammer over sequence `l`, either a list 
or a vector (1-dimensional np.array). Outputs contiguous
subsequences from `l` of length `n`, reshaped into the 
two-dimensional shape `dim2`. 
"""
class NGrammer2D:

    def __init__(self,l,n,dim2):
        assert type(l) == list or is_vector(l) 
        assert len(l) > 1 
        self.l = l 

        self.check_length(n)
        self.n = n 

        if type(dim2) == type(None): 
            dim2 = (1,n) 
        self.check_dim2(dim2) 
        self.dim2 = dim2 

        self.ref_index = 0 

    def check_length(self,n):
        assert type(n) in {int,np.int32,np.int64}
        assert n > 0 

    def check_dim2(self,dim2): 
        assert_nm(dim2) 
        assert dim2[0] * dim2[1] == self.n

    def one_cycle(self,ref_index=0): 
        assert type(ref_index) in {int,np.int32,np.int64} 
        assert 0 <= ref_index < self.n
        self.ref_index = ref_index 

        c = 0
        while self.ref_index != ref_index or c < 1:
            if self.ref_index == ref_index: 
                c += 1 

                if c == 1:
                    yield self.__next__() 
                else: 
                    continue 
            
            yield self.__next__()
        return 

    def __next__(self):
        sv = subvec(self.l,self.ref_index,self.n)
        self.ref_index = (self.ref_index + 1) % len(self.l) 

        return np.array(sv).reshape(self.dim2) 

    def reset_spec(self,new_n=None,new_dim2=None,new_ref_index=0):
        if type(new_n) != type(None):
            self.check_length(new_n)
            self.n = new_n

        if type(new_dim2) != type(None): 
            self.check_dim2(new_dim2) 

        new_ref_index %= len(self.l) 
        self.ref_index = new_ref_index