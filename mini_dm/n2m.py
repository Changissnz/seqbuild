from morebs2.matrix_methods import is_number,vector_to_string
from intigers.extraneous import to_trinary_relation
from collections import defaultdict
import numpy as np 

class N2MAutocorrelator:

    def __init__(self,nm):
        self.nm = None 
        self.set_nm(nm) 
        # subset of input indices -> 
        # sign-change vector for output -> 
        # number of times of occurrence given 
        #     sample inputs 
        self.ftable = defaultdict(None)
        return 
    
    def set_nm(self,nm): 
        assert type(nm) == tuple and len(nm) == 2 
        assert type(nm[0]) in {int,np.int32,np.int64}
        assert type(nm[1]) in {int,np.int32,np.int64}
        assert min(nm) >= 0 
        self.nm = nm 

    """
    adds a pair of (x_i,e_i) to memory; e_i the error term for input x_i. 
    """
    def add(self,x0,x1,e0,e1):

        # get the trinary vector for 
        r0 = to_trinary_relation(x1,x0)
        r1 = to_trinary_relation(e1,e0)  

        s0 = vector_to_string(r0,int)
        s1 = vector_to_string(r1,int)

        if s0 not in self.ftable: 
            self.ftable[s0] = defaultdict(int) 
        self.ftable[s0][s1] += 1
        return

"""
def f(vec,index):
    return 
"""
class N2MVectorFunction:

    def __init__(self):
        return 