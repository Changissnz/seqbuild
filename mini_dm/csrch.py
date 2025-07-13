#from morebs2.aprng_gauge import * 
from itertools import combinations
from morebs2.matrix_methods import is_vector
from types import MethodType,FunctionType
import numpy as np

"""
a "one-searcher" algorithm. Given a sequence 
S, with the first k elements of S as S' being 
the target combination, algorithm terminates 
search upon reaching a combination C (ideally, 
this would be S) that satisfies conditional 
function `cfunc`. The elements of S are ordered 
in decreasing priority for being considered a 
candidate for the "best" combination in the search. 

A visualization: 

s_1 s_2 ... s_k |||||||| s_(k+1) ... s_(|S|). 

The tail of the sequence, starting from s_(k+1), 
is de-prioritized elements that may be in the 
"best" combination (output value) if and only if
the head of the sequence does not output True from 
`cfunc`. 
"""
class ClosestCombinationSearch:

    def __init__(self,k,seq,cfunc,num_attempts = 10 ** 5):
        assert type(k) in {int,np.int32,np.int64}
        assert k > 0
        assert is_vector(seq) or type(seq) == list
        assert k <= len(seq) 
        assert type(cfunc) in {FunctionType,MethodType}
        self.k = k 
        self.seq = seq 
        self.cfunc = cfunc 
        self.num_attempts = num_attempts 
        self.num_attempts_ = num_attempts
        self.bx = [self.one_C(self.k)]

        self.j = 0
        self.l = 0 
        self.l_term = False
        return 
    
    def find(self): 
        stat = True
        q,stat3 = None,False 
        while stat: 
            if self.l_term: 
                stat = False 
                continue 

            q,stat2 = next(self) 

            if type(q) == type(None) and \
                not stat2:  
                continue 

            if stat2: 
                stat3 = stat2
                stat = False 
                continue 
        return q,stat3 
        
    def __next__(self):
        if self.num_attempts_ <= 0: 
            return None,False   

        if len(self.bx) == 0: 
            if self.l_term: 
                return None,False   
            
            self.l += 1
            bs1 = self.one_C(self.k + self.l)
            bs2 = self.one_C(self.k - self.l)

            if type(bs1) != type(None): 
                self.bx.append(bs1)
            if type(bs2) != type(None): 
                self.bx.append(bs2) 

            if len(self.bx) == 0: 
                self.l_term = True 
            return None, False   
        
        self.j = self.j % len(self.bx) 
        bx_ = self.bx[self.j]

        try: 
            qi = next(bx_)
        except: 
            self.bx.pop(self.j)
            return None,False

        vx = list(qi)  
        self.num_attempts_ -= 1 
        stat = self.cfunc(vx)
        self.j = (self.j + 1) % len(self.bx) 
        return vx,stat 
    
    def one_C(self,sz):
        if sz >= len(self.seq): return None
        if sz <= 0: return None 
        return combinations(self.seq, sz)