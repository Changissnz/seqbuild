"""
file contains methods to aid in calculating analytical values on integer vectors
"""
from mini_dm.iseq import * 

def affine_fit_for_pair__multiple(i1,i2): 
    if i1 == 0: return np.nan 
    return np.int32(ceil(i2/i1))

"""
Outputs a sequence of (multiple,additive) pairs for every 
contiguous pair in integer sequence `l`. 

The maximum multiple m to consider for each contiguous pair 
is: 

MAX abs([ceil(l[i]/l[i-1])]); |l| > i >= 1. 

If `exclude_neg` set to False, considers the span [-m,m]\\{0}.  
"""
class AffineFitCandidates: 

    def __init__(self,l,exclude_neg:bool=True):
        self.l = l 
        self.exclude_neg = exclude_neg
        self.load()
        self.m = self.max_multiple()
        self.i = 1 

    def load(self):
        self.l = np.array(self.l,"int32")
        assert len(self.l.shape) == 1 and len(self.l) > 1 
        return 

    def max_multiple(self): 
        q = 0 
        for i in range(1,len(self.l)): 
            m = affine_fit_for_pair__multiple(self.l[i-1],self.l[i]) 
            if np.isnan(m): 
                continue 

            lm = [q,m]
            j = 0 if abs(q) > abs(m) else 1  
            q = lm[j]
        return abs(q) 
    
    def candidates_at_index(self,i,exclude_neg:bool=True): 
        assert i >= 1 and i < len(self.l) 
        self.exclude_neg = exclude_neg

        ds = [] 
        v1,v2 = self.l[i-1],self.l[i]

        if self.m == 0: 
            ds.append((i,int(v2 - (v1 * 1)))) 
            return ds 

        if not self.exclude_neg: 
            mx = -self.m 
            for i in range(mx,0): 
                ds.append((i,int(v2 - (v1 * i))))
        for i in range(1,self.m+1): 
            ds.append((i,int(v2 - (v1 * i))))
        
        return ds 

    def next_candidate_set(self,exclude_neg:bool=True): 
        if self.i >= len(self.l): return None 
        q = self.candidates_at_index(self.i,exclude_neg) 
        self.i += 1
        return q 