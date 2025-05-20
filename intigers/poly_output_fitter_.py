from .process_seq import * 
from copy import deepcopy
from morebs2 import aprng_gauge

"""
Finds an n'th degree polynomial P with all coefficients 
variable except for the n'th power, set to argument 
`coeff`, s.t. P(x1) = P(x2). 

If `prng` is not None, during the search process, algorithm 
chooses a pseudo-random candidate coefficient for every power 
except for 1 and n.
"""
class PolyOutputFitterVar2:

    def __init__(self,n,x1,x2,coeff=1,prng=None):
        for q in [n,x1,x2]: 
            assert type(q) in {int,np.int64} 
        assert n > 1 
        assert x1 != x2 
        self.n = n 
        self.x1 = np.int64(x1)
        self.x2 = np.int64(x2) 
        self.ref = None 
        self.set_poly(coeff) 
        self.prng = prng 

        self.powerdiff_vec() 
        q = np.zeros((2,),dtype=np.int64) 
        self.running_diff = deepcopy(q) 
        q[0] = self.apply(self.x1) 
        q[1] = self.apply(self.x2) 
        self.update_runningdiff(q) 

        self.stat  = True 
        if not self.is_solvable(): 
            print("[??] cannot compute...") 
            self.stat = False 

    def set_poly(self,coeff:int=1): 
        self.poly = np.zeros((self.n,),dtype=np.int64) 
        self.poly[0] = coeff
        self.index = 1  

    def powerdiff_vec(self): 
        self.pdvec = np.zeros((self.n,),dtype=np.int64)

        self.ref = self.x1 if self.x1 < self.x2 else self.x2 

        for i in range(1,self.n+1): 
            px1 = self.x1 ** i 
            px2 = self.x2 ** i
            self.pdvec[self.n - i] = abs(px1-px2) 

    def apply(self,x): 
        q = 0 
        for i in range(1,self.n+1): 
            j = self.n - i 
            q += self.poly[j] * x ** i  
        return q 

    def next_coeff_(self,i): 
        q = self.running_diff / self.pdvec[i]
        if q[0] == 0.0: return -q[1]
        return q[0]

    def next_coeff(self,i): 
        q = self.next_coeff_(i) 
        x = floor(q) 

        if self.n - 1 == i: 
            return np.int64(x)

        x2 = 1 if x > 0 else - 1 
        if type(self.prng) == type(None): 
            return np.int64(x + x2) 
        l = [x + x2]  
        c = 1 
        s2,s3 = l[-1],l[-1]
        target = 1 if x > 0 else -1 

        while s2 != target and s3 != target: 
            s2,s3 = s - c,s + c 
            l.extend([s2,s3])
            c += 1
        j = self.prng() % len(target)
        return np.int64(l[j])

    def solve(self): 
        while self.index < self.n: 
            q = self.next_coeff(self.index)
            j = self.n - self.index # ?  
            new_diff = np.array([self.x1**j,self.x2 ** j],dtype=np.int64) 
            new_diff = q * new_diff 
            self.update_runningdiff(new_diff) 
            self.poly[self.index] = q 
            self.index += 1 
        return "finnamoto" 

    def resolve(self,pwr,new_coeff): 
        assert pwr > 1 
        j = self.n - pwr 
        assert j >= 0 

        new_coeffvec = np.zeros((self.n,),dtype=np.int64) 
        new_coeffvec[:j] = self.poly[:j] 
        self.poly = new_coeffvec 
        self.poly[j] = new_coeff
        self.running_diff = np.zeros((2,),dtype=np.int64)  

        y1 = self.apply(self.x1)
        y2 = self.apply(self.x2) 
        q = np.array([y1,y2],dtype=np.int64)
        self.update_runningdiff(q) 
        self.index = j + 1 
        self.solve()
        return

    def is_solvable(self): 
        stat = False 
        targ = self.pdvec[0] 

        for i in range(1,self.n): 
            if targ // self.pdvec[i] == targ / self.pdvec[i]: 
                return True 
        return False 

    def update_runningdiff(self,new_diff): 
        self.running_diff = self.running_diff + new_diff 
        q = np.min(self.running_diff) 
        self.running_diff = self.running_diff - q 
