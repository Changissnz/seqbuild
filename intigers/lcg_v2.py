from .ag_ext import * 

class TFunc__NoRepeats: 

    def __init__(self):
        self.l = set() 

    def output(self,i,i2):
        stat = not (i in self.l) 
        self.l |= {i} 
        return stat 

# TODO: test
class LCGV2:

    def __init__(self,start,m,a,n0,n1,log_sc:int):
        assert n0 < n1
        assert not (m == 0 and a == 0)
        assert type(log_sc) == int 
        assert log_sc >= 0 

        self.s = start 
        self.m = m 
        self.a = a 
        self.r = [n0,n1]
        self.log_sc = [] 

        self.cycle = []

        self.fired = False
        self.cycled = False
        self.dir = None

    def __next__(self):
        if not self.fired:
            self.fired = not self.fired
            self.cycle.append(self.s)
            return self.s 

        s_ = self.s * self.m + self.a
        s_ = modulo_in_range(s_,self.r) 
        d = 1 if s_ > self.s else -1 

        if type(self.dir) != type(None): 
            if d != self.dir and self.cycled == False: 
                self.cycled = True
            else: 
                self.cycled = False
        else: 
            self.dir = d
        
        self.s = s_ 

        if self.cycled:
            self.cycle.clear()
        self.cycle.append(s_)
        return self.s 
