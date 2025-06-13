from .ag_ext import * 

class TFunc__NoRepeats: 

    def __init__(self,s:set):
        assert type(s) == set 
        self.l = s

    def output(self,i,i2):
        stat = not (i in self.l) 
        self.l |= {i} 
        return stat 

class LCGV2:

    def __init__(self,start,m,a,n0,n1,log_sc:int):
        assert n0 < n1
        assert not (m == 0 and a == 0)
        assert type(log_sc) == int 
        assert log_sc >= 0 

        self.s = start
        self.s_ = start  
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

    #-------------- functions to discover LCG attributes 

    def gen_dir(self): 
        if self.m >= 0 and self.a >= 0: 
            return 1
        elif self.m <= 0 and self.a <= 0:
            return -1
        return 0 

    def init_dir(self):
        q = self.s_
        qx = q * self.m + self.a 

        if qx > q:
            return 1
        elif qx < q: 
            return -1 
        return 0

    def coverage(self): 
        q = self.s 
        d,c,f = self.dir,self.cycled,self.fired 
        self.reset_vars()

        fx = self.__next__ 
        tf = TFunc__NoRepeats({self.s}) 
        ag = APRNGGauge(fx,self.r,0.5)
        cov,_ = ag.measure_cycle(max_size=self.r[1]-self.r[0],\
            term_func=tf.output)            
        self.s,self.dir,self.cycled,self.fired \
            = q,d,c,f 
        return cov 

    def reset_vars(self):
        self.s = self.s_  
        self.dir,self.cycled = None,False 
        self.fired = False 

    def cycle_length(self): 
        q = self.s
        d,c = self.dir,self.cycled 
        self.reset_vars()
        ct = 0 
        while not self.cycled: 
            self.__next__()
            ct += 1
        self.s,self.dir,self.cycled = q,d,c
        return ct 
        