from morebs2.graph_basics import * 
from .ag_ext import * 
from .tvec import * 

# no sub-cycle is continuous
CYCLE_CATEGORIES = {"closed","sub-cycle"}

def io_map__signchange_count(m):
    p,n = 0,0
    for k,v in m.items(): 
        if v > k: p += 1
        elif v < k: n += 1
        else: pass 
    return p,n 

def travel_io_map_till_repeat(m,k):
    q = [k]
    stat = True 
    while stat: 
        v = m[k]
        if v in q: break 
        q.append(v) 
        k = v 
    return q 

class TFunc__NoRepeats: 

    def __init__(self,s:set):
        assert type(s) == set 
        self.l = s

    def output(self,i,i2):
        stat = not (i in self.l) 
        self.l |= {i} 
        return stat 

class CycleDescriptor:

    def __init__(self):
        self.d = defaultdict(None)

    def __str__(self):
        s = "closed: " + str(self.d["closed"]) + "\n" 
        s += "sub-cycle heads: " + \
            str(self.d["sub-cycle"]) + "\n"
        return s 

    def update(self,k,v): 
        assert k in CYCLE_CATEGORIES
        if k == "closed":
            assert type(v) == bool 
        else:
            assert type(v) in {type(None),set}

        self.d[k] = v 
        return

    def is_closed(self):
        return self.d["closed"] 

    def is_continuous(self):
        return type(self.d["sub-cycle"]) == type(None)

class LCGV2:

    def __init__(self,start,m,a,n0,n1,sc_size:int):
        assert n0 < n1
        assert not (m == 0 and a == 0)
        assert type(sc_size) == int 
        assert sc_size >= 0 

        self.s = start
        self.s_ = start  
        self.m = m 
        self.a = a 
        self.r = [n0,n1]

        self.sc_size = sc_size 
        self.log_sc = [] 

        self.cycle = []

        self.fired = False
        self.cycled = False
        self.dir = None
        self.map_io = dict()
        self.gd = None 
        self.cycle_descriptors = [] 

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

    def io_map(self):
        self.map_io.clear()
        for x in range(self.r[0],self.r[1]): 
            y = modulo_in_range(x * self.m + \
                self.a,self.r) 
            self.map_io[x] = y 

    def io_map_partition(self): 
        qx = defaultdict(set)
        for k,v in self.map_io.items():
            qx[k] = set([v]) 

        self.gd = GraphComponentDecomposition(qx) 
        self.gd.decompose()

    def io_map_summary(self):
        for i in range(len(self.gd.components)):
            cd = self.component_index_summary(i)
            self.cycle_descriptors.append(cd) 

    def component_index_summary(self,i):
        q = self.gd.components[i] 
        q = flatten_setseq(q) 

        is_closed = True
        sub_cycle = set()
        for q_ in q:
            p = travel_io_map_till_repeat(self.map_io,q_)
            px = set(p)

            if px != q: 
                sub_cycle |= {q_} 

            if not is_closed: continue 

            if not px.issubset(q): 
                is_closed = False

        if len(sub_cycle) == 0: 
            sub_cycle = None

        cd = CycleDescriptor()
        cd.update("closed",is_closed) 
        cd.update("sub-cycle",sub_cycle)
        return cd 

        