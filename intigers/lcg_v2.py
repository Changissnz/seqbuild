from morebs2.graph_basics import * 
from morebs2.numerical_generator import modulo_in_range,prg__LCG
from .extraneous import to_trinary_relation
from .seq_struct import * 
from mini_dm.ag_ext import * 

def output_LCG_matrix(lcg,shape,range_outputter):
    assert type(shape) == tuple and len(shape) == 2
    x = np.zeros((shape[0],shape[1]))

    for i in range(shape[0]):
        for j in range(shape[1]):
            x[i,j] = modulo_in_range(lcg(),range_outputter())
    return x 

# no sub-cycle -> continuous
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


"""
Measures categorical entropy over rows or columns of matrix `m`. The variable 
`sl_info` is for the `seg_length` parameter of <APRNGGaugeV2.std_cat_entropy>. 
It is one of five types: None, int (denumerator for max pairwise distance), 
vector<int> (denumerator for each row or column), float (literal), and vector<float>. 
"""
def APRNGGaugeV2__matrix_cat_entropy(m,franges,is_rowwise:bool=True,is_local_frange:bool=True,\
    sl_info = None,count_type="absdiff",round_depth:int=5):
    assert is_2dmatrix(m)
    m_ = m.T if not is_rowwise else m 
    f = None 
    i = None  

    if type(franges) == type(None): 
        if is_local_frange: 
            franges = []
            for x in m_: 
                franges.append((np.min(x),np.max(x))) 
            franges = np.array(franges,dtype=np.int32) 
        else: 
            franges = (np.min(m_),np.max(m_))

    if is_bounds_vector(franges): 
        assert is_proper_bounds_vector(franges)
        assert franges.shape[0] == len(m_)
        f,i = franges[0],0 
    else: 
        assert type(franges) == tuple
        assert len(franges) == 2 
        f = franges 

    j = 0 if is_vector(sl_info) else None
    sl = sl_info  

    ag = APRNGGaugeV2(None,f,pradius=5)
    lx = [] 
    for x in m_: 
        iseq = IntSeq(x) 

        if type(i) != type(None): 
            f = franges[i] 
            ag.reload_var("frange",tuple(f) )

        if type(j) != type(None):
            sl = sl_info[j] 
            assert sl > 0

        sl_ = sl 
        if type(sl_) in {int,np.int32,np.int64}:
            sl_ = (f[1] - f[0]) / sl 

        ce = ag.std_cat_entropy(iseq,seg_length=sl_,start_value=f[0],\
            count_type=count_type)
        lx.append(ce) 

        if type(i) != type(None): 
            i += 1 

        if type(j) != type(None): j += 1 

    return np.round(np.array(lx),round_depth)


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

"""
an extension of the standard linear congruential generator (LCG). 
In addition to outputting values by the same mechanism as the original 
LCG, the <LCGV2> is able to discover sub-cycles in its input-output 
map. Algorithm is able to discover sub-cycles in the LCG by treating 
the LCG output values in the manner of nodes for a connected graph. 
Algorithm is based on graph decomposition. For instantiation of this 
class, if variable `preproc_gd` is set to True, class instance conducts 
a structural search for sub-cycles via graph decomposition algorithms 
(see method `gd_preproc`).
"""
class LCGV2:

    def __init__(self,start,m,a,n0,n1,sc_size:int,preproc_gd:bool=False):
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

        if preproc_gd: 
            self.gd_preproc() 

    """
    main method 
    """
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

        if self.sc_size > 0:
            if self.s == s_: 
                d = 0
            self.log_sc.append(d) 
            rx = len(self.log_sc) - self.sc_size 
            while rx > 0:
                self.log_sc.pop(0) 
                rx -= 1
        
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

    """
    main method #2 
    """
    def gd_preproc(self):
        self.io_map()
        self.io_map_partition()
        self.io_map_summary()

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

        