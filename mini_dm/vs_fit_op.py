from .vs_fit_mem import * 
from types import MethodType,FunctionType
from collections import defaultdict

"""
structure that aids in modifying and applying an 
<AffineDelta> instance. 
"""
class VSTransform:

    def __init__(self,ad:AffineDelta):
        assert type(ad) == AffineDelta
        self.ad = ad
        return 
    
    #------------------------ auto-correcting and processing functions 
    
    #------------------------ I/O comparison functions
    
    """
    outputs the two size differences 
    (0) diffabs (absolute size difference) 
    (1) diffsize (size of difference)
    """
    def cmp_ad(self,ad1):
        sz_diff = abs(self.ad.size() - ad1.size())
        diff_ad = (self.ad - ad1).size() 
        return sz_diff,diff_ad 
    
    """
    calculates the difference in expected 
    output value between <AffineDelta> `ad1` 
    and the class <AffineDelta> `ad` with `x` 
    as the input. 
    """
    def diff_ad(self,x,ad1,op_type="all",\
        dfunc=lambda x,x2: np.sum(np.abs(x - x2))): 

        tv = ad1.fit(x)
        return self.ad.expected_diff(tv,x,op_type,\
            dfunc)

    def diff_derivative(self,x,m_delta,a_delta,dfunc=lambda x,x2:x-x2):
        ad1 = AffineDelta(self.m + m_delta,self.a + a_delta,\
            self.ad.ma_order)
        return dfunc(ad1.fit(x),self.ad.fit(x))

    """
    calculates the difference in contribution vectors of input values 
    `x1`,`x2`. 
    """
    def cvec_diff_ad(self,x1,x2,dfunc=lambda x,x2:np.abs(x - x2)):
        return dfunc(self.ad.cvec(x1),self.ad.cvec(x2))

    """
    outputs a function, hypothesis output difference. 
    """
    def to_hypdiff_func(self): 

        def f(x,hyp_y):
            actual =  self.ad.fit(x) 
            return actual - hyp_y 
        return f 
    
    """
    outputs a function, multiple-additive (MA) difference. 
    """
    def to_ma_diff_func(self):
        
        def f(ad):
            assert type(ad) == AffineDelta
            return self.ad - ad 
        return f 


class IOFit:

    def __init__(self,x,y,unknown_func,hypdiff_func,madiff_func):
        assert type(unknown_func) in {type(None),FunctionType,MethodType}

        self.x = x
        self.y = y 
        # function w/ unknown parameter values 
        self.unknownf = unknown_func
        # hypothesis function for x-to-y mappings 
        self.hyp = None 
        # expected/actual output difference b/t `unknownf` and `hyp`
        self.hypdiff_func = hypdiff_func
        # MA difference b/t actual and hypothesis 
        self.madiff_func = madiff_func

        # used for forming hypotheses on MA 
        self.mahm = None # <> 
        self.d,self.cv,self.ma_dim,self.ma_order = \
            None,None,None,None 

    def ranged_process(self,prange): 
        return -1 
    
    """
    CAUTION: no type-checking of arguments 
    """
    def load_mahyp_auxvar(self,d,cv,ma_dim,ma_order):
        self.d = d 
        self.cv = cv 
        self.ma_dim = ma_dim 
        self.ma_order = ma_order 
        return
    
    def init_MAHypMach(self):
        dx = self.process_MAHypMach_on_auxvar()
        ks = sorted(dx.keys())

        xs = [self.x[k] for k in ks] 
        ys = [self.y[k] for k in ks] 
        d = np.array([self.d[k] for k in ks])
        mem_type = "MA"

        mhm = MAHypMach(xs,ys,d,mem_type=mem_type)
        mhm.load_ma_hyp_dict(dx,clear_mem=True)
        self.mahm = mhm 

    def process_MAHypMach_on_auxvar(self): 
        l = len(self.x) 
        dx = dict() 
        for i in range(l): 
            x,y=self.io_sample(i)
            dim0 = vs_dim(x) 
            dim1 = vs_dim(y) 

            if max([dim0,dim1]) != max(self.ma_dim): 
                print("not of same dim.")
                continue 
            
            d = get_vs_element(self.d,i)
            cv = get_vs_element(self.cv,i,cf=lambda x: type(x) == tuple) 
            ad = MAHypMach.io_to_AffineDelta(x,y,d,cv,self.ma_dim,self.ma_order)
            dx[i] = ad 
        return dx 

    
    """
    outputs a description of the x and y types of the class 
    data. The description for each of x,y is 
            u | {unique data element dimensions}
    """
    def type(self):
        xinf = self.stat_count("x")
        yinf = self.stat_count("y")

        x = "u" if type(xinf) == type(None) else \
            set(xinf.keys())
        y = "u" if type(yinf) == type(None) else \
            set(yinf.keys())
        return x,y 
    
    def stat_count(self,seq="x"):
        assert seq in {"x","y"} 
        q = getattr(self,seq) 

        if type(q) == type(None):
            return None
        
        dcount = defaultdict(int) 
        for q_ in q:
            if is_vector(q_):
                dcount[len(q_)] += 1
            else: 
                dcount[0] += 1 
        return dcount 

    def add_value(self,x_,y_):
        self.x.append(x_)

        if type(self.y) != type(None): 
            self.y.append(y_)
        return

    def io_stat(self):
        stat1 = type(self.x) != type(None)
        stat2 = type(self.y) != type(None)

        if stat1 and stat2: 
            assert len(self.x) == len(self.y)
        return stat1,stat2 
    
    def load_hyp(self,hyp,ma_dim): 
        assert type(hyp) == MADHyp
        self.hyp = hyp 
        self.hyp_madim = ma_dim
        self.hyp_proc(self.hyp_madim)

    def hyp_proc(self,ma_dim):
        assert type(self.hyp) == MADHyp

        # convert to <AffineDelta>. 
        ad = self.hyp.solve_into_AffineDelta(ma_dim)

        # convert to <VSTransform>
        vst = VSTransform(ad) 
        self.hfunc = vst.to_hypdiff_func()
        return self.hfunc 
    

    #-------------------------- expected/actual difference functions 

    """
    difference in output values of `x` b/t hypothesis and 
    function `U` (`unknownf`) [expected]: 

        [got w/ hypothesis] - [output from `U`] 
    """
    def hyp_diff(self,x):
        q = self.unknownf(x) 
        return self.hfunc(x,q)
    
    """
    outputs an <AffineDelta> instance ad2, 
    ad2 = (actual MA) - (MA for `ad`). 
    """
    def ma_diff(self,ad): 
        return self.madiff_func(ad) 
    
    def ma_output_diff(self,ad,i,dfunc=lambda x,x2: x - x2): 
        x,y = self.io_sample(i)
        return dfunc(y,ad.fit(x))


    """
    fetches the i'th sample from either `x` or `y`. 
    """
    def io_sample(self,i,default_source = "y"):
        assert default_source in {"y","unknown"} 
        q = self.x[i]

        q2 = None 
        if default_source == "y":
            if type(self.y) != type(None): 
                q2 = self.y[i]
        

        if type(q2) == type(None):
            q2 = self.unknownf(q)
        return (q,q2)
    
    @staticmethod
    def io_pointpair_to_AffineDelta(i0,i1,o0,o1,ma_order):
        d0 = vs_dim(i0)
        d1 = vs_dim(i1)

        assert vs_dim(o0) == vs_dim(o1) 

        M,A = None,None

        dx1 = i1 - i0 
        dx2 = o1 - o0 

        dx3 = safe_div(dx2,dx1)
        M = np.array(dx3)

        # case: multiplication is first 
        if ma_order == 0: 
            t = i0 * M 
            A = o0 - t

            t2 = i1 * M 
            A2 = o1 - t2  
        else: 
            t = np.array(safe_div(o0,M)) 
            A = t - i0

            t = np.array(safe_div(o1,M)) 
            A2 = t - i1  
        return AffineDelta(M,A,ma_order)  