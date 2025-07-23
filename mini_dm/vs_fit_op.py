from .vs_fit_mem import * 
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
        diff_ad0 = (self.ad - ad1).size() 
        return sz_diff,diff_ad0 
    
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


"""
container that holds an input-output pair of sequences `x` and `y`. 
Either the input `x` or output `y` can be unknown. 

NOTE: 
    implementation for that is lacking as of this point in development. 

The `unknown_func` is usually the target function (to be found) 
that outputs the correct y_i for an input x_i. Both the 
`hypdiff_func` and `madiff_func` can be used as metrical functions 
for solution improvement (reduction in error in x-to-y mapping).  
The `madiff_func` is used specifically for measuring the parametric 
differences between two affine functions, while `hypdiff_func` 
is applicable for the general case of any function intended to 
fit `x` to `y`'. Initializing the hypothesis container via a 
<HypMach> instance can proceed by one of two ways:
1. via the `init_HypMach` method, which uses the 
   `process_HypMach_on_auxvar` method to form hypotheses 
   based on the auxiliary variables `d`, `cv`, `ma_dim`, and 
   `ma_order` provided by program user. Every (x_i,y_i) pair
    is associated with one affine function as a hypothesis. 
2. via the `init_HypMach_v2` method, which uses the
   `superpart_by_AffineDelta` method to form hypotheses. 
   The super-partitioning is done by the condition of being 
   fitted by a common affine function for each of the sets. 

Container is integrated with affine functions as the primary 
function form to fit the points. As for other function forms, 
those have to be user-specified. 
"""
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
    
    def init_HypMach(self):
        dx = self.process_HypMach_on_auxvar()
        ks = sorted(dx.keys())

        xs = [self.x[k] for k in ks] 
        ys = [self.y[k] for k in ks] 
        d = np.array([self.d[k] for k in ks])
        mem_type = "MA"

        mhm = HypMach(xs,ys,d,mem_type=mem_type)
        mhm.load_ma_hyp_dict(dx,clear_mem=True)
        self.mahm = mhm 

    def init_HypMach_v2(self,ma_order,index0_vec=None,\
        index1_function=None): 
        hm = self.superpart_by_AffineDelta(ma_order,\
            index0_vec=index0_vec,index1_function=index1_function)
        
        mhm = HypMach(self.x,self.y,None,mem_type="MA") 


        mhm.mhm = hm 
        self.mahm = mhm 

    def process_HypMach_on_auxvar(self): 
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
            ad = HypMach.io_to_AffineDelta(x,y,d,cv,self.ma_dim,self.ma_order)
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

    """
    return: 
        (is input X null?,is input X partially null?),
        (is input Y null?,is input Y partially null?)
    """
    def io_stat(self):
        stat1 = type(self.x) != type(None)
        stat2 = type(self.y) != type(None)

        if stat1 and stat2: 
            assert len(self.x) == len(self.y)
        
        stat3,stat4 = True,True 
        if stat1:
            stat3 = None not in self.x 
        
        if stat2: 
            stat4 = None not in self.y 

        return (stat1,stat3),(stat2,stat4) 
    
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
    
    """
    given 2 points `(i0,o0)` and `(i1,o1)`, calculates 
    an <AffineDelta> that fits the two. 
    """
    @staticmethod
    def io_pointpair_to_AffineDelta(i0,i1,o0,o1,ma_order):

        assert vs_dim(o0) == vs_dim(o1) 

        M,A = None,None

        dx1 = i1 - i0 
        dx2 = o1 - o0 

        dx3 = safe_div(dx2,dx1)
        if not is_number(dx3,set()): 
            M = np.array(dx3)
        else: 
            M = dx3 

        # case: multiplication is first 
        if ma_order == 0: 
            t = i0 * M 
            A = o0 - t

            # ? 
            t2 = i1 * M 
            A2 = o1 - t2  
        else: 
            t = np.array(safe_div(o0,M)) 
            A = t - i0

            # ? 
            t = np.array(safe_div(o1,M)) 
            A2 = t - i1  
        return AffineDelta(M,A,ma_order)  

    """
    calculates a super-partition P, consisting of subsets of 
    indices, the indices corresponding to the (x_i,y_i) i/o 
    samples. Each subset of the super-partition is indices of 
    the elements (x_i,y_i) that are part of one <AffineDelta> 
    instance. 

    A super-partition contains all indices 0,1,...,(|x|-1), 
    and an index i_j can exist in more than one subset. 
    """
    def superpart_by_AffineDelta(self,ma_order,\
        index0_vec=None,index1_function=None): 

        l = len(self.x) 
        if type(index0_vec) == type(None): 
            index0_vec = [i for i in range(0,l)]  
            index1_function = lambda x: [i for i in range(x+1,l)] 

        hm_ = HypMem([],[],"MA")
        hm_.init_partition()
        for i in index0_vec: 
            x0,y0 = self.io_sample(i)
            rx = index1_function(i) 
            for j in rx: 
                x1,y1 = self.io_sample(j) 
                ad = IOFit.io_pointpair_to_AffineDelta(x0,x1,y0,y1,ma_order)
                hm_.add_to_partition((i,j),ad)
        return hm_ 

    """
    NOTE: method assumes all (x_i,y_i) pairs are of equal 
          dimension to one another. 
    """    
    def error_by_hyp(self,h0): 
        hm0 = HypMem(indices=[],info=[],mem_type="ERROR")
        l = len(self.x)
        
        for i in range(l): 
            x,y = self.io_sample(i)
            q0 = h0(x)
            dx = y - q0 
            hm0.add(i,dx) 
        return hm0 

    def cmp_two_hyp(self,h0,h1,cfunc1=default_cfunc2,cfunc2=default_cfunc2):
        hm0 = self.error_by_hyp(h0)
        hm1 = self.error_by_hyp(h1) 
        return hm0.cmp_error(hm1,cfunc1,cfunc2),[hm0,hm1]