from .vs_fit import * 
from types import MethodType,FunctionType

class MADHyp(MADescriptor):

    def __init__(self,rv_vec,rvt_vec,t_vec,s_vec,d_vec,ma_order=None):
        super().__init__(rv_vec,rvt_vec,t_vec,s_vec,d_vec,ma_order)
        return 
    


class VSTransform:

    def __init__(self,ad:AffineDelta):
        assert type(ad) == AffineDelta
        self.ad = ad
        return 
    
    #------------------------ I/O comparison functions
    
    def cmp_ad(self,ad1):
        sz_diff = abs(self.ad.size() - ad1.size())
        diff_ad = (self.ad - ad1).size() 
        return sz_diff,diff_ad 
    
    def diff_ad(self,x,ad1,op_type="all",\
        dfunc=lambda x,x2: np.sum(np.abs(x - x2))): 

        tv = ad1.fit(x)
        return self.ad.expected_diff(tv,x,op_type,\
            dfunc)

    def diff_derivative(self,x,m_delta,a_delta,dfunc=lambda x,x2:x-x2):
        ad1 = AffineDelta(self.m + m_delta,self.a + a_delta,\
            self.ad.ma_order)
        return dfunc(ad1.fit(x),self.ad.fit(x))

    def cvec_diff_ad(self,x1,x2,dfunc=lambda x,x2:np.abs(x - x2)):
        return dfunc(self.ad.cvec(x1),self.ad.cvec(x2))

    def to_hypdiff_func(self): 

        def f(x,hyp_y):
            actual =  self.ad.fit(x) 
            return actual - hyp_y 
        return f 


class IOFit:

    def __init__(self,x,y,unknown_func,hypdiff_func):
        assert type(unknown_func) in {type(None),FunctionType,MethodType}

        self.x = x
        self.y = y 
        self.unknownf = unknown_func
        self.hyp = None 
        self.hypdiff_func = hypdiff_func

    def io_stat(self):
        stat1 = type(self.x) != type(None)
        stat2 = type(self.y) != type(None)

        if stat1 and stat2: 
            assert len(self.x) == len(self.y)
        return stat1,stat2 
    
    def load_hyp(self,hyp): 
        assert type(hyp) == MADHyp
        self.hyp = hyp 

    def hyp_proc(self):
        assert type(self.hyp) == MADHyp
        return -1 

    def io_sample_proc(self,i,default_source = "y"):
        assert default_source in {"y","unknown"} 
        q = self.x[i]

        q2 = None 
        if default_source == "y":
            if type(self.y) != type(None): 
                q2 = self.y[i]
        

        if type(q2) == type(None):
            q2 = self.unknownf(q)
        return (q,q2)