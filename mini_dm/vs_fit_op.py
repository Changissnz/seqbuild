from .vs_fit import * 

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
