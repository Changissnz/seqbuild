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
        ##def expected_diff(self,target_value,x,op_type,dfunc=sub):