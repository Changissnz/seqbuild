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