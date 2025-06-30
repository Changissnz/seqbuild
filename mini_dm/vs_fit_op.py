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
    
    def 