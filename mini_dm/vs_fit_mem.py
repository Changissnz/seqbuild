from .vs_fit import * 

"""
a hypothesis container for an <MADescriptor> instance. 
"""
class MADHyp(MADescriptor):

    def __init__(self,rv_vec,rvt_vec,t_vec,s_vec,d_vec,ma_order=None):
        super().__init__(rv_vec,rvt_vec,t_vec,s_vec,d_vec,ma_order)
        return 
    
    @staticmethod
    def naive_hyp_for_MADescriptor(md):
        md2 = md.default_set_naive("rv_vec")
        md2 = md2.default_set_naive("rvt_vec")

        return MADHyp(md2.rv_vec,md2.rvt_vec,\
            md2.t_vec,md2.s_vec,md2.d_vec,\
            md2.ma_order)


class MAHypMem: 

    def __init__(self,indices,info,mem_type):
        assert mem_type in {"MA","VECQUAL"}
        self.indices = indices 
        self.info = info 
        self.mem_type = mem_type  
        for i in self.info: assert self.type_check(i)

    def type_check(self,element): 
        if self.mem_type == "MA": 
            return type(element) == AffineDelta
        else: 
            return type(element) == MADHyp
        
    def add(self,idn,info): 
        assert self.type_check(info) 
        self.indices.append(idn) 
        self.info.append(info)
        return

class MAHypMach:

    def __init__(self,x,y,d=None):
        #assert x and is_vector(y) 
        assert None not in x 
        assert None not in y
        assert len(x) == len(y) 
        self.x = x 
        self.y = y 
        self.d = None 
        self.load_distance_map(d) 
        return
    
    def load_distance_map(self,d): 
        if type(d) == type(None): return 

        assert len(d) == len(self.x) 
        self.d = d 

