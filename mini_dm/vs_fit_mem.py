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

    def __init__(self,indices=[],info=[],mem_type="MA"):
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

    def __init__(self,x,y,d=None,mem_type="MA"):
        #assert x and is_vector(y) 
        assert None not in x 
        assert None not in y
        assert len(x) == len(y) 
        self.x = x 
        self.y = y 
        self.d = None 
        self.mem_type = mem_type
        self.load_distance_vec(d) 
        self.load_mem() 
        return
    
    def load_distance_vec(self,d): 
        if type(d) == type(None): return 
        assert len(d) == len(self.x) 
        assert is_vector(d) 
        self.d = d 

    def load_mem(self,indices=[],info=[]): 
        self.mhm = MAHypMem(indices,info,mem_type=self.mem_type)

    #------------------- hypothesis generation 

    @staticmethod 
    def io_to_AffineDelta(x,y,d,cv,ma_dim,ma_order): 
        if ma_order == 1: 
            cv = cv[::-1] 


        q = np.abs(y - x)

        # case: M,A are singletons 
        if is_number(q,set()): 
            return -1 

        assert len(d) == len(q) 

        # set the initial MA pair 
        ma = [None,None]

        if ma_dim[0] == 0: 
            ma[0] = 0.0
        else: 
            ma[0] = np.zeros((len(q),),dtype=float)

        if ma_dim[1] == 0: 
            ma[1] = 0.0
        else: 
            ma[1] = np.zeros((len(q),),dtype=float)

        sz_ma = []
        for x in ma:
            if is_number(x,set()):
                sz_ma.append(0)
            else:
                sz_ma.append(len(x))

        if is_number(y): assert max(sz_ma) == 0
        else: 
            assert max(sz_ma) == len(y) 

        def place_element(ma_index,e_index,element):
            if is_number(ma[ma_index]):
                ma[ma_index] = element 
                return 
            ma[ma_index][e_index] = element
            return   

        for (i,q_) in enumerate(q): 
                if q_ < d[i]: 
                    other = d[i] - q_ 

                    if cv[0] < cv[1]: 
                        place_element(0,i,min([q_,other]))
                        place_element(1,i,max([q_,other]))
                    else: 
                        place_element(0,i,max([q_,other]))
                        place_element(0,i,min([q_,other]))
                elif q_ == d[i]: 
                    place_element(0,i,cv[0][i] * d[i])
                    place_element(1,i,cv[1][i] * d[i])
                else: 
                    raise ValueError("d[{}] has to be at least the absolute diff. of y and x.".format(i))

        return AffineDelta(ma[0],ma[1],ma_order) 

    def naive_MA_for_sample(x,y,ma_dim,ma_order):
        return 
    
    def naive_vecqual_for_sample(self,index):
        return -1 