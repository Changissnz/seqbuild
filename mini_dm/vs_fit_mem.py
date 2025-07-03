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
        q = np.abs(y - x)
        if not is_number(q,set()): 
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
        for x_ in ma:
            if is_number(x_,set()):
                sz_ma.append(0)
            else:
                sz_ma.append(len(x_))

        if is_number(y,set()): assert max(sz_ma) == 0
        else: 
            assert max(sz_ma) == len(y) 

        def place_element(ma_index,e_index,element):
            if is_number(ma[ma_index],set()):
                ma[ma_index] = element 
                return 
            ma[ma_index][e_index] = element
            return   
        
        def get_x(i): 
            if is_number(x,set()): return x 
            return x[i] 
        
        def get_y(i):
            if is_number(y,set()): return y 
            return y[i] 
        
        def get_d(i): 
            if is_number(d,set()): return d 
            return d[i]
        
        def get_cv(i1,i2):
            if is_number(cv[i1],set()): 
                return cv[i1] 
            return cv[i1][i2] 
        
        # case: M,A are singletons 
        if is_number(q,set()): 
            q = [q] 

        # iterate through and solve for every index 
        for (i,q_) in enumerate(q): 

            d_ = get_d(i)
            cv0_ = get_cv(0,i)
            cv1_ = get_cv(1,i)

            x_ = get_x(i) 
            y_ = get_y(i)
            if q_ < d_: 
                c1,c2 = cv0_ * d_,cv1_ * d_ 

                cx = [c1,c2]
                if y_ > x_: 
                    j = np.argmin(cx) 
                else: 
                    j = np.argmax(cx) 
                cx[j] *= -1 
                c1,c2 = cx[0],cx[1] 

                if ma_order == 0: 
                    mx = safe_div(x_ + c1,x_)
                    mx2 = c2 
                else: 
                    mx = c1 
                    mx2 = safe_div(y_,x_ + c1)  

                mxs = [mx,mx2]
                place_element(0,i,mxs[ma_order])
                place_element(1,i,mxs[(ma_order + 1) % 2])
            elif q_ == d_: 

                if y_ < x_: 
                    d_ = -d_ 
                q = cv0_ * d_

                if ma_order == 0:
                    mx = safe_div(x_ + q,x_)
                else: 
                    mx = q
                place_element(ma_order,i,mx) 

                q = cv1_ * d_
                if ma_order == 0: 
                    mx2 = q 
                else:
                    mx2 = safe_div(y_,x_ + (cv0_ * d_)) 
                place_element((ma_order + 1) % 2,i,mx2)
            else: 
                raise ValueError("d[{}] has to be at least the absolute diff. of y and x.".format(i))

        return AffineDelta(ma[0],ma[1],ma_order) 

    def naive_MA_for_sample(x,y,ma_dim,ma_order):
        return 
    
    def naive_vecqual_for_sample(self,index):
        return -1 