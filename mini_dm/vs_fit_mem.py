from .vs_fit import * 
from types import MethodType,FunctionType
from morebs2.matrix_methods import is_2dmatrix

def get_vs_element(L,i,cf=lambda x: is_number(x,set())): 
    if cf(L): return L 
    return L[i] 

def vs_dim(L): 
    if is_number(L,set()): return 0  
    return len(L)

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

"""
abstract representation of a general hypothesis, which 
consists of 1 function `h` and the capability to update 
`h` given 1 input value. 
"""
class GHyp: 

    def __init__(self,h,dx_h,vf,error_term=None):
        assert type(h) in {MethodType,FunctionType}
        assert type(dx_h) in {MethodType,FunctionType}
        assert type(vf) in {MethodType,FunctionType}
        assert is_number(error_term) or is_vector(error_term) \
            or type(error_term) == type(None)
        
        self.h = h 
        self.dx_h = dx_h 
        self.vf = vf 
        self.num_updates = 0 
        self.error_term = error_term

    def update(self,q,make_copy:bool=True): 
        g = self if not make_copy else deepcopy(self) 
        g.h,g.vf = g.dx_h.update(q) 
        g.num_updates += 1 
        return g 
    
    def load_error(self,et):
        assert is_number(et) or is_vector(et) 
        self.error_term = et 

    def vector_form(self): 
        return self.vf() 

def default_cfunc1(S): 
    S_ = None 
    try:
        S_ = np.array(S) 
    except:
        assert False, "sequence of elements is irregular in shape."

    if len(S_.shape) == 1:
        return np.mean(S_) 
    return np.mean(S_,axis=0)

def default_cfunc2(S): 
    S = np.abs(S) 
    return default_cfunc1(S) 

"""
memory structure to contain values for qualities of 
hypotheses on the mapping function between the input and 
output space. 
"""
class HypMem: 

    def __init__(self,indices=[],info=[],mem_type="MA"):
        assert mem_type in {"MA","VECQUAL","ERROR","GHYP"}
        self.indices = indices 
        self.info = info 

        self.mem_type = mem_type  
        
        for i in self.info: assert self.type_check(i)

        # vars. used for `mem_type` := ERROR 
        self.condensed_error1 = None 
        self.condensed_error2 = None 

        # vars used for partitioning of function set 
        self.partition = None 

    def type_check(self,element): 
        if self.mem_type == "MA": 
            return type(element) == AffineDelta
        elif self.mem_type == "VECQUAL": 
            return type(element) == MADHyp
        elif self.mem_type == "GHYP": 
            return type(element) == GHyp
        else: 
            return is_number(element,set()) or is_vector(element)
        
    def add(self,idn,info): 
        assert self.type_check(info) 
        if self.mem_type != "GHYP":
            self.indices.append(idn) 
            self.info.append(info)
        else:
            i = len(self.info) 
            for (j,i2) in enumerate(self.info): 
                if i2.error > info.error:
                    i = j
                    break 
            self.info.insert(i,info) 

        return

    def init_partition(self):
        self.partition = []  

    def add_to_partition(self,idn,info): 
        ix = None 
        for i,inf in enumerate(self.info): 
            if inf == info:
                ix = i 
                break 
        
        if type(ix) != type(None):
            self.partition[ix] |= {idn[0],idn[1]}
        else: 
            self.partition.append({idn[0],idn[1]})
            self.info.append(info) 

    def condense_error_term(self,cfunc1=default_cfunc2,cfunc2=default_cfunc2): 
        assert self.mem_type == "ERROR" 
        self.condensed_error1 = cfunc1(self.info) 
        if is_vector(self.condensed_error1):
            self.condensed_error2 = cfunc2(self.condensed_error1) 
    
    def c_error(self,i): 
        assert i in {1,2} 
        x = self.condensed_error1 if i == 1 else self.condensed_error2

        if type(x) == type(None): 
            self.condense_error_term() 
        return self.condensed_error1 if i == 1 else self.condensed_error2

    """
    compares error values with those of another <HypMem> `hm`. The functions 
    `cfunc1,cfunc2` are used to possibly condense the multi-dimensional form 
    of the initial error. 
    """
    def cmp_error(self,hm,cfunc1=default_cfunc2,cfunc2=default_cfunc2):
        assert type(hm) == HypMem
        assert hm.mem_type == "ERROR" 
        assert hm.mem_type == self.mem_type 

        self.condense_error_term(cfunc1,cfunc2) 
        hm.condense_error_term(cfunc1,cfunc2) 

        # get trinary relation 
        c0,c1 = None,None 

        c0 = to_trinary_relation_v2(self.condensed_error1,\
            hm.condensed_error1,zero_feature=False,abs_feature=True)
        
        if type(self.condensed_error2) != type(None) and \
            type(hm.condensed_error2) != type(None): 
            c1 = to_trinary_relation_v2(self.condensed_error2,\
                hm.condensed_error2,zero_feature=False,abs_feature=True) 
        return c0,c1 

    def clear(self): 
        self.indices.clear() 
        self.info.clear() 

class HypMach:

    def __init__(self,x,y,d=None,mem_type="MA"):
        assert not np.any(x == None)
        assert not np.any(y == None) 
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
        assert is_vector(d) or is_2dmatrix(d) or type(d) == list 
        self.d = d 

    def load_mem(self,indices=[],info=[]): 
        self.mhm = HypMem(indices,info,mem_type=self.mem_type)

    def load_ma_hyp_dict(self,mhd,clear_mem:bool=True):
        if clear_mem: self.mhm.clear() 

        for k,v in mhd.items(): 
            self.mhm.add(k,v) 

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

        """
        adds an element according to the index arguments. 
        """
        def place_element(ma_index,e_index,element):
            if is_number(ma[ma_index],set()):
                ma[ma_index] = element 
                return 
            ma[ma_index][e_index] = element
            return   
        
        def get_cv(i1,i2):
            if is_number(cv[i1],set()): 
                return cv[i1] 
            return cv[i1][i2] 
        
        # case: M,A are singletons 
        if is_number(q,set()): 
            q = [q] 

        # iterate through and solve for every index 
        for (i,q_) in enumerate(q): 
            d_ = get_vs_element(d,i) 
            cv0_ = get_cv(0,i)            
            cv1_ = get_cv(1,i)

            x_ = get_vs_element(x,i)
            y_ = get_vs_element(y,i)
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
    
    def load_delta(self,d): 
        return -1 