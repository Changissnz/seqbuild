from morebs2.matrix_methods import is_number,vector_to_string,\
    string_to_vector,euclidean_point_distance
from intigers.extraneous import to_trinary_relation,zero_div0,safe_div,\
    trinary_vector_invertible_difference
from collections import defaultdict
from copy import deepcopy
import numpy as np 

euclidean_point_distance__zero = lambda p : euclidean_point_distance(p,np.zeros((len(p),)))

def invertible_trinary_euclidean_distance(v1,v2,invertible_weight=0.5): 
    tvdiff = trinary_vector_invertible_difference(v1,v2,invertible_weight)
    return euclidean_point_distance__zero(tvdiff)

"""
NOTE: no argument check. 
"""
def N2M_AC__most_frequent_cfunc(x,assume_sorted:bool=True):
    assert len(x) > 0 

    if not assume_sorted: 
        x = sorted(x,lambda y: y[1]) 
    return x[-1]

def N2M_AC__weighted_average_cfunc(x):

    assert len(x) > 0 

    # iterate through and collect weighted mean 
    vx_ = []
    s = 0 
    for x_ in x:
        vx = string_to_vector(x_[0],castFunc=int) 
        vx2 = vx * x_[1] 
        vx_.append(vx2) 
        s += x_[1] 
    vx_ = np.array(vx_) 
    vx_ = np.sum(vx_,axis=0) 
    vx_ = safe_div(vx_,s)
    q = np.round(vx_,1)
    return np.array(q,dtype=int)

class N2MAutocorrelator:

    def __init__(self,nm):
        self.nm = None 
        self.set_nm(nm) 
        # sign-change vector of input -> 
        # sign-change vector of output -> 
        # frequency of occurrence given 
        #     sample inputs 
        self.ftable = defaultdict(None)
        self.seqc = defaultdict(list)
        return 
    
    def set_nm(self,nm): 
        assert type(nm) == tuple and len(nm) == 2 
        assert type(nm[0]) in {int,np.int32,np.int64}
        assert type(nm[1]) in {int,np.int32,np.int64}
        assert min(nm) >= 0 
        self.nm = nm 

    """
    adds a pair of (x_i,e_i) to memory; e_i the error term for input x_i. 

    Outputs a boolean for the status of normal relation b/t the 
    input variables. This status is calculated through mean-based 
    functions.     
    """
    def add(self,x0,x1,e0,e1):

        # get the trinary vector for 
        r0 = to_trinary_relation(x1,x0)
        r1 = to_trinary_relation(e1,e0)  

        s0 = vector_to_string(r0,int)
        s1 = vector_to_string(r1,int)

        if s0 not in self.ftable: 
            self.ftable[s0] = defaultdict(int) 
        self.ftable[s0][s1] += 1
        self.seqc[s0].append(s1)
        return
    
    """
    guesses the difference 
        `e1 -e0`, 
    `e1` error-term for `x1` and `e0` 
    for `x0`. Guess is based on the 
    """
    def induce_derivative(self,x0,x1,\
        cfunc1=N2M_AC__most_frequent_cfunc,\
        cfunc2=N2M_AC__weighted_average_cfunc,\
        frequency_type:str="absolute"):
        assert frequency_type in {"ratio","absolute"}

        r0 = to_trinary_relation(x1,x0)
        ks = self.closest_keyset(r0)
        
        qs = [] 
        for k in ks: 
            sk = self.summarize_key_relation(r0,k,frequency_type)
            #sk = self.summarize_key(k,frequency_type)
            o = cfunc1(sk) 
            qs.append(o) 
        return cfunc2(qs)
    
    """
    Calculates the set of closest elements to 

    d := stringized vector, difference vector between two inputs x1,x0 (x1-x0). 
    dfunc := distance function F(point1,point2); usually set to 
                `invertible_trinary_euclidean_distance` or 
                `euclidean_point_distance`. 
    """
    def closest_keyset(self,d,dfunc=invertible_trinary_euclidean_distance): 
        qx = list(self.ftable.keys()) 

        # case: empty memory 
        if len(qx) == 0: return None 

        qx = [string_to_vector(qx_,castFunc=int) for qx_ in qx] 
        qx = [(qx_,dfunc(qx_,d)) for qx_ in qx]
        qx2 = sorted(qx,key=lambda x:x[1]) 

        rx = qx2.pop(0) 
        lx = [rx[0]] 

        stat = True 
        while stat: 
            if len(qx2) == 0: 
                stat = False 
                continue 
            
            rx2 = qx2.pop(0)
            if np.round(rx[1] - rx2[1],5) == 0.0: 
                lx.append(rx2[0]) 
            else: 
                stat = False 
         
        lx = [vector_to_string(lx_,castFunc=int) for lx_ in lx] 
        return lx 
    
    """
    summarizes the delta errors for an input-difference vector `d`. 
    The reference input-difference vector `rd` is used to calculate 
    a trinary vector that contains information on the direct or inverse 
    correlation between the reference and `d` on the matter of `d`'s 
    information, index correlations to the m-space (output space). 
    """
    def summarize_key_relation(self,rd,d,output_type:str="ratio"):
        q = safe_div(rd,d) 

        sk = self.summarize_key(d,output_type) 
        sk2 = []
        for sk_ in sk: 
            qv = string_to_vector(sk_[0]) 
            qv = qv * q 
            sk2.append((qv,sk_[1]))
        return sk2 
    
    """
    summarizes the m-space index correlations for 
    input-difference vector `d`. If `output_type` 
    is `absolute`, then output a sequence with 
    elements of the form 
        (error delta, frequency). 
    If mode `ratio` is used instead, then the 
    frequencies are normalized for a sum of 1.0. 
    """
    def summarize_key(self,d,output_type:str="ratio"):
        assert output_type in {"ratio","absolute"}

        if type(d) != str: 
            d = vector_to_string(d,int) 
        q = self.ftable[d]
        q = [(k,v) for k,v in q.items()]
        q = sorted(q,key=lambda x:x[1])

        if output_type == "absolute": 
            return q 
    
        s = sum([q_[1] for q_ in q]) 
        q = [(q_[0],zero_div0(q_[1],s)) for q_ in q] 
        return q 

    def cycle_patterns(self): 
        return -1 

"""
def f(vec,index):
    return 
"""
class N2MVectorFunction:

    def __init__(self):
        return 
    
