from morebs2.matrix_methods import vector_to_string,\
    string_to_vector,euclidean_point_distance,is_vector,is_number
from intigers.extraneous import trinary_diff,to_trinary_relation_v2,safe_div,\
    trinary_vector_invertible_difference,active_trinary_vector_indices,\
    trinary_vector_intersection,round_to_trinary_vector,round_trinary
from collections import defaultdict
from .n2m_index import * 

#---------------------------------- auxiliary methods with trinary 
#---------------------------------- vectors and n2m index-mapping 

euclidean_point_distance__zero = lambda p : euclidean_point_distance(p,np.zeros((len(p),)))

def invertible_trinary_euclidean_distance(v1,v2,invertible_weight=0.5): 
    tvdiff = trinary_vector_invertible_difference(v1,v2,invertible_weight)
    return euclidean_point_distance__zero(tvdiff)

def n2m_delta_correlate__orderedproportional(delta_x_relation,delta_err): 
    assert is_vector(delta_x_relation)
    assert is_vector(delta_err)

    part = paired_n2m_partition__orderedproportional(\
        len(delta_x_relation),len(delta_err)) 
    q = []
    for p in part: 
        x0,x1 = min(p[0]),max(p[0]) + 1
        qr = np.mean(delta_x_relation[x0:x1]) 
        t = round_trinary(qr) 
        q.append(t) 
    return np.array(q) * delta_err  

#----------------------------------------------------------------------------


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
        vx = x_[0] 
        vx2 = vx * x_[1] 
        vx_.append(vx2) 
        s += x_[1] 
    vx_ = np.array(vx_) 
    vx_ = np.sum(vx_,axis=0) 
    vx_ = safe_div(vx_,s)
    q = np.round(vx_,1)
    return np.array(q,dtype=int)

"""
function calculates an integer i that is to be a divisor 
(or a denumerator with 1 as numerator). Integer i is 
0 if `x0` and `x_ref` do not intersect in active indices, 
or the number `k` of active indices in `x0` not active in 
`x_ref` [divisor of `k + 1`]. 

The intent of this function is to penalize the fit score 
of `x0` with `x_ref` if `x0` has extra active indices in 
comparison with `x_ref`. 
"""
def N2M_AC__setdiff_weighted_support_(x0,x_ref): 

    tv = trinary_vector_intersection(x0,x_ref)

    if np.all(tv == 0): 
        return 0 

    ai0 = active_trinary_vector_indices(x0)
    ai1 = active_trinary_vector_indices(x_ref) 
    sz = len(ai0 - ai1) 
    return sz + 1 

def N2M_AC__setdiff_weighted_support_function(coeff=2.0): 
    assert coeff > 0.0

    def f(x0_,x_ref_):
        q = N2M_AC__setdiff_weighted_support_(x0_,x_ref_) 
        if q == 1: return q 
        return q * coeff 

    return f  

class N2MAutocorrelator:

    def __init__(self,nm,is_active_seqc:bool=True):
        assert type(is_active_seqc) == bool 

        self.nm = None 
        self.set_nm(nm) 
        # sign-change vector of input -> 
        # sign-change vector of output -> 
        # frequency of occurrence given 
        #     sample inputs 
        self.ftable = defaultdict(None)
        self.seqc = defaultdict(list)
        self.seq_stat = is_active_seqc
        return 
    
    def set_nm(self,nm): 
        assert_nm(nm) 
        self.nm = nm 

    """
    adds a pair of (x_i,e_i) to memory; e_i the error term for input x_i. 

    Outputs a boolean for the status of normal relation b/t the 
    input variables. This status is calculated through mean-based 
    functions.     
    """
    def add(self,x0,x1,e0,e1):
        assert len(x0) == len(x1) 
        assert len(x0) == self.nm[0] 
        assert len(e0) == len(e1) 
        assert len(e0) == self.nm[1] 

        # get the trinary vector for 
        r1 = to_trinary_relation_v2(e1,e0,True,False)  
        self.add_v2(x0,x1,r1)

        return
    
    def add_v2(self,x0,x1,dx_err): 
        r0 = to_trinary_relation_v2(x1,x0,True,False)
        r1 = dx_err 
        if is_number(r1):
            r1 = [r1] 

        s0 = vector_to_string(r0,int)
        s1 = vector_to_string(r1,int)

        if s0 not in self.ftable: 
            self.ftable[s0] = defaultdict(int) 
        self.ftable[s0][s1] += 1

        if self.seq_stat: 
            self.seqc[s0].append(s1)

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

        r0 = to_trinary_relation_v2(x1,x0,True,False)
        ks = self.closest_keyset(r0)
        qs = [] 
        for k in ks: 
            k_ = string_to_vector(k,int)
            sk = self.summarize_key_relation(r0,k_,frequency_type)
            o = cfunc1(sk) 
            qs.append(o) 
        return cfunc2(qs)
    
    def induce_derivative_v2(self,x0,x1):
        ## NOTE: v v 
        r0 = to_trinary_relation_v2(x1,x0,True,False)
        q = np.zeros((self.nm[1],))

        c = 0 
        m = 0.0 
        for k in self.ftable.keys():
            v = self.sample_support_for_dvec(k,r0)
            if is_number(v): v = [v] 
            m = max(np.abs(v + [m])) 
            q += v 
            c += 1 
        q = safe_div(q,m)
        return round_to_trinary_vector(q) 
        
    def sample_support_for_dvec(self,skey,d,coefficient = 1.0): 
        vx = string_to_vector(skey,int) 
        f = N2M_AC__setdiff_weighted_support_function(coeff=coefficient) 
        q = f(vx,d) 
        if q == 0: 
            return 0  

        qs = self.summarize_key_relation(d,vx,output_type="absolute")
        if len(qs) == 0: 
            assert False 
        qs2 = [qs_[0] * qs_[1] for qs_ in qs]
        qs2 = np.array(qs2) 
        qs2 = np.sum(qs2,axis=0) 
        return safe_div(qs2,q)

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
    
    # TODO: test this. 
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
            qv_ = string_to_vector(sk_[0]) 
            # NOTE: v v 
            qv = n2m_delta_correlate__orderedproportional(np.array(q),qv_) 
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
        q = [(q_[0],safe_div(q_[1],s)) for q_ in q] 
        return q 

    def cycle_patterns(self): 
        return -1 
