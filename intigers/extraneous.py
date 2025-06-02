import numpy as np 
from morebs2.matrix_methods import is_vector
from morebs2.measures import zero_div 
from morebs2.numerical_generator import modulo_in_range

zero_div0 = lambda num,denum: zero_div(num,denum,0)

def safe_div(V1,V2):
    stat1 = is_vector(V1)
    stat2 = is_vector(V2) 

    if not stat1 and not stat2: 
        return zero_div0(V1,V2)

    if stat1 and stat2:
        assert len(V1) == len(V2)
        q = []
        for x,x2 in zip(V1,V2):
            q2 = zero_div0(x,x2)
            q.append(q2)
        return np.array(q)

    i = 0
    VX = None
    f = None 

    if stat1:
        VX = V1
        f = lambda x: zero_div0(x,V2)
    else: 
        VX = V2
        f = lambda x: zero_div0(V1,x) 

    q = [] 
    for x in VX:
        q.append(f(x))
    return q 
    
def safe_npint32_value(v):
    r = (np.iinfo(np.int32).min,np.iinfo(np.int32).max)
    if v >= r[0] and v <= r[1]: return np.int32(v) 
    v_ = modulo_in_range(v,r) 
    return np.int32(v_)

def safe_npint32_vec(V):
    return np.array([safe_npint32_value(v_) for \
        v_ in V],dtype=np.int32) 

def safe_npint32__prg_vec(prg,sz):
    v = np.zeros((sz,),dtype=np.int32) 
    for i in range(sz):
        v[i] = safe_npint32_value(prg()) 
    return v 

def modulated_vec_op(v1,v2,op):

    V,V2 = None,None 
    if len(v1) > len(v2): 
        V,V2 = v1,v2
    else:
        V,V2 = v2,v1 

    q = []
    for (i,v) in enumerate(V): 
        i2 = i % len(V2)
        v_ = V2[i2] 
        q.append(op(v,v_))
    return np.array(q) 

def modulated_vecdot(v1,v2,op1,op2):
    V = modulated_vec_op(v1,v2,op1)
    return op2(V)