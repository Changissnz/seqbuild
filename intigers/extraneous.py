import numpy as np 
from morebs2.matrix_methods import is_vector
from morebs2.measures import zero_div 
from morebs2.numerical_generator import modulo_in_range
from math import ceil

zero_div0 = lambda num,denum: zero_div(num,denum,0)

#------------------------- safe division 

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

#----------------------- vector-to-vector operation 

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

#------------------------ operations related to LCM and GCD 

def lcm_times(x0,x1,m): 

    def f(): 
        return np.lcm(x0,x1) * m 
    return f 

"""
Generates a sequence of sequences (called `sets` in function name).

The cumulative number of integers from the sequences is aimed for `n`.
There are `c` centers. Each center is the base for variable multiple 
values provided by pseudo-random number generator `prg`. The ratio `r` 
specifies the number of elements in (n-c) that are not a multiple of any
center, except for {-2,-1,0,1,2} since those numbers are relatively 
common multiples. The `prg` draws values for the center from the range 
`crange` and values for the multiples from `mrange`. For the elements 
that are not a multiple of any center (except for {-2,-1,0,1,2}), the 
integer `num_attempts_per_nc` specifies the number of attempts the `prg` 
makes to output a satisfying integer value. If the number of attempts 
fails for any of these integers, algorithm terminates.

NOTE: there is a non-null probability the algorithm does not produce the 
      wanted `n` values. Use with discretion and awareness of the parameters.
"""
def prg__integer_sets__mult(n,c,r,crange,mrange,prg,\
    num_attempts_per_nc):
    assert type(n) == int and type(r) == float 
    assert type(c) in {int,list} 

    lc = c if type(c) == int else len(c) 
    assert n >= lc and lc > 0
    assert r >= 0.0 and r <= 1.0 

    if type(c) == int:
        assert len(crange) == 2
        assert crange[0] < crange[1] 
        
    assert mrange[0] < mrange[1]

    def not_divisible_by_centers(ix,C):

        for c in C:
            if abs(c) in {0,1,2}: 
                continue
            if ix / c == ix // c: return False
        return True 

    def attempt_p_times(C,p,rnge):
        while p > 0: 
            ix = modulo_in_range(prg(),rnge)
            if not_divisible_by_centers(ix,C):
                return ix
            p -= 1 
        return None

    qx = []
    centers = []

    minnie = float('inf')
    maxie = -float('inf')

    # declare the centers
    if type(c) == int: 
        for _ in range(c):
            c_ = modulo_in_range(prg(),crange)
            qx.append([c_])
            centers.append(c_)
            minnie = min(minnie,c_)
            maxie = max(maxie,c_) 
    else: 
        qx = [[c_] for c_ in c]
        
    # get ratio of remaining elements that will not be a multiple of
    # any center 
    l = int(ceil((n - lc) * r))

    # assign values from `prg` to the centers
    rl = n - lc - l
    for _ in range(rl):
        ci = prg() % len(centers)
        ctr = centers[ci]
        m = modulo_in_range(prg(),mrange)

        q = ctr * m 
        qx[ci].append(q) 

        minnie = min(minnie,q)
        maxie = max(maxie,q) 

    # assign values that are not divisible by any centers
    # except for 1 and 2
    rx = []
    for _ in range(l):
        ix = attempt_p_times(centers,num_attempts_per_nc,[minnie,maxie])
        if type(ix) == type(None): 
            break 
        rx.append(ix)

    qx.append(rx)
    return qx 

"""
an overhead version of the method `prg__integer_sets__mult`. 
Outputs two values: 
- a flattened version (one-dimensional) of the sequence 
  of sequences,
- ?number of integer output values equals `n`?
"""
def prg__integer_seq__mult(n,c,r,crange,mrange,prg,\
    num_attempts_per_nc):

    qx = prg__integer_sets__mult(n,c,r,crange,mrange,prg,\
        num_attempts_per_nc)

    qx2 = []
    for qx_ in qx: qx2.extend(qx_)

    return qx2,len(qx2) == n