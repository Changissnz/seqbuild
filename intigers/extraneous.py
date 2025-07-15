import numpy as np 
from morebs2.matrix_methods import is_vector,is_number,is_valid_range 
from morebs2.measures import zero_div 
from morebs2.numerical_generator import modulo_in_range,\
    prg__LCG,euclidean_point_distance
from math import ceil

zero_div0 = lambda num,denum: zero_div(num,denum,0)

def to_trinary_relation(v1,v2):
    if v1 == v2: return 0 
    if v1 > v2: return 1 
    return -1 

def trinary_vector_to_indexvalue_map(tv):
    q = dict() 
    for (i,t) in enumerate(tv): 
        assert t in {-1,0,1} 
        if t != 0:
            q[i] = t 
    return q 

def trinary_diff(t0,t1,invertible_weight=0.5):
    assert t0 in {-1,1,0} and t1 in {-1,1,0}
    if set([t0,t1]) == {-1,1}: 
        return invertible_weight
    return int(abs(t0-t1))

def trinary_vector_invertible_difference(v1,v2,invertible_weight): 
    assert is_vector(v1) and is_vector(v2) 
    assert len(v1) == len(v2) 
    if is_vector(invertible_weight): 
        assert len(invertible_weight) == len(v1)

    def iw(i):
        if is_number(invertible_weight,set()): 
            return invertible_weight
        return invertible_weight[i] 

    vx = []
    for (i,(v1_,v2_)) in enumerate(zip(v1,v2)):
        iw_ = iw(i)
        vx.append(trinary_diff(v1_,v2_,iw_)) 
    return np.array(vx)  

"""
vector-input version of method<to_trinary_relation>; 

`zero_feature` results in non-absolute comparison for cases of (v1,v2) 
pairs that have 0.0 in them. 
`abs_feature` results in absolute comparison. 
To use these 2 features, set at most one of them to True. 
"""
def to_trinary_relation_v2(v1,v2,zero_feature:bool=False,abs_feature:bool=True):

    stat1 = is_vector(v1)
    stat2 = is_vector(v2) 

    def next_index(i):
        v1_ = v1 if not stat1 else v1[i]
        v2_ = v2 if not stat2 else v2[i]
        return v1_,v2_ 

    l = None 
    if stat1 and stat2: 
        l = len(v1) 
        assert l == len(v2)
    elif stat1:
        l = len(v1) 
    elif stat2:
        l = len(v2) 
    else:
        pass

    if type(l) == type(None):
        if zero_feature: 
            if v1 == 0.0 or v2 == 0.0:
                return to_trinary_relation(v1,v2)
        if abs_feature: 
            v1,v2 = np.abs(v1),np.abs(v2) 
        return to_trinary_relation(v1,v2) 
    
    i = 0 
    lx = [] 
    while i < l:
        x1,x2 = next_index(i)

        if zero_feature and (x1 == 0.0 or x2 == 0.0): 
            q = to_trinary_relation(x1,x2)
        elif abs_feature: 
            q = to_trinary_relation(np.abs(x1),np.abs(x2))
        else: 
            q = to_trinary_relation(x1,x2)

        lx.append(q) 
        i += 1
    return np.array(lx) 

def trinary_vector_intersection(v1,v2):
    assert is_vector(v1) and is_vector(v2) 
    assert len(v1) == len(v2) 
    return v1 * v2 

def active_trinary_vector_indices(v1,keys=(-1,1)): 
    return set([i for (i,v1_) in \
        enumerate(v1) if v1_ in keys])

def round_trinary(t,is_distance_roundtype:bool=True):
    assert is_number(t,set())

    if is_distance_roundtype: 
        X = [-1,0,1]
        i = np.argmin([abs(x - t) for x in X])
        return X[i] 
    if t < 0: return -1
    if t > 0: return 1 
    return 0 


def round_to_trinary_vector(V,is_distance_roundtype:bool=True):
    if is_number(V,set()): return round_trinary(V,is_distance_roundtype) 
    return np.array([round_trinary(v_,is_distance_roundtype) \
        for v_ in V]) 

#------------------------- safe division 

def safe_div(V1,V2):

    if type(V1) == list: 
        V1 = np.array(V1) 
    if type(V2) == list: 
        V2 = np.array(V2) 

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

    def value_at_index(i,is_x):
        V = V1 if is_x else V2

        if is_number(V,set()): 
            return V 
        return V[i]
    
    l = len(V1) if stat1 else len(V2) 
    q = [] 
    for i in range(l): 
        x = value_at_index(i,True)
        y = value_at_index(i,False) 
        q.append(zero_div0(x,y)) 
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
        centers = c 
        minnie = min(centers)
        maxie = max(centers)
        
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

    maxie = maxie + 1 if maxie == minnie else maxie 

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

def prg__single_to_range_outputter(lcg):
    def f(): 
        q1 = lcg() 
        q2 = lcg() 
        if q1 < q2:
            return (q1,q2)
        q2 = (q1 + 1) + q1 - q2 
        return (q1,q2) 
    return f 

def prg__single_to_ndim_index_outputter(lcg,dim):

    def f():
        x = np.zeros((len(dim),),dtype=np.int32) 
        for i in range(len(dim)): 
            x[i] = int(lcg()) % dim[i] 
        return x 
    return f 

def prg__single_to_nvec(prg,n):

    def f():
        q = np.zeros((n,))
        for i in range(n): 
            q[i] = prg()
        return q 
    return f 


#------------------------------- operations related to euclid's distance

"""
d := int, positive, dimension of points
r := iterable, length 2, ordered range for values of interest
inclusive := bool, 
"""
def euclid_pd_over_range__int(d,r,inclusive=False,\
    output_type=0):
    assert type(d) == int and d > 0 
    assert r[0] < r[1] 
    assert type(r[0]) == type(r[1])
    assert type(r[0]) == int 
    assert len(r) == 2 
    assert output_type in {0,1} 

    rx0,rx1 = r[0],r[1] 

    if inclusive: 
        rx1 += 1
    
    p = np.zeros((d,),dtype=np.float32)
    dx = []
    for r0 in range(rx0,rx1):
        p2 = np.ones((d,),dtype=np.float32)
        p2 = p2 * r0
        pd = euclidean_point_distance(p,p2)

        x = pd if output_type == 0 else (r0,pd)
        dx.append(x)

    return dx 

def multiple_sqrt_seq(x,mrange,output_type=0):
    assert is_valid_range(mrange,True,False) 
    assert output_type in {0,1}

    qx = []
    for i in range(mrange[0],mrange[1]): 
        dx = np.sqrt(x * i)
        y = (i,dx) if output_type else dx 
        qx.append(y) 
    return qx 