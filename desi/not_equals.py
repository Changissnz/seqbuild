from intigers.extraneous import * 
from morebs2.fit_2n2 import * 
from mini_dm.minmax_freq import  vec_to_frequency_map
from types import MethodType,FunctionType

DEFAULT_ADDITIVE = 10 ** -6 
DEFAULT_SUBMAT_TYPES = ["L+U","L+L","R+U","R+L"]

def DEFAULT_NOTEQUALS_ADDITIVE(d,m): 
    assert round(m - 0.0,5) > 0.0

    V = None 
    if type(d) in {int,np.int32,np.int64}: 
        V = np.arange(1,d+1)
        V = V * m 
        return V 

    assert type(d) == tuple and len(d) == 2 
    V = np.arange(1,(d[0] * d[1]) + 1)
    V = V * m 
    V = np.reshape(V,d)
    return V 

"""
compares two vectors V1,V2 of equal lengths. 
"""
def not_equals__pairvec(V1,V2,prg,indices=None):
    assert is_vector(V1) 
    assert is_vector(V2) 
    assert len(V1) == len(V2) 

    V1 = np.array(V1)
    V2 = np.array(V2)

    qr = np.ones((len(V2),),dtype=np.float64)
    vx = np.ones((len(V2),),dtype=np.float64)
    vx = vx * DEFAULT_ADDITIVE

    ix = set([i for i in range(len(V2))]) 

    if type(indices) == type(None):
        indices = sorted(ix)

    complement = sorted(set(ix) - set(indices)) 
    vx[complement] = 0.0
    qr[complement] = 0.0 

    if type(prg) in {MethodType,FunctionType}:
        p = prg()
        vx = safe_div(V2,p) 

    while np.any(np.round(V1 - V2,5) == 0.0):
        V2 = V2 + vx + qr 

    return V1,V2 

def not_equals__matrix_whole(M,prg,submat_type):
    assert is_2dmatrix(M) 
    assert type(prg) in {FunctionType,MethodType}
    assert submat_type in DEFAULT_SUBMAT_TYPES

    cx0,cx1 = M.shape[0],M.shape[1] 
    x0,x1 = 0,0 
    p = None
    while x1 < cx1:
        p = (x0,x1)

        q = submatrix__2d(M,p,submat_type)

        # check for uniqueness
        ucheck = q.flatten()
        ulen = len(ucheck)
        ustat = ulen == len(np.unique(ucheck)) 

        if not ustat: 
            mx = modulo_in_range(prg(),[0.05,1.0]) 
            adder = DEFAULT_NOTEQUALS_ADDITIVE(q.shape,mx)

            while not ustat: 
                q += adder 
                ustat = ulen == len(np.unique(q.flatten()))

        x1 += 1
        x0 = int(round(x1 * (cx0 / cx1))) 
        x0 = min([M.shape[0] - 1,x0])

    return M