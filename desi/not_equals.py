from intigers.extraneous import * 
from morebs2.fit_2n2 import * 
from types import MethodType,FunctionType

DEFAULT_ADDITIVE = 10 ** -6 

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

def not_equals__matrix(M,prg,modulo_range= (0.,1.),axes = {0,1}):
    assert axes.issubset({0,1})
    return -1 