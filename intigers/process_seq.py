"""
file contains methods to aid in calculating analytical values on integer vectors
"""

import numpy as np

def diffvec(l,cast_type=np.int32):
    assert cast_type in [np.int32,np.float32] 
    assert type(l) == np.ndarray and len(l.shape) == 1
    d = []
    for i in range(1,len(l)):
        d.append(l[i] - l[i-1])
    return np.array(d,dtype=cast_type) 

def gleqvec(l,rounding_depth=5): 
    assert type(l) == np.ndarray and len(l.shape) == 1
    d = []
    for i in range(1,len(l)):
        l1,l2 = round(l[i-1],rounding_depth), round(l[i],rounding_depth)
        equals = l1 == l2 
        if equals: 
            d.append(0)
            continue 
        d.append(-1) if l1 > l2 else d.append(1)
    return np.array(d,dtype='int32')

def affine_fit_on_vec():

    return -1  