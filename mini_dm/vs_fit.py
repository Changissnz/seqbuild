"""
Planned file that works on vector-singleton fits between vectors. 
"""

from morebs2.matrix_methods import is_number

"""
most under or least over; 

takes a delta function to increment point in search 
"""
def mulo___inc1(p_current,p_ref,diff_func,delta_func,term_func,max_iterations:int):
    d = diff_func(p_ref,p_current)
    q = delta_func(p_current)
    stat = term_func(q) 

    return d,q,stat

class AffineDelta: 

    def __init__(self,m,a,ma_order:int): 
        assert is_number(m) 
        assert is_number(a) 
        assert ma_order in {0,1}

        self.m = m
        self.a = a 
        self.ma_order = ma_order
        return
    
    @staticmethod
    def prg__affine_delta(struct_dim,prg):
        return -1 