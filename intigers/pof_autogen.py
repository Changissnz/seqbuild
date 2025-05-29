"""
file helps with generating instances of classes from file<poly_output_fitter> 
without having to declare as many variable values for instantiation. 

Uses basic data structures such as <IntSeq> and <LCG> to aid in these requirements. 
"""

from .seq_struct import * 
from .poly_output_fitter_ import * 
from morebs2.numerical_generator import prg__n_ary_alternator

# num of options 
DEFAULT_NUM_POLYSIBLING_RANGE = [1,10]

def lcm_times(x0,x1,m): 

    def f(): 
        return np.lcm(x0,x1) * m 
    return f 

'''
'''
class POFV2ConditionAutoGen: 

    def __init__(self,prg):
        self.prg = prg 
        return 

    """
    for 2 integers i1,i2, generates q <PolyOutputFitterVar2> instances 
    according to the method parameters. 
    """
    def integerpair_op(self,i1,i2,\
        sibling_range=DEFAULT_NUM_POLYSIBLING_RANGE,coeff_range=DEFAULT_COEFF_RANGE,\
        power_range = DEFAULT_POWER_RANGE,deepcopy_prng:bool=False):
        assert type(i1) in {int,np.int32,np.int64} and \
            type(i2) in {int,np.int32,np.int64}

        q = modulo_in_range(self.prg(),sibling_range)
        for i in range(q):
            #
            pofv = self.one_new_POFV2(deepcopy(i1),deepcopy(i2),\
                coeff_range,power_range,deepcopy_prng)
            yield pofv

    def one_new_POFV2(self,n0,n1,\
        coeff_range=DEFAULT_COEFF_RANGE,\
        power_range=DEFAULT_POWER_RANGE,\
        deepcopy_prng:bool=True,order_pair:bool=True): 
        
        coeff = modulo_in_range(self.prg(),coeff_range)
        pwr = modulo_in_range(self.prg(),power_range)
        prg = self.prg if not deepcopy_prng else deepcopy(self.prg)
        pofv = PolyOutputFitterVar2(pwr,n0,n1,coeff,prng=self.prg,\
            default_sizemod=False,order_pair=order_pair)
        pofv.solve() 
        return pofv

    # NOTE: caution required for large integers. The exponential
    #       values are not suited for program. 
    def POFV2_to_POFV1_siblings(self,pofv2,sibling_integers): 

        for s in sibling_integers: assert type(s) in {int,np.int32,np.int64} 

        q = [] 
        c = pofv2.apply(pofv2.x1)
        for s in sibling_integers:
            # search for the largest power that base s can use 
            mp = DEFAULT_MAXPOW4BASE(s)
            pwrange = [2,None]
            if mp < 2: 
                print("[!!] large value. be warned.")
            pwrange[1] = min(mp,DEFAULT_POWER_RANGE[1])
            if pwrange[0] == pwrange[1]: pwrange[1] += 1 
            n = modulo_in_range(self.prg(),pwrange)

            pofv1 = PolyOutputFitterVar1(n,s,c,self.prg,default_sizemod=False)
            pofv1.solve() 
            q.append(pofv1) 
        return q 

class UDLSSAutoGen: 

    def __init__(self): 
        return -1 
