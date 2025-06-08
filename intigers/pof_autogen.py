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

def safe_power_range_for_base(b,power_range=DEFAULT_POWER_RANGE):
    assert power_range[0] > 0 
    assert power_range[1] > power_range[0] 

    mp = DEFAULT_MAXPOW4BASE(b,power_range[1])
    pwrange = [power_range[0],None]
    if mp < power_range[0]: 
        print("[!!] large value. be warned.")
        return None 

    pwrange[1] = min(mp,DEFAULT_POWER_RANGE[1])
    return [pwrange[0],max(pwrange[1],pwrange[0]+1)]

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
        
        nx = max([n0,n1])
        pwrange = safe_power_range_for_base(nx,power_range)

        coeff = modulo_in_range(self.prg(),coeff_range)
        pwr = modulo_in_range(self.prg(),pwrange)
        prg = self.prg if not deepcopy_prng else deepcopy(self.prg)
        pofv = PolyOutputFitterVar2(pwr,n0,n1,coeff,prng=self.prg,\
            default_sizemod=False,order_pair=order_pair)
        pofv.solve() 
        return pofv

    # NOTE: caution required for large integers. Their exponential
    #       values are not suited for program. 
    def POFV2_to_POFV1_siblings(self,pofv2,sibling_integers,adjust_excess=False): 

        for s in sibling_integers: assert type(s) in {int,np.int32,np.int64} 

        q = [] 
        solvestat = [] 
        c = pofv2.apply(pofv2.x1)
        for s in sibling_integers:
            # search for the largest power that base s can use 
            pwrange = safe_power_range_for_base(s)

            n = modulo_in_range(self.prg(),pwrange)
            pofv1 = PolyOutputFitterVar1(n,s,c,self.prg,default_sizemod=False,\
                adjust_excess=adjust_excess)
            pofv1.solve() 
            q.append(pofv1)
            solvestat.append(pofv1.is_solved()) 
        return q,solvestat 

class UDLSSAutoGen: 

    def __init__(self): 
        return -1 
