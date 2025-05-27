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

class PolyEqCondition:

    def __init__(self,cep,target_value:np.int64): 
        assert type(cep) == CEPoly 
        self.cep = cep 
        self.tv = target_value 

    @staticmethod
    def from_POFV2(pofv,x=None):
        x = pofv.x1 if type(x) == type(None) else x 
        assert type(x) in {int,np.int32,np.int64}

        cep = pofv.to_CEPoly() 
        tv = pofv.apply(x)
        return PolyEqCondition(cep,tv) 

    def output(self,x): 
        return self.cep.apply(x) == self.tv 

'''
'''
class POFV2ConditionAutoGen: 

    def __init__(self,prg):
        self.prg = prg 
        return 

    def integerpair_op(self,i1,i2,\
        sibling_range=DEFAULT_NUM_POLYSIBLING_RANGE,coeff_range=DEFAULT_COEFF_RANGE,\
        power_range = DEFAULT_POWER_RANGE,deepcopy_prng:bool=False):
        assert type(i1) in {int,np.int32,np.int64} and \
            type(i2) in {int,np.int32,np.int64}

        q = modulo_in_range(self.prg(),sibling_range)
        for i in range(q):
            #
            pofv = self.one_new_POFV2(i1,i2,\
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
            default_sizemod=True,order_pair=order_pair)
        pofv.solve() 
        return pofv

class UDLSSAutoGen: 

    def __init__(self): 
        return -1 
