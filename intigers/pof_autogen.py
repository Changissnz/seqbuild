"""
file helps with generating instances of classes from file<poly_output_fitter> 
without having to declare as many variable values for instantiation. 

Uses basic data structures such as <IntSeq> and <LCG> to aid in these requirements. 
"""

from .seq_struct import * 
from .poly_output_fitter_ import * 

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

    def __init__(self,intseq,prg):
        assert type(intseq) == IntSeq
        self.intseq = intseq 
        self.prg = prg 
        return 

    def integer_pair(self): 
        return -1

    def one_new_POFV2(self,n0,n1): 
        # PolyOutputFitterVar2(n,x1,x2,coeff=1,prng=None)
        return -1 

class UDLSSAutoGen: 

    def __init__(self): 
        return -1 
