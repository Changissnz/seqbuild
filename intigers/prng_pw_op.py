from types import MethodType,FunctionType
from operator import add,sub,mul
from .extraneous import safe_div 

DEFAULT_PAIRWISE_OPS = [add,sub,mul,safe_div] 

"""
for every pair of input integers (x1,x2), chooses 
an operator O using `prg` and outputs O(x1,x2). 
"""
class PRNGPairwiseOp:

    def __init__(self,operators,prg):
        assert type(operators) == list 
        for x in operators: 
            assert type(x) in {MethodType,FunctionType}
        self.operators = operators
        self.prg = prg

    def __add__(self,x):
        assert type(x) in {MethodType,FunctionType}
        self.operators.append(x)
        return self 

    def apply(self,x,x2):
        i = self.prg() % len(self.operators)
        q = self.operators[i]
        return q(x,x2) 
    
"""
return: 
- a function built from the base function, `pairwise_op`, and 
the operator of the float `weight`, `weight_op`. 
"""
def one_weighted_pairwise_operator(pairwise_op,weight_op,weight,order=0):
    assert order in {0,1,2}  

    if order == 0:
        return lambda x,x2: pairwise_op(weight_op(weight,x),x2)
    elif order == 1:
        return lambda x,x2: pairwise_op(x,weight_op(weight,x2))
    else:  # order == 2
        return lambda x,x2: weight_op(weight,pairwise_op(x,x2))
    
def prg__one_weighted_pairwise_operator(prg,base_op_seq,weight_op_seq):
    weight = prg() 
    op1 = int(prg() % len(base_op_seq))
    op2 = int(prg() % len(weight_op_seq))

    op1 = base_op_seq[op1]
    op2 = weight_op_seq[op2] 

    weight_order = int(prg() % 3) 
    return one_weighted_pairwise_operator(op1,op2,weight,weight_order)