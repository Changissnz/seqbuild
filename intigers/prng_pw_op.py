from types import MethodType,FunctionType

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
    