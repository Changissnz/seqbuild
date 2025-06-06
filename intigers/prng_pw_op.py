from types import MethodType,FunctionType

class PRNGPairwiseOp:

    def __init__(self,operators,prg):
        assert type(operators) == list 
        self.operators = operators
        self.prg = prg

    def apply(self,x,x2):
        i = self.prg() % len(self.operators)
        q = self.operators[i]
        return q(x,x2) 