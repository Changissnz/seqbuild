from types import MethodType,FunctionType

"""
outputs a function for every `next` call, according to 
the counter `i` modulated by |prngs|. 
"""
class ModPRNGOutputter:

    def __init__(self,prngs):
        for x in prngs: 
            assert type(x) in {FunctionType,MethodType}

        self.prngs = prngs
        self.i = 0

    def __next__(self):
        if self.i >= len(self.prngs):
            self.i = self.i % len(self.prngs)
        q = self.prngs[self.i]
        self.i += 1
        return q 