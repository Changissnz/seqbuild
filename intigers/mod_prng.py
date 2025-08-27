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

class BasicIterableContainer:

    def __init__(self,L,cycle_on:bool=True):
        assert type(L) == list
        assert len(L) > 0 
        self.l = L
        self.cycle_on = cycle_on
        self.i = 0 

    def __next__(self):
        if self.i >= len(self.l): 
            return None

        q = self.l[self.i] 
        self.i = self.i + 1
        
        if self.cycle_on: 
            self.i = self.i % len(self.l)
        return q 

def prg__iterable(L,cycle_on:bool=True):
    bic = BasicIterableContainer(L,cycle_on)
    return bic.__next__ 