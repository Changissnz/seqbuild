from types import MethodType,FunctionType
from numpy import add,ufunc,array 

"""
outputs a function for every `next` call, according to 
the counter `i` modulated by |prngs|. 
"""
class ModPRNGOutputter:

    def __init__(self,prngs):
        for x in prngs: 
            assert type(x) in {FunctionType,MethodType}, "got {},{}".format(x,type(x))
        
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

"""
Extension of class<BasicIterableContainer>. After every full iteration 
of `l`, the base vector, updates `l` to 
`l` := f(`l`,v). 
"""
class NoiseAddedIterableContainer(BasicIterableContainer): 

    def __init__(self,L,v,f=add):
        super().__init__(L,cycle_on=True)
        assert type(f) in {FunctionType,MethodType,ufunc}

        self.l = array(L)
        self.v = v 
        self.f = f  
        self.started = False 
        return 

    def __next__(self): 
        if self.started and self.i == 0: 
            self.l = self.f(self.l,self.v)

        self.started = True 
        return super().__next__() 

def prg__iterable(L,cycle_on:bool=True):
    bic = BasicIterableContainer(L,cycle_on)
    return bic.__next__ 