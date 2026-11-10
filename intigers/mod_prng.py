from types import MethodType,FunctionType
from numpy import add,ufunc,array 
from copy import deepcopy

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

    def __init__(self,L,v,f=add,op_first:bool=True):
        super().__init__(L,cycle_on=True)
        assert type(f) in {FunctionType,MethodType,ufunc}
        assert type(op_first) == bool

        self.l = array(L)
        self.l_ = deepcopy(self.l)
        self.v = v 
        self.f = f  
        # makes copy of function, in case function is 
        # not constant. 
        self.f_ = deepcopy(f) 

        self.started = False 
        if op_first: self.started = True 
        return 

    def __next__(self): 
        if self.started and self.i == 0: 
            self.l = self.f(self.l,self.v)

        self.started = True 
        return super().__next__() 

    def reset(self,reset_index:bool=True,v=None): 
        assert type(reset_index) == bool 
        if reset_index:
            self.l = deepcopy(self.l_)
            self.i = 0 
            self.f = deepcopy(self.f_)

        if type(v) != type(None):
            self.v = v 
        
        if not reset_index and type(v) != type(None):
            self.l = self.f(self.l,self.v)


def prg__iterable(L,cycle_on:bool=True):
    bic = BasicIterableContainer(L,cycle_on)
    return bic.__next__ 