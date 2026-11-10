from types import MethodType,FunctionType,BuiltinFunctionType
from numpy import add,ufunc,array 
from copy import deepcopy

from morebs2.numerical_generator import safe_modulo_in_range
from morebs2.matrix_methods import is_valid_range

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

    def __len__(self): 
        return len(self.l)

"""
Extension of class<BasicIterableContainer>. After every full iteration 
of `l`, the base vector, updates `l` to 
`l` := f(`l`,v). 
"""
class NoiseAddedIterableContainer(BasicIterableContainer): 

    def __init__(self,L,v,f=add,op_first:bool=True):
        super().__init__(L,cycle_on=True)
        assert type(f) in {FunctionType,MethodType,ufunc,BuiltinFunctionType}
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

"""
Extension of class<NoiseAddedIterableContainer>. Uses some k pairwise functions,
    f_1,f2,...,f_k, 
in cyclical order, to update the base vector `l`. The new base vector is of length 
w in the range of `l_range`. 
"""
class NoiseAddedIterableContainerV2(BasicIterableContainer): 

    def __init__(self,L,l_range,v,function_vec,op_first:bool=True):

        super().__init__(L,cycle_on=True)
        assert is_valid_range(l_range,True,True) and l_range[0] > 0 
        assert type(function_vec) == list 
        for f in function_vec: 
            assert type(f) in {FunctionType,MethodType,ufunc,BuiltinFunctionType}, "got {}".format(type(f)) 

        assert type(op_first) == bool

        self.l = array(L)
        self.l_ = deepcopy(self.l)
        self.l_range = l_range 

        self.v = v 
        self.func_vec = function_vec

        # makes copy of function_vec, in case function is 
        # not constant. 
        self.func_vec_ = deepcopy(function_vec) 

        self.started = False 
        if op_first: self.started = True 

        # used to store the most previous output 
        self.prev_output = self.l_[-1]
        return 

    def __next__(self): 
        if self.started and self.i == 0: 
            self.update_vector() 

        self.started = True 
        self.prev_output = super().__next__() 
        return self.prev_output

    def reset(self,reset_index:bool=True,v=None): 
        assert type(reset_index) == bool 
        if reset_index:
            self.l = deepcopy(self.l_)
            self.i = 0 
            self.func_vec = deepcopy(self.func_vec_)

            self.prev_output = self.l_[-1]

        if type(v) != type(None):
            self.v = v 

        if not reset_index and type(v) != type(None):
            self.update_vector() 

    def update_vector(self): 
        new_length = safe_modulo_in_range(int(self.prev_output),self.l_range)

        l2 = [] 
        for i in range(new_length): 
            f = self.func_vec[i % len(self.func_vec)] 

            j = i % len(self.l)
            l0 = self.l[j]
            l2.append(f(l0,self.v)) 

        self.l = l2 


def prg__iterable(L,cycle_on:bool=True):
    bic = BasicIterableContainer(L,cycle_on)
    return bic.__next__ 