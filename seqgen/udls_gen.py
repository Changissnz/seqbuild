from intigers.poly_output_fitter_ import * 
from types import MethodType,FunctionType
from morebs2.numerical_generator import wrap_sp_modulo_over_generator,prg__single_to_nvec,\
    prg__constant,prg_sequence_to_unique_sequence,sign_preserving_modulo
from morebs2.matrix_methods import is_valid_range
from collections import deque

DEFAULT_UDLSGEN_COEFF_ABSMAX = 556 #2048 
DEFAULT_UDLSGEN_CACHE_CAPACITY = 5555
# range of output value size before <UDLSGen> generates another under-determined linear 
# system of equations to solve. 
DEFAULT_UDLSGEN_NEW_UDLS_RATE1_RANGE = (51,576)
# range of output value size before <UDLSGen> assigns new values for the free variables 
# of its under-determined linear system of equations. 
DEFAULT_UDLSGEN_NEW_UDLS_RATE2_RANGE = (5,18) 

"""
(U)nder-(D)etermined (L)inear (S)ystem (Gen)erator. 
"""
class UDLSGen: 

    def __init__(self,prg,varsize_range,allow_rate1_change:bool,allow_rate2_change:bool,\
        coeff_absmax=DEFAULT_UDLSGEN_COEFF_ABSMAX,draw_values_from_cache:bool=True,\
        verbose:bool=False): 

        assert type(prg) in {MethodType,FunctionType}
        assert is_valid_range(varsize_range,True,False) 

        # NOTE: max variable size is 8, min is 2 
        if not 2 < varsize_range[0] < varsize_range[1] <= 8: 
            start = modulo_in_range(varsize_range[0],[2,9]) 
            end = modulo_in_range(varsize_range[1],[2,9]) 
            if end <= start: 
                end = start + 1 
            varsize_range = (start,end) 

        assert coeff_absmax > 0 and type(coeff_absmax) == int 
        assert type(allow_rate1_change) == type(allow_rate2_change) == bool == \
            type(draw_values_from_cache) == type(verbose) 

        self.prg = prg__single_to_int(prg) 
        self.vs_range = varsize_range
        self.coeff_absmax = coeff_absmax
        self.allow_rate1_change = allow_rate1_change 
        self.allow_rate2_change = allow_rate2_change
        self.draw_values_from_cache = draw_values_from_cache
        self.verbose = verbose 

        self.udls = None 

        self.cache = deque()  

        self.next_shift1,self.next_shift2 = None,None 
        self.c1 = 0 
        self.c2 = 0 

        self.load_udls(False) 
        self.set_next_shift(1) 
        self.set_next_shift(2) 

    def __next__(self):
        V = self.one_input_vector()
        x = self.udls.apply(V)

        self.cache.extend(V) 
        self.cache.append(x) 

        self.c1 += 1 
        self.c2 += 1 

        self.update_by_shift_number(1) 
        self.update_by_shift_number(2) 

        return x 

    def one_input_vector(self): 
        prg_ = wrap_sp_modulo_over_generator(self.prg,self.coeff_absmax) 
        V = prg__single_to_nvec(prg_,self.udls.m.shape[1])()
        return V 

    def load_udls(self,draw_from_cache:bool=False,num_attempts_remaining=3): 

        X,Y = self.new_XY_dataset(draw_from_cache)
        
        if self.verbose: 
            print("\t\t** new dataset")
            print("-- input")
            print(X)
            print()
            print("-- output")
            print(Y) 
            print()
            print("-----")
        
        try: 
            udls = UDLinSysSolver(X,Y) 
            udls.solve()

            self.udls = udls 
            self.assign_values_to_freevars() 
        except: 
            if num_attempts_remaining > 0: 
                return self.load_udls(draw_from_cache,num_attempts_remaining - 1) 
            
        return 

    def new_XY_dataset(self,draw_from_cache:bool=False): 

        num_vars = modulo_in_range(self.prg(),self.vs_range) 
        if num_vars == 2: num_vars = 3 

        num_samples = modulo_in_range(self.prg(),[2,num_vars]) 
        num_x = num_vars * num_samples 
        prg_ = wrap_sp_modulo_over_generator(self.prg,self.coeff_absmax)

        X,Y = None,None 
        if not draw_from_cache:             

            X = prg_unique_sequence(prg_,num_x) 
            X = np.reshape(X,(num_samples,num_vars)) 
            X = np.array(X,dtype=int) 

            Y = prg__single_to_nvec(prg_,num_samples)() 
            Y = np.array(Y,dtype=int)
        else: 
            X,Y = [],[] 

            r = num_x 
            while len(self.cache) > 0 and r > 0: 
                q = self.cache.popleft() 
                X.append(sign_preserving_modulo(q,self.coeff_absmax)) 
                r -= 1 

            while r > 0: 
                X.append(sign_preserving_modulo(prg_(),self.coeff_absmax)) 
                r -= 1 

            r = num_samples 
            while len(self.cache) > 0 and r > 0: 
                q = self.cache.popleft() 
                Y.append(sign_preserving_modulo(q,self.coeff_absmax))  
                r -= 1 

            while r > 0: 
                Y.append(sign_preserving_modulo(prg_(),self.coeff_absmax))  

            prg2 = prg__constant(1)
            X = prg_sequence_to_unique_sequence(prg2,X) 

            X = np.reshape(X,(num_samples,num_vars))
            X = np.array(X,dtype=int) 

            Y = np.array(Y,dtype=int) 

        return X,Y 

    def assign_values_to_freevars(self): 

        fvars = sorted(self.udls.fvars)
        prg_ = wrap_sp_modulo_over_generator(self.prg,self.coeff_absmax) 
        fvar2value_map = {f:prg_() for f in fvars} 
        self.udls.set_freevar_values(fvar2value_map) 
    
    def update_by_shift_number(self,shift_number:int): 

        assert shift_number in {1,2} 

        if shift_number == 1: 
            if self.c1 >= self.next_shift1: 
                self.load_udls(self.draw_values_from_cache)  
                self.c1 = 0 
                self.set_next_shift(1) 
            return 
    
        if self.c2 >= self.next_shift2: 
            self.assign_values_to_freevars() 
            self.c2 = 0 
            self.set_next_shift(2) 

    def set_next_shift(self,shift_number:int): 
        assert shift_number in {1,2} 

        if shift_number == 1: 
            self.next_shift1 = modulo_in_range(self.prg(),DEFAULT_UDLSGEN_NEW_UDLS_RATE1_RANGE) \
                if self.allow_rate1_change else float('inf')
            return 

        self.next_shift2 = modulo_in_range(self.prg(),DEFAULT_UDLSGEN_NEW_UDLS_RATE2_RANGE) \
            if self.allow_rate2_change else float('inf') 
