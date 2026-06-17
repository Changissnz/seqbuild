"""
Self-referential algorithms to derive integers from base floats. 
"""
import numpy as np 
from morebs2.matrix_methods import is_number
from intigers.prng_pw_op import * 
from intigers.mod_prng import * 
from copy import deepcopy 
from collections import defaultdict 
from operator import add,sub,mul


DEFAULT_PDD_EXCLUDED_NUMERICAL_CHARS = {"e","E",".","-","+"} 
DEFAULT_PDD_EXCLUDED_OPERATING_NUMBER_TYPES = {float,np.float16,\
    np.float32,np.float64,complex,np.complex64,np.complex128}
DEFAULT_PDD_OPERATING_FLOAT_RANGE = [-10**32,10**32]

#------------------------------------------------------------------------

def float_to_str__type_exclude_EDOT(f): 

    s = str(f) 
    t = "" 
    for s_ in s: 
        if s_ not in DEFAULT_PDD_EXCLUDED_NUMERICAL_CHARS: 
            t += s_ 
    return str(int(t))

#------------------------ methods to select substrings in integers 

def int_substr_at_index_(I_:str,current_index,length): 
    assert type(I_) == str 
    assert type(int(I_)) == int 
    assert type(current_index) == type(length) == int 

    i0 = current_index 
    i1 = (current_index + length) % len(I_) 

    # case: end index cycles back to the prefix 
    if i1 < i0: 
        p0 = I_[i0:] 
        p1 = I_[:i1] 
        return int(p0 + p1), i1  

    # case: start equals end 
    if i0 == i1: 
        return int(I_[i0]), i1 

    # case: proper substring 
    return int(I_[i0:i1]), i1 

def int_split_next__type_sequential_ref(I,current_index,ref_index):  
    assert is_number(I,DEFAULT_PDD_EXCLUDED_OPERATING_NUMBER_TYPES)
    I_ = str(I)
    x = int(I_[ref_index]) 
    return int_substr_at_index_(I_,current_index,x) 

def int_split__type_current_ref(I,current_index):  
    I_ = str(I) 
    x = int(I_[current_index]) 
    return int_substr_at_index_(I_,current_index,x) 

#--------------------------- methods to derive new integers from an integer 

"""
auxiliary method for method<integer_derivation__type_PWO> 
"""
def integer_derivation__type_PWO_selector_(I,current_index,ref_index=None): 

    first = None 
    ##print("LEN {} {} {}".format(len(str(I)),current_index,ref_index)) 
    if type(ref_index) == type(None): 
        ##print("CURRENT INDEX: ",type(current_index))
        first,current_index = int_split__type_current_ref(I,current_index) 
    else: 
        first,current_index = int_split_next__type_sequential_ref(I,current_index,ref_index)
        ref_index = (ref_index + 1) % len(str(I))  

    return first,current_index,ref_index

"""
fetches two integer substrings i0 and i1 from I, and performs 
op(i0,i1) to produce a new integer
"""
# NOTE: does not check if op(i0,i1) can execute
def integer_derivation__type_PWO(I,op,current_index,ref_index=None,skip_counts=0): 
    assert skip_counts >= 0 

    first,current_index,ref_index = integer_derivation__type_PWO_selector_(I,current_index,ref_index)

    for _ in range(skip_counts): 
        _,current_index,ref_index = integer_derivation__type_PWO_selector_(I,current_index,ref_index)

    second,current_index,ref_index = integer_derivation__type_PWO_selector_(I,current_index,ref_index) 
    return round(op(first,second)),current_index,ref_index 

#--------------------------------------------------------------------------

# NOTE: this is a custom safe division function that defaults to 1 + 42/43 when 
#       Python cannot handle the float division of num / denum. 
def safe_div_v2(num,denum): 
    try: 
        return safe_div(num,denum)
    except: 
        return 1 + 42/43 

# +,-,/,* and weighted custom 
DEFAULT_PDD_PAIRWISE_OPS = [add,sub,mul,safe_div_v2,"w.c."] 


DEFAULT_PDD_PRIMARY_MOVES = ["cache op",\
    integer_derivation__type_PWO_selector_,\
    integer_derivation__type_PWO] 

DEFAULT_CACHE_OP_DRAW_SIZE_RANGE = [1,6]  
DEFAULT_PDD_PAIRWISE_OP_WEIGHT_RANGE = (0.,5+54/55-16/17)

DEFAULT_PDD_SELECTOR_SKIP_RANGE = [0,6] 

#--------------------------------------------------------------------------

def integer_to_prg_iterable(I):  
    assert is_number(I,DEFAULT_PDD_EXCLUDED_OPERATING_NUMBER_TYPES)
    S = list(str(I)) 
    S = [int(s) for s in S] 
    return prg__iterable(S)

#--------------------------------------------------------------------------

"""
Iterator takes an initial float f, and outputs floats based on 
selector+delta operations it conducts revolving around f. 

The procedure for outputting floats: 
- choose one of the moves in `DEFAULT_PDD_PRIMARY_MOVES`.
- 
"""
class PRNGDecimalDelta:   

    def __init__(self,f,cache_size=200,enable_base_delta:bool=False):  
        s = str(f) 
        assert len(s) >= 4, "float must have at least four values" 
        assert type(cache_size) == int and cache_size >= 5 
        assert type(enable_base_delta) == bool 

        self.f = f 
        self.cache_size = cache_size
        self.enable_base_delta = enable_base_delta
        self.s = float_to_str__type_exclude_EDOT(f) 

        self.primary_move_list = deepcopy(DEFAULT_PDD_PRIMARY_MOVES) 
        self.pw_op_list = deepcopy(DEFAULT_PDD_PAIRWISE_OPS)

        self.queue = [] 
        self.cindex = 0 
        self.rindex = 0 
        self.start()    

        return

    def choose_pairwise_operator(self): 
        x = DEFAULT_PDD_PAIRWISE_OPS[self.rindex % len(DEFAULT_PDD_PAIRWISE_OPS)] 
        self.rindex = (self.rindex + 1) % len(self.s)  

        if x == "w.c.": 
            return self.generate_weighted_pairwise_operator() 
        return x  

    def generate_weighted_pairwise_operator(self): 
        # fetch four numbers 
        N = [] 
        for i in range(4): 
            d = self.deciding_int() 
            N.append(d) 
            self.rindex = self.rindex + 1 
        prg = prg__iterable(N)
        pw_op = prg__one_weighted_pairwise_operator(prg,DEFAULT_PAIRWISE_OPS,\
            DEFAULT_PAIRWISE_OPS,weight_range= DEFAULT_PDD_PAIRWISE_OP_WEIGHT_RANGE)

        return pw_op 
    
    def __next__(self): 
        self.cindex = self.cindex % len(self.s) 

        action = self.choose_act() 
        i = self.deciding_int() 
        i2 = None 

        if action == integer_derivation__type_PWO_selector_: 
            i2 = self.fetch_substr(i)
        
        elif action == integer_derivation__type_PWO: 
            i2 = self.pairwise_op_with_skipped_selections(i) 
        else: 
            i2 = self.pairwise_ops_on_queue_elements(i) 

        self.queue.append(i2) 
        self.s = self.s[::-1] 
        self.update_queue()     
        return i2 

    #--------------------------------- primary move: substring selection of integerized float 

    """
    obtains a substring from `s` 
    """
    def fetch_substr(self,dint:int): 
        assert is_number(dint,DEFAULT_PDD_EXCLUDED_OPERATING_NUMBER_TYPES)

        add_ref = bool(dint % 2) 
        
        I = int(self.s) 
        I_ = None 
        if add_ref: 
            I_,self.cindex,_ = integer_derivation__type_PWO_selector_(\
                I,self.cindex,self.rindex % len(self.s))  
        else: 
            I_,self.cindex,_ = integer_derivation__type_PWO_selector_(\
                I,self.cindex,None) 
        self.rindex = self.rindex + 1

        return I_

    #--------------------------------- primary move: pairwise operation on two substrings of integerized float 

    def pairwise_op_with_skipped_selections(self,i): 
        skip_counts = modulo_in_range(i,DEFAULT_PDD_SELECTOR_SKIP_RANGE)

        I = int(self.s) 
        pw_op = self.choose_pairwise_operator() 

        ref_index = self.rindex % len(self.s) if I % 2 else None 

        i2,self.cindex,ref_index = integer_derivation__type_PWO(\
            I,pw_op,self.cindex,ref_index=ref_index,skip_counts=skip_counts)

        if type(ref_index) != type(None): 
            self.rindex = ref_index 

        return i2 
        
    #--------------------------------- primary move: q pairwise operations with integer parameter and q elements from queue 

    def pairwise_ops_on_queue_elements(self,i): 
        x = modulo_in_range(i,DEFAULT_CACHE_OP_DRAW_SIZE_RANGE) 
        C = []
        for _ in range(x): 
            d = self.deciding_int() 
            pw0 = self.choose_pairwise_operator()
            C.append((d,pw0)) 
        i2 = i 

        for x0,pw in C: 
            i2 = pw(i2,x0)  
        return i2 

    #------------------------------------- post-calculation of next float, updating queue 

    def update_queue(self): 
        q = -1 
        f2 = 0 
        while len(self.queue) > self.cache_size: 
            r = self.queue.pop(0) 

            if self.enable_base_delta: 
                # approach 1 
                f2 += modulo_in_range(r * q,DEFAULT_PDD_OPERATING_FLOAT_RANGE)   
                q = -1 * q

        if self.enable_base_delta:
            self.update_float(int(f2))

    #--------------------------------- auxiliary methods for calculating the next float 

    def start(self): 
        q = int(self.s[-1]) 

        if q % 2: 
            self.s = self.s[::-1] 
        
        add_ref = int(self.s[-1])
        I = self.fetch_substr(add_ref)
        self.queue.append(I) 

        return 

    def choose_act(self): 
        action = self.primary_move_list[self.rindex % len(self.primary_move_list)] 

        # shift move list by 1 
        Q = self.primary_move_list[1:] 
        self.primary_move_list = Q + [self.primary_move_list[0]]
        return action 

    def deciding_int(self): 
        if len(self.queue) > 0: 
            i = self.rindex % len(self.queue)
            return round(self.queue[i]) 
        else: 
            self.start() 
            return self.deciding_int() 

    def update_float(self,delta): 
        assert is_number(delta,DEFAULT_PDD_EXCLUDED_OPERATING_NUMBER_TYPES) 

        self.f = self.f + delta 
        self.s = float_to_str__type_exclude_EDOT(self.f) 
        self.cindex = self.cindex % len(self.s)