import numpy as np 
from morebs2.numerical_generator import prg_seqsort
from operator import add,sub,mul,truediv,floordiv
from math import ceil 
import random


def intlist_no_dups_no_zero(il):
    return sorted(list(set(il) - {0}))

"""
prg is used to decide between neg. and pos. 
in the case of 
"""
def intlist_no_dups_no_zero_abs(il,prg,exclude_one:bool=True): 
    il2 = intlist_no_dups_no_zero(il)
    negs = set(il2) - set(np.abs(il2))

    def ex_one(I):
        if exclude_one: 
            I = set(I) - {1,-1}
        return prg_seqsort(list(I),prg)

    if len(negs) == 0: 
        return ex_one(il2)

    il2 = set(il2)
    for n in negs:
        q = prg() % 2

        if q: 
            il2 -= {n}
        else: 
            il2 -= {-n}

    return ex_one(il2)

#--------------------------------------------------------------------

def prgen(s1,s2):
    def prgen_():  
        return random.randint(s1,s2) 
    return prgen_

def stdop_vec(l,operation,cast_type=np.int32):
    ##assert operation in {add,sub,mul,truediv,floordiv}
    assert cast_type in {np.int32,np.float32} 
    assert type(l) == np.ndarray and len(l.shape) == 1
    d = []
    for i in range(1,len(l)):
        try: 
            q = operation(l[i],l[i-1])
        except: 
            q = np.nan   
        d.append(q if not np.isinf(q) else np.nan)
    return np.array(d,dtype=cast_type)

def diffvec(l,cast_type=np.int32): 
    return stdop_vec(l,sub,cast_type)

def divvec(l,div_type=truediv,cast_type=np.float32):
    return stdop_vec(l,div_type,cast_type)

def gleqvec(l,rounding_depth=5): 
    assert type(l) == np.ndarray and len(l.shape) == 1
    d = []
    for i in range(1,len(l)):
        l1,l2 = round(l[i-1],rounding_depth), round(l[i],rounding_depth)
        equals = l1 == l2 
        if equals: 
            d.append(0)
            continue 
        d.append(-1) if l1 > l2 else d.append(1)
    return np.array(d,dtype='int32')

#---------------------------------------------------------------------
#----------------------- vector-to-vector operation 

def modulated_vec_op(v1,v2,op):

    V,V2 = None,None 
    if len(v1) > len(v2): 
        V,V2 = v1,v2
    else:
        V,V2 = v2,v1 

    q = []
    for (i,v) in enumerate(V): 
        i2 = i % len(V2)
        v_ = V2[i2] 
        q.append(op(v,v_))
    return np.array(q) 

def modulated_vecdot(v1,v2,op1,op2):
    V = modulated_vec_op(v1,v2,op1)
    return op2(V)

#---------------------------------------------------------------------

class IntSeq:

    def __init__(self,l):
        self.l = l 
        self.load()

    #------------------- dunder methods 

    def __iter__(self):
        self.index = 0
        return self

    def __next__(self):
        if self.index < len(self.l):
            x = self.l[self.index]
            self.index += 1
            return x
        raise StopIteration

    def __str__(self):
        return str(self.l)

    def __min__(self):
        return min(self.l)
    
    def __max__(self):
        return max(self.l)

    def __add__(self,x):
        assert type(x) in {np.ndarray,np.int32,int}
        L2 = self.l + x 
        return IntSeq(L2)

    def __sub__(self,x):
        return self.__add__(-x)

    def __abs__(self): 
        return IntSeq(np.abs(self.l)) 

    def __len__(self): 
        return len(self.l)

    def __getitem__(self,i):
        if isinstance(i, slice): 
            return self.l.__getitem__(i)    

        if type(i) in {list,np.ndarray}:
            qx = []
            for i_ in i:
                qx.append(self.__getitem__(i_))
            return qx 

        assert i < len(self) 
        return self.l[i]  

    def __eq__(self,l): 
        assert type(l) == np.ndarray
        if len(l) != len(self): return False
        return np.all(l == self.l) 

    #---------------------------------- 


    def load(self):
        if type(self.l) == type(None): 
            self.l = np.array([],"int32")
        else: 
            self.l = np.array(self.l,"int32")
            assert len(self.l.shape) == 1 
        return 

    # TODO: test 
    def element_indices(self,elements):
        elements = set(elements) 
        i = set()
        for x in elements: 
            i1 = set(np.where(self.l == x)[0])
            i |= i1 
        return i 

    def remove_elements(self,elements):
        indices = self.element_indices(elements)
        self.remove_element_indices(indices)
        return 

    def remove_element_indices(self,indices): 
        l2 = [l_ for (i,l_) in enumerate(self.l) if i not in indices] 
        self.l = np.array(l2,dtype=np.int32)

    def append(self,q): 
        if type(q) in [int,np.int32]: 
            q = np.int32(q)
            self.l = np.append(self.l,q)

        for q_ in q: 
            self.l = np.append(self.l,np.int32(q_))
        return

    """
    upper-right hand triangle of operation-output values 
    """
    def optri(self,operation,cast_type):

        # get the first 
        dvec = stdop_vec(self.l,operation,cast_type) 
        l = len(dvec)
        if l == 0: 
            print("[!] NONE.")
            return

        x = np.zeros((l,l),dtype=cast_type) 
        x[0] = np.copy(dvec)
        lx = l - 1 

        while lx > 0: 
            dvec = stdop_vec(dvec,operation,cast_type) 
            dx = l - len(dvec) 
            nx = np.zeros((dx,))
            dvec2 = np.append(nx,dvec)
            x[dx] = np.copy(dvec2) 
            lx -= 1
        return x 

    '''
    difference triangle 
    '''
    def difftri(self,cast_type=np.int32): 
        return self.optri(sub,cast_type)

    '''
    dividor (multiple) triangle 
    '''
    def divtri(self,div_type=truediv,cast_type=np.float32): 
        return self.optri(div_type,cast_type)

    # TODO: test 
    def match_map(self,is2,match_func): 
        assert type(is2) == IntSeq 

        D = dict() 
        for l in self.l:
            D[l] = match_func(l,is2.l)
        return D 

    def diffcat_vec(self,seg_length:float,start_value = None):
        mx0 = self.__min__() if type(start_value) \
            == type(None) else start_value

        lx = []
        for x in self:
            q = abs(x - mx0)
            l = int(ceil(q / seg_length))
            lx.append(l) 
        return lx 

