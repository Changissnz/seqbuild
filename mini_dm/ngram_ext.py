from .ngram import * 
from morebs2.matrix_methods import vector_to_string,string_to_vector,cr 
from collections import Counter 

DEFAULT_CYMETRIC_SUFFIX_MEASURE_TYPES = ["seq","set","valfreq","l-seq"]

"""
output sequence is by sorted unique values: 
    [v0,f0,v1,f1,...]
"""
# NOTE: duplicate? 
def seq_to_value_frequency_seq(s): 
    q = sorted(set(s))
    C = Counter(s)
    S = [] 
    for q_ in q: 
        S.extend((q_,C[q_]))
    return S 
    

"""
Built on top of n-gram partitioning. A pre/post map for a
one-dimensional sequence S connects every n0-gram to the 
next disjoint n1-gram, n0 the size of the prefix and n1 
that of the suffix. 

NOTE: the n-gram partitioning goes exactly to the end, without 
    cycling (first portion at the end of S, second portion at the 
    beginning). 

A pre/post map has the (key,value) variable of  
        prefix -> start index of occurrence -> postfix. 

"""
class SequencePrePostMap: 

    def __init__(self,S,pre_size,post_size): 
        assert type(pre_size) == type(post_size) == int 
        assert pre_size > 0 and post_size > 0 
        assert pre_size + post_size <= len(S), "pre + post sizes cannot exceed sequence length"

        self.S = S 
        self.pre_size = pre_size 
        self.post_size = post_size 

        # main variable 
        #       prefix -> start index of occurrence -> postfix. 
        self.M = dict() 

        self.c = 0 

        self.preproc() 

        self.fin_stat = False 

    def __getitem__(self,k):
        if isinstance(k, slice): 
            assert False, "invalid key"

        if k not in self.M: 
            return None 
        return self.M[k] 

    def preproc(self): 
        self.grammer = NGrammer2D(self.S,self.pre_size + self.post_size,None) 
        self.end_index = len(self.S) - (self.pre_size + self.post_size) 
        return 

    def process(self): 
        while not self.fin_stat: 
            next(self) 

    def __next__(self): 
        if self.fin_stat: 
            return 

        if self.c > self.end_index:
            self.fin_stat = True 
            return 

        x = next(self.grammer).flatten() 
        pref,suff = x[:self.pre_size],x[self.pre_size:] 

        pref_ = vector_to_string(pref,cr)
        if pref_ not in self.M: 
            self.M[pref_] = dict() 

        self.M[pref_][self.c] = suff 
        self.c += 1 
        return 

    def key_to_seq_by_type(self,k,t): 
        assert t in DEFAULT_CYMETRIC_SUFFIX_MEASURE_TYPES

        if t == "seq": 
            return self.key_to_seq(k) 
        
        if t == "valfreq": 
            return self.key_to_seq__valuefreq_form(k,True)  

        return self.key_to_seq__set_form(k) 

    def key_to_seq_(self,k,f): 
        if k not in self.M: return []

        L = self.M[k] 
        T = sorted(L.keys()) 
        S = [] 
        for t in T: 
            x = f(L[t]) 
            if type(x) in {np.ndarray,list}: 
                S.extend(x) 
            else: 
                S.append(x)  
        return S 

    def key_to_seq(self,k): 
        return self.key_to_seq_(k,lambda x: x) 

    def key_to_seq__valuefreq_form(self,k,label_suffixes:bool=False): 

        if k not in self.M: return []

        L = self.M[k] 
        T = sorted(L.keys()) 

        S = [] 
        suffixes = sorted([vector_to_string(v,cr) for v in L.values()]) 
        suffixes_ = sorted(set(suffixes)) 

        if label_suffixes: 
            lmap = {s:i for (i,s) in enumerate(suffixes_)}  

            def g(suffix): 
                return lmap[suffix]

        else: 
            def g(suffix): 
                return suffix 

        vf_seq = seq_to_value_frequency_seq(suffixes) 
        vf_seq_ = [] 

        for (i,v) in enumerate(vf_seq): 

            if i % 2 == 0: 
                v = g(v) 
                vf_seq_.append(v) 
            else: 
                vf_seq_.append(v)
        return vf_seq_ 

    def key_to_seq__set_form(self,k): 
        return self.key_to_seq_(k,lambda x: sorted(set(x))) 

    def key_to_seq__labelled(self,k): 

        if k not in self.M: return []

        L = self.M[k] 
        T = sorted(L.keys()) 

        S = [] 
        suffixes = sorted(set([vector_to_string(v,cr) for v in L.values()])) 
        lmap = {s:i for (i,s) in enumerate(suffixes)} 

        def g(suffix): 
            return lmap[vector_to_string(suffix,cr)] 
        return self.key_to_seq_(k,g) 