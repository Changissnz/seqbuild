"""
file contains code to operate super-partitions
"""
from morebs2.graph_basics import flatten_setseq
from morebs2.numerical_generator import prg_seqsort
from copy import deepcopy

"""

A super-partition S for a sequence of elements (e_1,e_2,...e_n)
consists of unique subsets of these elements, such that S 
contains all elements (e_1,...,e_n) and an element e_i can 
exist in more than 1 subset. 

Structure is an operator that uses a PRNG to sort a super-partition
`sp`. Then it calculates the quickest full partition from `sp` under 
that under a new ordering. 
"""
class SuperPartitionOp:

    def __init__(self,sp,prg):
        assert type(sp) == list
        self.sp = sp 
        self.indices = [i for i in range(len(self.sp))] 
        self.unique_elements = flatten_setseq(self.sp) 

        self.prg = prg 
        return
    
    """
    main method 
    """
    def one_partition(self):
        self.resort() 
        return self.quickest_partition() 
    
    def resort(self):
        self.indices = prg_seqsort(self.indices,self.prg)
        self.sp = [self.sp[i] for i in self.indices]
        return

    def quickest_partition(self):
        q = deepcopy(self.unique_elements)
        q_ = set() 
        px = [] 
        ix = [] 
        for (i,s) in enumerate(self.sp):
            if q_ == q: break 

            s_ = s - q_
            if len(s_) == 0: continue
            px.append(s_) 
            ix.append(self.indices[i]) 
            q_ = q_ | s_ 
        return px 