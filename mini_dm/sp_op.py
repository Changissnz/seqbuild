"""
file contains code to operate super-partitions
"""
from morebs2.graph_basics import flatten_setseq
from morebs2.numerical_generator import prg_seqsort,prg_seqsort_ties
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
    
    def resort_v2(self,reversed:bool=False):
        indexset_pairs = [(i,self.sp[i]) for i in self.indices] 
        vf = lambda x: len(x[1]) 
        nx = prg_seqsort_ties(indexset_pairs,self.prg,vf)

        if reversed: 
            nx = nx[::-1] 
        self.indices = [nx_[0] for nx_ in nx]
        self.sp = [nx_[1] for nx_ in nx] 
    
    """
    `quickest` is equivalent to greedy search, in this case. 
    """
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
        return px,ix 
    
    """
    attempts to calculate a partition with the least number of sets 
    from the super-partition `sp`, sorted by subset size. 
    """
    def smallest_partition(self): 
        self.resort_v2(reversed=True)
        return self.quickest_partition() 
    
    """
    attempts to calculate a partition with the greatest number of sets 
    from the super-partition `sp`, sorted by subset size. 
    """
    def largest_partition(self): 
        self.resort_v2(reversed=False)
        return self.quickest_partition() 