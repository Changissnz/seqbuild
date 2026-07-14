from morebs2.numerical_generator import prg__single_to_int,\
    prg_partition_for_sz__n_rounds,prg_seqsort,prg_decimal,\
    prg_partition_for_float,wrap_ranged_modulo_over_generator,\
    modulo_in_range
from morebs2.matrix_methods import is_valid_range
from types import MethodType,FunctionType
import numpy as np 

DEFAULT_PARTMULTIGEN_SEGSIZE_RANGE = [21,213]

class PartitionedMultiGen: 

    def __init__(self,prg_seq,super_range,segment_size_range=DEFAULT_PARTMULTIGEN_SEGSIZE_RANGE):
        assert len(prg_seq) > 1 
        for x in prg_seq: assert type(x) in {MethodType,FunctionType} 

        super_range = [int(x) for x in super_range]
        assert is_valid_range(super_range,True,False) #or is_valid_range(super_range,False,False) 
        assert is_valid_range(segment_size_range,True,False) 

        segment_size_range = list(segment_size_range)
        if segment_size_range[0] < len(prg_seq): 
            segment_size_range[0] = len(prg_seq) * 2 
        if segment_size_range[1] <= segment_size_range[0]: 
            segment_size_range[1] = segment_size_range[0] + 1 

        self.prg_seq = prg_seq
        # version of `prg_seq` where every prg is constrained by its assigned ranged modulo 
        self.prg_seq_ = [] 

        self.super_range = super_range
        self.segsize_range = segment_size_range
        self.pindex = 0  

        self.c = 0 
        # number of partitioning changes 
        self.pc = 0 
        self.current_segsize = None 
        self.si2gi_map = dict()
        self.assign_partition() 
        return

    def __next__(self): 

        if self.c >= self.current_segsize: 
            self.assign_partition() 

        gi = self.si2gi_map[self.c] 
        G = self.prg_seq_[gi]

        self.c += 1 
        return G() 

    def assign_partition(self): 
        self.pc += 1 

        G = self.prg_seq[self.pindex]
        self.pindex = (self.pindex + 1) % len(self.prg_seq)

        G_ = prg__single_to_int(G) 

        # get segment size and partition it into |prg_seq| parts.
        self.c = 0 
        self.current_segsize = modulo_in_range(G_(),self.segsize_range) 
        g2size_vec = prg_partition_for_sz__n_rounds(self.current_segsize,len(self.prg_seq),G_,0.5,53) 

        # assign every index in the segment to a PRNG index 
        self.si2gi_map.clear()
        seg_indices = prg_seqsort([i for i in range(self.current_segsize)],G_)

        gvec = np.cumsum(g2size_vec) 
        start_index = 0 
        qi = 0 

        for (i,s) in enumerate(seg_indices):
            end = gvec[qi] 

            if i >= end: 
                qi += 1  
            self.si2gi_map[s] = qi 

        # assign every PRNG a non-overlapping sub-range of the super-range 
        F = self.super_range[1] - self.super_range[0]
        var = prg_decimal(G,[0.1,0.85])
        P = list(prg_partition_for_float(F,len(self.prg_seq),G_,var,n=321,rounding_depth=5))
        P.insert(0,self.super_range[0]) 
        P = np.cumsum(P) 

        self.prg_seq_.clear()
        for i in range(len(P) - 1): 
            r = list(P[i:i+2]) 
            G2 = self.prg_seq[i]
            G2_ = wrap_ranged_modulo_over_generator(G2,r)
            self.prg_seq_.append(G2_)

        return 