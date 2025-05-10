from intigers.seq_struct import * 

class OpTri45N90Split: 

    def __init__(self,split):
        assert type(split) is tuple and len(split) == 2 
        assert len(split[0]) > 0 
        self.split = split 
        self.derivative_seqs = dict() 

    def j_derivative_seq(self,j,store_seq:bool=True): 
        assert j < len(self.split[0]) + 1 

        # case: already stored in memory 
        if j in self.derivative_seqs: 
            return self.derivative_seqs[j]

        k = len(self.split[0]) - 1
        start_seq = deepcopy(self.split[1])
        start_seq.insert(0,self.split[0][-1])

        if store_seq:
            self.derivative_seqs[k+1] = deepcopy(start_seq)

        while k >= j: 
            s = self.split[0][k -1] 
            seq = [s]
            for ss in start_seq: 
                seq.append(seq[-1] + ss) 
            
            if store_seq: 
                self.derivative_seqs[k] = deepcopy(seq) 
            start_seq = seq 
            k -= 1 

        return deepcopy(start_seq)

