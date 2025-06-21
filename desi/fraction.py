from intigers.seq_struct import * 
from intigers.extraneous import zero_div0  
from morebs2.numerical_generator import modulo_in_range
from morebs2.matrix_methods import is_valid_range 
from types import MethodType,FunctionType

"""
Outputs float values that originate from quotient outputs of pairs in `intseq`.

index_selector := function|method, parameter-less. Outputs pairs of numbers that act 
                  as 2-dimensional indices. 
length_outputter := function|method, parameter-less. Outputs positive integers. 
range_outputter := function|method, parameter-less. Outputs ordered pairs of 
                    numbers. 
"""
class QValueOutputter: 

    def __init__(self,intseq,index_selector,length_outputter,range_outputter,\
            adjustment_type:int):
        assert type(intseq) == IntSeq 
        assert len(intseq) >= 2 
        assert type(index_selector) in {MethodType,FunctionType}
        assert type(length_outputter) in {MethodType,FunctionType}
        assert type(range_outputter) in {MethodType,FunctionType}
        assert adjustment_type in {1,2}

        self.intseq = intseq
        self.iselector = index_selector
        self.l_out = length_outputter
        self.r_out = range_outputter
        self.adj_type = adjustment_type
        self.s = "" 
        return 

    def change_adj_type(self):
        self.adj_type = 1 if self.adj_type == 2 else 2 

    def __next__(self):
        l = modulo_in_range(self.l_out(),[1,18]) 
        assert type(l) in {int,np.int32,np.int64}

        s_ = self.s[:l]
        self.s = self.s[l:]

        diff = l - len(s_) 
        if diff > 0:
            self.set_next_q()
            s_ = s_ + self.s[:diff]
            self.s = self.s[diff:]

        return self.adjust_integer(np.int64(s_)) 
    
    def set_next_q(self):

        i1,i2 = self.iselector() 
        i1 = i1 % len(self.intseq)
        i2 = i2 % len(self.intseq)
        
        q = str(zero_div0(self.intseq[i1],self.intseq[i2])) 
        x1 = q.find(".")

        if x1 != -1: 
            q = q[x1+1:] 

        x2 = q.find("E")

        if x2 != -1: 
            q = q[:x2]
        self.s = q 
        return q 
    
    def adjust_integer(self,i):
        assert type(i) in {int,np.int32,np.int64}

        mrange = self.r_out()
        assert is_valid_range(mrange,False,False)

        if self.adj_type == 1:
            return self.adjustment_type1(i,mrange)
        return self.adjustment_type2(i,mrange)
    
    def adjustment_type1(self,i,mrange): 
        return modulo_in_range(i,mrange)
    
    def adjustment_type2(self,i,mrange): 
        diff = mrange[1] - mrange[0]
        i = abs(i) 

        while i >= diff:
            i = i / 10
        return mrange[0] + i