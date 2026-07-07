
from desi.approximators import * 

DEFAULT_LPSGEN_ADJTYPE_RATE_CHANGE_RANGE = [10,91] 

"""
Subclass of class<LPSValueOutputter>. Initialized with 1-3 PRNGs, all 
of them outputting single real numbers at every call. Reformatted 
to output lengths and ranges, before instantiating class<LPSValueOutputter>. 
"""
class LPSGen(LPSValueOutputter):

    def __init__(self,prg1,prg2,prg3,adjustment_type_rate_change_range=DEFAULT_LPSGEN_ADJTYPE_RATE_CHANGE_RANGE): 
        assert type(prg1) in {MethodType,FunctionType} 
        assert is_valid_range(adjustment_type_rate_change_range,True,False)
        assert adjustment_type_rate_change_range[0] > 0

        prg1 = prg__single_to_nvec(prg1,2) 

        if type(prg3) in {MethodType,FunctionType}: 
            prg3 = prg__single_to_range_outputter(prg3) 
    
        super().__init__(prg1,prg2,prg3,\
            adjustment_type=1)

        self.atr_crange = adjustment_type_rate_change_range
        self.change_threshold = self.atr_crange[0]  
        self.c = 0 

    def next_change_threshold(self,v): 
        self.change_threshold = modulo_in_range(int(v),self.atr_crange) 

    def __next__(self):

        x = super().__next__() 

        self.c += 1 
        if self.c >= self.change_threshold: 
            self.c = 0 
            self.adj_type = modulo_in_range((self.adj_type + 1),[1,3])
            self.next_change_threshold(x) 

        return x 