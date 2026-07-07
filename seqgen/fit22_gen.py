from desi.approximators import * 

DEFAULT_FIT22GEN_ADJTYPE_RATE_CHANGE = [13,69] 
DEFAULT_FIT22GEN_ADJTYPE_RATE_CHANGE_2 = [9,95] 


class Fit22Gen(Fit22ValueOutputter): 

    def __init__(self,prg1,prg2,prg3,prg4,adjtype1_rate_change_range = DEFAULT_FIT22GEN_ADJTYPE_RATE_CHANGE,\
        adjtype2_rate_change_range= DEFAULT_FIT22GEN_ADJTYPE_RATE_CHANGE_2): 
        assert type(prg1) in {MethodType,FunctionType} 

        assert is_valid_range(adjtype1_rate_change_range,True,False)
        assert adjtype1_rate_change_range[0] > 0

        assert is_valid_range(adjtype2_rate_change_range,True,False)
        assert adjtype2_rate_change_range[0] > 0

        prg1 = prg__single_to_nvec(prg1,2) 

        if type(prg3) in {MethodType,FunctionType}: 
            prg3 = prg__single_to_range_outputter(prg3) 

        super().__init__(prg1,prg2,prg3,prg4,\
            point_conn_type=1,adjustment_type=1)

        self.atr_crange = adjtype1_rate_change_range
        self.atr2_crange = adjtype2_rate_change_range
        self.change1_threshold = self.atr_crange[0]  
        self.change2_threshold = self.atr2_crange[0]  

        self.c1 = 0 
        self.c2 = 0 

    def next_change_threshold(self,v,type_number):
        assert type_number in {1,2}

        if type_number == 1: 
            self.change1_threshold = modulo_in_range(int(v),self.atr_crange)
        else: 
            self.change2_threshold = modulo_in_range(int(v),self.atr2_crange)

    def switch_mode(self,x,type_number):
        if type_number == 1:  
            if self.c1 >= self.change1_threshold: 
                self.c1 = 0 
                self.adj_type = modulo_in_range((self.adj_type + 1),[1,3])
                self.next_change_threshold(x,type_number) 
            return 

        if self.c2 >= self.change2_threshold: 
            self.c2 = 0 
            self.point_conn_type = modulo_in_range((self.point_conn_type + 1),[1,3])
            self.next_change_threshold(x,type_number) 

    def __next__(self):

        x = super().__next__() 

        self.c1 += 1 
        self.c2 += 1 

        self.switch_mode(x,1)
        self.switch_mode(x,2)
        return x 