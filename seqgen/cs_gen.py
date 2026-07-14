from intigers.modulo_ops import * 

DEFAULT_CSHIFTER_GEN_FRANGE = [4,40] 

class CongruenceShifterGen: 

    def __init__(self,prg,super_range,con_shift_frequency_range=DEFAULT_CSHIFTER_GEN_FRANGE,\
        allow_stretch_and_shrink:bool=False): 
        assert type(prg) in {MethodType,FunctionType}
        assert is_valid_range(super_range,True,False) or  is_valid_range(super_range,False,False)
        assert is_valid_range(con_shift_frequency_range,True,False)
        assert type(allow_stretch_and_shrink) == bool 

        self.prg = prg 
        self.super_range = super_range 
        self.cshift_frange = con_shift_frequency_range
        self.allow_sns = allow_stretch_and_shrink
        self.c = 0 
        self.c2 = 0 # for stretch-n-shrink

        # congruence shift 
        self.next_shift = None
        # stretch-n-shrink shift 
        self.next_shift_sns = float('inf')  
        self.set_next_shift() 
        self.prev = 0 
        return

    def __next__(self): 

        x_ = self.prg() 
        x = modulo_in_range(x_,self.super_range) 
        self.c += 1 
        self.c2 += 1 

        self.adjust_modrange__cshift(x_) 
        self.adjust_modrange__sns()

        self.prev = x 
        return x  

    def adjust_modrange__cshift(self,x_):
        if self.c == self.next_shift: 
            if x_ != self.prev: 
                self.super_range = modrange_for_congruence(x_,self.prev,self.super_range)
                if self.super_range[0] >= self.super_range[1]:
                    print("...strange...")
                    self.super_range = sorted(self.super_range)

            self.c = 0
            self.set_next_shift() 

    def adjust_modrange__sns(self):
        if self.c2 == self.next_shift_sns:
            d = prg_decimal(self.prg,[0.25,0.9]) 

            l = (self.super_range[1] - self.super_range[0]) * d
            l = round(l,5)

            s0,s1 = 1,1 

            # stretch 
            if prg_decimal(self.prg,[0.,1.]) >= 0.5:
                s0 = -1 
            # shrink 
            else:
                s1 = -1 

            start = self.super_range[0] + (s0 * l / 2) 
            end = self.super_range[1] + (s1 * l / 2) 
            ##print('before: ',self.super_range)
            self.super_range = (start,end) 
            ##print('after: ',self.super_range)


    def set_next_shift(self):
        self.next_shift = modulo_in_range(int(self.prg()),self.cshift_frange) 

        if self.allow_sns:
            self.next_shift_sns = ceil(modulo_in_range(int(self.prg()),self.cshift_frange) / 2)