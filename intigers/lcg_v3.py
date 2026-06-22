from .modulo_ops import * 

DEFAULT_LCGV3_TERNARY_SIZE_RANGE = [3,21] 
DEFAULT_LCGV3_TERNARY_DELTA_TIMESTAMP_RANGE = [6,42] 

class LCGV3(LCGV2): 

    def __init__(self,start,m,a,n0,n1,prg,super_range,ternary_size_range=DEFAULT_LCGV3_TERNARY_SIZE_RANGE,\
        ternary_delta_timestamp_range = DEFAULT_LCGV3_TERNARY_DELTA_TIMESTAMP_RANGE,\
        add_noise:bool=False,verbose=False):

        assert n0 < n1 
        assert is_valid_range(ternary_size_range,True,False)
        assert ternary_size_range[0] >= 2 

        assert is_valid_range(ternary_delta_timestamp_range,True,False) 
        assert ternary_delta_timestamp_range[0] >= 2 
        assert type(add_noise) == bool 
        assert type(prg) in {MethodType,FunctionType} 

        if type(super_range) != type(None): 
            assert is_valid_range(super_range,False,False) or \
                is_valid_range(super_range,True,False)
        else: 
            super_range = (n0,n1) 

        self.ts_range = ternary_size_range
        self.tdt_range = ternary_delta_timestamp_range

        sc_size = modulo_in_range(int(start),self.ts_range)
        super().__init__(start,m,a,n0,n1,sc_size,preproc_gd=False)

        self.prg = prg 
        self.super_range = super_range 
        self.add_noise = add_noise

        self.activated_ternary = None 
        self.current_tsize = None  

        self.tdelta =  None 
        self.tdelta_counter = 0
        self.num_ternaries = 0 

        self.new_ternary()
        self.prev = self.s 

    def __next__(self): 

        self.prev = self.s
        q = self.s #super().__next__() 
        q = self.ternary_adjustment(q) 
        q = self.apply_noise(q) 
        self.update_ternary() 
        self.s = q 
        return q 

    def apply_noise(self,q): 
        if not self.add_noise:
            return q 

        return modulo_in_range(q + self.prg(),self.r) 

    def update_ternary(self): 
        self.tdelta_counter += 1 

        # case: no new ternary 
        if self.tdelta_counter < self.tdelta: 
            return 

        # case: new ternary 
        self.new_ternary() 
        self.tdelta_counter = 0 

    def new_ternary(self):

        self.update_vars() 
        q = int(self.s) % 3 

        if q == 0: 
            base_value = modulo_in_range(int(self.prg()),[-1,2]) 
            cx = prg_decimal(self.prg,[0.,1.])
            tv = TrinaryVec.one_instance__v1(base_value,self.current_tsize,cx,self.prg) 
        elif q == 1: 
            k0 = prg_decimal(self.prg,[0.,1]) 
            k1 = prg_decimal(self.prg,[0.,1]) 
            k2 = prg_decimal(self.prg,[0.,1]) 
            s = k0 + k1 + k2 

            k0,k1,k2 = zero_div(k0,s,1/3),\
                zero_div(k1,s,1/3),zero_div(k2,s,1/3)
            
            fm = {-1:k0,0:k1,1:k2}

            tv = TrinaryVec.one_instance__v2(self.current_tsize,fm,self.prg) 
        else: 
            tv = TrinaryVec.one_instance__v3(self.current_tsize,self.prg,self.super_range) 
        
        self.num_ternaries += 1 
        self.activated_ternary = tv 

    def update_vars(self): 
        self.current_tsize = modulo_in_range(int(self.s),self.ts_range)
        self.tdelta = modulo_in_range(int(self.s),self.tdt_range) 
        
        r0 = modulo_in_range(self.prg(),self.super_range)
        r1 = modulo_in_range(self.prg(),self.super_range)
        R = sorted([r0,r1]) 
        
        # case: start equals ends, default start to super-range start  
        if R[0] == R[1]: 
            R[1] = R[1] + (self.super_range[1] - self.super_range[0]) / 2.0 
            R[1] = modulo_in_range(R[1],self.super_range) 
            R = sorted(R)

        #assert R[0] != R[1] 
        if R[0] == R[1]: 
            R[1] += 1 
        self.r = R 

    def ternary_adjustment(self,q): 
        t = next(self.activated_ternary) 
        if t == 0: 
            self.s = self.prev 
            return self.prev 

        if t == -1: 
            d = q - self.r[0] 
            if d == 0: 
                d = self.r[1] - self.r[0] 
        else: 
            d = self.r[1] - q

        x = self.prg() % d 
        o = q + (x * t)

        if not self.r[0] <= o < self.r[1]: 
            o = modulo_in_range(o,self.r) 

        return o 

