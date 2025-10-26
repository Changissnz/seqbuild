from .lcg_v2 import * 
from .tvec import * 
from .extraneous import prg__single_to_int

"""
"""
def is_value_below_modulo_range_length(v,mr):
    assert is_valid_range(mr,True,False) or is_valid_range(mr,False,False)
    return abs(v) <= mr[1] - mr[0]

# NOTE: not the most efficient approach
"""
`set_limit` sets limits for the search (respectively for integer and float 
types). Variable is either int|float and guarantees that scanning algorithm 
does not become an infinite loop. 
"""
def iterative_adjustment_for_direction__modrange(pv,av,sign,modrange,\
    modindex,moddir,set_limit=int):  
    assert is_valid_range(modrange,True,False) or is_valid_range(modrange,False,False)

    # set limits 
    lim = [-float('inf'),float('inf')]
    if type(set_limit) != type(None):
        if set_limit in {int,np.int32,np.int64}: 
            lim[0] = np.iinfo(np.int32).min 
            lim[1] = np.iinfo(np.int32).max 
        elif set_limit in {float,np.float32,np.float64}: 
            lim[0] = np.iinfo(np.float64).min 
            lim[1] = np.iinfo(np.float64).max 

    x1 = modulo_in_range(av,modrange)
    
    # case: do nothing
    if to_trinary_relation(x1,pv) == sign: 
        return modrange

    newrange = [modrange[0],modrange[1]]

    stat = True
    while stat:
        newrange[modindex] += moddir 
        x1 = modulo_in_range(av,newrange)
        stat = to_trinary_relation(x1,pv) != sign
        
        stat2 = lim[0] <= newrange[modindex] <= lim[1] 
        stat = stat and stat2 

    return newrange 

"""
NOTE: Method is specifically for `av` that does not satisfy the 
`is_value_below_modulo_range_length(av,modrange)` relation. 

Adjusts a modulo range `modrange` to M_r so that the output from 
`modulo_in_range(av,M_r)` satisfies the `sign`, -1 for `pv > av` or 
1 for `pv < av`. 

pv := "parent" value, also a "reference" for the binary relation with 
      `av`.
av := actual value, v = m * pv + a. 
sign := wanted binary relation between `av` and `pv`. 
modrange := base modulo range to adjust. 
"""
def multimodular_number__modulo_range_adjustment(pv,av,sign,modrange):
    assert sign in {-1,1}

    if sign == 1: 
        md,mi = 1,1
    else: 
        md,mi = -1,0
    return iterative_adjustment_for_direction__modrange(pv,av,sign,\
        modrange,mi,md,int) 

"""
NOTE: Method is specifically for `av` that satisfies the 
`is_value_below_modulo_range_length(av,modrange)` relation. 

wv := wanted value
av := actual value
cv := current value, the output from `modulo_in_range(av,modrange)`
modrange := base modulo range to adjust. 
"""
def unimodular_number__modulo_range_adjustment(wv,av,cv,modrange):

    sign = to_trinary_relation(wv,cv) 
    assert sign != 0 

    # adjust modrange if necessary
    cv_ = cv 
    mr = [modrange[0],modrange[1]] 
    if sign == 1: 
        if mr[1] < wv:
            mr[1] = wv
            cv_ = modulo_in_range(av,mr)
            mr[1] += wv - cv_ 
        else: 
            if mr[1] == wv: 
                mr[1] = wv + 1 
            d = wv - av 
            mr[0] = d 
    else:
        if mr[0] > wv: 
            mr[0] = wv
            cv_ = modulo_in_range(av,mr)
            mr[1] = mr[1] - cv + wv 
        elif mr[0] == wv: 
            mr[1] = mr[1] - (cv_ - mr[0])
        else: 
            mr[0] -= (cv_ - wv) 
    return mr 

"""
Calculates a new modrange M_r such that for the actual 
value `av`, `modulo_in_range(av,M_r) == pv`. 
"""
def modrange_for_congruence(pv,av,modrange):
    assert type(pv) in {int,np.int32,np.int64}
    assert type(av) in {int,np.int32,np.int64}
    assert is_valid_range(modrange,True,False) or is_valid_range(modrange,False,False)

    if not is_value_below_modulo_range_length(av,modrange): 
        if pv >= 0: 
            mr2 = [0,abs(av)] 

            m = 1 if av < 0 else -1 
            mr2[1] += m * pv 
        else:            
            av_ = -av if av > 0 else av 

            pvs = pv >= 0 
            avs = av >= 0 
            s = -1 if pvs == avs else 1 
            mr2 = [av_ + 1,0] 
            mr2[0] = mr2[0] + s * pv - 1
        return mr2 

    cv = modulo_in_range(av,modrange)
    difference = cv - pv 
    mr2 = [modrange[0],modrange[1]]
    
    if av >= 0: 
        mr2[0] -= difference 
    else: 
        mr2[1] -= difference 
    return mr2 

#--------------------------------------------------------------------- 

DEFAULT_LCGV3_AUTO_TRINARYDELTA_RANGE = (3,9) 
DEFAULT_RANGEDMOD_ABSMAX = 10 ** 6 / 2 - 1

#DEFAULT_LCGV3_RMOD_RANGE = [2,19] 

"""
class is extension of <LCGV2>, a data structure programmed to be able to recognize 
sub-cycles in its input-output map. In addition to this capability in awareness, 
<LCGV3> can `fit` its output sequence according to a trinary vector, set at 
`tv`. The trinary vector V_t denotes the contiguous pairwise relations for the output 
sequence S. If V_t is of length l, then S is of length (l + 1). There are modes to 
the usage of this class. Broadly, there are two modes: `trinary-guided` and `reflective 
modification`. 

In the basic use case, the method `autoset_tv` sets a trinary vector to `tv` by generative 
scheme type 1 or 2. Then the <LCGV3> will output values that satisfy `tv`. When the end of 
`tv` is reached, algorithm re-iterates back at index 0. And algorithm continues doing this 
unless set with another trinary vector for guidance. 

In trinary-guided mode, <LCGV3> automatically changes its trinary vector `tv` at pre-set 
points in its generation of values. Using two positive integers `delta_one` and `delta_two`, 
<LCGV3> additionally uses a counter to keep track of the number of values produced that were 
fitted over the trinary vector. 

Trinary-guided mode allows <LCGV3> to produce output values that are varyingly different 
than its parent class <LCGV2> even when the two have the same parameter values (where that 
applies). 

In reflective modification mode, <LCGV3> applies changes to its pre-output values at the 
same rate specified by `delta_one` and `delta_two`. However, the changes take place via a 
different methodology, described now. It first "shadows" over `trinary_length` 
number of output values from its base LCG. It then calculates a trinary relation vector 
V_t for the vector of these output values. It applies changes to V_t using its `prg`. The 
resulting vector is V_t', set to the class variable `tv`. The <LCGV3> outputs values 
according to this new vector `tv`. 

In some more brief terms, <LCGV3> uses its capability, having an awareness of certain 
trinary relation patterns from its base LCG, to apply changes to its pre-output values to 
be of different quality in trinary relations than what the base LCG produces.
"""
class LCGV3(LCGV2): 

    def __init__(self,start,m,a,n0,n1,sc_size:int,preproc_gd:bool=False,\
        prg=None,super_range=None,exclude_zero__auto_td:bool=False,\
        is_rmod:bool=False,verbose=False):
        super().__init__(start,m,a,n0,n1,sc_size,preproc_gd)
        if type(super_range) != type(None): 
            assert is_valid_range(super_range,False,False) or \
                is_valid_range(super_range,True,False)

        self.tv = None
        self.tvi = None 
        self.r_ = [self.r[0],self.r[1]]
        self.prg = prg
        self.super_range = super_range

        # variables for use with class feature 
        # `auto trinary-vector delta`. 
        self.auto_trinarydelta = False 
        self.gen_type = None 
        self.trinary_length = None 
        self.delta_one = None
        self.delta_two = None
        self.is_delta2_mutable = None 
        self.delta_counter = 0 
        self.delta_counter2 = 0 
            # reserved used for auto trinary delta,rmod feature 
        self.exclude_zero__auto_td = exclude_zero__auto_td

        # variables for use with class feature 
        # `reflective mod`. 
        self.is_rmod = is_rmod

        self.stat__new_trinary = False 
        self.verbose = verbose 

    #------------------------------ main methods, __next__ methods for output values 

    """
    method is reserved for only automatic trinary-vector 
    mode. Repeatedly calls `__next__` method, collecting 
    those output values into a list container B, until 
    tinary-vector is updated (in proportion to 
    `trinary_length * delta_one`). Outputs the sequence B, 
    called a batch, of output values. 
    """
    def next_batch__auto_td(self): 
        assert self.auto_trinarydelta

        bx = []
        stat = True 
        while stat: 
            q = next(self) 
            bx.append(q)
            stat = not self.stat__new_trinary 

        return bx 

    def __next__(self):
        self.stat__new_trinary = False 
        did_fire = self.fired
        s_ = self.s 
        q = super().__next__()

        # used to trim down runtime from exploding values 
        if q > 0:
            q = modulo_in_range(q,[0,99990])
        else: 
            q = modulo_in_range(q,[-99999,0])

        # case: first element, no trinary vec comparison 
        if not did_fire:
            if type(self.super_range) != type(None):
                return modulo_in_range(q,self.super_range)
            return q

        # resets `r` to a smaller modrange if it exceeds large value 
        self.reset_big_modrange() 

        if type(self.tv) == TrinaryVec: 
            assert type(self.prg) != type(None)
            sc = self.tv[self.tvi]
            try: 
                mr = self.adjust_modulo_range(s_,q,sc,prg__single_to_int(self.prg))
            except: 
                print("...strange#2...")
                self.r = deepcopy(self.super_range)
                return modulo_in_range(q,self.super_range)  

            if type(mr) != type(None): 
                self.r = sorted([int(mr[0]),int(mr[1])])
                q = s_ * self.m + self.a
                q = modulo_in_range(q,self.r) 

            self.tvi = self.tvi + 1 

            # case: auto-delta for trinary vector feature of this class 
            if self.tvi == len(self.tv): 
                if self.auto_trinarydelta:
                    self.delta_counter += 1 
                    self.auto_td__update_trinary_vec() 
                if self.is_rmod: 
                    self.delta_counter += 1
                    self.reflective_mod() 
            self.tvi = self.tvi % len(self.tv) 
            self.s = q 

        if type(self.super_range) != type(None):
            q = modulo_in_range(q,self.super_range)

        return q 

    #-------------------------------------- auxiliary methods used for 
    #-------------------------------------- `reflective modification, a technique 
    #-------------------------------------- used to alter output values. 

    """
    sets a new trinary vector for trinary-guided mode. The 
    new trinary vector is a modification of a trinary pattern 
    for the output sequence from the superclass of this instance. 
    The function is fittingly called `reflective modification` 
    because it takes an actual trinary pattern from its original 
    output values, modifies the pattern, and produces values to 
    fit the new trinary vector pattern.  
    """
    def reflective_mod(self): 
        
        # get the next sequence of values. 
        # sequence is to be used as reference 
        # for generator output values.         
        rx = self.trinary_length + 1
        qs = [] 
        for _ in range(rx): 
            q = super().__next__()
            qs.append(q)
        qs = np.array(qs) 

        # convert sequence to trinary relation (to reference 0). 
        dvq = diffvec(qs,cast_type=np.float32)
        rx = to_trinary_relation_v2(dvq,None,zero_feature=True,abs_feature=False)

        # modify the trinary relation 
        DEFAULT_RMOD_SRANGE = (-5000,5000.) 
        tv2 = TrinaryVec.one_instance__v3(len(rx),self.prg,DEFAULT_RMOD_SRANGE) 

        rx3 = rx + tv2.l 
        rx3 = [modulo_in_range(rx3_,(-1,1)) for rx3_ in rx3] 
        rx3 = np.array(rx3,dtype=int) 
        tv = TrinaryVec(rx3) 

        # set as new trinary vector 
        if self.exclude_zero__auto_td: 
            tv = TrinaryVec.omit_value(tv,0,self.prg) 

        self.set_tv(tv)
        self.stat__new_trinary = True
        return 

    #-------------------------------------- auxiliary methods used for 
    #-------------------------------------- cases in trinary-guided mode 

    def adjust_modulo_range(self,reference,new_value,sign_change,prg2):
        # case: sign already satisfied
        if to_trinary_relation(new_value,reference) == sign_change: 
            return None 

        actual_value = self.m * reference + self.a 

        # case: congruence 
        if sign_change == 0: 
            r_ = modrange_for_congruence(reference,actual_value,self.r)
            return r_ 

        ################################################################### 
        
        stat = is_value_below_modulo_range_length(actual_value,self.r) 
        if stat: 
            wv = reference 
            d = modulo_in_range(prg2(),[0,self.r[1] - self.r[0]])
            wv += (d * sign_change)

            try: 
                return unimodular_number__modulo_range_adjustment(wv,actual_value,\
                    new_value,self.r)
            except: 
                print("...strange...")
                self.r = deepcopy(self.super_range)
                return self.r   

        return multimodular_number__modulo_range_adjustment(reference,\
            actual_value,sign_change,self.r)
        
    def reset_big_modrange(self):
        if np.max(np.abs(self.r)) > DEFAULT_RANGEDMOD_ABSMAX:
            self.r = [-int(DEFAULT_RANGEDMOD_ABSMAX / 100000 * 49),\
                int(DEFAULT_RANGEDMOD_ABSMAX / 100000 * 49)] 
        return

    def set_tv(self,tv):
        assert type(tv) == TrinaryVec
        self.tv = tv 
        self.tvi = 0 
        return

    #-------------------------------- methods to automatically update trinary vector 

    def auto_td__update_trinary_vec(self): 
        if self.delta_counter2 >= self.delta_two: 
            if self.is_delta2_mutable: 
                l = modulo_in_range(int(self.prg()),DEFAULT_LCGV3_AUTO_TRINARYDELTA_RANGE)
                self.static_autoset(self.gen_type,l,self.is_delta2_mutable)
            self.delta_counter2 = 0

        if self.delta_counter >= self.delta_one: 
            self.autoset_tv(self.gen_type,self.trinary_length,ext_prg=self.prg)
            self.delta_counter2 += 1
            self.delta_counter = 0 

    def autoset_tv(self,gen_type,l,ext_prg=None):
        assert gen_type in {1,2}
        qx = ext_prg if type(ext_prg) != \
            type(None) else self.__next__ 

        tv = None
        if gen_type == 1:
            L = [-1,0,1]
            bv = L[qx() % 3]
            cr = (qx() % 1000) / 1000.0 
            tv = TrinaryVec.one_instance__v1(bv,l,cr,qx) 
        else: 
            k0 = (qx() % 1000) / 1000.0
            k1 = (qx() % 1000) / 1000.0
            k2 = (qx() % 1000) / 1000.0
            s = k0 + k1 + k2 
            k0,k1,k2 = zero_div0(k0,s),\
                zero_div0(k1,s),zero_div0(k2,s)
            fm = {-1:k0,0:k1,1:k2}
            tv = TrinaryVec.one_instance__v2(l,fm,qx)

        if self.exclude_zero__auto_td: 
            tv = TrinaryVec.omit_value(tv,0,self.prg) 
        self.set_tv(tv)
        self.stat__new_trinary = True 

    # NOTE: requires non-null `prg`. 
    """
    sets variables to static values for `auto_trinarydelta` mode. 
    """ 
    def static_autoset(self,gen_type,l,is_delta2_mutable:bool=False):
        self.auto_trinarydelta = True 

        self.gen_type = gen_type 
        self.trinary_length = l 

        self.delta_one = modulo_in_range(int(self.prg()),DEFAULT_LCGV3_AUTO_TRINARYDELTA_RANGE) 
        self.delta_two = modulo_in_range(int(self.prg()),DEFAULT_LCGV3_AUTO_TRINARYDELTA_RANGE) 
        self.is_delta2_mutable = is_delta2_mutable 
        self.delta_counter = 0 

        self.autoset_tv(self.gen_type,self.trinary_length,self.prg) 
    
    def set_prg(self,prg):
        assert type(prg) in {FunctionType,MethodType}
        self.prg = prg 