from .process_seq import * 
from copy import deepcopy
from mini_dm.matfactor_eval import * 
from morebs2 import aprng_gauge
from morebs2.poly_struct import CEPoly 
from morebs2.numerical_generator import modulo_in_range
from .extraneous import safe_npint32_value
from math import floor 

# function used to help ensure variable does not exceed size expectations
def BASE_CHECK_FOR_EXP(b,exp): 
    try: 
        x = abs(1000 * (b ** exp))
        return x <= abs(np.iinfo(np.int32).min) and x <= np.iinfo(np.int32).max 
    except:
        return False 


def min_base4pow_geq(f,p): 
    assert abs(f) > 0
    assert p > 0 

    f = abs(f) 
    i = 2 
    stat = i ** p < f 

    while stat: 
        i += 1 
        j = i ** p 
        stat = j < f 
    return i 


NPINT32_MAX = max(np.abs([np.iinfo(np.int32).min,np.iinfo(np.int32).max])) 
DEFAULT_MAXBASE4POW = lambda p: min_base4pow_geq(NPINT32_MAX,p) 
DEFAULT_COEFF_RANGE = [-996,1001]
DEFAULT_POWER_RANGE = [2,10]

def DEFAULT_MAXPOW4BASE(b,max_pow=DEFAULT_POWER_RANGE[1]): 
    assert max_pow > 0 
    if b == 1: return max_pow

    i = 2
    while True: 
        if b >= DEFAULT_MAXBASE4POW(i): 
            i -= 1 
            break
        i += 1
        if i > max_pow: 
            return max_pow 
    return i 

class PolyOutputFitterVar1:

    def __init__(self,n,x1,c,prg,default_sizemod:bool=True,adjust_excess:bool=False):
        for q in [n,x1,c]: assert type(q) in {int,np.int64,np.int32}
        assert n >= 1 
        assert x1 != 0 
        if default_sizemod: 
            q = DEFAULT_MAXBASE4POW(n) // 1.5
            x1,x2 = x1 % q + 1,c % q + 1 

        self.n = n 
        self.x1 = np.int32(x1) 
        self.c = np.int32(c)
        self.prg = prg 
        self.adjust_excess = adjust_excess
        self.coeff = np.zeros((self.n + 1,),dtype=np.int32)
        self.index = 0
        self.rem = c 

    """
    main method 
    """
    def solve(self): 
        while self.one_coeff(): 
            continue 

    def is_solved(self): 
        try: 
            return self.apply(self.x1) == self.c
        except: return False 


    def one_coeff(self): 
        if self.n + 1 <= self.index: return False 
        pwr = self.n - self.index
        if self.rem == 0: return False 

        if pwr == 0:
            if abs(self.rem) >= NPINT32_MAX:
                self.coeff[-1] = safe_npint32_value(self.rem)
                if self.adjust_excess:
                    self.c = self.apply(self.x1) 
            else: 
                self.coeff[-1] = self.rem 
            self.rem -= self.rem 
            return True 

        q = ceil(abs(self.rem / (self.x1 ** pwr)))
         
        coeff_range = [-NPINT32_MAX//2,NPINT32_MAX//2] if \
            q >= NPINT32_MAX else [-q,q]
        md = modulo_in_range(self.prg(),coeff_range)
        self.coeff[self.index] = md
        self.rem = self.rem - (md * self.x1 ** pwr)
        self.index += 1 
        return True 

    def apply(self,x): 
        q = 0 
        for i in range(0,self.n+1): 
            j = self.n - i 
            q += self.coeff[j] * x ** i  
        return q 

    def to_CEPoly(self): 
        return PolyOutputFitterVar1.to_CEPoly_(self) 
    
    @staticmethod
    def to_CEPoly_(pofv1): 

        m = np.zeros((0,2),dtype=np.int32) 
        for i in range(0,pofv1.n+1): 
            j = pofv1.n - i 
            q = np.array([pofv1.coeff[i],j])
            m = np.vstack((m,q)) 
        return CEPoly(m) 

"""
Finds an n'th degree polynomial P with all coefficients 
variable except for the n'th power, set to argument 
`coeff`, s.t. P(x1) = P(x2). 

If `prng` is not None, during the search process, algorithm 
chooses a pseudo-random candidate coefficient for every power 
except for 1 and n.
"""
class PolyOutputFitterVar2:

    def __init__(self,n,x1,x2,coeff=1,prng=None,default_sizemod:bool=False,\
        order_pair:bool=True):
        for q in [n,x1,x2]: 
            assert type(q) in {int,np.int32,np.int64} 
        assert n > 1 
        assert x1 != x2 

        if default_sizemod: 
            q = DEFAULT_MAXBASE4POW(n) // 1.5
            x1,x2 = x1 % q + 1,x2 % q + 1
            if x1 == x2: x2 += 6

        if order_pair and x1 > x2: x1,x2 = x2,x1 

        self.n = n 
        self.x1 = np.int64(x1)
        self.x2 = np.int64(x2) 
        assert self.x1 != 0 and self.x2 != 0 
        self.ref = None 
        self.set_poly(coeff) 
        self.prng = prng 

        self.powerdiff_vec() 
        q = np.zeros((2,),dtype=np.int64) 
        self.running_diff = deepcopy(q) 
        q[0] = self.apply(self.x1) 
        q[1] = self.apply(self.x2) 
        self.update_runningdiff(q) 

        self.stat  = True 
        if not self.is_solvable(): 
            print("[??] cannot compute...") 
            self.stat = False 

    def to_CEPoly(self): 
        m = np.zeros((0,2),dtype=np.int32) 
        l = len(self.poly) 

        for (i,x) in enumerate(self.poly): 
            if x != 0: 
                q = np.array([x,l-i],dtype=np.int32) 
                m = np.vstack((m,q)) 
        return CEPoly(m) 

    def set_poly(self,coeff:int=1): 
        self.poly = np.zeros((self.n,),dtype=np.int64) 
        self.poly[0] = coeff
        self.index = 1  

    def powerdiff_vec(self): 
        self.pdvec = np.zeros((self.n,),dtype=np.int64)

        self.ref = self.x1 if self.x1 < self.x2 else self.x2 

        for i in range(1,self.n+1): 
            px1 = self.x1 ** i 
            px2 = self.x2 ** i
            self.pdvec[self.n - i] = abs(px1-px2) 

    def apply(self,x): 
        q = 0 
        for i in range(1,self.n+1): 
            j = self.n - i 
            q += self.poly[j] * x ** i  
        return q 

    def is_solved(self): 
        return self.apply(self.x1) == self.apply(self.x2) 

    def next_coeff_(self,i): 
        #if self.pdvec[i] == 0: return 0 
        q = self.running_diff / self.pdvec[i]
        if q[0] == 0.0: return -q[1]
        return q[0]

    def next_coeff(self,i): 
        q = self.next_coeff_(i) 
        x = floor(q) 

        if self.n - 1 == i: 
            return np.int64(x)

        x2 = 1 if x > 0 else - 1 
        if type(self.prng) == type(None): 
            return np.int64(x + x2) 

        # NOTE: arbitrary range setting 
        if x == 0: x = 1
        if self.prng() % 2: 
            x2 = x2 * -1
        return np.int64(x + (x2 * self.prng() % x))

    def solve(self): 
        while self.index < self.n: 
            q = self.next_coeff(self.index)
            j = self.n - self.index # ?  
            new_diff = np.array([self.x1**j,self.x2 ** j],dtype=np.int64) 
            new_diff = q * new_diff 
            self.update_runningdiff(new_diff) 
            self.poly[self.index] = q 
            self.index += 1 
        return "finnamoto" 

    def resolve(self,pwr,new_coeff): 
        assert pwr > 1 
        j = self.n - pwr 
        assert j >= 0 

        new_coeffvec = np.zeros((self.n,),dtype=np.int64) 
        new_coeffvec[:j] = self.poly[:j] 
        self.poly = new_coeffvec 
        self.poly[j] = new_coeff
        self.running_diff = np.zeros((2,),dtype=np.int64)  

        y1 = self.apply(self.x1)
        y2 = self.apply(self.x2) 
        q = np.array([y1,y2],dtype=np.int64)
        self.update_runningdiff(q) 
        self.index = j + 1 
        self.solve()
        return

    def is_solvable(self): 
        stat = False 
        targ = self.pdvec[0] 

        for i in range(1,self.n): 
            if targ // self.pdvec[i] == targ / self.pdvec[i]: 
                return True 
        return False 

    def update_runningdiff(self,new_diff): 
        self.running_diff = self.running_diff + new_diff 
        q = np.min(self.running_diff) 
        self.running_diff = self.running_diff - q 

"""
Solves under-determined linear system of equations. 
"""
class UDLinSysSolver:

    def __init__(self,m,y): 
        assert type(m) == np.ndarray and type(y) == np.ndarray 
        assert m.ndim == 2 and m.shape[0] != m.shape[1] 
        assert m.shape[0] > 1 and m.shape[1] > 1
        assert m.shape[0] < m.shape[1] 
        assert y.ndim == 1 and y.shape[0] == m.shape[0] 
        self.m = m 
        self.M = None 
        self.y = y 
        self.Y = None 

        self.colstat = np.zeros((self.m.shape[1],),dtype=np.int32) 

        # index -> (varvec,constant) 
        self.eqtn = dict() 
        self.fvars = None 
        self.fvars_ = None 
        self.varvec = None 

        self.rep_order = None 
        self.cancel_order = None 

        self.constat = True 
        self.inconsistent = None 
        self.consistency_check()

        self.mreps = None 

    """
    main method
    """
    def solve(self): 
        if not self.constat: 
            print("[!] inconsistent. program terminated.")
            return 

        self.initial_eval()
        self.cancel() 
        self.postcancel_solve()
        self.missing_reps() 

    """
    pre-solve evaluation of variables (column-wise values). 
    """
    def initial_eval(self): 
        mfe = MatFactorEval(self.m)
        q1,q2 = mfe.identity_eval(rank_type=1)

        self.rep_order = sorted(q1,key=lambda x:q1[x])
        q = [i for i in range(self.m.shape[1])] 
        for q_ in q: 
            if q_ not in self.rep_order: 
                self.rep_order.insert(0,q_) 

        self.cancel_order = self.rep_order[:self.m.shape[0]]
        self.idn_colsets = q2 
        self.colstat = mfe.id_col

    def consistency_check(self):  
        cfc = MatrixConsistencyCheck(self.m,self.y) 
        cfc.check() 
        self.inconsistent = cfc.inconsistent
        self.constat = len(self.inconsistent) == 0 
        return

    #------------- cancellation methods

    def cancel(self): 
        # for rows 
        bis = aprng_gauge.BatchIncrStruct(self.m.shape[0],False,False,\
            subset_size=self.m.shape[0] - 1)
        self.M = deepcopy(self.m)
        self.Y = deepcopy(self.y)

        # iterate through the columns and cancel 
        for c in self.cancel_order: 
            q = np.asarray(next(bis),dtype=int) 
            self.cancel_var_for_rows(self.M,self.Y,c,q)

    def cancel_var_for_rows(self,M,Y,c,rows):
        failed = [] 
        for r in rows:   
            row,y = self.cancel_one_index(M,Y,r,c)
            if type(row) == type(None): 
                failed.append(r) 
                continue 

            M[r] = row 
            Y[r] = y 
        return failed 

    def post_cancelvar_adjust(self,M,Y): 
        for (j,r) in enumerate(M): 
            y_,i,b = self.solve_constant(r,Y[j]) 
            if b != False: 
                self.eqtn[i] = (np.zeros((M.shape[1],)),y_) 
                self.adjust_matrix_with_constant(M,Y,i) 
                continue 

    def cancel_one_index(self,M,Y,r,ci): 
        # iterate through rows and choose a row 
        # with identical 0-indices and 

        # get 0-indices 
        i0 = set(np.where(M[r] == 0)[0])

        lowest_excess_cancel = float('inf')
        r2,y = None,None 
        for r_ in range(M.shape[0]): 
            r2_,y_ = None,None 
            if r_ == r: continue 
            if M[r_,ci] == 0: continue 

            v = np.lcm(M[r,ci],M[r_,ci])
            m0,m1 = v / M[r,ci], v / M[r_,ci]

            r0 = m0 * M[r] 
            r1 = m1 * M[r_] 
            r2_ = r0 - r1 
            y_ = m0 * Y[r] - m1 * Y[r_] 

            if np.isnan(y_) or np.any(np.isnan(r2_)): 
                continue 

            excess,excess2 = self.check_cancel(M[r],r2_,ci) 

            excess = excess | excess2  

            if len(excess) < lowest_excess_cancel: 
                lowest_excess_cancel = len(excess)
                r2 = r2_
                y = y_  
             
        return r2,y 

    def check_cancel(self,r,r1,ci): 
        old0 = set(np.where(r == 0)[0])
        new0 = set(np.where(r1 == 0)[0])
        excess_delta = new0 - old0 - {ci} 
        excess_delta2 = old0 - new0 - {ci} 
        return excess_delta,excess_delta2

    #-------------------- post-cancellation solve 

    def postcancel_solve(self): 
        # iterate through all rows and solve for constants 
        self.apply_constants() 
        fv = self.freevars() 
        self.fvars = fv 
        for ri in range(self.M.shape[0]): 
            self.solve_row_with_freevar(ri,fv) 
        return

    def solve_row_with_freevar(self,ri,freevar): 
        indices = set(np.where(self.M[ri] != 0)[0])
        indices = indices - freevar 
        if len(indices) == 0: 
            return False 
        if len(indices) > 1: return False 

        q = indices.pop()
        rw = self.M[ri] 
        x = rw[q] 
        rw = -1 * rw 
        rw[q] = 0 
        rw = rw / x 
        self.eqtn[q] = (rw,self.Y[ri] / x) 

    def freevars(self): 
        q = set() 
        for rx in self.M: 
            i0 = set(np.where(rx != 0)[0])
            if len(i0) == 0: continue 

            if len(q) == 0: 
                q = i0 
            else: 
                q = q & i0 
        return q 

    def apply_constants(self): 
        for (i,r) in enumerate(self.M): 
            y = self.Y[i]
            y_,c,stat = self.solve_constant(r,y) 
            if stat: 
                self.eqtn[c] = (np.zeros((self.M.shape[1],)),y_) 
                self.adjust_matrix_with_constant(self.M,self.Y,c)  

    def solve_constant(self,r,y): 
        i0 = np.where(r != 0)[0]
        if len(i0) != 1: 
            return None,None,False 
        i0 = i0[0]
        y = y / r[i0]
        return y, i0,True 

    def adjust_matrix_with_constant(self,M,Y,constant_index): 
        assert constant_index in self.eqtn 
        assert np.all(self.eqtn[constant_index][0] == 0) 
        x = self.eqtn[constant_index][1] 

        cx = M[:,constant_index] * x 
        for (i,c) in enumerate(cx): 
            Y[i] = Y[i] - c 
        M[:,constant_index] = 0 

    #------------------------- post-solve plug and output 

    def value_map(self,fvmap): 
        varvec = indexvalue_map_to_vector(fvmap,self.M.shape[1]) 

        for k,v in self.eqtn.items(): 
            # case: constant 
            if np.all(v[0] == 0): 
                fvmap[k] = v[1] 
            else: 
                s1 = np.sum(varvec * v[0])
                s1 = s1 + v[1] 
                fvmap[k] = s1 
        return fvmap 

    def solve_missing_reps(self,varvec): 
        self.M = deepcopy(self.m) 
        self.Y = deepcopy(self.y) 

        eqtn = deepcopy(self.eqtn) 
        self.eqtn.clear() 
        for i,v in enumerate(varvec): 
            if v == 0: continue 
            self.eqtn[i] = (np.zeros((self.M.shape[1],)),v) 
            self.adjust_matrix_with_constant(self.M,self.Y,i)
        self.eqtn = eqtn 

        for r in self.mreps: 
            r2,y,freevars = self.solve_for_rep(r)
            if type(r2) != type(None): 
                self.fvars_ = freevars 
                self.eqtn[r] = (r2,y)
                return 

    # TODO: not fully tested. 
    def solve_for_rep(self,rep):
        most_reps = 0 
        r_,y,new_freevars = None,None,None 
        for (i,r) in enumerate(self.M): 
            x = r[rep]
            if x != 0: 
                xi = set(np.where(r != 0)[0]) 
                q = xi - {rep}
                if len(q) < most_reps: 
                    continue 

                r_ = - r / x
                r_[rep] = 0 
                y = self.Y[i] / x 
                most_reps = len(q)
                new_freevars = q  
        return r_,y,new_freevars
            
    def set_freevar_values(self,fvmap,map_type=1):
        if map_type == 1:  
            assert set(fvmap.keys()) == self.fvars 
            fvmap = self.value_map(fvmap) 
            self.varvec = indexvalue_map_to_vector(fvmap,self.M.shape[1])
        else: 
            assert set(fvmap.keys()) == self.fvars_
            for k,v in fvmap.items(): 
                self.varvec[k] = v 
            fvmap = vector_to_indexvalue_map(self.varvec)
            fvmap = self.value_map(fvmap) 

            self.varvec = indexvalue_map_to_vector(fvmap,self.M.shape[1])  
            for (i,v) in enumerate(self.m): 
                if self.apply(v) != np.round(self.y[i],5): 
                    self.constat = False 
            return 

    def missing_reps(self): 
        col = set([i for i in range(self.M.shape[1])])  
        self.mreps = col - self.fvars 
        return

    def apply(self,vec): 
        assert type(self.varvec) != type(None) 
        return round(np.sum(self.varvec*vec),5)

class LinCombo: 

    def __init__(self,x): 
        assert type(x) == np.ndarray
        self.x = x 

    def __str__(self):
        s = ""
        L = len(self.x) 

        for i in range(len(self.x) - 1): 
            j = L - i
            x = str(self.x[i]) + "X_" + str(j) + " + "
            s += x

        s += str(self.x[-1])
        return s 

    def degree(self): 
        return len(self.x) - 1 

    def apply(self,X):
        assert type(X) in {np.ndarray,list} 
        assert len(X) == len(self.x) - 1, "GOT {} {}".format(len(X),len(self.x))
        return np.dot(X,self.x[:-1]) + self.x[-1]