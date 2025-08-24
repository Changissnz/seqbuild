from intigers.mdr_v2 import * 
from intigers.tvec import * 
from intigers.intfactor import * 
from mini_dm.iseq import * 
from intigers.poly_output_fitter_ import * 

DEFAULT_SHADOW_FITTERS = {"mdr","mdrv2","tvec","fvec","optri","pofv1","pofv2"}

class QualVec: 

    def __init__(self,vec,qual,qual_op=sub):
        assert len(vec) >= 2 
        assert is_vector(vec)
        assert qual in {"tvec","fvec","optri"} 

        self.vec = vec 
        self.qual = qual 
        self.qual_op = qual_op 
        self.qvec = None 
        return 

    def qualvec(self): 
        if self.qual == "tvec": 
            V0 = stdop_vec(self.vec,sub,float)
            V = to_trinary_relation_v2(V0,None,zero_feature=True,abs_feature=False)
            self.qvec = TrinaryVec(V) 
        elif self.qual == "fvec":
            isfso = ISFactorSetOps(np.array(self.vec,dtype=int),int_limit=DEFAULT_INT_MAX_THRESHOLD,str_mode_full=True) 
            isfso.factor_count_() 
            factors = isfso.all_factors()
            self.qvec = np.array(sorted(factors),dtype=int) if len(factors) != 0 else None 
            self.qvec_ = deepcopy(self.qvec)
            self.qvec2 = isfso.factors 
        else: 
            self.qvec = stdop_vec(self.vec,self.qual_op,cast_type = float)
        return

    def apply_noise(self,prg): 
        # case: fvec 

        for i in range(len(self.qvec)): 
            self.qvec[i] = self.qvec[i] + prg() 

        if self.qual == "tvec": 
            self.qvec = round_to_trinary_vector(self.qvec,is_distance_roundtype=True)
        return

    def refit__tvec(self):
        rx = [self.vec[0]]  
        for i in range(1,len(self.vec)): 
            ij = to_trinary_relation(self.vec[i],self.vec[i-1])
            if ij != self.qvec[i-1]: 
                if self.qvec[i-1] == 0: 
                    rx.append(rx[-1]) 
                    continue 
                diff = abs(self.vec[i] - self.vec[i-1]) 
                rx.append(rx[-1] + diff if self.qvec[i-1] == 1 else -diff) 
                continue 
            rx.append(self.vec[i])
        return rx 

    def refit__fvec(self): 

        def fit_for_one(x_):
            fs = self.qvec2[x_] 
            fs2 = [] 
            mfs = []
            for f in fs: 
                i = np.where(self.qvec == f)[0] 
                assert len(i) == 1 
                i = i[0] 
                f2 = self.qvec_[i] 
                fs2.append(f2) 

                mfs_ = zero_div0(x_,f) 
                mfs.append(mfs_) 

            fs2,mfs = np.array(fs2),np.array(mfs) 
            return np.mean(mfs * fs2)

        v2 = [] 
        for x in self.vec: 
            v2.append(fit_for_one(x))
        return np.array(v2)

    def refit__optri(self): 
        return -1 

class ShadowFitter:

    def __init__(self,fitter_idn,prg): 
        assert fitter_idn in DEFAULT_SHADOW_FITTERS
        self.fitter_idn = fitter_idn 
        self.fitter = None 
        return 

    def apply_delta(self):
        return 

    def refit_sequence(self):
        return 

class ShadowGen: 

    def __init__(self,fitting_struct):
        return