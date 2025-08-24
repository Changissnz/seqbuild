from intigers.mdr_v2 import * 
from intigers.tvec import * 
from intigers.intfactor import * 
from mini_dm.iseq import * 
from mini_dm.nsfr import * 
from intigers.poly_output_fitter_ import * 
from intigers.extraneous import to_trinary_relation,to_trinary_relation_v2,zero_div0


DEFAULT_SHADOW_FITTERS = {"mdr","mdrv2","tvec","fvec","optri","pofv1","pofv2"}

class QualVec: 

    def __init__(self,vec,qual,qual_op=sub,inverted_qual_op=add):
        assert len(vec) >= 2 
        assert is_vector(vec)
        assert qual in {"tvec","fvec","optri"} 

        self.vec = vec 
        self.qual = qual 
        self.qual_op = qual_op 
        self.inverted_qual_op = inverted_qual_op
        self.qvec,self.qvec_,self.qvec2 = None,None,None
        self.qualvec() 
        return 

    def qualvec(self): 
        if self.qual == "tvec": 
            V0 = stdop_vec(self.vec,sub,np.float32)
            self.qvec = to_trinary_relation_v2(V0,None,zero_feature=True,abs_feature=False)
            #self.qvec = TrinaryVec(V) 
        elif self.qual == "fvec":
            isfso = ISFactorSetOps(np.array(self.vec,dtype=int),int_limit=DEFAULT_INT_MAX_THRESHOLD,str_mode_full=True) 
            isfso.factor_count_() 
            factors = isfso.all_factors()
            self.qvec = np.array(sorted(factors),dtype=int) if len(factors) != 0 else None 
            self.qvec_ = deepcopy(self.qvec)
            self.qvec2 = isfso.factor_map() 
        else: 
            self.qvec = stdop_vec(self.vec,self.qual_op,cast_type = float)
        return

    def apply_noise(self,prg): 

        for i in range(len(self.qvec)): 
            self.qvec[i] = self.qvec[i] + prg() 

        # case: tvec 
        if self.qual == "tvec": 
            self.qvec = round_to_trinary_vector(self.qvec,is_distance_roundtype=True)
            return self.refit__tvec() 
        # case: fvec 
        if self.qual == "fvec":
            return self.refit__fvec() 
        return self.refit__optri()

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
                i = np.where(self.qvec_ == f)[0] 
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
        v2 = [self.vec[0]] 
        for v in self.qvec: 
            v2.append(self.inverted_qual_op(v2[-1],v)) 
        return v2 

DEFAULT_SHADOWGEN_BASESEQ_LENGTH_RANGE = [7,91]

class ShadowGen(NSFileReader): 

    def __init__(self,prg,file_path,fitting_struct,cast_func = cr):
        assert type(prg) in {MethodType,FunctionType}
        assert os.path.exists(file_path)

        self.prg = prg 
        self.file_path = file_path

        file_obj = open(self.file_path,'r') 
        is_periodic = True 
        self.nsfr = NSFileReader(file_obj,cast_func,is_periodic)

        self.queue = [] 
        return

    def __next__(self): 
        if len(self.queue) == 0: 
            self.load_next_fitter()
        return self.queue.pop(0)

    def load_next_fitter(self): 
        L = int(round(modulo_in_range(self.prg(),DEFAULT_SHADOWGEN_BASESEQ_LENGTH_RANGE)))
        S = [] 
        for _ in range(L): 
            S.append(next(self.nsfr)) 
        self.fit_sequence(S)

        q = self.mod_sequence()
        self.queue.extend(q) 

    def fit_sequence(self,S): 
        if self.fitting_struct[:3] == "mdr": 
            t = 1 
            if self.fitting_struct == "mdrv2":  
                md = ModuloDecompV2(s,False) 
                t = 2 
            else: 
                md = ModuloDecomp(S) 
                md.merge(False)  
            self.fitter = ModuloDecompRepr(md,t)
        else: 
            self.fitter = QualVec(deepcopy(S),self.fitting_struct,qual_op=sub)
        return 

    def mod_sequence(self): 

        if type(self.fitter) == ModuloDecompRepr:
            d = int(round(self.prg() % 2)) 

            # case: shift the partition 
            if d: 
                lx = len(self.fitter.afs_prt) 
                lx = lx ** 2 
                s = int(round(self.prg() % lx))
                self.fitter.shift_afs_prt(s) 
            # case: add noise to partition 
            else: 
                self.fitter.noise_to_afs_prt(self.prg,True)
            return self.fitter.reconstruct() 
        return self.fitter.apply_noise(self.prg) 