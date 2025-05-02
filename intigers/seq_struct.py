from .process_seq import * 
from mini_dm.minmax_freq import * 
from copy import deepcopy 

class IntSeq:

    def __init__(self,l):
        self.l = l 
        self.load()

    def load(self):
        if type(self.l) == type(None): 
            self.l = np.array([],"int32")
        else: 
            self.l = np.array(self.l,"int32")
            assert len(self.l.shape) == 1 
        return 

    def __len__(self): 
        return len(self.l)

    def append(self,q): 
        if type(q) in [int,np.int32]: 
            q = np.int32(q)
            self.l = np.append(self.l,q)

        for q_ in q: 
            self.l = np.append(self.l,np.int32(q_))
        return

    """
    upper-right hand triangle of operation-output values 
    """
    def optri(self,operation,cast_type):

        # get the first 
        dvec = stdop_vec(self.l,operation,cast_type) 
        l = len(dvec)
        if l == 0: 
            print("[!] NONE.")
            return

        x = np.zeros((l,l),dtype=cast_type) 
        x[0] = np.copy(dvec)
        lx = l - 1 

        while lx > 0: 
            dvec = stdop_vec(dvec,operation,cast_type) 
            dx = l - len(dvec) 
            nx = np.zeros((dx,))
            dvec2 = np.append(nx,dvec)
            x[dx] = np.copy(dvec2) 
            lx -= 1
        return x 

    '''
    difference triangle 
    '''
    def difftri(self,cast_type=np.int32): 
        return self.optri(sub,cast_type)

    '''
    dividor (multiple) triangle 
    '''
    def divtri(self,div_type=truediv,cast_type=np.float32): 
        return self.optri(div_type,cast_type)

"""
runs frequency count of (multiple,additive) pairs on contiguous pairs in 
integer sequence. 

Decomposes integer sequence into format: ((multiple,additive),index list). 
"""
class AffineFitSearch: 

    def __init__(self,l,exclude_neg:bool,log_revd:bool=False):
        self.afc = AffineFitCandidates(l)
        self.exclude_neg = exclude_neg
        self.log_revd = log_revd 
        self.d = defaultdict(list) 
        self.mmf = None 
        return 

    """
    preproc method 
    """
    def load_all_candidates(self): 
        i = self.afc.i
        
        while True: 
            sx = self.afc.next_candidate_set(self.exclude_neg)
            #print("SX: ", sx) 
            if type(sx) == type(None): 
                break 
            self.d[i] = sx 
            #print("D: ",self.d)
            i = self.afc.i 
        self.load_mmf()
        
    def load_mmf(self):
        self.mmf = MinMaxFreq(self.d,self.log_revd)

    """
    main method
    """
    def count(self): 
        while not self.mmf.fin:
            self.mmf.count_one()
        self.mmf.finalize_count() 

    def frequency_count(self): 
        return self.mmf.sorted_counts

    """
    main method #2
    """
    def default_affine_decomp(self): 
        assert self.log_revd         
        assert len(self.mmf.revd) > 0 

        l = set([i for i in range(1,len(self.afc.l))])

        x = -1 
        sol = [] 
        while len(l) > 0: 
            q = self.mmf.sorted_counts[x]
            x -= 1 
            indices = set(self.mmf.revd[q[0]])
            new_indices = set([i for i in indices if i in l])
            if len(new_indices) == 0: continue 
            l = l - new_indices
            sol.append((q[0],new_indices))
        return sol 

    """
    arguments: 
    - decomp := list, element is ((multiple,additive),index set) 

    return: 
    - list, element is {(multiple,additive),contiguous span [start index,end index)}. 

    EX: 
        input: [((2, 1), {1, 2, 5}), ((3, -25), {4}), ((3, -20), {3})]
        output: [[(2, 1), [1, 2]], [(3, -20), [3, 3]], [(3, -25), [4, 4]], [(2, 1), [5, 5]]]
    """
    @staticmethod 
    def decomp_to_span_fmt(decomp): 
        revd = dict() 

        for x in decomp: 
            for y in x[1]: 
                revd[y] = x[0]
        if len(revd) == 0: return None 

        mx = max(revd)
        spand = []  
        prev = revd[1] 
        spand.append([prev,[1,1]]) 
        for i in range(2,mx+1): 
            now = revd[i]
            if now == prev:
                spand[-1][1][1] += 1 
            else: 
                i1,i2 = i,i
                spand.append([now,[i1,i2]])
            prev = now 
        return spand 

"""
contains a decomposition scheme for integer sequences, based on 
the modulo operation. 
"""
class ModuloDecomp: 

    def __init__(self,l): 
        assert type(l) == IntSeq 
        self.l = l 
        self.gvec = None 
        self.gleqvec_prt = self.gleqvec_partition()
        self.subvec_fit = dict()
        self.afs_prt = [] 
        self.afs_prt_mod = [] 

    """
    partitions integer sequence based on (greater|lesser|equal)-sign change. 
    """
    def gleqvec_partition(self): 
        self.gvec = gleqvec(self.l.l,rounding_depth=5)
        ilist = []
        
        for i in range(1,len(self.gvec)): 
            if self.gvec[i] * -1 == self.gvec[i-1]: 
                if len(ilist) > 0: 
                    if ilist[-1] != i -1: 
                        ilist.append(i)
                else: ilist.append(i)

        ilist.append(len(self.gvec)) 
        return ilist 

    """
    runs AffineFitSearch for i'th subsequence  
    """
    def afs_on_subsequence_(self,i,exclude_neg:bool=True): 
        assert i >= 0 and i < len(self.gleqvec_prt) 

        prev = 0 if i == 0 else self.gleqvec_prt[i-1] + 1 
        now = self.gleqvec_prt[i] + 1
        chunk = self.l.l[prev:now]
        if len(chunk) < 2: 
            return ((prev,now),[])

        afs = AffineFitSearch(chunk,exclude_neg,log_revd=True)
        afs.load_all_candidates()
        afs.count()
        q = afs.default_affine_decomp()
        a = ((prev,now),q)
        return a 

    def afs_on_subsequence(self,i,exclude_neg:bool=True):
        # case: (i > 0)-index 
        if i > 0: 
            mod_val = self.subsequence_modulo_connect(i)
            if type(mod_val) == type(None): 
                return 
            self.afs_prt_mod.append(mod_val)

        a_ = self.afs_on_subsequence_(i,exclude_neg)
        if len(a_[1]) == 0: 
            self.afs_prt.append(a_) 
            return 

        b_ = AffineFitSearch.decomp_to_span_fmt(a_[1])
        for b2 in b_: 
            b2[1][0] += a_[0][0]
            b2[1][1] += a_[0][0]
        a = (a_[0],b_)
        self.afs_prt.append(a)

    """
    calculates the integer value for modulo operation connecting 
    the last element of previous subsequence (i-1) to the first 
    element of subsequence i. 
    """
    def subsequence_modulo_connect(self,i): 
        assert i >= 1 
        if i >= len(self.gleqvec_prt): 
            return None  

        r = self.afs_prt[i-1][1]
        ma = r[-1] 
        m,a,j = ma[0][0],ma[0][1],ma[1][1] + 1 
        prev,now = self.l.l[j-1],self.l.l[j]
        value = prev * m + a 
        mod_val = value - now
        return mod_val

    """
    main method
    """
    def merge(self,exclude_neg):
        self.continuous_merge(exclude_neg) 
        i = 1 
        while i < len(self.gleqvec_prt): 
            stat = self.premerge_contiguous(i,exclude_neg) 
            if not stat: 
                i += 1 

    def continuous_merge(self,exclude_neg:bool=True):
        for i in range(len(self.gleqvec_prt)): 
            self.afs_on_subsequence(i,exclude_neg)
        return

    def premerge_contiguous(self,i,exclude_neg:bool=True): 
        assert i >= 1 and i < len(self.gleqvec_prt) 

        gp = deepcopy(self.gleqvec_prt)
        apm = deepcopy(self.afs_prt_mod) 
        
        cost = 1 + len(self.afs_prt[i-1][1]) + len(self.afs_prt[i][1])# one modulo value 
        self.gleqvec_prt.pop(i-1)
        self.afs_prt_mod.pop(i-1)
        a_ = self.afs_on_subsequence_(i-1,exclude_neg=exclude_neg)  
        b_ = AffineFitSearch.decomp_to_span_fmt(a_[1])
        for b2 in b_: 
            b2[1][0] += a_[0][0]
            b2[1][1] += a_[0][0]
        cost2 = len(b_)
        a = (a_[0],b_)
        if cost2 < cost:
            self.afs_prt[i] = a 
            self.afs_prt.pop(i-1) 
            
            if i - 1 < len(self.afs_prt_mod): 
                mod_val = self.subsequence_modulo_connect(i)
                if type(mod_val) == type(None): 
                    return 
                self.afs_prt_mod.pop(i-1)
                self.afs_prt_mod.insert(i - 1,mod_val)
             
            return True 
        else: 
            self.gleqvec_prt = gp 
            self.afs_prt_mod = apm 
            return False 