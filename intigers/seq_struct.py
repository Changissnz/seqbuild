from .process_seq import * 
from mini_dm.minmax_freq import * 

from copy import deepcopy 

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
            if type(sx) == type(None): 
                break 
            self.d[i] = sx 
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
the modulo operation. Scheme splits an integer sequence into 
contiguous subsequences based on the condition of the (i +1)'th 
value being less than the i'th value. 

`gleqvec_prt`: partition of the indices of integer sequence `l` based 
            on the (i + 1)'th value of `l` being less than the i'th 
            value. 
`afs_prt`: corresponding index spans and (multiple,additive) for `gleqvec_prt1`. 
           The form is 
    [0] (start range, end range) 
    [1] [[(m_0, a_0), [i_00, i_01]],...,[(m_j,a_j),[i_j0,i_j1]]; 
        m_0 := start range, 
        a_j := end range. 
`afs_prt_mod`: corresponding modulus separating every subsequence of 
            `gleqvec_prt`. Number of modulus is |`afs_prt`| - 1. 
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
        self.fin = False 

    def __str__(self):
        s = "\t\tmodulo decomp of {}".format(self.l) + "\n"
        s += "\t* GPRT" + "\n"
        s += str(self.gleqvec_prt) + "\n"
        s += "\t* APRT" + "\n"
        s += str(self.afs_prt) + "\n"
        s += "\t* APRTMOD" + "\n"
        s += str(self.afs_prt_mod) + "\n"
        return s 

    def __eq__(self,md):
        assert type(md) == ModuloDecomp
        return self.afs_prt == md.afs_prt and \
            self.afs_prt_mod == md.afs_prt_mod 

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

        ## ?? 
        if now < 0 or prev < 0:  
            mod_val *= -1 
        if now > 0 and mod_val < 0:
            mod_val *= -1
        elif now < 0 and mod_val > 0: 
            mod_val *= -1 
        ## ?? 

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
        self.fin = True 

    def continuous_merge(self,exclude_neg:bool=True):
        for i in range(len(self.gleqvec_prt)): 
            self.afs_on_subsequence(i,exclude_neg) 
        return

    def premerge_contiguous(self,i,exclude_neg:bool=True): 
        assert i >= 1 and i < len(self.gleqvec_prt) 
        
        cost = 1 + len(self.afs_prt[i-1][1]) + len(self.afs_prt[i][1])# one modulo value 
        q1 = self.gleqvec_prt.pop(i-1)
        q2 = self.afs_prt_mod.pop(i-1)
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
            self.gleqvec_prt.insert(i-1,q1)
            self.afs_prt_mod.insert(i-1,q2)
            return False 

class ModuloDecompRepr: 

    def __init__(self,md,reconstruct_type:int=1):
        assert reconstruct_type in {1,2}  
        assert issubclass(type(md),ModuloDecomp)
        self.load(md)
        self.rtype = reconstruct_type 
        return 
    
    @staticmethod 
    def from_list(l,exclude_neg:bool=True,md_type:int=1):
        assert md_type in {1,2}

        seq = IntSeq(l) 

        if md_type == 1: 
            md = ModuloDecomp(seq)
            md.merge(exclude_neg)
        else: 
            md = ModuloDecompV2(seq,exclude_neg)
        mdr = ModuloDecompRepr(md,md_type) 
        return mdr 

    def load(self,md):
        self.gleqvec_prt = md.gleqvec_prt
        self.afs_prt = md.afs_prt 
        self.afs_prt_mod = md.afs_prt_mod 
        self.first = md.l.l[0] 
        return

    def reconstruct_(self,first): 
        i = 0 
        l = [first]
        prev_ma = None  
        while i < len(self.gleqvec_prt): 
            q = self.afs_prt[i]
            if i > 0: 
                assert type(prev_ma) != type(None) 
                mod_val = self.afs_prt_mod[i-1] 
                val = (l[-1] * prev_ma[0] + prev_ma[1]) % mod_val 
                l.append(val)
            
            for x in q[1]:
                ma = x[0]
                for j in range(x[1][0],x[1][1] + 1): 
                    val = l[-1] * ma[0] + ma[1] 
                    l.append(val)
            i += 1 
            if len(q[1]) > 0: 
                prev_ma = q[1][-1][0] 

        return l

    def reconstruct_v2_(self,first):
        i = 0 
        l = [first]
        prev_ma = None  
        while i < len(self.gleqvec_prt): 
            q = self.afs_prt[i]            
            for x in q[1]:
                ma = x[0]
                for j in range(x[1][0],x[1][1] + 1): 
                    val = l[-1] * ma[0] + ma[1] 
                    l.append(val)

            if i < len(self.gleqvec_prt) - 1:  
                mod_val = self.afs_prt_mod[i] 
                l[-1] = l[-1] % mod_val  

            i += 1 
            if len(q[1]) > 0: 
                prev_ma = q[1][-1][0] 

        return l

    def reconstruct(self): 
        if self.rtype == 1:
            return self.reconstruct_(self.first)
        return self.reconstruct_v2_(self.first) 

    def reset_first(self,f): 
        assert type(f) in {int,np.int32}
        self.first = np.int32(f) 