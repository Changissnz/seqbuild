from .process_seq import * 
from mini_dm.minmax_freq import * 

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
            indices = set(self.mmf.revd[q[0]])
            new_indices = set([i for i in indices if i in l])
            if len(new_indices) == 0: continue 
            l = l - new_indices
            sol.append((q[0],new_indices))
        return sol 