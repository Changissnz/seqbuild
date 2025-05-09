from intigers.seq_struct import * 

class MDRGen: 

    """
    mdr := ModuloDecompRepr
    prgen := function, pseudo-random integer generator. Call with `prgen()`.  
    """
    def __init__(self,mdr,prgen,exclude_neg:bool):
        assert type(mdr) == ModuloDecompRepr
        self.mdr = mdr
        self.mdr2 = None 
        self.prg = prgen 
        self.exclude_neg = exclude_neg
        self.cache = [] 
        self.seed2seq = defaultdict(np.ndarray) 

    def __next__(self): 
        if len(self.cache) == 0: 
            r = self.mdr.reconstruct() 
            self.mdr.first = np.int32(self.prg())
            self.cache.extend(r) 
        q = self.cache.pop(0) 
        return q 

    def load_int_seed(self,i): 
        assert type(i) in {int,np.int32}
        self.mdr.first = np.int32(i) 
        
    """
    generates sequences until 
    """
    def generate_sequence__type_novelgen(self,intseed_cycle,max_sequences:int,\
        clear_prev_info:bool=True):
        assert type(intseed_cycle) == list and \
            len(intseed_cycle) > 0
        assert len(set(intseed_cycle)) == len(intseed_cycle)
        assert self.mdr.first not in intseed_cycle 
        assert type(max_sequences) == int and max_sequences > 0 

        if clear_prev_info: 
            self.seed2seq.clear() 

        self.mdr2 = deepcopy(self.mdr) 
        terminated_seeds = set() 
        intseed_cycle.insert(0,self.mdr2.first) 
        i = 0 
        c = 0 
        while len(terminated_seeds) < len(intseed_cycle) and \
            c < max_sequences: 
            r, stat = self.add_sequence__type_novelgen_(intseed_cycle[i])

            if not stat: 
                terminated_seeds = terminated_seeds | {intseed_cycle[i]}
            else: 
                c += 1 

            intsq = IntSeq(r) 
            md = ModuloDecomp(intsq) 
            md.merge(self.exclude_neg)
            new_mdr = ModuloDecompRepr(md) 
            self.mdr2 = new_mdr 
            i = (i + 1) % len(intseed_cycle) 
        return 

    def add_sequence__type_novelgen_(self,intseed): 
        self.mdr2.reset_first(intseed) 
        r = self.mdr2.reconstruct()
        q1,q2 = r[0],np.array(r[1:],dtype=np.int32)

        if q1 not in self.seed2seq:
            x = np.empty((0,len(q2)),dtype=np.int32)
        else: 
            x = self.seed2seq[q1]

        stat = np.any(np.all(x == q2, axis=1))
        if stat: 
            return r, False 

        x = np.vstack((x,q2))
        self.seed2seq[q1] = x 
        return r, True 

    def display_s2s(self):
        for k,v in self.seed2seq.items(): 
            print("seed: ",k)
            print("sequences: ")
            print(v) 
            print() 
    