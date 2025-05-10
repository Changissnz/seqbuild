from intigers.seq_struct import * 

class MDRGen: 

    """
    Integer generator that uses an instance of a <ModuloDecompRepr>. 
    Two modes are available. See the methods `next_gentype1` and `next_gentype2`.

    mdr := ModuloDecompRepr
    prgen := function, pseudo-random integer generator. Call with `prgen()`.
    exclude_neg := bool, include negatives in `AffineFitSearch`. 
    gentype := 1|2
    gt2_rcswitch := bool, (see comment for variable)  
    gt2_sel1 := bool, (see comment for variable)  
    gt2_sel2 := bool, (see comment for variable)  
    gt2_seed_in_output := bool, (see comment for variable)  
    """
    def __init__(self,mdr,prgen,exclude_neg:bool,gentype=1,\
        gt2_rcswitch:bool=False,gt2_sel1:bool=False,gt2_sel2:bool=False,\
        gt2_sel3:bool=True,gt2_seed_in_output:bool=True):

        assert type(mdr) == ModuloDecompRepr
        assert type(exclude_neg) == bool 
        assert gentype in {1,2} 

        self.mdr = mdr
        self.mdr2 = None 
        self.prg = prgen 
        self.exclude_neg = exclude_neg
        self.gentype = gentype

        # generator type 1 vars 
        self.cache = [] 

        # generator type 2 vars 
        self.seed2seq = defaultdict(np.ndarray) 
            # can switch between row and column (default is row) 
        self.gentype2_rcswitch = gt2_rcswitch
            # chooses index of row or column 
        self.gentype2_prindex_selector1 = gt2_sel1 
            # possible to switch midway during drawing 
            # sequences from an integer seed 
        self.gentype2_prindex_selector2 = gt2_sel2 
            # sequential selection of integer seed in seed cycle 
        self.gentype2_prindex_selector3 = gt2_sel3 
            # seed value is in output 
        self.gentype2_seed_in_output = gt2_seed_in_output
            # seed cycle used to make `seed2seq`
        self.gentype2_seedcycle = [] 
            # container to draw values for output 
        self.gentype2_cache = [] 
            # indices 
        self.gentype2_intseed_draw = None 
        self.gentype2_isrow_draw = True 
        self.gentype2_rc_index = 0 


    def __next__(self): 
        if self.gentype == 1: 
            return self.next_gentype1() 
        return self.next_gentype2() 

    """
    outputs the next element in `cache`. If `cache` is 
    empty, generates another sequence with integer seed 
    i, the next output from `prg`. 
    """
    def next_gentype1(self): 
        if len(self.cache) == 0: 
            r = self.mdr.reconstruct() 
            self.mdr.reset_first(np.int32(self.prg()))
            self.cache.extend(r) 
        q = self.cache.pop(0) 
        return q 

    #--------------- methods for generator type 2 -----------------# 
    
    # TODO: test 
    """
    """
    def next_gentype2(self): 
        #print("intseed draw: ",self.gentype2_intseed_draw)
        #print("rc index: ",self.gentype2_rc_index)

        if len(self.gentype2_cache) == 0: 
            if len(self.seed2seq) == 0: 
                return None 

            self.gentype2_selector_delta() 
            q = self.gentype2_seedcycle[self.gentype2_intseed_draw]
            
            qx = None 
            if self.gentype2_isrow_draw:
                qx = np.copy(self.seed2seq[q][self.gentype2_rc_index,:])
                self.seed2seq[q] = np.delete(self.seed2seq[q],\
                    self.gentype2_rc_index,axis=0)
            else: 
                qx = np.copy(self.seed2seq[q][:,self.gentype2_rc_index])
                self.seed2seq[q] = np.delete(self.seed2seq[q],\
                    self.gentype2_rc_index,axis=1)

            if len(self.seed2seq[q]) == 0: 
                del self.seed2seq[q] 
            
            if self.gentype2_seed_in_output:
                qx =  np.insert(qx,0,q) 
            self.gentype2_cache = qx 

        value = self.gentype2_cache[0] 
        self.gentype2_cache = np.delete(self.gentype2_cache,0)
        return value 

    def gentype2_selector_delta(self): 
        assert len(self.gentype2_cache) == 0 

        # case: start 
        if type(self.gentype2_intseed_draw) == type(None): 
            # not sequential, choose index based on `prg`
            if not self.gentype2_prindex_selector3: 
                q = self.prg()
                q = q % len(self.gentype2_seedcycle)
                self.gentype2_intseed_draw = q 
            else: 
                self.gentype2_intseed_draw = 0 
        else: 
            # subcase: choose another seed only if the current one 
            #           is exhausted. 
            q = self.gentype2_seedcycle[self.gentype2_intseed_draw]

                # not exhausted 
            if q in self.seed2seq: 
                    # subsubcase: choose another intseed 
                if self.gentype2_prindex_selector2: 
                    q = self.prg() % len(self.gentype2_seedcycle) 
                    self.gentype2_intseed_draw = q 
                # exhausted, choose another intseed 
            else:
                    # delete exhausted intseed  
                self.gentype2_seedcycle.pop(self.gentype2_intseed_draw) 
                if not self.gentype2_prindex_selector3: 
                    q = self.prg() % len(self.gentype2_seedcycle) 
                else: 
                    q = 0 
                self.gentype2_intseed_draw = q 

        if self.gentype2_rcswitch: 
            stat = bool(self.prg() % 2) 
            self.gentype2_isrow_draw = stat 
        
        if self.gentype2_prindex_selector1: 
            dx = 0 if self.gentype2_isrow_draw else 1 
            intseed = self.gentype2_seedcycle[self.gentype2_intseed_draw]
            mx = self.seed2seq[intseed].shape[dx]
            self.gentype2_rc_index = self.prg() % mx 
        return

    def load_int_seed(self,i): 
        assert type(i) in {int,np.int32}
        if self.gentype == 1: 
            self.mdr.reset_first(np.int32(i))
        else: 
            self.mdr2.reset_first(np.int32(i))
        
    """
    generates sequences until no new sequence for any integer seed is outputted. 
    """
    def generate_sequence__type_novelgen(self,intseed_cycle,max_sequences:int,\
        clear_prev_info:bool=True):
        assert type(intseed_cycle) == list and \
            len(intseed_cycle) > 0
        assert len(set(intseed_cycle)) == len(intseed_cycle)
        assert self.mdr.first not in intseed_cycle, "seed {}, have {}".format(self.mdr.first,intseed_cycle)
        assert type(max_sequences) == int or type(max_sequences) == type(None) 
        max_sequences = max_sequences if type(max_sequences) != type(None) else np.inf 
        assert max_sequences > 0 

        if clear_prev_info: 
            self.seed2seq.clear() 
            self.gentype2_seedcycle.clear() 

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
        
        self.gentype2_seedcycle = list(np.array(intseed_cycle,dtype=np.int32)) 
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

    def sequence_size(self): 
        s = 0 
        for v in self.seed2seq.values(): 
            s += v.shape[0]
        return s 

    def display_s2s(self):
        for k,v in self.seed2seq.items(): 
            print("seed: ",k)
            print("sequences: ")
            print(v) 
            print() 
    