from intigers.seq_struct import * 

class OpTri45N90Split: 

    def __init__(self,split,opfunc=add):
        assert type(split) is tuple and len(split) == 2 
        assert len(split[0]) > 0 
        self.split = split 
        self.opfunc = add 
        self.derivative_seqs = dict() 

    def degree(self): 
        return len(self.split[0])

    '''
    main method 
    '''
    def j_derivative_seq(self,j,store_seq:bool=True): 
        assert j < self.degree() + 1 

        # case: already stored in memory 
        if j in self.derivative_seqs: 
            return self.derivative_seqs[j]

        k = len(self.split[0]) - 1
        start_seq = deepcopy(self.split[1])
        start_seq.insert(0,self.split[0][-1])

        if store_seq:
            self.derivative_seqs[k+1] = deepcopy(start_seq)

        while k >= j: 
            s = self.split[0][k -1] 
            seq = [s]
            for ss in start_seq: 
                seq.append(self.opfunc(seq[-1],ss))
            
            if store_seq: 
                self.derivative_seqs[k] = deepcopy(seq) 
            start_seq = seq 
            k -= 1 

        return deepcopy(start_seq)

    # TODO: test 
    def reproduce(self,p45_indices,p90_indices):
        p45,p90 = [],[] 

        for i in p45_indices: 
            assert i in self.derivative_seqs
            p45.extend(self.derivative_seqs[i])
        for i in p90_indices: 
            assert i in self.derivative_seqs
            p90.extend(self.derivative_seqs[i])
        return OpTri45N90Split((p45,p90),self.opfunc)

# TODO: test 
class OpTriFlipDerivation:

    def __init__(self,m,intseed,opfunc,axis:int): 
        assert type(m) == np.ndarray 
        assert m.ndim == 2 
        assert m.shape[0] == m.shape[1]
        assert type(intseed) in {int,np.int32}
        # 0 -> lower right triangle 
        # 1 -> upper left triangle 
        assert axis in {0,1}

        self.mf = np.flip(m,axis=axis)
        self.m_ = np.zeros((self.mf.shape[0],self.mf.shape[0]),dtype=np.int32) 
        self.intseed = intseed 
        self.opfunc = opfunc 
        self.axis = axis 

    def construct_(self,i,intseed):  
        assert i >= 0 and i < self.mf.shape[0]

        if self.axis == 0: 
            subrow = self.mf[0,-(i+1):]
        else: 
            subrow = self.mf[0,i:] 

        q = [intseed] 
        for v in subrow: 
            q.append(self.opfunc(q[-1],v)) 

        if self.axis == 0: 
            i2 = self.mf.shape[0] - 1 - i 
        else:
            q = q[::-1]  
            i2 = i
        self.m_[i2,i2:] = np.array(q)             
        return q 

    def construct(self): 
        intseed = self.intseed 

        for i in range(self.mf.shape[0]): 
            q = self.construct_(i,intseed)
            intseed = q[-1] 

# TODO: test 
class OpTriGen:

    def __init__(self,intseed,m,prg,gentype:int,\
        forward_func=add,backward_func=sub,add_noise:bool=False):
        assert type(intseed) in {int,np.int32}
        assert type(m) == np.ndarray 
        assert m.ndim == 2 
        assert m.shape[0] == m.shape[1]
        assert gentype in {1,2,3} 

        self.intseed = np.int32(intseed) 
        self.m = m 
        self.prg = prg 
        self.gentype = gentype 
        self.ffunc = forward_func
        self.bfunc = backward_func
        self.ots = None 
        self.cache = [] 

    def jagged_split_(self,row_indices,split_index): 
        l = []
        for (i,r) in enumerate(row_indices): 
            assert r >= i 
            l.append(self.m[r,i]) 
        
        p45,p90 = l[:split_index],l[split_index:]
        return (p45,p90)

    def jagged_split_prng(self):

        row_indices = [] 
        for i in range(self.m.shape[0]): 
            j = self.prg() % (i +1) 
            row_indices.append(j)

        split_index = self.prg() % (len(row_indices) +1)
        return self.jagged_split_(row_indices,split_index) 

    def set_jagged_split(self): 
        js = self.jagged_split_prng()
        self.ots = OpTri45N90Split(js,self.ffunc)

    def store_ots_row_(self,i,store_seq=True): 
        q = self.ots.j_derivative_seq(i,store_seq)
        self.cache.extend(q) 
        return

    def store_ots_row(self,store_seq=True): 
        s = (self.prg() % self.ots.degree()) + 1
        self.store_ots_row_(s,store_seq)  