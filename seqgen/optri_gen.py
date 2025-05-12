from intigers.seq_struct import * 

"""
collates values from circumference of an upper-right hand 
triangular matrix. 

00 01 02 03 
   21    10
      20 11
         12
-- clockwise 

00 21 20 12 
   01    11
      02 10
         03
-- counter-clockwise 


          
Argument `start_edge` is the starting edge for the circumference 
values. 
"""
def triangular_matrix_circumference(m,start_edge=0,is_clockwise=True): 
    assert type(m) is np.ndarray 
    assert m.ndim == 2 
    assert m.shape[0] == m.shape[1] 
    assert m.shape[0] > 0 
    assert start_edge in {0,1,2} 

    if m.shape[0] == 1: return np.array(m[0,0],dtype=np.int32)

    if is_clockwise: 
        edge0 = m[0,:] 
        edge1 = m[:,-1] 
        edge2 = np.diag(m)[::-1]
    else: 
        edge0 = np.diag(m) 
        edge1 = m[:,-1][::-1]
        edge2 = m[0][::-1] 

    if start_edge == 0: 
        edges = [edge0,edge1,edge2] 
    elif start_edge == 1: 
        edges = [edge1,edge2,edge0] 
    else:
        edges = [edge2,edge0,edge1] 

    q = edges[0]
    q = np.append(q,edges[1][1:]) 
    q = np.append(q,edges[2][1:-1]) 
    return q 
    
"""
The `split` represents a partial derivation of j derivative sequences from 
an <OpTri>, a triangular matrix that contains m derivative sequences for a 
vector of length (m+1). The `split` is a 2-tuple of the form
    (P45: 45-degree part,P90: 90-degree part). 
In an upper-right hand triangle (see diagram),

    y x x x x 
    0 y x x x 
    0 0 y x x 
    0 0 0 y x 
    0 0 0 0 y   , 

the 45-degree part is a contiguous subsequence belonging to the diagonal 
(denoted y) starting at index (0,0) [uppermost leftist]. The 90-degree part 
is the i'th subrow, i being |P45|-1. 

The number of derivative sequences, otherwise known as the degree, from the 
`split` is j := |P45| <= m. 

The 1st derivative sequence is of the greatest length: |P45| + |P90|. 
For a j'th derivative sequence, the length is |P45| + |P90| - (j - 1). 

Derivative sequences are calculated starting from the |P45|'th one. 

                     
y ? ? ? ? ? ? ? ?  [derivative sequence 1]
  y ? ? ? ? ? ? ?  [derivative sequence 2]
    y ? ? ? ? ? ?  [derivative sequence 3]
      y x x x x x  [derivative sequence 4] 

*   Illustration of a 45-90 split for an <OpTri>; `y` is of P45 and `x` is 
    of P90. 

The |P45|'th derivative sequence is simply the (|P45| - 1)'th row of the <OpTri>, 
in other words,
    <P45[-1],P90[0],P90[1],...,P90[-1]>. 

For all sequences S' of arbitrary derivative order j, the derivative sequence is 
<x_i : x_i := P45[j-1] if i == 0 else x_i := opfunc(x_(i-1), S[i-1]))>, 
i is the index of `x_i` in S' and sequence S is of derivative order (j+1). 

The pairwise `opfunc` is typically set to the addition operation. 
"""
class OpTri45N90Split: 

    def __init__(self,split,opfunc=add):
        assert type(split) is tuple and len(split) == 2 
        assert len(split[0]) > 0 
        self.split = split 
        self.opfunc = add 
        self.derivative_seqs = dict() 
        self.m_ = None 

    def degree(self): 
        return len(self.split[0])

    def size(self): 
        r = len(self.split[1]) + 1
        s = r
        for i in range(1,len(self.split[0])): 
            s += r + i 
        return s 

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
    def to_triangular_matrix(self): 
        d = self.degree() 
        self.m_ = np.zeros((d,d),dtype=np.int32)
        for i in range(1,d+1): 
            if i not in self.derivative_seqs:
                self.j_derivative_seq(i)
            self.m_[i-1,i-1:] = np.array(\
                self.derivative_seqs[i],dtype=np.int32) 
        return deepcopy(self.m_)

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
"""
Given an <OpTri> upper-right hand triangular matrix `m`, 
flip `m` by the 0-axis or the 1-axis. 

`m` :=
    r x x x 
    0 s x x 
    0 0 t x 
    0 0 0 u , 

`m` w/ 0-axis flip := 

    0 0 0 u
    0 0 t x 
    0 s x x 
    r x x x , 

`m` w/ 1-axis flip := 

    x x x r
    x x s 0
    x t 0 0 
    u 0 0 0 . 

    ** Brief ** 
The flipped <OpTri> is `mf`. And the sequences for the derivative of 
`mf` is first calculated by iterating through the rows of `mf`, starting 
with index 0. If index is 0, the starting value for the pairwise opfunc 
operation is `intseed`. Otherwise, the starting value is the last element 
of the calculated sequence pertaining to the previous row. 

The sequences are either reversed w.r.t. themselves (if 1-axis flip) or 
reversed w.r.t. their index ordering amongst themselves. The result is 
an upper-right hand triangular matrix `m_`, the flip-derivative of the 
original matrix `m`. 
"""
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
        self.m2_ = None 
        self.intseed = intseed 
        self.opfunc = opfunc 
        self.axis = axis 

    def size(self): 
        r = self.m_.shape[0] 
        return np.sum([_ for _ in range(1,r+1)])
    
    def reset_axis(self,axis):
        assert axis in {0,1} 
        self.mf = np.flip(self.mf,self.axis)
        self.axis = axis 
        self.mf = np.flip(self.mf,self.axis)

    def construct_(self,i,intseed):  
        assert i >= 0 and i < self.mf.shape[0]

        if self.axis == 0: 
            subrow = self.mf[i,-(i+1):]
        else: 
            subrow = self.mf[i,:(self.mf.shape[0] -i)] 

        q = [intseed] 

        for v in subrow: 
            q.append(self.opfunc(q[-1],v)) 
        q = q[1:]

        if self.axis == 0: 
            i2 = self.mf.shape[0] - 1 - i 
        else:
            q = q[::-1]  
            i2 = i

        self.m_[i2,i2:] = np.array(q)             
        return q 

    def construct(self):
        self.m2_ = deepcopy(self.m_) 
        intseed = self.intseed 

        j = 0 if self.axis == 1 else -1 
        for i in range(self.mf.shape[0]): 
            q = self.construct_(i,intseed)
            intseed = q[j] 

    def source_seq(self,func): 
        q = np.array([self.intseed],dtype=np.int32)
        for r in self.m_[0]: 
            q = np.append(q, func(r,q[-1])) 
        return q 

    def revert(self): 
        axis = (self.axis + 1) % 2 
        self.m_,self.m2_ = self.m2_,None 
        self.reset_axis(axis)


    def reproduce(self,backward_func):
        source = self.source_seq(backward_func) 
        seq = IntSeq(source) 
        ot = seq.optri(backward_func,np.int32)

        otfd = OpTriFlipDerivation(ot,self.intseed,\
            self.opfunc,self.axis)
        return otfd 
     
# TODO: test 
class OpTriGen:

    def __init__(self,intseed,m,prg,gentype:int,\
        forward_func=add,backward_func=sub,add_noise:bool=False):
        assert type(intseed) in {int,np.int32}
        assert type(m) == np.ndarray 
        assert m.ndim == 2 
        assert m.shape[0] == m.shape[1]
        assert gentype in {1,2} 

        self.intseed = np.int32(intseed) 
        self.m = m 
        self.prg = prg 
        self.gentype = gentype 
        self.ffunc = forward_func
        self.bfunc = backward_func
        self.cache = [] 

        # variables for gentype #1 
        self.ots = None 
        self.ots_available_rows = None 

        # variables for gentype #2 
        self.otfd = None 
        self.otfd_available_rows = None 
        self.otfd_prev_axis = None 

        # is a new triangular matrix set for drawing values from? 
        self.reloaded = False 

        # variables used for drawing values in cases of `prg` == None. 
            # variables for gentype #1. 
        self.base_dim = m.shape[0] 
            # variables for gentype #2. 
        self.binary_alternator = 0 

    def __next__(self): 
        if len(self.cache) == 0: 
            self.load_by_gentype()
        else: 
            self.reloaded = False  
        q = self.cache.pop(0)
        return q 

    def batch_size(self): 
        if self.gentype == 1: 
            return self.ots.size() 
        return self.otfd.size() 

    def load_by_gentype(self): 
        if self.gentype == 1: 
            self.load_from_jagged45N90() 
            return 
        self.load_from_flipderivation() 

    def load_data_struct(self,datastruct):
        if type(datastruct) == OpTri45N90Split:
            self.ots = datastruct 
            self.ots_available_rows = [i for i in range(1,len(self.ots.split[0]) + 1)] 
        elif type(datastruct) == OpTriFlipDerivation:
            self.otfd = datastruct 
            self.otfd.construct() 
            self.otfd_available_rows = [i for i in range(self.otfd.m_.shape[0])] 
            self.otfd_prev_axis = [self.otfd.axis]  
        else:
            assert False 
        

    #-------------------- generator type #1: Jagged 45-90 Split 
    # TODO: test this section. 

    def jagged_split_(self,row_indices,split_index): 
        l = []
        for (i,r) in enumerate(row_indices): 
            assert r <= i 
            l.append(self.m[r,i]) 
        
        p45,p90 = l[:split_index],l[split_index:]
        return (p45,p90)

    def jagged_split_prng(self):

        row_indices = [] 
        for i in range(self.m.shape[0]): 
            j = self.prg() % (i +1) 
            row_indices.append(j)

        split_index = self.prg() % len(row_indices) +1 
        return self.jagged_split_(row_indices,split_index) 

    def set_jagged_split(self): 
        # case: PRNG-less 
        if type(self.prg) == type(None): 
            js = (np.diag(self.m),[])
        else: 
            js = self.jagged_split_prng()
        ots = OpTri45N90Split(js,self.ffunc)
        self.load_data_struct(ots)

    def store_ots_row_(self,i,store_seq=True): 
        q = self.ots.j_derivative_seq(i,store_seq)
        self.cache.extend(q) 
        return

    def store_ots_row(self,store_seq=True): 
        if len(self.ots_available_rows) == 0: 
            if type(self.prg) == type(None): 
                self.new_default_45N90_split()
            else: 
                self.reproduce_OpTri_jagged45N90()
                self.reloaded = True 
        else: 
            self.reloaded = False 

        if type(self.prg) == type(None): 
            i = 0 
        else: 
            i = self.prg() % len(self.ots_available_rows)
        s = self.ots_available_rows.pop(i)
        self.store_ots_row_(s,store_seq)  

    def reproduce_OpTri_jagged45N90(self): 
        p45_seqsize = (self.prg() % self.ots.degree()) + 1 
        p90_seqsize = (self.prg() % self.ots.degree()) + 1 

        p45_indices,p90_indices = [],[] 

        for _ in range(p45_seqsize): 
            s = (self.prg() % self.ots.degree()) + 1 
            p45_indices.append(s)
        for _ in range(p90_seqsize): 
            s = (self.prg() % self.ots.degree()) + 1 
            p90_indices.append(s)
        ots = self.ots.reproduce(p45_indices,p90_indices)
        self.load_data_struct(ots)

    def load_from_jagged45N90(self): 
        stat = None 
        if type(self.ots_available_rows) == type(None): 
            self.set_jagged_split() 
            stat = True 
        self.store_ots_row()

        if stat: 
            self.reloaded = True 
        return

    #------------------------- generator type #2: Flip-Derivation 
    # TODO: test this section. 

    def set_otfd(self): 
        if type(self.prg) == type(None): 
            axis = 0 
        else: 
            axis = self.prg() % 2 

        otfd = OpTriFlipDerivation(self.m,self.intseed,self.ffunc,axis)
        self.load_data_struct(otfd)

    def store_otfd_row(self): 

        if len(self.otfd_available_rows) == 0: 
            # case: PRNG-less 
            if type(self.prg) == type(None): 
                self.new_default_FlipDerivation()
            # case: flip onto another axis
            elif len(self.otfd_prev_axis) < 2 and \
                self.prg() % 2 == 1: 
                self.otherflip_OpTriFlipDerivation()
            # case: reproduce from previous <OpTriFlipDerivation>
            else: 
                self.reproduce_OpTriFlipDerivation()
            self.reloaded = True 
        else: 
            self.reloaded = False 

        if type(self.prg) == type(None): 
            i = 0 
        else: 
            i = self.prg() % len(self.otfd_available_rows)
        
        s = self.otfd_available_rows.pop(i)
        seq = self.otfd.m_[s,s:]
        self.cache.extend(seq) 

    def otherflip_OpTriFlipDerivation(self): 
        na = (self.otfd_prev_axis[-1] + 1) % 2
        self.otfd_prev_axis.append(na) 
        self.otfd.reset_axis(na)
        self.otfd.construct() 
        self.otfd_available_rows = [i for i in range(self.otfd.m_.shape[0])] 

    def reproduce_OpTriFlipDerivation(self): 
        self.otfd = self.otfd.reproduce(self.bfunc)
        axis = self.prg() % 2 if type(self.prg) != type(None) else 0 
        if self.otfd.axis != axis: 
            self.otfd.reset_axis(axis)  
        self.otfd.construct() 
        self.otfd_available_rows = [i for i in range(self.otfd.m_.shape[0])] 
        self.otfd_prev_axis = [self.otfd.axis]
        
    def load_from_flipderivation(self): 
        stat = None 
        if type(self.otfd_available_rows) == type(None): 
            self.set_otfd() 
            stat = True 

        self.store_otfd_row()
        if stat: 
            self.reloaded = True 
        return

    #------------------------- default methods for PRNG-less mode 

    def new_default_45N90_split(self):
        new_seq = None 
        self.ots.to_triangular_matrix()
        d = self.ots.degree() 

        if d == self.base_dim: 
            new_seq = triangular_matrix_circumference(self.ots.m_,\
                start_edge=0,is_clockwise=True)            
        else: 
            new_seq = triangular_matrix_circumference(self.ots.m_,\
                start_edge=0,is_clockwise=False)
            new_seq = new_seq[:self.base_dim]

        split = (list(new_seq),[]) 
        ots = OpTri45N90Split(split,opfunc=self.ffunc)
        self.load_data_struct(ots)
        self.reloaded = True 
        return 

    def new_default_FlipDerivation(self): 
        
        if len(self.otfd_prev_axis) < 2: 
            self.otherflip_OpTriFlipDerivation()
            self.reloaded = True 
            return 

        if self.binary_alternator == 0:
            self.otfd.revert() 
        self.binary_alternator = (self.binary_alternator + 1) % 2
        self.reproduce_OpTriFlipDerivation()
        self.reloaded = True 