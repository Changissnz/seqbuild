from mini_dm.affine_fit_search import * 
from desi.nvec_gen import * 
from collections import deque 

DEFAULT_AFSGEN_NUM_CENTROIDS = 5 

DEFAULT_AFSGEN_CENTROID_UPDATE_RANGE = [2,13] 

"""
(A)ffine (F)it (S)earch (Gen)erator.

Built on top of class<AffineFitSearch__TypeN2NA>. Starts by initializing an 
<AffineFitSearch__TypeN2NA> instance on some q I/O vector pairs, generated 
by <PointSetGen__TypeAffine>. 

Output values are drawn from `queue`. When `queue` is empty, algorithm populates 
`queue` by first outputting one value v from `prg`. Value v is then used as input 
into each of the <AffineFitSearch__TypeN2NA>'s top <AffineDelta> solutions. The 
resulting difference is an (n x m) matrix M, flattened into a 1-dimensional vector 
and appended to the queue. 

If `diff_mode` is set to True, `prg` is used to generate an (n x m) matrix M_2, and 
the flattened (M - M_2) vector is instead appended to the queue.

prg := PRNG, single output 
prg2 := PRNG, range output XOR None  
prg3 := PRNG, range output XOR None  
prg4 := PRNG, range output XOR None  
num_vectors := x, positive integer specifying the number of I/O pairs the <AffineFitSearch__TypeN2NA> 
                trains over. 
diff_mode := ~ 
"""
class AFSGen: 

    def __init__(self,prg,prg2,prg3,prg4,num_vectors,diff_mode:bool=True):
        assert type(prg) in {MethodType,FunctionType}
        assert type(prg2) in {MethodType,FunctionType,type(None)}
        assert type(prg3) in {MethodType,FunctionType,type(None)}
        assert type(prg4) in {MethodType,FunctionType,type(None)}
        
        assert type(num_vectors) == int and num_vectors > 0 
        assert type(diff_mode) == bool 

        self.num_vectors = num_vectors
        self.diff_mode = diff_mode 

        self.prg = prg 
        self.prg2 = prg2 
        self.prg3 = prg3 
        self.prg4 = prg4 

        if type(self.prg2) != type(None): self.prg2 = prg__single_to_range_outputter(self.prg2) 
        if type(self.prg3) != type(None): self.prg3 = prg__single_to_range_outputter(self.prg3) 
        if type(self.prg4) != type(None): self.prg4 = prg__single_to_range_outputter(self.prg4) 

        self.afsearch = None 
        self.instantiate_AFSearch()

        self.queue = deque() 
        return 

    def __next__(self): 

        if len(self.queue) == 0: 
            q = self.prg() 
            x = self.afsearch.map_input(q)
            x = self.output_difference(x) 
            v = list(x.flatten()) 
            v = prg_seqsort(v,self.prg)

            self.queue.extend(v) 
            self.adjust_centroids() 
        return self.queue.popleft() 

    def output_difference(self,q):
        if not self.diff_mode: return q 

        d = q.shape[0] * q.shape[1] 
        L = np.array([self.prg() for _ in range(d)]) 
        L = np.reshape(L,q.shape)
        return q - L 

    def adjust_centroids(self): 
        x = modulo_in_range(int(self.prg()),DEFAULT_AFSGEN_CENTROID_UPDATE_RANGE) 
        for _ in range(x): 
            next(self.afsearch)
            
    def instantiate_AFSearch(self): 

        psg = PointSetGen__TypeAffine(self.num_vectors,self.prg,self.prg2,self.prg3,self.prg4)
        psg.generate_points(False,True) 

        x,y = np.array(psg.input_seq),np.array(psg.point_seq) 
        k = DEFAULT_AFSGEN_NUM_CENTROIDS

        ma_dim = [y.shape[1],y.shape[1]]

        if prg_decimal(self.prg,[0.,1.]) > 0.5: 
            ma_dim[0] = 0 
        
        if ma_dim[0] != 0: 
            if prg_decimal(self.prg,[0.,1.]) < 0.5:
                ma_dim[1] = 0 

        ma_order = 0 if prg_decimal(self.prg,[0.,1.]) > 0.5 else 1 
        initial_cv = "random" 

        is_bfs_queue = prg_decimal(self.prg,[0.,1.]) > 0.5 

        score_improvement_type = "stochastic" 

        self.afsearch = AffineFitSearch__TypeN2NA(x,y,self.prg,k,ma_dim,\
            ma_order,initial_cv,is_bfs_queue=is_bfs_queue,\
            score_improvement_type = score_improvement_type,\
            score_improvement_reference = "local")

        for _ in range(75): next(self.afsearch) 
        self.afsearch.score_imp_type = "random"
