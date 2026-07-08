from mini_dm.affine_fit_search import * 
from desi.nvec_gen import * 
from collections import deque 

DEFAULT_AFSGEN_MAX_VECTOR_SIZE_RANGE = [75,150]
DEFAULT_AFSGEN_NUM_CENTROIDS = 5 

DEFAULT_AFSGEN_CENTROID_UPDATE_RANGE = [2,13] 

class AFSGen: 

    def __init__(self,prg,prg2,prg3,prg4,num_vectors):
        assert type(prg) in {MethodType,FunctionType}
        assert type(prg2) in {MethodType,FunctionType,type(None)}
        assert type(prg3) in {MethodType,FunctionType,type(None)}
        assert type(prg4) in {MethodType,FunctionType,type(None)}
        
        assert type(num_vectors) == int and num_vectors > 0 
        #self.num_vectors = modulo_in_range(num_vectors,DEFAULT_AFSGEN_MAX_VECTOR_SIZE_RANGE)
        self.num_vectors = num_vectors

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
            v = list(x.flatten()) 
            v = prg_seqsort(v,self.prg)

            self.queue.extend(v) 
            self.adjust_centroids() 
        return self.queue.popleft() 

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
