from desi.multi_metric import * 

class AGV2GuidedGen: 

    def __init__(self,base_prg,aux_prg,\
        base_output_span,density_measure):

        self.base_prg = base_prg 
        self.aux_prg = aux_prg 

        self.bspan = None 
        self.base_seq = [] 
        self.bs_summary = None 
        self.guided_rep = None        
        self.density = None
        self.set_span(base_output_span)
        self.set_density(density) 

        self.output_mode = "base" 
        return 

    def set_span(self,span):
        assert is_number(span,{float,np.float16,np.float32,np.float64}) 
        assert span > 1 
        self.bspan = span 

    def set_density(self,density):
        assert type(density) == np.ndarray
        if is_vector(density):
            d_ = np.zeros((2,2),dtype=int) 
            d_[0],d_[1] = density,density
            density = d_
         
        assert density.shape == (2,2)
        assert_nm(tuple(density[0])) 
        assert_nm(tuple(density[1])) 
        self.density = density
        return

    def next__base(self):
        assert self.output_mode == "base" 
        
        if len(self.base_seq) >= self.bspan:
            return None

        x = self.base_prg()
        self.base_seq.append(x)
        return x 

    def next__guidedrepl(self): 
        assert len(self.base_seq) > 0
        return -1 

    def inspect_base_seq(self,ngram=None):
        if type(ngram) == type(None): 
            ngram = default_ngram_length(len(self.base_seq))  
        mm = MultiMetric(self.base_seq)
        smry = mm.summarize(ngram,False)

        mm.load_mc_map()
        mcm = mm.modcomplex_map 
        self.bs_summary = [smry,mcm]
        return self.bs_summary 