from desi.multi_metric import * 

def mod_uwpd_of_sequence(S,m):
    assert m != 0 
    S_ = np.array([s_ % m for s_ in S]) 
    return uwpd(S_,pairwise_op=lambda x1,x2: np.abs(x2 - x1),\
        accum_op=lambda x1,x2: x1 + x2)

"""
map with 
key : factor : float 
value : (uwpd of modulated sequence,coverage of modulated sequence) 
"""
def factorseq_to_uwpdcov_map(seq,fseq):
    uwpdcov = dict()
    for m in fseq: 
        pdseq = mod_uwpd_of_sequence(seq,m)
        mseq = np.array([s_ % m for s_ in seq]) 
        covseq = coverage_of_sequence(mseq,[0,m],max_radius=0.5)
        uwpdcov[m] = (pdseq,covseq)
    return uwpdcov

"""
container to hold 
    (coverage,uwpd,%f_i) 
scores over the course of some j iterations. 

Used to calculate density measures 
"""
class AGV2DensityLog:

    def __init__(self,cat_sz:int=10): 
        self.cat_sz = None 
        self.reload_cat_sz(cat_sz) 

        self.measures = dict() 

        # density measures 
            # coverage 
        self.dmap0 = defaultdict(int) 
            # uwpd
        self.dmap1 = defaultdict(int) 
            # specific factor `refvar`
        self.dmap2 = defaultdict(int) 

        self.covuwpd_log = [] 
        self.factormap_log = defaultdict(list) 
        self.refvar = None
        return 

    """
    calculates a map with 
    key: category 
    value: frequency of category. 

    There are `cat_sz` number of categories.
    """
    def catfreq_map(self,density_measure,lx = None):
        q = self.fetch_vec(density_measure)

        clog = None 
        if type(lx) != type(None): 
            assert type(lx) in {int,np.int32,np.int64} 
            clog = q[-lx:]
        else: 
            clog = q 

        scat_vec = stdcat_vec(np.array(clog),self.cat_sz,0.0) 
        fmap = vec_to_frequency_map(scat_vec) 
        return fmap 

    def fetch_vec(self,density_measure):
        assert type(density_measure) ! = type(None) 

        q = None 
        if density_measure == "cov": 
            q = [c[0] for c in self.covuwpd_log]
        elif density_measure == "uwpd": 
            q = [c[1] for c in self.covuwpd_log] 
        else:
            assert density_measure in self.factormap_log
            q = [c[0] for c in self.factormap_log[density_measure]]

        return q 

    def update_covuwpd(self,cov,uwpd): 
        self.covuwpd_log.append((cov,uwpd)) 

    def update_factormap(self,d):
        for k,v in d.items():
            self.factormap_log[k].append(v)

    def reload_cat_sz(self,cat_sz):
        assert type(cat_sz) in {int,np.int32,np.int64}
        assert cat_sz > 0 
        self.cat_sz = cat_sz 
        return

    def reload_refvar(self,refvar): 
        if refvar == "cov":
            self.refvar = refvar 
        elif refvar == "uwpd": 
            self.refvar = refvar
        else:
            assert is_number(refvar)
            self.refvar = refvar

        self.dmap0.clear() 
        self.dmap1.clear() 
        self.dmap2.clear() 
        return

    @staticmethod 
    def label_value(v,dx):
        assert is_number(v) 
        seqcat_length = 1.0 / dx 
        return int(round(v/seqcat_length,0))

    def density_count(self):
        return -1 

    def update_density_count(self,v): 
        return -1 

    def label_sample_vector(self,seq):
        return -1 
        
class AGV2GuidedGen: 

    def __init__(self,base_prg,aux_prg,\
        base_output_span,density_measure,ngram_length=None):

        self.base_prg = base_prg 
        self.aux_prg = aux_prg 

        self.bspan = None 
        self.base_seq = [] 
        self.output_queue = [] 

        self.bs_summary = None 
        self.guided_rep = None        
        self.density = None
        self.set_span(base_output_span)
        self.set_density(density_measure) 

        self.ngram_length = ngram_length
        self.output_mode = "base" 

        self.init_density_log() 
        self.inspect_base_seq()
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

        self.base_seq.clear() 
        for _ in range(self.bspan): 
            x = self.base_prg()
            self.base_seq.append(x)
        return x 

    def next__guidedrepl(self): 
        assert len(self.base_seq) > 0
        return -1 

    def inspect_base_seq(self):

        self.next__base()
        q = self.seq_summary(self.base_seq,True)  
        self.bs_summary = q 
        return self.bs_summary 

    def init_density_log(self):
        self.agd_log = AGV2DensityLog() 
        return

    def seq_summary(self,seq,full_output:bool=False):
        mm = MultiMetric(seq)
        smry = mm.summarize(self.ngram_length,False)

        if not full_output:
            return [smry,None] 

        mm.load_mc_map()
        mcm = mm.modcomplex_map 
        return [smry,mcm]

    def reload_mcm(self):
        return -1 

    def unguided_iter(self,num_vectors:int): 
        assert type(num_vectors) in {int,np.int32,np.int64} 
        assert num_vectors > 0 

        for i in range(num_vectors): 
            self.next__base()
            self.log_seq(self.base_seq)

    def factor_kcomplexity(self,seq):
        qmcm = sorted(list(self.bs_summary[1].keys()))
        return factorseq_to_uwpdcov_map(self.base_seq,qmcm)
    
    def log_seq(self,seq): 
        summry = self.seq_summary(seq,False)[0] 
        fmap = self.factor_kcomplexity(seq)

        cov = np.mean(summry[0])
        uwpd_ = np.mean(summry[1]) 

        self.agd_log.update_covuwpd(cov,uwpd_) 
        self.agd_log.update_factormap(fmap)