from desi.multi_metric import * 
from desi.seqcov_perm import * 

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

Used to adjust output distribution from <AGV2GuidedGen>. 
There are 2 ways to use this container: 
1. `refvar_catvec` and `refvar_catmap`,
2. `covuwpd_log` and `factormap_log`. 

In the first approach, the method <log_sample_cat> is called 
to record a quality category for a sample (a vector). Then the 
method <refvar_frequency_map> can be used to calculate categorical 
frequencies for the last `previous_iter` iterations. The first 
approach observes the reference variable, set at `refvar`. The 
reference variable is one of ('cov','uwpd',modulus:float). 

In the second approach, the method <catfreq_map> is used to 
calculate the categorical frequencies of any quality 
('cov','uwpd',modulus:float). The calculation draws information 
from `covuwpd_log` and `factormap_log`, maps that are updated 
from calling method<update_one_element>. 

NOTE: 
The measure of density is not fully implemented. Instead, class 
essentially serves as frequency counter. 
"""
class AGV2DensityLog:

    def __init__(self,cat_sz:int=10): 
        self.cat_sz = None 
        self.reload_cat_sz(cat_sz) 

        self.covuwpd_log = [] 
        self.factormap_log = defaultdict(list) 
        self.refvar = None
        self.refvar_catvec = []
        self.refvar_catmap = defaultdict(int) # category -> frequency 
        return 

    def log_sample_cat(self,sample_cat): 
        self.refvar_catvec.append(sample_cat) 
        self.refvar_catmap[sample_cat] += 1 
    
    """
    main method #1.
    """
    def refvar_frequency_map(self,previous_iter:int): 
        assert type(previous_iter) in {int,np.int32,np.int64} 
        assert previous_iter > 0 
        fmap = defaultdict(int)

        for x in self.refvar_catvec[-previous_iter:]:
            fmap[x] += 1 
        return fmap 

    """
    main method #2. 

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
        assert type(density_measure) != type(None) 

        q = None 
        if density_measure == "cov": 
            q = [c[0] for c in self.covuwpd_log]
        elif density_measure == "uwpd": 
            q = [c[1] for c in self.covuwpd_log] 
        else:
            assert density_measure in self.factormap_log
            q = [c[0] for c in self.factormap_log[density_measure]]

        return q 

    def update_one_element(self,cov,uwpd,d): 
        self.update_covuwpd(cov,uwpd) 
        self.update_factormap(d)

    """
    method used to help with handling density measure requirements; 
    method acts as filtration against unwanted output sequences (from 
    <AGV2GuidedGen>).
    """
    def classify_value(self,vx,l,rv):
        if self.refvar == "cov": 
            v0 = std_classify_one_value(vx,1/self.cat_sz,0.0)

        elif self.refvar == "uwpd": 
            max_value = max_float_uwpd(l,rv)
            v0 = std_classify_one_value(vx,max_value/self.cat_sz,0.0)
        else: 
            max_value = max_float_uwpd(l,[0,self.refvar])
            v0 = std_classify_one_value(vx,max_value/self.cat_sz,0.0)
        return v0 

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
        return

    @staticmethod 
    def label_value(v,dx):
        assert is_number(v) 
        seqcat_length = 1.0 / dx 
        return int(round(v/seqcat_length,0))
        
"""
generator is a wrapper around a base generator `base_prg`. 
The <AGV2GuidedGen> is a generator that uses a numerical 
sequence of length `base_output_span` from `base_prg` as 
a template. Then the generator outputs 'sibling' sequences 
from this template sequence. These 'sibling' sequences have 
qualities of coverage,uwpd, and modular uwpd that are prioritized 
by the algorithm due to low or non-existing frequency (0-frequency 
is new). In more words, this generator adjusts base sequences 
from `base_prg` in order for the output sequences to satisfy 
uniform density distribution of some quality ('cov','uwpd', modular 
uwpd). 

Generator should be considered as more of a wrapper than a 
generator in and of itself. The functions of this class are 
designed specifically to modulate base sequences for the output 
of alternative base sequences. 
"""
class AGV2GuidedGen: 

    def __init__(self,base_prg,aux_prg,\
        base_output_span,density_measure,ngram_length=None,\
        is_context_fed:bool=False):

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
        self.is_context_fed = is_context_fed

        self.output_mode = "base" 

        self.init_density_log() 
        self.inspect_base_seq()

        self.permuter = None 
        return 

    #------------------------ setter functions 

    def set_span(self,span):
        assert is_number(span,{float,np.float16,np.float32,np.float64}) 
        assert span > 1 
        self.bspan = span 

    def set_refvar(self,refvar):
        self.agd_log.reload_refvar(refvar) 

    def set_density(self,density):
        assert type(density) in {int,np.int32,np.int64} 
        assert density > 0 
        self.density = density 
        return 

    def set_permuter(self,epsilon=None): 
        if type(epsilon) != type(None):
            assert type(epsilon) in {float,np.float64,np.float32}

        q = self.agd_log.refvar
        super_range = [min(self.base_seq),max(self.base_seq)]

        if type(epsilon) == type(None):
            delta = (self.aux_prg() % 10000.) / 10000.
        else: 
            delta = epsilon 

        if q == "cov":
            max_radius = (super_range[1] - super_range[0]) / len(self.base_seq) 
            p = SeqCoveragePermuter(np.array(self.base_seq),delta,max_radius / 100,super_range,self.aux_prg)
            p.set_partition(9)
        elif q == "uwpd": 
            mfpd = max_float_uwpd(len(self.base_seq),super_range)
            mfpd = round(mfpd * delta,5) 

            p = SeqUWPDPermuter(np.array(self.base_seq),delta,super_range,self.aux_prg)
        else: 
            super_range = [0,q] 
            p = SeqUWPDPermuter(np.array(self.base_seq),delta,super_range,self.aux_prg,q)
        self.permuter = p 

    def init_density_log(self):
        self.agd_log = AGV2DensityLog() 
        return

    #---------------------------- next functions 

    def next__base(self):
        assert self.output_mode == "base" 

        self.base_seq.clear() 
        for _ in range(self.bspan): 
            x = self.base_prg()
            self.base_seq.append(x)
        return x 

    def inspect_base_seq(self):

        self.next__base()
        q = self.seq_summary(self.base_seq,True)  
        self.bs_summary = q 
        return self.bs_summary 

    def unguided_iter(self,num_vectors:int): 
        assert type(num_vectors) in {int,np.int32,np.int64} 
        assert num_vectors > 0 

        for i in range(num_vectors): 
            self.next__base()
            self.log_seq(self.base_seq)

    def next__guidedrepl(self,update_base_seq:bool=True,num_attempts:int=10):
        if update_base_seq:  
            q = self.base_seq 
            self.next__base() 

        if not self.is_context_fed: 
            qs,cl = self.highest_scoring_permutation(num_attempts)
        else: 
            qs,cl = self.context_fed_permutation() 
        self.agd_log.log_sample_cat(cl)
        return qs 

    def highest_scoring_permutation(self,num_attempts=10):
        v,cl = None,None 
        score = float('inf')
        fmap = self.agd_log.refvar_frequency_map(self.density)

        while num_attempts > 0: 
            self.set_permuter() 
            qs = self.permuter.apply()
            cl = self.classify_seq(qs) 

            if cl in fmap:
                if fmap[cl] < score: 
                    v = qs 
                    score = fmap[cl] 
            else:
                return qs,cl 
            num_attempts -= 1         
        return v,cl 

    def context_fed_permutation(self): 
        # choose a category 
        fmap = self.agd_log.refvar_frequency_map(self.density)

        csize = self.agd_log.cat_sz
        qx = set([_ for _ in range(csize)]) - set(fmap.keys())

        # case: all categories exist 
        if len(qx) == 0:
            return self.highest_scoring_permutation(num_attempts=10) 

        # case: choose a category 
        qx = sorted(qx) 
        i = int(self.base_prg()) % len(qx) 
        cat = qx[i] 

            # convert it to percentage 
        s0 = 1 / csize * cat
        s1 = s0 + 1 / csize 
        s0 = (s0 + s1) / 2.0 

        bstat = int(self.aux_prg()) % 2 
        if bstat: 
            s0 = s0 * -1 

        self.set_permuter(epsilon=s0)
        qs = self.permuter.apply()
        cl = self.classify_seq(qs)

        return qs,cl 

    #------------------------------- sequence summarization 

    def classify_seq(self,seq): 
        # classify value 
        qscore = self.quality_score_for_sequence(seq) 

        rv = None
        if self.agd_log.refvar == "cov": 
            rv = (0.,1.)
        elif self.agd_log.refvar == "uwpd":
            rv = (min(seq),max(seq))
        else: 
            rv = (0.,self.agd_log.refvar)

        l = len(seq)
        return self.agd_log.classify_value(qscore,l,rv)

    def quality_score_for_sequence(self,seq):
        qs = self.seq_summary(seq,False) 

        if self.agd_log.refvar in {"cov","uwpd"}:
            summry = qs[0] 
            i = 0 if self.agd_log.refvar == "cov" else 1 
            qx = np.mean(summry[i]) 
            return qx

        return self.factor_kcomplexity(seq)[self.agd_log.refvar][0]


    #------------------------------- log sequence into memory 

    def seq_summary(self,seq,full_output:bool=False):
        mm = MultiMetric(seq)
        smry = mm.summarize(self.ngram_length,False)

        if not full_output:
            return [smry,None] 

        mm.load_mc_map()
        mcm = mm.modcomplex_map 
        return [smry,mcm]

    def factor_kcomplexity(self,seq):
        qmcm = sorted(list(self.bs_summary[1].keys()))
        return factorseq_to_uwpdcov_map(self.base_seq,qmcm)

    def log_seq(self,seq): 
        summry = self.seq_summary(seq,False)[0] 
        fmap = self.factor_kcomplexity(seq)

        cov = np.mean(summry[0])
        uwpd_ = np.mean(summry[1]) 

        self.agd_log.update_one_element(cov,uwpd_,fmap)