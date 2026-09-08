from .prng_stats import * 
from .make_cmd import MAIN_method_for_object

class ChainedPRNGComparator: 

    def __init__(self,prng_struct_seq,prg_varname,ngram_info,partition_size,\
        output_prg_range=None,do_embed_generators:bool=True): 

        if type(output_prg_range) != type(None): 
            assert is_valid_range(output_prg_range,True,False) or is_valid_range(output_prg_range,False,False)
        assert len(prng_struct_seq) > 1 
        assert type(do_embed_generators) == bool 

        self.prng_struct_seq = prng_struct_seq
        self.prg_varname = prg_varname 
        self.ngram_info = ngram_info
        self.psize = partition_size
        self.output_prg_range = output_prg_range 
        self.do_embed = do_embed_generators

        self.pstat_seq = None 
        self.preproc() 
        self.processed = False 
        return 

    def preproc(self): 

        # set every prng struct's parameter `prg_varname` to the previous prng. 
        if self.do_embed: 
            for i in range(1,len(self.prng_struct_seq)): 
                q0 = deepcopy(self.prng_struct_seq[i-1]) 
                prg = MAIN_method_for_object(q0) 
                setattr(self.prng_struct_seq[i],self.prg_varname,prg) 

        prgs = [MAIN_method_for_object(p) for p in self.prng_struct_seq] 

        # apply modulii range if not None 
        if type(self.output_prg_range) != type(None): 
            prgs = [wrap_ranged_modulo_over_generator(p,self.output_prg_range) for p in \
                prgs] 

        # initialize PRNGStats for every generator 
        self.pstat_seq = [PRNGStats(p,self.ngram_info,self.psize) for p in prgs] 

    def run(self): 
        for (i,x) in enumerate(self.pstat_seq):
            print("* running input prng #{}".format(i))
            x.run() 
        self.processed = True 

    # NOTE: duplicate code with <PRNGIOStats.cmp_matrix>. 
    def cmp_matrix(self,metric_numbers:set = {1,2,3,4,5,6,7}): 
        if not self.processed: 
            print("[!] comparator has not processed all PRNGs") 
            return 

        M = np.zeros((len(self.prng_struct_seq),len(self.prng_struct_seq))) + np.nan

        for i in range(len(self.pstat_seq) - 1): 
            for j in range(i+1,len(self.pstat_seq)): 
                b = self.boolcmp_stats(j,i,metric_numbers)
                M[i,j] = b 
        return M 


    # NOTE: duplicate code of method<PRNGIOStats.boolcmp_stats> 
    """
    return: 
    1(more random)|-1(less random)|0(equally random)
    """
    def boolcmp_stats(self,i0,i1,metric_numbers:set):

        G0,G1 = self.pstat_seq[i0],self.pstat_seq[i1] 

        q0,c0 = G0.latest_results() 
        q1,c1 = G1.latest_results() 

        return is_PRNG_more_QC_random((q0,c0),(q1,c1),metric_numbers)

"""
uses a unique copy of `prg` to generate every struct in `struct_name_list` 
(no embedded generators). Then uses a <ChainedPRNGComparator> to produce a 
cumulative comparison matrix over `num_iter` segments. 
"""
def prg_to_chained_PRNG_cmp_matrix(prg,struct_name_list,ngram_info=[8],partition_size=5,\
    output_prg_range=[-50000,50000],num_iter=1,include_source_prg_in_cmp:bool=False): 

    a = int(include_source_prg_in_cmp)
    M = np.zeros((len(struct_name_list) + a,len(struct_name_list) + a))
    S = [one_P_random_struct(l,deepcopy(prg)) for l in struct_name_list]  
    if include_source_prg_in_cmp: S.insert(0,prg) 

    X = ChainedPRNGComparator(S,"prg",ngram_info,partition_size,output_prg_range,do_embed_generators = False) 

    for _ in range(num_iter):  
        X.run() 
        q = X.cmp_matrix() 
        print("{}: {}".format(_,q))
        M = M + q  
    return M 