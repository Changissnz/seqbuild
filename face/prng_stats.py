"""
data structures used to simplify QUALTEST and CHAINTEST b/t and w/in generators. 
"""
from .prng_bool_cmp import * 

DEFAULT_QUALTEST_SEGMENT_LENGTH = 200 
DEFAULT_CHAINTEST_SEGMENT_LENGTH = 1000 

"""
Method<run> conducts QUALTEST and CHAINTEST for their respective 
iteration numbers a.k.a segments (see the default variables in this file). 
The segments are split into `partition_size` chunks, and `ngram_info` 
specifies the n-gram length/s for the measurement processes. 
"""
class PRNGStats: 

    def __init__(self,prg,ngram_info,partition_size): 
        assert type(prg) in {MethodType,FunctionType} 
        assert type(partition_size) == int and partition_size > 1 

        if type(ngram_info) != int: 
            assert type(ngram_info) == list or is_vector(ngram_info) 
            for x in ngram_info: 
                assert type(x) in {int,np.int32,np.int64} 
                assert x > 1 
            self.deg_vec = ngram_info
            self.gauge_depth = None 
        else: 
            assert ngram_info > 1 
            self.gauge_depth = ngram_info 
            self.deg_vec = None 

        self.prg = prg
        self.psize = partition_size

        self.cqueue = deque() 
        self.qqueue = deque() 

        self.qresults = [] 
        self.qdiff = [] 

        self.cresults = [] 
        self.cdiff = [] 
        return  

    def run(self):
        self.qtest() 
        self.ctest() 
        print("============================================================")
        return

    def latest_results(self):
        if len(self.qresults) == 0: return None,None 

        return self.qresults[-1],self.cresults[-1] 

    def qtest(self): 
        X = self.qtest_() 

        self.qresults.append(X) 
        
        if len(self.qresults) > 1: 
            qdiff = cmp_generators__MultiMetric_(X,self.qresults[-1]) 
            self.qdiff.append(qdiff) 
        return 

    def qtest_(self): 
        self.populate_queue(False) 
        I = [self.qqueue.popleft() for _ in range(DEFAULT_QUALTEST_SEGMENT_LENGTH)] 
        G = prg__iterable(I) 

        return gauge_generator__MultiMetric(G,DEFAULT_QUALTEST_SEGMENT_LENGTH,self.gauge_depth,\
            self.deg_vec,set_frange=True,condense_ngram_output=True,full_output=True)

    def ctest(self): 
        X = self.ctest_() 

        self.cresults.append(X) 
        
        if len(self.cresults) > 1: 
            cdiff = X - self.cresults[-1]
            self.cdiff.append(cdiff) 
        return

    def ctest_(self): 
        self.populate_queue(True) 
        I = [self.cqueue.popleft() for _ in range(DEFAULT_CHAINTEST_SEGMENT_LENGTH)] 
        G = prg__iterable(I) 
        ag2 = APRNGGaugeV2(G,(0.,1.),0.5) 
        return ag2.chaintest(DEFAULT_CHAINTEST_SEGMENT_LENGTH,self.psize) 
        

    def populate_queue(self,is_cqueue:bool):
 
        if is_cqueue:
            Q = self.cqueue     
            Q2 = self.qqueue 
            s = DEFAULT_CHAINTEST_SEGMENT_LENGTH
        else:
            Q = self.qqueue
            Q2 = self.cqueue
            s = DEFAULT_QUALTEST_SEGMENT_LENGTH

        print("populating {}->{}".format(len(Q),s))

        while len(Q) < s: 
            x = self.prg() 
            Q.append(x) 
            Q2.append(x) 
        print("done populating")

'''
A comparator class for input PRNGs, of `input_prg_list`, when they are each used as a 
PRNG variable in an output (secondary) PRNG, `output_prg_struct`. The `output_prg_struct` 
is a Python class, such that its `__next__` method is a PRNG function. 

This algorithm starts off by creating |`input_prg_list`| versions of the `output_prg_struct`, 
each i'th one fitted its `input_prg_varname` the i'th input PRNG. 

Class is built on top of class<PRNGStats>, serving as a comparator of QUALTEST and CHAINTEST 
scores. These are the four primary comparisons: 
- pairwise comparison between input PRNGs. 
- comparison between input PRNG and its corresponding output PRNG. 
- pairwise comparison between output PRNGs. 
- pairwise comparison between two (input,output) PRNG differences. 

If the method variable<metric_numbers> is non-null (a set of integers specifying the 
boolean metrics from file<face.prng_bool_cmp> to use), these comparisons output 
1 (more random)|0 (equally random)|-1 (less random).
'''
class PRNGIOStats: 

    def __init__(self,input_prg_list,input_prg_varname,output_prg_struct,ngram_info,\
        partition_size,absdiff:bool=True,output_prg_range=None,collect_variance_over_runs:bool=False): 
        
        assert len(input_prg_list) > 0 
        
        for x in input_prg_list: 
            assert type(x) in {MethodType,FunctionType}
        assert type(absdiff) == bool 

        if type(output_prg_range) != type(None): 
            assert is_valid_range(output_prg_range,True,False) or \
            is_valid_range(output_prg_range,False,False)

        assert type(collect_variance_over_runs) == bool 

        self.input_prg_list = input_prg_list
        self.input_prg_varname = input_prg_varname
        self.output_prg_struct = output_prg_struct
        self.cvar_over_runs = collect_variance_over_runs
        self.ngram_info = ngram_info 
        self.psize = partition_size
        self.absdiff = absdiff 

        self.input_prng_stats = [] 
        self.output_prng_stats = []  
        self.output_prg_range = output_prg_range

        # index -> score term -> list of values over runs 
        self.input_prng_summary = defaultdict(None) 
        self.output_prng_summary = defaultdict(None) 
        self.preproc() 

        return  

    def preproc(self):
        for (i,x) in enumerate(self.input_prg_list):

            G = deepcopy(self.output_prg_struct) 
            setattr(G,self.input_prg_varname,deepcopy(x))

            x_ = x 
            Q = deepcopy(G.__next__)
            if type(self.output_prg_range) != type(None): 
                x_ = wrap_ranged_modulo_over_generator(x_,self.output_prg_range)
                Q = wrap_ranged_modulo_over_generator(Q,self.output_prg_range)

            G0 = PRNGStats(x_,deepcopy(self.ngram_info),self.psize)
            self.input_prng_stats.append(G0)

            G1 = PRNGStats(Q,deepcopy(self.ngram_info),self.psize)
            self.output_prng_stats.append(G1)

            self.input_prng_summary[i] = defaultdict(list)
            self.output_prng_summary[i] = defaultdict(list)

    def update_summary(self,index,is_input_prg:bool): 
        assert type(is_input_prg) == bool 
        assert 0 <= index < len(self.input_prg_list) 

        x = self.input_prng_stats[index] if is_input_prg else \
            self.output_prng_stats[index] 

        q0,q1 = x.latest_results()
        qc = QCInfo((q0,q1)) 

        D = qc.retrieve_all() 
        
        D2 = self.input_prng_summary[index] if is_input_prg else \
            self.output_prng_summary[index] 

        for k,v in D.items(): 
            D2[k].append(v) 

    def run(self): 
        for (i,x) in enumerate(self.input_prng_stats):
            print("* running input prng #{}".format(i))
            x.run() 
            
            if self.cvar_over_runs: self.update_summary(i,True)

        for (i,y) in enumerate(self.output_prng_stats):
            print("* running output prng #{}".format(i))
            y.run()

            if self.cvar_over_runs: self.update_summary(i,False)

    """
    Compares the intra-variance between two PRNGs. 

    index_type_pair0 := (PRNG index, ?is input PRNG?)
    index_type_pair1 := (PRNG index, ?is input PRNG?)
    """
    def cmp_var_over_runs(self,index_type_pair0,index_type_pair1): 
        assert len(index_type_pair0) == len(index_type_pair1) == 2 
 
        G0 = self.input_prng_summary[index_type_pair0[0]] if \
            index_type_pair0[1] else self.output_prng_summary[index_type_pair0[0]]
        G1 = self.input_prng_summary[index_type_pair1[0]] if \
            index_type_pair1[1] else self.output_prng_summary[index_type_pair1[0]]

        var0 = {k:np.var(v) for k,v in G0.items()}
        var1 = {k:np.var(v) for k,v in G1.items()} 

        D = dict()  

        for k,v in var0.items(): 
            r = to_trinary_relation(v,var1[k]) 
            D[k] = r 
        
        print("** variance comparison")
        Q = sorted(SCORE_TERMS) 
        for q in Q: 
            print("{} -> {}".format(q,D[q])) 

        return to_trinary_relation_v2(sum(list(D.values())),None,True,False)

    # TODO: test 
    def cmp_matrix(self,cmp_type:str,metric_numbers:set=None): 
        assert cmp_type in {"i","o","io"} 

        l = len(self.input_prg_list)
        M = np.zeros((l,l)) + np.nan 

        if cmp_type == "i": 
            F = self.input_cmp 
        elif cmp_type == "o": 
            F = self.output_cmp
        else: 
            F = self.pairwise_io_cmp 

        for i in range(l - 1): 
            for j in range(i+1,l): 
                c = F(j,i,metric_numbers) 
                M[i,j] = c 
        return M 

    """
    input PRNG#0 - input PRNG#1 

        XOR [if `metric_numbers` != NULL]

    is input PRNG#0 more random? 
    """
    def input_cmp(self,index0,index1,metric_numbers:set=None):
        G0 = self.input_prng_stats[index0]
        G1 = self.input_prng_stats[index1]

        if type(metric_numbers) == set: 
            return self.boolcmp_stats(G0,G1,metric_numbers)
        
        assert type(metric_numbers) == type(None)
        return self.cmp_stats(G0,G1)

    """
    output PRNG#0 - output PRNG#1 

        XOR [if `metric_numbers` != NULL]

    is output PRNG#0 more random? 
    """
    def output_cmp(self,index0,index1,metric_numbers:set=None):
        G0 = self.output_prng_stats[index0]
        G1 = self.output_prng_stats[index1]

        if type(metric_numbers) == set: 
            return self.boolcmp_stats(G0,G1,metric_numbers)
        
        assert type(metric_numbers) == type(None)
        return self.cmp_stats(G0,G1)
        
    """
    output PRNG - input PRNG 

        XOR [if `metric_numbers` != NULL]

    is output PRNG more random? 
    """
    def io_cmp(self,input_index,metric_numbers:set=None): 
        G0 = self.input_prng_stats[input_index]
        G1 = self.output_prng_stats[input_index]

        if type(metric_numbers) == set: 
            return self.boolcmp_stats(G1,G0,metric_numbers)
        
        assert type(metric_numbers) == type(None)
        return self.cmp_stats(G1,G0)


    """
    (output PRNG#0 - input PRNG#0) - (output PRNG#1 - input PRNG#1)

        XOR [if `metric_numbers` != NULL]

    is (output PRNG#0 - input PRNG#0) more random? 
    """
    def pairwise_io_cmp(self,index0,index1,metric_numbers:set=None):

        q0,c0 = self.io_cmp(index0)
        q1,c1 = self.io_cmp(index1) 

        x0,x1 = q0[0],q1[0] 
        x0 = (x0[:3],x0[3:5],x0[5]) 
        x1 = (x1[:3],x1[3:5],x1[5]) 

        q0 = (x0,q0[1],q0[2]) 
        q1 = (x1,q1[1],q1[2]) 

        if type(metric_numbers) == set: 
            return is_PRNG_more_QC_random((q0,c0),(q1,c1),metric_numbers)

        D0 = cmp_generators__MultiMetric_(q0,q1,self.absdiff)
        D1 = c0 - c1 
    
        return D0,D1 

    def cmp_stats(self,G0,G1):
        assert type(G0) == type(G1) == PRNGStats 

        q0,c0 = G0.latest_results() 
        q1,c1 = G1.latest_results() 

        D0 = cmp_generators__MultiMetric_(q0,q1,self.absdiff)
        D1 = c0 - c1 

        return D0,D1 

    """
    return: 
    1(more random)|-1(less random)|0(equally random)
    """
    def boolcmp_stats(self,G0,G1,metric_numbers:set):
        assert type(G0) == type(G1) == PRNGStats 
        q0,c0 = G0.latest_results() 
        q1,c1 = G1.latest_results() 

        return is_PRNG_more_QC_random((q0,c0),(q1,c1),metric_numbers)

def rank_prng_cmp_matrix(M,verbose=True): 
    assert len(M.shape) == 2 
    assert M.shape[0] == M.shape[1] 
    
    M_ = deepcopy(M)
    M_[np.isnan(M_)] = 0. 

    D = [] 

    for i in range(M_.shape[0]): 
        s0 = np.sum(M_[:,i]) 
        s1 = np.sum(M_[i,:]) 
        D.append((i,s0 - s1)) 

    D = sorted(D,key=lambda x:x[1],reverse = True) 
    I = [d[0] for d in D] 

    for i in range(M_.shape[0]): 
        x = I.index(i) 
        y = D[x][1] 
        if verbose: print("index {} -> {} -> {}".format(i,x,y))

    return D  