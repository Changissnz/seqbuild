"""
data structures used to simplify QUALTEST and CHAINTEST b/t and w/in generators. 
"""

from desi.multi_metric import * 
from collections import deque 
from intigers.mod_prng import prg__iterable 

DEFAULT_QUALTEST_SEGMENT_LENGTH = 200 
DEFAULT_CHAINTEST_SEGMENT_LENGTH = 1000 
DEFAULT_CHAINTEST_RADIUS = 0.5 

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

        while len(Q) < s: 
            x = self.prg() 
            Q.append(x) 
            Q2.append(x) 

class PRNGIOStats: 

    '''
    ''' 
    def __init__(self,input_prg_list,input_prg_varname,output_prg_struct,ngram_info,\
        partition_size): 
        
        assert len(input_prg_list) > 0 
        
        for x in input_prg_list: 
            assert type(x) in {MethodType,FunctionType}

        self.input_prg_list = input_prg_list
        self.input_prg_varname = input_prg_varname
        self.output_prg_struct = output_prg_struct
        self.ngram_info = ngram_info 
        self.psize = partition_size

        self.input_prng_stats = [] 
        self.output_prng_stats = []  
        self.preproc() 
        return  

    def preproc(self):
        for x in self.input_prg_list:
            G0 = PRNGStats(x,deepcopy(self.ngram_info),self.psize)
            self.input_prng_stats.append(G0)

            G = deepcopy(self.output_prg_struct) 
            setattr(G,self.input_prg_varname,deepcopy(x))
            G1 = PRNGStats(G.__next__,deepcopy(self.ngram_info),self.psize)
            self.output_prng_stats.append(G1) 

    def run(self): 
        for (i,x) in enumerate(self.input_prng_stats):
            print("* running input prng #{}".format(i))
            x.run() 

        for (i,y) in enumerate(self.output_prng_stats):
            print("* running output prng #{}".format(i))
            y.run() 

    def input_cmp(self,index0,index1):
        G0 = self.input_prng_stats[index0]
        G1 = self.input_prng_stats[index1]

        return self.cmp_stats(G0,G1)
        
    def io_cmp(self,input_index): 
        G0 = self.input_prng_stats[input_index]
        G1 = self.output_prng_stats[input_index]
        return self.cmp_stats(G1,G0) 

    def pairwise_io_cmp(self,index0,index1):

        q0,c0 = self.io_cmp(index0)
        q1,c1 = self.io_cmp(index1) 

        x0,x1 = q0[0],q1[0] 
        x0 = (x0[:3],x0[3:5],x0[5]) 
        x1 = (x1[:3],x1[3:5],x1[5]) 

        q0 = (x0,q0[1],q0[2]) 
        q1 = (x1,q1[1],q1[2]) 

        D0 = cmp_generators__MultiMetric_(q0,q1)
        D1 = c0 - c1 
    
        return D0,D1 

    def cmp_stats(self,G0,G1):
        assert type(G0) == type(G1) == PRNGStats 

        q0,c0 = G0.latest_results() 
        q1,c1 = G1.latest_results() 

        D0 = cmp_generators__MultiMetric_(q0,q1)
        D1 = c0 - c1 

        return D0,D1 