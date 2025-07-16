
import numpy as np  
from morebs2.numerical_generator import * 
from morebs2.search_space_iterator import * 

"""
Partitioned Uniform Crawler through a solution space. 
Crawler takes a positive real number `unit`, and searches 
the vicinity of vector `v`. This search is comprised of 
a sequence of vectors. Each element in this sequence is 
formatted in the manner of `partition`, a variable that 
specifies output dimensionality, before being outputted 
by the `next(PUCrawler)` call. 
"""
class PUCrawler: 

    def __init__(self,v,unit:float,partition): 
        assert is_vector(v) 
        assert unit > 0. 
        assert type(unit) in {float,np.float32,np.float64} 
        assert type(partition) in {tuple,list} 
        if type(partition) == tuple:
            assert len(partition) == 2 
            assert partition[0] * partition[1] == len(v) 
        else: 
            assert sum(partition) == len(v) 

        self.v = v 
        self.unit = unit 
        self.partition = partition 
        self.ssi = None 
        self.init_ssi() 
        return 
    
    def init_ssi(self):
        bx = self.v - self.unit 
        bx1 = self.v + self.unit
        bx2 = bx1 + self.unit 
        bounds = np.array([bx,bx2]).T 
        startPoint = bounds[:,0] 
        columnOrder = [i for i in range(len(self.v))]  
        self.ssi = SearchSpaceIterator(bounds,startPoint,\
            columnOrder,SSIHop=3,cycleOn=False,\
            cycleIs=0)
        
    def __next__(self):
        if self.ssi.reached_end(): return None
        v = next(self.ssi) 
        pv = PUCrawler.vec_to_partition_(v,self.partition) 
        return pv 

    @staticmethod 
    def vec_to_partition_(v,p): 

        if type(p) == tuple: 
            return v.reshape(p) 

        i = 0 
        q = []
        for p_ in p:
            v_ = v[i:i+p_] 
            q.append(deepcopy(v_)) 
            i = i + p_ 
        return q 