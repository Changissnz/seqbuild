"Partitioned Uniform Crawler through a solution space"

import numpy as np  
from morebs2.numerical_generator import * 

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

        self.unit = unit 
        self.partition = partition 
        return 