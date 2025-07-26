from intigers.mdr_v2 import * 
from intigers.extraneous import * 
from intigers.tvec import *
from mini_dm.ag_ext import *  
from mini_dm.ngram import * 

class MultiMetric:

    def __init__(self,l):
        self.l = np.array(l)
        assert is_vector(self.l) 

    def summarize(self):
        return -1 