"""
file contains code to extend the <RChainHead> class from the 
library <morebs2>
"""

class MutableRInstFunction:

    def __init__(self): 
        return -1 

class RCHAccuGen: 

    def __init__(self,rch,perm_prg,acc_queue=[]):
            assert type(acc_queue) == list 
            self.rch = rch 
            self.perm_prg = perm_prg 
            self.acc_queue = [] 
            return 

        def fetch_varlist(self):
            return -1 

        def apply(self,x): 
            self.rch.apply(x)
            vx = deepcopy(self.rch.vpath)
            self.acc_queue.extend(vx)
            self.update()
            return 

        def update(self): 
            return -1 

        def __next__(self):
            return -1