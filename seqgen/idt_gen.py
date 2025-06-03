from intigers.idt_proc import * 

class IDecForest: 

    def __init__(self,s):
        assert type(s) == IntSeq
        self.S = [s] 
    
    def one_tree(self): 
        return -1 


IntSeq2Tree: 

    def __init__(self,intseq,l:int,d:int,prg,verbose:bool=False):