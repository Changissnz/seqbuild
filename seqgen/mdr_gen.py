from intigers.seq_struct import * 

class MDRGen: 

    """
    mdr := ModuloDecompRepr
    prgen := function, pseudo-random integer generator. Call with `prgen()`.  
    """
    def __init__(self,mdr,prgen):
        assert type(mdr) == ModuloDecompRepr
        self.mdr = mdr
        self.prg = prgen 
        self.cache = [] 

    def __next__(self): 
        if len(self.cache) == 0: 
            r = self.mdr.reconstruct() 
            self.mdr.first = np.int32(self.prg())
            self.cache.extend(r) 
        q = self.cache.pop(0) 
        return q 
