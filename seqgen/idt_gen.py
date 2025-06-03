from intigers.idt_proc import * 

class IDecForest: 

    def __init__(self,s,prg,prg2=None):
        assert type(s) == IntSeq
        assert len(s) >= 2 
        self.SX = [] 
        self.TS = [] 
        self.S = s 
        self.T = None 
        self.prg = prg 
        self.prg2 = prg2

    def one_tree(self): 
        I = self.s 
        
        l,d = None,None 
        rnge = [2,ceil(len(I) / 2)] 
        if rnge[0] == rnge[1]: 
            rnge[1] += 1 
        q = modulo_in_range(self.prg(),rnge)
        
        if self.prg() % 2: 
            l = q
        else: 
            d = q     

        prgx = self.prg2 if type(self.prg2) != type(None) else self.prg 
        self.T = IntSeq2Tree(I,l,q,prgx)
        return self.T