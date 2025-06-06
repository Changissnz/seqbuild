from intigers.idt_proc import * 

class IDecForest: 

    def __init__(self,s,prg,prg2=None):
        assert type(s) == IntSeq
        assert len(s) >= 2 
        self.ST = []
        self.S = s 
        self.T = None 
        self.prg = prg 
        self.prg2 = prg2

    """
    outputs n values
    """
    def output_n(self,n):  

        return -1

    #---------------------- tree and integer sequence creation

    def one_tree(self): 
        I = self.S
        
        l,d = None,None 
        rnge = [2,ceil(len(I) / 2)] 
        if rnge[0] == rnge[1]: rnge[1] += 1 

        q = modulo_in_range(self.prg(),rnge)
        if self.prg() % 2: 
            l = q
        else: 
            d = q     

        prgx = self.prg2 if type(self.prg2) != type(None) else self.prg 
        cnvrt = IntSeq2Tree(I,l,q,prgx)
        cnvrt.convert()
        self.T = cnvrt.root 
        return self.T

    def one_new_IntSeq(self,gentype:int):
        assert gentype in {1,2,3} 
        return -1

    # choose seq for tree 
    def process_seq_at_tree(self,tindex:int,S):
        return -1 

    #------------------------- 