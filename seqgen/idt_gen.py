from intigers.mod_prng import *
from intigers.idt_proc import * 

DEFAULT_FOREST_NEWSEQ_NUMCENTERS = [2,15]
DEFAULT_FOREST_NEWSEQ_MULTRANGE = [-20,20]

class IDecForest: 

    def __init__(self,s,prng_outputter,cache_size,prg,prg2=None):
        assert type(s) == IntSeq
        assert len(s) >= 2 
        assert type(prng_outputter) == ModPRNGOutputter
        self.ST = []
        self.S = s 
        self.T = None 
        self.mpo = prng_outputter
        self.cache_size = cache_size
        self.prg = prg 
        self.prg2 = prg2

    """
    outputs n values
    """
    def output_n(self,n):  
        return -1

    #---------------------- tree and integer sequence creation

    def one_tree(self): 
        self.one_tree_()
        self.ST.append((self.S,self.T))

    def one_tree_(self): 
        # fetch arguments for <IntSeq2Tree> conversion 
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

        # set the `entryf` functions for the new tree 
        self.T = cnvrt.root 
        TNode.dfs(self.T,False,False,True,set_attr=('entryf',self.mpo.__next__))
        return self.T

    def one_new_IntSeq(self,gentype:int):
        assert gentype in {1,2,3} 
        return -1

    """

    """
    @staticmethod
    def derive_new_IntSeq__type_fm(S,prg):
        assert type(S) == IntSeq

        # set number of centers
        ncr = (DEFAULT_FOREST_NEWSEQ_NUMCENTERS[0],\
            min([len(S),DEFAULT_FOREST_NEWSEQ_NUMCENTERS[1]]))
        c_ = modulo_in_range(prg(),ncr)

        # fetch the centers, each a factor
        isfso = ISFactorSetOps(S.l,int_limit=NPINT32_MAX)
        isfso.factor_count_()
        q = isfso.dsort(keys=None)

        cx = []
        while c_ > 0:
            i = prg() % len(q)
            cx.append(q.pop(i)[0])
            c_ -= 1 

        # generate the sequence 
        nrange = (int(ceil(len(S) / 2)), len(S) * 2)
        n = modulo_in_range(prg(),nrange)

        rx = modulo_in_range(prg(),[0,1001])
        rx = rx / 1000.0 

        iseq = prg__integer_seq__mult(n,cx,rx,None,\
            DEFAULT_FOREST_NEWSEQ_MULTRANGE,prg,\
            num_attempts_per_nc=150)
        iseq = intlist_no_dups_no_zero(iseq)
        return IntSeq(iseq)

    #------------------------ integer generation

    # choose seq for tree 
    def process_seq_at_tree__sequential(self,tindex:int,S):
        return -1 

    def process_seq_at_tree__inflow(self,tindex:int,S):
        return -1

    def process_seq_at_tree__splat(self,tindex:int,S):
        return -1 

    #------------------------- 