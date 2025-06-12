
from seqgen.idt_gen import *
from seqgen.rch_gen import * 

def IDecForest_sampleZ(L=None,make_first_tree:bool=True,prg1=None):
    if type(L) == type(None):
        L = [3, 10, 395, 273, 275, 404, 278, 280, 154, \
            416, 33, 162, 38, 300, 44, 174, 304, 180, 182, \
            54, 314, 60, 61, 62, 63, 65, 72, 336, 96, 356, \
            360, 362, 365, 238, 114, 242, 252]
    else: 
        assert type(L) in {list,np.ndarray}
        assert len(L) >= 2 

    argx = RCHAccuGen_argseq__caseX() 
    s = IntSeq(L) 

    rgs = RCHAccuGen.generate_n_instances(argx[0],argx[1],argx[2],\
        argx[3],argx[4],argx[5])
    rgsf = [rgs_.apply for rgs_ in rgs]

    mpo = ModPRNGOutputter(rgsf)

    if type(prg1) == type(None): 
        prg1 = prg__LCG(13,87,211,5000)
    #prg2 = None 
    #prg3 = prg__LCG(31,83,311,1500)

    idf = IDecForest(s,mpo,100,[15,50],15,prg1,prg2=None)

    x3,x4,x5 = None,None,None 
    if make_first_tree: 
        idf.one_tree() 
        x3 = idf.one_new_IntSeq(0)
        x4 = idf.one_new_IntSeq(1)
        x5 = idf.one_new_IntSeq(2)
    return idf,(x3,x4,x5)