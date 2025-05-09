from seqgen.mdr_gen import * 
import unittest,random

random.seed(100) 

def prgen(s1,s2):
    def prgen_():  
        return random.randint(s1,s2) 
    return prgen_

### lone file test 
"""
python3 -m tests.test_mdr_gen
"""
###
class MDRGenMethods(unittest.TestCase):

    def test__MDRGen__generate_sequence__type_novelgen(self):
        prg = prgen(2,100)

        l = [2,5,11,4,14,44,6,27,3,15] # unstable 2 6 
        intsq = IntSeq(l) 
        md = ModuloDecomp(intsq)
        md.merge(False)
        mdr = ModuloDecompRepr(md)

        mdrg = MDRGen(mdr,prg,False) 

        intseed_cycle = [3,-2]
        mdrg.generate_sequence__type_novelgen(intseed_cycle,\
            max_sequences=50,clear_prev_info=True)
        assert mdrg.sequence_size() == 50 

        intseed_cycle = [5]
        mdrg.generate_sequence__type_novelgen(intseed_cycle,\
            max_sequences=100,clear_prev_info=True)
        assert mdrg.sequence_size() == 8 

    def test__MDRGen__next_gentype2__case1(self):
        prg = prgen(2,100)
        l = [2,5,11,4,14,44,6,27,3,15] # unstable 2 6 
        intsq = IntSeq(l) 
        md = ModuloDecomp(intsq)
        md.merge(False)
        mdr = ModuloDecompRepr(md)

        mdrg = MDRGen(mdr,prg,False,gentype=2)

        intseed_cycle = [5]
        mdrg.generate_sequence__type_novelgen(intseed_cycle,\
            max_sequences=100,clear_prev_info=True)

        q = mdrg.seed2seq[np.int32(2)][[0,1],:]
        q = np.insert(q,0,[2,2],axis=1) 
        q = q.flatten() 
        for i in range(20): 
            q2 = next(mdrg)
            assert q[i] == q2 

    def test__MDRGen__next_gentype2__case2(self):
        
        xt = [2,27,672,3232,3232,2,44,215,215,215,2,5,11,23,71,647,2741,13705,5,49]
        #random.seed(100)
        prg = prgen(2,100)
        l = [2,5,11,4,14,44,6,27,3,15] # unstable 2 6 
        intsq = IntSeq(l) 
        md = ModuloDecomp(intsq)
        md.merge(False)
        mdr = ModuloDecompRepr(md)
        mdrg = MDRGen(mdr,prg,False,gentype=2,gt2_rcswitch=True,\
            gt2_sel1=True,gt2_sel2=True,gt2_sel3=False)

        intseed_cycle = [5]
        mdrg.generate_sequence__type_novelgen(intseed_cycle,\
            max_sequences=100,clear_prev_info=True)
        q = mdrg.seed2seq[np.int32(2)][[0,1],:]
        q = np.insert(q,0,[2,2],axis=1) 
        q = q.flatten() 
        for i in range(20): 
            q2 = next(mdrg)
            assert xt[i] == q2, "have {}, want {}".format(q2,xt[i]) 

if __name__ == '__main__':
    unittest.main()