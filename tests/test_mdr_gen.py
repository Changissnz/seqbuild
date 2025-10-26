from seqgen.mdr_gen import * 
from intigers.mdr_v2 import * 
from morebs2.numerical_generator import prg__n_ary_alternator

import unittest
import time 

random.seed(100) 

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
        assert mdrg.sequence_size() == 50, "got {}".format(mdrg.sequence_size())

        intseed_cycle = [5]
        mdrg.generate_sequence__type_novelgen(intseed_cycle,\
            max_sequences=50,clear_prev_info=True)
        assert mdrg.sequence_size() == 50, "got {}".format(mdrg.sequence_size())

    def test__MDRGen__next_gentype2__case1(self): 
       
        prg = prgen(2,100)
        l = [5, 11, 23, 9, 29, 89, 13, 62, 49, 245, 5, 11, 23, 47, 143, 431, 271, 686, -2717, 10895]
        
        intsq = IntSeq(l) 
        md = ModuloDecomp(intsq)
        md.merge(False)
        mdr = ModuloDecompRepr(md)

        mdrg = MDRGen(mdr,prg,False,gentype=2,preproc=True) 
        q = [int(next(mdrg)) for _ in range(30)]

        ans = [5, 11, 23, 9, 29, 89, 13, 62, 49, 245, 5, 11, 23, 47, 143, 431, 271, 686, -2717, 10895, 5, 11, 23, 47, 95, 191, -581, -5794, 2063, 9442]        
        assert q == ans 

    def test__MDRGen__next_gentype2__case2(self):
        
        xt = [5, 23, 23, 2, 15, 15, 24303, 2, 5, 11, 23, 71, 215, 135, 672, 3, 2, 6, 647, 5]       
        prg = prgen(2,100)
        l = [2,5,11,4,14,44,6,27,3,15] # unstable 2 6 
        intsq = IntSeq(l) 
        md = ModuloDecomp(intsq)
        md.merge(False)
        mdr = ModuloDecompRepr(md)
        mdrg = MDRGen(mdr,prg,False,gentype=2,gt2_rcswitch=True,\
            gt2_sel1=True,gt2_sel2=True,gt2_sel3=False)

        q = [int(next(mdrg)) for _ in range(30)]
        ans = [60, 12, 487, 31, 2019, 60, 94, 11703, 31, 63, 127, 255, \
            511, -913, -3327, -6441, -12669, 28, 57, 115, 231, 463, 927, \
            2777, 155, 3401, 4121, 60, 121, 243]
        assert q == ans 

    def test__MDRGen__next__speedtest(self):
        lx = [14,1400,15,-43,-410,32,108,90,934]
        mdx = ModuloDecompV2(IntSeq(lx),False,DEFAULT_MDRGEN_MAXMULT)
        mdr = ModuloDecompRepr(mdx,reconstruct_type=2)

        prg = prg__n_ary_alternator(-0,50,0) 

        t = time.time()
        print("MDRG#1")
        mdrg = MDRGen(mdr,prg,exclude_neg=False,gentype=1) 
        q = np.array([next(mdrg) for _ in range(500)]) 

        print("MDRG#2")
        mdrg2 = MDRGen(mdr,prg,exclude_neg=False,gentype=2,\
            gt2_rcswitch=False,gt2_sel1=True,gt2_sel2=True,\
            gt2_sel3=True,gt2_seed_in_output=True)
        q2 = np.array([next(mdrg2) for _ in range(500)]) 

        print("MDRG#3")
        mdrg3 = MDRGen(mdr,prg,exclude_neg=False,gentype=2,\
            gt2_rcswitch=True,gt2_sel1=False,gt2_sel2=True,\
            gt2_sel3=True,gt2_seed_in_output=True)
        q3 = np.array([next(mdrg3) for _ in range(500)]) 

        print("MDRG#4")
        mdrg4 = MDRGen(mdr,prg,exclude_neg=False,gentype=2,\
            gt2_rcswitch=True,gt2_sel1=True,gt2_sel2=False,\
            gt2_sel3=True,gt2_seed_in_output=True)
        q4 = np.array([next(mdrg4) for _ in range(500)]) 

        print("MDRG#5")
        mdrg5 = MDRGen(mdr,prg,exclude_neg=False,gentype=2,\
            gt2_rcswitch=True,gt2_sel1=True,gt2_sel2=True,\
            gt2_sel3=False,gt2_seed_in_output=False)
        q5 = np.array([next(mdrg5) for _ in range(500)]) 

        print("MDRG#6")
        mdrg6 = MDRGen(mdr,prg,exclude_neg=False,gentype=2,\
            gt2_rcswitch=True,gt2_sel1=True,gt2_sel2=True,\
            gt2_sel3=True,gt2_seed_in_output=False)
        q6 = np.array([next(mdrg6) for _ in range(500)]) 

        print("MDRG#7")
        mdrg7 = MDRGen(mdr,prg,exclude_neg=False,gentype=2,\
            gt2_rcswitch=True,gt2_sel1=False,gt2_sel2=False,\
            gt2_sel3=True,gt2_seed_in_output=False)
        q7 = np.array([next(mdrg7) for _ in range(500)]) 

        t2 = time.time() 

        print("\t\tTIME: ",t2 -t)
        assert t2 - t < 15 

if __name__ == '__main__':
    unittest.main()