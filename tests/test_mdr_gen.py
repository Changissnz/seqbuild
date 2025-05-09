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

if __name__ == '__main__':
    unittest.main()