from seqgen.gg_gen import * 
import unittest

### lone file test 
"""
python -m tests.test_gg_gen
"""
###
class AGV2GuidedGenMethods(unittest.TestCase):

    def test__AGV2GuidedGen__next__guidedrepl__case1(self):
                
        base_prg = prg__LCG(677,251,734,41242)
        aux_prg = prg__LCG(78,63,27,10000) 
        base_output_span = 10 
        density_measure = np.array(((5,30),(2,1)))

        agg = AGV2GuidedGen(base_prg,aux_prg,base_output_span,density_measure,\
                ngram_length=4)
        agg.set_refvar("cov")
        agg.set_permuter() 
        qr = agg.next__guidedrepl()
        qs = []

        for _ in range(15): 
            qr = agg.next__guidedrepl()
            assert not np.all(qr == agg.base_seq)
            #print("QR")
            #print(qr) 
            #print("BASE")
            #print(agg.base_seq)


if __name__ == '__main__':
    unittest.main()