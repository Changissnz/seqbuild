from seqgen.gg_gen import * 
import unittest

def AGV2GuidedGen__sample_z(): 
    base_prg = prg__LCG(677,251,734,41242)
    aux_prg = prg__LCG(78,63,27,10000) 
    base_output_span = 10 
    density_measure = 30 

    agg = AGV2GuidedGen(base_prg,aux_prg,base_output_span,density_measure,\
            ngram_length=4)
    return agg 

### lone file test 
"""
python -m tests.test_gg_gen
"""
###
class AGV2GuidedGenMethods(unittest.TestCase):

    def test__AGV2GuidedGen__next__guidedrepl__case1(self):
                
        agg = AGV2GuidedGen__sample_z() 
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

    def test__AGV2GuidedGen__set_permuter__case1(self):

        agg = AGV2GuidedGen__sample_z()
        agg.set_refvar("cov")

        qvec = []
        for _ in range(5): 
            agg.set_permuter() 
            qs = agg.permuter.apply()
            qvec.append(qs)

        for (i,q) in enumerate(qvec):
            for (j,q2) in enumerate(qvec):
                if i == j: 
                    continue 
                assert not np.any(q == q2) 

    def test__AGV2GuidedGen__next__guidedrepl__case2(self):
        base_prg = prg__LCG(677,251,734,41242)
        aux_prg = prg__LCG(78,63,27,10000) 
        base_output_span = 10 
        density_measure = np.array(((5,30),(2,1)))
        density_measure = 30 

        agg = AGV2GuidedGen(base_prg,aux_prg,base_output_span,density_measure,\
                ngram_length=7,is_context_fed=True)
        agg.set_refvar(5)

        qs = []

        for _ in range(100): 
            qr = agg.next__guidedrepl(update_base_seq=True)
            agg.set_permuter() 
            qs.append(qr) 

        #qs_sol = [4, 3, 5, 3, 4, 3, 4, 4, 3, 3, 4, 4, \
        #        4, 2, 3, 4, 4, 4, 4, 4, 4, 3, 4, 3, 4, \
        #        4, 4, 4, 4, 3, 4, 4, 3, 2, 4, 4, 4, 4, \
        #        4, 3, 4, 4, 4, 3, 4, 3, 2, 3, 3, 4, 5, \
        #        4, 5, 3, 3, 4, 4, 4, 3, 4, 4, 4, 3, 4, \
        #        4, 4, 2, 3, 3, 3, 4, 3, 4, 4, 3, 3, 2, \
        #        4, 3, 4, 4, 4, 3, 3, 4, 3, 4, 3, 3, 4, \
        #        4, 4, 2, 4, 4, 3, 3, 3, 3, 4]
        qs_sol = [5, 4, 5, 4, 4, 4, 4, 5, 5, 4, 5, 4, 4, \
            4, 4, 4, 5, 5, 5, 4, 5, 4, 4, 4, 4, 5, 5, 4, \
            5, 4, 5, 5, 4, 4, 5, 5, 5, 4, 4, 4, 5, 5, 4, \
            4, 5, 4, 4, 4, 4, 4, 5, 5, 5, 5, 4, 5, 4, 4, \
            5, 5, 4, 5, 4, 5, 4, 5, 3, 4, 5, 4, 4, 4, 4, \
            5, 4, 4, 4, 5, 4, 4, 5, 5, 3, 4, 4, 4, 4, 3, \
            3, 5, 4, 4, 4, 4, 5, 4, 4, 4, 4, 4]


        assert agg.agd_log.refvar_catvec == qs_sol,"got {}".format(agg.agd_log.refvar_catvec)

if __name__ == '__main__':
    unittest.main()