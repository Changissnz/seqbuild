from mini_dm.affine_fit_search import * 
from morebs2.numerical_generator import * 
from intigers.extraneous import * 
from desi.nvec_gen import * 
import unittest 


def sample_affine_dataset_NX(prg,num_points): 

    ro_prg_ = prg__LCG(3,-13,42,990)
    ro_prg_ = prg__LCG(4,-13,42,990)
    ro_prg = prg__single_to_range_outputter(ro_prg_) 

    ro_prg2 = prg__LCG(37,-131,4211,9800)
    ro_prg2 = prg__single_to_range_outputter(ro_prg2) 

    psg = PointSetGen__TypeAffine(num_points,prg,ro_prg,ro_prg2)

    X = psg.generate_points(is_ordered=True,clear_data=False)

    I = np.array(psg.input_seq)
    O = np.array(psg.point_seq)

    return I,O 

def AffineFitSearch__sample_TRI(idn,I,O,k,prg,ma_dim_,score_improvement_reference,spread_per_guess=DEFAULT_AFFINE_FIT_SPG): 
    assert idn in {0,1,2} 

    ma_dim = list(ma_dim_) 
    n = O.shape[1] 
    if ma_dim[0] == 1: 
        ma_dim[0] = n 
    if ma_dim[1] == 1: 
        ma_dim[1] = n 
    
    if idn == 0: 
        afst = AffineFitSearch__TypeN2NA(I,O,prg,k,ma_dim,0,"random",is_bfs_queue=True,\
            spread_per_guess = spread_per_guess,score_improvement_reference = score_improvement_reference,\
            score_improvement_type="greedy")
    elif idn == 1: 
        afst = AffineFitSearch__TypeN2NA(I,O,prg,k,ma_dim,0,(0.5,0.5),1.0,is_bfs_queue=False,\
            spread_per_guess= spread_per_guess,score_improvement_reference=score_improvement_reference,\
            score_improvement_type="greedy")
    elif idn == 2: 
        afst = AffineFitSearch__TypeN2NA(I,O,prg,k,ma_dim,1,(0.5,0.5),0.2,is_bfs_queue=False,\
            spread_per_guess= spread_per_guess,score_improvement_reference=score_improvement_reference,\
            score_improvement_type="stochastic")

    return afst 


### lone file test 
"""
py -m tests.test_affine_fit_search 
"""
###
class AffineFitSearch__TypeN2NAMethods(unittest.TestCase):

    # greedy search
    # dim(X,Y) = (0,9)
    def test__AffineFitSearch__TypeN2NA__case_1(self): 
        print("Case #1: Greedy")

        prg = prg__LCG(6,3,2,350)
        I,O = sample_affine_dataset_NX(prg,1200) 
        k = 5 

        afst = AffineFitSearch__sample_TRI(0,I,O,k,prg,(1,1),"local")

        V0 = afst.soln_error_vec() 

        for i in range(1000): 
            if afst.finstat: break 
            ##print("{}: {}: {}".format(i,default_cfunc2(afst.best__global_error()),len(afst.queue)))
            next(afst)

        # termination @ 73 
        assert i == 73 

        V1 = afst.soln_error_vec() 

        # compare errors of pre-search to post-search solution 
        assert np.sum(V1) < np.sum(V0) 

        assert np.all(V0 == np.array([ 59248.07389,63508.16637,\
            85991.31259,110407.033,124146.6061])) 

        assert np.all(V1 == np.array([59248.07389,61552.60723,\
            63508.16637,63857.14056,65812.6997]))

    # stochastic 
    # same (X,Y) dataset as in Case 1 
    def test__AffineFitSearch__TypeN2NA__case_2(self): 
        print("Case #2: Stochastic")

        prg = prg__LCG(6,3,2,350)
        I,O = sample_affine_dataset_NX(prg,1200) 
        k = 5 

        afst = AffineFitSearch__sample_TRI(2,I,O,k,prg,(0,1),"local")


        V0 = afst.soln_error_vec() 

        for i in range(1000): 
            if afst.finstat: break 
            #print("{}: {}: {}".format(i,default_cfunc2(afst.best__global_error()),len(afst.queue)))
            next(afst)

        assert i == 999  

        V1 = afst.soln_error_vec() 

        assert np.sum(V1) < np.sum(V0) 

        assert np.all(V0 == np.array([81065.64347,82807.87267,107855.81141,\
            112996.34028,122200.55412])) 

        assert np.all(V1 == np.array([81065.64347,82807.87267,92135.1689,\
            92161.46216,92266.63519]))

    # greedy  
    # dim(X,Y) = (5,5) 
    def test__AffineFitSearch__TypeN2NA__case_3(self): 

        print("Case #3: Greedy")

        # (5,5)
        prg = prg__LCG(2,-112,271,3502)
        I,O = sample_affine_dataset_NX(prg,1200) 
        k = 5 

        afst = AffineFitSearch__sample_TRI(0,I,O,k,prg,(1,1),"local")

        V0 = afst.soln_error_vec() 

        for i in range(200): 
            if afst.finstat: break 
            #print("{}: {}: {}".format(i,default_cfunc2(afst.best__global_error()),len(afst.queue)))
            next(afst)

        assert i == 199 

        V1 = afst.soln_error_vec() 

        assert np.sum(V1) < np.sum(V0) 

        assert np.all(V0 == np.array([533234.97053, 554229.18479, 555355.71567, \
            613552.69112,615846.7739 ])) 

        assert np.all(V1 == np.array([306622.55567, 335887.73373, 336509.49373, \
            351619.85621, 362826.82882])) 

    # stochastic for first 20 rounds, then greedy for next 200 rounds 
    # same dataset as in case 3 
    def test__AffineFitSearch__TypeN2NA__case_4(self): 

        print("Case #4: Stochastic+Greedy")

        # (5,5)
        prg = prg__LCG(2,-112,271,3502)
        I,O = sample_affine_dataset_NX(prg,1200) 
        k = 5 

        afst = AffineFitSearch__sample_TRI(0,I,O,k,prg,(0,1),"local")

        V0 = afst.soln_error_vec() 

        afst.score_imp_type = "stochastic"
        for i in range(20): 
            if afst.finstat: break 

            ##print("{}: {}: {}".format(i,default_cfunc2(afst.best__global_error()),len(afst.queue)))
            next(afst)

        afst.score_imp_type = "greedy" 
        for i in range(200): 
            if afst.finstat: break 

            ##print("{}: {}: {}".format(i,default_cfunc2(afst.best__global_error()),len(afst.queue)))
            next(afst)

        assert i == 199

        V1 = afst.soln_error_vec() 

        assert np.sum(V1) < np.sum(V0) 

        assert np.all(V0 == np.array(([585842.31002, 588256.39759, \
            593173.24262,609567.5662,609781.41277]))) 

        assert np.all(V1 == np.array([278667.78373, 278691.78373, 321468.74373,\
            321516.74373,321516.74373]))



if __name__ == '__main__':
    unittest.main()