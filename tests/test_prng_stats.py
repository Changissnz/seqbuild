from morebs2.numerical_generator import prg__constant
from face.prng_stats import * 
from seqgen.rch_gen import * 
from intigers.lcg_v3 import * 
from seqgen.optri_gen import * 
import unittest

def PRNGIOStats__template_Q(prg_struct):
    prg_seq = [prg__constant(10),prg__LCG(1,2,3,20),prg__single_to_int(prg__LCG(56,-142,566,2112))]

    return PRNGIOStats(prg_seq,"prg",prg_struct,ngram_info=[3,4,5],\
        partition_size=5,absdiff=True) 


### lone file test 
"""
py -m tests.test_prng_stats 
"""
###
class PRNGIOStatsMethods(unittest.TestCase):

    """
    demonstrates relative randomness comparisons on a nilpotent <RCHAccuGen> 
    """
    def test__PRNGIOStats__cmp__case1(self):

        prg = prg__LCG(56,-142,566,2112)
        num_nodes = 5
        dim_range = [4,8]
        ufreq_range = [2,6] 

        rch = RCHAccuGen.one_new_RCHAccuGen__v1(num_nodes,dim_range,prg,\
            ufreq_range,mutrate=0.5,queue_capacity=1000)

        pios = PRNGIOStats__template_Q(rch)
        pios.run()    

        # check that input PRNG relative randomness checks out 
        assert pios.input_cmp(0,1,{1,2,3,4,5,6,7}) == -1 
        assert pios.input_cmp(0,2,{1,2,3,4,5,6,7}) == -1 
        assert pios.input_cmp(1,2,{1,2,3,4,5,6,7}) == -1 

        # compare output randomness 
        assert pios.output_cmp(0,1,{1,2,3,4,5,6,7}) == 0 
        assert pios.output_cmp(0,2,{1,2,3,4,5,6,7}) == 0 
        assert pios.output_cmp(1,2,{1,2,3,4,5,6,7}) == 0 

        # compare I/O
        assert pios.io_cmp(0,{1,2,3,4,5,6,7}) == 0
        assert pios.io_cmp(1,{1,2,3,4,5,6,7}) == -1 
        assert pios.io_cmp(2,{1,2,3,4,5,6,7}) == -1 

        # compaire pairwise I/O
        assert pios.pairwise_io_cmp(0,1,{1,2,3,4,5,6,7}) == 1 
        assert pios.pairwise_io_cmp(0,2,{1,2,3,4,5,6,7}) == 1
        assert pios.pairwise_io_cmp(1,2,{1,2,3,4,5,6,7}) == -1 

    """
    demonstrates relative randomness comparisons on a case of <RCHAccuGen> where 
    input generator #2 results in a more random output generator than with generator #3. 
    """
    def test__PRNGIOStats__cmp__case2(self):

        prg = prg__LCG(451.11,-155.44,546.222,2001.004) 
        num_nodes = 5
        dim_range = [4,8]
        ufreq_range = [2,6] 

        rch = RCHAccuGen.one_new_RCHAccuGen__v1(num_nodes,dim_range,prg,\
            ufreq_range,mutrate=0.5,queue_capacity=1000)

        pios = PRNGIOStats__template_Q(rch)
        pios.run()  

        assert pios.output_cmp(0,1,{1,2,3,4,5,6,7})  == -1 
        assert pios.output_cmp(0,2,{1,2,3,4,5,6,7})  == -1 
        assert pios.output_cmp(1,2,{1,2,3,4,5,6,7}) == 1 

        assert pios.io_cmp(0,{1,2,3,4,5,6,7}) == 0 
        assert pios.io_cmp(1,{1,2,3,4,5,6,7}) == 1 
        assert pios.io_cmp(2,{1,2,3,4,5,6,7}) == -1

        assert pios.pairwise_io_cmp(0,1,{1,2,3,4,5,6,7}) == -1 
        assert pios.pairwise_io_cmp(1,2,{1,2,3,4,5,6,7}) == 1 
        assert pios.pairwise_io_cmp(0,1,{1,2,3,4,5,6,7}) == -1

    """
    case demonstrates how LCGV3, with its objective of satisfying its PRNG-assigned 
    trinary vector, is measured as not as random as its input PRNGs. 
    """
    def test__PRNGIOStats__cmp__case3(self):

        prg = prg__LCG(451.11,-155.44,546.222,2001.004) # 0,-1,1 
        lv3 = LCGV3(0,3,17,-997,1112,prg,super_range=[-1500,1500])

        pios = PRNGIOStats__template_Q(lv3)
        pios.run()    

        assert pios.output_cmp(0,1,{1,2,3,4,5,6,7}) == -1 
        assert pios.output_cmp(0,2,{1,2,3,4,5,6,7}) == -1 
        assert pios.output_cmp(1,2,{1,2,3,4,5,6,7}) == -1 

        assert pios.io_cmp(0,{1,2,3,4,5,6,7}) == 1 
        assert pios.io_cmp(1,{1,2,3,4,5,6,7}) == -1 
        assert pios.io_cmp(2,{1,2,3,4,5,6,7}) == -1 

        assert pios.pairwise_io_cmp(0,1,{1,2,3,4,5,6,7}) == 1 
        assert pios.pairwise_io_cmp(0,2,{1,2,3,4,5,6,7}) == 1 
        assert pios.pairwise_io_cmp(1,2,{1,2,3,4,5,6,7}) == -1 

    """
    demonstrates another curious case, an <OpTriGen> where output PRNG is not necessarily measured 
    as more rnadom than input PRNG
    """
    def test__PRNGIOStats__cmp__case4(self):
        
        prg = prg__LCG(451.11,-155.44,546.222,2001.004) # 0,-1,1 

        l = [2,6,4,3,8,10]
        seq = IntSeq(l) 
        ot = seq.optri(sub,np.int32)

        random.seed(100)
        prg = prgen(0,101)

        otg = OpTriGen(2,ot,prg,2,\
            forward_func=add,backward_func=sub,add_noise=False)

        pios = PRNGIOStats__template_Q(otg)
        pios.run()    

        assert pios.output_cmp(0,1,{1,2,3,4,5,6,7}) == 1 
        assert pios.output_cmp(0,2,{1,2,3,4,5,6,7}) == 1 
        assert pios.output_cmp(1,2,{1,2,3,4,5,6,7}) == -1 

        assert pios.io_cmp(0,{1,2,3,4,5,6,7}) == 1 
        assert pios.io_cmp(1,{1,2,3,4,5,6,7}) == 1 
        assert pios.io_cmp(2,{1,2,3,4,5,6,7}) == 0  

        assert pios.pairwise_io_cmp(0,1,{1,2,3,4,5,6,7}) == 1
        assert pios.pairwise_io_cmp(0,2,{1,2,3,4,5,6,7}) == 1
        assert pios.pairwise_io_cmp(1,2,{1,2,3,4,5,6,7}) == 1 

    """
    same OpTriGen> as in case 4, except `add_noise` is set to True 
    """
    def test__PRNGIOStats__cmp__case5(self):
        
        prg = prg__LCG(451.11,-155.44,546.222,2001.004) 

        l = [2,6,4,3,8,10]
        seq = IntSeq(l) 
        ot = seq.optri(sub,np.int32)

        random.seed(100)
        prg = prgen(0,101)

        otg = OpTriGen(2,ot,prg,2,\
            forward_func=add,backward_func=sub,add_noise=True)

        pios = PRNGIOStats__template_Q(otg)
        pios.run()    

        assert pios.output_cmp(0,1,{1,2,3,4,5,6,7}) == 1 
        assert pios.output_cmp(0,2,{1,2,3,4,5,6,7}) == 0 
        assert pios.output_cmp(1,2,{1,2,3,4,5,6,7}) == -1 

        assert pios.io_cmp(0,{1,2,3,4,5,6,7}) == 1 
        assert pios.io_cmp(1,{1,2,3,4,5,6,7}) == 1 
        assert pios.io_cmp(2,{1,2,3,4,5,6,7}) == 0  

        assert pios.pairwise_io_cmp(0,1,{1,2,3,4,5,6,7}) == 1
        assert pios.pairwise_io_cmp(0,2,{1,2,3,4,5,6,7}) == 1
        assert pios.pairwise_io_cmp(1,2,{1,2,3,4,5,6,7}) == 1 


if __name__ == '__main__':
    unittest.main()