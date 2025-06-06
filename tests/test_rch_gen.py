from intigers.mod_prng import * 
from seqgen.rch_gen import * 
from morebs2.numerical_generator import prg__LCG,prg__n_ary_alternator

import unittest

def RCHAccuGen_argseq__caseX():
    n = 10
    prg3 = prg__LCG(11,31,17,1500) 
    num_nodes_range = [2,13] 
    dim_range = [3,12]
    ufreq_range = [2,11]
    qcap_range = [100,1001]

    return n,prg3,num_nodes_range,\
        dim_range,ufreq_range,qcap_range

### lone file test 
"""
python -m tests.test_rch_gen
"""
###
class RCHAccuGenMethods(unittest.TestCase):

    def test__RCHAccuGen__one_new_RCHAccuGen__v1__case1(self):
        num_nodes = 2
        dim_range = [3,6]
        prg = prg__LCG(3,4,5,33)
        ufreq_range = [2,3]
        mutrate = 1.0 
        queue_capacity = 1000 

        rg = RCHAccuGen.one_new_RCHAccuGen__v1(num_nodes,dim_range,prg,\
                ufreq_range,mutrate,queue_capacity)

        # size check with |<RChainHead>|
        rg.apply(16) 
        assert len(rg.acc_queue) == 3 

        rg.apply(30)
        rg.apply(233) 
        rg.apply(52) 

        # counter check for update 
        D =  {0: {0: 2}, 1: {0: 2}}
        assert rg.update_log == D,"GOT : {}".format(rg.update_log)

        # size check for queue
        for _ in range(1200): 
            rg.apply(prg())
        assert len(rg.acc_queue) == 1000 

    def test__RCHAccuGen__one_new_RCHAccuGen__v1__case2(self):

        num_nodes = 1
        dim_range = [3,6]
        prg = prg__LCG(3,4,5,33)
        ufreq_range = [2,3]
        mutrate = 1.0 
        queue_capacity = 1000 

        rg = RCHAccuGen.one_new_RCHAccuGen__v1(num_nodes,dim_range,prg,\
                ufreq_range,mutrate,queue_capacity)
        rch1 = deepcopy(rg.rch)

        x1 = rg.apply(16) 
        assert len(rg.acc_queue) == 2 

        # test for correct change of function 
        q1 = rg.rch.s[0].f(34)

        x2 = rg.apply(16) # change function here 

        q2 = rg.rch.s[0].f(34)

        x3 = rg.apply(16)

        assert q1 != q2
        assert x1 == x2
        assert x2 != x3 

    def test__RCHAccuGen__one_new_RCHAccuGen__v1__case3(self):
        prg3 = prg__LCG(1,4,7,150) 

        num_nodes = 8
        dim_range = [3,10]
        prg = prg__LCG(3,4,5,33)
        ufreq_range = [2,10]
        mutrate = 0.3#1.0 #0.3 
        queue_capacity = 1000 

        rg = RCHAccuGen.one_new_RCHAccuGen__v1(num_nodes,dim_range,prg,\
                ufreq_range,mutrate,queue_capacity)

        # run some numbers 
        for i in range(20): 
            x = prg3() 
            y = rg.apply(x)
            #print("X: ",x, " Y: ",y)

        # check for correct log of updates 
        D = {4: {0: 5}, 1: {0: 3}, 5: {0: 2}}
        assert rg.update_log == D 

        # check for correct number of <MutableRInstFunction>
        c = len([m for m in rg.mutgen if m != set()])
        assert c == 3

    """
    inconclusive test, mainly for demonstrating output values. 
    """
    def test__RCHAccuGen__generate_n_instances__case1(self):

        # set each as RCHAccuGen 
        argx = RCHAccuGen_argseq__caseX()
        rgs = RCHAccuGen.generate_n_instances(argx[0],argx[1],argx[2],\
            argx[3],argx[4],argx[5])

        prg2 = prg__n_ary_alternator(1,15,3)

        for _ in range(10):
            x = prg2()
            print("input: ",x*2) 
            for (i,r) in enumerate(rgs):
                y = r.apply(x*2)
                print("\t ",y)
            print() 

    """
    uses a <ModPRNGOutputter> to test out output from 10 
    <RCHAccuGen.apply> instances. 
    """
    def test__RCHAccuGen__generate_n_instances__case2(self):
        argx = RCHAccuGen_argseq__caseX()
        rgs = RCHAccuGen.generate_n_instances(argx[0],argx[1],argx[2],\
            argx[3],argx[4],argx[5])
        rgsf = [rgs_.apply for rgs_ in rgs]

        prg = prg__n_ary_alternator(0,len(rgs),0)
        mpo = ModPRNGOutputter(rgsf)

        prg2 = prg__LCG(4,-3,11,59)
        for i in range(10):
            q = next(mpo)
            rx = q(prg2()) 
            assert type(rx) in {int,np.int32} 

if __name__ == '__main__':
    unittest.main()