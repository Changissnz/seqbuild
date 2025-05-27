from intigers.pof_autogen import * 
import unittest
from morebs2.numerical_generator import prg__constant,prg__n_ary_alternator,LCG
import random 

### lone file test 
"""
python -m tests.test_pof_autogen
"""
###
class LCGMethodTests(unittest.TestCase):

    def test___LCG_next(self): 
        print("GuudFudamis")
        random.seed(1202)
        prg = LCG(5,390)
        prg1 = LCG(32,3901)
        prg2 = LCG(75,391014456456)
        prg3 = LCG(518,1473719390)

        D = {0: [122, 184], 1: [3431, 3371], \
        2: [358024985255, 24069917960],\
        3: [453575051, 1425679607]} 


        for i in range(10):
            for (j,x) in enumerate([prg,prg1,prg2,prg3]):
                #print(next(x),x.multiplier, x.increment)
                next(x) 
                assert x.multiplier == D[j][0]
                assert x.increment == D[j][1]



class POFAutogenMethods(unittest.TestCase):

    def test__POFV2ConditionAutoGen__integerpair_op__case1(self):
        prg = prg__n_ary_alternator(s0=-5,s1=29,start=0)
        pofv2_ca = POFV2ConditionAutoGen(prg) 
        i1,i2 = 5000,300
        q = pofv2_ca.integerpair_op(i1,i2)
        c = 0 
        for q_ in q: 
            assert q_.is_solved() 
            c += 1 
        assert c == 5 
        return


if __name__ == '__main__':
    unittest.main()