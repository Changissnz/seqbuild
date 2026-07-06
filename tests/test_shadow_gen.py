from seqgen.shadow_gen import * 
from morebs2.numerical_generator import * 
from intigers.mod_prng import * 
from intigers.extraneous import * 
import unittest

'''
def ShadowGen__sample_X(fitting_struct): 
    prg = prg__LCG(87,51,-311,40040)
    file_path = "dummy_file.txt"
    sg = ShadowGen(prg,file_path,fitting_struct)
    return sg 
'''

### lone file test 
"""
py -m tests.test_shadow_gen 
"""
###

class ShadowGenMethods(unittest.TestCase):

    """
    1 unique value in kernel file, 
    4 unique values in cyclical aux. PRNG 
    """
    def test__ShadowGen__next__case1(self):

        prg = prg__iterable([5,6,10,12])  

        file_path = "dummy_file4.txt"

        # subcase 1: 98 unique values over 2000 iterations 
        fitting_struct = "fvec" 
        sg = ShadowGen(prg,file_path,fitting_struct)

        X = [] 
        for _ in range(2000): 
            x = next(sg) 
            X.append(x) 
            if _ % 100 == 0: print("#1.1-{}: {}".format(_,x)) 

        assert len(set(X)) == 98
        sg.close() 

        # subcase 2: 180 unique values over 2000 iterations 
        fitting_struct = "tvec"
        sg = ShadowGen(prg,file_path,fitting_struct)

        X = [] 
        for _ in range(2000): 
            x = next(sg) 
            X.append(x) 
            if _ % 100 == 0: print("#1.2-{}: {}".format(_,x)) 

        assert len(set(X)) == 180 
        sg.close() 

        # subcase 3: 1392 unique values over 2000 iterations 
        fitting_struct = "optri"
        sg = ShadowGen(prg,file_path,fitting_struct)

        X = [] 
        for _ in range(2000): 
            x = next(sg) 
            X.append(x) 
            if _ % 100 == 0: print("#1.3-{}: {}".format(_,x)) 

        assert len(set(X)) == 1334, "got {}".format(len(set(X)))
        sg.close() 

        # subcase 4: 320 unique values over 2000 iterations 
        prg = prg__iterable([5,6,10,12])  

        fitting_struct = "mdrv2"
        sg = ShadowGen(prg,file_path,fitting_struct)

        X = [] 
        for _ in range(2000): 
            x = next(sg) 
            X.append(x) 
            if _ % 100 == 0: print("#1.3-{}: {}".format(_,x)) 

        assert len(set(X)) == 320, "got {}".format(len(set(X)))
        sg.close() 

        # subcase 5: 320 unique values over 2000 iterations 
        #           check for last 10 different from that of `mdrv2`
        prg = prg__iterable([5,6,10,12])  

        fitting_struct = "mdr"
        sg = ShadowGen(prg,file_path,fitting_struct)

        X2 = [] 
        for _ in range(2000): 
            x = next(sg) 
            X2.append(x) 
            if _ % 100 == 0: print("#1.3-{}: {}".format(_,x)) 

        assert len(set(X2)) == 320, "got {}".format(len(set(X)))
        sg.close() 

        D = np.round(np.array(X2[-10:]) - np.array(X[-10:]),5) 
        assert list(D).count(0) == 1, "got {}".format(D) 

    """
    20 unique values in kernel file, 
    30 unique values in cyclical aux. PRNG 
    """
    def test__ShadowGen__next__case2(self): 
        prg = prg__n_ary_alternator(-15,15,-15)

        file_path = "dummy_file.txt"
        
        # case 1: 
        fitting_struct = "fvec"
        #fitting_struct = "mdr" 

        sg = ShadowGen(prg,file_path,fitting_struct)

        #Q = np.array([next(sg) for _ in range(2000)] )

        X = [] 
        for _ in range(2000): 
            x = next(sg) 
            X.append(x) 
            if _ % 100 == 0: print("#2.1-{}: {}".format(_,x)) 
        assert len(set(X)) == 1764
        sg.close() 

        # case 2: 
        fitting_struct = "tvec"
        prg = prg__n_ary_alternator(-15,15,-15)
        sg = ShadowGen(prg,file_path,fitting_struct)

        X = [] 
        for _ in range(2000): 
            x = next(sg) 
            X.append(x) 
            if _ % 100 == 0: print("#2.2-{}: {}".format(_,x)) 
        assert len(set(X)) == 1315
        sg.close() 

        # case 3: 
        fitting_struct = "optri"
        prg = prg__n_ary_alternator(-15,15,-15)
        sg = ShadowGen(prg,file_path,fitting_struct)

        X = [] 
        for _ in range(2000): 
            x = next(sg) 
            X.append(x) 
            if _ % 100 == 0: print("#2.3-{}: {}".format(_,x)) 
        assert len(set(X)) == 175 
        sg.close() 

        # case 4: 
        fitting_struct = "mdrv2" 
        prg = prg__n_ary_alternator(-15,15,-15)
        sg = ShadowGen(prg,file_path,fitting_struct)

        X = [] 
        for _ in range(2000): 
            x = next(sg) 
            X.append(x) 
            if _ % 100 == 0: print("#2.4-{}: {}".format(_,x)) 
        assert len(set(X)) == 1998 
        sg.close() 

if __name__ == '__main__':
    unittest.main()
