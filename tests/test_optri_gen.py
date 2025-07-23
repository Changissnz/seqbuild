from seqgen.optri_gen import * 
import unittest


### lone file test 
"""
python3 -m tests.test_optri_gen
"""
###

class OpTriGenMethods(unittest.TestCase):

    def test__triangular_matrix_perimeter(self):
        l = [2,5,1,8,13] 
        seq = IntSeq(l) 
        ot = seq.optri(sub,np.int32)
        print("OT")
        print(ot)

        c1 = triangular_matrix_perimeter(ot,0,True)
        c2 = triangular_matrix_perimeter(ot,1,True)
        c3 = triangular_matrix_perimeter(ot,2,True)

        sol1 = np.array([3,  -4,   7,   5,\
            -2, -13, -31,  18,  -7], dtype=np.int32)
        sol2 = np.array([  5,  -2, -13, -31,\
            18,  -7,   3,  -4,   7], dtype=np.int32)
        sol3 = np.array([-31,  18,  -7,   3,\
            -4,   7,   5,  -2, -13], dtype=np.int32)
        assert np.all(c1==sol1)
        assert np.all(c2==sol2)
        assert np.all(c3==sol3)

        c4 = triangular_matrix_perimeter(ot,0,False)
        c5 = triangular_matrix_perimeter(ot,1,False)
        c6 = triangular_matrix_perimeter(ot,2,False)

        sol4 = np.array([  3,  -7,  18, -31,\
            -13,  -2,   5,   7,  -4], dtype=np.int32)
        sol5 = np.array([-31, -13,  -2,   5,\
            7,  -4,   3,  -7,  18], dtype=np.int32)
        sol6 = np.array([  5,   7,  -4,   3,\
            -7,  18, -31, -13,  -2], dtype=np.int32)

        assert np.all(c4==sol4)
        assert np.all(c5==sol5)
        assert np.all(c6==sol6)
        return 

    def test__OpTri45N90Split__j_derivative_seq(self):
        split = ([12,13,14],[15,16,17]) 
        ots = OpTri45N90Split(split)
        s = ots.j_derivative_seq(1)

        assert s == [12, 25, 52, 94, 152, 227]
        assert ots.derivative_seqs == \
            {3: [14, 15, 16, 17],\
            2: [13, 27, 42, 58, 75],\
            1: [12, 25, 52, 94, 152, 227]}

    def test__OpTriFlipDerivation__construct(self):
        l = [2,6,4,3,8,10]
        seq = IntSeq(l) 
        ot = seq.optri(sub,np.int32)
        otfd = OpTriFlipDerivation(ot,2,add,0)

        otfd.construct() 
        sol = np.array([[-21, -23, -24, -19, -17],\
            [  0, -29, -28, -22, -25],\
            [  0,   0, -19, -14, -23],\
            [  0,   0,   0, -12, -26],\
            [  0,   0,   0,   0, -10]], dtype=np.int32)
        assert np.all(otfd.m_ == sol) 

        otfd.reset_axis(1)
        otfd.construct() 
        sol = np.array([[ 10,   6,   8,   9,   4],\
            [  0,   8,  14,  13,   7],\
            [  0,   0,  11,   4,  -1],\
            [  0,   0,   0,  -5,  -3],\
            [  0,   0,   0,   0, -17]], dtype=np.int32)
        assert np.all(otfd.m_ == sol)

    '''
    case 1 is of generator type #2 (flip-derivation)
    '''
    def test__OpTriGen__next__case1(self):
                
        l = [2,6,4,3,8,10]
        seq = IntSeq(l) 
        ot = seq.optri(sub,np.int32)

        random.seed(100)
        prg = prgen(0,101)

        otg = OpTriGen(2,ot,prg,2,\
            forward_func=add,backward_func=sub,add_noise=False)

        l = [] 

        q = next(otg)
        assert otg.reloaded
        l.append(q) 

        q = next(otg)
        assert not otg.reloaded
        l.append(q) 

        m1 = np.array([\
            [-21, -23, -24, -19, -17],\
            [  0, -29, -28, -22, -25],\
            [  0,   0, -19, -14, -23],\
            [  0,   0,   0, -12, -26],\
            [  0,   0,   0,   0, -10]], dtype=np.int32)

        assert np.all(otg.otfd.m_ == m1)

        sz = otg.batch_size()
        s = 2 

        mx = otg.otfd.mf 
        while not otg.reloaded: 
            q = next(otg) 
            l.append(q) 
            s += 1 

        assert s - 1 == sz 

        q = next(otg)
        s += 1 

        solx = np.array([\
            [ -22,    3,  -20,    4,  -25],\
            [   0,  -24,  -72,  -25,  -78],\
            [   0,    0, -128,  -33, -133],\
            [   0,    0,    0, -142, -337],\
            [   0,    0,    0,    0, -546]], dtype=np.int32)
        assert np.all(solx == otg.otfd.m_) 

        q = next(otg)
        l.append(q) 
        s += 1 
        
        sol_mx = np.array([\
            [ -27,   29,  -24,   23,  -25],\
            [ -56,   53,  -47,   48,    0],\
            [-109,  100,  -95,    0,    0],\
            [-209,  195,    0,    0,    0],\
            [-404,    0,    0,    0,    0]], dtype=np.int32)
        assert np.all(otg.otfd.mf == sol_mx) 

        sol_mx2 = np.array([\
            [ -22,    3,  -20,    4,  -25],\
            [   0,  -24,  -72,  -25,  -78],\
            [   0,    0, -128,  -33, -133],\
            [   0,    0,    0, -142, -337],\
            [   0,    0,    0,    0, -546]], dtype=np.int32)
        assert np.all(otg.otfd.m_ == sol_mx2) 

        while not otg.reloaded: 
            q = next(otg) 
            l.append(q) 
            s += 1 

        assert s - 1 == sz * 2 

        sol_m = np.array([\
            [-547, -524, -548, -519, -546],\
            [   0, -472, -519, -466, -522],\
            [   0,    0, -511, -411, -520],\
            [   0,    0,    0, -207, -416],\
            [   0,    0,    0,    0, -402]], dtype=np.int32)
        assert np.all(sol_m == otg.otfd.m_)
        assert np.all(sol_mx2 == otg.otfd.m2_)


    '''
    case 2 is of generator type #1 (jagged 45-90 split)
    '''
    def test__OpTriGen__next__case2(self):
        random.seed(100) 
        l = [2,6,4,3,8,10]
        seq = IntSeq(l) 
        ot = seq.optri(sub,np.int32)
        prg = prgen(0,101)
        otg = OpTriGen(2,ot,prg,1,\
            forward_func=add,backward_func=sub,add_noise=False)

        l = [] 
        s = 0
        q = next(otg)
        s += 1 
        assert otg.reloaded
        l.append(q) 

        q = next(otg)
        s += 1
        assert not otg.reloaded
        l.append(q) 

        bs1 = otg.batch_size()

        while not otg.reloaded: 
            q = next(otg) 
            s += 1 
            l.append(q)

        assert s - 1 == bs1 

        bs2 = otg.batch_size()
        q = next(otg)
        s += 1

        assert not otg.reloaded

        while not otg.reloaded: 
            q = next(otg) 
            s += 1 
            l.append(q)
        l = np.array(l,dtype=np.int32)

        assert s - 1 == bs1 + bs2 

        sol = np.array([ 4, -2,  1,  5, -9, -9, \
            -2,  1,  5, -9,  4,  2,  1,  6, 13, \
            22, 31, 39, 50, 59, -2, -1,  5,  7, \
            9,  9,  8, 11,  9,  1,  6,  2,  2, \
            0, -1,  3, -2,  5, -4,  0, -2, -1, \
            4, -5,  7], dtype=np.int32)

        assert np.all(sol == l)

#############################################################
    '''
    case 3 is of generator type #2 (flip-derivation). 
    Case tests PRNG-less mode.
    '''
    def test__OpTriGen__next__case3(self):

        random.seed(100)
        l = [2,6,4,3,8,10]

        seq = IntSeq(l) 
        ot = seq.optri(sub,np.int32)
        prg = prgen(0,101)
        otg = OpTriGen(2,ot,None,2,\
            forward_func=add,backward_func=sub,add_noise=False)

        l = [] 
        q = next(otg)
        l.append(q)

        otfd_ = deepcopy(otg.otfd)
        otfd2 = otfd_.reproduce(otg.bfunc) 
        otfd2.construct() 

        for i in range(29): 
            q = next(otg)

        if i == 14: assert otg.reloaded 
        else: assert not otg.reloaded 

        assert not otg.binary_alternator
        l.append(q) 

        q = next(otg)
        l.append(q) 
        assert otg.binary_alternator

        x = [] 
        for (i,r) in enumerate(otg.otfd.m_): 
            x.extend(r[i:])

        assert np.all(otg.otfd.m_ == otfd2.m_)

        for i in range(14):
            q = next(otg)
            l.append(q) 

        assert l[-15:] == x 

        for i in range(15):
            q = next(otg)
            l.append(q) 

        qt = otg.otfd.m_
        qt2 = otg.otfd.m2_ 

        q = np.array([otg.intseed],dtype=np.int32)
        for r in qt[0]: 
            q = np.append(q, sub(r,q[-1])) 
        source_seq = IntSeq(q)
        optri = source_seq.optri(sub,np.int32)

        otfd4 = OpTriFlipDerivation(optri,otg.intseed,add,0)
        otfd4.construct()  

        q = next(otg)
        assert np.all(otg.otfd.m_ == otfd4.m_)

    '''
    case 4 is of generator type #1 (jagged 45-90 split). 
    Case tests PRNG-less mode.
    '''
    def test__OpTriGen__next__case4(self):

        l = [2,6,4,3,8,10]
        seq = IntSeq(l) 
        ot = seq.optri(sub,np.int32)

        random.seed(100)
        prg = prgen(0,101)

        otg = OpTriGen(2,ot,prg,1,\
            forward_func=add,backward_func=sub,add_noise=False)

        l = [] 
        q = next(otg)
        l.append(q)

        for i in range(16): 
            q = next(otg)


        sol = np.array([\
            [ 4,  2,  1,  6, 13, 22, 31, 39, 50, 59],\
            [ 0, -2, -1,  5,  7,  9,  9,  8, 11,  9],\
            [ 0,  0,  1,  6,  2,  2,  0, -1,  3, -2],\
            [ 0,  0,  0,  5, -4,  0, -2, -1,  4, -5],\
            [ 0,  0,  0,  0, -9,  4, -2,  1,  5, -9]], dtype=np.int32)

        assert np.all(otg.ots.to_matrix() == sol)

        if i == 14: assert otg.reloaded 
        else: assert not otg.reloaded 


if __name__ == '__main__':
    unittest.main()