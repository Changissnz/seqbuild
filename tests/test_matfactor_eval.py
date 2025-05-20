from mini_dm.matfactor_eval import * 
import unittest

### lone file test 
"""
python3 -m tests.test_matfactor_eval 
"""
###

mfe_testvar_1 = [({0, 2, 3, 4}, np.float64(2.0)), ({0, 2}, np.float64(2.0)), \
            ({1, 3}, np.float64(3.0)), ({0, 2, 3, 4}, np.float64(2.0)), \
            ({0, 1, 2}, np.float64(1.0)), ({0, 2}, np.float64(2.0)), \
            ({1, 3}, np.float64(3.0)), ({0, 1, 2}, np.float64(1.0))]

class MatFactorEvalMethods(unittest.TestCase):

    def test__is_self_factorable_list(self):
        L = [1, 3, 12, 15]
        L2 = [2,4,12,36,144]

        stat1 = is_self_factorable_list(L,True)
        stat2 = is_self_factorable_list(L2,True)
        assert not stat1 and stat2 

        stat1 = is_self_factorable_list(L)
        stat2 = is_self_factorable_list(L2)

        assert stat1 == [3] and stat2 == []

    def test__MatFactorEval__column_factor_map(self):
        M = np.array([\
            [1,1,1,5,7],\
            [2,3,2,10,14],\
            [2,3,2,15,20]])

        mfe = MatFactorEval(M)

        mfe.column_factor_map()

        assert mfe.cfm == \
        {0: {1: ([{0, 2, 3, 4}], [np.float64(2.0)]), \
            2: ([{0, 2}, {1, 3}], [np.float64(2.0), np.float64(3.0)])}, \
            1: {0: ([{0, 2, 3, 4}], [np.float64(2.0)]), \
            2: ([{0, 1, 2}], [np.float64(1.0)])}, \
            2: {0: ([{0, 2}, {1, 3}], [np.float64(2.0), np.float64(3.0)]), \
            1: ([{0, 1, 2}], [np.float64(1.0)])}}

        q = columnfactormap_to_sequence(mfe.cfm)
        assert q == mfe_testvar_1

    def test__columnfactor_rank(self): 
        q2 = columnfactor_rank(mfe_testvar_1) 

        ans1 = {0: np.float64(10.0), \
            2: np.float64(10.0), \
            3: np.float64(10.0), \
            4: np.float64(4.0), \
            1: np.float64(8.0)}

        assert q2 == ans1 

    def test__columnfactor_identity(self): 

        M = np.array([\
            [1,1,1,2,3],\
            [2,2,2,4,6],\
            [4,4,4,8,12]])

        mfe = MatFactorEval(M)
        mfe.column_factor_map()
        cfs = columnfactormap_to_sequence(mfe.cfm)
        q4 = columnfactor_identity(cfs,3) 
        assert q4 == [{0, 1, 2, 3, 4}]

        M = np.array([\
            [1,1,1],\
            [2,2,2]])
        mfe = MatFactorEval(M)
        mfe.column_factor_map()
        cfs = columnfactormap_to_sequence(mfe.cfm)
        q4 = columnfactor_identity(cfs,2) 
        assert q4 == [{0, 1, 2}]

    def test__MatFactorEval__identity_eval(self):
        M = np.array([\
            [1,1,1],\
            [2,2,2]])

        mfe = MatFactorEval(M)
        q1,q2 = mfe.identity_eval(rank_type=2) 
        assert q1 == {0: np.float64(4.0), 1: np.float64(4.0), 2: np.float64(4.0)}
        assert q2 == [{0, 1, 2}]
        assert mfe.id_col == [] 
        
        M = np.array([\
            [1,1,1,2,3],\
            [2,2,2,4,6],\
            [4,4,4,8,12]])
        mfe = MatFactorEval(M)
        q1,q2 = mfe.identity_eval(rank_type=2) 
        assert q1 == {0: np.float64(16.0), 1: np.float64(16.0), \
            2: np.float64(16.0), 3: np.float64(16.0), 4: np.float64(16.0)}
        assert q2 == [{0, 1, 2, 3, 4}]
        assert mfe.id_col == [] 

class MatrixConsistencyCheckMethods(unittest.TestCase):

    def test__MatrixConsistencyCheck__identity_eval(self):

        M = np.array([\
                    [1,1,1],\
                    [2,2,2]])

        Y = np.array([2,6])

        mfe = MatrixConsistencyCheck(M,Y)
        mfe.check() 
        assert mfe.inconsistent == [(0, 1)]
        M = np.array([\
                    [1,1,1],\
                    [2,2,2]])

        Y = np.array([2,4])

        mfe = MatrixConsistencyCheck(M,Y)
        mfe.check() 
        assert mfe.inconsistent == [] 

        M = np.array([\
                    [1,1,1],\
                    [2,2,2],\
                    [3,3,3]])

        Y = np.array([2,4,5])

        mfe = MatrixConsistencyCheck(M,Y)
        mfe.check() 
        assert mfe.inconsistent == [(0, 2), (1, 2)]


        M = np.array([\
                    [1,1,1],\
                    [2,2,2],\
                    [3,3,3]])

        mfe = MatrixConsistencyCheck(M,Y)
        mfe.check() 
        assert mfe.inconsistent == [(0, 2), (1, 2)]


        M = np.array([\
                    [1,0,1],\
                    [2,2,2],\
                    [3,3,3]])

        Y = np.array([2,4,5])

        mfe = MatrixConsistencyCheck(M,Y)
        mfe.check() 
        assert mfe.inconsistent == [(1, 2)]

        M = np.array([\
                    [1,0,1],\
                    [2,0,2],\
                    [3,3,3]])

        Y = np.array([2,4,5])

        mfe = MatrixConsistencyCheck(M,Y)
        mfe.check() 
        assert mfe.inconsistent == []


        M = np.array([\
                    [1,0,0],\
                    [0,0,0],\
                    [0,0,0]])

        Y = np.array([2,4,5])

        mfe = MatrixConsistencyCheck(M,Y)
        mfe.check() 
        assert mfe.inconsistent == []

        M = np.array([\
                    [1,0,0],\
                    [1,0,0],\
                    [0,0,0]])

        Y = np.array([2,4,5])

        mfe = MatrixConsistencyCheck(M,Y)
        mfe.check() 
        assert mfe.inconsistent == [(0, 1)]

if __name__ == '__main__':
    unittest.main()