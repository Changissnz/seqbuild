from face.comm_lang import * 
import unittest

### lone file test 
"""
python -m tests.test_comm_lang 
"""
###
class CommLangMethods(unittest.TestCase):

    def test__MAKE_proc__case1(self):
        cmd = "make lcg with 400,532,31,4577"
        splitstr_cmd = cmd.split(" ")
        var_map = dict() 

        prg = MAKE_proc(splitstr_cmd,var_map)

        lx = []

        for _ in range(12): 
            lx.append(prg()) 

        assert equal_iterables(lx, \
            [400.0, 2289.0, 297.0, 2417.0, 4315.0, \
            2534.0, 2481.0, 1747.0, 304.0, 1564.0, \
            3642.0, 1504.0])

    def test__OPEN_proc__case1(self):

        fp = "face/sample_script/commondeth.txt"
        cmd2 = "open file {}".format(fp)
        splitstr_cmd2 = cmd2.split(" ")
        fi = OPEN_proc(splitstr_cmd2)
        assert type(fi) == io.BufferedWriter 
        fi.close()

    def test__CommLangParser__process_command__case1(self): 
        clp = CommLangParser("face/sample_script/commond_one.txt") 

        clp.load_next_command()
        #'set G = make lcg with 400,532,31,4577.'
        q = clp.process_command()
        assert q[0] == "G" 
        assert type(q[1]) in {MethodType,FunctionType}

        clp.close() 

    """
    tests for `lcg` and `multimetric` commands. 
    """
    def test__CommLangParser__process_command__case2(self):

        clp = CommLangParser("face/sample_script/commond_one.txt") 

        clp.load_next_command()
        #'set G = make lcg with 400,532,31,4577.'
        q = clp.process_command()
        assert q[0] == "G" 
        assert type(q[1]) in {MethodType,FunctionType}

        clp.load_next_command() 
        q2 = clp.process_command() 

        clp.load_next_command() 
        q3 = clp.process_command() 

        assert "G" in clp.vartable
        assert "V" in clp.vartable

        clp.load_next_command() 
        q3 = clp.process_command() 

        clp.load_next_command() 
        q4 = clp.process_command() 

        q4sol0 = np.array([0.77917 , 0.34337 , 0.31333]) 
        assert np.all(np.round(q4[1][0] - q4sol0,5) == 0) 

        q4sol1 = (137, 153)
        assert q4sol1 == q4[1][1] 

        q4sol2 = np.float64(0.705)
        assert q4sol2 == q4[1][2] 

        clp.close() 

    """
    tests for `lcg` and `mdr` commands. 
    """
    def test__CommLangParser__process_command__case3(self):

        clp = CommLangParser("face/sample_script/commond_two.txt") 

        clp.load_next_command()
        q = clp.process_command()
        assert q[0] == "G" 
        assert type(q[1]) in {MethodType,FunctionType}

        clp.load_next_command() 
        q2 = clp.process_command() 

        clp.load_next_command() 
        q3 = clp.process_command() 

        assert "G" in clp.vartable
        assert "V" in clp.vartable

        clp.load_next_command() 
        q3 = clp.process_command() 

        clp.load_next_command() 
        q4 = clp.process_command() 

        q3sol = [np.int32(300), np.int32(1495), \
            np.int32(8876), np.int32(882), np.int32(1707), np.int32(303), \
            np.int32(731), np.int32(1006), np.int32(8029), np.int32(3321), \
            np.int32(9825), np.int32(35), np.int32(14), np.int32(10), np.int32(6), \
            np.int32(-198), np.int32(20), np.int32(-495), np.int32(680), np.int32(486)]
        q4sol = [np.int32(47), np.int32(-23), np.int32(-232), np.int32(1128), \
            np.int32(2199), np.int32(269), np.int32(629), np.int32(700), np.int32(5581), \
            np.int32(2473), np.int32(7281), np.int32(292), np.int32(271), np.int32(267), \
            np.int32(263), np.int32(59), np.int32(277), np.int32(-238), np.int32(92), \
            np.int32(-102)]

        assert q3[1] == q3sol 
        assert q4[1] == q4sol 

        clp.close() 

    """
    tests for `lcgv3` command. 
    """
    def test__CommLangParser__process_file__case1(self): 

        clp = CommLangParser("face/sample_script/commond_three.txt") 
        clp.process_file()

        assert "V" in clp.vartable
        assert "V2" in clp.vartable
        assert "V3" in clp.vartable
        assert "V4" in clp.vartable
        assert "V5" in clp.vartable

        V = np.array(clp.vartable["V"])
        V2 = np.array(clp.vartable["V2"])
        V3 = np.array(clp.vartable["V3"])
        V4 = np.array(clp.vartable["V4"])
        V5 = np.array(clp.vartable["V5"])

        ix = np.where(V == V2)[0] 
        ix2 = np.where(V == V3)[0] 
        ix3 = np.where(V == V4)[0] 
        ix4 = np.where(V2 == V3)[0]  
        ix5 = np.where(V2 == V4)[0]  
        ix6 = np.where(V3 == V4)[0]  
        ix7 = np.where(V3 == V5)[0]  

        assert len(ix) == 0 and len(ix2) == 0 and len(ix3) == 0 
        assert len(ix4) == 1000 and len(ix5) == 1000 and len(ix6) == 1000 
        assert len(ix7) == 0

        clp.close() 

    """
    tests for keyword `convert`
    """
    def test__CommLangParser__process_file__case2(self): 
        clp = CommLangParser("face/sample_script/commond_four.txt") 
        clp.process_file()

        R = clp.vartable['R']
        for r in R: 
            assert is_valid_range(r,False,False) 

        R2 = clp.vartable['R2']
        R2sol = np.array([0, 3, 8, 5]) 
        R2sol2 = np.array([1, 1, 5, 9]) 
        assert np.all(R2[-2] == R2sol)
        assert np.all(R2[-1] == R2sol2) 
        for r in R2: 
            assert len(r) == 4 

        R3 = clp.vartable['R3']
        R3sol = np.array([67843., 18047.,  4907., 81609., 32115., 62765.])
        assert np.all(R3[0] == R3sol) 
        for r in R3: 
            assert len(r) == 6 

        R4 = clp.vartable['R4']
        R4sol = np.array([-1,  0, -1,  0,  0, -1,  0])
        R4sol2 = np.array([ 0, -1, -1,  1, -1, -1,  1])
        assert np.all(R4[0] == R4sol)
        assert np.all(R4[1] == R4sol2)
        for r in R4: 
            assert len(r) == 7 

        clp.close() 

    """
    tests for `qval` command. 
    """
    def test__CommLangParser__process_file__case3(self): 
        clp = CommLangParser("face/sample_script/commond_five.txt") 
        clp.process_file()

        q = clp.vartable["Q"]

        qsol = [429,39797,8817,19535,65501] 

        qs = []
        for _ in range(5): 
            qs.append(int(round(next(q)))) 

        assert qs == qsol 
        clp.close()

    """
    tests for `lcgv2` and `mdrv2` commands. 
    """
    def test__CommLangParser__process_file__case4(self): 
        clp = CommLangParser("face/sample_script/commond_six.txt") 
        clp.process_file()

        mv1 = clp.vartable['MV1']
        mv2 = clp.vartable['MV2']
        assert mv1 == mv2

        mv11 = clp.vartable['MV11']
        mv12 = clp.vartable['MV12']
        assert mv11 == mv12

        mv21 = clp.vartable['MV21']
        mv22 = clp.vartable['MV22']
        assert mv21 == mv22

        clp.close()

    """
    tests for `mdrgen` command. 
    """
    def test__CommLangParser__process_file__case5(self): 
        clp = CommLangParser("face/sample_script/commond_seven.txt") 
        clp.process_file()

        V2 = clp.vartable['V2']

        V2sol = [67,7110,6701,5968,2217,5630]
        assert len(V2) == 150 
        assert V2[:6] == V2sol

        V3 = clp.vartable['V3']
        V3sol = [6388,683457,744880,744147,740396]

        assert V3[:5] == V3sol 

        clp.close() 

    """
    tests for `optri` command. 
    """
    def test__CommLangParser__process_file__case6(self): 

        clp = CommLangParser("face/sample_script/commond_eight.txt") 
        clp.process_file()

        OT = clp.vartable['OT']
        V2 = clp.vartable['V2']
        V3 = clp.vartable['V3']
        V4 = clp.vartable['V4']

        X2 = clp.vartable['X2']
        X3 = clp.vartable['X3']
        X4 = clp.vartable['X4']

        assert not np.any(np.round(np.abs(X2[0] - X3[0]),5) == 0.)  
        assert not np.any(np.round(np.abs(X2[0] - X4[0]),5) == 0.)  
        assert not np.any(np.round(np.abs(X3[0] - X4[0]),5) == 0.)  

        clp.close() 

    """
    tests for `pid` command. 
    """
    def test__CommLangParser__process_file__case7(self): 
        clp = CommLangParser("face/sample_script/commond_nine.txt") 
        clp.process_file()
        V = clp.vartable['V']

        qm = np.mean(V)
        assert qm == 49712.56 
        clp.close() 

    """
    tests for `write` keyword. 
    """
    def test__CommLangParser__process_file__case8(self): 
        clp = CommLangParser("face/sample_script/commond_ten.txt") 
        clp.process_file()

        p1 = "command_the_fourth.txt" 
        p2 = "command_the_force.txt"
        assert os.path.exists(p1)
        assert os.path.exists(p2)

        clp.close() 

    """
    tests for command 
        `set <generator> for <start range,end range>.`. 
    """
    def test__CommLangParser__process_file__case9(self): 
        clp = CommLangParser("face/sample_script/commond_13.txt") 
        clp.process_file()
        V = np.array(clp.vartable['V'])
        V2 = np.array(clp.vartable['V2'])
        assert not np.any(V==V2)

        for v in V2: 
            assert v >= -200. and v < 6000.01

        clp.close() 

    """
    tests for processing erroneous file 
    """
    def test__CommLangParser__process_file__case10(self): 
        clp = CommLangParser("face/sample_script/commond_12error.txt") 
        clp.process_file()
        clp.close() 

    def test__CommLangParser__process_file__case11(self):
        clp = CommLangParser("face/sample_script/commond_14.txt") 
        clp.process_file() 

        for i in range(1,6): 
            s = "O" + str(i) 
            assert s in clp.vartable
        clp.close() 
        return

    def test__CommLangParser__process_file__case12(self):
        clp = CommLangParser("face/sample_script/commond_16.txt") 
        clp.process_file() 

        assert "X" in clp.vartable 

        X = clp.vartable['X'] 

        Xsol = [3.,22.,98.,402.,1618.,1484.,\
            948.,470.,224.,906.,302.,1218.]

        assert X == Xsol

        clp.close() 
        return  

    def test__CommLangParser__process_file__case13(self): 
        clp = CommLangParser("face/sample_script/commond_18.txt") 
        clp.process_file()

        q = clp.vartable['q'] 
        assert q[2] == 0 

        clp.close() 

    def test__CommLangParser__process_file__case14(self): 
        clp = CommLangParser("face/sample_script/commond_21.txt") 
        clp.process_file()

        v2 = clp.vartable['V2']
        v3 = clp.vartable['V3']
        assert None in v2 and None not in v3 
        clp.close() 

    def test__CommLangParser__process_file__case15(self):
        clp = CommLangParser("face/sample_script/commond_22.txt")
        clp.process_file() 

        vsol = [36,156,260,285,958,434,45,750,384,653,821,185]
        V = clp.vartable['V']
        assert V == vsol 

if __name__ == '__main__':
    unittest.main()
