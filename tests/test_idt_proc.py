from seqgen.rch_gen import * 
from intigers.idt_proc import * 
from morebs2.numerical_generator import prg__LCG 

import unittest

def IntSeq2Tree__caseX():
    prg = prg__LCG(3,7,9,33)
    prg2 = prg__n_ary_alternator(-20,30,-19) 

    q = [] 
    for i in range(50): 
        q_ = prg()
        q2_ = prg2() 
        q.append((q_ * q2_) % 422)
    q = list(set(q) - {0})

    ns = IntSeq(q) 
    l = 4 
    d = None  
    is2t = IntSeq2Tree(ns,l,d,prg,verbose=False) 
    is2t.convert()
    return is2t 

def integer_sequence__caseY():

    prg = prg__LCG(3,7,9,33)
    prg2 = prg__n_ary_alternator(-20,30,-19) 

    q = [] 

    # 
    for i in range(50): 
        q_ = prg()
        q2_ = prg2() 
        q.append((q_ * q2_) % 422)
    return list(set(q) - {0})

### lone file test 
"""
python -m tests.test_idt_proc
"""
###
class IDTProcMethods(unittest.TestCase):

    def test__IDTProc__travel_path__case1(self):
        prg3 = prg__LCG(1,4,7,150) 

        num_nodes = 8
        dim_range = [3,10]
        prg = prg__LCG(3,4,5,33)
        ufreq_range = [2,10]
        mutrate = 0.3
        queue_capacity = 1000 

        rg = RCHAccuGen.one_new_RCHAccuGen__v1(num_nodes,dim_range,prg,\
                ufreq_range,mutrate,queue_capacity)

        def set_v():
            return rg.apply

        is2t = IntSeq2Tree__caseX()
        qx1,qx2,_ = TNode.dfs(is2t.root,False,True,True,set_attr=('entryf',set_v))

        qx = IDTProc(is2t.root)
        
        # iso-output test 
        p,_ = qx.process_value(300)
        p2 = qx.travel_path(300,p[1:])
        p3 = qx.travel_path(240,p[1:])
        p4 = qx.travel_path(240,p[1:])

        assert p2 != p3 
        assert p3 != p4

    def test__IDTProc__travel_path__case2(self):
        prg3 = prg__LCG(1,4,7,150) 

        num_nodes = 2
        dim_range = [3,10]
        prg = prg__LCG(3,4,5,33)
        ufreq_range = [20,30]
        mutrate = 1.0 
        queue_capacity = 1000 

        rg = RCHAccuGen.one_new_RCHAccuGen__v1(num_nodes,dim_range,prg,\
                ufreq_range,mutrate,queue_capacity)

        def set_v():
            return rg.apply


        ###
        is2t = IntSeq2Tree__caseX()
        qx1,qx2,_ = TNode.dfs(is2t.root,False,True,True,set_attr=('entryf',set_v))

        qx = IDTProc(is2t.root)

        p,_ = qx.process_value(300)
        p2 = qx.travel_path(300,p[1:])
        p3 = qx.travel_path(240,p[1:])
        p4 = qx.travel_path(240,p[1:])
        assert p3 == p4 

    def test__IDTProc__next__case1(self):
        q = integer_sequence__caseY()
        ns = IntSeq(q) 
        l = None#4 
        d = 6 

        prg = prg__LCG(3,7,9,33)
        is2t = IntSeq2Tree(ns,l,d,prg,verbose=False) 
        is2t.convert()

        times2 = lambda x: 2 * x 
        def f_output():
            return times2

        TNode.dfs(is2t.root,False,True,True,set_attr=('entryf',f_output)) 

        q0 = q[:9]
        itp = IDTProc(is2t.root)
        itp.inflow_set(set(q0))

        qr = None 
        qr = next(itp)
        qr2 = next(itp)

        assert len(qr) == 1 
        assert len(qr2) == 4 

        cx = sum([len(v) for v in qr.values()])
        cx2 = sum([len(v) for v in qr2.values()])
        assert cx == cx2
        assert cx == len(q0)

        vx1 = flatten_setseq([set(v) for v in qr.values()])
        vx2 = flatten_setseq([set(v) for v in qr2.values()])
        assert vx1 == vx2

        q0x = set([2 * q for q in q0])
        assert vx1 == set(q0x)," got {} actual {}".format(vx1,set(q0))

if __name__ == '__main__':
    unittest.main()