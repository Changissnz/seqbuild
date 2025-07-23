from mini_dm.n2m_index import * 
from morebs2.numerical_generator import prg__LCG

import unittest

### lone file test 
"""
python -m tests.test_n2m_index
"""
###
class N2MIndexMapGenMethods(unittest.TestCase):

    def test__N2MIndexMapGen__one_new_relation__case1(self):
        nm = (5,9)
        index_degree_range = [3,5]
        ns_range = [2,4]
        ms_range = [3,6]
        prg = prg__LCG(43,100,31,511)
        index_degree_is_geq = False 
        nimg = N2MIndexMapGen(nm,index_degree_range,\
                ns_range,ms_range,prg,\
                index_degree_is_geq)

        stat = True 
        while stat: 
            r = nimg.one_new_relation() 
            if type(r) == type(None):
                stat = False 

        for _ in range(9): 
            q = nimg.n2m_imap.degree_of_mindex(0)
            assert q >= index_degree_range[0] and q <= index_degree_range[1]

        for q in nimg.n2m_imap.n2m_map:
            sz = len(string_to_vector(q[0])) 
            assert sz >= ns_range[0] and sz <= ns_range[1]

    def test__N2MIndexMapGen__make__case1(self):

        nm = (15,92)
        index_degree_range = [7,15] 
        ns_range = [2,14]
        ms_range = [3,16]
        prg = prg__LCG(431,10120,3101,511000)
        index_degree_is_geq = False 
        nimg = N2MIndexMapGen(nm,index_degree_range,\
                ns_range,ms_range,prg,\
                index_degree_is_geq=False)

        nimg.make() 
        dm = nimg.n2m_imap.mindex_degree_map()
        tm = nimg.target_index_degree_map

        c = 0 
        c2 = 0 
        for k,v in dm.items():
            if v == tm[k]: 
                c += 1
            
            if v >= ms_range[0] and v <= ms_range[1]:
                c2 += 1 

        assert c == 91#84, "got {}".format(c)
        assert c2 == 92, "got {}".format(c2)
        assert len(dm) == 92, "got {}".format(len(dm))

if __name__ == '__main__':
    unittest.main()


