from face.cl_guide_parser import *
import unittest

### lone file test 
"""
python -m tests.test_cl_guide_parser
"""
###
class CLGuideParserMethods(unittest.TestCase):

    def test__CLGuideParser__process_structures(self):
        cgp = CLGuideParser() 
        cgp.process_structures()
        assert len(cgp.structures) >= 10

        cgp.close()

    def test__CLGuideParser__process_keywords(self): 
        cgp = CLGuideParser() 
        cgp.process_structures()
        cgp.process_keywords() 

        keywords_sol = {\
            0: {'make', 'set', 'open', 'write', 'run', 'convert'}, \
            1: {'for', 'to', 'with', 'iter'}, \
            2: {'open': ['file', 'seq', 'obj'], \
                'convert': ['range', 'ndim', 'nvec', 'tvec']}}
        assert cgp.keywords == keywords_sol
        cgp.close()

    def test__CLGuideParser__process_command_forms(self): 
        cgp = CLGuideParser() 
        cgp.process_structures()
        cgp.process_keywords() 
        cgp.process_command_forms()
        
        keys = cgp.command_forms.keys()
        sol_keys = set(['make', 'run', 'set', 'open', 'write', 'convert', 'show']) 
        assert sol_keys.issubset(keys)

        cgp.close() 

if __name__ == '__main__':
    unittest.main()