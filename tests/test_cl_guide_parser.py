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

if __name__ == '__main__':
    unittest.main()