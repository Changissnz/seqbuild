from intigers.mdr_v2 import * 
from intigers.tvec import * 
from intigers.intfactor import * 
from mini_dm.iseq import * 
from intigers.poly_output_fitter_ import * 

DEFAULT_SHADOW_FITTERS = {"mdr","mdrv2","tvec","fvec","optri","pofv1","pofv2"}

class ShadowFitter:

    def __init__(self,fitting_struct,prg): 
        assert fitting_struct in DEFAULT_SHADOW_FITTERS
        self.fitting_struct = fitting_struct 
        return 

    def fit_sequence(self,S): 
        if self.fitting_struct == "mdr": 

        return 

    def apply_delta(self):
        return 

    def refit_sequence(self):
        return 

class ShadowGen: 

    def __init__(self,fitting_struct):
        return