from .n2m import * 
from .vs_fit_op import * 
from intigers.mod_prng import * 
from morebs2.numerical_generator import prg__constant

DEFAULT_VFO_AFFINEHYP_POINTSIZE_RANGE = [10,26]

class VFOSearch(IOFit):

    def __init__(self,x,y,unknown_func,hypdiff_func,madiff_func,\
        prg):
        super().__init__(x,y,unknown_func,hypdiff_func,madiff_func)
        assert type(prg) in {FunctionType,MethodType}
        self.prg = prg 
        self.n2mac = None
        self.init_n2m_ac() 

    def init_n2m_ac(self):
        q0,q1 = self.type()

        assert len(q0) == 1
        assert len(q1) == 1

        k0 = q0.pop() 
        k1 = q1.pop() 
        self.n2mac = N2MAutocorrelator((k0,k1),False)
        return
    
    """
    hypotheses of <AffineDelta> instances 
    """
    def preproc(self,ma_order=0): 
        # calculating initial affine hypotheses 
        l = len(self.x) 

        if l < DEFAULT_VFO_AFFINEHYP_POINTSIZE_RANGE[1]: 
            self.init_HypMach_v2(ma_order,None,\
                None)
            return 
        
        q = modulo_in_range(self.prg(),\
            DEFAULT_VFO_AFFINEHYP_POINTSIZE_RANGE)
        q = min([q,l // 2])
        assert q > 0 

        index0_vec = []
        index1_vec = []
        indices = [i for i in range(l)] 
        
        index0_vec = prg_choose_n(indices,q,self.prg,is_unique_picker=True) 
        index1_vec = prg_choose_n(indices,q,self.prg,is_unique_picker=True) 

        cx = prg__constant(index1_vec)  
        def index1_function(x): 
            return cx() 

        self.init_HypMach_v2(ma_order,index0_vec,\
            index1_function)
        
        return 