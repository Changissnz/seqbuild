from .n2m import * 
from .vs_fit_op import * 
from intigers.mod_prng import * 
from morebs2.numerical_generator import prg__constant

DEFAULT_VFO_AFFINEHYP_POINTSIZE_RANGE = [10,26]
DEFAULT_DX_MULTIPLE_RANGE = [1.0,4]

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
    
    def preproc_v2(self,ma_dim,ma_order): 
        cv = (0.5,0.5) 
        dx = self.default_dx_() 
        self.load_mahyp_auxvar(dx,cv,ma_dim,ma_order) 
        self.init_HypMach() 
    
    def default_dx_(self):
        l = len(self.x)
        dx_vec = [] 
        for i in range(l): 
            x,y = self.io_sample(i)
            dx = y - x 

            if is_number(dx):
                r = [int(DEFAULT_DX_MULTIPLE_RANGE[0] * dx),\
                     int(DEFAULT_DX_MULTIPLE_RANGE[1] * dx) + 1] 
                x = modulo_in_range(self.prg(),r) 
                dx_vec.append(x) 
                continue

            r0 = np.array(DEFAULT_DX_MULTIPLE_RANGE[0] * dx,\
                dtype=int) 
            r1 = np.array(DEFAULT_DX_MULTIPLE_RANGE[1] * dx + 1,\
                dtype= int)
            r = np.array([r0,r1]).T 
            dx = [modulo_in_range(self.prg(),r_) for r_ in r] 
            dx_vec.append(dx) 
        return dx_vec 