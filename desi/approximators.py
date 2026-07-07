from .not_equals import * 
from .fraction import * 
from morebs2.poly_interpolation import * 
from morebs2.numerical_generator import prg__constant

DEFAULT_FIT22_PARTITION_SIZE_RANGE = (5,103) 

DEFAULT_POLYNOMIAL_PARTITION_SIZE_RANGE = (50,1000)
DEFAULT_POLYNOMIAL_POINT_SIZE_RANGE = (5,21) 

def default_length_outputter(prg_2vec): 
    x = list(prg_2vec()) + list(prg_2vec()) 
    q = prg__LCG(x[0],x[1],x[2],x[3])
    return prg__single_to_int(q) 

def default_range_outputter(prg_2vec): 
    x = list(prg_2vec()) + list(prg_2vec()) 
    q = prg__LCG(x[0],x[1],x[2],x[3])
    return prg__single_to_range_outputter(q)

def default_bool_outputter(prg_2vec): 
    return wrap_ranged_modulo_over_generator(\
        default_length_outputter(prg_2vec),[0,2])

def prg__nvec_to_int_type(px): 

    def f(): 
        return np.array(px(),dtype=int)
    return f  

class Fit22ValueOutputter(GenericIntSeqOp): 

    def __init__(self,px,length_outputter,range_outputter,bool_outputter,\
            point_conn_type:int,adjustment_type:int=1):  
        assert type(px) in {MethodType,FunctionType}
        assert type(length_outputter) in {type(None),MethodType,FunctionType}
        assert type(range_outputter) in {type(None),MethodType,FunctionType}
        assert type(bool_outputter) in {type(None),MethodType,FunctionType}

        assert point_conn_type in {1,2}

        self.px = px  
        self.point_conn_type = point_conn_type

        self.p0 = None
        self.p1 = None
        self.H = None

        self.x_ = None
        self.x2_ = None
        self.a_ = None
        self.ax = None 
        
        if type(length_outputter) != type(None): 
            length_outputter = prg__single_to_int(length_outputter)
        else: 
            length_outputter = default_length_outputter(px) 

        if type(range_outputter) == type(None): 
            range_outputter = default_range_outputter(px)  

        if type(bool_outputter) == type(None): 
            bool_outputter = default_bool_outputter(px) 
        else: 
            bool_outputter = wrap_ranged_modulo_over_generator(\
                prg__single_to_int(bool_outputter),[0,2]) 

        self.b_out = bool_outputter 
        super().__init__(length_outputter,range_outputter,adjustment_type)
        return
    
    def set_next_value(self):
        if type(self.x_) == type(None):
            self.update_points() 
        elif self.x_ > self.x2_: 
            self.update_points() 
        else: pass 

        # do check on input before function apply 
        A = (self.ax + 1) % 2 
        qmin = np.min(self.H.ps[:,A]) 
        qmax = np.max(self.H.ps[:,A]) 
        x_ = modulo_in_range(self.x_,[qmin,qmax]) 

        QF = self.H.f if self.ax else self.H.g 
        try: 
            sx = QF(self.x_)
        except: 
            sx = QF(x_)

        self.x_ += self.a_ 

        sx_ = float_to_string(sx,True,True)
        self.s = sx_ 
        return 
    
    def update_points(self):
        # case: start 
        if type(self.p0) == type(None):
            self.p0 = self.px() 
            self.p1 = self.px() 
        # output 1 or 2 new points 
        elif self.point_conn_type == 1:
            self.p0 = self.p1
            self.p1 = self.px()
        else:
            self.p0,self.p1 = self.px(),self.px() 

        self.fix_equal_points()

        # instantiate either <LogFit22> or <Exp2Fit22> 
        mx = np.array([self.p0,self.p1]) if \
            self.p0[1] < self.p1[1] else \
            np.array([self.p1,self.p0])
        lORe = self.b_out() % 2

        if lORe:
            self.H = LogFit22(mx,direction=[0,1]) 
        else: 
            self.H = Exp2Fit22(mx,direction=[0,1])

        # set delta variables for exploring the fitter 
        self.ax = self.b_out() % 2 

        ax_ = (self.ax + 1) % 2        
        self.x_ = min(self.H.ps[:,ax_])
        self.x2_ = max(self.H.ps[:,ax_]) 

        ap = modulo_in_range(self.l_out(),DEFAULT_FIT22_PARTITION_SIZE_RANGE) 
        self.a_ = (self.x2_ - self.x_) / ap

    
    def fix_equal_points(self):

        self.p0 = np.array(self.p0)
        self.p1 = np.array(self.p1)

        self.p0,self.p1 = \
            not_equals__pairvec(self.p0,self.p1,self.px,indices=None)



"""
LPS stands for Lagrange Polynomial Solver, interpolation of a 
polynomial's span of points using the Lagrange basis.  
"""
class LPSValueOutputter(GenericIntSeqOp):

    def __init__(self,px,length_outputter,range_outputter,\
        adjustment_type:int=1): 
        
        assert type(px) in {MethodType,FunctionType}
        assert type(length_outputter) in {type(None),MethodType,FunctionType}
        assert type(range_outputter) in {type(None),MethodType,FunctionType}
        self.px = px 

        self.p0,self.p1 = None,None 
        self.adder_end,self.adder,self.adder_i = None,None,None

        if type(length_outputter) != type(None): 
            length_outputter = prg__single_to_int(length_outputter)
        else: 
            length_outputter = default_length_outputter(px) 

        if type(range_outputter) == type(None): 
            range_outputter = default_range_outputter(px)  

        super().__init__(length_outputter,range_outputter,adjustment_type)

    #--------------------------------------------------------------------------

    def set_next_value(self):
        if type(self.p0) == type(None):
            self.update_points() 
        elif self.adder_i > self.adder_end: 
            self.update_points() 
        else: pass 

        x = self.p0 + (self.adder_i * self.adder)
        sx = self.lps.output_by_lagrange_basis(x)
        sx_ = float_to_string(sx,True,True)
        self.s = sx_ 

        self.adder_i += 1 

        return 
    
    def update_points(self):

        self.q = LPSValueOutputter.new_pointset(self.px,self.l_out) 

        # assumes q is a proper set of points for 
        # fitting via <LagrangePolySolver>
        self.lps = LagrangePolySolver(self.q, prefetch = True)
        self.p0,self.p1 = self.lps.minumum,self.lps.maximum 

        self.adder_end = modulo_in_range(self.l_out(),DEFAULT_POLYNOMIAL_PARTITION_SIZE_RANGE)
        self.adder = (self.p1 - self.p0) / self.adder_end 
        self.adder_i = 0   

    """

    """
    @staticmethod 
    def new_pointset(px,l_out):

        num_points = modulo_in_range(\
            l_out(),DEFAULT_POLYNOMIAL_POINT_SIZE_RANGE) 

        M = []
        for _ in range(num_points):
            x = px() 
            assert type(x) in {list,tuple} or is_vector(x) 
            assert len(x) == 2 
            M.append(x) 

        # NOTE: this is additional code. Guarantees matrix 
        #       contains only unique numbers. 
        M = np.round(np.array(M),5) 
        q = np.array(M.flatten()  * 10 ** 5,dtype=int) 

        x = prg__constant(2.5)
        q = intlist_no_dups_no_zero_abs__maintain_size(q,x,False) 
        q = np.round(np.array(q) / (10 ** 5),5) 
        M = q.reshape((num_points,2)) 

        submat_type = l_out() % 4 
        t = DEFAULT_SUBMAT_TYPES[submat_type]

        I = np.argsort(M[:,1])
        M = M[I] 
        
        # NOTE: first method may yield crashes (equal values retained)
        return not_equals__matrix_whole(M,l_out,t)

