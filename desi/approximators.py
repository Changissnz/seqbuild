from morebs2.fit_2n2 import * 
from intigers.extraneous import safe_div 
from .fraction import * 

DEFAULT_FIT22_PARTITION_SIZE_RANGE = [5,103] 

class Fit22ValueOutputter(GenericIntSeqOp): 

    def __init__(self,px,length_outputter,range_outputter,bool_outputter,\
            point_conn_type:int,adjustment_type:int=1):  
        assert type(px) in {MethodType,FunctionType}
        assert point_conn_type in {1,2}
        assert type(bool_outputter) in {MethodType,FunctionType}

        self.px = px 
        self.point_conn_type = point_conn_type
        self.b_out = bool_outputter

        self.p0 = None
        self.p1 = None
        self.f = None

        self.x_ = None
        self.x2_ = None
        self.a_ = None
        self.ax = None 
        
        super().__init__(length_outputter,range_outputter,adjustment_type)
        return
    
    def set_next_value(self):
        if type(self.x_) == type(None):
            self.update_points() 
        elif self.x_ > self.x2_: 
            self.update_points() 
        else: pass 

        if self.ax: 
            sx = self.f.f(self.x_)
        else: 
            sx = self.f.g(self.x_) 
        self.x_ += self.a_ 

        sx_ = float_to_string(sx,True,True)
        self.s = sx_ 
        return 
    
    def update_points(self):
        # case: start 
        if type(self.p0) == type(None):
            self.p0 = self.px() 
            self.p1 = self.px()
            return 
        
        # output 1 or 2 new points 
        if self.point_conn_type == 1:
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

        if round(self.p0[0] - self.p1[0],5) != 0.0 and \
            round(self.p0[1] - self.p1[1],5) != 0.0:
            return 
        
        p = self.px() 
        qr = safe_div(self.p1,p)

        while round(self.p0[0] - self.p1[0],5) == 0.0 or \
            round(self.p0[1] - self.p1[1],5) == 0.0:
            self.p1 = self.p1 + qr 

"""
LPS stands for Lagrange Polynomial Solver, interpolation of a 
polynomial's span of points using the Lagrange basis.  
"""
class LPSValueOutputter(GenericIntSeqOp):

    def __init__(self): 
        return 