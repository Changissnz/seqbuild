from .approximators import * 

class PointInterpolationContainer:

    def __init__(self,pts,starting_fvec):
        assert is_2dmatrix(pts)
        pts = pts[np.argsort(pts[:,0])]
        self.pts = pts
        self.fvec = starting_fvec
        self.declare_structs() 

    def declare_structs(self):
        self.lps = LagrangePolySolver(self.pts, prefetch = False)
        self.cdc = ChainedDCurves(self.pts,self.fvec,0)
        self.cdc.draw() 
        self.x = self.cdc.amin 
        self.x2 = self.cdc.amax 

    def differential_at(self,x,lps_first:bool=True):
        q = self.lps.output_by_lagrange_basis(x) 
        q2 = self.cdc.fit(x) 

        if lps_first: return q - q2 
        return q2 - q 
    
    def modulate_dcurve_of_point(self,x): 
        index = self.cdc.point_to_dcurve(x)
        assert index != -1
        self.cdc.dcurve_seq[index].modulate_fit() 


"""
PID stands for Point Interpolation Differential. 
"""
class PIDValueOutputter(GenericIntSeqOp):

    def __init__(self,px,frequency_outputter,\
        length_outputter,range_outputter,adjustment_type):
        assert type(px) in {MethodType,FunctionType}
        
        assert type(frequency_outputter) in {MethodType,FunctionType}

        self.px = px
        self.pts = None

        self.f_out = frequency_outputter
        super().__init__(length_outputter,range_outputter,adjustment_type)

        self.adder_end,self.adder,self.adder_i = -1,None,0
        self.pic = None 
        return 
    
    def set_next_value(self):

        if type(self.pic) == type(None) or self.adder_i > self.adder_end:
            self.pts = LPSValueOutputter.new_pointset(self.px,self.l_out)
            fvec = [self.f_out() % 2 for  _ in range(len(self.pts) - 1)] 
            self.pic = PointInterpolationContainer(self.pts,fvec) 

            self.adder_end = modulo_in_range(self.f_out(),\
                DEFAULT_POLYNOMIAL_PARTITION_SIZE_RANGE)
            self.adder = (self.pic.x2 - self.pic.x) / self.adder_end
            self.adder_i = 0
            self.fvec_frequency = modulo_in_range(\
                self.f_out(),[1,max([2,ceil(self.adder_end / 5)])])
        
        x =  self.pic.x + (self.adder_i * self.adder)
        stat = (self.adder_i % self.fvec_frequency) == 0

        if stat:
            self.pic.modulate_dcurve_of_point(x) 

        q = self.pic.differential_at(x,lps_first=True)
        self.adder_i += 1 

        q_ = int(q) 
        q -= q_

        self.s += str(abs(q_)) 
        self.s += float_to_string(abs(q),True,True)
        return