from .approximators import * 

class PointInterpolationContainer:

    def __init__(self,pts,starting_fvec):
        assert is_2dmatrix(pts)
        pts = pts[np.argsort(pts[:,0])]
        self.pts = pts
        self.fvec = starting_fvec
        self.declare_structs() 

    def declare_structs(self):
        self.lps = LagrangePolySolver(self.pts, prefetch = True)
        self.cdc = ChainedDCurves(self.pts,self.fvec,0)
        self.x = self.cdc.amin 
        self.x2 = self.cdc.amax 

    def differential_at(self,x,lps_first:bool=True):
        q = self.lps.output_by_lagrange_basis(x) 
        q2 = self.cdc.fit(x) 

        if lps_first: return q - q2 
        return q2 - q 

"""
PID stands for Point Interpolation Differential. 
"""
class PIDValueOutputter(GenericIntSeqOp):

    def __init__(self,px,pts,frequency_outputter,\
        length_outputter,range_outputter,adjustment_type):
        assert type(px) in {MethodType,FunctionType}
        
        assert type(frequency_outputter) in {MethodType,FunctionType}

        self.px = px
        self.pts = np.array(self.px())
        assert is_2dmatrix(self.pts)
        assert self.pts.shape[1] == 2

        self.f_out = frequency_outputter
        super().__init__(length_outputter,range_outputter,adjustment_type)

        self.adder_end,self.adder,self.adder_i = None,None,None
        return 
    
    def __next__(self): 

        if type(self.adder_i) == type(None): 
            self.set_next_value()
            return self.__next__() 
        
        if self.adder_i > self.adder_end: 
            self.set_next_value() 
            return self.__next__()
        
        x =  self.pic.x + (self.adder_i * self.adder)

        self.adder_i += 1 
        return x 

    def set_next_value(self):
        starting_fvec = [self.f_out() % 2 for  _ in range(len(self.pts) - 1)]

        if type(self.pic) != type(None): 
            q = np.array(self.px())
            self.pts = q 
            self.pic = PointInterpolationContainer(self.pts,starting_fvec) 
        else: 
            self.pic = PointInterpolationContainer(self.pts,starting_fvec)

        self.adder_end = modulo_in_range(self.f_out(),DEFAULT_POLYNOMIAL_PARTITION_SIZE_RANGE) 
        self.adder = (self.pic.x2 - self.pic.x) / self.adder_end
        self.adder_i = 0 

        # update after 
        self.pts = np.array(self.px())