from .affine_ops import * 
from collections import defaultdict
from morebs2.matrix_methods import is_2dmatrix
from morebs2.measures import normalize_vector
from types import FunctionType,MethodType

"""
container to store <AffineDelta> hypotheses for best fit for an I/O 
dataset (X,Y). One of the following dimensions for input-output must apply 
for every sample in the dataset.  
- 1-to-n,
- n-to-n,
- 1-to-1. 

All samples of X must be of the same dimension to one another, and likewise 
for Y-space. 
"""
class IOAffineFit: 

    def __init__(self,X,Y,prg): 
        assert type(prg) in {MethodType,FunctionType}

        self.X = X 
        self.Y = Y 
        self.input_dim = None 
        self.output_dim = None 
        self.prg = prg 

        assert self.data_check() 

        # these variables are used for constructing 
        # AffineDelta hypotheses 
        self.D = None 
        self.cv = None 
        self.ma_dim = None 
        self.ma_order = None 

        self.i2f_map = dict()  

    def data_check(self): 

        assert is_2dmatrix(self.X) or is_vector(self.X)
        assert is_2dmatrix(self.Y) or is_vector(self.Y)

        xshape = self.X.shape 
        yshape = self.Y.shape 

        if len(xshape) == len(yshape): 
            if xshape != yshape: 
                return False 
        else: 
            if len(self.X) != len(self.Y): 
                return False 
             
        self.input_dim = 0 if is_vector(self.X) else xshape[1] 
        self.output_dim = 0 if is_vector(self.Y) else yshape[1] 

        if self.input_dim > 0: 
            if self.output_dim == 0: 
                return False 

        return True 

    def io_sample(self,i): 
        assert 0 <= i < len(self.X) 
        return (self.X[i],self.Y[i]) 

    """
    """
    def init_hypotheses(self,ma_dim,ma_order,cv): 
        ##assert dim2_to_output_dim(ma_dim) == self.output_dim 
        assert ma_order in {0,1}

        self.i2f_map.clear() 

        if cv != "random": 
            #if type(cv) == tuple:
            assert len(cv) == 2 
            assert min(cv) >= 0.0 and sum(cv) == 1.0 

        self.ma_dim = ma_dim 
        self.ma_order = ma_order 

        for i in range(self.X.shape[0]): 
            ad = self.hypothesis_at_index(i,cv)
            self.i2f_map[i] = ad 

        return  

    def hypothesis_at_index(self,i,cv):

        if cv == "random": 
            x0 = abs(self.prg()) 
            x1 = abs(self.prg()) 
            X = np.array([x0,x1]) 
            cv = normalize_vector(X,5)

        x,y=self.io_sample(i)

        ad = AffineDelta.io_to_AffineDelta(x,y,cv,self.ma_dim,self.ma_order)
        return ad 

    #--------------------------------------------------------------