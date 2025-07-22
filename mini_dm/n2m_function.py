from .n2m_index import * 
from .iseq import modulated_vec_op,modulated_vecdot 
from morebs2.matrix_methods import is_number 
from copy import deepcopy 

"""
class that operates an n-to-m vector function to map 
any n-vector to an m-vector, both of real numbers. 
"""
class N2MVectorFunction:

    def __init__(self,nm,indexset_pair_seq,function_map,mode="replace"):
        assert_nm(nm)  
        self.nm = nm 

        self.check_args(indexset_pair_seq,function_map)
        assert mode in {"replace","accumulate"}


        self.indexset_pair_seq = indexset_pair_seq
        self.function_map = function_map
        self.mode = mode 
        return 

    def check_args(self,ip_seq,f_map): 
        assert len(ip_seq) > 0
        assert type(ip_seq) == list 
        assert len(f_map) > 0 
        assert type(f_map) == dict 
        # check the index set pairs 
        for x in ip_seq: 
            p0 = string_to_vector(x[0]) 
            assert min(p0) >= 0 and max(p0) < self.nm[0]

            p1 = string_to_vector(x[1])
            assert min(p1) >= 0 and max(p1) < self.nm[1]
        
        # check the function map
        qx = list(f_map.keys())
        qx = set([type(k) in {MethodType,FunctionType} for k in qx]) 
        assert qx == {True}

        qx2 = list(f_map.values()) 
        qx2_ = []
        for qx22 in qx2: qx2_.extend(qx22) 
        qx2_ = sorted(qx2_)
        qx2 = np.array(qx2_,dtype=np.int32)
        assert np.all(qx2 == np.arange(len(ip_seq),dtype=np.int32)) 
        return

    """
    main method 
    """
    def apply(self,x,default_value=0.0): 
        assert is_vector(x) 
        assert len(x) == self.nm[0]
        assert is_number(default_value)

        mvec = np.ones((self.nm[1],),dtype=np.float64)
        mvec = mvec * default_value
        for i in range(len(self.indexset_pair_seq)):
            mvec = self.apply_function(x,mvec,i)
        return mvec 

    def apply_function(self,nvec,mvec,ip_index:int): 
        q = self.indexset_pair_seq[ip_index] 
        indices0 = string_to_vector(q[0])
        indices1 = string_to_vector(q[1])

        vx = deepcopy(nvec[indices0])

        f = self.ip_index_to_function(ip_index)
        assert type(f) != type(None)

        vx2 = f(vx)

        if self.mode == "replace":
            mvec[indices1] = vx2
        else: 
            mvec[indices1] += vx2 
        return mvec 

    def ip_index_to_function(self,ip_index:int):
        for k,v in self.function_map.items():
            if ip_index in v: 
                return k 
        return None