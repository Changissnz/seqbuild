from mini_dm.vs_fit_op import * 
from math import ceil 
from morebs2.numerical_generator import prg_partition_for_sz__n_rounds

DEFAULT_POINTSET_PARTSIZE_RATIO_RANGE = [0.15,0.23] 

"""
generates points using <AffineDelta> instances, which are in turn 
created through `prg`. Ordering of points can correspond to that of 
the <AffineDelta> sequence, or be unordered. 

Main method is `generate_points`. Main variables are 
`input_seq` and `point_seq`. 
"""
class PointSetGen__TypeAffine:

    def __init__(self,num_points,prg,ro_prg,ro_prg2,ro_prg3=None):
        assert type(prg) in {FunctionType,MethodType}
        assert type(ro_prg) in {FunctionType,MethodType}
        assert type(ro_prg2) in {FunctionType,MethodType}
        assert type(ro_prg3) in {FunctionType,MethodType,type(None)}

        self.num_points = None 
        self.set_num_points(num_points)
        # outputs single integers 
        self.prg = prg 
        # used for var m 
        self.ro_prg = ro_prg
        # used for var a 
        self.ro_prg2 = ro_prg2
        # spare, used for output 
        self.ro_prg3 = ro_prg3

        # output sequence 
        self.point_seq = []
        # input sequence 
        self.input_seq = [] 

        # information pertaining to `point_seq`
        # every element is (affine delta, list of indices)
        self.ad_seq = []
        self.ad_indices_seq = [] 

        self.prt = None 

        # to set parameter values for 
        #       [dim_range,ma_order,input_dim]
        self.gen_ad_parameters = [None,None,None]
        self.set_new_gen_ad_parameters()

    def generate_n_partition(self,n):
        # calculate the number of sets for partition 
        r = modulo_in_range(self.prg(),\
            DEFAULT_POINTSET_PARTSIZE_RATIO_RANGE)        
        rx = int(ceil(n * r))
        var = 0.5 
        num_rounds = 5 
        prt = prg_partition_for_sz__n_rounds(n,rx,self.prg,var,num_rounds)
        return prt 
    
    def set_num_points(self,n): 
        assert type(n) in {int,np.int32,np.int64} 
        assert n > 0 
        self.num_points = n 

    def set_new_gen_ad_parameters(self): 
        dim = int(modulo_in_range(self.prg(),DEFAULT_AFFINEVEC_DIMRANGE))
        dim_range = [dim,dim+1]
        ma_order = int((self.prg() % 2)) 

        if int(self.prg() % 2): 
            input_dim = dim 
        else: 
            input_dim = 0 

        self.gen_ad_parameters = [dim_range,ma_order,input_dim]

    def clear(self): 
        self.point_seq.clear()
        self.ad_seq.clear() 
        self.ad_indices_seq.clear() 
        self.gen_ad_parameters = [None,None,None]

    def one_new_AffineDelta(self):
        return AffineDelta.one_instance_(self.prg,\
            self.ro_prg,self.ro_prg2,dim_range=self.gen_ad_parameters[0],\
            ma_order=self.gen_ad_parameters[1])
    
    def one_pointset(self,size_range):
        assert is_valid_range(size_range,inclusive=False)   
        
        ad = self.one_new_AffineDelta()
        num_points = modulo_in_range(self.prg(), size_range)
        L = []
        L0 = [] 
        for _ in range(num_points):
            x,q = self.one_io_pair(ad) 
            q = self.modulate_point(q) 
            L0.append(x) 
            L.append(q) 
        return ad,np.array(L0),np.array(L) 
    
    def one_io_pair(self,ad):
        x,y = None,None 
        if self.gen_ad_parameters[2] == 0: 
            x = self.prg() 
            y = ad.fit(x) 
        else: 
            x = [] 
            for _ in range(self.gen_ad_parameters[2]): 
                x.append(self.prg()) 
            x = np.array(x) 
            y = ad.fit(x) 
        return x,y 
         
    """
    main method 
    """
    def generate_points(self,is_ordered:bool,clear_data:bool=True):

        if clear_data: 
            self.clear() 
            self.set_new_gen_ad_parameters()

        self.prt = self.generate_n_partition(self.num_points)

        px = None 
        if is_ordered:
            px = self.generate_points_ordered()
        else: 
            px = self.generate_points__unordered()


        return px 

    """
    generates points that are not ordered according to the 
    order of <AffineDelta>s in `ad_seq`. For every point to 
    be generated, uses `prg` to select an <AffineDelta> A, and 
    A is used to calculate the point. 
    """
    def generate_points__unordered(self):

        s = len(self.ad_indices_seq)         
        l = len(self.prt)

        for _ in range(l): 
            self.ad_seq.append(self.one_new_AffineDelta())
        qx = [[] for _ in range(l)]
        self.ad_indices_seq.extend(qx) 

        # 
        index_range = [s,len(self.ad_indices_seq)]

        lx = len(self.point_seq) 
        for i in range(self.num_points): 
            adi = modulo_in_range(int(self.prg()),index_range)
            p = self.ad_seq[adi].fit(self.prg())
            p = self.modulate_point(p) 
            self.point_seq.append(p) 
            self.ad_indices_seq[adi].append(i+lx)
            
        return
    
    """
    generates points that are ordered according to the 
    order of <AffineDelta>s in `ad_seq`. Ordering of 
    <AffineDelta>s (in `ad_seq`) is spatially correlated 
    to the ordering of elements in `input_seq` and `point_seq`.
    """
    def generate_points_ordered(self):
        q = 0 
        lx = len(self.point_seq) 
        for p in self.prt:
            ad,L0,L = self.one_pointset([p,p+1]) 
            self.ad_seq.append(ad) 
            self.input_seq.extend(L0) 
            self.point_seq.extend(L) 
            self.ad_indices_seq.append([lx + x for x in range(q,q+p)]) 
            q = q + p 

    def modulate_point(self,p): 
        if type(self.ro_prg3) == type(None):
            r1 = self.ro_prg() 
            r2 = self.ro_prg2() 
            rx = ((r1[0] + r2[0]) / 2.0, (r1[1] + r2[1]) / 2.0)
        else: 
            rx = self.ro_prg3()
        return modulo_in_range(p,rx) 

        
