from .seq_struct import * 
from morebs2.numerical_generator import prg_choose_n

class TrinaryVec(IntSeq):

    def __init__(self,l):
        assert set(l).issubset({0,1,-1})
        super().__init__(l)
        
    # TODO: test these two methods
    @staticmethod
    def one_instance__v1(base_value,length,change_ratio,prg): 
        assert base_value in {-1,0,1}
        assert change_ratio >= 0.0 and change_ratio <= 1.0 

        vec = np.zeros((length,),dtype=np.int32) \
            + base_value 

        n = int(ceil(change_ratio * length))
        Q = [i for i in range(length)]
        ix = prg_choose_n(Q,n,prg,is_unique_picker=True)

        L = {-1,0,1} - set([base_value]) 
        L = sorted(L) 

        for i in ix:
            d = L[prg() % 2]
            vec[i] = d 
        return TrinaryVec(vec)

    @staticmethod
    def one_instance__v2(length,frequency_map,prg): 
        assert 1.0 - sum(frequency_map.values()) < 10 ** -5 
        assert set(frequency_map.keys()).issubset({-1,0,1})
        rx = [(k,int(round(v * length))) for k,v in frequency_map.items()]
        rx = np.array(rx,dtype=np.int32) 
        sx = length - sum(rx[:,1])

        # choose random index to add calibration to 
        if sx != 0:
            rp = prg() % len(rx) 
            rx[rp,1] = rx[rp,1] + sx 

        Q = [i for i in range(length)]
        vec = np.zeros((length,),dtype=np.int32)

        for q in rx:
            k,v = q[0],q[1] 
            q2 = prg_choose_n(Q,v,prg,is_unique_picker=True)
            vec[q2] += k

        return TrinaryVec(vec)

    """
    calculates rightward contiguous range for l[i]. 
    """
    def contiguous_range(self,i):
        q = self.l[i]
        start,end = i,None
        for j in range(i+1,len(self.l)): 
            if self.l[j] != q: 
                end = j
                break 

        if type(end) == type(None):
            end = len(self.l) 
        return [start,end] 
