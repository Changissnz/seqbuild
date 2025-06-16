

class TrinaryVec:

    def __init__(self,l):
        assert set(l).issubset({0,1,-1})
        self.l = l  

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
