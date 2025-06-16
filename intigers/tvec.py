

class TrinaryVec:

    def __init__(self,l):
        assert set(l).issubset({0,1,-1})
        self.l = l  