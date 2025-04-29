import numpy as np

def diffvec(l,cast_type=np.int32):
    assert cast_type in [np.int32,np.float32] 
    assert type(l) == np.ndarray and len(l.shape) == 1
    d = []
    for i in range(1,len(l)):
        d.append(l[i] - l[i-1])
    return np.array(d,dtype=cast_type) 

def gleqvec(l,rounding_depth=5): 
    assert type(l) == np.ndarray and len(l.shape) == 1
    d = []
    for i in range(1,len(l)):
        l1,l2 = round(l[i-1],rounding_depth), round(l[i],rounding_depth)
        equals = l1 == l2 
        if equals: 
            d.append(0)
            continue 
        d.append(-1) if l1 > l2 else d.append(1)
    return np.array(d,dtype='int32')

def affine_fit_on_vec():

    return -1  


class IntSeq:

    def __init__(self,l):
        self.l = l 
        self.load()

    def load(self):
        if type(self.l) == type(None): 
            self.l = np.array([],"int32")
        else: 
            self.l = np.array(self.l,"int32")
            assert len(self.l.shape) == 1 
        return 

    def append(self,q): 
        if type(q) in [int,np.int32]: 
            q = np.int32(q)
            self.l = np.append(self.l,q)

        for q_ in q: 
            self.l = np.append(self.l,np.int32(q_))
        return

    '''
    difference triangle 
    '''
    def difftri(self):
        # get the first 
        dvec = diffvec(self.l) 
        l = len(dvec)
        if l == 0: 
            print("[!] NONE.")
            return

        x = np.zeros((l,l),dtype='int32') 
        x[0] = np.copy(dvec)
        lx = l - 1 

        while lx > 0: 
            dvec = diffvec(dvec) 
            dx = l - len(dvec) 
            nx = np.zeros((dx,))
            dvec2 = np.append(nx,dvec)
            x[dx] = np.copy(dvec2) 
            lx -= 1
        return x 