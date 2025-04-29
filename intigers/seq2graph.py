from .process_seq import * 

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

    def optri(self,operation,cast_type):

        # get the first 
        dvec = stdop_vec(self.l,operation,cast_type) 
        l = len(dvec)
        if l == 0: 
            print("[!] NONE.")
            return

        x = np.zeros((l,l),dtype='int32') 
        x[0] = np.copy(dvec)
        lx = l - 1 

        while lx > 0: 
            dvec = stdop_vec(dvec,operation,cast_type) 
            dx = l - len(dvec) 
            nx = np.zeros((dx,))
            dvec2 = np.append(nx,dvec)
            x[dx] = np.copy(dvec2) 
            lx -= 1
        return x 

    '''
    difference triangle 
    '''
    def difftri(self,cast_type=np.int32): 
        return self.optri(sub,cast_type)