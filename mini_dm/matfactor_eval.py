from collections import defaultdict 
import numpy as np 

"""
In a self-factorable list, every element is a multiple of all 
elements less than it. 
"""
def is_self_factorable_list(l,output_bool:bool=False):

    l = np.array(l) 

    not_factorable = []
    for (i,l_) in enumerate(l): 
        xs = np.where(l < l_)[0]
        if len(xs) == 0: continue 

        q = l[xs] 
        if not np.all(l_ / q == np.asarray(l_ / q,dtype=np.int32)): 
            if output_bool: return False 
            not_factorable.append(i)

    if output_bool: return True 
    return not_factorable

"""
A column-factor map is 

row index -> <column set>,<factor for column set>. 

Outputs a sequence with elements 
(column set, factor for column set). 
"""
def columnfactormap_to_sequence(cfd): 
    s = [] 
    v_ = list(cfd.values())
    for v2_ in v_: 
        q = list(v2_.values())
        for q_ in q: 
            s1,s2 = q_[0],q_[1] 
            for s1_,s2_ in zip(s1,s2):  
                s.append((s1_,s2_))

    return s 

"""
The factor rank of a column is the summation 
of the factors corresponding to the column sets 
that the column is in. 

seq := list, sequence with elements 
(column set, factor for column set).
"""
def columnfactor_rank(seq,rank_type:int=2): 
    assert rank_type in {1,2}
    d = defaultdict(float) 

    for x0,x1 in seq: 
        for x in x0: 
            if rank_type == 1: 
                d[x] += 1 
            else: 
                d[x] += x1 
    return d 

"""
Converts a sequence with elements 
    (column set, factor for column set)
into a sequence of 
    (column set, <factors for column set>). 
"""
def columnfactorpair_to_seq(seq): 

    def indexia(q):
        for (i,s) in enumerate(sts2): 
            if s[0] == q: return i 
        return -1 


    # collect all sets first 
    sts = [] 
    for (v1,_) in seq: 
        if v1 not in sts: 
            sts.append(v1)
    
    sts2 = [] 
    for v1,v2 in seq: 
        i = indexia(v1) 
        if i == -1: 
            sts2.append((v1,[v2]))
        else: 
            sts2[i][1].append(v2) 
    return sts2 

"""
Determines each column set of a sequence with elements 
(column set, factor for column set) 
that has identity columns in it based on the factor count 
w.r.t. the number of rows belonging to the source matrix. 
"""
def columnfactor_identity(seq,num_rows): 
    sts2 = columnfactorpair_to_seq(seq) 
    q = [] 
    nr = np.prod([i for i in range(1,num_rows+1)])
    for (r1,r2) in sts2:  
        if len(r2) != nr: continue 
        b = is_self_factorable_list(r2,output_bool=True)
        if b: 
            q.append(r1)
    return q 

def indexvalue_map_to_vector(m,sz): 
    l = np.zeros((sz,))
    for k,v in m.items(): 
        l[k] = v 
    return l 

def vector_to_indexvalue_map(v): 
    d = {}
    for (i,v_) in enumerate(v): 
        d[i] = v_
    return d 

"""
Evaluates factors between rows of a matrix for 
column-wise identities. 
"""
class MatFactorEval: 

    def __init__(self,M): 
        assert type(M) == np.ndarray 
        assert M.ndim == 2 
        assert M.shape[0] > 1 and M.shape[1] > 1 
        self.M = M 
        self.id_col = None 
        self.cfm = None 

    #----------------- preliminary: compares literal matrix values 
    def identity_columns_(self): 
        col = [] 
        for c in range(self.M.shape[1]): 
            if self.is_identity_column(c):
                col.append(c) 
        return col 

    def is_identity_column(self,c):  
        return np.all(self.M[0,c] == self.M[:,c])

    #----------------- factor evaluation: comparisons with application of 
    #                  multiples. 

    def identityvar_in_two_rows(self,r1,r2): 
        x = r1 / r2
        q = []
        r = [] 
        for i in range(len(x) -1): 
            for j in range(i + 1,len(x)): 
                if x[j] == x[i]: 
                    stat = False 
                    for q_ in q: 
                        if i in q_: 
                            q_ |= {j} 
                            stat = True 
                            break  
                    if not stat: 
                        q.append({i,j}) 
                        if x[j] < 1: 
                            r.append(r2[i] / r1[i])
                        else: 
                            r.append(x[j])
        return q,r 

    def identity_var_for_row(self,r): 
        d = {} 
        row = self.M[r] 

        for i in range(self.M.shape[0]): 
            if i == r: continue 
            iv = self.identityvar_in_two_rows(row,self.M[i]) 
            d[i] = iv  
        return d 

    def column_factor_map(self): 
        d = {} 
        for i in range(self.M.shape[0]): 
            # ineff: recount 
            q = self.identity_var_for_row(i)
            d[i] = q 
        self.cfm = d 

    def identity_eval(self,rank_type=1): 
        self.id_col = self.identity_columns_() 
        self.column_factor_map() 
        seq = columnfactormap_to_sequence(self.cfm)
        q1 = columnfactor_rank(seq,rank_type)
        q2 = columnfactor_identity(seq,self.M.shape[0]) 
        return q1,q2 

class MatrixConsistencyCheck:

    def __init__(self,M,Y): 
        assert type(M) == np.ndarray 
        assert M.ndim == 2 
        assert M.shape[0] >= 1 and M.shape[1] >= 1 
        assert type(Y) == np.ndarray 
        assert Y.ndim == 1 and len(Y) == M.shape[0] 
        self.M = M 
        self.Y = Y
        self.inconsistent = None 

    def check(self): 
        self.inconsistent = [] 
        for r in range(self.M.shape[0]-1): 
            for r2 in range(r + 1, self.M.shape[0]): 
                stat = self.check_row_pairs(r,r2)
                if not stat: 
                    self.inconsistent.append((r,r2))
        return self.inconsistent

    def check_row_pairs(self,r1,r2): 
        x = self.M[r1] / self.M[r2] 
        if np.any(np.isinf(x)): return True

        xu = np.unique(x)
        if len(xu) > 2: return True 

        if len(xu) == 2: 
            if not np.any(np.isnan(xu)): return True 

        xu = set(xu)
        n = xu.pop() 
        if np.isnan(n): 
            if len(xu) == 0: return True 
            n = xu.pop() 
        q = round(self.Y[r2] * n - self.Y[r1],5) 
        return q == 0.0 
        
