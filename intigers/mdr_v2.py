from .seq_struct import * 
from .extraneous import safe_div 

def afs_prt__size(afs_prt):
    c = 0 
    for x in afs_prt:  
        c += len(x[1])
    return c 

def afs_prt__density(afs_prt): 
    if len(afs_prt) == 0: return 0.0 

    start_index = afs_prt[0][0][0]
    end_index = afs_prt[-1][0][1] 

    # get active ranges 
    c = afs_prt__size(afs_prt) 
    return safe_div(end_index-start_index,c)

"""
transforms a list `afs_prt`, representing a partition 
from an affine-fit search over an integer sequence, from 
a feed-backward to a feed-forward indexing. There is one 
key difference between these two indexing schemes. 

Consider the format for one element in `afs_prt`: 

[0] (start range, end range) 
[1] [[(m_0, a_0), [i_00, i_01]],...,[(m_j,a_j),[i_j0,i_j1]]; 
    i_00 := start range, 
    i_j1 := end range. 

In a feed-backward indexing, the output value at the `i_x0`'th index 
uses as input the (`i_x0` - 1)'th output value. 

In a feed-forward indexing, the output value at the (`i_x0` + 1)'th index 
uses as input the `i_x0`'th output value. 
"""
def afs_prt__feedbw2feedfw(afs_prt):
    afs_prt_ = deepcopy(afs_prt)
    for i in range(len(afs_prt_)):
        q = afs_prt_[i][1] 
        for j,q_ in enumerate(q): 
            q2_ = q_
            q2_[1][0] -= 1
            q2_[1][1] -= 1 
            q[j] = q2_ 

    return afs_prt_ 

"""
version is an improvement over the original <ModuloDecomp>, 
which suffered from a deficit due to constraints of feed-backward 
processing, a trait indicated by the indexing style of its 
variable `afs_prt`. 

When given two contiguous subsequences from an `afs_prt`, 
S0 = S[i] and S1 = S[i+1], S0 and S1 subsequences of `l`,
and the last integer v of S0 and first integer v2 of S1 result 
in 
    difference = (v * m + a) - v2 such that 
    |difference| < v2, 
the deficit of feed-backward processing does not permit v 
to be calculated into v2 via its (m,a) and modulus (of `afs_prt_mod`), 
the modulus being `difference`. 

<ModuloDecompV2> uses a feed-forward indexing scheme to overcome 
this deficit of its superclass. 
"""
class ModuloDecompV2(ModuloDecomp):

    def __init__(self,l,exclude_neg:bool=True):
        super().__init__(l)  
        self.exclude_neg = exclude_neg 
        self.merge(exclude_neg) 
        self.afs_prt = afs_prt__feedbw2feedfw(self.afs_prt)

        for i in range(len(self.afs_prt) - 1):
            self.fix_subend(i)

    def size(self):
        return afs_prt__size(self.afs_prt)

    def density(self): 
        return  afs_prt__density(self.afs_prt)

    def check_subend(self,ap_index): 
        ix = self.gleqvec_prt[ap_index] 
        qseq = self.afs_prt[ap_index][1] 
        m,a = qseq[-1][0]

        prev,now = self.l.l[ix],self.l.l[ix+1] 
        x = m * prev + a 

        mx = self.afs_prt_mod[ap_index]
        if x % mx == now: 
            return True 

        diff = x - now
        return not diff < abs(now)

    def fix_subend(self,ap_index): 
        stat = self.check_subend(ap_index) 
        qseq = self.afs_prt[ap_index][1] 
        if stat: 
            x = qseq[-1]
            x[1][1] += 1
            return 
        
        ix = self.gleqvec_prt[ap_index] 
        qseq = self.afs_prt[ap_index][1] 
        m,a = qseq[-1][0]
        prev,now = self.l.l[ix],self.l.l[ix+1] 

        MA = [m,a]

        def add_one(index): 
            MA[index] += 1 

        x = m * prev + a 
        diff = x - now
        index = 0 if prev != 0 else 1 
        while abs(diff) < abs(now): 
            add_one(index) 
            got = prev * MA[0] + MA[1] 
            diff = got - now 
        qseq.append([tuple(MA),[ix,ix]])

        got = MA[0] * prev + MA[1] 
        diff = got - now
        self.afs_prt_mod[ap_index] = diff 