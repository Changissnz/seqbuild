from .lcg_v2 import * 
from .tvec import * 
from .extraneous import prg__single_to_int
from morebs2.numerical_generator import prg_decimal

"""
"""
def is_value_below_modulo_range_length(v,mr):
    assert is_valid_range(mr,True,False) or is_valid_range(mr,False,False)
    return abs(v) <= mr[1] - mr[0]

# NOTE: not the most efficient approach
"""
`set_limit` sets limits for the search (respectively for integer and float 
types). Variable is either int|float and guarantees that scanning algorithm 
does not become an infinite loop. 
"""
def iterative_adjustment_for_direction__modrange(pv,av,sign,modrange,\
    modindex,moddir,set_limit=int):  
    assert is_valid_range(modrange,True,False) or is_valid_range(modrange,False,False)

    # set limits 
    lim = [-float('inf'),float('inf')]
    if type(set_limit) != type(None):
        if set_limit in {int,np.int32,np.int64}: 
            lim[0] = np.iinfo(np.int32).min 
            lim[1] = np.iinfo(np.int32).max 
        elif set_limit in {float,np.float32,np.float64}: 
            lim[0] = np.iinfo(np.float64).min 
            lim[1] = np.iinfo(np.float64).max 

    x1 = modulo_in_range(av,modrange)
    
    # case: do nothing
    if to_trinary_relation(x1,pv) == sign: 
        return modrange

    newrange = [modrange[0],modrange[1]]

    stat = True
    while stat:
        newrange[modindex] += moddir 
        x1 = modulo_in_range(av,newrange)
        stat = to_trinary_relation(x1,pv) != sign
        
        stat2 = lim[0] <= newrange[modindex] <= lim[1] 
        stat = stat and stat2 

    return newrange 

"""
NOTE: Method is specifically for `av` that does not satisfy the 
`is_value_below_modulo_range_length(av,modrange)` relation. 

Adjusts a modulo range `modrange` to M_r so that the output from 
`modulo_in_range(av,M_r)` satisfies the `sign`, -1 for `pv > av` or 
1 for `pv < av`. 

pv := "parent" value, also a "reference" for the binary relation with 
      `av`.
av := actual value, v = m * pv + a. 
sign := wanted binary relation between `av` and `pv`. 
modrange := base modulo range to adjust. 
"""
def multimodular_number__modulo_range_adjustment(pv,av,sign,modrange):
    assert sign in {-1,1}

    if sign == 1: 
        md,mi = 1,1
    else: 
        md,mi = -1,0
    return iterative_adjustment_for_direction__modrange(pv,av,sign,\
        modrange,mi,md,int) 

"""
NOTE: Method is specifically for `av` that satisfies the 
`is_value_below_modulo_range_length(av,modrange)` relation. 

wv := wanted value
av := actual value
cv := current value, the output from `modulo_in_range(av,modrange)`
modrange := base modulo range to adjust. 
"""
def unimodular_number__modulo_range_adjustment(wv,av,cv,modrange):

    sign = to_trinary_relation(wv,cv) 
    assert sign != 0 

    # adjust modrange if necessary
    cv_ = cv 
    mr = [modrange[0],modrange[1]] 
    if sign == 1: 
        if mr[1] < wv:
            mr[1] = wv
            cv_ = modulo_in_range(av,mr)
            mr[1] += wv - cv_ 
        else: 
            if mr[1] == wv: 
                mr[1] = wv + 1 
            d = wv - av 
            mr[0] = d 
    else:
        if mr[0] > wv: 
            mr[0] = wv
            cv_ = modulo_in_range(av,mr)
            mr[1] = mr[1] - cv + wv 
        elif mr[0] == wv: 
            mr[1] = mr[1] - (cv_ - mr[0])
        else: 
            mr[0] -= (cv_ - wv) 
    return mr 

"""
Calculates a new modrange M_r such that for the actual 
value `av`, `modulo_in_range(av,M_r) == pv`. 
"""
# NOTE: some cases yield wrong modranges 
def modrange_for_congruence(pv,av,modrange):
    assert is_number(pv,{complex,np.complex128})
    assert is_number(av,{complex,np.complex128})
    assert is_valid_range(modrange,True,False) or is_valid_range(modrange,False,False)

    if not is_value_below_modulo_range_length(av,modrange): 
        if pv >= 0: 
            mr2 = [0,abs(av)] 

            m = 1 if av < 0 else -1 
            mr2[1] += m * pv 
        else:            
            av_ = -av if av > 0 else av 

            pvs = pv >= 0 
            avs = av >= 0 
            s = -1 if pvs == avs else 1 
            mr2 = [av_ + 1,0] 
            mr2[0] = mr2[0] + s * pv - 1
        return mr2 
    elif not is_value_below_modulo_range_length(pv,modrange): 
        start = (pv - av)
        end = start + modrange[1] - modrange[0]
        return (start,end) 

    cv = modulo_in_range(av,modrange)
    difference = cv - pv 
    mr2 = [modrange[0],modrange[1]]
    
    if av >= 0: 
        mr2[0] -= difference 
    else: 
        mr2[1] -= difference 
    return mr2 

#--------------------------------------------------------------------- 