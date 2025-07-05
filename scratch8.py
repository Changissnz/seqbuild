from mini_dm.vs_fit_op import * 
from morebs2.numerical_generator import prg__LCG,prg_partition_for_sz,\
    prg__n_ary_alternator,prg__constant,prg_seqsort

m_type = "vec"
a_type = "vec" 
prg = prg__LCG(71,688,31,900) 
r_out1 = prg__constant((50.,5000.))
r_out2 = prg__constant((50.,5000.))
ma_order = 0

ad = AffineDelta.one_instance(m_type,a_type,prg,r_out1,r_out2,dim_range=None,ma_order=ma_order)

x1,x2 = 56,400 
y1 = ad.fit(x1)
y2 = ad.fit(x2)

ad2 = IOFit.io_pointpair_to_AffineDelta(x1,x2,y1,y2,ma_order=0)

ad.ma_order = 1
y1 = ad.fit(x1) 
y2 = ad.fit(x2) 
ad3 = IOFit.io_pointpair_to_AffineDelta(x1,x2,y1,y2,ma_order=1)

ad2.ma_order = 1 
assert ad == ad2 

assert ad2 == ad3 