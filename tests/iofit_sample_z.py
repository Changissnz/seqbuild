from morebs2.search_space_iterator import * 
from mini_dm.vs_fit_op import * 

def IOFit_sample_z(): 

    m = np.array([4,5,16,27,54,65]) 
    a = 43.0 
    ad = AffineDelta(m,a,0) 

    bounds = np.array([[-1.0,1.0],\
                    [2.0,7.0],\
                        [14.0,22.0],\
                        [-20,-2],\
                        [55.,105.0],\
                        [200.,490.]])

    start_point = deepcopy(bounds[:,0])
    column_order = np.arange(6)
    ssi_hop = np.array([2,5,8,18,50,290]) 
    cycle_on = True 
    cycle_is = 0 

    ssi = SearchSpaceIterator(bounds,start_point,\
        column_order,ssi_hop,cycle_on,cycle_is)
    xs,ys = [],[]

    for i in range(100):
        x = next(ssi) 
        xs.append(x) 
        ys.append(ad.fit(x))

    iof = IOFit(xs,ys,None,None,None)
    return iof 