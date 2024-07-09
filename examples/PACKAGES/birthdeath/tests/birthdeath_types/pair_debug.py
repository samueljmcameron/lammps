import numpy as np
import matplotlib.pyplot as plt

from generalised_forces import division_offset,death_offset


if __name__ == "__main__":

    from lammpstools import DumpLoader
    import sys


    dumpfile = sys.argv[1]
    b0 = 1.0
    d0 = 3.0
    cutoff = 1.5
    sigma = 0.9
    cnum = 6.
    width = 0.001
    
    dl = DumpLoader(dumpfile,integerquantities=['id','type'])


    N = dl.data[-1]['N']
    xs = dl.data[-1]['x']
    ys = dl.data[-1]['y']
    zs = dl.data[-1]['z']
    divisions = dl.data[-1]['division']
    deaths = dl.data[-1]['death']
    types = dl.data[-1]['type']


    expected_divisions = np.zeros([N],float)
    expected_deaths = np.zeros([N],float)

    for i in range(N):
        
        if types[i] == 1:
            expected_divisions[i] += b0
        for j in range(i+1,N):

            rij = np.array([xs[i] - xs[j],
                            ys[i] - ys[j],
                            zs[i] - zs[j]])

            rdist = np.linalg.norm(rij)
            print(rdist)

            if types[i] == 1 and types[j] == 1:
                righttype = True
            else:
                righttype = False

            if rdist <= cutoff and righttype:
            
                expected_divisions[i] -= division_offset(rdist,b0,sigma,width,cnum)
                expected_deaths[i] += death_offset(rdist,d0,sigma,width,cnum)
                expected_divisions[j] -= division_offset(rdist,b0,sigma,width,cnum)
                expected_deaths[j] += death_offset(rdist,d0,sigma,width,cnum)



    for i in range(N):
        print(xs[i],ys[i],zs[i],expected_divisions[i],divisions[i])
        print(xs[i],ys[i],zs[i],expected_deaths[i],deaths[i])
