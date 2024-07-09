import numpy as np


def division_offset(rdist,b0,sigma,width,cnum):

    return b0/(1+np.exp((rdist-sigma)/width))/cnum


def death_offset(rdist,d0,sigma,width,cnum):

    return d0/(1+np.exp((rdist-sigma)/width))/cnum

