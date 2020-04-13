#!/usr/bin/env python3

import numpy as np

# ---------------------------- Processing functions ---------------------------
def agc(ts, k, agctype):
    # Auto-gain normalization using a running window
    # k is the window length
    # agctype is either mean, rms, median 
    n = len(ts)

    k2 = int(k/2) # This will floor the number if k is even; (k-1)/2 
    if np.mod(k, 2) == 0: # even numbers need will not have a centered window
        k = int( k + 1) 
    
    stat = np.ones([n])
    # Normalize
    if agctype == "std":
        for i in range(0, k2):
            stat[i] = np.std( abs( ts[0:(i+k2)] ) )
            stat[ (n-k2+i) ] = np.std( abs(ts[ (n-2*k2+i):n] ) )
        for i in range(k2,n-k2):
            stat[i] = np.std( abs( ts[ (i-k2):(i+k2) ] ) )
    elif agctype == "mean":
        for i in range(0, k2):
            stat[i] = np.mean( abs( ts[0:(i+k2)] ) )
            stat[ (n-k2+i) ] = np.mean( abs(ts[ (n-2*k2+i):n] ) )
        for i in range(k2,n-k2):
            stat[i] = np.mean( abs( ts[ (i-k2):(i+k2) ] ) )
    else:
        for i in range(0, k2):
            stat[i] = np.std( ts[i:(i+k2)] )
            stat[ (n-k2+i) ] = np.std( ts[ (n-2*k2+i):n] )
        for i in range(k2,n-k2):
            stat[i] = np.std( ts[ (i-k2):(i+k2) ] )
    
    stat[stat == 0] = 1
    ts = ts/stat 
    return ts

