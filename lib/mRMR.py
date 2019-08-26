#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 16:18:11 2019

@author: edison
"""

def mRMR(X,y):
    import os
    import pandas
    import tempfile
    import subprocess
    
    import numpy as np
    
    dir_path = os.path.dirname( os.path.realpath( __file__ ) )
    mrmr_bin = os.path.join( dir_path, "mrmr" )
    
    colNames = ["class"]+["V"+str(a) for a in range(1,X.shape[1]+1)]
    df = pandas.DataFrame(np.hstack((np.reshape(y,(y.size,1)),X)), columns=colNames)
    
    with tempfile.NamedTemporaryFile( 'w', dir = dir_path ) as tmpIn:
        df.to_csv( tmpIn.name, index=False )
        result = subprocess.check_output( [mrmr_bin, "-t", "1", "-i", tmpIn.name, "-n", "220", "-s", str(X.shape[0])] )
        tmp = {}
        for row in result.decode().split( '***' )[4].strip().split( '\n' )[1:]:
            tokens = row.split( '\t' )
            tmp[int(tokens[1].strip())] = float(tokens[-1].strip())
        scores = []
        for i in range(1,X.shape[1]+1):
            if i in tmp:
                scores.append( tmp[i] )
            else:
                scores.append( np.min(list(tmp.values())) - 0.001 )
        scores = np.array( scores, dtype=float )
        tmpIn.close()
    
    if len( scores ) < X.shape[1]:
        scores = np.concatenate( (scores, [scores[-1] - 0.001] * (X.shape[1] - len(scores))) )
    
    return scores