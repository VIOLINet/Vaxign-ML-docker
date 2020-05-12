#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 16:18:11 2019

@author: edison
"""

def mRMR(X,y):
    import os
    import pandas
    import hashlib
    import tempfile
    import subprocess
    
    import numpy as np
    
    from pathlib import Path
    
    dir_path = os.path.dirname( os.path.realpath( __file__ ) )
    mrmr_bin = os.path.join( dir_path, "mrmr" )
    computed_dir_path = os.path.join( dir_path, "mrmr.computed" )
    if not os.path.exists( computed_dir_path ):
        os.mkdir( computed_dir_path )
    computed_files = os.listdir( computed_dir_path )
    md5s = {}
    for file in computed_files:
        md5s[Path(file).stem] = file
    
    colNames = ["class"]+["V"+str(a) for a in range(1,X.shape[1]+1)]
    df = pandas.DataFrame(np.hstack((np.reshape(y,(y.size,1)),X)), columns=colNames)
    
    with tempfile.NamedTemporaryFile( 'w', dir = dir_path ) as tmpIn:
        df.to_csv( tmpIn.name, index=False )
        md5 = hashlib.md5( open( tmpIn.name ).read().encode() ).hexdigest()
        if ( md5 in md5s.keys() ):
            result = open( os.path.join( computed_dir_path, md5s[md5] ) ).read()
            tmp = {}
            for row in result.split( '***' )[4].strip().split( '\n' )[1:]:
                tokens = row.split( '\t' )
                tmp[int(tokens[1].strip())] = float(tokens[-1].strip())
            scores = []
            for i in range(1,X.shape[1]+1):
                if i in tmp:
                    scores.append( tmp[i] )
                else:
                    scores.append( np.min(list(tmp.values())) - 0.001 )
            scores = np.array( scores, dtype=float )
        else:
            result = subprocess.check_output( [mrmr_bin, "-t", "1", "-i", tmpIn.name, "-n", "220", "-s", str(X.shape[0])] )
            open( os.path.join( computed_dir_path, "%s.out" % md5 ), 'w' ).write( result.decode() )
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