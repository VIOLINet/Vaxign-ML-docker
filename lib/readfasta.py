#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 15:18:43 2018

@author: edison
"""

def readFasta( fileName, key="full", strip = True ):
    fasta = {}
    curID = ''
    curSequence = ''
    for line in open( fileName ).readlines():
        if line[0] == '>':
            if strip:
                line = line.strip().strip( '>' )
            if curID != '':
                fasta[curID] = curSequence
                curID = ''
                curSequence = ''
            if key == "short":
                curID = line.split( ' ' )[0]
            elif key == "vi":
                tokens = line.split( '|' )
                for ( i, token ) in enumerate( tokens ):
                    if token == "vi" or token == ">vi":
                        curID = tokens[i+1]
                        break
            else:
                curID = line
        else:
            if strip:
                line = line.strip()
            curSequence = curSequence + line
    fasta[curID] = curSequence
    return fasta

def readFastaDesc( fileName, key="full", strip=True ):
    desc = []
    for line in open( fileName ).readlines():
        if line[0] == '>':
            if strip:
                line = line.strip().strip( '>' )
            if key == "short":
                desc.append( line.split( ' ' )[0] )
            elif key == "vi":
                tokens = line.split( '|' )
                for ( i, token ) in enumerate( tokens ):
                    if token == "vi" or token == ">vi":
                        desc.append( tokens[i+1] )
                        break
            else:
                desc.append( line )
    return desc