#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 17:45:42 2019

@author: edison
"""

import os
import sys
import shutil
import argparse
import multiprocessing

import numpy as np

from Bio import SeqIO
from pathlib import Path

from sklearn.externals import joblib

from mRMR import mRMR
from xgboost import XGBClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import MinMaxScaler
from sklearn.feature_selection import SelectKBest
from sklearn.model_selection import GridSearchCV, StratifiedShuffleSplit

LIB_PATH = os.path.dirname( os.path.realpath( __file__ ) )

sys.path.append( LIB_PATH )
import readfasta
from feature import Feature

class Train:
    
    def __init__( self ):
        self.version = "1.0"
    
    def train( self, positiveFasta, negativeFasta, outputDir, organism, incFeatures, multiFlag, process ):
        featureDir = os.path.join( outputDir, "_FEATURE" )
        masterLabels = ["ID", "Label", "Gram"]
        masterData = {}
        
        for (g, inputFasta) in enumerate( [negativeFasta, positiveFasta] ):
            sequenceIDs = readfasta.readFastaDesc( inputFasta, key="full" )
            if organism.lower() in ["gram+","g+"]:
                for fastaID in sequenceIDs:
                    masterData[fastaID] = [fastaID, str(g), "1"]
                for method in incFeatures:
                    tsvFile = os.path.join( featureDir, method.upper(), "%s.%s.tsv" % ( Path( inputFasta ).stem, method ) )
                    for (i, line) in enumerate( open( tsvFile ).read().splitlines() ):
                        tokens = line.split( '\t' )
                        if i == 0:
                            if method == "psortb":
                                masterLabels += tokens[2:]
                            else:
                                masterLabels += tokens[1:]
                        else:
                            if method == "psortb":
                                masterData[tokens[0]] += tokens[2:]
                            else:
                                masterData[tokens[0]] += tokens[1:]
            elif organism.lower() in ["gram-","g-"]:
                for fastaID in sequenceIDs:
                    masterData[fastaID] = [fastaID, str(g), "0"]
                for method in incFeatures:
                    tsvFile = os.path.join( featureDir, method.upper(), "%s.%s.tsv" % ( Path( inputFasta ).stem, method ) )
                    for (i, line) in enumerate( open( tsvFile ).read().splitlines() ):
                        tokens = line.split( '\t' )
                        if i == 0:
                            if method == "psortb":
                                masterLabels += tokens[2:]
                            else:
                                masterLabels += tokens[1:]
                        else:
                            if method == "psortb":
                                masterData[tokens[0]] += tokens[2:]
                            else:
                                masterData[tokens[0]] += tokens[1:]

        rows = []
        samples = []
        groups = []
        for fastaID in masterData.keys():
            tokens = masterData[fastaID]
            samples.append( fastaID )
            groups.append( int( tokens[1] ) )
            rows.append( tokens[2:] )
        X = np.array( rows )
        X = X.astype( float )
        y = np.array( groups )
        
        scaler = MinMaxScaler( copy=False )
        scaler.fit( X )
        joblib.dump( scaler, os.path.join( outputDir, "Scaler.sav" ) )
        
        X = scaler.transform( X )
        est = XGBClassifier( objective='binary:logistic', silent=True, nthread=1, 
                    eval_metric='auc', random_state=26 )

        estPipe = Pipeline( [('feature_selection', SelectKBest( mRMR )),('classification', est)] )
        grid = [{
            "feature_selection__k": list(range(20, X.shape[1], 20))[:10],
            'classification__learning_rate': [0.3, 0.1],
            
            'classification__n_estimators': [60, 80, 100, 120, 140, 160],
            
            'classification__max_depth': [3, 6, 9],
            'classification__min_child_weight': [1, 3],
            
            'classification__scale_pos_weight': [1, 6],
            
            'classification__max_delta_step': [0, 3],
        }]
    
        cv = StratifiedShuffleSplit( n_splits=3, random_state=6 )
        if multiFlag.lower() not in ['t','true']:
            xgb = GridSearchCV( estimator=estPipe, param_grid=grid, cv=cv, iid=False,
                               verbose=1, n_jobs=1 )
        else:
            xgb = GridSearchCV( estimator=estPipe, param_grid=grid, cv=cv, iid=False,
                               verbose=1, n_jobs=process )
        xgb.fit( X, y )
        joblib.dump( xgb, os.path.join( outputDir, "VaxignML.sav" ) )
        
        y_prob = xgb.predict_proba(X)
        joblib.dump( y_prob, os.path.join( outputDir, "VaxignML.scores" ) )
    
    def check_args( self, args ):
        if not os.path.exists( args.positiveFasta ):
            sys.stderr.write( "Positive sample fasta not found! Please check your input file.\n" )
            exit(1)
        elif not any( SeqIO.parse( open( args.positiveFasta, 'r' ), "fasta" ) ):
            sys.stderr.write( "Invalid positive sample fasta! Please check your input file is in fasta format.\n" )
            exit(1)
        if not os.path.exists( args.negativeFasta ):
            sys.stderr.write( "Negative sample fasta not found! Please check your input file.\n" )
            exit(1)
        elif not any( SeqIO.parse( open( args.negativeFasta, 'r' ), "fasta" ) ):
            sys.stderr.write( "Invalid negative sample fasta! Please check your input file is in fasta format.\n" )
            exit(1)
        if not os.path.exists( args.outputDir ):
            os.mkdir( args.outputDir )
        if args.organism.lower() not in ["gram+","g+","gram-","g-"]:
            sys.stderr.write( "Incorrect organism! Please choose from: [gram+,gram-]\n" )
            exit(1)
        if args.multiFlag.lower() not in ['t','true','f','false']:
            sys.stderr.write( "Incorrect input for multi-processes! Please choose from: [T,F]\n" )
            exit(1)
        if args.process >= multiprocessing.cpu_count():
            sys.stderr.write( "Requested processes exceed limit!\n" )
            exit(1)
        if args.rawFlag.lower() not in ['t','true','f','false']:
            sys.stderr.write( "Incorrect input for keeping raw output file! Please choose from: [T,F]\n" )
            exit(1)
    
    def main( self ):
        try:
            parser = argparse.ArgumentParser( description="Train Vaxign-ML model." )
            parser.add_argument( "--positive", '-p', dest='positiveFasta', required=True )
            parser.add_argument( "--negative", '-n', dest='negativeFasta', required=True )
            parser.add_argument( "--output", '-o', dest='outputDir', required=True )
            parser.add_argument( "--organismtype", '-t', dest='organism', required=True)
            parser.add_argument( "--multi", '-m', dest='multiFlag', default="True" )
            parser.add_argument( "--cpuprocess",'-c', dest='process', type=int, default=int(multiprocessing.cpu_count()/2) )
            parser.add_argument( "--raw", '-r', dest='rawFlag', default="False" )
            args = parser.parse_args()
            self.check_args( args )
            
            featureDir = os.path.join( args.outputDir, "_FEATURE" )
            if not os.path.exists( featureDir ):
                os.mkdir( featureDir )
            
            incFeatures = []
            incFeatures.append( "psortb" )
            Feature.run_psortb( args.positiveFasta, featureDir, args.organism, args.multiFlag, args.process, args.rawFlag )
            Feature.run_psortb( args.negativeFasta, featureDir, args.organism, args.multiFlag, args.process, args.rawFlag )
            incFeatures.append( "spaan" )
            Feature.run_spaan( args.positiveFasta, featureDir, args.rawFlag )
            Feature.run_spaan( args.negativeFasta, featureDir, args.rawFlag )
            incFeatures.append( "signalp" )
            Feature.run_signalp( args.positiveFasta, featureDir, args.organism, args.multiFlag, args.process, args.rawFlag )
            Feature.run_signalp( args.negativeFasta, featureDir, args.organism, args.multiFlag, args.process, args.rawFlag )
            incFeatures.append( "tmhmm" )
            Feature.run_tmhmm( args.positiveFasta, featureDir, args.multiFlag, args.process, args.rawFlag )
            Feature.run_tmhmm( args.negativeFasta, featureDir, args.multiFlag, args.process, args.rawFlag )
            incFeatures.append( "imgen" )
            Feature.run_immugen( args.positiveFasta, featureDir, args.rawFlag )
            Feature.run_immugen( args.negativeFasta, featureDir, args.rawFlag )
            incFeatures.append( "mdesc" )
            Feature.run_descriptor( args.positiveFasta, featureDir, args.rawFlag )
            Feature.run_descriptor( args.negativeFasta, featureDir, args.rawFlag )
            
            self.train( args.positiveFasta, args.negativeFasta, args.outputDir, args.organism, incFeatures, args.multiFlag, args.process )
            
            if args.rawFlag.lower() in ['f','false']:
                shutil.rmtree( featureDir )
        except:
            print( sys.exc_info() )
        

if __name__ == "__main__":
    Train().main()