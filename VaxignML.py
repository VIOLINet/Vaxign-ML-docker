#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 11:27:40 2019

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
from scipy.stats import percentileofscore

SCRIPT_PATH = os.path.dirname( os.path.realpath( __file__ ) )
LIB_PATH = os.path.join( SCRIPT_PATH, "lib" )
MODEL_PATH = os.path.join( SCRIPT_PATH, "model" )

sys.path.append( SCRIPT_PATH )
sys.path.append( LIB_PATH )
import readfasta
from feature import Feature

class VaxignML:

    def makeInput( self, inputFasta, outputDir, organism, incFeatures ):
        featureDir = os.path.join( outputDir, "_FEATURE" )
        if organism.lower() in ["gram+","g+","gram-","g-"]:
            masterLabels = ["ID", "Gram"]
            masterData = {}
    
            sequenceIDs = readfasta.readFastaDesc( inputFasta, key="full" )
            if organism.lower() in ["gram+","g+"]:
                for fastaID in sequenceIDs:
                    masterData[fastaID] = [fastaID, "1"]
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
                    masterData[fastaID] = [fastaID,"0"]
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
        else:
            masterLabels = ["ID"]
            masterData = {}
    
            sequenceIDs = readfasta.readFastaDesc( inputFasta, key="full" )
            for fastaID in sequenceIDs:
                masterData[fastaID] = [fastaID]
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

        output = ["\t".join( masterLabels )]
        for fastaID in masterData.keys():
            output.append( "\t".join( masterData[fastaID] ) )
        open( os.path.join( outputDir, "%s.input.tsv" % Path( inputFasta ).stem  ), 'w' ).write( '\n'.join( output ) )
        
    def predict( self, inputFasta, outputDir, modelDir, model='virus'):
        inputFile = os.path.join( outputDir, "%s.input.tsv" % Path( inputFasta ).stem  )
        
        if (sys.version_info > (3, 0)):
            scaler = joblib.load( os.path.join( modelDir, "Scaler_%s.sav" % model ) )
            vaxignML = joblib.load( os.path.join( modelDir, "VaxignML_%s.sav" % model ) )
            scores = joblib.load( os.path.join( modelDir, "VaxignML_%s.scores" % model ) )
        else:
            scaler = joblib.load( os.path.join( modelDir, "Scaler_%s.sav.2" % model ) )
            vaxignML = joblib.load( os.path.join( modelDir, "VaxignML_%s.sav.2" % model ) )
            scores = joblib.load( os.path.join( modelDir, "VaxignML_%s.scores.2" % model ) )
        
        if model == 'bacteria':
            labels = open( inputFile ).read().splitlines()[0].split( '\t' )[1:]
            samples = []
            X = np.array( [] ).reshape( 0, len( labels ) )
            groups = []
            for line in open( inputFile ).read().splitlines()[1:]:
                tokens = line.split( '\t' )
                fastaID = tokens[0]
                samples.append( fastaID )
                groups.append( int( tokens[1] ) )
                value = np.array( tokens[1:] ).reshape( 1, len( labels ) )
                X = np.concatenate( ( X, value ), axis = 0 )
            X = X.astype( float )
            X = scaler.transform( X )
        else:
            labels = open( inputFile ).read().splitlines()[0].split( '\t' )[1:]
            samples = []
            X = np.array( [] ).reshape( 0, len( labels ) )
            for line in open( inputFile ).read().splitlines()[1:]:
                tokens = line.split( '\t' )
                fastaID = tokens[0]
                samples.append( fastaID )
                value = np.array( tokens[1:] ).reshape( 1, len( labels ) )
                X = np.concatenate( ( X, value ), axis = 0 )
            X = X.astype( float )
            X = scaler.transform( X )
        
        y_pred = vaxignML.predict( X )
        y_prob = vaxignML.predict_proba( X )
        output = ["sample\tprediction\tprotegenicity"]
        for (i, fastaID) in enumerate( samples ):
            output.append( "%s\t%f\t%f" % ( fastaID, y_pred[i], percentileofscore(scores[:,1], y_prob[i,1]) ) )
        open( os.path.join( outputDir, "%s.result.tsv" % Path( inputFasta ).stem  ), 'w' ).write( '\n'.join( output ) )

    def check_args( self, args ):
        if not os.path.exists( args.inputFasta ):
            sys.stderr.write( "Input fasta not found! Please check your input file.\n" )
            exit(1)
        elif not any( SeqIO.parse( open( args.inputFasta, 'r' ), "fasta" ) ):
            sys.stderr.write( "Invalid input fasta! Please check your input file is in fasta format.\n" )
            exit(1)
        if not os.path.exists( args.outputDir ):
            os.mkdir( args.outputDir )
        if not ( os.path.exists( args.savedModel ) 
            and os.path.exists( os.path.join( MODEL_PATH, "Scaler_bacteria.sav" ) )
            and os.path.exists( os.path.join( MODEL_PATH, "VaxignML_bacteria.sav" ) )
            and os.path.exists( os.path.join( MODEL_PATH, "VaxignML_bacteria.scores" ) )
            and os.path.exists( os.path.join( MODEL_PATH, "Scaler_virus.sav" ) )
            and os.path.exists( os.path.join( MODEL_PATH, "VaxignML_virus.sav" ) )
            and os.path.exists( os.path.join( MODEL_PATH, "VaxignML_virus.scores" ) )
            ):
            sys.stderr.write( "Unable to find previous train model! Please check your input.\n" )
            exit(1)
        if args.organism.lower() not in ["gram+","g+","gram-","g-", "virus", "v"]:
            sys.stderr.write( "Incorrect organism! Please choose from: [gram+,gram-,virus]\n" )
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
            parser = argparse.ArgumentParser( description="Vaxign-ML predicts bacterial protective antigens." )
            parser.add_argument( "--input", '-i', dest='inputFasta', required=True )
            parser.add_argument( "--output", '-o', dest='outputDir', required=True )
            parser.add_argument( "--organismtype", '-t', dest='organism', required=True)
            parser.add_argument( "--savedModel", '-s', dest='savedModel', default=MODEL_PATH )
            parser.add_argument( "--multi", '-m', dest='multiFlag', default="True" )
            parser.add_argument( "--process", '-p', dest='process', type=int, default=int(multiprocessing.cpu_count()/2) )
            parser.add_argument( "--raw", '-r', dest='rawFlag', default="False" )
            args = parser.parse_args()
            self.check_args( args )
            
            featureDir = os.path.join( args.outputDir, "_FEATURE" )
            if not os.path.exists( featureDir ):
                os.mkdir( featureDir )
            if args.organism.lower() in ["gram+","g+","gram-","g-"]:
                incFeatures = []
                incFeatures.append( "psortb" )
                Feature.run_psortb( args.inputFasta, featureDir, args.organism, args.multiFlag, args.process, args.rawFlag )
                incFeatures.append( "spaan" )
                Feature.run_spaan( args.inputFasta, featureDir, args.rawFlag )
                incFeatures.append( "signalp" )
                Feature.run_signalp( args.inputFasta, featureDir, args.organism, args.multiFlag, args.process, args.rawFlag )
                incFeatures.append( "tmhmm" )
                Feature.run_tmhmm( args.inputFasta, featureDir, args.multiFlag, args.process, args.rawFlag )
                incFeatures.append( "imgen" )
                Feature.run_immugen( args.inputFasta, featureDir, args.rawFlag )
                incFeatures.append( "mdesc" )
                Feature.run_descriptor( args.inputFasta, featureDir, args.rawFlag )
                
                self.makeInput( args.inputFasta, args.outputDir, args.organism, incFeatures )
                
                if args.rawFlag.lower() in ['f','false']:
                    shutil.rmtree( featureDir )
                
                self.predict( args.inputFasta, args.outputDir, args.savedModel, model="bacteria" )
            elif args.organism.lower() in ["virus","v"]:
                incFeatures = []
                incFeatures.append( "spaan" )
                Feature.run_spaan( args.inputFasta, featureDir, args.rawFlag )
                incFeatures.append( "tmhmm" )
                Feature.run_tmhmm( args.inputFasta, featureDir, args.multiFlag, args.process, args.rawFlag )
                incFeatures.append( "imgen" )
                Feature.run_immugen( args.inputFasta, featureDir, args.rawFlag )
                incFeatures.append( "mdesc" )
                Feature.run_descriptor( args.inputFasta, featureDir, args.rawFlag )
                
                self.makeInput( args.inputFasta, args.outputDir, args.organism, incFeatures )
                
                if args.rawFlag.lower() in ['f','false']:
                    shutil.rmtree( featureDir )
                
                self.predict( args.inputFasta, args.outputDir, args.savedModel, model="virus" )
            
        except:
            print( sys.exc_info() )
        

if __name__ == "__main__":
    VaxignML().main()