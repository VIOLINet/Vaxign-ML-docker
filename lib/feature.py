#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 18:56:49 2018

@author: edison
"""

import os
import re
import sys
import glob
import math
import shutil
import argparse
import subprocess
import multiprocessing

from Bio import SeqIO
from pathlib import Path

SCRIPT_PATH = os.path.dirname( os.path.realpath( __file__ ) )
PSORTB_PATH = os.path.join( "/usr", "local", "psortb", "bin", "psort" )
SIGNALP_PATH = os.path.join( SCRIPT_PATH, "signalp", "4.1", "signalp" )
TMHMM_PATH = os.path.join( SCRIPT_PATH, "tmhmm", "2.0c", "bin", "tmhmm" )
SPAAN_PATH = os.path.join( SCRIPT_PATH, "spaan", "SPAAN", "askquery.pl" )
IMGEN_PATH = os.path.join( SCRIPT_PATH, "iedb", "immunogenicity", "predict_immunogenicity.py" )

SPLIT_LIMIT = 1000.0
MIN_PEPTIDE_LENGTH = 16

sys.path.append( SCRIPT_PATH )
import readfasta
from propy.PyPro import GetProDes

class Feature:

    def __init__( self ):
        self.version = "1.0"
    
    @staticmethod
    def split_files( inputFasta, tmpDir ):
        fasta = readfasta.readFasta( inputFasta, key="full", strip=False )
        sequenceIDs = readfasta.readFastaDesc( inputFasta, key="full" )
        inFiles = []
        outFiles = []
        split = int( math.ceil( len( sequenceIDs ) ) / SPLIT_LIMIT )
        size = int( math.ceil( len( sequenceIDs ) / float( split+1 ) ) )
        for i in range( split+1 ):
            textBuffer = ''
            for fastaID in list(fasta.keys())[i*size:(i+1)*size]:
                textBuffer = textBuffer + fastaID + fasta[fastaID]
            inFile = os.path.join( tmpDir, "input.fasta.%i" % i )
            open( inFile,'w' ).write( textBuffer )
            inFiles.append( inFile )
            outFile = os.path.join( tmpDir, "output.raw.%i" % i )
            outFiles.append( outFile )
        return ( inFiles, outFiles )
    
    @staticmethod
    def combine_and_cleanup_files( inFiles, outFiles, finalFile ):
        open( finalFile, 'w' ).write( '' )
        for outFile in outFiles:
            open( finalFile, 'a' ).write( open( outFile, 'r' ).read() )
            os.remove( outFile )
        for inFile in inFiles:
            os.remove( inFile )
        
    @staticmethod
    def run_psortb( inputFasta, outputDir, organism, multiFlag, process, rawFlag ):
        command = PSORTB_PATH
        rawOutput = os.path.join( outputDir, "PSORTB" )
        if not os.path.exists( rawOutput ):
            os.mkdir( rawOutput )
        
        sequenceIDs = readfasta.readFastaDesc( inputFasta, key="full" )
        
        for file in os.listdir( rawOutput ):
                if 'psortb_grampos.txt' in file:
                    os.remove( os.path.join( rawOutput, file ) )
        
        if organism.lower() in ["gram+","g+"]:
            os.system( "%s -p -i %s" % (command, inputFasta) )
            for file in os.listdir( rawOutput ):
                if 'psortb_grampos.txt' in file:
                    rawFile = file
                    break
        elif organism.lower() in ["gram-","g-"]:
            os.system( "%s -n -i %s" % (command, inputFasta) )
            for file in os.listdir( rawOutput ):
                if 'psortb_gramneg.txt' in file:
                    rawFile = file
                    break
        
        output = ['\t'.join( ["ID","SubcellularLocation","Extracellular_Probability", "CytoplasmicMembrane_Probability", "Cytoplasmic_Probability",
                              "Cellwall_Probability","Periplasmic_Probability","OuterMembrane_Probability"] )]
        locs = ["Extracellular", "CytoplasmicMembrane", "Cytoplasmic","Cellwall","Periplasmic","OuterMembrane"]
        values = {}
        for entry in re.split( '[-]{79}',open( os.path.join( rawOutput, rawFile ) ).read() )[:-1]:
            fastaID = entry.strip().split( '\n' )[0][7:].strip()
            value = [""]+["0.0"]*(len(locs))
            value[0] = re.split( '[ ]+', entry[entry.find( "Final Prediction:" ):].split( '\n' )[1] )[1]
            for row in entry[entry.find( "Localization Scores:" ):entry.find( "Final Prediction:" )].split( '\n' )[1:-1]:
                tokens = re.split( '[ ]+', row )
                value[locs.index(tokens[1])+1] = str( float( tokens[2] ) / 10.0 )
            values[fastaID] = value
        for fastaID in sequenceIDs:
            output.append( '\t'.join( [fastaID] + values[fastaID] ) )
        open( os.path.join( rawOutput, "%s.psortb.tsv" % Path( inputFasta ).stem ), 'w' ).write( '\n'.join( output ) )
        if not rawFlag:
            os.remove( os.path.join( rawOutput, rawFile ) )
    
    @staticmethod
    def run_spaan( inputFasta, outputDir, rawFlag ):
        command = SPAAN_PATH
        rawOutput = os.path.join( outputDir, "SPAAN" )
        if not os.path.exists( rawOutput ):
            os.mkdir( rawOutput )
        
        sequenceIDs = readfasta.readFastaDesc( inputFasta, key="full" )
        
        os.system( "%s %s %s" % (command, inputFasta, os.path.join( rawOutput, "%s.output" % Path( inputFasta ).stem ) ) )
        
        output = ['ID\tSPAAN_Score']
        values = {}
        for row in open( os.path.join( rawOutput, "%s.output" % Path( inputFasta ).stem ) ).read().split( '\n' )[1:]:
            if row == '':
                continue
            tokens = row.split( '\t' )
            values[tokens[2].strip( '>' )] = tokens[1]
        for fastaID in sequenceIDs:
            if fastaID in values.keys():
                output.append( '\t'.join( [fastaID, values[fastaID]] ) )
            else:
                output.append( fastaID + '\t0.0' )
        open( os.path.join( rawOutput, "%s.spaan.tsv" % Path( inputFasta ).stem ), 'w' ).write( '\n'.join( output ) )
        if not rawFlag:
            os.remove( os.path.join( rawOutput, "%s.output" % Path( inputFasta ).stem ) )
    
    @staticmethod
    def run_signalp( inputFasta, outputDir, organism, multiFlag, process, rawFlag ):
        command = SIGNALP_PATH
        rawOutput = os.path.join( outputDir, "SIGNALP" )
        if not os.path.exists( rawOutput ):
            os.mkdir( rawOutput )
        
        sequenceIDs = readfasta.readFastaDesc( inputFasta, key="full" )
        
        if multiFlag.lower() in ['t','true'] and len(sequenceIDs) > SPLIT_LIMIT:
            ( inFiles, outFiles ) = Feature.split_files( inputFasta, rawOutput )
        else:
            inFiles = [inputFasta]
            outFiles = [os.path.join( rawOutput, "%s.output" % Path( inputFasta ).stem )]
        
        clines = []
        for ( f, inFile ) in enumerate( inFiles ):
            outFile = outFiles[f]
            if organism.lower() in ["gram+","g+"]:
                cline = "%s -t gram+ %s > %s" % (command, inFile, outFile)
            elif organism.lower() in ["gram-","g-"]:
                cline = "%s -t gram- %s > %s" % (command, inFile, outFile)
            clines.append( cline )
        Feature.mp_run( clines, process )
        
        if multiFlag.lower() in ['t','true'] and len(sequenceIDs) > SPLIT_LIMIT:
            Feature.combine_and_cleanup_files( inFiles, outFiles, os.path.join( rawOutput, "%s.output" % Path( inputFasta ).stem ) )
        
        fasta = readfasta.readFasta( inputFasta, key="full", strip=False )
        sequenceIDs = {}
        for fastaID in fasta.keys():
            shortID = fastaID.strip().strip( '>' ).split( ' ' )[0]
            sequenceIDs[shortID] = fastaID.strip().strip( '>' )
        values = {}
        for ( i, row ) in enumerate( open( os.path.join( rawOutput, "%s.output" % Path( inputFasta ).stem ) ).readlines() ):
            if row == '' or row.startswith( '#' ):
                continue
            tokens = re.split( '[ ]+', row )
            if tokens[0] in values.keys():
                print( "Duplicate SingalP Result: %s" % tokens[0] )
            values[tokens[0]] = tokens[8]
        
        output = ["ID\tSignalP_DScore"]
        for shortID in sequenceIDs.keys():
            if shortID in values.keys():
                output.append( '\t'.join( [sequenceIDs[shortID], values[shortID]] ) )
            else:
                output.append( '\t'.join( [sequenceIDs[shortID], '0.000'] ) )
        open( os.path.join( rawOutput, "%s.signalp.tsv" % Path( inputFasta ).stem ), 'w' ).write( '\n'.join( output ) )
        if not rawFlag:
            os.remove( os.path.join( rawOutput, "%s.output" % Path( inputFasta ).stem ) )
    
    @staticmethod
    def run_tmhmm( inputFasta, outputDir, multiFlag, process, rawFlag ):
        command = TMHMM_PATH
        rawOutput = os.path.join( outputDir, "TMHMM" )
        if not os.path.exists( rawOutput ):
            os.mkdir( rawOutput )
        
        sequenceIDs = readfasta.readFastaDesc( inputFasta, key="full" )
        
        if multiFlag.lower() in ['t','true'] and len(sequenceIDs) > SPLIT_LIMIT:
            ( inFiles, outFiles ) = Feature.split_files( inputFasta, rawOutput )
        else:
            inFiles = [inputFasta]
            outFiles = [os.path.join( rawOutput, "%s.output" % Path( inputFasta ).stem )]
        
        clines = []
        for ( f, inFile ) in enumerate( inFiles ):
            outFile = outFiles[f]
            cline = "%s %s > %s" % (command, inFile, outFile)
            clines.append( cline )
        Feature.mp_run( clines, process )
        
        if multiFlag.lower() in ['t','true'] and len(sequenceIDs) > SPLIT_LIMIT:
            Feature.combine_and_cleanup_files( inFiles, outFiles, os.path.join( rawOutput, "%s.output" % Path( inputFasta ).stem ) )
        
        fasta = readfasta.readFasta( inputFasta, key="full", strip=False )
        sequenceIDs = {}
        for fastaID in fasta.keys():
            shortID = fastaID.strip().strip( '>' ).split( ' ' )[0]
            sequenceIDs[shortID] = fastaID.strip().strip( '>' )
        values = {}
        for row in open( os.path.join( rawOutput, "%s.output" % Path( inputFasta ).stem ) ).read().split( '\n' ):
            if row.startswith( '#' ):
                tokens = row.split( ' ' )
                if tokens[1] not in values.keys():
                    values[tokens[1]] = [""] * 4
                if "Number of predicted TMHs:" in row:
                    values[tokens[1]][0] = tokens[-1].strip()
                if "Exp number of AAs in TMHs:" in row:
                    values[tokens[1]][1] = tokens[-1].strip()
                if "Exp number, first 60 AAs:" in row:
                    values[tokens[1]][2] = tokens[-1].strip()
                if "Total prob of N-in:" in row:
                    values[tokens[1]][3] = tokens[-1].strip()
        
        output = ["ID\tPredicted_TMH#\tExp_AAs#\tExp_first_60_AAs#\tTotal_N-in_prob"]
        for shortID in sequenceIDs.keys():
            if shortID in values.keys():
                output.append( '\t'.join( [sequenceIDs[shortID]] + values[shortID] ) )
            else:
                output.append( '\t'.join( [sequenceIDs[shortID]] + [""] * 4 ) )
        open( os.path.join( rawOutput, "%s.tmhmm.tsv" % Path( inputFasta ).stem ), 'w' ).write( '\n'.join( output ) )
        if not rawFlag:
            os.remove( os.path.join( rawOutput, "%s.output" % Path( inputFasta ).stem ) )
        for tmpFile in glob.glob( os.path.join( os.getcwd(), "TMHMM_*" ) ):
            shutil.rmtree( tmpFile )
    
    @staticmethod
    def run_descriptor( inputFasta, outputDir, rawFlag ):
        rawOutput = os.path.join( outputDir, "MDESC" )
        if not os.path.exists( rawOutput ):
            os.mkdir( rawOutput )
        
        fasta = readfasta.readFasta( inputFasta, key="full" )
        sequenceIDs = readfasta.readFastaDesc( inputFasta, key="full" )
        
        aacomp = []
        ctd = []
        seqord = []
        autocor = []
        output = []
        for fastaID in sequenceIDs:
            sequence = fasta[fastaID]
            Des = GetProDes(sequence)
            features = {}
            dict.update( features, Des.GetAAComp() )
            keys = list( features.keys() )
            if len(aacomp) == 0:
                aacomp.append( '\t'.join( ['ID'] + keys ) )
                output.append( '\t'.join( ['ID'] + keys ) )
            tmp = []
            for key in keys:
                tmp.append( str( features[key] ) )
            aacomp.append( '\t'.join( [fastaID] + tmp ) )
            output.append( '\t'.join( [fastaID] + tmp ) )
            
            features = Des.GetCTD()
            keys = list( features.keys() )
            if len( ctd ) == 0:
                ctd.append( '\t'.join( ['ID'] + keys ) )
                output[0] += '\t%s' % '\t'.join( keys )
            tmp = []
            for key in keys:
                tmp.append( str( features[key] ) )
            ctd.append( '\t'.join( [fastaID] + tmp ) )
            output[-1] += '\t%s' % '\t'.join( tmp )
            
            features = {}
            dict.update( features, Des.GetQSO( maxlag = MIN_PEPTIDE_LENGTH ) )
            keys = list( features.keys() )
            if len( seqord ) == 0:
                seqord.append( '\t'.join( ['ID'] + keys ) )
                output[0] += '\t%s' % '\t'.join( keys )
            tmp = []
            for key in keys:
                tmp.append( str( features[key] ) )
            seqord.append( '\t'.join( [fastaID] + tmp ) )
            output[-1] += '\t%s' % '\t'.join( tmp )
            
            features = {}
            dict.update( features, Des.GetMoreauBrotoAuto( maxlag = MIN_PEPTIDE_LENGTH ) )
            dict.update( features, Des.GetGearyAuto( maxlag = MIN_PEPTIDE_LENGTH ) )
            keys = list( features.keys() )
            if len( autocor ) == 0:
                autocor.append( '\t'.join( ['ID'] + keys ) )
                output[0] += '\t%s' % '\t'.join( keys )
            tmp = []
            for key in keys:
                tmp.append( str( features[key] ) )
            autocor.append( '\t'.join( [fastaID] + tmp ) )
            output[-1] += '\t%s' % '\t'.join( tmp )
        
        if rawFlag:
            open( os.path.join( rawOutput, "%s.aacomp.tsv" % Path( inputFasta ).stem ), 'w' ).write( '\n'.join( aacomp ) )
            open( os.path.join( rawOutput, "%s.ctd.tsv" % Path( inputFasta ).stem ), 'w' ).write( '\n'.join( ctd ) )
            open( os.path.join( rawOutput, "%s.seqord.tsv" % Path( inputFasta ).stem ), 'w' ).write( '\n'.join( seqord ) )
            open( os.path.join( rawOutput, "%s.autocor.tsv" % Path( inputFasta ).stem ), 'w' ).write( '\n'.join( autocor ) )
        open( os.path.join( rawOutput, "%s.mdesc.tsv" % Path( inputFasta).stem ), 'w' ).write( '\n'.join( output ) )
        
        for fastaID in sequenceIDs:
            sequence = fasta[fastaID]
    
    @staticmethod
    def run_immugen( inputFasta, outputDir, rawFlag ):
        command = IMGEN_PATH
        rawOutput = os.path.join( outputDir, "IMGEN" )
        if not os.path.exists( rawOutput ):
            os.mkdir( rawOutput )
        
        fasta = readfasta.readFasta( inputFasta, key="full" )
        sequenceIDs = readfasta.readFastaDesc( inputFasta, key="full" )
        
        tmpInput = os.path.join( rawOutput, "sequence.tmp" )
        open( tmpInput, 'w' ).write( '' )
        for seq in fasta.values():
            open( tmpInput, 'a' ).write( "%s\n" % seq )
        
        os.system( "python2.7 %s %s > %s" % (command, tmpInput, os.path.join( rawOutput, "%s.output" % Path( inputFasta ).stem ) ) )
        
        output = ['ID\tImmunogenicity_Score']
        values = {}
        for row in open( os.path.join( rawOutput, "%s.output" % Path( inputFasta ).stem ) ).read().split( '\n' )[4:]:
            if len( row ) == 0:
                continue
            tokens = row.split( ',' )
            seq = tokens[0]
            fastaID = list( fasta.keys() )[list( fasta.values() ).index( seq )]
            values[fastaID] = tokens[2]
        for fastaID in sequenceIDs:
            if fastaID in values.keys():
                output.append( '\t'.join( [fastaID, values[fastaID]] ) )
            else:
                output.append( fastaID + '\t0.0' )
        open( os.path.join( rawOutput, "%s.imgen.tsv" % Path( inputFasta ).stem ), 'w' ).write( '\n'.join( output ) )
        if not rawFlag:
            os.remove( os.path.join( rawOutput, "%s.output" % Path( inputFasta ).stem ) )
        os.remove( tmpInput )
    
    def check_args( self, args ):
        if args.method.lower() not in ["all","psortb","spaan","signalp","tmhmm","mdesc", "imgen"]:
            sys.stderr.write( "Incorrect feature method! Please choose from: [psortb,spaan,signalp,tmhmm,mdesc,imgen]\n" )
            exit(1)
        if not os.path.exists( args.inputFasta ):
            sys.stderr.write( "Input fasta not found! Please check your input file.\n" )
            exit(1)
        elif not any( SeqIO.parse( open( args.inputFasta, 'r' ), "fasta" ) ):
            sys.stderr.write( "Invalid input fasta! Please check your input file is in fasta format.\n" )
            exit(1)
        if not os.path.exists( args.outputDir ):
            os.mkdir( args.outputDir )
        if args.organism.lower() not in ["gram+","g+"]:
            sys.stderr.write( "Incorrect organism! Please choose from: [gram+]\n" )
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
    
    @staticmethod
    def mp_run( clines, poolsize ):
        # (c) The Scottish Crop Research Institute 2010
        # Author: Leighton Pritchard
        #
        # Contact:
        # leighton.pritchard@hutton.ac.uk
        #
        # Leighton Pritchard,
        # Information and Computing Sciences,
        # James Hutton Institute,
        # Errol Road,
        # Invergowrie,
        # Dundee,
        # DD6 9LH,
        # Scotland,
        # UK
        #
        # The MIT License
        #
        # Copyright (c) 2010-2014 The James Hutton Institute
        #
        # Permission is hereby granted, free of charge, to any person obtaining a copy
        # of this software and associated documentation files (the "Software"), to deal
        # in the Software without restriction, including without limitation the rights
        # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
        # copies of the Software, and to permit persons to whom the Software is
        # furnished to do so, subject to the following conditions:
        #
        # The above copyright notice and this permission notice shall be included in
        # all copies or substantial portions of the Software.
        #
        # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
        # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
        # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
        # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
        # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
        # OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
        # THE SOFTWARE.
        pool = multiprocessing.Pool(processes=poolsize)
        complete = []
        pool_output = [pool.apply_async(subprocess.call,
                                        (str(cline), ),
                                        {'stderr': subprocess.PIPE,
                                         'shell': sys.platform != "win32"})
                       for cline in clines]
        pool.close()
        pool.join()

    def main( self ):
        try:
            parser = argparse.ArgumentParser( description="Compute feature from sequences in fasta format for Vaxign Prediction." )
            parser.add_argument( "--feature", '-f', dest='method', required=True )
            parser.add_argument( "--input", '-i', dest='inputFasta', required=True )
            parser.add_argument( "--output", '-o', dest='outputDir', required=True )
            parser.add_argument( "--organismtype", '-t', dest='organism', required=True)
            parser.add_argument( "--multi", '-m', dest='multiFlag', default="True" )
            parser.add_argument( "--process", '-p', dest='process', type=int, default=int(multiprocessing.cpu_count()/2) )
            parser.add_argument( "--raw", '-r', dest='rawFlag', default="False" )
            args = parser.parse_args()
            
            self.check_args( args )
            
            if args.method.lower() == "all" or args.method.lower() == "psortb":
                Feature.run_psortb( args.inputFasta, args.outputDir, args.organism, args.multiFlag, args.process, args.rawFlag )
            if args.method.lower() == "all" or args.method.lower() == "spaan":
                Feature.run_spaan( args.inputFasta, args.outputDir, args.rawFlag )
            if args.method.lower() == "all" or args.method.lower() == "signalp":
                Feature.run_signalp( args.inputFasta, args.outputDir, args.organism, args.multiFlag, args.process, args.rawFlag )
            if args.method.lower() == "all" or args.method.lower() == "tmhmm":
                Feature.run_tmhmm( args.inputFasta, args.outputDir, args.multiFlag, args.process, args.rawFlag )
            if args.method.lower() == "all" or args.method.lower() == "mdesc":
                Feature.run_descriptor( args.inputFasta, args.outputDir, args.rawFlag )
            if args.method.lower() == "all" or args.method.lower() == "imgen":
                Feature.run_immugen( args.inputFasta, args.outputDir, args.rawFlag )
        except:
            print( sys.exc_info() )
        

if __name__ == "__main__":
    Feature().main()