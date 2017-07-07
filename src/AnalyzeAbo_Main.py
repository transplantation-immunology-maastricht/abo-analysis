# This file is part of abo-analysis.
#
# abo-analysis is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# abo-analysis is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with abo-analysis. If not, see <http://www.gnu.org/licenses/>.

SoftwareVersion = "abo-analysis Version 1.0"

import os

import sys
import getopt

from AnalyzeAbo import *


def usage():
    print("usage:\n" + 
    "\tThis script is written for python 2.7.11\n" + 
    "\tI haven't written the usage tutorial yet.  Oops.  Do this now please."
    )      
    
    
# Read Commandline Arguments.  Return true if everything looks okay for read extraction.
def readArgs():
    # Default to None.  So I can easily check if they were not passed in.
    global referenceFileName    
    global allelesFileName
    global allelesAFileName
    global allelesBFileName
    global allelesOFileName
    global outputDirectoryName

    referenceFileName       = None
    allelesFileName         = None
    allelesAFileName        = None
    allelesBFileName        = None
    allelesOFileName        = None
    outputDirectoryName     = None

    if(len(sys.argv) < 3):
        print ('I don\'t think you have enough arguments.\n')
        usage()
        raise()
        return False    

    # getopt.getopt(..) is a function for parsing args the way smart people do it
    # For More info, use google or 
    # https://www.tutorialspoint.com/python/python_command_line_arguments.htm
    try:
        opts, args = getopt.getopt(sys.argv[1:]
            ,"hvr:a:o:A:B:O:"
            ,["help", "version", "reference=","alleles=","output=","alleles-a=","alleles-b=","alleles-o="])

        for opt, arg in opts:

            if opt in ('-h', '--help'):
                print (SoftwareVersion)
                usage()
                return False

            elif opt in ('-v', '--version'):
                print (SoftwareVersion)
                return False

            elif opt in ("-i", "--alleles"):
                allelesFileName = arg
                
            elif opt in ("-i", "--reference"):
                referenceFileName = arg
                
            elif opt in ("-o", "--output"):
                outputDirectoryName = arg

            elif opt in ("-A", "--alleles-a"):
                allelesAFileName = arg
                
            elif opt in ("-B", "--alleles-b"):
                allelesBFileName = arg
                
            elif opt in ("-O", "--alleles-o"):
                allelesOFileName = arg

            else:
                print('Unknown commandline option: ' + opt)
                raise()

    except getopt.GetoptError, errorMessage:
        print ('Something seems wrong with your commandline parameters.')
        print (errorMessage)
        usage()
        return False

    print('Reference File:' + str(referenceFileName))
    print('Alleles File:' + str(allelesFileName))
    print('Output Directory:' + str(outputDirectoryName))

    # Quick sanity check.
    if(len(referenceFileName) < 4):
        print('referenceFileName is too short:' + str(referenceFileName))
        return False
    
    #if(len(allelesFileName) < 4):
    #    print('allelesFileName is too short:' + str(allelesFileName))
    #    return False
    if(len(outputDirectoryName) < 4):
        print('Output directory is too short:' + str(outputDirectoryName))
        return False
        
    if not os.path.isdir(outputDirectoryName):
        os.mkdir(outputDirectoryName)

    return True
    

if __name__=='__main__':

    try:    
        if(readArgs()):
            
            print('Starting to Analyze ABO')
            
            #compareIndividualABOAllelesToReference(referenceFileName, allelesFileName, outputDirectoryName)
            batchABOAllelesAgainstReference(referenceFileName, allelesFileName, outputDirectoryName)
            
           
            print ('Done analyzing ABO.')    
        else:
            print('\nI\'m giving up because I was not satisfied with your commandline arguments.')  
            
    except Exception:
        # Top Level exception handling like a pro.
        # This is not really doing anything.
        print 'Fatal problem during read extraction:'
        print sys.exc_info()
        raise


