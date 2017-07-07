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

import pysam

from AnalyzeAbo_Common import *
from PairwiseAlignmentResult import *

from os.path import join, isdir
from os import makedirs

def compareIndividualABOAllelesToReference(referenceFileName, allelesFileName, outputDirectoryName):
    print('I load the ABO Reference first.')
    
    alignmentInfoOutputFile = createOutputFile(join(outputDirectoryName, 'AlleleAlignments.txt'))
    alignmentDataSpreadsheet = createOutputFile(join(outputDirectoryName, 'AlleleData.csv'))
    # Headers on my spreadsheet
    alignmentDataSpreadsheet.write('AlleleName,AlignmentType,Begin,End,Length\n')
    
    for recordIndex, record in loadInputRecords(referenceFileName):
        referenceSequenceID = str(record.id)
        referenceSequence = str(record.seq)
        
    print('The reference sequence is called:' + str(referenceSequenceID))
    
    print('Then load the ABO allele sequences.')
    alleleSequences = loadInputRecords(allelesFileName)
    
    #print ('Guess how many alleles, too late it is this many:' + str(len(alleleSequences)))
    
    
    for index, record in alleleSequences:
        currentReadID = str(record.id)
        currentSequence = str(record.seq)

            
        sequenceLength = len(currentSequence)
     

        aboPhenotype=getPhenotype(currentReadID)
        
        print ('Analyzing allele:' + currentReadID)
        alignmentInfoOutputFile.write('>' + str(currentReadID) + '\n')
        alignmentInfoOutputFile.write('Phenotype=' + str(aboPhenotype) + '\n')

        # Align, find interesting polymorphisms.
        alleleAlignment = PairwiseAlignmentResult(currentReadID, referenceSequence , currentSequence)
        
        alignmentInfoOutputFile.write('Alignment Score:' + str(alleleAlignment.alignmentScore) + '\n')
        alignmentInfoOutputFile.write(str(alleleAlignment.sequence1Aligned) + '\n')
        alignmentInfoOutputFile.write(str(alleleAlignment.sequence2Aligned) + '\n\n')
        
        alignmentDataSpreadsheet.write(str(currentReadID) + ': ' + str(sequenceLength) + '\n')
        
        # loop through alignment tuples, print em to a spreadsheet.
        for alignmentTuple in alleleAlignment.allTuples:
            regionLength = int(alignmentTuple[2]) - int(alignmentTuple[1]) + 1 
            if (regionLength > 0):
                alignmentDataSpreadsheet.write(',' + str(alignmentTuple[0])
                    + ',' + str(alignmentTuple[1])
                    + ',' + str(alignmentTuple[2])
                    + ',' + str(regionLength)             
                    + '\n')
        

    alignmentInfoOutputFile.close()
    alignmentDataSpreadsheet.close()
    
    #return sequenceList
    
# Perform BW Alignment.  Align all reads against the Reference.

def batchABOAllelesAgainstReference(referenceFileName, allelesFileName, outputDirectoryName):

    
    groupwiseAlignmentResultsFile = createOutputFile(join(outputDirectoryName, 'GroupwiseABOAlignmentResults.txt'))
    alignmentDataSpreadsheet = createOutputFile(join(outputDirectoryName, 'GroupwiseABOAlignmentData.csv'))
    
    alignmentSubdir = join(outputDirectoryName, 'AllAlleleAlignment')
    alignReads(referenceFileName, allelesFileName, alignmentSubdir)
    analyzeAlignment(alignmentSubdir)
    
    groupwiseAlignmentResultsFile.close()
    alignmentDataSpreadsheet.close()
    


def alignReads(referenceLocation, alleleFileLocation, resultsOutputDir):
    # Align Alleles against Reference    
    if not isdir(resultsOutputDir):
        makedirs(resultsOutputDir)
    
    # Part 1 Index the Reference        
    try:
        # Copy the reference sequence to the alignment directory. This is a complicated way to do it.
        newReferenceLocation = join(resultsOutputDir,'AlignmentReference.fasta')
        refSequence = list(SeqIO.parse(referenceLocation, 'fasta'))[0]
        refSequence.id = 'AlignmentReference'
        sequenceWriter = createOutputFile(newReferenceLocation)
        SeqIO.write([refSequence], sequenceWriter, 'fasta')
        sequenceWriter.close()
                    
        # Index The Reference
        cmd = ("bwa index " + newReferenceLocation)
        os.system(cmd)
        
    except Exception:
        print ('Exception indexing alignment reference. Is bwa installed? folder writing permission issue?')                  
        raise 
    
    # Part 2 Align
    try:
        # align | sam->bam | sort
        tempAlignmentName = join(resultsOutputDir,'alignment')
        bwaMemArgs = "-t 4 -x ont2d"
        cmd = ("bwa mem " + 
            bwaMemArgs + " " +  
            newReferenceLocation + " " +
            alleleFileLocation + 
            " | samtools view  -Sb - | samtools sort - "
            + tempAlignmentName)
        #print ('alignment command:\n' + cmd)
        os.system(cmd)
        alignmentOutputName = tempAlignmentName + '.bam'
        
    except Exception:
        print ('Exception aligning reads against reference. Are bwa and samtools installed?')                  
        raise 
    
    # Part 3 Index Alignment
    try:
        cmd = ("samtools index " + alignmentOutputName)
        #print ('alignment index command:\n' + cmd)
        os.system(cmd)
        #print ('index command:\n' + cmd)
    except Exception:
        print ('Exception indexing alignment reference. Is bwa installed?')                  
        raise 



def analyzeAlignment(alignmentOutputDirectory):#, totalReadCount):
    
    # Ideas regarding ABO:
    # Parse the entire length of the alignment
    # Huge spreadsheet, headers:
    # [RefPosition], [RefNucleotide] [AlignedSequenceCount] [MostCommonABase] [MostCommonBBase] [MostCommonOBase]...
    # [AMatchPercent] [AMismatchPercent] [AInsertionPercent] [ADeletionPercent]
    # [BMatchPercent] [BMismatchPercent] [BInsertionPercent] [BDeletionPercent]
    # [OMatchPercent] [OMismatchPercent] [OInsertionPercent] [ODeletionPercent]#
    # 17 columns so far. I can also summarize this but better start with the raw.

    print ('\nStep 2.) Parse the alignment and create a new consensus sequence.')
    
    # Load up the Alignment Reference file, we'll need it.
    alignmentReferenceFileName = join(alignmentOutputDirectory,'AlignmentReference.fasta')
    alignmentRef = list(SeqIO.parse(alignmentReferenceFileName, 'fasta'))[0]
    
    # Count the reads in the input file
    #totalReadCount = len(list(SeqIO.parse(self.readInput, self.readInputFormat)))
    #self.readInputFormat
    #self.readInput
            
    # We generate a new consensus sequence from the alignment results.
    #newConsensusSequence = ""
    
    # Open the bam file
    bamfile = pysam.AlignmentFile(join(alignmentOutputDirectory,'alignment.bam'), 'rb')  
    
    # Open alignment analysis text file
    alignmentSummaryFile = createOutputFile(join(alignmentOutputDirectory,'AlignmentSpreadsheet.csv')) 
    alignmentSummaryFile.write('Ref_Position,Ref_Base,etc\n')
    
    # A smaller log. I will provide human-readable descriptions of the
    # bases that were adjusted in the new consensus sequence.
    # TODO: Provide surrounding sequence as well, maybe it's a repeat region....
    # Acutally NAH, I want to just put it in the wrangler log. 
    #adjustedBasesSummaryFile = createOutputFile(join(alignmentOutputDirectory,'AdjustedBases.txt')) 
    
    # Keep a running total of adjustments made to the reference.
    # If this total is 0, then theoretically the consensus matches the alignment reference, and we're done.
    #totalSequenceAdjustments = 0
    
    # Iterate the reference sequence column by column.
    pileupIterator = bamfile.pileup(alignmentRef.id)
    # TODO: Check where this pileup iterator starts. Are there reads mapped before or after the reference begins/ends?
    for pileupColumn in pileupIterator:
        
        #
        referencePosition = 0
        referenceBase = ''
        referenceAdjustment = '?'
        alignedCount = 0
        #unalignedCount = 0
        #matchCount = 0
        #mismatchCount = 0
        #inCount = 0
        #delCount = 0
        #aCount = 0
        #gCount = 0
        #cCount = 0
        #tCount = 0
        
        referencePosition = pileupColumn.reference_pos
        referenceBase = alignmentRef[pileupColumn.reference_pos].upper()
        alignedCount = pileupColumn.nsegments
        #unalignedCount = totalReadCount - alignedCount
        
        # Iterate the Reads at this position           
        for pileupRead in pileupColumn.pileups:
            
            # wait, what is the sequene name here?
            queryAlleleName = pileupRead.alignment.query_name
            
            currentBase = pileupRead.alignment.query_sequence[pileupRead.query_position].upper()
            bloodGroup = getPhenotype(queryAlleleName)
            
            print ('Index,AlleleName,Bloodgroup,RefBase,CurrentBase:' 
                   + str(referencePosition)
                   + ',' + queryAlleleName 
                   + ',' + bloodGroup
                   + ',' + referenceBase 
                   + ',' + currentBase
                   
                   
                   
                   
                   )
            
            # If this read is a deletion
            #if(pileupRead.is_del == 1):
            #    pass
            #    delCount += 1
            # else if this read is an insertion
            #elif(pileupRead.indel > 0):
            #    pass
                #print ('INSERTION DETECTED, INDEL=' + str(pileupRead.indel))  
            #    inCount += 1                   
            # Else if it is a refskip (TODO What does this mean? no read aligned? Count these?)
            if(pileupRead.is_refskip):
                print('This read is a refskip, i dont know what that means:' + pileupRead.alignment.query_name)
                raise Exception('This read is a refskip, i dont know what that means:' + pileupRead.alignment.query_name)
            # else this means we have a base aligned at this position for this read.
            #else:    
            #    currentBase = pileupRead.alignment.query_sequence[pileupRead.query_position].upper()                    
                #print('Reference,Current:' + referenceBase + ',' + currentBase)
                #print('Curr')
            #    if(currentBase == referenceBase):
            #        matchCount += 1
             #   else:
            #        mismatchCount += 1

            # Write it to my csv file
            #[RefPosition], [RefNucleotide] [AlignedSequenceCount] [MostCommonABase] [MostCommonBBase] [MostCommonOBase]...
    # 
        alignmentSummaryFile.write(
            str(referencePosition )
            + ',' + str(referenceBase)
            + ',' + str(alignedCount)

            + '\n')
            
            
    alignmentSummaryFile.close()


"""               
            # Count the nucleotide 
            #if (currentBase == 'A'):
            #    aCount += 1
            #elif (currentBase == 'G'):
            #    gCount += 1
            #elif (currentBase == 'C'):
            #    cCount += 1
            #elif (currentBase == 'T'):
             #   tCount += 1
            #else:
            #    print('Unknown Base found in Alignment at position ' + str(referencePosition) + ':' + currentBase)
            #    raise Exception('Unknown Base in Alignment')
            
            
            # TODO: What if the query insertion sequence is longer than one base?
            # Maybe I can only adjust one base per iteration, is that okay? Probably for the Best, actually..
            # Don't worry bout it for now.
        
        # Calculate highest frequency base
        # I hope this algorithm makes sense, probably there is a smarter way to do it.
        #if(aCount >= gCount and aCount >= cCount and aCount >= tCount):
        #    mostFrequentBase = 'A'
        #    mostFrequentBaseCount = aCount
        #elif(gCount >= cCount and gCount >= tCount):
        #    mostFrequentBase = 'G'
        #    mostFrequentBaseCount = gCount
        #elif(cCount >= tCount):
        #    mostFrequentBase = 'C'
        #    mostFrequentBaseCount = cCount
        #else:
        #    mostFrequentBase = 'T'
        #    mostFrequentBaseCount = tCount
        
        #TODO: Detect heterozygosity here
        # Do the base frequencies look "normal"?
        # High proportion of Inserts or Deletions?
        # Maybe I don't care, because I want to build a read-clusterer tool.
        
        
        # Add the next base to the new consensus sequence            
        #if (matchCount >= mismatchCount and matchCount >= inCount and matchCount >= delCount):
            # Aligned bases match the reference, add reference base to the consensus.
        #    referenceAdjustment='-'
        #    newConsensusSequence += referenceBase
            
        #elif (inCount >= mismatchCount and inCount >= delCount):
            # Aligned bases show an insertion.
            # Add the Reference Base and the Insertion Base to the consensus.  
        #    totalSequenceAdjustments += 1 
        #    referenceAdjustment='I'  
        #    newConsensusSequence += referenceBase + mostFrequentBase         
            
        #    self.wranglerLog.write(str(referencePosition) + ':Insertion' +
        #        '\n(' + str(inCount) + '/' + str(alignedCount) + ') = ' + str((100.0 * inCount) / alignedCount) + '% of aligned reads'
        #        '\n(' + referenceBase + ' > ' + referenceBase + mostFrequentBase + ')' +
        #        '\n')
            
            #TODO: I need to insert multiple bases, if that is waht the alignment suggests.

        #elif (delCount >= mismatchCount):
            # Reads show a deletion.
            # Don't add anything to the consensus.
        #    totalSequenceAdjustments += 1
        #    referenceAdjustment='D'
            
         ##   self.wranglerLog.write(str(referencePosition) + ':Deletion' +
          #      '\n(' + str(delCount) + '/' + str(alignedCount) + ') = ' + str((100.0 * delCount) / alignedCount) + '% of aligned reads'
          #      '\n(' + referenceBase + ' > _)' +
          #      '\n')
            
        #else:
            # Mismatch base.
            # Add the highest read count base to the reference.
            # It might actually be the same base as the reference,
            # Because this just means there are more mismatches than matches.
            # Problematic base, at least we'll notice here.
            # TODO: What to do with highly heterozygous Positions?
            # I should report those that look particularly heterozygous, somewhere.
            newConsensusSequence += mostFrequentBase 
            totalSequenceAdjustments += 1     
            referenceAdjustment='M'   
            
            self.wranglerLog.write(str(referencePosition) + ':Mismatch' +
                '\n(' + str(mostFrequentBaseCount) + '/' + str(alignedCount) + ') = ' + str((100.0 * mostFrequentBaseCount) / alignedCount) + '% of aligned reads'
                '\n(' + referenceBase + ' > ' + mostFrequentBase + ')' +
                '\n')
          

        # Write a line to the alignment Summary 
        alignmentSummaryFile.write(str(referencePosition) + 
            ',' + str(referenceBase) +
            ',' + str(referenceAdjustment) + 
            ',' + str(alignedCount) + 
            ',' + str(unalignedCount) + 
            ',' + str(matchCount) + 
            ',' + str(mismatchCount) + 
            ',' + str(inCount) + 
            ',' + str(delCount) + 
            ',' + str(aCount) + 
            ',' + str(gCount) + 
            ',' + str(cCount) + 
            ',' + str(tCount) +
            '\n')
        
    print('\nTotal Sequence Adjustments:' + str(totalSequenceAdjustments) + ' (How many bases the consensus differs from the reference.)\n')    
    
    # Write the newly constructed consensus sequence.
    currentConsensusSequenceFileName = join(alignmentOutputDirectory, 'Consensus.fasta')        
    consensusWriter = createOutputFile(currentConsensusSequenceFileName)          
       
    SeqIO.write([SeqRecord(Seq(newConsensusSequence,
        IUPAC.unambiguous_dna),
        id="GeneratedConsensusSequence|Coverage=GarbageInformation", description="") ], consensusWriter, 'fasta')
    consensusWriter.close()
        
    self.wranglerLog.write('Total Sequence Adjustments:' + str(totalSequenceAdjustments) + '\n')
        """
    # Close Summary Files
    
    
    #adjustedBasesSummaryFile.close()
    
    #return totalSequenceAdjustments    

    
    
    
            
def getPhenotype(currentSeqID):

    if(currentSeqID.startswith('ABO')):
        #print('This read starts with ABO.')
        pass
    else:
        print('doesnt start with abo.')
        raise Exception ('I want the sequence IDs to start with ABO')

    currentPhenotype = '??'

    dividerCharacter = currentSeqID[3]
    #print('Divider = ' + dividerCharacter)
    
    # Strange logic to decide the ABO types
    if(currentSeqID[4] == 'O'):
        currentPhenotype = 'O'
        
    elif(currentSeqID[4] == 'A'):
        currentPhenotype = 'A'
        
    elif(currentSeqID[4] == 'B'):
        currentPhenotype = 'B'
        
    else:
        print('I don\'t know what phenotype is this allele:' + str(currentSeqID))
        raise Exception ('I want the sequence IDs to start with ABO')

    return currentPhenotype
    
def scaleValues(inputValues):
    # scale from 0:1
    scaledData = []
    
    oldMinimum = min(inputValues)
    oldMaximum = max(inputValues)

    # Calculated based on a formula from here:
    # http://stats.stackexchange.com/questions/25894/changing-the-scale-of-a-variable-to-0-100

    if(oldMaximum == oldMinimum):
        return inputValues
    else:
        for currentValue in inputValues:
            #print('x is:' + str(x))
            scaledData.append(((1 - 0)/(oldMaximum - oldMinimum)) * (currentValue - oldMaximum) + 1)
        
    return scaledData

