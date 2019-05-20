#!/usr/local/bin/python2.7.3


from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import numpy as np
import re
import csv
import os
import sys
import os.path
import tempfile
import cgi
import cgitb
import time
import copy

"""Converters and Checkers"""
def ConvertInt(Chr):
    if Chr == "X":
        Chr = 0
    if Chr == "Y":
        Chr = -1
    Chr = int(Chr)
    return Chr

def ConvertStr(Chr):
    if Chr == 0:
        Chr = "X"
    if Chr == -1:
        Chr = "Y"
    Chr = str(Chr)
    return Chr

def NumAlphaConvert(n): # converts numeral to alpha
    out = []
    if n == 0:
        out = 'A'
    elif n == 1:
        out = 'B'
    elif n == 2:
        out = 'C'
    elif n == 3:
        out = 'D'
    elif n == 4:
        out = 'E'
    elif n == 5:
        out = 'F'
    elif n == 6:
        out = 'G'
    elif n == 7:
        out = 'H'
    elif n == 8:
        out = 'I'
    elif n == 9:
        out = 'J'
    elif n == 10:
        out = 'K'
    elif n == 11:
        out = 'L'
    elif n == 12:
        out = 'M'
    elif n == 13:
        out = 'N'
    elif n == 14:
        out = 'O'
    elif n == 15:
        out = 'P'   
    else:
        out = 'error: too many rearrangements'
    return out

def CheckFASTAString(strg, search=re.compile(r'[^a,^t,^c,^g,^A,^T,^C,^G,^n,^N,^-]').search):
    return not bool(search(strg))

"""End Converters and Checkers"""

"""Genomic Data"""
assembly = "GRCh37"
def DefCentromeres(): # retrieve the centromere position
    #centromere_file = os.path.join(GetAppDir(),"centromere2.txt")
    centromere_file = "centromere_%s.txt" % assembly
    f = open(centromere_file, "r") # currently a static file that needs to be updated with each build
    centromeres = np.loadtxt(centromere_file, usecols = (2,3))
    f.close()
    centromeres = centromeres.tolist()
    return centromeres

def DefBands(): #defines breakpoints cytogenetic bands
    Mbn_file = os.path.join(GetAppDir(),"Mbn.csv")
    f = open(Mbn_file, "rb") #currently a static file that needs to be updated with each build
    reader = csv.reader(f)
    bands = []
    for row in reader:
        bands.append(row)
    f.close()
    return bands
"""End Genomic Data"""

"""Parameters"""
def EndFilter():
    endfilter = 5 #ensures alignment is near the end of the FASTA rather than a homologous region in center
    return endfilter

def AssemblyFilter():
    assemblyfilter = "Primary Assembly"
    return assemblyfilter

def EValueThreshold():
    evaluethreshold = 0.05
    return evaluethreshold

def OverlapCheck():
    overlapcheck = 10
    return overlapcheck

def Identity():
    identity = 0.90
    return identity
"""End Parameters"""

"""Retrieving Input"""
def GetAppDir():
    return os.path.dirname(os.path.realpath(__file__))

def read_fasta_input(fasta_input_file):
    arg_len = len(fasta_input_file)
    if arg_len < 2:
        print 'Usage: python(2.7+)', sys.argv[0], '[path_to_fasta_file]'
        sys.exit()
    path_to_fasta_file = fasta_input_file
    if not os.path.isfile(path_to_fasta_file):
        print 'Fasta input file "', path_to_fasta_file, '" does not exist.'
        sys.exit()
    return path_to_fasta_file


def read_fasta_input_from_web(form):
    fasta_input_file = tempfile.NamedTemporaryFile(delete=True).name
    open(fasta_input_file,'w').write(form['fasta_input'].value)
    return fasta_input_file

def GetInputFASTA(fasta_input_file):
    FASTAStartlines = []
    TempInputSeqs = []
    InputSeqs = []
    FullSeq = []
    InputError = []
    FASTASeqs = open(fasta_input_file,'r')
    inputlines = FASTASeqs.readlines()
    print "here next"
    for i in range(len(inputlines)):
        if re.search('>',inputlines[i]) != None:
            FASTAStartlines.append(i)
        print i
    for i in range(len(FASTAStartlines)):
        if i == (len(FASTAStartlines)-1):
            elem = FASTAStartlines[i]
            stripelem = inputlines[elem].strip('\n')
            InputSeqs.append(stripelem)
            FullSeq = []
            for j in range(FASTAStartlines[i],len(inputlines)):
                if j == (len(inputlines)-1):
                    FullSeq = ''.join(FullSeq)
                else:
                    FullSeq += inputlines[j+1]
                    if re.search('\n', FullSeq[-1]) != None:
                        FullSeq.remove('\n')
                    if re.search('\r', FullSeq[-1]) != None:
                        FullSeq.remove('\r')
            InputSeqs.append(FullSeq)
        else:
            elem = FASTAStartlines[i]
            stripelem = inputlines[elem].strip('\n')
            InputSeqs.append(stripelem)
            FullSeq = []
            for j in range(FASTAStartlines[i],FASTAStartlines[i+1]):
                if j == (FASTAStartlines[i+1]-1):
                    FullSeq = ''.join(FullSeq)
                else:
                    FullSeq += inputlines[j+1]
                    if re.search('\n', FullSeq[-1]) != None:
                        FullSeq.remove('\n')
                    if re.search('\r', FullSeq[-1]) != None:
                        FullSeq.remove('\r')
            InputSeqs.append(FullSeq)
    if InputSeqs == []:
        InputError.append("Input Error IE001 - No FASTA input has been entered. Please check formating. Note FASTA should start with \">\"." )
    for i in range(len(InputSeqs)):
        if i % 2 == 0:
            pass
        else:
            if CheckFASTAString(InputSeqs[i]) == False:
                InputError.append("Input Error IE002 - Fasta Input:\" %s\" - FASTA sequence contains non-FASTA characters. Please check input." % InputSeqs[i-1])
    print "here third"
    print "Inputs", InputSeqs, InputError
    return InputSeqs, InputError
"""End Retrieving Input"""


""" BLAST I/O and BioPython Parsing"""
def GetXmlBlastRecords(fasta_file): #module to handle I/O to NCBI Blast and return XML output
    fasta_input_file = fasta_file
    XmlBlastRecords = tempfile.NamedTemporaryFile(delete=True).name
    db_file = os.path.join(GetAppDir(),"human_genomic.02")
    blastcmd = os.path.join(GetAppDir(), "blastn")
    blastn_cline = NcbiblastnCommandline(cmd=blastcmd, query=fasta_input_file, db=db_file, outfmt=5, evalue=0.001, out=XmlBlastRecords, task="megablast", max_target_seqs=25, best_hit_overhang=0.1)
    blastn_cline()
    return XmlBlastRecords

def ParseBlastXml(fasta_file): #parses xml data into individual records
    XmlBlastRecords = GetXmlBlastRecords(fasta_file)
    result_handle = open(XmlBlastRecords) #temporay file but will return XmlBlastRecords
    blast_records = NCBIXML.parse(result_handle)
    #blast_record = list(blast_records)
    #os.remove(XmlBlastRecords)
    return blast_records
""" End BLAST I/O and BioPython Parsing"""

"""Retrieve Necessary Data and Direction for Analysis from BioPython Parsed Records"""
def BlastOutput(blast_records): #retrieve necessary data from BioPython defined records   
    blast_output_error = []    
    blast_output= []
    QxyD = []
    quality_data=[]
    FASTAEnds = []
    endfilter = EndFilter()
    assemblyfilter = AssemblyFilter()
    evaluethreshold = EValueThreshold()
    identity = Identity()
    overlapcheck = OverlapCheck()
    for iteration in blast_records:
        SxC = SxBP = SxBPO = SxS = SxE = SxD = QxD = QxS = QxE = QxBP = QxBPO = xExpect = xScore = xIdentity = xSeq = []
        SyC = SyBP = SyBPO = SyS = SyE = SyD = QyD = QyS = QyE = QxBP = QxBPO = yExpect = yScore = yIdentity = ySeq = []
        endfilter = EndFilter()
        IQLength = iteration.query_length
        QScore = 1                    # QScore is score filter (bits) for BLAST Sequence - reset to zero for each alignment 
        for alignment in iteration.alignments:    # find best left(x) BLAST result per submission
            for hsp in alignment.hsps:
                if assemblyfilter in alignment.title:    # filters BLAST results to Primary Assembly
                    if hsp.query_start < endfilter:            # filters for left(x) BLAST results to only those that have a good reads to start of queary
                        QScoreTemp = hsp.expect           # set QScore to current record
                        if hsp.identities/float(hsp.align_length) > identity and QScore > QScoreTemp:            # check to see if current record the biggest QScore
                            QScore = QScoreTemp                #Set current record as biggest QScore
                            D = hsp.frame                #Build Direction vector D = (Q,S) e.g. (1,-1) == (+,-)
                            SxD    = D[1]                    #Set subject direction
                            QxD    = D[0]                    #set query direction
                            QxS = hsp.query_start            #set query start
                            SxS = hsp.sbjct_start            #set subject start
                            SxC = str(re.findall('chromosome [\dXY]+', alignment.title)) #find chromosome #
                            SxC = str(re.findall('[\dXY]+',    SxC)[0])        #refine chromosome number
                            if D == (1,1):                #define BLAST result end (breakpoint)
                                QxE = QxS + len(hsp.query) - 1
                                SxE = SxS + len(hsp.sbjct) - 1
                            elif D == (1,-1):
                                QxE = QxS + len(hsp.query) - 1
                                SxE = SxS - len(hsp.sbjct) + 1
                            elif D == (-1,1):
                                QxE = QxS - len(hsp.query) + 1
                                SxE = SxS + len(hsp.sbjct) - 1
                            elif D == (-1,-1):
                                QxE = QxS - len(hsp.query) + 1
                                SxE = SxS - len(hsp.sbjct) + 1
                            SxBP = SxE                              #making a BP definition for reference
                            SxBPO = SxS
                            QxBP = QxE
                            QxBPO = QxS
                            xExpect = hsp.expect            #adding additional information (not for filter)
                            xIdentity = (float(hsp.identities) / float(hsp.align_length))*100
                            xScore = hsp.bits
                            xSeq = str(hsp.query)
        QScore = 1   # QScore is score filter (bits) for BLAST Sequence - reset to zero for each alignment
        for alignment in iteration.alignments:    # find best right(y) BLAST result per submission
            for hsp in alignment.hsps:
                if assemblyfilter in alignment.title:    # filters BLAST results to Primary Assembly  
                    if hsp.query_end > (IQLength - endfilter):  # filters for right(y) BLAST results to only those that have a good reads to end of query
                        QScoreTemp = hsp.expect           # set QScore to current record
                        if hsp.identities/float(hsp.align_length) > identity and QScore > QScoreTemp:            # check to see if current record the biggest QScore
                            QScore = QScoreTemp                #Set current record as biggest QScore
                            D = hsp.frame                #Build Direction vector D = (Q,S) e.g. (1,-1) == (+,-)
                            SyD    = D[1]                    #Set subject direction
                            QyD    = D[0]                    #set query direction
                            QyS = hsp.query_start            #set query start
                            SyS = hsp.sbjct_start            #set subject start
                            SyC = str(re.findall('chromosome [\dXY]+', alignment.title)) #find chromosome #
                            SyC = str(re.findall('[\dXY]+',    SyC)[0])        #refine chromosome number
                            if D == (1,1):                #define BLAST result end (breakpoint)
                                QyE = QyS + len(hsp.query) - 1
                                SyE = SyS + len(hsp.sbjct) - 1
                            elif D == (1,-1):
                                QyE = QyS + len(hsp.query) - 1
                                SyE = SyS - len(hsp.sbjct) + 1
                            elif D == (-1,1):
                                QyE = QyS - len(hsp.query) + 1
                                SyE = SyS + len(hsp.sbjct) - 1
                            elif D == (-1,-1):
                                QyE = QyS - len(hsp.query) + 1
                                SyE = SyS - len(hsp.sbjct) + 1
                            SyBP = SyS                              #making a BP definition for reference
                            SyBPO = SyE
                            QyBP = QyS
                            QyBPO = QyE
                            yExpect = hsp.expect            #adding additional information (not for filter)
                            yIdentity = (float(hsp.identities) / float(hsp.align_length))*100
                            yScore = hsp.bits
                            ySeq = str(hsp.query)
        if SxC == [] and SyC == []:    #Error Handling for poor alignments
            blast_output_error.append(" BLAST Ouput Error BE003 - Fasta Input:\" %s\" - neither end aligned well to reference. Change filter parameter, remove poor alignment from ends (e.g., \"N \") and/or insert more flanking nucleotides." % iteration.query)   
        elif SxC == [] and SyC != []:
            blast_output_error.append(" BLAST Ouput Error BE004 - Fasta Input:\" %s\" - left end did not align well to reference. Change filter parameter, remove poor alignment from ends (e.g., \"N \") and/or insert more flanking nucleotides." % iteration.query) 
        elif SxC != [] and SyC == []:
            blast_output_error.append(" BLAST Ouput Error BE005 - Fasta Input:\" %s\" right end did not align well to reference. Change filter parameter, remove poor alignment from ends (e.g., \"N \") and/or insert more flanking nucleotides." % iteration.query) 
        elif SxS == SxBP and abs(SxE-SyBP) < overlapcheck: #Error Handling of non-breakpoint sequences
            blast_output_error.append(" BLAST Ouput Error BE001 - Fasta Input:\" %s\" - did not find breakpoint. Insert more flanking nucleotides." % iteration.query)
        elif SxE == SxBP and abs(SxS-SyBP) < overlapcheck:
            blast_output_error.append(" BLAST Ouput Error BE002 - Fasta Input:\" %s\" - did not find breakpoint. Insert more flanking nucleotides." % iteration.query)
        elif SxC != [] and SyC != []: #load data into string for each iteration
            blast_output.append([SxC, SxBP, SxBPO, SxS, SxE, SxD])
            blast_output.append([SyC, SyBP, SyBPO, SyS, SyE, SyD])
            QxyD.append([QxD, QxS, QxE,QxBP, QxBPO])
            QxyD.append([QyD, QyS, QyE,QyBP, QyBPO])
            quality_data.append([xExpect, xScore, xIdentity])
            quality_data.append([yExpect, yScore, yIdentity])
            FASTAEnds.append([xSeq])
            FASTAEnds.append([ySeq])
    return blast_output, QxyD, quality_data, FASTAEnds, blast_output_error

def GetArm(blast_output): #determines chromosome arm of breakpoints (i.e. p or q)
    blast_arm = range(0, len(blast_output))
    centromeres = DefCentromeres()
    for i in range(0, len(blast_output)):
        ChrNum = ConvertInt(blast_output[i][0])
        if blast_output[i][2] > centromeres[ChrNum + 1][1]:
            SxyA = 'q'
        elif blast_output[i][2] < centromeres[ChrNum + 1][1]:
            SxyA = 'p'
        blast_arm[i] = SxyA
    return blast_arm


def GetBands(blast_output):#determines cytogenetic band of breakpoints
    BANDxy = []
    bands = DefBands()
    for o in range(len(blast_output)):
        for i in range(len(bands)):
            if bands[i][0] == str(blast_output[o][0]) and blast_output[o][2] >= int(bands[i][1]) and blast_output[o][2] < int(bands[i][2]):
                blast_output[o][3]
                BANDxy.append(bands[i][3])
    return BANDxy
"""End Retrieve Necessary Data and Direction for Analysis"""

"""Combine Essential Data for BP Analysis"""
def BPEssentials(blast_output, blast_arm, QxyD, BANDxy, quality_data, FASTAEnds, M1Raw = "Nothing"): #creates a list for each breakpoint end
    blast_output2 = map(list, blast_output)
    blast_arm2 = map(list, blast_arm)
    M1Raw = list(range(0,len(blast_output2)))
    for i in range(0,len(blast_output2)):
        M1Raw[i] = blast_output2[i] + blast_arm2[i] + QxyD[i] + [BANDxy[i]] + quality_data[i] + FASTAEnds[i]
    for i in range(0,len(blast_output2)/2):
        tempx = [''.join(['I',str(i),'x'])]
        tempy = [''.join(['I',str(i),'y'])]
        M1Raw[i*2] = list(tempx) + M1Raw[i*2]
        M1Raw[i*2 + 1] = list(tempy) + M1Raw[i*2 + 1]
    return M1Raw
"""End Combine Essential Data for BP Analysis"""




"""Building Break Points"""
def D_Picker(M1Raw): #pick derivative chromosome
    TempDnC = []
    DnC = list(range(0,len(M1Raw)))
    for i in range(0, len(M1Raw)):
        DnC[i] = M1Raw[i][1]

        DnC[i] = ConvertInt(DnC[i])
    DnC = list(set(DnC))
    DnC.sort()
    for elem in DnC:
        elem = ConvertStr(elem)
        TempDnC.append(elem)   
    DnC = TempDnC
    return DnC

def O_Picker(M1Raw,InputSeqs):
    DnO = []
    M1Final = copy.deepcopy(M1Raw)
    for i in range(len(M1Raw)/2):   #Find overlap
        if M1Raw[i*2][8] == 1 and M1Raw[i*2 + 1][8] == 1:
            DnO.append([M1Raw[i*2][10] - M1Raw[i*2 + 1][9]])
        elif M1Raw[i*2][8] == 1 and M1Raw[i*2 + 1][7] == -1:
            DnO.append([M1Raw[i*2][10] - M1Raw[i*2 + 1][10]])
        elif M1Raw[i*2][8] == -1 and M1Raw[i*2 + 1][7] == 1:
            DnO.append([M1Raw[i*2][9] - M1Raw[i*2 + 1][9]])
        elif M1Raw[i*2][8] == -1 and M1Raw[i*2 + 1][7] == -1:
            DnO.append([M1Raw[i*2][9] - M1Raw[i*2 + 1][10]])
    for i in range(len(M1Raw)/2):    #Find BP overlap or insertion or normal
        if DnO[i][0] == -1:
            DnO[i].extend(['Balanced'])
            M1Final[i*2].extend([''])
            M1Final[i*2+1].extend([''])
        elif DnO[i][0] < -1:
            DnO[i].extend(['Gain'])
            if M1Raw[i*2][8] == 1:
                DnO[i].extend([M1Raw[i*2][10]])
            if M1Raw[i*2][8] == -1:
                DnO[i].extend([M1Raw[i*2][9]])
            if M1Raw[i*2 + 1][8] == 1:
                DnO[i].extend([M1Raw[i*2 + 1][9]])
            if M1Raw[i*2 + 1][8] == -1:
                DnO[i].extend([M1Raw[i*2 + 1][10]]) 
            temp = list(InputSeqs[i+1])
            m = ''.join(temp[DnO[i][2]:DnO[i][3] - 1])
            DnO[i].extend([m])
            M1Final[i*2].extend([m])
            M1Final[i*2+1].extend([m])
        elif DnO[i][0] > -1:
            DnO[i].extend(['Loss'])
            if M1Raw[i*2][6] == 1:
                DnO[i].extend([M1Raw[i*2][2] - DnO[i][0], M1Raw[i*2][2]])
            elif M1Raw[i*2][6] == -1:
                DnO[i].extend([M1Raw[i*2][2] + DnO[i][0], M1Raw[i*2][2]])
            if M1Raw[i*2 + 1][6] == 1:
                DnO[i].extend([M1Raw[i*2 + 1][2], M1Raw[i*2 + 1][2] + DnO[i][0]])
            elif M1Raw[i*2 + 1][6] == -1:
                DnO[i].extend([M1Raw[i*2 + 1][2], M1Raw[i*2 + 1][2] - DnO[i][0]])
            if list(str(DnO[i][2])) == list(str(DnO[i][3])):         # what to write in for MH x portion
                temp = list(str(DnO[i][2]))
                temp[len(list(str(DnO[i][2]))) - 1] = ''.join(['{',temp[len(list(str(DnO[i][2]))) - 1],'}'])    
                temp = ''.join(temp)
                DnO[i].append(temp)
            else:
                counter = 0
                BPs = list(str(DnO[i][2]))
                BPe = list(str(DnO[i][3]))
                for l in range(len(BPs)):
                    if BPs[l] != BPe[l]:
                        counter += 1
                temp = len(BPs) - counter
                finalBPs = []
                finalBPe = []
                for l in range(temp , len(BPs)):
                    finalBPs.extend(BPs[l])
                    finalBPe.extend(BPe[l])
                finalBPs = ''.join(finalBPs)
                finalBPe = ''.join(finalBPe)
                finalBP = ''.join(['{', finalBPs, '-', finalBPe, '}'])
                final = BPs
                for j in range(counter):
                    del final[-1]
                final.append(finalBP)
                final = ''.join(final)
                DnO[i].append(final)
            if list(str(DnO[i][4])) == list(str(DnO[i][5])):         # what to write in for MH y portion
                temp = list(str(DnO[i][4]))
                temp[len(list(str(DnO[i][4]))) - 1] = ''.join(['{',temp[len(list(str(DnO[i][4]))) - 1],'}'])    
                temp = ''.join(temp)
                DnO[i].append(temp)
            else:
                counter = 0
                BPs = list(str(DnO[i][4]))
                BPe = list(str(DnO[i][5]))
                for l in range(len(BPs)):
                    if BPs[l] != BPe[l]:
                        counter += 1
                temp = len(BPs) - counter
                finalBPs = []
                finalBPe = []
                for l in range(temp , len(BPs)):
                    finalBPs.extend(BPs[l])
                    finalBPe.extend(BPe[l])
                finalBPs = ''.join(finalBPs)
                finalBPe = ''.join(finalBPe)
                finalBP = ''.join(['{', finalBPs, '-', finalBPe, '}'])
                final = BPs
                for k in range(counter):
                    del final[-1]
                final.append(finalBP)
                final = ''.join(final)
                DnO[i].append(final)
            M1Final[i*2][2] = DnO[i][6]
            M1Final[i*2].extend([''])
            M1Final[i*2+1][2] = DnO[i][7]
            M1Final[i*2+1].extend([''])
        InputSeqsName = ''.join(InputSeqs[i*2])
        M1Final[i*2].extend([InputSeqsName])
        M1Final[i*2+1].extend([InputSeqsName])
    for l in range(len(DnO)):
        if DnO[l][1] == 'OVERLAP':
            M1Final[l*2][2] = DnO[l][6]
            M1Final[l*2 + 1][2] = DnO[l][7]
    return DnO, M1Final

def Rearrangement_Pieces(M1Raw, DnO):
    RP = []
    for i in range(len(M1Raw)/2):
        temp1 = ''.join([str(M1Raw[i*2][1]), str(M1Raw[i*2][13])])
        if M1Raw[i*2][6] == 1:
            temp2 = '(+)'
        elif M1Raw[i*2][6] == -1:
            temp2 = '(-)'
        if DnO[i][1] == 'Loss':
            temp3 = ''.join(['(', str(DnO[i][6]), ')'])
        else:
            temp3 = ''.join(['(', str(M1Raw[i*2][2]), ')'])
        if DnO[i][1] == 'Loss' or DnO[i][1] == 'Balanced' :
            temp4 = ''.join(['::'])
        elif DnO[i][1] == 'Gain' :
            temp4 = ''.join(['::',DnO[i][4],'::'])
        temp5 = ''.join([str(M1Raw[i*2 + 1][1]), str(M1Raw[i*2 + 1][13])])
        if M1Raw[i*2 + 1][6] == 1:
            temp6 = '(+)'
        elif M1Raw[i*2 + 1][6] == -1:
            temp6 = '(-)'
        if DnO[i][1] == 'Loss':
            temp7 = ''.join(['(', str(DnO[i][7]), ')'])
        else:
            temp7 = ''.join(['(', str(M1Raw[i*2 + 1][2]), ')'])
        RP.extend([''.join([temp1,temp2,temp3,temp4,temp5,temp6,temp7])])
    return RP 
"""End Building Break Points"""

"""Full Nomenclature"""
def Derivative_Definition(M1Raw, M1Final, DnC):
    D = []
    whilecounter = 0
    centromeres = DefCentromeres()
    UsageCheck=[]
    nomenclature_error = []
    loop_error = []
    for DC in DnC:               # loop through each D chr
        Rp=DnR0pB=DnR0pC=DnR0pD=DnR0pA=DnR0pN=DnR0pFASTA=DnR0pBO=[]
        Rq=DnR0qB=DnR0qC=DnR0qD=DnR0qA=DnR0qN=DnR0qFASTA=DnR0qBO=[]
        D.append(DC)
        ########################################### P PART
        NC = []
        for j in range(0, len(M1Raw)):                # loop through each Ix & Iy data set   D chr
            if M1Raw[j][1] == DC and M1Raw[j][7] == 'p': #find any p side recombinations
                DC = ConvertInt(DC)
                NC.extend([M1Raw[j][0],abs(M1Raw[j][4] - centromeres[DC+1][0]),M1Raw[j][0],abs(M1Raw[j][5] - centromeres[DC+1][0])]) #find distance from centromere for query start and end
                DC = ConvertStr(DC)
        if NC == []:            #if no p side found than set Pname to None
            Pname = 'Empty'
        elif NC != []:
            NCmini = NC.index(min(NC))        #find p recombination closest to centromere
            for i in range(len(M1Raw)):        #loop through M1Raw
                if re.search('y',NC[NCmini-1]) != None and M1Raw[i][0] == NC[NCmini-1] and M1Raw[i][6] == 1:    # find right (y) BLAST with NC value, with the same name (e.g., IOy) and (+) direction - this is the normal part of the chromosome
                    DnR0pC = str(M1Raw[i-1][1])
                    NR0pC = M1Raw[i][1]
                    if M1Final[i][18] == '':
                        BPU = '::'
                    else:
                        BPU =''.join( ['::',M1Final[i][18],"::"])
                    DnR0pN = ''.join(["->",str(DnR0pC),str(M1Raw[i-1][13]),"(",str(M1Final[i-1][2]),")",BPU,str(NR0pC),str(M1Raw[i][13]),"(",str(M1Final[i][2]),")->"])         # add x pside BP information
                    DnR0pD = M1Raw[i-1][6]
                    DnR0pA = M1Raw[i-1][7]
                    DnR0pB = M1Raw[i-1][2]
                    DnR0pBO = M1Raw[i-1][3]
                    DnR0pFASTA = M1Final[i][19]  
                elif re.search('y',NC[NCmini-1]) != None and M1Raw[i][0] == NC[NCmini-1] and M1Raw[i][6] == -1:
                    nomenclature_error.append("Nomenclature Error NE001 - Right side of FASTA Input:\" %0s\" - Sequence is inverted and cannot be closest to the centromere. Please check for other breakpoints on the %1s side of chromosome %2s." % (M1Final[i][17], "pter", DC))
                elif re.search('x',NC[NCmini-1]) != None and M1Raw[i][0] == NC[NCmini-1] and M1Raw[i][6] == -1:  # find left (x) BLAST with NC value, with the same name (e.g., IOx) and (-) direction - this is the normal part of the chromosome
                    DnR0pC = M1Raw[i+1][1]
                    NR0pC = M1Raw[i][1]
                    if M1Final[i][18] == '':
                        BPU = '::'
                    else:
                        BPU =''.join( ['::',M1Final[i][18],"::"]) 
                    DnR0pN = ''.join(["->",str(DnR0pC),str(M1Raw[i+1][13]),"(",str(M1Final[i+1][2]),")",BPU,str(NR0pC),str(M1Raw[i][13]),"(",str(M1Final[i][2]),")->"])    # add y pside BP information
                    DnR0pD = M1Raw[i+1][6]
                    DnR0pA = M1Raw[i+1][7]
                    DnR0pB = M1Raw[i+1][2]
                    DnR0pBO = M1Raw[i+1][3]
                    DnR0pFASTA = M1Final[i][19]  
                elif re.search('x',NC[NCmini-1]) != None and M1Raw[i][0] == NC[NCmini-1] and M1Raw[i][6] == 1:
                    nomenclature_error.append("Nomenclature Error NE002 - Left side of FASTA Input:\" %0s\" - Sequence is inverted and cannot be closest to the centromere. Please check for other breakpoints on the %1s side of chromosome %2s." % (M1Final[i][17], "pter", DC))
        if DnR0pN != []:
            UsageCheck.append(DnR0pFASTA)
            Rp = [[DnR0pB, DnR0pC, DnR0pD, DnR0pA, DnR0pN,'pside',DnR0pFASTA, DnR0pBO]]
            Pname = None
            while Pname == None:
                whilecounter = whilecounter + 1
                if whilecounter > (len(M1Raw)/2):
                    loop_error.append("Loop Error LE001 - Breakpoints create loop, check more multiple homologous chromosome involvement and/or additional breakpoints.")
                    break
                elif whilecounter <= (len(M1Raw)/2):
                    Rp, Pname, DnRnpFASTA = P_Extender(M1Raw, M1Final, Rp)
                    if DnRnpFASTA != []:
                        UsageCheck.append(DnRnpFASTA) 
            D.append(Rp)
        elif DnR0pN == []:
            Rp = None
                
        ########################################### Q PART
        NC = []
        if Rp != None:
            arrow = ""
        elif Rp == None:
            arrow = "->"
        for j in range(0, len(M1Raw)):                # loop through each Ix & Iy data set   D chr
            if str(M1Raw[j][1]) == str(DC) and str(M1Raw[j][7]) == 'q': #find any q side recombinations
                DC = ConvertInt(DC)
                NC.extend([str(M1Raw[j][0]),str(abs(M1Raw[j][4] - centromeres[DC+1][1])),str(M1Raw[j][0]),str(abs(M1Raw[j][5] - centromeres[DC+1][1]))]) #find distance from centromere for query start and end
                DC = ConvertStr(DC)
        if NC == []:            #if no q side found than set Qname to Empty
            Qname = 'Empty'
        elif NC != []:
            NCmini = NC.index(min(NC))        #find q recombination closest to centromere
            for i in range(len(M1Raw)):        #loop through M1Raw
                if re.search('y',NC[NCmini-1]) != None and M1Raw[i][0] == NC[NCmini-1] and M1Raw[i][6] == -1:    # find right (y) BLAST with NC value, with the same name (e.g., IOy) and (-) direction - this is the normal part of the chromosome
                    DnR0qC = M1Raw[i-1][1]
                    NR0qC = M1Raw[i][1]
                    if M1Final[i][18] == '':
                        BPU = '::'
                    else:
                        BPU =''.join( ['::',M1Final[i][18],"::"])                  
                    DnR0qN = ''.join([arrow,str(NR0qC),str(M1Raw[i][13]),"(",str(M1Final[i][2]),")",BPU,str(DnR0qC),str(M1Raw[i-1][13]),"(",str(M1Final[i-1][2]),")->"])         # add x qside BP information
                    DnR0qD = M1Raw[i-1][6]
                    DnR0qA = M1Raw[i-1][7]
                    DnR0qB = M1Raw[i-1][2]
                    DnR0qBO = M1Raw[i-1][3]
                    DnR0qFASTA = M1Final[i][19]  
                elif re.search('y',NC[NCmini-1]) != None and M1Raw[i][0] == NC[NCmini-1] and M1Raw[i][6] == 1:
                    nomenclature_error.append("Nomenclature Error NE003 - Right side of FASTA Input:\" %0s\" - Sequence is inverted and cannot be closest to the centromere. Please check for other breakpoints on the %1s side of chromosome %2s." % (M1Final[i][17], "qter", DC))
                elif re.search('x',NC[NCmini-1]) != None and M1Raw[i][0] == NC[NCmini-1] and M1Raw[i][6] == 1:  # find left (x) BLAST with NC value, with the same name (e.g., IOx) and (+) direction - this is the normal part of the chromosome
                    DnR0qC = M1Raw[i+1][1]
                    NR0qC = M1Raw[i][1]
                    if M1Final[i][18] == '':
                        BPU = '::'
                    else:
                        BPU =''.join( ['::',M1Final[i][18],"::"])
                    DnR0qN = ''.join([arrow,str(NR0qC),str(M1Raw[i][13]),"(",str(M1Final[i][2]),")",BPU,str(DnR0qC),str(M1Raw[i+1][13]),"(",str(M1Final[i+1][2]),")->"])    # add y qside BP information
                    DnR0qD = M1Raw[i+1][6]
                    DnR0qA = M1Raw[i+1][7]
                    DnR0qB = M1Raw[i+1][2]
                    DnR0qBO = M1Raw[i+1][3]
                    DnR0qFASTA = M1Final[i][19]
                elif re.search('x',NC[NCmini-1]) != None and M1Raw[i][0] == NC[NCmini-1] and M1Raw[i][6] == -1:
                    nomenclature_error.append("Nomenclature Error NE004 - Left side of FASTA Input:\" %0s\" - Sequence is inverted and cannot be closest to the centromere. Please check for other breakpoints on the %1s side of chromosome %2s." % (M1Final[i][17], "pter", DC))
        if DnR0qN == []:
            Rq = None
        elif DnR0qN != []:
            if Rp != None:
                D.append(DC)
            Rq = [[DnR0qB, DnR0qC, DnR0qD, DnR0qA, DnR0qN,'qside', DnR0qFASTA,DnR0qBO]]
            UsageCheck.append(DnR0qFASTA)
            Qname = None
            while Qname == None:
                whilecounter = whilecounter + 1
                if whilecounter > (len(M1Raw)/2):
                    loop_error.append("Loop Error LE002 - Breakpoints create loop, check more multiple homologous chromosome involvement and/or additional breakpoints.")
                    break
                elif whilecounter <= (len(M1Raw)/2):
                    Rq, Qname, DnRnqFASTA  = Q_Extender(M1Raw,M1Final, Rq)
                    if DnRnqFASTA != []:
                        UsageCheck.append(DnRnqFASTA) 
            D.append(Rq)
            Rp=DnR0pB=DnR0pC=DnR0pD=DnR0pA=DnR0pN=DnR0pFASTA=DnR0pBO=[]
            Rq=DnR0qB=DnR0qC=DnR0qD=DnR0qA=DnR0qN=DnR0qFASTA=DnR0qBO=[]
    UsageCounter = []
    for i in range(len(UsageCheck)):
        counter = 0
        UsageCounter.append([0])
        elem = UsageCheck[i]
        for j in range(i,len(UsageCheck)):
            if UsageCheck[j].count(elem) > 0:
                counter += 1
                UsageCounter[i] = counter
    for i in range(len(UsageCounter)):
        if UsageCounter[i] > 1:
            nomenclature_error.append("Nomenclature Error NE005 - FASTA Input:\" %0s\" - FASTA is used multiple times. Check for inversions on either side of breakpoints and/or copy number variations." % UsageCheck[i])
    return D, nomenclature_error, loop_error

def P_Extender(M1Raw, M1Final, Rp, Pname = None):
    temp=DnRnpB=DnRnpC=DnRnpD=DnRnpA=DnRnpN=DnRnpFASTA=DnRnpBO=[]
    for i in range(len(M1Raw)):
        if Rp[len(Rp)-1][1] == M1Raw[i][1]:          # len(Rp)-1 because index starts at 0
            if Rp[len(Rp)-1][0] < Rp[len(Rp)-1][7]:   
                if Rp[len(Rp)-1][7] <  M1Raw[i][3]:
                    temp.extend([M1Raw[i][0],abs(M1Raw[i][4] -  Rp[len(Rp)-1][0])])
                    temp.extend([M1Raw[i][0],abs(M1Raw[i][5] -  Rp[len(Rp)-1][0])])
                elif Rp[len(Rp)-1][7] >=  M1Raw[i][3]:
                    pass
            elif Rp[len(Rp)-1][0] > Rp[len(Rp)-1][7]:   
                if Rp[len(Rp)-1][7] >  M1Raw[i][3]:
                    temp.extend([M1Raw[i][0],abs(M1Raw[i][4] -  Rp[len(Rp)-1][0])])
                    temp.extend([M1Raw[i][0],abs(M1Raw[i][5] -  Rp[len(Rp)-1][0])])
                elif M1Raw[i][4] <= Rp[len(Rp)-1][7] and M1Raw[i][5] <= Rp[len(Rp)-1][7]:
                    pass
    if temp == []:
        Pname = "join(DaRn-1pC,DaRn-1pA,'ter')"
        Rp[len(Rp) - 1][4] = ''.join([str(Rp[len(Rp) - 1][1]),str(Rp[len(Rp) - 1][3]),'ter',str(Rp[len(Rp) - 1][4])])
        return Rp, Pname, DnRnpFASTA
    else:
        tempmini = temp.index(min(temp))
        for i in range(len(M1Raw)):
            if re.search('y',temp[tempmini-1]) != None and M1Raw[i][0] == temp[tempmini-1]:                                 # case: SyC, y == +
                DnRnqC = M1Raw[i-1][1]
                NRnqC = M1Raw[i][1]
                if M1Final[i][18] == '':
                    BPU = '::'
                else:
                    BPU =''.join( ['::',M1Final[i][18],"::"])                  
                DnRnpN = ''.join([str(NRnqC),str(M1Raw[i][13]),"(",str(M1Final[i][2]),")",BPU,str(DnRnqC),str(M1Raw[i-1][13]),"(",str(M1Final[i-1][2]),")->"])         # add x qside BP information
                DnRnpD = M1Raw[i-1][6]
                DnRnpA = M1Raw[i-1][7]
                DnRnpB = M1Raw[i-1][2]
                DnRnpBO = M1Raw[i-1][3]
                DnRnpFASTA = M1Final[i][19]  
            elif re.search('x',temp[tempmini-1]) != None and M1Raw[i][0] == temp[tempmini-1]:                                # case: SxC, x == -   
                DnRnqC = M1Raw[i+1][1]
                NRnqC = M1Raw[i][1]
                if M1Final[i][18] == '':
                    BPU = '::'
                else:
                    BPU =''.join( ['::',M1Final[i][18],"::"])                  
                DnRnpN = ''.join([str(NRnqC),str(M1Raw[i][13]),"(",str(M1Final[i][2]),")",BPU,str(DnRnqC),str(M1Raw[i+1][13]),"(",str(M1Final[i+1][2]),")->"])         # add x qside BP information
                DnRnpD = M1Raw[i+1][6]
                DnRnpA = M1Raw[i+1][7]
                DnRnpB = M1Raw[i+1][2]
                DnRnpBO = M1Raw[i+1][3]
                DnRnpFASTA = M1Final[i][19]  
        Rp.append([DnRnpB, DnRnpC, DnRnpD, DnRnpA, DnRnpN,'pside',DnRnpFASTA,DnRnpBO])
        return Rp, Pname, DnRnpFASTA

def Q_Extender(M1Raw, M1Final, Rq, Qname = None):
    temp = DnRnqC = NRnqC = DnRnqN = DnRnqD = DnRnqA = DnRnqB = DnRnqFASTA = DnRnqBO = []
    for i in range(len(M1Raw)):
        if Rq[len(Rq)-1][1] == M1Raw[i][1]:          # len(Rp)-1 because index starts at 0
            if Rq[len(Rq)-1][0] < Rq[len(Rq)-1][7]:   
                if Rq[len(Rq)-1][7] <  M1Raw[i][3]:
                    temp.extend([M1Raw[i][0],abs(M1Raw[i][4] -  Rq[len(Rq)-1][0])])
                    temp.extend([M1Raw[i][0],abs(M1Raw[i][5] -  Rq[len(Rq)-1][0])])
                elif Rq[len(Rq)-1][7] >=  M1Raw[i][3]:
                    pass
            elif Rq[len(Rq)-1][0] > Rq[len(Rq)-1][7]: 
                if Rq[len(Rq)-1][7] >  M1Raw[i][3]:
                    temp.extend([M1Raw[i][0],abs(M1Raw[i][4] -  Rq[len(Rq)-1][0])])
                    temp.extend([M1Raw[i][0],abs(M1Raw[i][5] -  Rq[len(Rq)-1][0])])
                elif M1Raw[i][4] <= Rq[len(Rq)-1][7] and M1Raw[i][5] <= Rq[len(Rq)-1][7]:
                    pass
    if temp == []:
        Qname = "join(DaRn-1pC,DaRn-1pA,'ter')"
        Rq[len(Rq) - 1][4] = ''.join([str(Rq[len(Rq) - 1][4]),str(Rq[len(Rq) - 1][1]),str(Rq[len(Rq) - 1][3]),'ter'])
        return Rq, Qname, DnRnqFASTA
    else:
        tempmini = temp.index(min(temp))
        for i in range(len(M1Raw)):
            if re.search('y',temp[tempmini-1]) != None and M1Raw[i][0] == temp[tempmini-1]:                                     # case: SyC, y == +
                DnRnqC = M1Raw[i-1][1]
                NRnqC = M1Raw[i][1]
                if M1Final[i][18] == '':
                    BPU = '::'
                else:
                    BPU =''.join( ['::',M1Final[i][18],"::"])                  
                DnRnqN = ''.join([str(NRnqC),str(M1Raw[i][13]),"(",str(M1Final[i][2]),")",BPU,str(DnRnqC),str(M1Raw[i-1][13]),"(",str(M1Final[i-1][2]),")->"])         # add x qside BP information
                DnRnqD = M1Raw[i-1][6]
                DnRnqA = M1Raw[i-1][7]
                DnRnqB = M1Raw[i-1][2]
                DnRnqBO = M1Raw[i-1][3]
                DnRnqFASTA = M1Final[i][19]  
            elif re.search('x',temp[tempmini-1]) != None and M1Raw[i][0] == temp[tempmini-1]:                                    # case: SxC, x == -   
                DnRnqC = M1Raw[i+1][1]
                NRnqC = M1Raw[i][1]
                if M1Final[i][18] == '':
                    BPU = '::'
                else:
                    BPU =''.join( ['::',M1Final[i][18],"::"])                  
                DnRnqN = ''.join([str(NRnqC),str(M1Raw[i][13]),"(",str(M1Final[i][2]),")",BPU,str(DnRnqC),str(M1Raw[i+1][13]),"(",str(M1Final[i+1][2]),")->"])         # add x qside BP information
                DnRnqD = M1Raw[i+1][6]
                DnRnqA = M1Raw[i+1][7]
                DnRnqB = M1Raw[i+1][2]
                DnRnqBO = M1Raw[i+1][3]
                DnRnqFASTA = M1Final[i][19]  
        Rq.append([DnRnqB, DnRnqC, DnRnqD, DnRnqA, DnRnqN,'qside', DnRnqFASTA,DnRnqBO])
        return Rq, Qname, DnRnqFASTA

def Naming(D, DnC):
    Kseq = []
    Kseq.append([''.join(['seq[',assembly,'] '])])
    for DC in DnC:                  # pick corresponding D chr
        finalp = []
        finalq = [] 
        for i in range(len(D)):     # pick corresponding D chr
            if D[i] == DC:             # pick corresponding D chr
                for o in range(len(D[i + 1])):
                    if D[i + 1][o][5] == 'pside':
                        finalp = finalp + [D[i + 1][o][4]]
                        finalp.reverse()              # MAKE SURE REVERSING WITHIN LOOP CORRECT
                    elif D[i + 1][o][5] == 'qside':
                        finalq = finalq + [D[i + 1][o][4]]
        strp = ''.join(finalp)           
        strq = ''.join(finalq)             # empty ones == ''
        if DC == -1:
            DC = "X"        
        if DC == 0:
            DC = "Y"        
        if finalp != [] and finalq != []:
            Kseq.append([''.join(['der(',str(DC),')(',strp,strq,')'])])
        elif finalp == [] and finalq != []:
            Kseq.append([''.join(['der(',str(DC),')(',str(DC),'pter',strp,strq,')'])])
        elif finalp != [] and finalq == []:
            Kseq.append([''.join(['der(',str(DC),')(',strp,strq,str(DC),'qter',')'])])
        elif finalp == [] and finalq == []:
            Kseq.append([''.join(['der(',str(DC),')(',strp,strq,'error unable name centromere',')'])])
    Kseq2 = []
    for i in range(len(Kseq)):
        if i == 0 or i == 1:                    #add "," between der - not working
            Kseq2 = Kseq2 + Kseq[i]
        if i > 1:
            Kseq2 = Kseq2 + [','] + Kseq[i]
    Kseq2 = ''.join(Kseq2)
    return Kseq2
"""End Full Nomenclature"""

"""Finishing"""
def Tables_Small(M1Raw, i):
# <th>Rearrangement </th>
    print "All data within tables are acquired from BLAST"
    #<th>Rearrangement </th>
    print """
        <table border="1">
        <tr>
        <th>Query Start</th>
        <th>Query End</th>
        <th>Query Length</th>
        <th>Query Strand</th>
        <th>Identity</th>
        <th>Sbjct Chromosome</th>
        <th>Sbjct Band</th>
        <th>Sbjct Strand</th>
        <th>Sbjct Start</th>
        <th>Sbjct End</th>
        <th>Max Score</th>
        <th>e-value</th>
        </tr>
    """
    ## ALL THE x stuff here
    print "<tr>"
    #print "<td>%s</td>" % str(i)  # input/seq #
    print "<td>%s</td>" % M1Raw[i*2][9] # query start
    print "<td>%s</td>" % M1Raw[i*2][10] # query end
    print "<td>%s</td>" % str(abs(M1Raw[i*2+1][9] - M1Raw[i*2+1][10]) + 1) # query length
    if M1Raw[i*2][8] == 1:
        print "<td> (+) </td>"  # query (+) strand
    elif M1Raw[i*2][8] == -1:
        print "<td> (-) </td>"  # query (-) strand
    print "<td>%s</td>" % ''.join([str(M1Raw[i*2][16]), '%']) # identity
    print "<td>%s</td>" % M1Raw[i*2][1] # sbjct chromosome
    print "<td>%s</td>" % M1Raw[i*2][13] # sbjct band
    if M1Raw[i*2][6] == 1:
        print "<td> (+) </td>"  # sbjct direction
    elif M1Raw[i*2][6] == -1:
        print "<td> (-) </td>"  # sbjct direction
    print "<td>%s</td>" % M1Raw[i*2][4] # sbjct start
    print "<td>%s</td>" % M1Raw[i*2][5] # sbjct end
    print "<td>%s</td>" % M1Raw[i*2][15] # blast Score
    print "<td>%s</td>" % M1Raw[i*2][14] # blast e-value 
    print "</tr>"

    ## ALL THE y stuff here
    print "<tr>"
    #print "<td>%s</td>" % str(i)  # input/seq #
    print "<td>%s</td>" % M1Raw[i*2+1][9] # query start
    print "<td>%s</td>" % M1Raw[i*2+1][10] # query end
    print "<td>%s</td>" % str(abs(M1Raw[i*2+1][9] - M1Raw[i*2+1][10]) + 1) # query length
    if M1Raw[i*2+1][8] == 1:
        print "<td> (+) </td>"  # query (+) strand
    elif M1Raw[i*2+1][8] == -1:
        print "<td> (-) </td>"  # query (-) strand
    print "<td>%s</td>" % ''.join([str(M1Raw[i*2+1][16]), '%']) # identity
    print "<td>%s</td>" % M1Raw[i*2+1][1] # sbjct chromosome
    print "<td>%s</td>" % M1Raw[i*2+1][13] # sbjct band
    if M1Raw[i*2+1][6] == 1:
        print "<td> (+) </td>"  # sbjct direction
    elif M1Raw[i*2+1][6] == -1:
        print "<td> (-) </td>"  # sbjct direction
    print "<td>%s</td>" % M1Raw[i*2+1][4] # sbjct start
    print "<td>%s</td>" % M1Raw[i*2+1][5] # sbjct end
    print "<td>%s</td>" % M1Raw[i*2+1][15] # blast Score
    print "<td>%s</td>" % M1Raw[i*2+1][14] # blast e-value 
    print "</tr>"
    print "</table>"
    
def HTML_Header():
    print """
    <html>
    <head>
      <title>BOSToN Output</title>
    </head>
    <body>
    <big><span style="font-weight: bold;">BOSToN Output:</span></big><br>
    """ 

def HTML_Footer():
    print """
    </p>
    </body>
    </html>
    """    
def SingleBreakPointTables(M1Raw,i):
    BPTable = []
    QStart = QEnd = QLength = QStrand = Identity = SChr = SBand = SStrand = SStart = SEnd = Score = Evalue = []
    
    ## All the x stuff here
    QStart = str(M1Raw[i*2][9]) # query start
    QEnd = str(M1Raw[i*2][10]) # query end
    QLength = str(abs(M1Raw[i*2][9] - M1Raw[i*2][10]) + 1) # query length
    if M1Raw[i*2][8] == 1:
        QStrand = "(+)"  # query (+) strand
    elif M1Raw[i*2][8] == -1:
        QStrand = "(-)"  # query (-) strand
    Identity =  str(M1Raw[i*2][16]) # identity
    SChr = str(M1Raw[i*2][1]) # sbjct chromosome
    SBand = str(M1Raw[i*2][13]) # sbjct band
    if M1Raw[i*2][6] == 1:
        SStrand = "(+)"  # sbjct direction
    elif M1Raw[i*2][6] == -1:
        SStrand = "(-)"  # sbjct direction
    SStart = str(M1Raw[i*2][4]) # sbjct start
    SEnd = str(M1Raw[i*2][5]) # sbjct end
    Score = str(M1Raw[i*2][15]) # blast Score
    Evalue = str(M1Raw[i*2][14]) # blast e-value
    BPTable.append([QStart, QEnd, QLength, QStrand, Identity, SChr, SBand, SStrand, SStart, SEnd, Score, Evalue]) 
  
    ## ALL THE y stuff here
    QStart = str(M1Raw[i*2+1][9]) # query start
    QEnd = str(M1Raw[i*2+1][10]) # query end
    QLength = str(abs(M1Raw[i*2+1][9] - M1Raw[i*2+1][10]) + 1) # query length
    if M1Raw[i*2+1][8] == 1:
        QStrand = "(+)"  # query (+) strand
    elif M1Raw[i*2+1][8] == -1:
        QStrand = "(-)"  # query (-) strand
    Identity =  str(M1Raw[i*2+1][16]) # identity
    SChr = str(M1Raw[i*2+1][1]) # sbjct chromosome
    SBand = str(M1Raw[i*2+1][13]) # sbjct band
    if M1Raw[i*2+1][6] == 1:
        SStrand = "(+)"  # sbjct direction
    elif M1Raw[i*2+1][6] == -1:
        SStrand = "(-)"  # sbjct direction
    SStart = str(M1Raw[i*2+1][4]) # sbjct start
    SEnd = str(M1Raw[i*2+1][5]) # sbjct end
    Score = str(M1Raw[i*2+1][15]) # blast Score
    Evalue = str(M1Raw[i*2+1][14]) # blast e-value
    BPTable.append([QStart, QEnd, QLength, QStrand, Identity, SChr, SBand, SStrand, SStart, SEnd, Score, Evalue])
    return BPTable
"""End Finishing"""

"""Essential Executes for M1 Matrix"""
def BPExecute(fasta_file):
    #seq = GetXmlBlastRecords(FASTA)
    blast_records = ParseBlastXml(fasta_file)
    blast_output, QxyD, quality_data, FASTAEnds, blast_output_error = BlastOutput(blast_records)
    blast_arm = GetArm(blast_output)
    BANDxy = GetBands(blast_output)
    M1Raw = BPEssentials(blast_output, blast_arm, QxyD, BANDxy, quality_data, FASTAEnds)
    return M1Raw #, seq    
"""End Essential Executes for M1 Matrix"""


"""Code Exectution"""
def WebExecute():
    cgitb.enable() # enable debugging
    start = time.clock()
    print "Content-Type: text/html;charset=utf-8"
    print # needed to separate heading
    HTML_Header()
    InputError = blast_output_error = M1Final = RP = nomenclature_error = nomenclature = [] 
    fasta_input_file = read_fasta_input_from_web(cgi.FieldStorage())
    InputSeqs, InputError = GetInputFASTA(fasta_input_file)
    if InputError != []:
        print '<br><b>Input Errors detected: </b></br>'
        print '</br>.'.join(InputError)
    else:
        blast_records = ParseBlastXml(fasta_input_file)
        blast_output, QxyD, quality_data, FASTAEnds, blast_output_error = BlastOutput(blast_records)
        if blast_output_error != []:
            print '<br><b>Ouput Errors detected: </b></br>'
            print '</br>.'.join(blast_output_error)
        else:
            blast_arm = GetArm(blast_output)
            BANDxy = GetBands(blast_output)
            M1Raw = BPEssentials(blast_output, blast_arm, QxyD, BANDxy, quality_data, FASTAEnds)
            print """<p style="font-family: Courier New,Courier,monospace;">"""
            with open(fasta_input_file, 'r') as fasta_input:
                f_input = fasta_input.read()
                f_input = f_input.replace("\r\n", "<br />")
                print f_input
            print "</p>"
            print "<p>"
            if InputError != []:
                print '<br>Input Errors detected: </br>'
                print '</br>.'.join(InputError)
            elif blast_output_error != []:
                print '</br>'.join(blast_output_error)
            
            elif range(len(M1Raw)) == []:
                print "<br> No Breakpoint Detected </br>"
            
            elif range(len(M1Raw)) == [0,1]:
                print "<br> Single Breakpoint Detected </br>"
                DnO, M1Final = O_Picker(M1Raw,InputSeqs)
                RP = Rearrangement_Pieces(M1Raw, DnO)
                for i in range(len(M1Raw)/2):
                    print """<br><span style="font-weight: bold;">Rearrangement_%s </span></br>""" % NumAlphaConvert(i)
                    Tables_Small(M1Raw, i)
                    print "BLA(S)T Output:  %s " % RP[i]
                    print "<br></br>"
            
            else:
                print "<br>Mulitple Breakpoints Detected</br>"
                DnO, M1Final = O_Picker(M1Raw,InputSeqs)
                RP = Rearrangement_Pieces(M1Raw, DnO)
                for i in range(len(M1Raw)/2):
                    print """<br><span style="font-weight: bold;">Rearrangement_%s </span></br>""" % NumAlphaConvert(i)
                    Tables_Small(M1Raw, i)
                    print "BLA(S)T Output:  %s " % RP[i]
                    print "<br></br>"
                DnC = D_Picker(M1Raw)
                D, nomenclature_error, loop_error = Derivative_Definition(M1Raw, M1Final, DnC)
                if nomenclature_error != [] and loop_error == []:
                    print '</br>'.join(nomenclature_error)
                elif loop_error != []:
                    print '</br>'.join(nomenclature_error)
                    print '</br>'.join(loop_error)
                    print "<br> Nomenclature up to loop: </br>"
                    nomenclature = Naming(D, DnC)
                    print nomenclature
                else:
                    nomenclature = Naming(D, DnC)
                    print nomenclature
    HTML_Footer()
    Errors = InputError, blast_output_error, nomenclature_error, loop_error 
    return Errors, M1Final, RP, nomenclature, start

def PythonExecuteFull(): 
    InputError = blast_output_error = M1Final = RP = nomenclature_error = nomenclature = [] 
    fasta_input_file = raw_input("Enter FASTA.xml file: ")
    start = time.clock()
    read_fasta_input(fasta_input_file)
    InputSeqs, InputError = GetInputFASTA(fasta_input_file)
    if InputError != []:
        print 'Input Errors detected:'
        print '/n'.join(InputError)
    else:    
        print "Executing BLAST. Please wait."
        blast_records = ParseBlastXml(fasta_input_file)
        blast_output, QxyD, quality_data, FASTAEnds, blast_output_error = BlastOutput(blast_records)
        if blast_output_error != []:
            print 'Ouput Errors detected:'
            print '/n.'.join(blast_output_error)
        else:
            blast_arm = GetArm(blast_output)
            BANDxy = GetBands(blast_output)
            M1Raw = BPEssentials(blast_output, blast_arm, QxyD, BANDxy, quality_data, FASTAEnds)
            if range(len(M1Raw)) == []:
                print "No Breakpoint Detected"
                for i in range(len(InputSeqs)):
                    print InputSeqs[i]
            elif range(len(M1Raw)) == [0,1]:
                print "Single Breakpoint Detected"
                DnO, M1Final = O_Picker(M1Raw,InputSeqs)
                RP = Rearrangement_Pieces(M1Raw, DnO)
                for i in range(len(M1Raw)/2):
                    print ">Rearrangement_", NumAlphaConvert(i)
                    BPTableHeader = ["Q_Start","Q_End","Q_Length","Q_Strand","Identity","S_Chr","S_Band","S_Strand","S_Start","S_End","Score", "E-value"]
                    print ', '.join(BPTableHeader)
                    BPTable = SingleBreakPointTables(M1Raw, i)
                    for j in range(len(BPTable)):
                        print ', '.join(BPTable[j])
                    print RP[i]
            else:
                print "Mulitple Breakpoints Detected"
                DnC = D_Picker(M1Raw)
                DnO, M1Final = O_Picker(M1Raw,InputSeqs)
                RP = Rearrangement_Pieces(M1Raw, DnO)
                for i in range(len(M1Raw)/2):
                    print ">Rearrangement_", NumAlphaConvert(i)
                    BPTableHeader = ["Q_Start","Q_End","Q_Length","Q_Strand","Identity","S_Chr","S_Band","S_Strand","S_Start","S_End","Score", "E-value"]
                    print ', '.join(BPTableHeader)
                    BPTable = SingleBreakPointTables(M1Raw, i)
                    for j in range(len(BPTable)):
                        print ', '.join(BPTable[j])
                    print RP[i]
                D, nomenclature_error, loop_error = Derivative_Definition(M1Raw, M1Final, DnC)
                if nomenclature_error != [] and loop_error == []:
                    print '/n'.join(nomenclature_error)
                elif loop_error != []:
                    print '/n'.join(nomenclature_error)
                    print '/n'.join(loop_error)
                    print "Nomenclature at loop:"
                    nomenclature = Naming(D, DnC)
                    print nomenclature
                else:
                    nomenclature = Naming(D, DnC)
                    print nomenclature
    Errors = InputError, blast_output_error, nomenclature_error, loop_error 
    return Errors, M1Final, RP, nomenclature, start

def PythonExecuteReturn(): 
    InputError = blast_output_error = M1Final = RP = nomenclature_error = nomenclature = [] 
    fasta_input_file = raw_input("Enter FASTA.xml file: ")
    read_fasta_input(fasta_input_file)
    InputSeqs, InputError = GetInputFASTA(fasta_input_file)
    if InputError != []:
        ','.join(InputError)
    else:    
        blast_records = ParseBlastXml(fasta_input_file)
        blast_output, QxyD, quality_data, FASTAEnds, blast_output_error = BlastOutput(blast_records)
        if blast_output_error != []:
            '.'.join(blast_output_error)
        else:
            blast_arm = GetArm(blast_output)
            BANDxy = GetBands(blast_output)
            M1Raw = BPEssentials(blast_output, blast_arm, QxyD, BANDxy, quality_data, FASTAEnds)
            if range(len(M1Raw)) == []:
                pass
            elif range(len(M1Raw)) == [0,1]:
                #Single Breakpoint
                DnO, M1Final = O_Picker(M1Raw,InputSeqs)
                RP = Rearrangement_Pieces(M1Raw, DnO)
            else:
                #Multiple Breakpoints
                DnO, M1Final = O_Picker(M1Raw,InputSeqs)
                RP = Rearrangement_Pieces(M1Raw, DnO)
                DnC = D_Picker(M1Raw)
                D, nomenclature_error, loop_error = Derivative_Definition(M1Raw, M1Final, DnC)
                if nomenclature_error != [] and loop_error == []:
                    print '/n'.join(nomenclature_error)
                elif loop_error != []:
                    print '/n'.join(nomenclature_error)
                    print '/n'.join(loop_error)
                    print "Nomenclature up to loop:"
                    nomenclature = Naming(D, DnC)
                else:
                    nomenclature = Naming(D, DnC)
    Errors = InputError, blast_output_error, nomenclature_error, loop_error 
    return Errors, M1Final, RP, nomenclature
"""Code Exectution"""

Errors, M1Final, RP, nomenclature, start = WebExecute()
#Errors, M1Final, RP, nomenclature, start = PythonExecuteFull()
#Errors, M1Final, RP, nomenclature = PythonExecuteReturn()

elapsed = (time.clock() - start)
print elapsed
print(sys.version)
#raw_input("press <enter> to end")
