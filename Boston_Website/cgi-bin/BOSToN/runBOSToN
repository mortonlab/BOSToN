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

### basic converters for nomenclature
def commify(n): #adds commas to an integer n
    if n is None: return None
    n = str(n)
    if '.' in n:
        dollars, cents = n.split('.')
    else:
        dollars, cents = n, None

    r = []
    for i, c in enumerate(str(dollars)[::-1]):
        if i and (not (i % 3)):
            r.insert(0, ',')
        r.insert(0, c)
    out = ''.join(r)
    if cents:
        out += '.' + cents
    return out
def convrearrg(n): # converts numeral to alpha
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
    else:
	out = 'error: too many rearrangements'
    return out

### end basic converters for nomenclature

### Static Definitions
Rtemp = "GRCh37/hg19"
def DefCentromeres(): # retrieve the centromere position
    centromere_file = os.path.join(GetAppDir(),"centromere2.txt")
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
### End of Static Definitions

### Parameters
def GetParameters():
    endfilter = 5 #ensures alignment is near the end of the FASTA rather than a homologous region in center
    return endfilter
### End Parameters

### Get the directory where the application is installed
def GetAppDir():
    return os.path.dirname(os.path.realpath(__file__))
###

### Read input file
def read_fasta_input():
    arg_len = len(sys.argv)
    if arg_len < 2:
        print 'Usage: python(2.7+)', sys.argv[0], '[path_to_fasta_file]'
        sys.exit()
    path_to_fasta_file = sys.argv[1]
    if not os.path.isfile(path_to_fasta_file):
        print 'Fasta input file "', path_to_fasta_file, '" does not exist.'
        sys.exit()
    return path_to_fasta_file
###

### Read web input
def read_fasta_input_from_web(form):
	# save form input to a temp input file
	fasta_input_file = tempfile.NamedTemporaryFile(delete=True).name
	open(fasta_input_file,'w').write(form['fasta_input'].value)
	return fasta_input_file
###

### BLAST I/O
def GetXmlBlastRecords(fasta_file): #module to handle I/O to NCBI Blast and return XML output
    fasta_input_file = fasta_file
    XmlBlastRecords = tempfile.NamedTemporaryFile(delete=True).name
    db_file = os.path.join(GetAppDir(),"human_genomic.02")
    blastcmd = os.path.join(GetAppDir(), "blastn")
    #print db_file
    blastn_cline = NcbiblastnCommandline(cmd=blastcmd, query=fasta_input_file, db=db_file, outfmt=5, out=XmlBlastRecords, task="megablast")
    #print("Running NCBI BLASTn ...")
    blastn_cline()
    return XmlBlastRecords
### End of BLAST I/O

### Parsing XML using BioPython
def ParseBlastXml(fasta_file): #parses xml data into individual records
    XmlBlastRecords = GetXmlBlastRecords(fasta_file)
    result_handle = open(XmlBlastRecords) #temporay file but will return XmlBlastRecords
    blast_records = NCBIXML.parse(result_handle)
    #blast_record = list(blast_records)
    os.remove(XmlBlastRecords)
    return blast_records

def BlastOutput(blast_records): #retrieve necessary data for analysis
    blast_output= []
    blast_output_error = []
    QxyD = []
    quality_data=[]
    SxC = SxBP = SxS = SxE = SxD = QxD = QxS = QxE = xExpect = xScore = xIdentity = []
    SyC = SyBP = SyS = SyE = SyD = QyD = QyS = QyE = xExpect = yScore = yIdentity = []
    for iteration in blast_records:
	endfilter = GetParameters()
	IQLength = iteration.query_length
	for alignment in iteration.alignments:	# find best left(x) BLAST result per submission
	    QScore = 0				    # QScore is score filter (bits) for BLAST Sequence - reset to zero for each alignment 
	    for hsp in alignment.hsps:
		if "Primary Assembly" in alignment.title:	# filters BLAST results to Primary Assembly
		    if hsp.query_start < endfilter:		    # filters for left(x) BLAST results to only those that have a good reads to start of queary
			QScoreTemp = hsp.bits			# set QScore to current record
			if QScore <	QScoreTemp:			# check to see if current record the biggest QScore
				QScore = QScoreTemp				#Set current record as biggest QScore
				D = hsp.frame				#Build Direction vector D = (Q,S) e.g. (1,-1) == (+,-)
				SxD	= D[1]					#Set subject direction
				QxD	= D[0]					#set query direction
				QxS = hsp.query_start			#set query start
				SxS = hsp.sbjct_start			#set subject start
				SxC = str(re.findall('chromosome [\dXY]+', alignment.title)) #find chromosome #
				SxC = int(re.findall('[\dXY]+',	SxC)[0])		#refine chromosome number
				if D == (1,1):				#define BLAST result end (breakpoint)
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
				QxBP = QxE
				xExpect = hsp.expect		    #adding additional information (not for filter)
				xIdentity = (float(hsp.identities) / float(hsp.align_length))*100
				xScore = hsp.bits
	for alignment in iteration.alignments:	# find best right(y) BLAST result per submission
	    QScore = 0				    # QScore is score filter (bits) for BLAST Sequence - reset to zero for each alignment 
	    for hsp in alignment.hsps:
		if "Primary Assembly" in alignment.title:	# filters BLAST results to Primary Assembly		
		    if hsp.query_end > (IQLength - endfilter):  # filters for right(y) BLAST results to only those that have a good reads to end of query
			QScoreTemp = hsp.bits			# set QScore to current record
			if QScore <	QScoreTemp:			# check to see if current record the biggest QScore
				QScore = QScoreTemp				#Set current record as biggest QScore
				D = hsp.frame				#Build Direction vector D = (Q,S) e.g. (1,-1) == (+,-)
				SyD	= D[1]					#Set subject direction
				QyD	= D[0]					#set query direction
				QyS = hsp.query_start			#set query start
				SyS = hsp.sbjct_start			#set subject start
				SyC = str(re.findall('chromosome [\dXY]+', alignment.title)) #find chromosome #
				SyC = re.findall('[\dXY]+',	SyC)[0]		#refine chromosome number
				if SyC == "X":
				    SyC = -1
				if SyC == "Y":
				    SyC = 0
				SyC = int(SyC)
				if D == (1,1):				#define BLAST result end (breakpoint)
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
				QyBP = QyS
				yExpect = hsp.expect		    #adding additional information (not for filter)
				yIdentity = (float(hsp.identities) / float(hsp.align_length))*100
				yScore = hsp.bits  
	if SxC == [] and SyC == []:	#Error Handling for poor alignments
	    blast_output_error += " BLAST Output Error 0001 - Fasta Input: ",iteration.query," did not align well to reference."    
	elif SxC == [] and SyC != []:
	    blast_output_error += " BLAST Output Error 0002 - Start of Fasta Input: ",iteration.query," did not align well to reference."
	elif SxC != [] and SyC == []:
	    blast_output_error += " BLAST Output Error 0003 - End of Fasta Input: ",iteration.query," did not align well to reference."
	elif SxC != [] and SyC != []: #load data into string for each iteration
	    blast_output.append([SxC, SxBP, SxS, SxE, SxD])
	    blast_output.append([SyC, SyBP, SyS, SyE, SyD])
	    QxyD.append([QxD, QxS, QxE])
	    QxyD.append([QyD, QyS, QyE])
	    quality_data.append([xExpect, xScore, xIdentity])
	    quality_data.append([yExpect, yScore, yIdentity])	
    return blast_output, QxyD, quality_data, blast_output_error
### End of Parsing XML using BioPython

### Defining Additional Data Required for BP Ends
def GetArm(blast_output): #determines chromosome arm of breakpoints (i.e. p or q)
	blast_arm = range(0, len(blast_output))
	centromeres = DefCentromeres()
	for i in range(0, len(blast_output)):
		#x chr
		if blast_output[i][0] == 'X':
			if blast_output[i][2] > centromeres[0][1]:
				SxyA = 'q'
			elif blast_output[i][2] < centromeres[0][1]:
				SxyA = 'p'
		#y chr
		elif blast_output[i][0] == 'Y':
			if blast_output[i][2] > centromeres[1][1]:
				SxyA = 'q'
			elif blast_output[i][2] < centromeres[1][1]:
				SxyA = 'p'
		#remaining chromosomes arm definition [SxC + 1] needed for chr1 is in 3rd position == [2]
		else:
			if blast_output[i][2] > centromeres[int(blast_output[i][0]) + 1][1]:
				SxyA = 'q'
			elif blast_output[i][2] < centromeres[int(blast_output[i][0]) + 1][1]:
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
    
### End of Defining Additional Data Required for BP Ends

### Combine Essential Data for BP Analysis
def BPEssentials(blast_output, blast_arm, QxyD, BANDxy, quality_data, M1 = "Nothing"): #creates a list for each breakpoint end
	blast_output2 = map(list, blast_output)
	blast_arm2 = map(list, blast_arm)
	M1 = list(range(0,len(blast_output2)))
	for i in range(0,len(blast_output2)):
	    M1[i] = blast_output2[i] + blast_arm2[i] + QxyD[i] + [BANDxy[i]] + quality_data[i]
		# print blast_output2[i] + blast_arm2[i]
	for i in range(0,len(blast_output2)/2):
		tempx = [''.join(['I',str(i),'x'])]
		tempy = [''.join(['I',str(i),'y'])]
		M1[i*2] = list(tempx) + M1[i*2]
		M1[i*2 + 1] = list(tempy) + M1[i*2 + 1]
		# print i
		# print i*2
		# print i*2 + 1
	return M1




def BPExecute(fasta_file):
    #seq = GetXmlBlastRecords(FASTA)
    blast_records = ParseBlastXml(fasta_file)
    blast_output, QxyD, quality_data, blast_output_error = BlastOutput(blast_records)
    blast_arm = GetArm(blast_output)
    BANDxy = GetBands(blast_output)
    M1 = BPEssentials(blast_output, blast_arm, QxyD, BANDxy, quality_data)
    return M1 #, seq    
### End Combine Essential Data for BP Analysis

### Build Chromosomes
def D_Picker(M1): #pick derivative chromosome
	DnC = list(range(0,len(M1)))
	for i in range(0, len(M1)):
		DnC[i] = M1[i][1]
	DnC = list(set(DnC))
	for i in range(len(DnC)):
		if DnC[i] == 'X':
			DnC[i] = -1
		elif DnC[i] == 'Y':
			DnC[i] = 0
		else:
			DnC[i] = int(DnC[i])
	DnC.sort()
	for i in range(len(M1)):
		if M1[i][1] != 'X' and M1[i][1] != 'Y':  
			M1[i][1] = int(M1[i][1])
	return DnC

def O_Picker(M1, seq):
	DnO = []
	for i in range(len(M1)/2):   #Find overlap
		if M1[i*2][7] == 1 and M1[i*2 + 1][7] == 1:
			DnO.append([M1[i*2][9] - M1[i*2 + 1][8]])
		elif M1[i*2][7] == 1 and M1[i*2 + 1][7] == -1:
			DnO.append([M1[i*2][9] - M1[i*2 + 1][9]])
		elif M1[i*2][7] == -1 and M1[i*2 + 1][7] == 1:
			DnO.append([M1[i*2][8] - M1[i*2 + 1][8]])
		elif M1[i*2][7] == -1 and M1[i*2 + 1][7] == -1:
			DnO.append([M1[i*2][8] - M1[i*2 + 1][9]])
	for i in range(len(M1)/2):    #Find BP overlap or insertion or normal
		if DnO[i][0] == -1:
			DnO[i].extend(['NORMAL'])
		elif DnO[i][0] < -1:
			DnO[i].extend(['INSERTION'])
			if M1[i*2][7] == 1:
				DnO[i].extend([M1[i*2][9]])
			if M1[i*2][7] == -1:
				DnO[i].extend([M1[i*2][8]])
			if M1[i*2 + 1][7] == 1:
				DnO[i].extend([M1[i*2 + 1][8]])
			if M1[i*2 + 1][7] == -1:
				DnO[i].extend([M1[i*2 + 1][9]]) 
			temp = list(seq[i][0])
			m = ''.join(temp[DnO[i][2]:DnO[i][3] - 1])
			DnO[i].append(m)
		elif DnO[i][0] > -1:
			DnO[i].extend(['OVERLAP'])
			if M1[i*2][5] == 1:
				DnO[i].extend([M1[i*2][2] - DnO[i][0], M1[i*2][2]])
			elif M1[i*2][5] == -1:
				DnO[i].extend([M1[i*2][2] + DnO[i][0], M1[i*2][2]])
			if M1[i*2 + 1][5] == 1:
				DnO[i].extend([M1[i*2 + 1][2], M1[i*2 + 1][2] + DnO[i][0]])
			elif M1[i*2 + 1][5] == -1:
				DnO[i].extend([M1[i*2 + 1][2], M1[i*2 + 1][2] - DnO[i][0]])
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
				final[len(final) - 1] = finalBP
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
				final[len(final) - 1] = finalBP
				final = ''.join(final)
				DnO[i].append(final)
	for l in range(len(DnO)):
		if DnO[l][1] == 'OVERLAP':
			M1[l*2][2] = DnO[l][6]
			M1[l*2 + 1][2] = DnO[l][7]
	return DnO

def Q_Extender(M1, Rq, Qname = None):
    if Rq == []:
	Rq = Qname = []
    elif Rq != []:
	temp = []
	for i in range(len(M1)):
		if Rq[len(Rq)-1][1] == M1[i][1]:          # len(Rq)-1 because index starts at 0
			if Rq[len(Rq)-1][2] == 1:   
				if M1[i][3] > Rq[len(Rq)-1][0] and M1[i][4] > Rq[len(Rq)-1][0]:
					temp.extend([M1[i][0],abs(M1[i][3] -  Rq[len(Rq)-1][0])])
					temp.extend([M1[i][0],abs(M1[i][4] -  Rq[len(Rq)-1][0])])
				elif M1[i][3] <= Rq[len(Rq)-1][0] and M1[i][4] <= Rq[len(Rq)-1][0]:
					pass
			elif Rq[len(Rq)-1][2] == -1:   
				if M1[i][3] < Rq[len(Rq)-1][0] and M1[i][4] < Rq[len(Rq)-1][0]:
					temp.extend([M1[i][0],abs(M1[i][3] -  Rq[len(Rq)-1][0])])
					temp.extend([M1[i][0],abs(M1[i][4] -  Rq[len(Rq)-1][0])])
				elif M1[i][3] >= Rq[len(Rq)-1][0] and M1[i][4] >= Rq[len(Rq)-1][0]:
					pass
	if temp == []:
		Qname = "join(DaRn-1pC,DaRn-1pA,'ter')"
		Rq[len(Rq) - 1][4] = ''.join([str(Rq[len(Rq) - 1][4]),str(Rq[len(Rq) - 1][1]),str(Rq[len(Rq) - 1][3]),'ter'])
		return Rq, Qname
	else:
		tempmini = temp.index(min(temp))
		for i in range(len(M1)):
			if re.search('y',temp[tempmini-1]) != None and M1[i][0] == temp[tempmini-1]:                                     # case: SyC, y == +
			    if DnRnqC == -1:
				DnRnqC = "X"	    
			    if DnRnqC == 0:
				DnRnqC = "Y"
			    NRnqC = M1[i][1]
			    if NRnqC == -1:
				NRnqC = "X"	    
			    if NRnqC == 0:
				NRnqC = "Y"				
			    DnRnqN = ''.join([str(DnRnqC),str(M1[i][10]),'(',str(M1[i][2]),')::',str(NRnqC),str(M1[i-1][10]),'(',str(M1[i-1][2]),')->'])                                 
			    DnRnqD = M1[i-1][5]
			    DnRnqA = M1[i-1][6]
			    if M1[i-1][5] == 1:
				    DnRnqB = M1[i-1][4]                  # 3 start  
			    elif M1[i-1][5] == -1:
				    DnRnqB = M1[i-1][3]                  # 4 end  
			elif re.search('x',temp[tempmini-1]) != None and M1[i][0] == temp[tempmini-1]:                                    # case: SxC, x == -   
			    if DnRnqC == -1:
				DnRnqC = "X"	    
			    if DnRnqC == 0:
				DnRnqC = "Y"
			    NRnqC = M1[i][1]
			    if NRnqC == -1:
				NRnqC = "X"	    
			    if NRnqC == 0:
				NRnqC = "Y"					
			    DnRnqN = ''.join([str(DnRnqC),str(M1[i][10]),'(',str(M1[i][2]),')::',str(NRnqC),str(M1[i+1][10]),'(',str(M1[i+1][2]),')->'])
			    DnRnqC = M1[i+1][1]
			    DnRnqD = M1[i+1][5]
			    DnRnqA = M1[i+1][6]
			    if M1[i+1][5] == 1:
				    DnRnqB = M1[i+1][4]                  # 3 start  
			    elif M1[i+1][5] == -1:
				    DnRnqB = M1[i+1][3]                  # 4 end  
		Rq.append([DnRnqB, DnRnqC, DnRnqD, DnRnqA, DnRnqN,'qside'])
		return Rq, Qname
def P_Extender(M1, Rp, Pname = None):
	temp = []
	for i in range(len(M1)):
		if Rp[len(Rp)-1][1] == M1[i][1]:          # len(Rp)-1 because index starts at 0
			if Rp[len(Rp)-1][2] == 1:   
				if M1[i][3] < Rp[len(Rp)-1][0] and M1[i][4] < Rp[len(Rp)-1][0]:
					temp.extend([M1[i][0],abs(M1[i][3] -  Rp[len(Rp)-1][0])])
					temp.extend([M1[i][0],abs(M1[i][4] -  Rp[len(Rp)-1][0])])
				elif M1[i][3] >= Rp[len(Rp)-1][0] and M1[i][4] >= Rp[len(Rp)-1][0]:
					pass
			elif Rp[len(Rp)-1][2] == -1:   
				if M1[i][3] > Rp[len(Rp)-1][0] and M1[i][4] > Rp[len(Rp)-1][0]:
					temp.extend([M1[i][0],abs(M1[i][3] -  Rp[len(Rp)-1][0])])
					temp.extend([M1[i][0],abs(M1[i][4] -  Rp[len(Rp)-1][0])])
				elif M1[i][3] <= Rp[len(Rp)-1][0] and M1[i][4] <= Rp[len(Rp)-1][0]:
					pass
	if temp == []:
		Pname = "join(DaRn-1pC,DaRn-1pA,'ter')"
		Rp[len(Rp) - 1][4] = ''.join([str(Rp[len(Rp) - 1][1]),str(Rp[len(Rp) - 1][3]),'ter',str(Rp[len(Rp) - 1][4])])
		return Rp, Pname
	else:
		tempmini = temp.index(min(temp))
		for i in range(len(M1)):
			if re.search('y',temp[tempmini-1]) != None and M1[i][0] == temp[tempmini-1]:                                 # case: SyC, y == +
			    DnRnpC = M1[i-1][1]
			    if DnRnpC == -1:
				DnRnpC = "X"	    
			    if DnRnpC == 0:
				DnRnpC = "Y"
			    NRnpC = M1[i][1]
			    if NRnpC == -1:
				NRnpC = "X"	    
			    if NRnpC == 0:
				NRnpC = "Y"				
			    DnRnpN = ''.join([str(DnRnpC),str(M1[i-1][10]),'(',str(M1[i-1][2]),')::',str(NRnpC),str(M1[i][10]),'(',str(M1[i][2]),')->'])          
			    DnRnpD = M1[i-1][5]
			    DnRnpA = M1[i-1][6]
			    if M1[i-1][5] == 1:
				    DnRnpB = M1[i-1][3]                  # 3 start  
			    elif M1[i-1][5] == -1:
				    DnRnpB = M1[i-1][4]                  # 4 end  
			elif re.search('x',temp[tempmini-1]) != None and M1[i][0] == temp[tempmini-1]:                                # case: SxC, x == -   
			    DnRnpC = M1[i-1][1]
			    if DnRnpC == -1:
				DnRnpC = "X"	    
			    if DnRnpC == 0:
				DnRnpC = "Y"
			    NRnpC = M1[i][1]
			    if NRnpC == -1:
				NRnpC = "X"	    
			    if NRnpC == 0:
				NRnpC = "Y"				
			    DnRnpN = ''.join([str(DnRnpC),str(M1[i+1][10]),'(',str(M1[i+1][2]),')::',str(NRnpC),str(M1[i][10]),'(',str(M1[i][2]),')->'])
			    DnRnpD = M1[i+1][5]
			    DnRnpA = M1[i+1][6]
			    if M1[i+1][5] == 1:
				    DnRnpB = M1[i+1][3]                  # 3 start  
			    elif M1[i+1][5] == -1:
				    DnRnpB = M1[i+1][4]                  # 4 end  
		Rp.append([DnRnpB, DnRnpC, DnRnpD, DnRnpA, DnRnpN,'pside'])
		return Rp, Pname


def Derivative_Definition(M1, DnC):
    D = []
    centromeres = DefCentromeres()
    Rp=DnR0pB=DnR0pC=DnR0pD=DnR0pA=DnR0pN=[]
    Rq=DnR0qB=DnR0qC=DnR0qD=DnR0qA=DnR0qN=[]
    for DC in DnC:               # loop through each D chr
	NC = []      
	D.append(DC)
	########################################### P PART
	for i in range(0, len(M1)):                # loop through each Ix & Iy data set   D chr
	    if M1[i][1] == DC and M1[i][6] == 'p':          #find any p side recombinations
		NC.extend([M1[i][0],abs(M1[i][3] - centromeres[DC+1][0]),M1[i][0],abs(M1[i][4] - centromeres[DC+1][0])]) #find distance from centromere for query start and end
	if NC == []:			#if no p side found than set Pname to None
	    Pname = 'Empty'
	elif NC != []:
	    NCmini = NC.index(min(NC))	    #find p recombination closest to centromere
	    for i in range(len(M1)):	    #loop through M1
		if re.search('y',NC[NCmini-1]) != None and M1[i][0] == NC[NCmini-1] and M1[i][5] == 1:	# find right (y) BLAST with NC value, with the same name (e.g., IOy) and (+) direction - this is the normal part of the chromosome
		    DnR0pC = M1[i-1][1]
		    if DnR0pC == -1:
			DnR0pC = "X"	    
		    if DnR0pC == 0:
			DnR0pC = "Y"
		    NR0pC = M1[i][1]
		    if NR0pC == -1:
			NR0pC = "X"	    
		    if NR0pC == 0:
			NR0pC = "Y"		    
		    DnR0pN = ''.join(["->",str(DnR0pC),str(M1[i-1][10]),"(",str(M1[i-1][2]),")::",str(NR0pC),str(M1[i][10]),"(",str(M1[i][2]),")->"])         # add x pside BP information
		    DnR0pD = M1[i-1][5]
		    DnR0pA = M1[i-1][6]
		    if M1[i-1][5] == 1:
			DnR0pB = M1[i-1][3]                  # 3 start  
		    elif M1[i-1][5] == -1:
			DnR0pB = M1[i-1][4]                  # 4 end  
		elif re.search('x',NC[NCmini-1]) != None and M1[i][0] == NC[NCmini-1] and M1[i][5] == -1:  # find left (x) BLAST with NC value, with the same name (e.g., IOx) and (-) direction - this is the normal part of the chromosome
		    DnR0pC = M1[i+1][1]
		    if DnR0pC == -1:
			DnR0pC = "X"	    
		    if DnR0pC == 0:
			DnR0pC = "Y"
		    NR0pC = M1[i][1]
		    if NR0pC == -1:
			NR0pC = "X"	    
		    if NR0pC == 0:
			NR0pC = "Y"		    
		    DnR0pN = ''.join(["->",str(DnR0pC),str(M1[i+1][10]),"(",str(M1[i+1][2]),")::",str(NR0pC),str(M1[i][10]),"(",str(M1[i][2]),")->"])	# add y pside BP information
		    DnR0pD = M1[i+1][5]
		    DnR0pA = M1[i+1][6]
		    if M1[i+1][5] == 1:
			DnR0pB = M1[i+1][3]                  # 3 start  
		    elif M1[i+1][5] == -1:
			DnR0pB = M1[i+1][4]                  # 4 end  		
	if DnR0pN != []:
	    Rp = [[DnR0pB, DnR0pC, DnR0pD, DnR0pA, DnR0pN,'pside']]
	    D.append(Rp)
	    Pname = None
	    while Pname == None:
		Rp, Pname = P_Extender(M1, Rp)	
	    Rp=DnR0pB=DnR0pC=DnR0pD=DnR0pA=DnR0pN=[]
	elif DnR0pN == []:
	    Rp = None
    
	########################################### Q PART
	for i in range(0, len(M1)):                # loop through each Ix & Iy data set   D chr
	    if M1[i][1] == DC and M1[i][6] == 'q':          #find any q side recombinations
		NC.extend([M1[i][0],abs(M1[i][3] - centromeres[DC+1][0]),M1[i][0],abs(M1[i][4] - centromeres[DC+1][0])]) #find distance from centromere for query start and end
	if NC == []:			#if no q side found than set Qname to Empty
	    Qname = 'Empty'
	elif NC != []:
	    NCmini = NC.index(min(NC))	    #find q recombination closest to centromere
	    for i in range(len(M1)):	    #loop through M1
		if re.search('y',NC[NCmini-1]) != None and M1[i][0] == NC[NCmini-1] and M1[i][5] == -1:	# find right (y) BLAST with NC value, with the same name (e.g., IOy) and (-) direction - this is the normal part of the chromosome
		    DnR0qC = M1[i-1][1]
		    if DnR0qC == -1:
			DnR0qC = "X"	    
		    if DnR0qC == 0:
			DnR0qC = "Y"
		    NR0qC = M1[i][1]
		    if NR0qC == -1:
			NR0qC = "X"	    
		    if NR0qC == 0:
			NR0qC = "Y"		    
		    DnR0qN = ''.join(["->",str(NR0qC),str(M1[i][10]),"(",str(M1[i][2]),")::",str(DnR0qC),str(M1[i-1][10]),"(",str(M1[i-1][2]),")->"])         # add x qside BP information
		    DnR0qD = M1[i-1][5]
		    DnR0qA = M1[i-1][6]
		    if M1[i-1][5] == 1:
			DnR0qB = M1[i-1][3]                  # 3 start  
		    elif M1[i-1][5] == -1:
			DnR0qB = M1[i-1][4]                  # 4 end  
		elif re.search('x',NC[NCmini-1]) != None and M1[i][0] == NC[NCmini-1] and M1[i][5] == 1:  # find left (x) BLAST with NC value, with the same name (e.g., IOx) and (+) direction - this is the normal part of the chromosome
		    DnR0qC = M1[i+1][1]
		    if DnR0qC == -1:
			DnR0qC = "X"	    
		    if DnR0qC == 0:
			DnR0qC = "Y"
		    NR0qC = M1[i][1]
		    if NR0qC == -1:
			NR0qC = "X"	    
		    if NR0qC == 0:
			NR0qC = "Y"			
		    DnR0qN = ''.join(["->",str(NR0qC),str(M1[i][10]),"(",str(M1[i][2]),")::",str(DnR0qC),str(M1[i+1][10]),"(",str(M1[i+1][2]),")->"])	# add y qside BP information
		    DnR0qD = M1[i+1][5]
		    DnR0qA = M1[i+1][6]
		    if M1[i+1][5] == 1:
			DnR0qB = M1[i+1][3]                  # 3 start  
		    elif M1[i+1][5] == -1:
			DnR0qB = M1[i+1][4]                  # 4 end  
	if DnR0qN != []:
	    Rq = [[DnR0qB, DnR0qC, DnR0qD, DnR0qA, DnR0qN,'qside']]
	    D.append(Rq)
	    Qname = None
	    while Qname == None:
		Rq, Qname = Q_Extender(M1,Rq)
	    Rq=DnR0qB=DnR0qC=DnR0qD=DnR0qA=DnR0qN=[]
	elif DnR0qN == []:
	    Rq = None
    return D


def Naming(D, DnC):
	Kseq = []
	Kseq.append([''.join(['seq[',Rtemp,'] '])])
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
			Kseq.append([''.join(['der(',str(DC),')(',strp,strq,'),'])])
		elif finalp == [] and finalq != []:
			Kseq.append([''.join(['der(',str(DC),')(',str(DC),'pter',strp,strq,'),'])])
		elif finalp != [] and finalq == []:
			Kseq.append([''.join(['der(',str(DC),')(',strp,strq,str(DC),'qter','),'])])
		elif finalp == [] and finalq == []:
			Kseq.append([''.join(['der(',str(DC),')(',strp,strq,'error unable name centromere',')'])])
	Kseq2 = []
	for i in range(len(Kseq)):
		if not i == 0:					#add "," between der - not working
			Kseq.append(',')
		Kseq2 = Kseq2 + Kseq[i]
	Kseq2 = ''.join(Kseq2)
	return Kseq2



### End Build Chromosome

### Display Information
def Tables(M1 = None):
	if M1 == None:
		pass
	elif M1 != None:
		print """
			<table border="1">
			<tr>
			<th>Rearrangement </th>
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
		
		for i in range(len(M1)/2):
			## ALL THE x stuff here
			print "<tr>"
			print "<td>%s</td>" % convrearrg(i) #input/seq #
			print "<td>%s</td>" % M1[i*2][8] # query start
			print "<td>%s</td>" % M1[i*2][9] # query end
			print "<td>%s</td>" % str(abs(M1[i*2][8] - M1[i*2][9]) + 1) # query length
			if M1[i*2][7] == 1:
				print "<td> (+) </td>"  # query (+) strand
			elif M1[i*2][7] == -1:
				print "<td> (-) </td>"  # query (-) strand
			print "<td>%s</td>" % ''.join([str(M1[i*2][13]), '%']) # identity
			print "<td>%s</td>" % M1[i*2][1] # sbjct chromosome
			print "<td>%s</td>" % M1[i*2][10] # sbjct band
			if M1[i*2][5] == 1:
				print "<td> (+) </td>"  # sbjct direction
			elif M1[i*2][5] == -1:
				print "<td> (-) </td>"  # sbjct direction
			print "<td>%s</td>" % M1[i*2][3] # sbjct start
			print "<td>%s</td>" % M1[i*2][4] # sbjct end
			print "<td>%s</td>" % M1[i*2][12] # blast Score
			print "<td>%s</td>" % M1[i*2][11] # blast e-value 
			print "</tr>"
			
			## ALL THE y stuff here
			print "<tr>"
			print "<td>%s</td>" % convrearrg(i)  # input/seq #
			print "<td>%s</td>" % M1[i*2+1][8] # query start
			print "<td>%s</td>" % M1[i*2+1][9] # query end
			print "<td>%s</td>" % str(abs(M1[i*2+1][8] - M1[i*2+1][9]) + 1) # query length
			if M1[i*2+1][7] == 1:
				print "<td> (+) </td>"  # query (+) strand
			elif M1[i*2+1][7] == -1:
				print "<td> (-) </td>"  # query (-) strand
			print "<td>%s</td>" % ''.join([str(M1[i*2+1][13]), '%']) # identity
			print "<td>%s</td>" % M1[i*2+1][1] # sbjct chromosome
			print "<td>%s</td>" % M1[i*2+1][10] # sbjct band
			if M1[i*2+1][5] == 1:
				print "<td> (+) </td>"  # sbjct direction
			elif M1[i*2+1][5] == -1:
				print "<td> (-) </td>"  # sbjct direction
			print "<td>%s</td>" % M1[i*2+1][3] # sbjct start
			print "<td>%s</td>" % M1[i*2+1][4] # sbjct end
			print "<td>%s</td>" % M1[i*2+1][12] # blast Score
			print "<td>%s</td>" % M1[i*2+1][11] # blast e-value 
			print "</tr>"
		print "</table>"

def Tables_Small(M1, i):
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
	print "<td>%s</td>" % M1[i*2][8] # query start
	print "<td>%s</td>" % M1[i*2][9] # query end
	print "<td>%s</td>" % str(abs(M1[i*2+1][8] - M1[i*2+1][9]) + 1) # query length
	if M1[i*2][7] == 1:
		print "<td> (+) </td>"  # query (+) strand
	elif M1[i*2][7] == -1:
		print "<td> (-) </td>"  # query (-) strand
	print "<td>%s</td>" % ''.join([str(M1[i*2][13]), '%']) # identity
	print "<td>%s</td>" % M1[i*2][1] # sbjct chromosome
	print "<td>%s</td>" % M1[i*2][10] # sbjct band
	if M1[i*2][5] == 1:
		print "<td> (+) </td>"  # sbjct direction
	elif M1[i*2][5] == -1:
		print "<td> (-) </td>"  # sbjct direction
	print "<td>%s</td>" % M1[i*2][3] # sbjct start
	print "<td>%s</td>" % M1[i*2][4] # sbjct end
	print "<td>%s</td>" % M1[i*2][12] # blast Score
	print "<td>%s</td>" % M1[i*2][11] # blast e-value 
	print "</tr>"
			
	## ALL THE y stuff here
	print "<tr>"
	#print "<td>%s</td>" % str(i)  # input/seq #
	print "<td>%s</td>" % M1[i*2+1][8] # query start
	print "<td>%s</td>" % M1[i*2+1][9] # query end
	print "<td>%s</td>" % str(abs(M1[i*2+1][8] - M1[i*2+1][9]) + 1) # query length
	if M1[i*2+1][7] == 1:
		print "<td> (+) </td>"  # query (+) strand
	elif M1[i*2+1][7] == -1:
		print "<td> (-) </td>"  # query (-) strand
	print "<td>%s</td>" % ''.join([str(M1[i*2+1][13]), '%']) # identity
	print "<td>%s</td>" % M1[i*2+1][1] # sbjct chromosome
	print "<td>%s</td>" % M1[i*2+1][10] # sbjct band
	if M1[i*2+1][5] == 1:
		print "<td> (+) </td>"  # sbjct direction
	elif M1[i*2+1][5] == -1:
		print "<td> (-) </td>"  # sbjct direction
	print "<td>%s</td>" % M1[i*2+1][3] # sbjct start
	print "<td>%s</td>" % M1[i*2+1][4] # sbjct end
	print "<td>%s</td>" % M1[i*2+1][12] # blast Score
	print "<td>%s</td>" % M1[i*2+1][11] # blast e-value 
	print "</tr>"
	print "</table>"
def Rearrangement_Pieces(M1, DnO):
	RP = []
	for i in range(len(M1)/2):
		temp1 = ''.join([str(M1[i*2][1]), str(M1[i*2][10])])
		if M1[i*2][5] == 1:
			temp2 = '(+)'
		elif M1[i*2][5] == -1:
			temp2 = '(-)'
		temp3 = ''.join(['(', str(M1[i*2][2]), ')'])
		if DnO[i][1] == 'OVERLAP' or DnO[i][1] == 'NORMAL' :
			temp4 = ''.join(['::'])
		elif DnO[i][1] == 'INSERTION' :
			temp4 = ''.join(['::',DnO[i][4],'::'])
		temp5 = ''.join([str(M1[i*2 + 1][1]), str(M1[i*2 + 1][10])])
		if M1[i*2 + 1][5] == 1:
			temp6 = '(+)'
		elif M1[i*2 + 1][5] == -1:
			temp6 = '(-)'
		temp7 = ''.join(['(', str(M1[i*2 + 1][2]), ')'])
		RP.extend([''.join([temp1,temp2,temp3,temp4,temp5,temp6,temp7])])
	return RP
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
### End Display Information

### Code Execution
cgitb.enable() # enable debugging

print "Content-Type: text/html;charset=utf-8"
print # needed to separate heading
fasta_input_file = read_fasta_input_from_web(cgi.FieldStorage())
M1 = BPExecute(fasta_input_file)
#print (M1)
DnC = D_Picker(M1)
#print (DnC)
D = Derivative_Definition(M1,DnC)
#print (D)
nomenclature = Naming(D, DnC)

HTML_Header()
print """<p style="font-family: Courier New,Courier,monospace;">"""
with open(fasta_input_file, 'r') as fasta_input:
    f_input = fasta_input.read()
    f_input = f_input.replace("\r\n", "<br />")
    print f_input
print "</p>"
print "<p>"
for i in range(len(M1)/2):
    print """<br><span style="font-weight: bold;">Rearrangement_%s </span></br>""" % convrearrg(i)
    Tables_Small(M1, i)
print "</p>"
print "<p>"
#print "test"
print (nomenclature)
HTML_Footer()

### End Code Execution
