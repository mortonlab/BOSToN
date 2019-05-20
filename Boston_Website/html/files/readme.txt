-----------------------------------------------------------------
BOSToN - BLA(S)T Ouput Sequence Tool of Nomenclature (Verion 1.0)
-----------------------------------------------------------------

I. Description
--------------
The BLA(S)T Ouput Sequence Tool of Nomenclature (BOSToN) is an algorithm designed to parse FASTA sequences,
covering multiple chromosomal breakpoints, in order to produce a clinically relevant nomenclature of nucleotide
resolved blanced rearrangements.

II. Copyright Notice
--------------------
© 2013 The Brigham and Women’s Hospital.  All Rights Reserved.  

III. Important Notices
--------------------
Information obtained by accessing the BOSToN Software (“Software”) under this Evaluation License should not be 
used for clinical diagnosis or patient care and management. The Software is not intended to render medical advice, 
or to replace or overrule a qualified, licensed health care provider's judgment or clinical diagnosis. Reports, 
reference sequence and variant computations and information generated during use of this Software must be reviewed 
and fully validated by appropriately qualified, trained and licensed clinical personnel. Any action taken or inferences 
made by user in response to the information provided herein is solely at the user's discretion. User assumes full 
responsibility for its use of information obtained from this Software.

IV. Evaluation License
-----------------------
The Brigham and Women’s Hospital (BWH) grants Company a limited, non-exclusive, non-transferable, 30-day license 
from the date of receipt of the Software  to use the Software for evaluation purposes only. Company also 
acknowledges that it is obtaining no title to or ownership of the Software or other materials received as a result 
of this License. Company shall use the Software solely for the purpose of evaluating whether to enter into a licensing 
relationship with BWH. Company shall not copy the Software and shall not transfer or disclose the Software to any 
other person or entity except for those employees of Company who require such knowledge of the Software in order to 
evaluate the Software. Company shall not reverse assemble or reverse compile any Software provided, in whole or in part, 
or permit any other person to do so. Nothing in this License shall be construed as granting or conferring, either express 
or implied, any rights by license or otherwise, under any patent, copyright, or other intellectual property rights 
relating to the Software with the exception of the Evaluation License described above.

V. BOSToN Versions
--------------------
BOSToN has been designed to run in three different modes:
1 - Standalone executable running in a Windows environment (described in this release)
2 - A callable package which returns error calls, breakpoint information, list of essential data and 
    final nomenclature (not supported in this release)
3 - An html package which can be run within a website environment and can be found at boston.bwh.harvard.edu
    (not supported in this release)

VI. Required Supporting Software
---------------------------------
BOSToN runs on Python v2.7.6 and uses both the NumPy v1.8.1 and BioPython supporting packages and must by setup and
installed before running BOSToN. This supporting software can be found at these locations:

For 32 bit operating systems -
Python v2.7.6 - https://www.python.org/downloads
NumPy v1.8.0 (x32) - http://www.lfd.uci.edu/~gohlke/pythonlibs *click on: numpy-MKL-1.8.1.win32-py2.7.exe*
BioPython (x32) - http://www.lfd.uci.edu/~gohlke/pythonlibs *click on: biopython-1.63.win32-py2.7.exe*

For 64 bit operating systems -
Python v2.7.6 - https://www.python.org/downloads *download the Nov. 10, 2013 release (new release may have registry errors with 64)*
NumPy v1.8.1 (x64) - http://www.lfd.uci.edu/~gohlke/pythonlibs *click on: numpy-MKL-1.8.1.win-amd64-py2.7.exe*
BioPython (x64) - http://www.lfd.uci.edu/~gohlke/pythonlibs *click on: biopython-1.63.win-amd64-py2.7.exe*

VII. Installing and Editing BOSToN
----------------------------------
BOSToN is stored and downloaded as a zip file and must be extracted and installed in the C:\Program directory for
links to operate properly. Once located in the C:\Program Files directory, BOSToN can be run in the Windows cmd line by
double left clicking on the BOSToN shortcut in the C:\Program Files\BOSToN_v1.0 folder. The BOSToN script can be run directly by 
double left clicking on the BOSToN.py executable in the C:\Program Files\BOSToN_v1.0\program folders or edit by right clicking
and opening using the Python IDLE. 

VIII. Instructions for Windows Standalone BOSToN
----------------------------------------------
After installing Python, NumPy and BioPython, the BOSToN script can be called in the BOSToN folder by clicking on the 
BOSToN icon. This opens a Window cmd line input which prompts user to enter enter FASTA sequences by calling a .txt 
file containing FASTA (and only FASTA) sequences. Some examples have been given for use and can be called by typing in
the following:
ex_pp.txt
ex_pq.txt
ex_qq.txt
ex_complex.txt

User FASTA inputs should be saved in a .txt file and called using its absolute or relative location. For example, if the FASTA
for the breakpoints is saved as "FASTA.txt" on your desktop, that FASTA input would be c:\Users\<username>\FASTA.txt, where 
<username> is the curent user's name (the Location of the file can be found by right clicking on the FASTA file on your desktop 
and selecting "Properties". 


Output will return the BLA(S)T Output Table for each FASTA input and the Next Generation Cyotgenomic Nomenclature (NGCN)
based on all the FASTA inputs (Note: non-FASTA formats and missing FASTA breakpoints may prevent the creation of the
BLA(S)T Ouput Tables and NGCN).
Example of BLA(S)T Output Table:
Query Start	Query End	Query Length	Query Strand	Identity	Sbjct Chromosome	Sbjct Band	Sbjct Strand	Sbjct Start	Sbjct End	Max Score	e-value
1		300		300	 	(+)		100.0%		17			q24.3	 	(-)		68248218	68247919	555.115		1.76031e-155
301		600		300		(+)		100.0%		3			p14.3	 	(+)		54949519	54949818	555.115		1.76031e-155
BLA(S)T Output: 17q24.3(-)(68247919)::3p14.3(+)(54949519)  

Example of NGCN:
seq[GRCh37] der(3)(17qter->17q24.3(68247919)::3p14.3(54949519)->3qter),der(17)(17pter->17q24.3(6824791{3-7})::3p14.3(549495{12-08})->3pter)

Program can be exited by clicking <return> key.

BOSToN can be re-initiated by clicking on the BOSToN icon in the BOSToN folder.

IX. Authors
----------
NGCN was developed by: 
Zehra Ordulu, M.D. (mordulu@partners.org)

BOSToN algorith was developed by:
Benjamin B. Currall, Ph.D. (bcurrall@partners.org)
Andrew Ivanov (andreyivanov88@gmail.com)
Lei Xia (leisummer@gmail.com)

Developing Laboratory:
Cynthia C. Morton, Ph.D. (cmorton@partners.org)

X. Release Notes
-----------------
NGNC - was published online on April 17, 2014 by Ordulu et al. in an article named "Describing Sequencing Results 
of Structural Chromosome Rearrangements with a Suggested Next-Generation Cytogenetic Nomenclature" in the 
American Journal of Human Genetics

BOSToN v1.0 - released April 28, 2014

XI. No Warranties
-----------------  
THE SOFTWARE IS PROVIDED “AS  IS” AND “AS  AVAILABLE”  TO THE FULLEST EXTENT PERMITTED BY LAW, BWH DISCLAIMS ALL 
WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO ANY IMPLIED WARRANTIES OF TITLE, 
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT. 

