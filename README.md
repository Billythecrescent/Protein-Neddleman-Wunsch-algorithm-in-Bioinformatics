# Protein-Neddleman-Wunsch-algorithm-in-Bioinformatics
This project is a Needleman-Wunsch algorithm practice in protein sequence alignment in Bioinformatics.
---
README.me for NW.cpp

Function: Align two protein sequences using Needleman-Wunsch algorithm.

Usage: After compiling, .exe file is available.
---
You have two options:

When you want to manually input the sequence:

'''
NW.exe -e sequence1 sequence2 matrix >> ../result/result.txt
'''

sequence1, sequence2 are the protein sequences to be aligned, while the matrix is the align matrix used, such as BLOSUM45.

>> ../result/result.txt means that the result is redirected to the result.txt in the result directory.

For example: 

'''
NW.exe -i -i KPL63462.1.FASTA KZL42154.1.FASTA BLOSUM45 >> ../result/result.txt
'''

When you want to put in two protein FASTA file:

'''
NW.exe -i file1 file2 matrix >> ../result/result.txt
'''

**Note: FASTA file must be put in the sequence directory before usage.**

**Muti-sequences in one FASTA file is supported**

Matrix avaliable: BLOSUM45, BLOSUM90, BLOSUM62, PAM100  >> ../result/result.txt

**Make sure you have the matrix file in the `lib` directory before you use the matrix name**

(PS: test.bat in "code" directory .bat file for test, result is in result.txt in "result" directory)
