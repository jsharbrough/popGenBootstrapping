# -*- coding: utf-8 -*-
import sys
import random
# -*- coding: utf-8 -*-
def polymorphismsSubstitutions(fofn,bootstrapNum):
    speciesDict = {'antipodarum':['>Ianthe','>Lady','>AlexMap','>Alexsex','>Yellow_Contig_56','>Ianthe_lane1_TAAGGCGA_trimmed_(paired)_contig_309','>Kaniere_1','>Rotoroa_1'],'estuarinus':['>Potamopyrgus_estuarinus']}
    referenceDict = {'antipodarum':">Yellow_Contig_56",'estuarinus':">Potamopyrgus_estuarinus"}
    speciesList = ['antipodarum','estuarinus']
    geneticCode = {'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'M', 'AGG': 'S', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'S', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'}
    alignmentFiles = open(fofn,'r')
    for alignment in alignmentFiles:
        seqDict = {}
        seqList = []
        currSeq = ''
        seqName = ''
        alignFile = open(alignment[0:-1],'r')
        for line in alignFile:
            if line[0] == '>':
                if seqName != '':
                    seqDict[seqName] = currSeq
                seqName = line
                while seqName[-1] == '\n' or seqName[-1] == '\t' or seqName[-1] == '\r':
                    seqName = seqName[0:-1]
       	        seqList.append(seqName)
       	        currSeq = ''
            else:	
                currSeq += line
                while currSeq[-1] == '\n' or currSeq[-1] == '\t' or currSeq[-1] == '\r':
                    currSeq = currSeq[0:-1]
        seqDict[seqName] = currSeq
        alignFile.close()
        currNum = 0
        sys.stdout.write('thetaS' + '\t' + 'thetaA' + '\t' + 'thetaC1' + '\t' + 'thetaR1' + '\t' + 'thetaC2' + '\t' + 'thetaR2' + '\t' + 'thetaC3' + '\t' + 'thetaR3' + '\t' + 'thetaC4' + '\t' + 'thetaR4' + '\t' + 'thetaC5' + '\t' + 'thetaR5' + '\t' + 'thetaC6' + '\t' + 'thetaR6' + '\t' + 'thetaC7' + '\t' + 'thetaR7' + '\t' + 'thetaMeanC' + '\t' + 'thetaMeanR' + '\t' + 'piS' + '\t' + 'piA' + '\t' + 'piC1' + '\t' + 'piR1' + '\t' + 'piC2' + '\t' + 'piR2' + '\t' + 'piC3' + '\t' + 'piR3' + '\t' + 'piC4' + '\t' + 'piR4' + '\t' + 'piC5' + '\t' + 'piR5' + '\t' + 'piC6' + '\t' + 'piR6' + '\t' + 'piC7' + '\t' + 'piR7' + '\t' + 'piMeanC' + '\t' + 'piMeanR' + '\n')
        codonDict,AADict = buildCodonDict(seqDict,seqList)
        while currNum < bootstrapNum:
            currSeqDict = seqDict
            #currSeqDict = bootstrapReplicate(codonDict,seqList)
            thetaS,thetaA,thetaC1,thetaR1,thetaC2,thetaR2,thetaC3,thetaR3,thetaC4,thetaR4,thetaC5,thetaR5,thetaC6,thetaR6,thetaC7,thetaR7,thetaMeanC,thetaMeanR,piS,piA,piC1,piR1,piC2,piR2,piC3,piR3,piC4,piR4,piC5,piR5,piC6,piR6,piC7,piR7,piMeanC,piMeanR = polSub(currSeqDict,seqList)
            if thetaS !=0 and piS != 0:
                sys.stdout.write(str(thetaS) + '\t' + str(thetaA) + '\t' + str(thetaC1) + '\t' + str(thetaR1) + '\t' + str(thetaC2) + '\t' + str(thetaR2) + '\t' + str(thetaC3) + '\t' + str(thetaR3) + '\t' + str(thetaC4) + '\t' + str(thetaR4) + '\t' + str(thetaC5) + '\t' + str(thetaR5) + '\t' + str(thetaC6) + '\t' + str(thetaR6) + '\t' + str(thetaC7) + '\t' + str(thetaR7) + '\t' + str(thetaMeanC) + '\t' + str(thetaMeanR) + '\t' + str(piS) + '\t' + str(piA) + '\t' + str(piC1) + '\t' + str(piR1) + '\t' + str(piC2) + '\t' + str(piR2) + '\t' + str(piC3) + '\t' + str(piR3) + '\t' + str(piC4) + '\t' + str(piR4) + '\t' + str(piC5) + '\t' + str(piR5) + '\t' + str(piC6) + '\t' + str(piR6) + '\t' + str(piC7) + '\t' + str(piR7) + '\t' + str(piMeanC) + '\t' + str(piMeanR) + '\n')
                currNum += 1

def uniquePolymorphismsSubstitutions(fofn,bootstrapNum):
    speciesDict = {'antipodarum':['>Ianthe','>Lady','>AlexMap','>Alexsex','>Yellow_Contig_56','>Ianthe_lane1_TAAGGCGA_trimmed_(paired)_contig_309','>Kaniere_1','>Rotoroa_1'],'estuarinus':['>Potamopyrgus_estuarinus']}
    referenceDict = {'antipodarum':">Yellow_Contig_56",'estuarinus':">Potamopyrgus_estuarinus"}
    speciesList = ['antipodarum','estuarinus']
    geneticCode = {'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'M', 'AGG': 'S', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'S', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'}
    alignmentFiles = open(fofn,'r')
    for alignment in alignmentFiles:
        seqDict = {}
        seqList = []
        currSeq = ''
        seqName = ''
        alignFile = open(alignment[0:-1],'r')
        for line in alignFile:
            if line[0] == '>':
                if seqName != '':
                    seqDict[seqName] = currSeq
                seqName = line
                while seqName[-1] == '\n' or seqName[-1] == '\t' or seqName[-1] == '\r':
                    seqName = seqName[0:-1]
       	        seqList.append(seqName)
       	        currSeq = ''
            else:	
                currSeq += line
                while currSeq[-1] == '\n' or currSeq[-1] == '\t' or currSeq[-1] == '\r':
                    currSeq = currSeq[0:-1]
        seqDict[seqName] = currSeq
        alignFile.close()
        currNum = 0
        sys.stdout.write('thetaUS' + '\t' + 'thetaUA' + '\t' + 'thetaUC1' + '\t' + 'thetaUR1' + '\t' + 'thetaUC2' + '\t' + 'thetaUR2' + '\t' + 'thetaUC3' + '\t' + 'thetaUR3' + '\t' + 'thetaUC4' + '\t' + 'thetaUR4' + '\t' + 'thetaUC5' + '\t' + 'thetaUR5' + '\t' + 'thetaUC6' + '\t' + 'thetaUR6' + '\t' + 'thetaUC7' + '\t' + 'thetaUR7' + '\t' + 'thetaUMeanC' + '\t' + 'thetaUMeanR' + '\n')
        codonDict,AADict = buildCodonDict(seqDict,seqList)
        while currNum < bootstrapNum:
            currSeqDict = bootstrapReplicate(codonDict,seqList)
            #currSeqDict = seqDict
            thetaUS,thetaUA,thetaUC1,thetaUR1,thetaUC2,thetaUR2,thetaUC3,thetaUR3,thetaUC4,thetaUR4,thetaUC5,thetaUR5,thetaUC6,thetaUR6,thetaUC7,thetaUR7,thetaUMeanC,thetaUMeanR = thetaU(currSeqDict,seqList)
            if thetaUS !=0:
                sys.stdout.write(str(thetaUS) + '\t' + str(thetaUA) + '\t' + str(thetaUC1) + '\t' + str(thetaUR1) + '\t' + str(thetaUC2) + '\t' + str(thetaUR2) + '\t' + str(thetaUC3) + '\t' + str(thetaUR3) + '\t' + str(thetaUC4) + '\t' + str(thetaUR4) + '\t' + str(thetaUC5) + '\t' + str(thetaUR5) + '\t' + str(thetaUC6) + '\t' + str(thetaUR6) + '\t' + str(thetaUC7) + '\t' + str(thetaUR7) + '\t' + str(thetaUMeanC) + '\t' + str(thetaUMeanR) + '\n')
                currNum += 1

def aN(N):
    aN = 0.0
    i = 1
    while i < N:
        aN += 1.0/i
        i += 1
    return aN
        
def polSub(seqDict,seqList,code='invertebrateMt'):
    geneticCodes = {'standard':{"TTT":"F",	"TTC":"F",	"TTA":"L",	"TTG":"L",	"TCT":"S",	"TCC":"S",	"TCA":"S",	"TCG":"S",	"TAT":"Y",	"TAC":"Y",	"TAA":"*",	"TAG":"*",	"TGT":"C",	"TGC":"C",	"TGA":"*",	"TGG":"W",	"CTT":"L",	"CTC":"L",	"CTA":"L",	"CTG":"L",	"CCT":"P",	"CCC":"P",	"CCA":"P",	"CCG":"P",	"CAT":"H",	"CAC":"H",	"CAA":"Q",	"CAG":"Q",	"CGT":"R",	"CGC":"R",	"CGA":"R",	"CGG":"R",	"ATT":"I",	"ATC":"I",	"ATA":"I",	"ATG":"M",	"ACT":"T",	"ACC":"T",	"ACA":"T",	"ACG":"T",	"AAT":"N",	"AAC":"N",	"AAA":"K",	"AAG":"K",	"AGT":"S",	"AGC":"S",	"AGA":"R",	"AGG":"R",	"GTT":"V",	"GTC":"V",	"GTA":"V",	"GTG":"V",	"GCT":"A",	"GCC":"A",	"GCA":"A",	"GCG":"A",	"GAT":"D",	"GAC":"D",	"GAA":"E",	"GAG":"E",	"GGT":"G",	"GGC":"G",	"GGA":"G",	"GGG":"G"},'invertebrateMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'M', 'AGG': 'S', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'S', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'vertebrateMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'M', 'AGG': '*', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': '*', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'}, 'yeastMt':{'CTT': 'T', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'M', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'T', 'CTA': 'T', 'CTC': 'T', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'}, 'coelenterateMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'ciliateNuc':{'CTT': 'L', 'TAG': 'Q', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': '*', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': 'Q', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'echinodermMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'S', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'S', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'N', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'}, 'euplotidNuc':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'C', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'}, 'bacterial':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': '*', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'yeastNuc':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'S', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': '*', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'}, 'ascidianMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'M', 'AGG': 'G', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'G', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'flatwormMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'S', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'S', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'N', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': 'Y', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'chlorophyceanMt':{'CTT': 'L', 'TAG': 'L', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': '*', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'trematodeMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'M', 'AGG': 'S', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'S', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'N', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'pterobranchiaMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'K', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'S', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'}}
    geneticCode = geneticCodes[code]
    startCodons = ['ATT','ATC','ATA','ATG','GTG'] #invertebrateMt code
    codonDict,AADict = buildCodonDict(seqDict,seqList)
    popList = []
    outList = [] 
    synL = 0.0
    nsynL = 0.0
    con1L = 0.0 
    rad1L = 0.0
    con2L = 0.0
    rad2L = 0.0
    con3L = 0.0
    rad3L = 0.0
    con4L = 0.0
    rad4L = 0.0
    con5L = 0.0
    rad5L = 0.0
    con6L = 0.0
    rad6L = 0.0
    con7L = 0.0
    rad7L = 0.0
    meanConL = 0.0
    meanRadL = 0.0
    for seq in seqList:
        currCodonList = codonDict[seq]
        totalSynSites,totalNonsynSites,totalC1Sites,totalR1Sites,totalC2Sites,totalR2Sites,totalC3Sites,totalR3Sites,totalC4Sites,totalR4Sites,totalC5Sites,totalR5Sites,totalC6Sites,totalR6Sites,totalC7Sites,totalR7Sites,totalMeanCSites,totalMeanRSites = countSites(currCodonList)
        if seq == '>Potamopyrgus_estuarinus':
            outList.append(seq)
        else:
            popList.append(seq)
            synL += totalSynSites
            nsynL += totalNonsynSites
            con1L += totalC1Sites
            rad1L += totalR1Sites
            con2L += totalC2Sites
            rad2L += totalR2Sites
            con3L += totalC3Sites
            rad3L += totalR3Sites
            con4L += totalC4Sites
            rad4L += totalR4Sites
            con5L += totalC5Sites
            rad5L += totalR5Sites
            con6L += totalC6Sites
            rad6L += totalR6Sites
            con7L += totalC7Sites
            rad7L += totalR7Sites
            meanConL += totalMeanCSites
            meanRadL += totalMeanRSites          
    N = len(popList)
    An = aN(N)
    synL /= N
    nsynL /= N
    con1L /= N 
    rad1L /= N
    con2L /= N
    rad2L /= N
    con3L /= N
    rad3L /= N
    con4L /= N
    rad4L /= N
    con5L /= N
    rad5L /= N
    con6L /= N
    rad6L /= N
    con7L /= N
    rad7L /= N
    meanConL /= N
    meanRadL /= N
    popPolymorphicSites = []
    i = 0
    sum2PQ_S = 0
    sum2PQ_N = 0
    sum2PQ_C1 = 0
    sum2PQ_C2 = 0
    sum2PQ_C3 = 0
    sum2PQ_C4 = 0
    sum2PQ_C5 = 0
    sum2PQ_C6 = 0
    sum2PQ_C7 = 0
    sum2PQ_R1 = 0
    sum2PQ_R2 = 0
    sum2PQ_R3 = 0
    sum2PQ_R4 = 0
    sum2PQ_R5 = 0
    sum2PQ_R6 = 0
    sum2PQ_R7 = 0
    sum2PQ_meanC = 0
    sum2PQ_meanR = 0
    synS = 0
    nsynS = 0
    con1S = 0
    con2S = 0
    con3S = 0
    con4S = 0
    con5S = 0
    con6S = 0
    con7S = 0
    meanConS = 0
    rad1S = 0
    rad2S = 0
    rad3S = 0
    rad4S = 0
    rad5S = 0
    rad6S = 0
    rad7S = 0
    meanRadS = 0
    synUS = 0.0
    nsynUS = 0.0
    con1US = 0.0
    rad1US = 0.0
    con2US = 0.0
    rad2US = 0.0
    con3US = 0.0
    rad3US = 0.0
    con4US = 0.0
    rad4US = 0.0
    con5US = 0.0
    rad5US = 0.0
    con6US = 0.0
    rad6US = 0.0
    con7US = 0.0
    rad7US = 0.0
    meanConUS = 0.0
    meanRadUS = 0.0
    outCodons = codonDict[outList[0]]
    while i < len(codonDict[seqList[0]]):
        currAlleleDict = {}
        currAlleleList = []
        currAADict = {}
        outCodon = outCodons[i]
        for seq in popList:
            currCodons = codonDict[seq]
            currCodon = currCodons[i]
            if currCodon not in currAlleleDict and 'N' not in currCodon and '-' not in currCodon:
                currAlleleDict[currCodon] = 1
                currAlleleList.append(currCodon)
            elif 'N' not in currCodon and '-' not in currCodon:
                currValue = currAlleleDict[currCodon]
                currValue += 1
                currAlleleDict[currCodon] = currValue
        if len(currAlleleDict) > 1:
            totalIndividuals = 0
            site1 = []
            site2 = []
            site3 = []
            for codon in currAlleleList:
                totalIndividuals += currAlleleDict[codon]
                if codon[0] not in site1:
                    site1.append(codon[0])
                if codon[1] not in site2:
                    site2.append(codon[1])
                if codon[2] not in site3:
                    site3.append(codon[2])
            currFreqDict = {}
            totalChanges = (len(site1) - 1) + (len(site2) - 1) + (len(site3) - 1)
            variableSites = []
            if len(site1) > 1:
                popPolymorphicSites.append(i*3)
                variableSites.append(i*3)
            if len(site2) > 1:
                popPolymorphicSites.append((i*3) + 1)
                variableSites.append((i*3) + 1)
            if len(site3) > 1:
                popPolymorphicSites.append((i*3) + 2)
                variableSites.append((i*3) + 1)
            aaList = []
            twoPQ = 2
            singleton = False
            for codon in currAlleleDict:
                freq = float(currAlleleDict[codon])/totalIndividuals
                currFreqDict[codon] = freq
                if currAlleleDict[codon] == 1:
                    singleton = True
                if i == 0 and codon in startCodons:
                    aa = 'M'
                else:
                    aa = geneticCode[codon]
                currAADict[codon] = aa
                if aa not in aaList:
                    aaList.append(aa)
            if totalChanges == 1:
                for codon in currAlleleDict:
                    freq = currFreqDict[codon]
                    twoPQ *= freq
                if len(aaList) == 1:
                    synS += 1
                    sum2PQ_S += twoPQ
                    if singleton == True:
                        synUS += 1
                else:
                    if singleton == True:
                        nsynUS += 1
                    nsynS += 1
                    sum2PQ_N += twoPQ
                    mutType = CRI(aaList) #[1,2,3,4,5,6,7,cri]
                    if mutType[0] == 0:
                        con1S += 1
                        sum2PQ_C1 += twoPQ
                        if singleton == True:
                            con1US += 1
                    else:
                        rad1S += 1
                        sum2PQ_R1 += twoPQ
                        if singleton == True:
                            rad1US += 1
                    if mutType[1] == 0:
                        con2S += 1
                        sum2PQ_C2 += twoPQ
                        if singleton == True:
                            con2US += 1
                    else:
                        rad2S += 1
                        sum2PQ_R2 += twoPQ
                        if singleton == True:
                            rad2US += 1
                    if mutType[2] == 0:
                        con3S += 1
                        sum2PQ_C3 += twoPQ
                        if singleton == True:
                            con3US += 1
                    else:
                        rad3S += 1
                        sum2PQ_R3 += twoPQ
                        if singleton == True:
                            rad3US += 1
                    if mutType[3] == 0:
                        con4S += 1
                        sum2PQ_C4 += twoPQ
                        if singleton == True:
                            con4US += 1
                    else:
                        rad4S += 1
                        sum2PQ_R4 += twoPQ
                        if singleton == True:
                            rad4US += 1
                    if mutType[4] == 0:
                        con5S += 1
                        sum2PQ_C5 += twoPQ
                        if singleton == True:
                            con5US += 1
                    else:
                        rad5S += 1
                        sum2PQ_R5 += twoPQ
                        if singleton == True:
                            rad5US += 1
                    if mutType[5] == 0:
                        con6S += 1
                        sum2PQ_C6 += twoPQ
                        if singleton == True:
                            con6US += 1
                    else:
                        rad6S += 1
                        sum2PQ_R6 += twoPQ
                        if singleton == True:
                            rad6US += 1
                    if mutType[6] == 0:
                        con7S += 1
                        sum2PQ_C7 += twoPQ
                        if singleton == True:
                            con7US += 1
                    else:
                        rad7S += 1
                        sum2PQ_R7 += twoPQ
                        if singleton == True:
                            rad7US += 1
                    if mutType[7] <= 0.5:
                        meanConS += 1
                        sum2PQ_meanC += twoPQ
                        if singleton == True:
                            meanConUS += 1
                    else:
                        meanRadS += 1
                        sum2PQ_meanR += twoPQ
                        if singleton == True:
                            meanRadUS += 1
            elif totalChanges == 2:
                if len(currAlleleDict) == 3:
                    ab = 0
                    ac = 0
                    bc = 0
                    codonA = currAlleleList[0]
                    codonB = currAlleleList[1]
                    codonC = currAlleleList[2]
                    if codonA[0] != codonB[0]:
                        ab += 1
                    if codonA[1] != codonB[1]:
                        ab += 1
                    if codonA[2] != codonB[2]:
                        ab += 1
                    if codonA[0] != codonC[0]:
                        ac += 1
                    if codonA[1] != codonC[1]:
                        ac += 1
                    if codonA[2] != codonC[2]:
                        ac += 1
                    if codonC[0] != codonB[0]:
                        bc += 1
                    if codonC[1] != codonB[1]:
                        bc += 1
                    if codonC[2] != codonB[2]:
                        bc += 1
                    if ab == ac and ac == bc:
                        if 'N' not in outCodon and '-' not in outCodon:
                            if outCodon == codonA:
                                aaList1 = [currAADict[codonA],currAADict[codonB]]
                                aaList2 = [currAADict[codonA],currAADict[codonC]]
                                codonList1 = [codonA,codonB]
                                codonList2 = [codonA,codonC]
                                if aaList1[0] == aaList1[1]:
                                    if aaList2[0] == aaList2[1]:
                                        synS += 2
                                        twoPQ = 4
                                        for allele in currFreqDict:
                                            twoPQ *= currFreqDict[allele]
                                        sum2PQ_S += twoPQ
                                        if singleton == True:
                                            synUS += 2
                                    else:
                                        twoPQ1 = 2*currFreqDict[codonList1[1]]*(1-currFreqDict[codonList1[1]]) #syn 
                                        twoPQ2 = 2*currFreqDict[codonList2[1]]*(1-currFreqDict[codonList2[1]]) #nsyn
                                        sum2PQ_S += twoPQ1
                                        sum2PQ_N += twoPQ2
                                        synS += 1
                                        nsynS += 1
                                        if singleton == True:
                                            synUS += 1
                                            nsynUS += 1
                                        mutType = CRI(aaList2) #[1,2,3,4,5,6,7,cri]
                                        if mutType[0] == 0:
                                            con1S += 1
                                            sum2PQ_C1 += twoPQ2
                                            if singleton == True:
                                                con1US += 1
                                        else:
                                            rad1S += 1
                                            sum2PQ_R1 += twoPQ2
                                            if singleton == True:
                                                rad1US += 1
                                        if mutType[1] == 0:
                                            con2S += 1
                                            sum2PQ_C2 += twoPQ2
                                            if singleton == True:
                                                con2US += 1
                                        else:
                                            rad2S += 1
                                            sum2PQ_R2 += twoPQ2
                                            if singleton == True:
                                                rad2US += 1
                                        if mutType[2] == 0:
                                            con3S += 1
                                            sum2PQ_C3 += twoPQ2
                                            if singleton == True:
                                                con3US += 1
                                        else:
                                            rad3S += 1
                                            sum2PQ_R3 += twoPQ2
                                            if singleton == True:
                                                rad3US += 1
                                        if mutType[3] == 0:
                                            con4S += 1
                                            sum2PQ_C4 += twoPQ2
                                            if singleton == True:
                                                con4US += 1
                                        else:
                                            rad4S += 1
                                            sum2PQ_R4 += twoPQ2
                                            if singleton == True:
                                                rad4US += 1
                                        if mutType[4] == 0:
                                            con5S += 1
                                            sum2PQ_C5 += twoPQ2
                                            if singleton == True:
                                                con5US += 1
                                        else:
                                            rad5S += 1
                                            sum2PQ_R5 += twoPQ2
                                            if singleton == True:
                                                rad5US += 1
                                        if mutType[5] == 0:
                                            con6S += 1
                                            sum2PQ_C6 += twoPQ2
                                            if singleton == True:
                                                con6US += 1
                                        else:
                                            rad6S += 1
                                            sum2PQ_R6 += twoPQ2
                                            if singleton == True:
                                                rad6US += 1
                                        if mutType[6] == 0:
                                            con7S += 1
                                            sum2PQ_C7 += twoPQ2
                                            if singleton == True:
                                                con7US += 1
                                        else:
                                            rad7S += 1
                                            sum2PQ_R7 += twoPQ2
                                            if singleton == True:
                                                rad7US += 1
                                        if mutType[7] <= 0.5:
                                            meanConS += 1
                                            sum2PQ_meanC += twoPQ2
                                            if singleton == True:
                                                meanConUS += 1
                                        else:
                                            meanRadS += 1
                                            sum2PQ_meanR += twoPQ2
                                            if singleton == True:
                                                meanRadUS += 1
                                elif aaList2[0] == aaList2[1]:
                                    nsynS += 1
                                    synS += 1
                                    if singleton == True:
                                        synUS += 1
                                        nsynUS += 1
                                    twoPQ1 = 2*currFreqDict[codonList1[1]]*(1-currFreqDict[codonList1[1]]) #nsyn 
                                    twoPQ2 = 2*currFreqDict[codonList2[1]]*(1-currFreqDict[codonList2[1]]) #syn
                                    sum2PQ_S += twoPQ2
                                    sum2PQ_N += twoPQ1
                                    mutType = CRI(aaList1) #[1,2,3,4,5,6,7,cri]
                                    if mutType[0] == 0:
                                        con1S += 1
                                        sum2PQ_C1 += twoPQ1
                                        if singleton == True:
                                            con1US += 1
                                    else:
                                        rad1S += 1
                                        sum2PQ_R1 += twoPQ1
                                        if singleton == True:
                                            rad1US += 1
                                    if mutType[1] == 0:
                                        con2S += 1
                                        sum2PQ_C2 += twoPQ1
                                        if singleton == True:
                                            con2US += 1
                                    else:
                                        rad2S += 1
                                        sum2PQ_R2 += twoPQ1
                                        if singleton == True:
                                            rad2US += 1
                                    if mutType[2] == 0:
                                        con3S += 1
                                        sum2PQ_C3 += twoPQ1
                                        if singleton == True:
                                            con3US += 1
                                    else:
                                        rad3S += 1
                                        sum2PQ_R3 += twoPQ1
                                        if singleton == True:
                                            rad3US += 1
                                    if mutType[3] == 0:
                                        con4S += 1
                                        sum2PQ_C4 += twoPQ1
                                        if singleton == True:
                                            con4US += 1
                                    else:
                                        rad4S += 1
                                        sum2PQ_R4 += twoPQ1
                                        if singleton == True:
                                            rad4US += 1
                                    if mutType[4] == 0:
                                        con5S += 1
                                        sum2PQ_C5 += twoPQ1
                                        if singleton == True:
                                            con5US += 1
                                    else:
                                        rad5S += 1
                                        sum2PQ_R5 += twoPQ1
                                        if singleton == True:
                                            rad5US += 1
                                    if mutType[5] == 0:
                                        con6S += 1
                                        sum2PQ_C6 += twoPQ1
                                        if singleton == True:
                                            con6US += 1
                                    else:
                                        rad6S += 1
                                        sum2PQ_R6 += twoPQ1
                                        if singleton == True:
                                            rad6US += 1
                                    if mutType[6] == 0:
                                        con7S += 1
                                        sum2PQ_C7 += twoPQ1
                                        if singleton == True:
                                            con7US += 1
                                    else:
                                        rad7S += 1
                                        sum2PQ_R7 += twoPQ1
                                        if singleton == True:
                                            rad7US += 1
                                    if mutType[7] <= 0.5:
                                        meanConS += 1
                                        sum2PQ_meanC += twoPQ1
                                        if singleton == True:
                                            meanConUS += 1
                                    else:
                                        meanRadS += 1
                                        sum2PQ_meanR += twoPQ1
                                        if singleton == True:
                                            meanRadUS += 1
                                else:
                                    nsynS += 2
                                    if singleton == True:
                                        nsynUS += 2
                                    twoPQ1 = 2*currFreqDict[codonList1[1]]*(1-currFreqDict[codonList1[1]]) #nsyn 
                                    twoPQ2 = 2*currFreqDict[codonList2[1]]*(1-currFreqDict[codonList2[1]]) #nsyn
                                    sum2PQ_N += twoPQ1 + twoPQ2
                                    mutType = CRI(aaList1) #[1,2,3,4,5,6,7,cri]
                                    if mutType[0] == 0:
                                        con1S += 1
                                        sum2PQ_C1 += twoPQ1
                                        if singleton == True:
                                            con1US += 1
                                    else:
                                        rad1S += 1
                                        sum2PQ_R1 += twoPQ1
                                        if singleton == True:
                                            rad1US += 1
                                    if mutType[1] == 0:
                                        con2S += 1
                                        sum2PQ_C2 += twoPQ1
                                        if singleton == True:
                                            con2US += 1
                                    else:
                                        rad2S += 1
                                        sum2PQ_R2 += twoPQ1
                                        if singleton == True:
                                            rad2US += 1
                                    if mutType[2] == 0:
                                        con3S += 1
                                        sum2PQ_C3 += twoPQ1
                                        if singleton == True:
                                            con3US += 1
                                    else:
                                        rad3S += 1
                                        sum2PQ_R3 += twoPQ1
                                        if singleton == True:
                                            rad3US += 1
                                    if mutType[3] == 0:
                                        con4S += 1
                                        sum2PQ_C4 += twoPQ1
                                        if singleton == True:
                                            con4US += 1
                                    else:
                                        rad4S += 1
                                        sum2PQ_R4 += twoPQ1
                                        if singleton == True:
                                            rad4US += 1
                                    if mutType[4] == 0:
                                        con5S += 1
                                        sum2PQ_C5 += twoPQ1
                                        if singleton == True:
                                            con5US += 1
                                    else:
                                        rad5S += 1
                                        sum2PQ_R5 += twoPQ1
                                        if singleton == True:
                                            rad5US += 1
                                    if mutType[5] == 0:
                                        con6S += 1
                                        sum2PQ_C6 += twoPQ1
                                        if singleton == True:
                                            con6US += 1
                                    else:
                                        rad6S += 1
                                        sum2PQ_R6 += twoPQ1
                                        if singleton == True:
                                            rad6US += 1
                                    if mutType[6] == 0:
                                        con7S += 1
                                        sum2PQ_C7 += twoPQ1
                                        if singleton == True:
                                            con7US += 1
                                    else:
                                        rad7S += 1
                                        sum2PQ_R7 += twoPQ1
                                        if singleton == True:
                                            rad7US += 1
                                    if mutType[7] <= 0.5:
                                        meanConS += 1
                                        sum2PQ_meanC += twoPQ1
                                        if singleton == True:
                                            meanConUS += 1
                                    else:
                                        meanRadS += 1
                                        sum2PQ_meanR += twoPQ1
                                        if singleton == True:
                                            meanRadUS += 1
                                    mutType = CRI(aaList2) #[1,2,3,4,5,6,7,cri]
                                    if mutType[0] == 0:
                                        con1S += 1
                                        sum2PQ_C1 += twoPQ2
                                        if singleton == True:
                                            con1US += 1
                                    else:
                                        rad1S += 1
                                        sum2PQ_R1 += twoPQ2
                                        if singleton == True:
                                            rad1US += 1
                                    if mutType[1] == 0:
                                        con2S += 1
                                        sum2PQ_C2 += twoPQ2
                                        if singleton == True:
                                            con2US += 1
                                    else:
                                        rad2S += 1
                                        sum2PQ_R2 += twoPQ2
                                        if singleton == True:
                                            rad2US += 1
                                    if mutType[2] == 0:
                                        con3S += 1
                                        sum2PQ_C3 += twoPQ2
                                        if singleton == True:
                                            con3US += 1
                                    else:
                                        rad3S += 1
                                        sum2PQ_R3 += twoPQ2
                                        if singleton == True:
                                            rad3US += 1
                                    if mutType[3] == 0:
                                        con4S += 1
                                        sum2PQ_C4 += twoPQ2
                                        if singleton == True:
                                            con4US += 1
                                    else:
                                        rad4S += 1
                                        sum2PQ_R4 += twoPQ2
                                        if singleton == True:
                                            rad4US += 1
                                    if mutType[4] == 0:
                                        con5S += 1
                                        sum2PQ_C5 += twoPQ2
                                        if singleton == True:
                                            con5US += 1
                                    else:
                                        rad5S += 1
                                        sum2PQ_R5 += twoPQ2
                                        if singleton == True:
                                            rad5US += 1
                                    if mutType[5] == 0:
                                        con6S += 1
                                        sum2PQ_C6 += twoPQ2
                                        if singleton == True:
                                            con6US += 1
                                    else:
                                        rad6S += 1
                                        sum2PQ_R6 += twoPQ2
                                        if singleton == True:
                                            rad6US += 1
                                    if mutType[6] == 0:
                                        con7S += 1
                                        sum2PQ_C7 += twoPQ2
                                        if singleton == True:
                                            con7US += 1
                                    else:
                                        rad7S += 1
                                        sum2PQ_R7 += twoPQ2
                                        if singleton == True:
                                            rad7US += 1
                                    if mutType[7] <= 0.5:
                                        meanConS += 1
                                        sum2PQ_meanC += twoPQ2
                                        if singleton == True:
                                            meanConUS += 1
                                    else:
                                        meanRadS += 1
                                        sum2PQ_meanR += twoPQ2
                                        if singleton == True:
                                            meanRadUS += 1
                            elif outCodon == codonB:
                                aaList1 = [currAADict[codonB],currAADict[codonA]]
                                aaList2 = [currAADict[codonB],currAADict[codonC]]
                                codonList1 = [codonB,codonA]
                                codonList2 = [codonB,codonC]
                                if aaList1[0] == aaList1[1]:
                                    if aaList2[0] == aaList2[1]:
                                        synS += 2
                                        if singleton == True:
                                            synUS += 1
                                        twoPQ = 4
                                        for allele in currFreqDict:
                                            twoPQ *= currFreqDict[allele]
                                        sum2PQ_S += twoPQ
                                    else:
                                        twoPQ1 = 2*currFreqDict[codonList1[1]]*(1-currFreqDict[codonList1[1]]) #syn 
                                        twoPQ2 = 2*currFreqDict[codonList2[1]]*(1-currFreqDict[codonList2[1]]) #nsyn
                                        sum2PQ_S += twoPQ1
                                        sum2PQ_N += twoPQ2
                                        synS += 1
                                        nsynS += 1
                                        if singleton == True:
                                            synUS += 1
                                            nsynUS += 1
                                        mutType = CRI(aaList2) #[1,2,3,4,5,6,7,cri]
                                        if mutType[0] == 0:
                                            con1S += 1
                                            sum2PQ_C1 += twoPQ2
                                            if singleton == True:
                                                con1US += 1
                                        else:
                                            rad1S += 1
                                            sum2PQ_R1 += twoPQ2
                                            if singleton == True:
                                                rad1US += 1
                                        if mutType[1] == 0:
                                            con2S += 1
                                            sum2PQ_C2 += twoPQ2
                                            if singleton == True:
                                                con2US += 1
                                        else:
                                            rad2S += 1
                                            sum2PQ_R2 += twoPQ2
                                            if singleton == True:
                                                rad2US += 1
                                        if mutType[2] == 0:
                                            con3S += 1
                                            sum2PQ_C3 += twoPQ2
                                            if singleton == True:
                                                con3US += 1
                                        else:
                                            rad3S += 1
                                            sum2PQ_R3 += twoPQ2
                                            if singleton == True:
                                                rad3US += 1
                                        if mutType[3] == 0:
                                            con4S += 1
                                            sum2PQ_C4 += twoPQ2
                                            if singleton == True:
                                                con4US += 1
                                        else:
                                            rad4S += 1
                                            sum2PQ_R4 += twoPQ2
                                            if singleton == True:
                                                rad4US += 1
                                        if mutType[4] == 0:
                                            con5S += 1
                                            sum2PQ_C5 += twoPQ2
                                            if singleton == True:
                                                con5US += 1
                                        else:
                                            rad5S += 1
                                            sum2PQ_R5 += twoPQ2
                                            if singleton == True:
                                                rad5US += 1
                                        if mutType[5] == 0:
                                            con6S += 1
                                            sum2PQ_C6 += twoPQ2
                                            if singleton == True:
                                                con6US += 1
                                        else:
                                            rad6S += 1
                                            sum2PQ_R6 += twoPQ2
                                            if singleton == True:
                                                rad6US += 1
                                        if mutType[6] == 0:
                                            con7S += 1
                                            sum2PQ_C7 += twoPQ2
                                            if singleton == True:
                                                con7US += 1
                                        else:
                                            rad7S += 1
                                            sum2PQ_R7 += twoPQ2
                                            if singleton == True:
                                                rad7US += 1
                                        if mutType[7] <= 0.5:
                                            meanConS += 1
                                            sum2PQ_meanC += twoPQ2
                                            if singleton == True:
                                                meanConUS += 1
                                        else:
                                            meanRadS += 1
                                            sum2PQ_meanR += twoPQ2
                                            if singleton == True:
                                                meanRadUS += 1
                                elif aaList2[0] == aaList2[1]:
                                    nsynS += 1
                                    synS += 1
                                    if singleton == True:
                                        synUS += 1
                                        nsynUS += 1
                                    twoPQ1 = 2*currFreqDict[codonList1[1]]*(1-currFreqDict[codonList1[1]]) #nsyn 
                                    twoPQ2 = 2*currFreqDict[codonList2[1]]*(1-currFreqDict[codonList2[1]]) #syn
                                    sum2PQ_S += twoPQ2
                                    sum2PQ_N += twoPQ1
                                    mutType = CRI(aaList1) #[1,2,3,4,5,6,7,cri]
                                    if mutType[0] == 0:
                                        con1S += 1
                                        sum2PQ_C1 += twoPQ1
                                        if singleton == True:
                                            con1US += 1
                                    else:
                                        rad1S += 1
                                        sum2PQ_R1 += twoPQ1
                                        if singleton == True:
                                            rad1US += 1
                                    if mutType[1] == 0:
                                        con2S += 1
                                        sum2PQ_C2 += twoPQ1
                                        if singleton == True:
                                            con2US += 1
                                    else:
                                        rad2S += 1
                                        sum2PQ_R2 += twoPQ1
                                        if singleton == True:
                                            rad2US += 1
                                    if mutType[2] == 0:
                                        con3S += 1
                                        sum2PQ_C3 += twoPQ1
                                        if singleton == True:
                                            con3US += 1
                                    else:
                                        rad3S += 1
                                        sum2PQ_R3 += twoPQ1
                                        if singleton == True:
                                            rad3US += 1
                                    if mutType[3] == 0:
                                        con4S += 1
                                        sum2PQ_C4 += twoPQ1
                                        if singleton == True:
                                            con4US += 1
                                    else:
                                        rad4S += 1
                                        sum2PQ_R4 += twoPQ1
                                        if singleton == True:
                                            rad4US += 1
                                    if mutType[4] == 0:
                                        con5S += 1
                                        sum2PQ_C5 += twoPQ1
                                        if singleton == True:
                                            con5US += 1
                                    else:
                                        rad5S += 1
                                        sum2PQ_R5 += twoPQ1
                                        if singleton == True:
                                            rad5US += 1
                                    if mutType[5] == 0:
                                        con6S += 1
                                        sum2PQ_C6 += twoPQ1
                                        if singleton == True:
                                            con6US += 1
                                    else:
                                        rad6S += 1
                                        sum2PQ_R6 += twoPQ1
                                        if singleton == True:
                                            rad6US += 1
                                    if mutType[6] == 0:
                                        con7S += 1
                                        sum2PQ_C7 += twoPQ1
                                        if singleton == True:
                                            con7US += 1
                                    else:
                                        rad7S += 1
                                        sum2PQ_R7 += twoPQ1
                                        if singleton == True:
                                            rad7US += 1
                                    if mutType[7] <= 0.5:
                                        meanConS += 1
                                        sum2PQ_meanC += twoPQ1
                                        if singleton == True:
                                            meanConUS += 1
                                    else:
                                        meanRadS += 1
                                        sum2PQ_meanR += twoPQ1
                                        if singleton == True:
                                            meanRadUS += 1
                                else:
                                    nsynS += 2
                                    if singleton == True:
                                        nsynUS += 2
                                    twoPQ1 = 2*currFreqDict[codonList1[1]]*(1-currFreqDict[codonList1[1]]) #nsyn 
                                    twoPQ2 = 2*currFreqDict[codonList2[1]]*(1-currFreqDict[codonList2[1]]) #nsyn
                                    sum2PQ_N += twoPQ1 + twoPQ2
                                    mutType = CRI(aaList1) #[1,2,3,4,5,6,7,cri]
                                    if mutType[0] == 0:
                                        con1S += 1
                                        sum2PQ_C1 += twoPQ1
                                        if singleton == True:
                                            con1US += 1
                                    else:
                                        rad1S += 1
                                        sum2PQ_R1 += twoPQ1
                                        if singleton == True:
                                            rad1US += 1
                                    if mutType[1] == 0:
                                        con2S += 1
                                        sum2PQ_C2 += twoPQ1
                                        if singleton == True:
                                            con2US += 1
                                    else:
                                        rad2S += 1
                                        sum2PQ_R2 += twoPQ1
                                        if singleton == True:
                                            rad2US += 1
                                    if mutType[2] == 0:
                                        con3S += 1
                                        sum2PQ_C3 += twoPQ1
                                        if singleton == True:
                                            con3US += 1
                                    else:
                                        rad3S += 1
                                        sum2PQ_R3 += twoPQ1
                                        if singleton == True:
                                            rad3US += 1
                                    if mutType[3] == 0:
                                        con4S += 1
                                        sum2PQ_C4 += twoPQ1
                                        if singleton == True:
                                            con4US += 1
                                    else:
                                        rad4S += 1
                                        sum2PQ_R4 += twoPQ1
                                        if singleton == True:
                                            rad4US += 1
                                    if mutType[4] == 0:
                                        con5S += 1
                                        sum2PQ_C5 += twoPQ1
                                        if singleton == True:
                                            con5US += 1
                                    else:
                                        rad5S += 1
                                        sum2PQ_R5 += twoPQ1
                                        if singleton == True:
                                            rad5US += 1
                                    if mutType[5] == 0:
                                        con6S += 1
                                        sum2PQ_C6 += twoPQ1
                                        if singleton == True:
                                            con6US += 1
                                    else:
                                        rad6S += 1
                                        sum2PQ_R6 += twoPQ1
                                        if singleton == True:
                                            rad6US += 1
                                    if mutType[6] == 0:
                                        con7S += 1
                                        sum2PQ_C7 += twoPQ1
                                        if singleton == True:
                                            con7US += 1
                                    else:
                                        rad7S += 1
                                        sum2PQ_R7 += twoPQ1
                                        if singleton == True:
                                            rad7US += 1
                                    if mutType[7] <= 0.5:
                                        meanConS += 1
                                        sum2PQ_meanC += twoPQ1
                                        if singleton == True:
                                            meanConUS += 1
                                    else:
                                        meanRadS += 1
                                        sum2PQ_meanR += twoPQ1
                                        if singleton == True:
                                            meanRadUS += 1
                                    mutType = CRI(aaList2) #[1,2,3,4,5,6,7,cri]
                                    if mutType[0] == 0:
                                        con1S += 1
                                        sum2PQ_C1 += twoPQ2
                                        if singleton == True:
                                            con1US += 1
                                    else:
                                        rad1S += 1
                                        sum2PQ_R1 += twoPQ2
                                        if singleton == True:
                                            rad1US += 1
                                    if mutType[1] == 0:
                                        con2S += 1
                                        sum2PQ_C2 += twoPQ2
                                        if singleton == True:
                                            con2US += 1
                                    else:
                                        rad2S += 1
                                        sum2PQ_R2 += twoPQ2
                                        if singleton == True:
                                            rad2US += 1
                                    if mutType[2] == 0:
                                        con3S += 1
                                        sum2PQ_C3 += twoPQ2
                                        if singleton == True:
                                            con3US += 1
                                    else:
                                        rad3S += 1
                                        sum2PQ_R3 += twoPQ2
                                        if singleton == True:
                                            rad3US += 1
                                    if mutType[3] == 0:
                                        con4S += 1
                                        sum2PQ_C4 += twoPQ2
                                        if singleton == True:
                                            con4US += 1
                                    else:
                                        rad4S += 1
                                        sum2PQ_R4 += twoPQ2
                                        if singleton == True:
                                            rad4US += 1
                                    if mutType[4] == 0:
                                        con5S += 1
                                        sum2PQ_C5 += twoPQ2
                                        if singleton == True:
                                            con5US += 1
                                    else:
                                        rad5S += 1
                                        sum2PQ_R5 += twoPQ2
                                        if singleton == True:
                                            rad5US += 1
                                    if mutType[5] == 0:
                                        con6S += 1
                                        sum2PQ_C6 += twoPQ2
                                        if singleton == True:
                                            con6US += 1
                                    else:
                                        rad6S += 1
                                        sum2PQ_R6 += twoPQ2
                                        if singleton == True:
                                            rad6US += 1
                                    if mutType[6] == 0:
                                        con7S += 1
                                        sum2PQ_C7 += twoPQ2
                                        if singleton == True:
                                            con7US += 1
                                    else:
                                        rad7S += 1
                                        sum2PQ_R7 += twoPQ2
                                        if singleton == True:
                                            rad7US += 1
                                    if mutType[7] <= 0.5:
                                        meanConS += 1
                                        sum2PQ_meanC += twoPQ2
                                        if singleton == True:
                                            meanConUS += 1
                                    else:
                                        meanRadS += 1
                                        sum2PQ_meanR += twoPQ2                                    
                                        if singleton == True:
                                            meanRadUS += 1
                            elif outCodon == codonC:
                                aaList1 = [currAADict[codonC],currAADict[codonA]]
                                aaList2 = [currAADict[codonC],currAADict[codonB]]
                                codonList1 = [codonA,codonB]
                                codonList2 = [codonA,codonC]
                                if aaList1[0] == aaList1[1]:
                                    if aaList2[0] == aaList2[1]:
                                        synS += 2
                                        if singleton == True:
                                            synUS += 2
                                        twoPQ = 4
                                        for allele in currFreqDict:
                                            twoPQ *= currFreqDict[allele]
                                        sum2PQ_S += twoPQ                                        
                                    else:
                                        twoPQ1 = 2*currFreqDict[codonList1[1]]*(1-currFreqDict[codonList1[1]]) #syn 
                                        twoPQ2 = 2*currFreqDict[codonList2[1]]*(1-currFreqDict[codonList2[1]]) #nsyn
                                        sum2PQ_S += twoPQ1
                                        sum2PQ_N += twoPQ2
                                        synS += 1
                                        nsynS += 1
                                        if singleton == True:
                                            synUS += 1
                                            nsynUS += 1
                                        mutType = CRI(aaList2) #[1,2,3,4,5,6,7,cri]
                                        if mutType[0] == 0:
                                            con1S += 1
                                            sum2PQ_C1 += twoPQ2
                                            if singleton == True:
                                                con1US += 1
                                        else:
                                            rad1S += 1
                                            sum2PQ_R1 += twoPQ2
                                            if singleton == True:
                                                rad1US += 1
                                        if mutType[1] == 0:
                                            con2S += 1
                                            sum2PQ_C2 += twoPQ2
                                            if singleton == True:
                                                con2US += 1
                                        else:
                                            rad2S += 1
                                            sum2PQ_R2 += twoPQ2
                                            if singleton == True:
                                                rad2US += 1
                                        if mutType[2] == 0:
                                            con3S += 1
                                            sum2PQ_C3 += twoPQ2
                                            if singleton == True:
                                                con3US += 1
                                        else:
                                            rad3S += 1
                                            sum2PQ_R3 += twoPQ2
                                            if singleton == True:
                                                rad3US += 1
                                        if mutType[3] == 0:
                                            con4S += 1
                                            sum2PQ_C4 += twoPQ2
                                            if singleton == True:
                                                con4US += 1
                                        else:
                                            rad4S += 1
                                            sum2PQ_R4 += twoPQ2
                                            if singleton == True:
                                                rad4US += 1
                                        if mutType[4] == 0:
                                            con5S += 1
                                            sum2PQ_C5 += twoPQ2
                                            if singleton == True:
                                                con5US += 1
                                        else:
                                            rad5S += 1
                                            sum2PQ_R5 += twoPQ2
                                            if singleton == True:
                                                rad5US += 1
                                        if mutType[5] == 0:
                                            con6S += 1
                                            sum2PQ_C6 += twoPQ2
                                            if singleton == True:
                                                con6US += 1
                                        else:
                                            rad6S += 1
                                            sum2PQ_R6 += twoPQ2
                                            if singleton == True:
                                                rad6US += 1
                                        if mutType[6] == 0:
                                            con7S += 1
                                            sum2PQ_C7 += twoPQ2
                                            if singleton == True:
                                                con7US += 1
                                        else:
                                            rad7S += 1
                                            sum2PQ_R7 += twoPQ2
                                            if singleton == True:
                                                rad7US += 1
                                        if mutType[7] <= 0.5:
                                            meanConS += 1
                                            sum2PQ_meanC += twoPQ2
                                            if singleton == True:
                                                meanConUS += 1
                                        else:
                                            meanRadS += 1
                                            sum2PQ_meanR += twoPQ2
                                            if singleton == True:
                                                meanRadUS += 1
                                elif aaList2[0] == aaList2[1]:
                                    nsynS += 1
                                    synS += 1
                                    if singleton == True:
                                        synUS += 1
                                        nsynUS += 1
                                    twoPQ1 = 2*currFreqDict[codonList1[1]]*(1-currFreqDict[codonList1[1]]) #nsyn 
                                    twoPQ2 = 2*currFreqDict[codonList2[1]]*(1-currFreqDict[codonList2[1]]) #syn
                                    sum2PQ_S += twoPQ2
                                    sum2PQ_N += twoPQ1
                                    mutType = CRI(aaList1) #[1,2,3,4,5,6,7,cri]
                                    if mutType[0] == 0:
                                        con1S += 1
                                        sum2PQ_C1 += twoPQ1
                                        if singleton == True:
                                            con1US += 1
                                    else:
                                        rad1S += 1
                                        sum2PQ_R1 += twoPQ1
                                        if singleton == True:
                                            rad1US += 1
                                    if mutType[1] == 0:
                                        con2S += 1
                                        sum2PQ_C2 += twoPQ1
                                        if singleton == True:
                                            con2US += 1
                                    else:
                                        rad2S += 1
                                        sum2PQ_R2 += twoPQ1
                                        if singleton == True:
                                            rad2US += 1
                                    if mutType[2] == 0:
                                        con3S += 1
                                        sum2PQ_C3 += twoPQ1
                                        if singleton == True:
                                            con3US += 1
                                    else:
                                        rad3S += 1
                                        sum2PQ_R3 += twoPQ1
                                        if singleton == True:
                                            rad3US += 1
                                    if mutType[3] == 0:
                                        con4S += 1
                                        sum2PQ_C4 += twoPQ1
                                        if singleton == True:
                                            con4US += 1
                                    else:
                                        rad4S += 1
                                        sum2PQ_R4 += twoPQ1
                                        if singleton == True:
                                            rad4US += 1
                                    if mutType[4] == 0:
                                        con5S += 1
                                        sum2PQ_C5 += twoPQ1
                                        if singleton == True:
                                            con5US += 1
                                    else:
                                        rad5S += 1
                                        sum2PQ_R5 += twoPQ1
                                        if singleton == True:
                                            rad5US += 1
                                    if mutType[5] == 0:
                                        con6S += 1
                                        sum2PQ_C6 += twoPQ1
                                        if singleton == True:
                                            con6US += 1
                                    else:
                                        rad6S += 1
                                        sum2PQ_R6 += twoPQ1
                                        if singleton == True:
                                            rad6US += 1
                                    if mutType[6] == 0:
                                        con7S += 1
                                        sum2PQ_C7 += twoPQ1
                                        if singleton == True:
                                            con7US += 1
                                    else:
                                        rad7S += 1
                                        sum2PQ_R7 += twoPQ1
                                        if singleton == True:
                                            rad7US += 1
                                    if mutType[7] <= 0.5:
                                        meanConS += 1
                                        sum2PQ_meanC += twoPQ1
                                        if singleton == True:
                                            meanConUS += 1
                                    else:
                                        meanRadS += 1
                                        sum2PQ_meanR += twoPQ1
                                        if singleton == True:
                                            meanRadUS += 1
                                else:
                                    nsynS += 2
                                    if singleton == True:
                                        nsynUS += 2
                                    twoPQ1 = 2*currFreqDict[codonList1[1]]*(1-currFreqDict[codonList1[1]]) #nsyn 
                                    twoPQ2 = 2*currFreqDict[codonList2[1]]*(1-currFreqDict[codonList2[1]]) #nsyn
                                    sum2PQ_N += twoPQ1 + twoPQ2
                                    mutType = CRI(aaList1) #[1,2,3,4,5,6,7,cri]
                                    if mutType[0] == 0:
                                        con1S += 1
                                        sum2PQ_C1 += twoPQ1
                                        if singleton == True:
                                            con1US += 1
                                    else:
                                        rad1S += 1
                                        sum2PQ_R1 += twoPQ1
                                        if singleton == True:
                                            rad1US += 1
                                    if mutType[1] == 0:
                                        con2S += 1
                                        sum2PQ_C2 += twoPQ1
                                        if singleton == True:
                                            con2US += 1
                                    else:
                                        rad2S += 1
                                        sum2PQ_R2 += twoPQ1
                                        if singleton == True:
                                            rad2US += 1
                                    if mutType[2] == 0:
                                        con3S += 1
                                        sum2PQ_C3 += twoPQ1
                                        if singleton == True:
                                            con3US += 1
                                    else:
                                        rad3S += 1
                                        sum2PQ_R3 += twoPQ1
                                        if singleton == True:
                                            rad3US += 1
                                    if mutType[3] == 0:
                                        con4S += 1
                                        sum2PQ_C4 += twoPQ1
                                        if singleton == True:
                                            con4US += 1
                                    else:
                                        rad4S += 1
                                        sum2PQ_R4 += twoPQ1
                                        if singleton == True:
                                            rad4US += 1
                                    if mutType[4] == 0:
                                        con5S += 1
                                        sum2PQ_C5 += twoPQ1
                                        if singleton == True:
                                            con5US += 1
                                    else:
                                        rad5S += 1
                                        sum2PQ_R5 += twoPQ1
                                        if singleton == True:
                                            rad5US += 1
                                    if mutType[5] == 0:
                                        con6S += 1
                                        sum2PQ_C6 += twoPQ1
                                        if singleton == True:
                                            con6US += 1
                                    else:
                                        rad6S += 1
                                        sum2PQ_R6 += twoPQ1
                                        if singleton == True:
                                            rad6US += 1
                                    if mutType[6] == 0:
                                        con7S += 1
                                        sum2PQ_C7 += twoPQ1
                                        if singleton == True:
                                            con7US += 1
                                    else:
                                        rad7S += 1
                                        sum2PQ_R7 += twoPQ1
                                        if singleton == True:
                                            rad7US += 1
                                    if mutType[7] <= 0.5:
                                        meanConS += 1
                                        sum2PQ_meanC += twoPQ1
                                        if singleton == True:
                                            meanConUS += 1
                                    else:
                                        meanRadS += 1
                                        sum2PQ_meanR += twoPQ1
                                        if singleton == True:
                                            meanRadUS += 1
                                    mutType = CRI(aaList2) #[1,2,3,4,5,6,7,cri]
                                    if mutType[0] == 0:
                                        con1S += 1
                                        sum2PQ_C1 += twoPQ2
                                        if singleton == True:
                                            con1US += 1
                                    else:
                                        rad1S += 1
                                        sum2PQ_R1 += twoPQ2
                                        if singleton == True:
                                            rad1US += 1
                                    if mutType[1] == 0:
                                        con2S += 1
                                        sum2PQ_C2 += twoPQ2
                                        if singleton == True:
                                            con2US += 1
                                    else:
                                        rad2S += 1
                                        sum2PQ_R2 += twoPQ2
                                        if singleton == True:
                                            rad2US += 1
                                    if mutType[2] == 0:
                                        con3S += 1
                                        sum2PQ_C3 += twoPQ2
                                        if singleton == True:
                                            con3US += 1
                                    else:
                                        rad3S += 1
                                        sum2PQ_R3 += twoPQ2
                                        if singleton == True:
                                            rad3US += 1
                                    if mutType[3] == 0:
                                        con4S += 1
                                        sum2PQ_C4 += twoPQ2
                                        if singleton == True:
                                            con4US += 1
                                    else:
                                        rad4S += 1
                                        sum2PQ_R4 += twoPQ2
                                        if singleton == True:
                                            rad4US += 1
                                    if mutType[4] == 0:
                                        con5S += 1
                                        sum2PQ_C5 += twoPQ2
                                        if singleton == True:
                                            con5US += 1
                                    else:
                                        rad5S += 1
                                        sum2PQ_R5 += twoPQ2
                                        if singleton == True:
                                            rad5US += 1
                                    if mutType[5] == 0:
                                        con6S += 1
                                        sum2PQ_C6 += twoPQ2
                                        if singleton == True:
                                            con6US += 1
                                    else:
                                        rad6S += 1
                                        sum2PQ_R6 += twoPQ2
                                        if singleton == True:
                                            rad6US += 1
                                    if mutType[6] == 0:
                                        con7S += 1
                                        sum2PQ_C7 += twoPQ2
                                        if singleton == True:
                                            con7US += 1
                                    else:
                                        rad7S += 1
                                        sum2PQ_R7 += twoPQ2
                                        if singleton == True:
                                            rad7US += 1
                                    if mutType[7] <= 0.5:
                                        meanConS += 1
                                        sum2PQ_meanC += twoPQ2
                                        if singleton == True:
                                            meanConUS += 1
                                    else:
                                        meanRadS += 1
                                        sum2PQ_meanR += twoPQ2                                    
                                        if singleton == True:
                                            meanRadUS += 1
                            else:
                                if len(aaList) > 1:
                                    mutType = CRI(aaList) #[1,2,3,4,5,6,7,cri]
                                else:
                                    mutType = ''                                
                        else:
                            if len(aaList) > 1:
                                mutType = CRI(aaList) #[1,2,3,4,5,6,7,cri]
                            else:
                                mutType = ''                             
                    else:
                        if ab > ac and ab > bc:
                            codonList1 = [codonC,codonB]
                            codonList2 = [codonC,codonA]
                        elif ac > ab and ac > bc:
                            codonList1 = [codonB,codonA]
                            codonList2 = [codonB,codonC]
                        elif bc > ab and bc > ac:
                            codonList1 = [codonA,codonB]
                            codonList2 = [codonA,codonC]
                        aaList1 = []
                        aaList2 = []
                        for comp in codonList1:
                            if i < 3:
                                if comp in startCodons:
                                    aaList1.append('M')
                                else:
                                    aaList1.append(geneticCode[comp])
                            else:
                                aaList1.append(geneticCode[comp])
                        for comp in codonList2:
                            if i < 3:
                                if comp in startCodons:
                                    aaList2.append('M')
                                else:
                                    aaList2.append(geneticCode[comp])
                            else:
                                aaList2.append(geneticCode[comp])
                        if aaList1[0] == aaList1[1]:
                            if aaList2[0] == aaList2[1]:
                                synS += 2
                                twoPQ = 4
                                if singleton == True:
                                    synUS += 2
                                for allele in currFreqDict:
                                    twoPQ *= currFreqDict[allele]
                                sum2PQ_S += twoPQ                                
                            else:
                                twoPQ1 = 2*currFreqDict[codonList1[1]]*(1-currFreqDict[codonList1[1]]) #syn 
                                twoPQ2 = 2*currFreqDict[codonList2[1]]*(1-currFreqDict[codonList2[1]]) #nsyn
                                sum2PQ_S += twoPQ1
                                sum2PQ_N += twoPQ2
                                synS += 1
                                nsynS += 1
                                if singleton == True:
                                    synUS += 1
                                    nsynUS += 1
                                mutType = CRI(aaList2) #[1,2,3,4,5,6,7,cri]
                                if mutType[0] == 0:
                                    con1S += 1
                                    sum2PQ_C1 += twoPQ2
                                    if singleton == True:
                                        con1US += 1
                                else:
                                    rad1S += 1
                                    sum2PQ_R1 += twoPQ2
                                    if singleton == True:
                                        rad1US += 1
                                if mutType[1] == 0:
                                    con2S += 1
                                    sum2PQ_C2 += twoPQ2
                                    if singleton == True:
                                        con2US += 1
                                else:
                                    rad2S += 1
                                    sum2PQ_R2 += twoPQ2
                                    if singleton == True:
                                        rad2US += 1
                                if mutType[2] == 0:
                                    con3S += 1
                                    sum2PQ_C3 += twoPQ2
                                    if singleton == True:
                                        con3US += 1
                                else:
                                    rad3S += 1
                                    sum2PQ_R3 += twoPQ2
                                    if singleton == True:
                                        rad3US += 1
                                if mutType[3] == 0:
                                    con4S += 1
                                    sum2PQ_C4 += twoPQ2
                                    if singleton == True:
                                        con4US += 1
                                else:
                                    rad4S += 1
                                    sum2PQ_R4 += twoPQ2
                                    if singleton == True:
                                        rad4US += 1
                                if mutType[4] == 0:
                                    con5S += 1
                                    sum2PQ_C5 += twoPQ2
                                    if singleton == True:
                                        con5US += 1
                                else:
                                    rad5S += 1
                                    sum2PQ_R5 += twoPQ2
                                    if singleton == True:
                                        rad5US += 1
                                if mutType[5] == 0:
                                    con6S += 1
                                    sum2PQ_C6 += twoPQ2
                                    if singleton == True:
                                        con6US += 1
                                else:
                                    rad6S += 1
                                    sum2PQ_R6 += twoPQ2
                                    if singleton == True:
                                        rad6US += 1
                                if mutType[6] == 0:
                                    con7S += 1
                                    sum2PQ_C7 += twoPQ2
                                    if singleton == True:
                                        con7US += 1
                                else:
                                    rad7S += 1
                                    sum2PQ_R7 += twoPQ2
                                    if singleton == True:
                                        rad7US += 1
                                if mutType[7] <= 0.5:
                                    meanConS += 1
                                    sum2PQ_meanC += twoPQ2
                                    if singleton == True:
                                        meanConUS += 1
                                else:
                                    meanRadS += 1
                                    sum2PQ_meanR += twoPQ2                                
                                    if singleton == True:
                                        meanRadUS += 1
                        elif aaList2[0] == aaList2[1]:
                            nsynS += 1
                            synS += 1
                            if singleton == True:
                                synUS += 1
                                nsynUS += 1
                            twoPQ1 = 2*currFreqDict[codonList1[1]]*(1-currFreqDict[codonList1[1]]) #nsyn 
                            twoPQ2 = 2*currFreqDict[codonList2[1]]*(1-currFreqDict[codonList2[1]]) #syn
                            sum2PQ_S += twoPQ2
                            sum2PQ_N += twoPQ1
                            mutType = CRI(aaList1) #[1,2,3,4,5,6,7,cri]
                            if mutType[0] == 0:
                                con1S += 1
                                sum2PQ_C1 += twoPQ1
                                if singleton == True:
                                    con1US += 1
                            else:
                                rad1S += 1
                                sum2PQ_R1 += twoPQ1
                                if singleton == True:
                                    rad1US += 1
                            if mutType[1] == 0:
                                con2S += 1
                                sum2PQ_C2 += twoPQ1
                                if singleton == True:
                                    con2US += 1
                            else:
                                rad2S += 1
                                sum2PQ_R2 += twoPQ1
                                if singleton == True:
                                    rad2US += 1
                            if mutType[2] == 0:
                                con3S += 1
                                sum2PQ_C3 += twoPQ1
                                if singleton == True:
                                    con3US += 1
                            else:
                                rad3S += 1
                                sum2PQ_R3 += twoPQ1
                                if singleton == True:
                                    rad3US += 1
                            if mutType[3] == 0:
                                con4S += 1
                                sum2PQ_C4 += twoPQ1
                                if singleton == True:
                                    con4US += 1
                            else:
                                rad4S += 1
                                sum2PQ_R4 += twoPQ1
                                if singleton == True:
                                    rad4US += 1
                            if mutType[4] == 0:
                                con5S += 1
                                sum2PQ_C5 += twoPQ1
                                if singleton == True:
                                    con5US += 1
                            else:
                                rad5S += 1
                                sum2PQ_R5 += twoPQ1
                                if singleton == True:
                                    rad5US += 1
                            if mutType[5] == 0:
                                con6S += 1
                                sum2PQ_C6 += twoPQ1
                                if singleton == True:
                                    con6US += 1
                            else:
                                rad6S += 1
                                sum2PQ_R6 += twoPQ1
                                if singleton == True:
                                    rad6US += 1
                            if mutType[6] == 0:
                                con7S += 1
                                sum2PQ_C7 += twoPQ1
                                if singleton == True:
                                    con7US += 1
                            else:
                                rad7S += 1
                                sum2PQ_R7 += twoPQ1
                                if singleton == True:
                                    rad7US += 1
                            if mutType[7] <= 0.5:
                                meanConS += 1
                                sum2PQ_meanC += twoPQ1
                                if singleton == True:
                                    meanConUS += 1
                            else:
                                meanRadS += 1
                                sum2PQ_meanR += twoPQ1
                                if singleton == True:
                                    meanRadUS += 1
                        else:
                            nsynS += 2
                            if singleton == True:
                                nsynUS += 2
                            twoPQ1 = 2*currFreqDict[codonList1[1]]*(1-currFreqDict[codonList1[1]]) #nsyn 
                            twoPQ2 = 2*currFreqDict[codonList2[1]]*(1-currFreqDict[codonList2[1]]) #nsyn
                            sum2PQ_N += twoPQ1 + twoPQ2
                            mutType = CRI(aaList1) #[1,2,3,4,5,6,7,cri]
                            if mutType[0] == 0:
                                con1S += 1
                                sum2PQ_C1 += twoPQ1
                                if singleton == True:
                                    con1US += 1
                            else:
                                rad1S += 1
                                sum2PQ_R1 += twoPQ1
                                if singleton == True:
                                    rad1US += 1
                            if mutType[1] == 0:
                                con2S += 1
                                sum2PQ_C2 += twoPQ1
                                if singleton == True:
                                    con2US += 1
                            else:
                                rad2S += 1
                                sum2PQ_R2 += twoPQ1
                                if singleton == True:
                                    rad2US += 1
                            if mutType[2] == 0:
                                con3S += 1
                                sum2PQ_C3 += twoPQ1
                                if singleton == True:
                                    con3US += 1
                            else:
                                rad3S += 1
                                sum2PQ_R3 += twoPQ1
                                if singleton == True:
                                    rad3US += 1
                            if mutType[3] == 0:
                                con4S += 1
                                sum2PQ_C4 += twoPQ1
                                if singleton == True:
                                    con4US += 1
                            else:
                                rad4S += 1
                                sum2PQ_R4 += twoPQ1
                                if singleton == True:
                                    rad4US += 1
                            if mutType[4] == 0:
                                con5S += 1
                                sum2PQ_C5 += twoPQ1
                                if singleton == True:
                                    con5US += 1
                            else:
                                rad5S += 1
                                sum2PQ_R5 += twoPQ1
                                if singleton == True:
                                    rad5US += 1
                            if mutType[5] == 0:
                                con6S += 1
                                sum2PQ_C6 += twoPQ1
                                if singleton == True:
                                    con6US += 1
                            else:
                                rad6S += 1
                                sum2PQ_R6 += twoPQ1
                                if singleton == True:
                                    rad6US += 1
                            if mutType[6] == 0:
                                con7S += 1
                                sum2PQ_C7 += twoPQ1
                                if singleton == True:
                                    con7US += 1
                            else:
                                rad7S += 1
                                sum2PQ_R7 += twoPQ1
                                if singleton == True:
                                    rad7US += 1
                            if mutType[7] <= 0.5:
                                meanConS += 1
                                sum2PQ_meanC += twoPQ1
                                if singleton == True:
                                    meanConUS += 1
                            else:
                                meanRadS += 1
                                sum2PQ_meanR += twoPQ1
                                if singleton == True:
                                    meanRadUS += 1
                            mutType = CRI(aaList2) #[1,2,3,4,5,6,7,cri]
                            if mutType[0] == 0:
                                con1S += 1
                                sum2PQ_C1 += twoPQ2
                                if singleton == True:
                                    con1US += 1
                            else:
                                rad1S += 1
                                sum2PQ_R1 += twoPQ2
                                if singleton == True:
                                    rad1US += 1
                            if mutType[1] == 0:
                                con2S += 1
                                sum2PQ_C2 += twoPQ2
                                if singleton == True:
                                    con2US += 1
                            else:
                                rad2S += 1
                                sum2PQ_R2 += twoPQ2
                                if singleton == True:
                                    rad2US += 1
                            if mutType[2] == 0:
                                con3S += 1
                                sum2PQ_C3 += twoPQ2
                                if singleton == True:
                                    con3US += 1
                            else:
                                rad3S += 1
                                sum2PQ_R3 += twoPQ2
                                if singleton == True:
                                    rad3US += 1
                            if mutType[3] == 0:
                                con4S += 1
                                sum2PQ_C4 += twoPQ2
                                if singleton == True:
                                    con4US += 1
                            else:
                                rad4S += 1
                                sum2PQ_R4 += twoPQ2
                                if singleton == True:
                                    rad4US += 1
                            if mutType[4] == 0:
                                con5S += 1
                                sum2PQ_C5 += twoPQ2
                                if singleton == True:
                                    con5US += 1
                            else:
                                rad5S += 1
                                sum2PQ_R5 += twoPQ2
                                if singleton == True:
                                    rad5US += 1
                            if mutType[5] == 0:
                                con6S += 1
                                sum2PQ_C6 += twoPQ2
                                if singleton == True:
                                    con6US += 1
                            else:
                                rad6S += 1
                                sum2PQ_R6 += twoPQ2
                                if singleton == True:
                                    rad6US += 1
                            if mutType[6] == 0:
                                con7S += 1
                                sum2PQ_C7 += twoPQ2
                                if singleton == True:
                                    con7US += 1
                            else:
                                rad7S += 1
                                sum2PQ_R7 += twoPQ2
                                if singleton == True:
                                    rad7US += 1
                            if mutType[7] <= 0.5:
                                meanConS += 1
                                sum2PQ_meanC += twoPQ2
                                if singleton == True:
                                    meanConUS += 1
                            else:
                                meanRadS += 1
                                sum2PQ_meanR += twoPQ2
                                if singleton == True:
                                    meanRadUS += 1
                elif len(currAlleleDict) == 2:
                    currFreqDict = {}
                    twoPQ = 2
                    for codon in currAlleleDict:
                        freq = float(currAlleleDict[codon])/totalIndividuals
                        twoPQ *= freq
                        currFreqDict[codon] = freq
                    if len(aaList) == 1:
                        synS += 2
                        if singleton == True:
                            synUS += 1
                        sum2PQ_S += (2*twoPQ)                        
        i += 1
    thetaS = synS/(synL*An)
    thetaA = nsynS/(nsynL*An)
    thetaC1 = con1S/(con1L*An) 
    thetaR1 = rad1S/(rad1L*An)
    thetaC2 = con2S/(con2L*An)
    thetaR2 = rad2S/(rad2L*An)
    thetaC3 = con3S/(con3L*An)
    thetaR3 = rad3S/(rad3L*An)
    thetaC4 = con4S/(con4L*An)
    thetaR4 = rad4S/(rad4L*An)
    thetaC5 = con5S/(con5L*An)
    thetaR5 = rad5S/(rad5L*An)
    thetaC6 = con6S/(con6L*An)
    thetaR6 = rad6S/(rad6L*An)
    thetaC7 = con7S/(con7L*An)
    thetaR7 = rad7S/(rad7L*An)
    thetaMeanC = meanConS/(meanConL*An)
    thetaMeanR = meanRadS/(meanRadL*An)
    piS = (31.0/30)*sum2PQ_S/synL
    piA = (31.0/30)*sum2PQ_N/nsynL
    piC1 = (31.0/30)*sum2PQ_C1/con1L
    piR1 = (31.0/30)*sum2PQ_R1/rad1L
    piC2 = (31.0/30)*sum2PQ_C2/con2L
    piR2 = (31.0/30)*sum2PQ_R2/rad2L
    piC3 = (31.0/30)*sum2PQ_C3/con3L
    piR3 = (31.0/30)*sum2PQ_R3/rad3L
    piC4 = (31.0/30)*sum2PQ_C4/con4L
    piR4 = (31.0/30)*sum2PQ_R4/rad4L
    piC5 = (31.0/30)*sum2PQ_C5/con5L
    piR5 = (31.0/30)*sum2PQ_R5/rad5L
    piC6 = (31.0/30)*sum2PQ_C6/con6L
    piR6 = (31.0/30)*sum2PQ_R6/rad6L
    piC7 = (31.0/30)*sum2PQ_C7/con7L
    piR7 = (31.0/30)*sum2PQ_R7/rad7L
    piMeanC = (31.0/30)*sum2PQ_meanC/meanConL
    piMeanR = (31.0/30)*sum2PQ_meanR/meanRadL
    sys.stdout.write(str(synS) + '\t' + str(nsynS) + '\t' + str(con1S) + '\t' + str(rad1S) + '\t' + str(synL) + '\t' + str(nsynL) + '\t' + str(con1L) + '\t' + str(rad1L) + '\n')
    return thetaS,thetaA,thetaC1,thetaR1,thetaC2,thetaR2,thetaC3,thetaR3,thetaC4,thetaR4,thetaC5,thetaR5,thetaC6,thetaR6,thetaC7,thetaR7,thetaMeanC,thetaMeanR,piS,piA,piC1,piR1,piC2,piR2,piC3,piR3,piC4,piR4,piC5,piR5,piC6,piR6,piC7,piR7,piMeanC,piMeanR
    

def thetaU(seqDict,seqList):
    geneticCodes = {'standard':{"TTT":"F",	"TTC":"F",	"TTA":"L",	"TTG":"L",	"TCT":"S",	"TCC":"S",	"TCA":"S",	"TCG":"S",	"TAT":"Y",	"TAC":"Y",	"TAA":"*",	"TAG":"*",	"TGT":"C",	"TGC":"C",	"TGA":"*",	"TGG":"W",	"CTT":"L",	"CTC":"L",	"CTA":"L",	"CTG":"L",	"CCT":"P",	"CCC":"P",	"CCA":"P",	"CCG":"P",	"CAT":"H",	"CAC":"H",	"CAA":"Q",	"CAG":"Q",	"CGT":"R",	"CGC":"R",	"CGA":"R",	"CGG":"R",	"ATT":"I",	"ATC":"I",	"ATA":"I",	"ATG":"M",	"ACT":"T",	"ACC":"T",	"ACA":"T",	"ACG":"T",	"AAT":"N",	"AAC":"N",	"AAA":"K",	"AAG":"K",	"AGT":"S",	"AGC":"S",	"AGA":"R",	"AGG":"R",	"GTT":"V",	"GTC":"V",	"GTA":"V",	"GTG":"V",	"GCT":"A",	"GCC":"A",	"GCA":"A",	"GCG":"A",	"GAT":"D",	"GAC":"D",	"GAA":"E",	"GAG":"E",	"GGT":"G",	"GGC":"G",	"GGA":"G",	"GGG":"G"},'invertebrateMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'M', 'AGG': 'S', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'S', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'vertebrateMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'M', 'AGG': '*', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': '*', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'}, 'yeastMt':{'CTT': 'T', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'M', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'T', 'CTA': 'T', 'CTC': 'T', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'}, 'coelenterateMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'ciliateNuc':{'CTT': 'L', 'TAG': 'Q', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': '*', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': 'Q', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'echinodermMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'S', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'S', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'N', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'}, 'euplotidNuc':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'C', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'}, 'bacterial':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': '*', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'yeastNuc':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'S', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': '*', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'}, 'ascidianMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'M', 'AGG': 'G', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'G', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'flatwormMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'S', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'S', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'N', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': 'Y', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'chlorophyceanMt':{'CTT': 'L', 'TAG': 'L', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': '*', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'trematodeMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'M', 'AGG': 'S', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'S', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'N', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'pterobranchiaMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'K', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'S', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'}}
    geneticCode = geneticCodes['invertebrateMt']
    startCodons = ['ATT','ATC','ATA','ATG','GTG'] #invertebrateMt code
    codonDict,AADict = buildCodonDict(seqDict,seqList)
    popList = []
    outList = [] 
    synL = 0.0
    nsynL = 0.0
    con1L = 0.0 
    rad1L = 0.0
    con2L = 0.0
    rad2L = 0.0
    con3L = 0.0
    rad3L = 0.0
    con4L = 0.0
    rad4L = 0.0
    con5L = 0.0
    rad5L = 0.0
    con6L = 0.0
    rad6L = 0.0
    con7L = 0.0
    rad7L = 0.0
    meanConL = 0.0
    meanRadL = 0.0
    for seq in seqList:
        currCodonList = codonDict[seq]
        totalSynSites,totalNonsynSites,totalC1Sites,totalR1Sites,totalC2Sites,totalR2Sites,totalC3Sites,totalR3Sites,totalC4Sites,totalR4Sites,totalC5Sites,totalR5Sites,totalC6Sites,totalR6Sites,totalC7Sites,totalR7Sites,totalMeanCSites,totalMeanRSites = countSites(currCodonList)
        if seq == '>Potamopyrgus_estuarinus':
            outList.append(seq)
        else:
            popList.append(seq)
            synL += totalSynSites
            nsynL += totalNonsynSites
            con1L += totalC1Sites
            rad1L += totalR1Sites
            con2L += totalC2Sites
            rad2L += totalR2Sites
            con3L += totalC3Sites
            rad3L += totalR3Sites
            con4L += totalC4Sites
            rad4L += totalR4Sites
            con5L += totalC5Sites
            rad5L += totalR5Sites
            con6L += totalC6Sites
            rad6L += totalR6Sites
            con7L += totalC7Sites
            rad7L += totalR7Sites
            meanConL += totalMeanCSites
            meanRadL += totalMeanRSites          
    N = len(popList)
    An = aN(N)
    synL /= N
    nsynL /= N
    con1L /= N 
    rad1L /= N
    con2L /= N
    rad2L /= N
    con3L /= N
    rad3L /= N
    con4L /= N
    rad4L /= N
    con5L /= N
    rad5L /= N
    con6L /= N
    rad6L /= N
    con7L /= N
    rad7L /= N
    meanConL /= N
    meanRadL /= N
    i = 0
    synUS = 0.0
    nsynUS = 0.0
    con1US = 0.0
    rad1US = 0.0
    con2US = 0.0
    rad2US = 0.0
    con3US = 0.0
    rad3US = 0.0
    con4US = 0.0
    rad4US = 0.0
    con5US = 0.0
    rad5US = 0.0
    con6US = 0.0
    rad6US = 0.0
    con7US = 0.0
    rad7US = 0.0
    meanConUS = 0.0
    meanRadUS = 0.0
    outCodons = codonDict[outList[0]]
    outSeq = seqDict[outList[0]]
    polymorphicCodons = []
    while i < len(seqDict[seqList[0]]):
        currAlleleDict = {}
        currAlleleList = []
        for seq in popList:
            currSeq = seqDict[seq]
            currNuc = currSeq[i]
            if currNuc not in currAlleleDict and currNuc != 'N' and '-' != currNuc:
                currAlleleDict[currNuc] = 1
                currAlleleList.append(currNuc)
            elif 'N' != currNuc and '-' != currNuc:
                currValue = currAlleleDict[currNuc]
                currValue += 1
                currAlleleDict[currNuc] = currValue
        singleton = False
        alleleNum = 0
        if len(currAlleleList) > 1:
            while alleleNum < len(currAlleleList) and singleton == False:
                if currAlleleDict[currAlleleList[alleleNum]] == 1:
                    singleton = True
                alleleNum += 1
        if singleton == True:
            polymorphicCodons.append(i/3)
        i += 1
    for codonNum in polymorphicCodons:
        currCodonDict = {}
        currCodonList = []
        outCodon = outCodons[codonNum]
        lineageDict = {}
        for seq in popList:
            currCodons = codonDict[seq]
            currCodon = currCodons[codonNum]
            if currCodon not in currCodonList and 'N' not in currCodon and '-' not in currCodon:
                currCodonDict[currCodon] = 1
                currCodonList.append(currCodon)
                lineageDict[currCodon] = [seq]
            elif 'N' not in currCodon and '-' not in currCodon:
                currValue = currCodonDict[currCodon]
                currValue += 1
                currCodonDict[currCodon] = currValue
                currList = lineageDict[currCodon]
                currList.append(seq)
                lineageDict[currCodon] = currList
        if len(currCodonList) == 2:
            aa1 = geneticCode[currCodonList[0]]
            aa2 = geneticCode[currCodonList[1]]
            if aa1 == aa2:
                synUS += 1
            else:
                nsynUS += 1
                cri = CRI([aa1,aa2])
                if cri[0] == 0:
                    con1US += 1
                else:
                    rad1US += 1
                if cri[1] == 0:
                    con2US += 1
                else:
                    rad2US += 1
                if cri[2] == 0:
                    con3US += 1
                else:
                    rad3US += 1
                if cri[3] == 0:
                    con4US += 1
                else:
                    rad4US += 1
                if cri[4] == 0:
                    con5US += 1
                else:
                    rad5US += 1
                if cri[5] == 0:
                    con6US += 1
                else:
                    rad6US += 1
                if cri[6] == 0:
                    con7US += 1
                else:
                    rad7US += 1
                if cri[7] <= 0.5:
                    meanConUS += 1
                else:
                    meanRadUS += 1
        else:
            minCodon = currCodonList[0]
            closestCodon = currCodonList[0]
            for codon in currCodonList:
                if currCodonDict[codon] == 1:
                    minCodon = codon
            for codon in currCodonList:
                numChanges = 0
                for nuc in codon:
                    for nuc2 in minCodon:
                        if nuc != nuc2:
                            numChanges += 1
                if numChanges == 1:
                    closestCodon = codon
            if minCodon != closestCodon:
                aa1 = geneticCode[closestCodon]
                aa2 = geneticCode[minCodon]
                if aa1 == aa2:
                    synUS += 1
                else:
                    nsynUS += 1
                    cri = CRI([aa1,aa2])
                    if cri[0] == 0:
                        con1US += 1
                    else:
                        rad1US += 1
                    if cri[1] == 0:
                        con2US += 1
                    else:
                        rad2US += 1
                    if cri[2] == 0:
                        con3US += 1
                    else:
                        rad3US += 1
                    if cri[3] == 0:
                        con4US += 1
                    else:
                        rad4US += 1
                    if cri[4] == 0:
                        con5US += 1
                    else:
                        rad5US += 1
                    if cri[5] == 0:
                        con6US += 1
                    else:
                        rad6US += 1
                    if cri[6] == 0:
                        con7US += 1
                    else:
                        rad7US += 1
                    if cri[7] <= 0.5:
                        meanConUS += 1
                    else:
                        meanRadUS += 1                   
    thetaUS = synUS/(An*synL)
    thetaUA = nsynUS/(An*nsynL)
    thetaUC1 = con1US/(An*con1L)
    thetaUC2 = con2US/(An*con2L)
    thetaUC3 = con3US/(An*con3L)
    thetaUC4 = con4US/(An*con4L)
    thetaUC5 = con5US/(An*con5L)
    thetaUC6 = con6US/(An*con6L)
    thetaUC7 = con7US/(An*con7L)
    thetaUR1 = rad1US/(An*rad1L)
    thetaUR2 = rad2US/(An*rad2L)
    thetaUR3 = rad3US/(An*rad3L)
    thetaUR4 = rad4US/(An*rad4L)
    thetaUR5 = rad5US/(An*rad5L)
    thetaUR6 = rad6US/(An*rad6L)
    thetaUR7 = rad7US/(An*rad7L)
    thetaUMeanC = meanConUS/(An*meanConL)
    thetaUMeanR = meanRadUS/(An*meanRadL)
    return thetaUS,thetaUA,thetaUC1,thetaUR1,thetaUC2,thetaUR2,thetaUC3,thetaUR3,thetaUC4,thetaUR4,thetaUC5,thetaUR5,thetaUC6,thetaUR6,thetaUC7,thetaUR7,thetaUMeanC,thetaUMeanR
                
                
        
        
    
        
def buildCodonDict(seqDict,seqList):
    geneticCodes = {'standard':{"TTT":"F",	"TTC":"F",	"TTA":"L",	"TTG":"L",	"TCT":"S",	"TCC":"S",	"TCA":"S",	"TCG":"S",	"TAT":"Y",	"TAC":"Y",	"TAA":"*",	"TAG":"*",	"TGT":"C",	"TGC":"C",	"TGA":"*",	"TGG":"W",	"CTT":"L",	"CTC":"L",	"CTA":"L",	"CTG":"L",	"CCT":"P",	"CCC":"P",	"CCA":"P",	"CCG":"P",	"CAT":"H",	"CAC":"H",	"CAA":"Q",	"CAG":"Q",	"CGT":"R",	"CGC":"R",	"CGA":"R",	"CGG":"R",	"ATT":"I",	"ATC":"I",	"ATA":"I",	"ATG":"M",	"ACT":"T",	"ACC":"T",	"ACA":"T",	"ACG":"T",	"AAT":"N",	"AAC":"N",	"AAA":"K",	"AAG":"K",	"AGT":"S",	"AGC":"S",	"AGA":"R",	"AGG":"R",	"GTT":"V",	"GTC":"V",	"GTA":"V",	"GTG":"V",	"GCT":"A",	"GCC":"A",	"GCA":"A",	"GCG":"A",	"GAT":"D",	"GAC":"D",	"GAA":"E",	"GAG":"E",	"GGT":"G",	"GGC":"G",	"GGA":"G",	"GGG":"G"},'invertebrateMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'M', 'AGG': 'S', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'S', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'vertebrateMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'M', 'AGG': '*', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': '*', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'}, 'yeastMt':{'CTT': 'T', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'M', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'T', 'CTA': 'T', 'CTC': 'T', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'}, 'coelenterateMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'ciliateNuc':{'CTT': 'L', 'TAG': 'Q', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': '*', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': 'Q', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'echinodermMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'S', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'S', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'N', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'}, 'euplotidNuc':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'C', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'}, 'bacterial':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': '*', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'yeastNuc':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'S', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': '*', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'}, 'ascidianMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'M', 'AGG': 'G', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'G', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'flatwormMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'S', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'S', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'N', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': 'Y', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'chlorophyceanMt':{'CTT': 'L', 'TAG': 'L', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': '*', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'trematodeMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'M', 'AGG': 'S', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'S', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'N', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'pterobranchiaMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'K', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'S', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'}}
    geneticCode = geneticCodes['invertebrateMt']
    startCodons = ['ATT','ATC','ATA','ATG','GTG','TTG']
    codonDict = {}
    AADict = {}
    for seq in seqList:
        nucleotideSeq = seqDict[seq]
        codonList = []
        i = 2
        while i < len(nucleotideSeq):
            currCodon = nucleotideSeq[i-2] + nucleotideSeq[i-1] + nucleotideSeq[i]
            codonList.append(currCodon)
            i += 3
        codonDict[seq] = codonList
        AAseq = ''
        codonNum = 1
        for codon in codonList:
            if codonNum == 1:
                if codon in startCodons:
                    aa = 'M'
                elif codon in geneticCode:
                    aa = geneticCode[codon]
                else:
                    aa = 'X'
            elif codon in geneticCode:
                aa = geneticCode[codon]
            else:
                aa = 'X'
            AAseq += aa
            codonNum += 1
        if AAseq[-1] == '*':
            AAseq = AAseq[0:-1]
        AADict[seq] = AAseq
    return codonDict,AADict
    

    
def reverseComplement(seq):
    seq_revc = ''
    for nuc in seq:
        if nuc == 'A':
            seq_revc = 'T' + seq_revc
        elif nuc == 'T':
            seq_revc = 'A' + seq_revc
        elif nuc == 'C':
            seq_revc = 'G' + seq_revc
        elif nuc == 'G':
            seq_revc = 'C' + seq_revc
        else:
            seq_revc = nuc + seq_revc
    return seq_revc
    
def seqDictGenerator(fasta):
    infile = open(fasta,'r')
    scaffoldDict = {}
    scaffoldList = []
    seqName = ''
    currSeq = ''
    for line in infile:
        if line[0] == '>':
            if seqName != '':
                scaffoldDict[seqName] = currSeq
            seqName = line
            while seqName[-1] == '\n' or seqName[-1] == '\t' or seqName[-1] == '\r':
                seqName = seqName[0:-1]
            scaffoldList.append(seqName)
            currSeq = ''
        else:
            currSeq += line
            while currSeq[-1] == '\n' or currSeq[-1] == '\t' or currSeq[-1] == '\r':
                currSeq = currSeq[0:-1]
    scaffoldDict[seqName] = currSeq 
    return scaffoldDict, scaffoldList

    
def bootstrapReplicate(codonDict,seqList):
    newSeqDict = {}
    seqLength = len(codonDict[seqList[0]])
    bootstrapNums = []
    while len(bootstrapNums) < seqLength:
        bootstrapNums.append(random.choice(range(seqLength)))
    for seq in seqList:
        oldSeq = codonDict[seq]
        newSeq = ''
        for site in bootstrapNums:
            newSeq += oldSeq[site]
        newSeqDict[seq] = newSeq
    return newSeqDict
     

def workingCode(fasta):
    infile = open(fasta,'r')
    seqDict = {}
    for line in infile:
        if line[0] == '>':
            seqName = line
            while seqName[-1] == '\n' or seqName[-1] == '\t' or seqName[-1] == '\r':
                seqName = seqName[0:-1]
            if seqName!= '>Potamopyrgus_estuarinus':
                seqDict[seqName] = 'antipodarum'
            else:
                seqDict[seqName] = 'estuarinus'
    return seqDict

def CRI(aaList):
    aaSchemeList = [1,2,3,4,5,6,7]
    aaSchemeDict = {1:{("R","H"):"C",("R","K"):"C",("R","D"):"R",("R","E"):"R",("R","A"):"R",("R","N"):"R",("R","C"):"R",("R","Q"):"R",("R","G"):"R",("R","I"):"R",("R","L"):"R",("R","M"):"R",("R","F"):"R",("R","P"):"R",("R","S"):"R",("R","T"):"R",("R","W"):"R",("R","Y"):"R",("R","V"):"R",("H","K"):"C",("H","D"):"R",("H","E"):"R",("H","A"):"R",("H","N"):"R",("H","C"):"R",("H","Q"):"R",("H","G"):"R",("H","I"):"R",("H","L"):"R",("H","M"):"R",("H","F"):"R",("H","P"):"R",("H","S"):"R",("H","T"):"R",("H","W"):"R",("H","Y"):"R",("H","V"):"R",("K","D"):"R",("K","E"):"R",("K","A"):"R",("K","N"):"R",("K","C"):"R",("K","Q"):"R",("K","G"):"R",("K","I"):"R",("K","L"):"R",("K","M"):"R",("K","F"):"R",("K","P"):"R",("K","S"):"R",("K","T"):"R",("K","W"):"R",("K","Y"):"R",("K","V"):"R",("D","E"):"C",("D","A"):"R",("D","N"):"R",("D","C"):"R",("D","Q"):"R",("D","G"):"R",("D","I"):"R",("D","L"):"R",("D","M"):"R",("D","F"):"R",("D","P"):"R",("D","S"):"R",("D","T"):"R",("D","W"):"R",("D","Y"):"R",("D","V"):"R",("E","A"):"R",("E","N"):"R",("E","C"):"R",("E","Q"):"R",("E","G"):"R",("E","I"):"R",("E","L"):"R",("E","M"):"R",("E","F"):"R",("E","P"):"R",("E","S"):"R",("E","T"):"R",("E","W"):"R",("E","Y"):"R",("E","V"):"R",("A","N"):"C",("A","C"):"C",("A","Q"):"C",("A","G"):"C",("A","I"):"C",("A","L"):"C",("A","M"):"C",("A","F"):"C",("A","P"):"C",("A","S"):"C",("A","T"):"C",("A","W"):"C",("A","Y"):"C",("A","V"):"C",("N","C"):"C",("N","Q"):"C",("N","G"):"C",("N","I"):"C",("N","L"):"C",("N","M"):"C",("N","F"):"C",("N","P"):"C",("N","S"):"C",("N","T"):"C",("N","W"):"C",("N","Y"):"C",("N","V"):"C",("C","Q"):"C",("C","G"):"C",("C","I"):"C",("C","L"):"C",("C","M"):"C",("C","F"):"C",("C","P"):"C",("C","S"):"C",("C","T"):"C",("C","W"):"C",("C","Y"):"C",("C","V"):"C",("Q","G"):"C",("Q","I"):"C",("Q","L"):"C",("Q","M"):"C",("Q","F"):"C",("Q","P"):"C",("Q","S"):"C",("Q","T"):"C",("Q","W"):"C",("Q","Y"):"C",("Q","V"):"C",("G","I"):"C",("G","L"):"C",("G","M"):"C",("G","F"):"C",("G","P"):"C",("G","S"):"C",("G","T"):"C",("G","W"):"C",("G","Y"):"C",("G","V"):"C",("I","L"):"C",("I","M"):"C",("I","F"):"C",("I","P"):"C",("I","S"):"C",("I","T"):"C",("I","W"):"C",("I","Y"):"C",("I","V"):"C",("L","M"):"C",("L","F"):"C",("L","P"):"C",("L","S"):"C",("L","T"):"C",("L","W"):"C",("L","Y"):"C",("L","V"):"C",("M","F"):"C",("M","P"):"C",("M","S"):"C",("M","T"):"C",("M","W"):"C",("M","Y"):"C",("M","V"):"C",("F","P"):"C",("F","S"):"C",("F","T"):"C",("F","W"):"C",("F","Y"):"C",("F","V"):"C",("P","S"):"C",("P","T"):"C",("P","W"):"C",("P","Y"):"C",("P","V"):"C",("S","T"):"C",("S","W"):"C",("S","Y"):"C",("S","V"):"C",("T","W"):"C",("T","Y"):"C",("T","V"):"C",("W","Y"):"C",("W","V"):"C",("Y","V"):"C",("H","R"):"C",("K","R"):"C",("D","R"):"R",("E","R"):"R",("A","R"):"R",("N","R"):"R",("C","R"):"R",("Q","R"):"R",("G","R"):"R",("I","R"):"R",("L","R"):"R",("M","R"):"R",("F","R"):"R",("P","R"):"R",("S","R"):"R",("T","R"):"R",("W","R"):"R",("Y","R"):"R",("V","R"):"R",("K","H"):"C",("D","H"):"R",("E","H"):"R",("A","H"):"R",("N","H"):"R",("C","H"):"R",("Q","H"):"R",("G","H"):"R",("I","H"):"R",("L","H"):"R",("M","H"):"R",("F","H"):"R",("P","H"):"R",("S","H"):"R",("T","H"):"R",("W","H"):"R",("Y","H"):"R",("V","H"):"R",("D","K"):"R",("E","K"):"R",("A","K"):"R",("N","K"):"R",("C","K"):"R",("Q","K"):"R",("G","K"):"R",("I","K"):"R",("L","K"):"R",("M","K"):"R",("F","K"):"R",("P","K"):"R",("S","K"):"R",("T","K"):"R",("W","K"):"R",("Y","K"):"R",("V","K"):"R",("E","D"):"C",("A","D"):"R",("N","D"):"R",("C","D"):"R",("Q","D"):"R",("G","D"):"R",("I","D"):"R",("L","D"):"R",("M","D"):"R",("F","D"):"R",("P","D"):"R",("S","D"):"R",("T","D"):"R",("W","D"):"R",("Y","D"):"R",("V","D"):"R",("A","E"):"R",("N","E"):"R",("C","E"):"R",("Q","E"):"R",("G","E"):"R",("I","E"):"R",("L","E"):"R",("M","E"):"R",("F","E"):"R",("P","E"):"R",("S","E"):"R",("T","E"):"R",("W","E"):"R",("Y","E"):"R",("V","E"):"R",("N","A"):"C",("C","A"):"C",("Q","A"):"C",("G","A"):"C",("I","A"):"C",("L","A"):"C",("M","A"):"C",("F","A"):"C",("P","A"):"C",("S","A"):"C",("T","A"):"C",("W","A"):"C",("Y","A"):"C",("V","A"):"C",("C","N"):"C",("Q","N"):"C",("G","N"):"C",("I","N"):"C",("L","N"):"C",("M","N"):"C",("F","N"):"C",("P","N"):"C",("S","N"):"C",("T","N"):"C",("W","N"):"C",("Y","N"):"C",("V","N"):"C",("Q","C"):"C",("G","C"):"C",("I","C"):"C",("L","C"):"C",("M","C"):"C",("F","C"):"C",("P","C"):"C",("S","C"):"C",("T","C"):"C",("W","C"):"C",("Y","C"):"C",("V","C"):"C",("G","Q"):"C",("I","Q"):"C",("L","Q"):"C",("M","Q"):"C",("F","Q"):"C",("P","Q"):"C",("S","Q"):"C",("T","Q"):"C",("W","Q"):"C",("Y","Q"):"C",("V","Q"):"C",("I","G"):"C",("L","G"):"C",("M","G"):"C",("F","G"):"C",("P","G"):"C",("S","G"):"C",("T","G"):"C",("W","G"):"C",("Y","G"):"C",("V","G"):"C",("L","I"):"C",("M","I"):"C",("F","I"):"C",("P","I"):"C",("S","I"):"C",("T","I"):"C",("W","I"):"C",("Y","I"):"C",("V","I"):"C",("M","L"):"C",("F","L"):"C",("P","L"):"C",("S","L"):"C",("T","L"):"C",("W","L"):"C",("Y","L"):"C",("V","L"):"C",("F","M"):"C",("P","M"):"C",("S","M"):"C",("T","M"):"C",("W","M"):"C",("Y","M"):"C",("V","M"):"C",("P","F"):"C",("S","F"):"C",("T","F"):"C",("W","F"):"C",("Y","F"):"C",("V","F"):"C",("S","P"):"C",("T","P"):"C",("W","P"):"C",("Y","P"):"C",("V","P"):"C",("T","S"):"C",("W","S"):"C",("Y","S"):"C",("V","S"):"C",("W","T"):"C",("Y","T"):"C",("V","T"):"C",("Y","W"):"C",("V","W"):"C",("V","Y"):"C",("R","*"):"R",("H","*"):"R",("K","*"):"R",("D","*"):"R",("E","*"):"R",("A","*"):"R",("N","*"):"R",("C","*"):"R",("Q","*"):"R",("G","*"):"R",("I","*"):"R",("L","*"):"R",("M","*"):"R",("F","*"):"R",("P","*"):"R",("S","*"):"R",("T","*"):"R",("W","*"):"R",("Y","*"):"R",("V","*"):"R",("*","R"):"R",("*","H"):"R",("*","K"):"R",("*","D"):"R",("*","E"):"R",("*","A"):"R",("*","N"):"R",("*","C"):"R",("*","Q"):"R",("*","G"):"R",("*","I"):"R",("*","L"):"R",("*","M"):"R",("*","F"):"R",("*","P"):"R",("*","S"):"R",("*","T"):"R",("*","W"):"R",("*","Y"):"R",("*","V"):"R"},2:{("R","H"):"C",("R","K"):"C",("R","D"):"C",("R","E"):"C",("R","A"):"R",("R","N"):"C",("R","C"):"C",("R","Q"):"C",("R","G"):"C",("R","I"):"R",("R","L"):"R",("R","M"):"R",("R","F"):"R",("R","P"):"R",("R","S"):"C",("R","T"):"C",("R","W"):"R",("R","Y"):"C",("R","V"):"R",("H","K"):"C",("H","D"):"C",("H","E"):"C",("H","A"):"R",("H","N"):"C",("H","C"):"C",("H","Q"):"C",("H","G"):"C",("H","I"):"R",("H","L"):"R",("H","M"):"R",("H","F"):"R",("H","P"):"R",("H","S"):"C",("H","T"):"C",("H","W"):"R",("H","Y"):"C",("H","V"):"R",("K","D"):"C",("K","E"):"C",("K","A"):"R",("K","N"):"C",("K","C"):"C",("K","Q"):"C",("K","G"):"C",("K","I"):"R",("K","L"):"R",("K","M"):"R",("K","F"):"R",("K","P"):"R",("K","S"):"C",("K","T"):"C",("K","W"):"R",("K","Y"):"C",("K","V"):"R",("D","E"):"C",("D","A"):"R",("D","N"):"C",("D","C"):"C",("D","Q"):"C",("D","G"):"C",("D","I"):"R",("D","L"):"R",("D","M"):"R",("D","F"):"R",("D","P"):"R",("D","S"):"C",("D","T"):"C",("D","W"):"R",("D","Y"):"C",("D","V"):"R",("E","A"):"R",("E","N"):"C",("E","C"):"C",("E","Q"):"C",("E","G"):"C",("E","I"):"R",("E","L"):"R",("E","M"):"R",("E","F"):"R",("E","P"):"R",("E","S"):"C",("E","T"):"C",("E","W"):"R",("E","Y"):"C",("E","V"):"R",("A","N"):"R",("A","C"):"R",("A","Q"):"R",("A","G"):"R",("A","I"):"C",("A","L"):"C",("A","M"):"C",("A","F"):"C",("A","P"):"C",("A","S"):"R",("A","T"):"R",("A","W"):"C",("A","Y"):"R",("A","V"):"C",("N","C"):"C",("N","Q"):"C",("N","G"):"C",("N","I"):"R",("N","L"):"R",("N","M"):"R",("N","F"):"R",("N","P"):"R",("N","S"):"C",("N","T"):"C",("N","W"):"R",("N","Y"):"C",("N","V"):"R",("C","Q"):"C",("C","G"):"C",("C","I"):"R",("C","L"):"R",("C","M"):"R",("C","F"):"R",("C","P"):"R",("C","S"):"C",("C","T"):"C",("C","W"):"R",("C","Y"):"C",("C","V"):"R",("Q","G"):"C",("Q","I"):"R",("Q","L"):"R",("Q","M"):"R",("Q","F"):"R",("Q","P"):"R",("Q","S"):"C",("Q","T"):"C",("Q","W"):"R",("Q","Y"):"C",("Q","V"):"R",("G","I"):"R",("G","L"):"R",("G","M"):"R",("G","F"):"R",("G","P"):"R",("G","S"):"C",("G","T"):"C",("G","W"):"R",("G","Y"):"C",("G","V"):"R",("I","L"):"C",("I","M"):"C",("I","F"):"C",("I","P"):"C",("I","S"):"R",("I","T"):"R",("I","W"):"C",("I","Y"):"R",("I","V"):"C",("L","M"):"C",("L","F"):"C",("L","P"):"C",("L","S"):"R",("L","T"):"R",("L","W"):"C",("L","Y"):"R",("L","V"):"C",("M","F"):"C",("M","P"):"C",("M","S"):"R",("M","T"):"R",("M","W"):"C",("M","Y"):"R",("M","V"):"C",("F","P"):"C",("F","S"):"R",("F","T"):"R",("F","W"):"C",("F","Y"):"R",("F","V"):"C",("P","S"):"R",("P","T"):"R",("P","W"):"C",("P","Y"):"R",("P","V"):"C",("S","T"):"C",("S","W"):"R",("S","Y"):"C",("S","V"):"R",("T","W"):"R",("T","Y"):"C",("T","V"):"R",("W","Y"):"R",("W","V"):"C",("Y","V"):"R",("H","R"):"C",("K","R"):"C",("D","R"):"C",("E","R"):"C",("A","R"):"R",("N","R"):"C",("C","R"):"C",("Q","R"):"C",("G","R"):"C",("I","R"):"R",("L","R"):"R",("M","R"):"R",("F","R"):"R",("P","R"):"R",("S","R"):"C",("T","R"):"C",("W","R"):"R",("Y","R"):"C",("V","R"):"R",("K","H"):"C",("D","H"):"C",("E","H"):"C",("A","H"):"R",("N","H"):"C",("C","H"):"C",("Q","H"):"C",("G","H"):"C",("I","H"):"R",("L","H"):"R",("M","H"):"R",("F","H"):"R",("P","H"):"R",("S","H"):"C",("T","H"):"C",("W","H"):"R",("Y","H"):"C",("V","H"):"R",("D","K"):"C",("E","K"):"C",("A","K"):"R",("N","K"):"C",("C","K"):"C",("Q","K"):"C",("G","K"):"C",("I","K"):"R",("L","K"):"R",("M","K"):"R",("F","K"):"R",("P","K"):"R",("S","K"):"C",("T","K"):"C",("W","K"):"R",("Y","K"):"C",("V","K"):"R",("E","D"):"C",("A","D"):"R",("N","D"):"C",("C","D"):"C",("Q","D"):"C",("G","D"):"C",("I","D"):"R",("L","D"):"R",("M","D"):"R",("F","D"):"R",("P","D"):"R",("S","D"):"C",("T","D"):"C",("W","D"):"R",("Y","D"):"C",("V","D"):"R",("A","E"):"R",("N","E"):"C",("C","E"):"C",("Q","E"):"C",("G","E"):"C",("I","E"):"R",("L","E"):"R",("M","E"):"R",("F","E"):"R",("P","E"):"R",("S","E"):"C",("T","E"):"C",("W","E"):"R",("Y","E"):"C",("V","E"):"R",("N","A"):"R",("C","A"):"R",("Q","A"):"R",("G","A"):"R",("I","A"):"C",("L","A"):"C",("M","A"):"C",("F","A"):"C",("P","A"):"C",("S","A"):"R",("T","A"):"R",("W","A"):"C",("Y","A"):"R",("V","A"):"C",("C","N"):"C",("Q","N"):"C",("G","N"):"C",("I","N"):"R",("L","N"):"R",("M","N"):"R",("F","N"):"R",("P","N"):"R",("S","N"):"C",("T","N"):"C",("W","N"):"R",("Y","N"):"C",("V","N"):"R",("Q","C"):"C",("G","C"):"C",("I","C"):"R",("L","C"):"R",("M","C"):"R",("F","C"):"R",("P","C"):"R",("S","C"):"C",("T","C"):"C",("W","C"):"R",("Y","C"):"C",("V","C"):"R",("G","Q"):"C",("I","Q"):"R",("L","Q"):"R",("M","Q"):"R",("F","Q"):"R",("P","Q"):"R",("S","Q"):"C",("T","Q"):"C",("W","Q"):"R",("Y","Q"):"C",("V","Q"):"R",("I","G"):"R",("L","G"):"R",("M","G"):"R",("F","G"):"R",("P","G"):"R",("S","G"):"C",("T","G"):"C",("W","G"):"R",("Y","G"):"C",("V","G"):"R",("L","I"):"C",("M","I"):"C",("F","I"):"C",("P","I"):"C",("S","I"):"R",("T","I"):"R",("W","I"):"C",("Y","I"):"R",("V","I"):"C",("M","L"):"C",("F","L"):"C",("P","L"):"C",("S","L"):"R",("T","L"):"R",("W","L"):"C",("Y","L"):"R",("V","L"):"C",("F","M"):"C",("P","M"):"C",("S","M"):"R",("T","M"):"R",("W","M"):"C",("Y","M"):"R",("V","M"):"C",("P","F"):"C",("S","F"):"R",("T","F"):"R",("W","F"):"C",("Y","F"):"R",("V","F"):"C",("S","P"):"R",("T","P"):"R",("W","P"):"C",("Y","P"):"R",("V","P"):"C",("T","S"):"C",("W","S"):"R",("Y","S"):"C",("V","S"):"R",("W","T"):"R",("Y","T"):"C",("V","T"):"R",("Y","W"):"R",("V","W"):"C",("V","Y"):"R",("R","*"):"R",("H","*"):"R",("K","*"):"R",("D","*"):"R",("E","*"):"R",("A","*"):"R",("N","*"):"R",("C","*"):"R",("Q","*"):"R",("G","*"):"R",("I","*"):"R",("L","*"):"R",("M","*"):"R",("F","*"):"R",("P","*"):"R",("S","*"):"R",("T","*"):"R",("W","*"):"R",("Y","*"):"R",("V","*"):"R",("*","R"):"R",("*","H"):"R",("*","K"):"R",("*","D"):"R",("*","E"):"R",("*","A"):"R",("*","N"):"R",("*","C"):"R",("*","Q"):"R",("*","G"):"R",("*","I"):"R",("*","L"):"R",("*","M"):"R",("*","F"):"R",("*","P"):"R",("*","S"):"R",("*","T"):"R",("*","W"):"R",("*","Y"):"R",("*","V"):"R"},3:{("R","H"):"C",("R","K"):"C",("R","D"):"R",("R","E"):"R",("R","A"):"R",("R","N"):"R",("R","C"):"R",("R","Q"):"R",("R","G"):"R",("R","I"):"R",("R","L"):"R",("R","M"):"R",("R","F"):"R",("R","P"):"R",("R","S"):"R",("R","T"):"R",("R","W"):"R",("R","Y"):"R",("R","V"):"R",("H","K"):"C",("H","D"):"R",("H","E"):"R",("H","A"):"R",("H","N"):"R",("H","C"):"R",("H","Q"):"R",("H","G"):"R",("H","I"):"R",("H","L"):"R",("H","M"):"R",("H","F"):"R",("H","P"):"R",("H","S"):"R",("H","T"):"R",("H","W"):"R",("H","Y"):"R",("H","V"):"R",("K","D"):"R",("K","E"):"R",("K","A"):"R",("K","N"):"R",("K","C"):"R",("K","Q"):"R",("K","G"):"R",("K","I"):"R",("K","L"):"R",("K","M"):"R",("K","F"):"R",("K","P"):"R",("K","S"):"R",("K","T"):"R",("K","W"):"R",("K","Y"):"R",("K","V"):"R",("D","E"):"C",("D","A"):"R",("D","N"):"C",("D","C"):"R",("D","Q"):"C",("D","G"):"R",("D","I"):"R",("D","L"):"R",("D","M"):"R",("D","F"):"R",("D","P"):"R",("D","S"):"R",("D","T"):"R",("D","W"):"R",("D","Y"):"R",("D","V"):"R",("E","A"):"R",("E","N"):"C",("E","C"):"R",("E","Q"):"C",("E","G"):"R",("E","I"):"R",("E","L"):"R",("E","M"):"R",("E","F"):"R",("E","P"):"R",("E","S"):"R",("E","T"):"R",("E","W"):"R",("E","Y"):"R",("E","V"):"R",("A","N"):"R",("A","C"):"R",("A","Q"):"R",("A","G"):"C",("A","I"):"R",("A","L"):"R",("A","M"):"R",("A","F"):"R",("A","P"):"C",("A","S"):"C",("A","T"):"C",("A","W"):"R",("A","Y"):"R",("A","V"):"R",("N","C"):"R",("N","Q"):"C",("N","G"):"R",("N","I"):"R",("N","L"):"R",("N","M"):"R",("N","F"):"R",("N","P"):"R",("N","S"):"R",("N","T"):"R",("N","W"):"R",("N","Y"):"R",("N","V"):"R",("C","Q"):"R",("C","G"):"R",("C","I"):"R",("C","L"):"R",("C","M"):"R",("C","F"):"R",("C","P"):"R",("C","S"):"R",("C","T"):"R",("C","W"):"R",("C","Y"):"R",("C","V"):"R",("Q","G"):"R",("Q","I"):"R",("Q","L"):"R",("Q","M"):"R",("Q","F"):"R",("Q","P"):"R",("Q","S"):"R",("Q","T"):"R",("Q","W"):"R",("Q","Y"):"R",("Q","V"):"R",("G","I"):"R",("G","L"):"R",("G","M"):"R",("G","F"):"R",("G","P"):"C",("G","S"):"C",("G","T"):"C",("G","W"):"R",("G","Y"):"R",("G","V"):"R",("I","L"):"C",("I","M"):"C",("I","F"):"R",("I","P"):"R",("I","S"):"R",("I","T"):"R",("I","W"):"R",("I","Y"):"R",("I","V"):"C",("L","M"):"C",("L","F"):"R",("L","P"):"R",("L","S"):"R",("L","T"):"R",("L","W"):"R",("L","Y"):"R",("L","V"):"C",("M","F"):"R",("M","P"):"R",("M","S"):"R",("M","T"):"R",("M","W"):"R",("M","Y"):"R",("M","V"):"C",("F","P"):"R",("F","S"):"C",("F","T"):"C",("F","W"):"R",("F","Y"):"R",("F","V"):"R",("P","S"):"C",("P","T"):"C",("P","W"):"R",("P","Y"):"R",("P","V"):"R",("S","T"):"C",("S","W"):"R",("S","Y"):"R",("S","V"):"R",("T","W"):"R",("T","Y"):"R",("T","V"):"R",("W","Y"):"C",("W","V"):"R",("Y","V"):"R",("H","R"):"C",("K","R"):"C",("D","R"):"R",("E","R"):"R",("A","R"):"R",("N","R"):"R",("C","R"):"R",("Q","R"):"R",("G","R"):"R",("I","R"):"R",("L","R"):"R",("M","R"):"R",("F","R"):"R",("P","R"):"R",("S","R"):"R",("T","R"):"R",("W","R"):"R",("Y","R"):"R",("V","R"):"R",("K","H"):"C",("D","H"):"R",("E","H"):"R",("A","H"):"R",("N","H"):"R",("C","H"):"R",("Q","H"):"R",("G","H"):"R",("I","H"):"R",("L","H"):"R",("M","H"):"R",("F","H"):"R",("P","H"):"R",("S","H"):"R",("T","H"):"R",("W","H"):"R",("Y","H"):"R",("V","H"):"R",("D","K"):"R",("E","K"):"R",("A","K"):"R",("N","K"):"R",("C","K"):"R",("Q","K"):"R",("G","K"):"R",("I","K"):"R",("L","K"):"R",("M","K"):"R",("F","K"):"R",("P","K"):"R",("S","K"):"R",("T","K"):"R",("W","K"):"R",("Y","K"):"R",("V","K"):"R",("E","D"):"C",("A","D"):"R",("N","D"):"C",("C","D"):"R",("Q","D"):"C",("G","D"):"R",("I","D"):"R",("L","D"):"R",("M","D"):"R",("F","D"):"R",("P","D"):"R",("S","D"):"R",("T","D"):"R",("W","D"):"R",("Y","D"):"R",("V","D"):"R",("A","E"):"R",("N","E"):"C",("C","E"):"R",("Q","E"):"C",("G","E"):"R",("I","E"):"R",("L","E"):"R",("M","E"):"R",("F","E"):"R",("P","E"):"R",("S","E"):"R",("T","E"):"R",("W","E"):"R",("Y","E"):"R",("V","E"):"R",("N","A"):"R",("C","A"):"R",("Q","A"):"R",("G","A"):"C",("I","A"):"R",("L","A"):"R",("M","A"):"R",("F","A"):"R",("P","A"):"C",("S","A"):"C",("T","A"):"C",("W","A"):"R",("Y","A"):"R",("V","A"):"R",("C","N"):"R",("Q","N"):"C",("G","N"):"R",("I","N"):"R",("L","N"):"R",("M","N"):"R",("F","N"):"R",("P","N"):"R",("S","N"):"R",("T","N"):"R",("W","N"):"R",("Y","N"):"R",("V","N"):"R",("Q","C"):"R",("G","C"):"R",("I","C"):"R",("L","C"):"R",("M","C"):"R",("F","C"):"R",("P","C"):"R",("S","C"):"R",("T","C"):"R",("W","C"):"R",("Y","C"):"R",("V","C"):"R",("G","Q"):"R",("I","Q"):"R",("L","Q"):"R",("M","Q"):"R",("F","Q"):"R",("P","Q"):"R",("S","Q"):"R",("T","Q"):"R",("W","Q"):"R",("Y","Q"):"R",("V","Q"):"R",("I","G"):"R",("L","G"):"R",("M","G"):"R",("F","G"):"R",("P","G"):"C",("S","G"):"C",("T","G"):"C",("W","G"):"R",("Y","G"):"R",("V","G"):"R",("L","I"):"C",("M","I"):"C",("F","I"):"R",("P","I"):"R",("S","I"):"R",("T","I"):"R",("W","I"):"R",("Y","I"):"R",("V","I"):"C",("M","L"):"C",("F","L"):"R",("P","L"):"R",("S","L"):"R",("T","L"):"R",("W","L"):"R",("Y","L"):"R",("V","L"):"C",("F","M"):"R",("P","M"):"R",("S","M"):"R",("T","M"):"R",("W","M"):"R",("Y","M"):"R",("V","M"):"C",("P","F"):"R",("S","F"):"C",("T","F"):"C",("W","F"):"R",("Y","F"):"R",("V","F"):"R",("S","P"):"C",("T","P"):"C",("W","P"):"R",("Y","P"):"R",("V","P"):"R",("T","S"):"C",("W","S"):"R",("Y","S"):"R",("V","S"):"R",("W","T"):"R",("Y","T"):"R",("V","T"):"R",("Y","W"):"C",("V","W"):"R",("V","Y"):"R",("R","*"):"R",("H","*"):"R",("K","*"):"R",("D","*"):"R",("E","*"):"R",("A","*"):"R",("N","*"):"R",("C","*"):"R",("Q","*"):"R",("G","*"):"R",("I","*"):"R",("L","*"):"R",("M","*"):"R",("F","*"):"R",("P","*"):"R",("S","*"):"R",("T","*"):"R",("W","*"):"R",("Y","*"):"R",("V","*"):"R",("*","R"):"R",("*","H"):"R",("*","K"):"R",("*","D"):"R",("*","E"):"R",("*","A"):"R",("*","N"):"R",("*","C"):"R",("*","Q"):"R",("*","G"):"R",("*","I"):"R",("*","L"):"R",("*","M"):"R",("*","F"):"R",("*","P"):"R",("*","S"):"R",("*","T"):"R",("*","W"):"R",("*","Y"):"R",("*","V"):"R"},4:{("R","H"):"C",("R","K"):"C",("R","D"):"R",("R","E"):"R",("R","A"):"R",("R","N"):"R",("R","C"):"R",("R","Q"):"C",("R","G"):"R",("R","I"):"R",("R","L"):"R",("R","M"):"R",("R","F"):"C",("R","P"):"R",("R","S"):"R",("R","T"):"R",("R","W"):"C",("R","Y"):"C",("R","V"):"R",("H","K"):"C",("H","D"):"R",("H","E"):"R",("H","A"):"R",("H","N"):"R",("H","C"):"R",("H","Q"):"C",("H","G"):"R",("H","I"):"R",("H","L"):"R",("H","M"):"R",("H","F"):"C",("H","P"):"R",("H","S"):"R",("H","T"):"R",("H","W"):"C",("H","Y"):"C",("H","V"):"R",("K","D"):"R",("K","E"):"R",("K","A"):"R",("K","N"):"R",("K","C"):"R",("K","Q"):"C",("K","G"):"R",("K","I"):"R",("K","L"):"R",("K","M"):"R",("K","F"):"C",("K","P"):"R",("K","S"):"R",("K","T"):"R",("K","W"):"C",("K","Y"):"C",("K","V"):"R",("D","E"):"C",("D","A"):"R",("D","N"):"R",("D","C"):"R",("D","Q"):"R",("D","G"):"R",("D","I"):"R",("D","L"):"R",("D","M"):"R",("D","F"):"R",("D","P"):"R",("D","S"):"R",("D","T"):"R",("D","W"):"R",("D","Y"):"R",("D","V"):"R",("E","A"):"R",("E","N"):"R",("E","C"):"R",("E","Q"):"R",("E","G"):"R",("E","I"):"R",("E","L"):"R",("E","M"):"R",("E","F"):"R",("E","P"):"R",("E","S"):"R",("E","T"):"R",("E","W"):"R",("E","Y"):"R",("E","V"):"R",("A","N"):"C",("A","C"):"C",("A","Q"):"R",("A","G"):"C",("A","I"):"R",("A","L"):"R",("A","M"):"R",("A","F"):"R",("A","P"):"C",("A","S"):"C",("A","T"):"C",("A","W"):"R",("A","Y"):"R",("A","V"):"R",("N","C"):"C",("N","Q"):"R",("N","G"):"C",("N","I"):"R",("N","L"):"R",("N","M"):"R",("N","F"):"R",("N","P"):"C",("N","S"):"C",("N","T"):"C",("N","W"):"R",("N","Y"):"R",("N","V"):"R",("C","Q"):"R",("C","G"):"C",("C","I"):"R",("C","L"):"R",("C","M"):"R",("C","F"):"R",("C","P"):"C",("C","S"):"C",("C","T"):"C",("C","W"):"R",("C","Y"):"R",("C","V"):"R",("Q","G"):"R",("Q","I"):"R",("Q","L"):"R",("Q","M"):"R",("Q","F"):"C",("Q","P"):"R",("Q","S"):"R",("Q","T"):"R",("Q","W"):"C",("Q","Y"):"C",("Q","V"):"R",("G","I"):"R",("G","L"):"R",("G","M"):"R",("G","F"):"R",("G","P"):"C",("G","S"):"C",("G","T"):"C",("G","W"):"R",("G","Y"):"R",("G","V"):"R",("I","L"):"C",("I","M"):"C",("I","F"):"R",("I","P"):"R",("I","S"):"R",("I","T"):"R",("I","W"):"R",("I","Y"):"R",("I","V"):"C",("L","M"):"C",("L","F"):"R",("L","P"):"R",("L","S"):"R",("L","T"):"R",("L","W"):"R",("L","Y"):"R",("L","V"):"C",("M","F"):"R",("M","P"):"R",("M","S"):"R",("M","T"):"R",("M","W"):"R",("M","Y"):"R",("M","V"):"C",("F","P"):"R",("F","S"):"R",("F","T"):"R",("F","W"):"C",("F","Y"):"C",("F","V"):"R",("P","S"):"C",("P","T"):"C",("P","W"):"R",("P","Y"):"R",("P","V"):"R",("S","T"):"C",("S","W"):"R",("S","Y"):"R",("S","V"):"R",("T","W"):"R",("T","Y"):"R",("T","V"):"R",("W","Y"):"C",("W","V"):"R",("Y","V"):"R",("H","R"):"C",("K","R"):"C",("D","R"):"R",("E","R"):"R",("A","R"):"R",("N","R"):"R",("C","R"):"R",("Q","R"):"C",("G","R"):"R",("I","R"):"R",("L","R"):"R",("M","R"):"R",("F","R"):"C",("P","R"):"R",("S","R"):"R",("T","R"):"R",("W","R"):"C",("Y","R"):"C",("V","R"):"R",("K","H"):"C",("D","H"):"R",("E","H"):"R",("A","H"):"R",("N","H"):"R",("C","H"):"R",("Q","H"):"C",("G","H"):"R",("I","H"):"R",("L","H"):"R",("M","H"):"R",("F","H"):"C",("P","H"):"R",("S","H"):"R",("T","H"):"R",("W","H"):"C",("Y","H"):"C",("V","H"):"R",("D","K"):"R",("E","K"):"R",("A","K"):"R",("N","K"):"R",("C","K"):"R",("Q","K"):"C",("G","K"):"R",("I","K"):"R",("L","K"):"R",("M","K"):"R",("F","K"):"C",("P","K"):"R",("S","K"):"R",("T","K"):"R",("W","K"):"C",("Y","K"):"C",("V","K"):"R",("E","D"):"C",("A","D"):"R",("N","D"):"R",("C","D"):"R",("Q","D"):"R",("G","D"):"R",("I","D"):"R",("L","D"):"R",("M","D"):"R",("F","D"):"R",("P","D"):"R",("S","D"):"R",("T","D"):"R",("W","D"):"R",("Y","D"):"R",("V","D"):"R",("A","E"):"R",("N","E"):"R",("C","E"):"R",("Q","E"):"R",("G","E"):"R",("I","E"):"R",("L","E"):"R",("M","E"):"R",("F","E"):"R",("P","E"):"R",("S","E"):"R",("T","E"):"R",("W","E"):"R",("Y","E"):"R",("V","E"):"R",("N","A"):"C",("C","A"):"C",("Q","A"):"R",("G","A"):"C",("I","A"):"R",("L","A"):"R",("M","A"):"R",("F","A"):"R",("P","A"):"C",("S","A"):"C",("T","A"):"C",("W","A"):"R",("Y","A"):"R",("V","A"):"R",("C","N"):"C",("Q","N"):"R",("G","N"):"C",("I","N"):"R",("L","N"):"R",("M","N"):"R",("F","N"):"R",("P","N"):"C",("S","N"):"C",("T","N"):"C",("W","N"):"R",("Y","N"):"R",("V","N"):"R",("Q","C"):"R",("G","C"):"C",("I","C"):"R",("L","C"):"R",("M","C"):"R",("F","C"):"R",("P","C"):"C",("S","C"):"C",("T","C"):"C",("W","C"):"R",("Y","C"):"R",("V","C"):"R",("G","Q"):"R",("I","Q"):"R",("L","Q"):"R",("M","Q"):"R",("F","Q"):"C",("P","Q"):"R",("S","Q"):"R",("T","Q"):"R",("W","Q"):"C",("Y","Q"):"C",("V","Q"):"R",("I","G"):"R",("L","G"):"R",("M","G"):"R",("F","G"):"R",("P","G"):"C",("S","G"):"C",("T","G"):"C",("W","G"):"R",("Y","G"):"R",("V","G"):"R",("L","I"):"C",("M","I"):"C",("F","I"):"R",("P","I"):"R",("S","I"):"R",("T","I"):"R",("W","I"):"R",("Y","I"):"R",("V","I"):"C",("M","L"):"C",("F","L"):"R",("P","L"):"R",("S","L"):"R",("T","L"):"R",("W","L"):"R",("Y","L"):"R",("V","L"):"C",("F","M"):"R",("P","M"):"R",("S","M"):"R",("T","M"):"R",("W","M"):"R",("Y","M"):"R",("V","M"):"C",("P","F"):"R",("S","F"):"R",("T","F"):"R",("W","F"):"C",("Y","F"):"C",("V","F"):"R",("S","P"):"C",("T","P"):"C",("W","P"):"R",("Y","P"):"R",("V","P"):"R",("T","S"):"C",("W","S"):"R",("Y","S"):"R",("V","S"):"R",("W","T"):"R",("Y","T"):"R",("V","T"):"R",("Y","W"):"C",("V","W"):"R",("V","Y"):"R",("R","*"):"R",("H","*"):"R",("K","*"):"R",("D","*"):"R",("E","*"):"R",("A","*"):"R",("N","*"):"R",("C","*"):"R",("Q","*"):"R",("G","*"):"R",("I","*"):"R",("L","*"):"R",("M","*"):"R",("F","*"):"R",("P","*"):"R",("S","*"):"R",("T","*"):"R",("W","*"):"R",("Y","*"):"R",("V","*"):"R",("*","R"):"R",("*","H"):"R",("*","K"):"R",("*","D"):"R",("*","E"):"R",("*","A"):"R",("*","N"):"R",("*","C"):"R",("*","Q"):"R",("*","G"):"R",("*","I"):"R",("*","L"):"R",("*","M"):"R",("*","F"):"R",("*","P"):"R",("*","S"):"R",("*","T"):"R",("*","W"):"R",("*","Y"):"R",("*","V"):"R"},5:{("R","H"):"C",("R","K"):"C",("R","D"):"R",("R","E"):"R",("R","A"):"R",("R","N"):"R",("R","C"):"R",("R","Q"):"R",("R","G"):"R",("R","I"):"R",("R","L"):"R",("R","M"):"R",("R","F"):"R",("R","P"):"R",("R","S"):"R",("R","T"):"R",("R","W"):"R",("R","Y"):"R",("R","V"):"R",("H","K"):"C",("H","D"):"R",("H","E"):"R",("H","A"):"R",("H","N"):"R",("H","C"):"R",("H","Q"):"R",("H","G"):"R",("H","I"):"R",("H","L"):"R",("H","M"):"R",("H","F"):"R",("H","P"):"R",("H","S"):"R",("H","T"):"R",("H","W"):"R",("H","Y"):"R",("H","V"):"R",("K","D"):"R",("K","E"):"R",("K","A"):"R",("K","N"):"R",("K","C"):"R",("K","Q"):"R",("K","G"):"R",("K","I"):"R",("K","L"):"R",("K","M"):"R",("K","F"):"R",("K","P"):"R",("K","S"):"R",("K","T"):"R",("K","W"):"R",("K","Y"):"R",("K","V"):"R",("D","E"):"C",("D","A"):"R",("D","N"):"R",("D","C"):"R",("D","Q"):"R",("D","G"):"R",("D","I"):"R",("D","L"):"R",("D","M"):"R",("D","F"):"R",("D","P"):"R",("D","S"):"R",("D","T"):"R",("D","W"):"R",("D","Y"):"R",("D","V"):"R",("E","A"):"R",("E","N"):"R",("E","C"):"R",("E","Q"):"R",("E","G"):"R",("E","I"):"R",("E","L"):"R",("E","M"):"R",("E","F"):"R",("E","P"):"R",("E","S"):"R",("E","T"):"R",("E","W"):"R",("E","Y"):"R",("E","V"):"R",("A","N"):"C",("A","C"):"C",("A","Q"):"C",("A","G"):"C",("A","I"):"C",("A","L"):"C",("A","M"):"C",("A","F"):"R",("A","P"):"C",("A","S"):"C",("A","T"):"C",("A","W"):"R",("A","Y"):"R",("A","V"):"C",("N","C"):"C",("N","Q"):"C",("N","G"):"C",("N","I"):"C",("N","L"):"C",("N","M"):"C",("N","F"):"R",("N","P"):"C",("N","S"):"C",("N","T"):"C",("N","W"):"R",("N","Y"):"R",("N","V"):"C",("C","Q"):"C",("C","G"):"C",("C","I"):"C",("C","L"):"C",("C","M"):"C",("C","F"):"R",("C","P"):"C",("C","S"):"C",("C","T"):"C",("C","W"):"R",("C","Y"):"R",("C","V"):"C",("Q","G"):"C",("Q","I"):"C",("Q","L"):"C",("Q","M"):"C",("Q","F"):"R",("Q","P"):"C",("Q","S"):"C",("Q","T"):"C",("Q","W"):"R",("Q","Y"):"R",("Q","V"):"C",("G","I"):"C",("G","L"):"C",("G","M"):"C",("G","F"):"R",("G","P"):"C",("G","S"):"C",("G","T"):"C",("G","W"):"R",("G","Y"):"R",("G","V"):"C",("I","L"):"C",("I","M"):"C",("I","F"):"R",("I","P"):"C",("I","S"):"C",("I","T"):"C",("I","W"):"R",("I","Y"):"R",("I","V"):"C",("L","M"):"C",("L","F"):"R",("L","P"):"C",("L","S"):"C",("L","T"):"C",("L","W"):"R",("L","Y"):"R",("L","V"):"C",("M","F"):"R",("M","P"):"C",("M","S"):"C",("M","T"):"C",("M","W"):"R",("M","Y"):"R",("M","V"):"C",("F","P"):"R",("F","S"):"R",("F","T"):"R",("F","W"):"C",("F","Y"):"C",("F","V"):"R",("P","S"):"C",("P","T"):"C",("P","W"):"R",("P","Y"):"R",("P","V"):"C",("S","T"):"C",("S","W"):"R",("S","Y"):"R",("S","V"):"R",("T","W"):"R",("T","Y"):"R",("T","V"):"C",("W","Y"):"C",("W","V"):"R",("Y","V"):"R",("H","R"):"C",("K","R"):"C",("D","R"):"R",("E","R"):"R",("A","R"):"R",("N","R"):"R",("C","R"):"R",("Q","R"):"R",("G","R"):"R",("I","R"):"R",("L","R"):"R",("M","R"):"R",("F","R"):"R",("P","R"):"R",("S","R"):"R",("T","R"):"R",("W","R"):"R",("Y","R"):"R",("V","R"):"R",("K","H"):"C",("D","H"):"R",("E","H"):"R",("A","H"):"R",("N","H"):"R",("C","H"):"R",("Q","H"):"R",("G","H"):"R",("I","H"):"R",("L","H"):"R",("M","H"):"R",("F","H"):"R",("P","H"):"R",("S","H"):"R",("T","H"):"R",("W","H"):"R",("Y","H"):"R",("V","H"):"R",("D","K"):"R",("E","K"):"R",("A","K"):"R",("N","K"):"R",("C","K"):"R",("Q","K"):"R",("G","K"):"R",("I","K"):"R",("L","K"):"R",("M","K"):"R",("F","K"):"R",("P","K"):"R",("S","K"):"R",("T","K"):"R",("W","K"):"R",("Y","K"):"R",("V","K"):"R",("E","D"):"C",("A","D"):"R",("N","D"):"R",("C","D"):"R",("Q","D"):"R",("G","D"):"R",("I","D"):"R",("L","D"):"R",("M","D"):"R",("F","D"):"R",("P","D"):"R",("S","D"):"R",("T","D"):"R",("W","D"):"R",("Y","D"):"R",("V","D"):"R",("A","E"):"R",("N","E"):"R",("C","E"):"R",("Q","E"):"R",("G","E"):"R",("I","E"):"R",("L","E"):"R",("M","E"):"R",("F","E"):"R",("P","E"):"R",("S","E"):"R",("T","E"):"R",("W","E"):"R",("Y","E"):"R",("V","E"):"R",("N","A"):"C",("C","A"):"C",("Q","A"):"C",("G","A"):"C",("I","A"):"C",("L","A"):"C",("M","A"):"C",("F","A"):"R",("P","A"):"C",("S","A"):"C",("T","A"):"C",("W","A"):"R",("Y","A"):"R",("V","A"):"C",("C","N"):"C",("Q","N"):"C",("G","N"):"C",("I","N"):"C",("L","N"):"C",("M","N"):"C",("F","N"):"R",("P","N"):"C",("S","N"):"C",("T","N"):"C",("W","N"):"R",("Y","N"):"R",("V","N"):"C",("Q","C"):"C",("G","C"):"C",("I","C"):"C",("L","C"):"C",("M","C"):"C",("F","C"):"R",("P","C"):"C",("S","C"):"C",("T","C"):"C",("W","C"):"R",("Y","C"):"R",("V","C"):"C",("G","Q"):"C",("I","Q"):"C",("L","Q"):"C",("M","Q"):"C",("F","Q"):"R",("P","Q"):"C",("S","Q"):"C",("T","Q"):"C",("W","Q"):"R",("Y","Q"):"R",("V","Q"):"C",("I","G"):"C",("L","G"):"C",("M","G"):"C",("F","G"):"R",("P","G"):"C",("S","G"):"C",("T","G"):"C",("W","G"):"R",("Y","G"):"R",("V","G"):"C",("L","I"):"C",("M","I"):"C",("F","I"):"R",("P","I"):"C",("S","I"):"C",("T","I"):"C",("W","I"):"R",("Y","I"):"R",("V","I"):"C",("M","L"):"C",("F","L"):"R",("P","L"):"C",("S","L"):"C",("T","L"):"C",("W","L"):"R",("Y","L"):"R",("V","L"):"C",("F","M"):"R",("P","M"):"C",("S","M"):"C",("T","M"):"C",("W","M"):"R",("Y","M"):"R",("V","M"):"C",("P","F"):"R",("S","F"):"R",("T","F"):"R",("W","F"):"C",("Y","F"):"C",("V","F"):"R",("S","P"):"C",("T","P"):"C",("W","P"):"R",("Y","P"):"R",("V","P"):"C",("T","S"):"C",("W","S"):"R",("Y","S"):"R",("V","S"):"R",("W","T"):"R",("Y","T"):"R",("V","T"):"C",("Y","W"):"C",("V","W"):"R",("V","Y"):"R",("R","*"):"R",("H","*"):"R",("K","*"):"R",("D","*"):"R",("E","*"):"R",("A","*"):"R",("N","*"):"R",("C","*"):"R",("Q","*"):"R",("G","*"):"R",("I","*"):"R",("L","*"):"R",("M","*"):"R",("F","*"):"R",("P","*"):"R",("S","*"):"R",("T","*"):"R",("W","*"):"R",("Y","*"):"R",("V","*"):"R",("*","R"):"R",("*","H"):"R",("*","K"):"R",("*","D"):"R",("*","E"):"R",("*","A"):"R",("*","N"):"R",("*","C"):"R",("*","Q"):"R",("*","G"):"R",("*","I"):"R",("*","L"):"R",("*","M"):"R",("*","F"):"R",("*","P"):"R",("*","S"):"R",("*","T"):"R",("*","W"):"R",("*","Y"):"R",("*","V"):"R"},6:{("R","H"):"C",("R","K"):"C",("R","D"):"R",("R","E"):"R",("R","A"):"R",("R","N"):"R",("R","C"):"R",("R","Q"):"R",("R","G"):"R",("R","I"):"R",("R","L"):"R",("R","M"):"R",("R","F"):"R",("R","P"):"R",("R","S"):"R",("R","T"):"R",("R","W"):"R",("R","Y"):"R",("R","V"):"R",("H","K"):"C",("H","D"):"R",("H","E"):"R",("H","A"):"R",("H","N"):"R",("H","C"):"R",("H","Q"):"R",("H","G"):"R",("H","I"):"R",("H","L"):"R",("H","M"):"R",("H","F"):"R",("H","P"):"R",("H","S"):"R",("H","T"):"R",("H","W"):"R",("H","Y"):"R",("H","V"):"R",("K","D"):"R",("K","E"):"R",("K","A"):"R",("K","N"):"R",("K","C"):"R",("K","Q"):"R",("K","G"):"R",("K","I"):"R",("K","L"):"R",("K","M"):"R",("K","F"):"R",("K","P"):"R",("K","S"):"R",("K","T"):"R",("K","W"):"R",("K","Y"):"R",("K","V"):"R",("D","E"):"C",("D","A"):"R",("D","N"):"R",("D","C"):"R",("D","Q"):"R",("D","G"):"R",("D","I"):"R",("D","L"):"R",("D","M"):"R",("D","F"):"R",("D","P"):"R",("D","S"):"R",("D","T"):"R",("D","W"):"R",("D","Y"):"R",("D","V"):"R",("E","A"):"R",("E","N"):"R",("E","C"):"R",("E","Q"):"R",("E","G"):"R",("E","I"):"R",("E","L"):"R",("E","M"):"R",("E","F"):"R",("E","P"):"R",("E","S"):"R",("E","T"):"R",("E","W"):"R",("E","Y"):"R",("E","V"):"R",("A","N"):"R",("A","C"):"R",("A","Q"):"R",("A","G"):"C",("A","I"):"C",("A","L"):"C",("A","M"):"C",("A","F"):"C",("A","P"):"C",("A","S"):"R",("A","T"):"R",("A","W"):"C",("A","Y"):"R",("A","V"):"C",("N","C"):"C",("N","Q"):"C",("N","G"):"R",("N","I"):"R",("N","L"):"R",("N","M"):"R",("N","F"):"R",("N","P"):"R",("N","S"):"C",("N","T"):"C",("N","W"):"R",("N","Y"):"C",("N","V"):"R",("C","Q"):"C",("C","G"):"R",("C","I"):"R",("C","L"):"R",("C","M"):"R",("C","F"):"R",("C","P"):"R",("C","S"):"C",("C","T"):"C",("C","W"):"R",("C","Y"):"C",("C","V"):"R",("Q","G"):"R",("Q","I"):"R",("Q","L"):"R",("Q","M"):"R",("Q","F"):"R",("Q","P"):"R",("Q","S"):"C",("Q","T"):"C",("Q","W"):"R",("Q","Y"):"C",("Q","V"):"R",("G","I"):"C",("G","L"):"C",("G","M"):"C",("G","F"):"C",("G","P"):"C",("G","S"):"R",("G","T"):"R",("G","W"):"C",("G","Y"):"R",("G","V"):"C",("I","L"):"C",("I","M"):"C",("I","F"):"C",("I","P"):"C",("I","S"):"R",("I","T"):"R",("I","W"):"C",("I","Y"):"R",("I","V"):"C",("L","M"):"C",("L","F"):"C",("L","P"):"C",("L","S"):"R",("L","T"):"R",("L","W"):"C",("L","Y"):"R",("L","V"):"C",("M","F"):"C",("M","P"):"C",("M","S"):"R",("M","T"):"R",("M","W"):"C",("M","Y"):"R",("M","V"):"C",("F","P"):"C",("F","S"):"R",("F","T"):"R",("F","W"):"C",("F","Y"):"R",("F","V"):"C",("P","S"):"R",("P","T"):"R",("P","W"):"C",("P","Y"):"R",("P","V"):"C",("S","T"):"C",("S","W"):"R",("S","Y"):"C",("S","V"):"R",("T","W"):"R",("T","Y"):"C",("T","V"):"R",("W","Y"):"R",("W","V"):"C",("Y","V"):"R",("H","R"):"C",("K","R"):"C",("D","R"):"R",("E","R"):"R",("A","R"):"R",("N","R"):"R",("C","R"):"R",("Q","R"):"R",("G","R"):"R",("I","R"):"R",("L","R"):"R",("M","R"):"R",("F","R"):"R",("P","R"):"R",("S","R"):"R",("T","R"):"R",("W","R"):"R",("Y","R"):"R",("V","R"):"R",("K","H"):"C",("D","H"):"R",("E","H"):"R",("A","H"):"R",("N","H"):"R",("C","H"):"R",("Q","H"):"R",("G","H"):"R",("I","H"):"R",("L","H"):"R",("M","H"):"R",("F","H"):"R",("P","H"):"R",("S","H"):"R",("T","H"):"R",("W","H"):"R",("Y","H"):"R",("V","H"):"R",("D","K"):"R",("E","K"):"R",("A","K"):"R",("N","K"):"R",("C","K"):"R",("Q","K"):"R",("G","K"):"R",("I","K"):"R",("L","K"):"R",("M","K"):"R",("F","K"):"R",("P","K"):"R",("S","K"):"R",("T","K"):"R",("W","K"):"R",("Y","K"):"R",("V","K"):"R",("E","D"):"C",("A","D"):"R",("N","D"):"R",("C","D"):"R",("Q","D"):"R",("G","D"):"R",("I","D"):"R",("L","D"):"R",("M","D"):"R",("F","D"):"R",("P","D"):"R",("S","D"):"R",("T","D"):"R",("W","D"):"R",("Y","D"):"R",("V","D"):"R",("A","E"):"R",("N","E"):"R",("C","E"):"R",("Q","E"):"R",("G","E"):"R",("I","E"):"R",("L","E"):"R",("M","E"):"R",("F","E"):"R",("P","E"):"R",("S","E"):"R",("T","E"):"R",("W","E"):"R",("Y","E"):"R",("V","E"):"R",("N","A"):"R",("C","A"):"R",("Q","A"):"R",("G","A"):"C",("I","A"):"C",("L","A"):"C",("M","A"):"C",("F","A"):"C",("P","A"):"C",("S","A"):"R",("T","A"):"R",("W","A"):"C",("Y","A"):"R",("V","A"):"C",("C","N"):"C",("Q","N"):"C",("G","N"):"R",("I","N"):"R",("L","N"):"R",("M","N"):"R",("F","N"):"R",("P","N"):"R",("S","N"):"C",("T","N"):"C",("W","N"):"R",("Y","N"):"C",("V","N"):"R",("Q","C"):"C",("G","C"):"R",("I","C"):"R",("L","C"):"R",("M","C"):"R",("F","C"):"R",("P","C"):"R",("S","C"):"C",("T","C"):"C",("W","C"):"R",("Y","C"):"C",("V","C"):"R",("G","Q"):"R",("I","Q"):"R",("L","Q"):"R",("M","Q"):"R",("F","Q"):"R",("P","Q"):"R",("S","Q"):"C",("T","Q"):"C",("W","Q"):"R",("Y","Q"):"C",("V","Q"):"R",("I","G"):"C",("L","G"):"C",("M","G"):"C",("F","G"):"C",("P","G"):"C",("S","G"):"R",("T","G"):"R",("W","G"):"C",("Y","G"):"R",("V","G"):"C",("L","I"):"C",("M","I"):"C",("F","I"):"C",("P","I"):"C",("S","I"):"R",("T","I"):"R",("W","I"):"C",("Y","I"):"R",("V","I"):"C",("M","L"):"C",("F","L"):"C",("P","L"):"C",("S","L"):"R",("T","L"):"R",("W","L"):"C",("Y","L"):"R",("V","L"):"C",("F","M"):"C",("P","M"):"C",("S","M"):"R",("T","M"):"R",("W","M"):"C",("Y","M"):"R",("V","M"):"C",("P","F"):"C",("S","F"):"R",("T","F"):"R",("W","F"):"C",("Y","F"):"R",("V","F"):"C",("S","P"):"R",("T","P"):"R",("W","P"):"C",("Y","P"):"R",("V","P"):"C",("T","S"):"C",("W","S"):"R",("Y","S"):"C",("V","S"):"R",("W","T"):"R",("Y","T"):"C",("V","T"):"R",("Y","W"):"R",("V","W"):"C",("V","Y"):"R",("R","*"):"R",("H","*"):"R",("K","*"):"R",("D","*"):"R",("E","*"):"R",("A","*"):"R",("N","*"):"R",("C","*"):"R",("Q","*"):"R",("G","*"):"R",("I","*"):"R",("L","*"):"R",("M","*"):"R",("F","*"):"R",("P","*"):"R",("S","*"):"R",("T","*"):"R",("W","*"):"R",("Y","*"):"R",("V","*"):"R",("*","R"):"R",("*","H"):"R",("*","K"):"R",("*","D"):"R",("*","E"):"R",("*","A"):"R",("*","N"):"R",("*","C"):"R",("*","Q"):"R",("*","G"):"R",("*","I"):"R",("*","L"):"R",("*","M"):"R",("*","F"):"R",("*","P"):"R",("*","S"):"R",("*","T"):"R",("*","W"):"R",("*","Y"):"R",("*","V"):"R"},7:{("R","H"):"C",("R","K"):"C",("R","D"):"R",("R","E"):"R",("R","A"):"R",("R","N"):"R",("R","C"):"R",("R","Q"):"R",("R","G"):"R",("R","I"):"R",("R","L"):"R",("R","M"):"R",("R","F"):"R",("R","P"):"R",("R","S"):"R",("R","T"):"R",("R","W"):"R",("R","Y"):"R",("R","V"):"R",("H","K"):"C",("H","D"):"R",("H","E"):"R",("H","A"):"R",("H","N"):"R",("H","C"):"R",("H","Q"):"R",("H","G"):"R",("H","I"):"R",("H","L"):"R",("H","M"):"R",("H","F"):"R",("H","P"):"R",("H","S"):"R",("H","T"):"R",("H","W"):"R",("H","Y"):"R",("H","V"):"R",("K","D"):"R",("K","E"):"R",("K","A"):"R",("K","N"):"R",("K","C"):"R",("K","Q"):"R",("K","G"):"R",("K","I"):"R",("K","L"):"R",("K","M"):"R",("K","F"):"R",("K","P"):"R",("K","S"):"R",("K","T"):"R",("K","W"):"R",("K","Y"):"R",("K","V"):"R",("D","E"):"C",("D","A"):"R",("D","N"):"R",("D","C"):"R",("D","Q"):"R",("D","G"):"R",("D","I"):"R",("D","L"):"R",("D","M"):"R",("D","F"):"R",("D","P"):"R",("D","S"):"R",("D","T"):"R",("D","W"):"R",("D","Y"):"R",("D","V"):"R",("E","A"):"R",("E","N"):"R",("E","C"):"R",("E","Q"):"R",("E","G"):"R",("E","I"):"R",("E","L"):"R",("E","M"):"R",("E","F"):"R",("E","P"):"R",("E","S"):"R",("E","T"):"R",("E","W"):"R",("E","Y"):"R",("E","V"):"R",("A","N"):"R",("A","C"):"R",("A","Q"):"R",("A","G"):"R",("A","I"):"C",("A","L"):"C",("A","M"):"C",("A","F"):"C",("A","P"):"C",("A","S"):"R",("A","T"):"R",("A","W"):"C",("A","Y"):"R",("A","V"):"C",("N","C"):"C",("N","Q"):"C",("N","G"):"C",("N","I"):"R",("N","L"):"R",("N","M"):"R",("N","F"):"R",("N","P"):"R",("N","S"):"C",("N","T"):"C",("N","W"):"R",("N","Y"):"C",("N","V"):"R",("C","Q"):"C",("C","G"):"C",("C","I"):"R",("C","L"):"R",("C","M"):"R",("C","F"):"R",("C","P"):"R",("C","S"):"C",("C","T"):"C",("C","W"):"R",("C","Y"):"C",("C","V"):"R",("Q","G"):"C",("Q","I"):"R",("Q","L"):"R",("Q","M"):"R",("Q","F"):"R",("Q","P"):"R",("Q","S"):"C",("Q","T"):"C",("Q","W"):"R",("Q","Y"):"C",("Q","V"):"R",("G","I"):"R",("G","L"):"R",("G","M"):"R",("G","F"):"R",("G","P"):"R",("G","S"):"C",("G","T"):"C",("G","W"):"R",("G","Y"):"C",("G","V"):"R",("I","L"):"C",("I","M"):"C",("I","F"):"C",("I","P"):"C",("I","S"):"R",("I","T"):"R",("I","W"):"C",("I","Y"):"R",("I","V"):"C",("L","M"):"C",("L","F"):"C",("L","P"):"C",("L","S"):"R",("L","T"):"R",("L","W"):"C",("L","Y"):"R",("L","V"):"C",("M","F"):"C",("M","P"):"C",("M","S"):"R",("M","T"):"R",("M","W"):"C",("M","Y"):"R",("M","V"):"C",("F","P"):"C",("F","S"):"R",("F","T"):"R",("F","W"):"C",("F","Y"):"R",("F","V"):"C",("P","S"):"R",("P","T"):"R",("P","W"):"C",("P","Y"):"R",("P","V"):"C",("S","T"):"C",("S","W"):"R",("S","Y"):"C",("S","V"):"R",("T","W"):"R",("T","Y"):"C",("T","V"):"R",("W","Y"):"R",("W","V"):"C",("Y","V"):"R",("H","R"):"C",("K","R"):"C",("D","R"):"R",("E","R"):"R",("A","R"):"R",("N","R"):"R",("C","R"):"R",("Q","R"):"R",("G","R"):"R",("I","R"):"R",("L","R"):"R",("M","R"):"R",("F","R"):"R",("P","R"):"R",("S","R"):"R",("T","R"):"R",("W","R"):"R",("Y","R"):"R",("V","R"):"R",("K","H"):"C",("D","H"):"R",("E","H"):"R",("A","H"):"R",("N","H"):"R",("C","H"):"R",("Q","H"):"R",("G","H"):"R",("I","H"):"R",("L","H"):"R",("M","H"):"R",("F","H"):"R",("P","H"):"R",("S","H"):"R",("T","H"):"R",("W","H"):"R",("Y","H"):"R",("V","H"):"R",("D","K"):"R",("E","K"):"R",("A","K"):"R",("N","K"):"R",("C","K"):"R",("Q","K"):"R",("G","K"):"R",("I","K"):"R",("L","K"):"R",("M","K"):"R",("F","K"):"R",("P","K"):"R",("S","K"):"R",("T","K"):"R",("W","K"):"R",("Y","K"):"R",("V","K"):"R",("E","D"):"C",("A","D"):"R",("N","D"):"R",("C","D"):"R",("Q","D"):"R",("G","D"):"R",("I","D"):"R",("L","D"):"R",("M","D"):"R",("F","D"):"R",("P","D"):"R",("S","D"):"R",("T","D"):"R",("W","D"):"R",("Y","D"):"R",("V","D"):"R",("A","E"):"R",("N","E"):"R",("C","E"):"R",("Q","E"):"R",("G","E"):"R",("I","E"):"R",("L","E"):"R",("M","E"):"R",("F","E"):"R",("P","E"):"R",("S","E"):"R",("T","E"):"R",("W","E"):"R",("Y","E"):"R",("V","E"):"R",("N","A"):"R",("C","A"):"R",("Q","A"):"R",("G","A"):"R",("I","A"):"C",("L","A"):"C",("M","A"):"C",("F","A"):"C",("P","A"):"C",("S","A"):"R",("T","A"):"R",("W","A"):"C",("Y","A"):"R",("V","A"):"C",("C","N"):"C",("Q","N"):"C",("G","N"):"C",("I","N"):"R",("L","N"):"R",("M","N"):"R",("F","N"):"R",("P","N"):"R",("S","N"):"C",("T","N"):"C",("W","N"):"R",("Y","N"):"C",("V","N"):"R",("Q","C"):"C",("G","C"):"C",("I","C"):"R",("L","C"):"R",("M","C"):"R",("F","C"):"R",("P","C"):"R",("S","C"):"C",("T","C"):"C",("W","C"):"R",("Y","C"):"C",("V","C"):"R",("G","Q"):"C",("I","Q"):"R",("L","Q"):"R",("M","Q"):"R",("F","Q"):"R",("P","Q"):"R",("S","Q"):"C",("T","Q"):"C",("W","Q"):"R",("Y","Q"):"C",("V","Q"):"R",("I","G"):"R",("L","G"):"R",("M","G"):"R",("F","G"):"R",("P","G"):"R",("S","G"):"C",("T","G"):"C",("W","G"):"R",("Y","G"):"C",("V","G"):"R",("L","I"):"C",("M","I"):"C",("F","I"):"C",("P","I"):"C",("S","I"):"R",("T","I"):"R",("W","I"):"C",("Y","I"):"R",("V","I"):"C",("M","L"):"C",("F","L"):"C",("P","L"):"C",("S","L"):"R",("T","L"):"R",("W","L"):"C",("Y","L"):"R",("V","L"):"C",("F","M"):"C",("P","M"):"C",("S","M"):"R",("T","M"):"R",("W","M"):"C",("Y","M"):"R",("V","M"):"C",("P","F"):"C",("S","F"):"R",("T","F"):"R",("W","F"):"C",("Y","F"):"R",("V","F"):"C",("S","P"):"R",("T","P"):"R",("W","P"):"C",("Y","P"):"R",("V","P"):"C",("T","S"):"C",("W","S"):"R",("Y","S"):"C",("V","S"):"R",("W","T"):"R",("Y","T"):"C",("V","T"):"R",("Y","W"):"R",("V","W"):"C",("V","Y"):"R",("R","*"):"R",("H","*"):"R",("K","*"):"R",("D","*"):"R",("E","*"):"R",("A","*"):"R",("N","*"):"R",("C","*"):"R",("Q","*"):"R",("G","*"):"R",("I","*"):"R",("L","*"):"R",("M","*"):"R",("F","*"):"R",("P","*"):"R",("S","*"):"R",("T","*"):"R",("W","*"):"R",("Y","*"):"R",("V","*"):"R",("*","R"):"R",("*","H"):"R",("*","K"):"R",("*","D"):"R",("*","E"):"R",("*","A"):"R",("*","N"):"R",("*","C"):"R",("*","Q"):"R",("*","G"):"R",("*","I"):"R",("*","L"):"R",("*","M"):"R",("*","F"):"R",("*","P"):"R",("*","S"):"R",("*","T"):"R",("*","W"):"R",("*","Y"):"R",("*","V"):"R"}}
    resultsList = []
    cri = 0
    for scheme in aaSchemeList:
        currScheme = aaSchemeDict[scheme]
        currValue = currScheme[(aaList[0],aaList[1])]
        if currValue == 'R':
            cri += 1
            resultsList.append(1)
        else:
            resultsList.append(0)
    cri = cri/7.0
    resultsList.append(cri)
    return resultsList

def countSites(codonList):
    code = 'invertebrateMt'
    geneticCodes = {'standard':{"TTT":"F",	"TTC":"F",	"TTA":"L",	"TTG":"L",	"TCT":"S",	"TCC":"S",	"TCA":"S",	"TCG":"S",	"TAT":"Y",	"TAC":"Y",	"TAA":"*",	"TAG":"*",	"TGT":"C",	"TGC":"C",	"TGA":"*",	"TGG":"W",	"CTT":"L",	"CTC":"L",	"CTA":"L",	"CTG":"L",	"CCT":"P",	"CCC":"P",	"CCA":"P",	"CCG":"P",	"CAT":"H",	"CAC":"H",	"CAA":"Q",	"CAG":"Q",	"CGT":"R",	"CGC":"R",	"CGA":"R",	"CGG":"R",	"ATT":"I",	"ATC":"I",	"ATA":"I",	"ATG":"M",	"ACT":"T",	"ACC":"T",	"ACA":"T",	"ACG":"T",	"AAT":"N",	"AAC":"N",	"AAA":"K",	"AAG":"K",	"AGT":"S",	"AGC":"S",	"AGA":"R",	"AGG":"R",	"GTT":"V",	"GTC":"V",	"GTA":"V",	"GTG":"V",	"GCT":"A",	"GCC":"A",	"GCA":"A",	"GCG":"A",	"GAT":"D",	"GAC":"D",	"GAA":"E",	"GAG":"E",	"GGT":"G",	"GGC":"G",	"GGA":"G",	"GGG":"G"},'invertebrateMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'M', 'AGG': 'S', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'S', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'vertebrateMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'M', 'AGG': '*', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': '*', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'}, 'yeastMt':{'CTT': 'T', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'M', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'T', 'CTA': 'T', 'CTC': 'T', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'}, 'coelenterateMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'ciliateNuc':{'CTT': 'L', 'TAG': 'Q', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': '*', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': 'Q', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'echinodermMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'S', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'S', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'N', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'}, 'euplotidNuc':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'C', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'}, 'bacterial':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': '*', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'yeastNuc':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'S', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': '*', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'}, 'ascidianMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'M', 'AGG': 'G', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'G', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'flatwormMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'S', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'S', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'N', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': 'Y', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'chlorophyceanMt':{'CTT': 'L', 'TAG': 'L', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': '*', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'trematodeMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'M', 'AGG': 'S', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'S', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'N', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'},'pterobranchiaMt':{'CTT': 'L', 'TAG': '*', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'K', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'S', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': 'W', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TAA': '*', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'}}
    geneticCode = geneticCodes[code]
    startCodons = ['ATT','ATC','ATA','ATG','GTG'] #invertebrateMt cod
    aaSchemeDict1 = {("R","H"):"C",("R","K"):"C",("R","D"):"R",("R","E"):"R",("R","A"):"R",("R","N"):"R",("R","C"):"R",("R","Q"):"R",("R","G"):"R",("R","I"):"R",("R","L"):"R",("R","M"):"R",("R","F"):"R",("R","P"):"R",("R","S"):"R",("R","T"):"R",("R","W"):"R",("R","Y"):"R",("R","V"):"R",("H","K"):"C",("H","D"):"R",("H","E"):"R",("H","A"):"R",("H","N"):"R",("H","C"):"R",("H","Q"):"R",("H","G"):"R",("H","I"):"R",("H","L"):"R",("H","M"):"R",("H","F"):"R",("H","P"):"R",("H","S"):"R",("H","T"):"R",("H","W"):"R",("H","Y"):"R",("H","V"):"R",("K","D"):"R",("K","E"):"R",("K","A"):"R",("K","N"):"R",("K","C"):"R",("K","Q"):"R",("K","G"):"R",("K","I"):"R",("K","L"):"R",("K","M"):"R",("K","F"):"R",("K","P"):"R",("K","S"):"R",("K","T"):"R",("K","W"):"R",("K","Y"):"R",("K","V"):"R",("D","E"):"C",("D","A"):"R",("D","N"):"R",("D","C"):"R",("D","Q"):"R",("D","G"):"R",("D","I"):"R",("D","L"):"R",("D","M"):"R",("D","F"):"R",("D","P"):"R",("D","S"):"R",("D","T"):"R",("D","W"):"R",("D","Y"):"R",("D","V"):"R",("E","A"):"R",("E","N"):"R",("E","C"):"R",("E","Q"):"R",("E","G"):"R",("E","I"):"R",("E","L"):"R",("E","M"):"R",("E","F"):"R",("E","P"):"R",("E","S"):"R",("E","T"):"R",("E","W"):"R",("E","Y"):"R",("E","V"):"R",("A","N"):"C",("A","C"):"C",("A","Q"):"C",("A","G"):"C",("A","I"):"C",("A","L"):"C",("A","M"):"C",("A","F"):"C",("A","P"):"C",("A","S"):"C",("A","T"):"C",("A","W"):"C",("A","Y"):"C",("A","V"):"C",("N","C"):"C",("N","Q"):"C",("N","G"):"C",("N","I"):"C",("N","L"):"C",("N","M"):"C",("N","F"):"C",("N","P"):"C",("N","S"):"C",("N","T"):"C",("N","W"):"C",("N","Y"):"C",("N","V"):"C",("C","Q"):"C",("C","G"):"C",("C","I"):"C",("C","L"):"C",("C","M"):"C",("C","F"):"C",("C","P"):"C",("C","S"):"C",("C","T"):"C",("C","W"):"C",("C","Y"):"C",("C","V"):"C",("Q","G"):"C",("Q","I"):"C",("Q","L"):"C",("Q","M"):"C",("Q","F"):"C",("Q","P"):"C",("Q","S"):"C",("Q","T"):"C",("Q","W"):"C",("Q","Y"):"C",("Q","V"):"C",("G","I"):"C",("G","L"):"C",("G","M"):"C",("G","F"):"C",("G","P"):"C",("G","S"):"C",("G","T"):"C",("G","W"):"C",("G","Y"):"C",("G","V"):"C",("I","L"):"C",("I","M"):"C",("I","F"):"C",("I","P"):"C",("I","S"):"C",("I","T"):"C",("I","W"):"C",("I","Y"):"C",("I","V"):"C",("L","M"):"C",("L","F"):"C",("L","P"):"C",("L","S"):"C",("L","T"):"C",("L","W"):"C",("L","Y"):"C",("L","V"):"C",("M","F"):"C",("M","P"):"C",("M","S"):"C",("M","T"):"C",("M","W"):"C",("M","Y"):"C",("M","V"):"C",("F","P"):"C",("F","S"):"C",("F","T"):"C",("F","W"):"C",("F","Y"):"C",("F","V"):"C",("P","S"):"C",("P","T"):"C",("P","W"):"C",("P","Y"):"C",("P","V"):"C",("S","T"):"C",("S","W"):"C",("S","Y"):"C",("S","V"):"C",("T","W"):"C",("T","Y"):"C",("T","V"):"C",("W","Y"):"C",("W","V"):"C",("Y","V"):"C",("H","R"):"C",("K","R"):"C",("D","R"):"R",("E","R"):"R",("A","R"):"R",("N","R"):"R",("C","R"):"R",("Q","R"):"R",("G","R"):"R",("I","R"):"R",("L","R"):"R",("M","R"):"R",("F","R"):"R",("P","R"):"R",("S","R"):"R",("T","R"):"R",("W","R"):"R",("Y","R"):"R",("V","R"):"R",("K","H"):"C",("D","H"):"R",("E","H"):"R",("A","H"):"R",("N","H"):"R",("C","H"):"R",("Q","H"):"R",("G","H"):"R",("I","H"):"R",("L","H"):"R",("M","H"):"R",("F","H"):"R",("P","H"):"R",("S","H"):"R",("T","H"):"R",("W","H"):"R",("Y","H"):"R",("V","H"):"R",("D","K"):"R",("E","K"):"R",("A","K"):"R",("N","K"):"R",("C","K"):"R",("Q","K"):"R",("G","K"):"R",("I","K"):"R",("L","K"):"R",("M","K"):"R",("F","K"):"R",("P","K"):"R",("S","K"):"R",("T","K"):"R",("W","K"):"R",("Y","K"):"R",("V","K"):"R",("E","D"):"C",("A","D"):"R",("N","D"):"R",("C","D"):"R",("Q","D"):"R",("G","D"):"R",("I","D"):"R",("L","D"):"R",("M","D"):"R",("F","D"):"R",("P","D"):"R",("S","D"):"R",("T","D"):"R",("W","D"):"R",("Y","D"):"R",("V","D"):"R",("A","E"):"R",("N","E"):"R",("C","E"):"R",("Q","E"):"R",("G","E"):"R",("I","E"):"R",("L","E"):"R",("M","E"):"R",("F","E"):"R",("P","E"):"R",("S","E"):"R",("T","E"):"R",("W","E"):"R",("Y","E"):"R",("V","E"):"R",("N","A"):"C",("C","A"):"C",("Q","A"):"C",("G","A"):"C",("I","A"):"C",("L","A"):"C",("M","A"):"C",("F","A"):"C",("P","A"):"C",("S","A"):"C",("T","A"):"C",("W","A"):"C",("Y","A"):"C",("V","A"):"C",("C","N"):"C",("Q","N"):"C",("G","N"):"C",("I","N"):"C",("L","N"):"C",("M","N"):"C",("F","N"):"C",("P","N"):"C",("S","N"):"C",("T","N"):"C",("W","N"):"C",("Y","N"):"C",("V","N"):"C",("Q","C"):"C",("G","C"):"C",("I","C"):"C",("L","C"):"C",("M","C"):"C",("F","C"):"C",("P","C"):"C",("S","C"):"C",("T","C"):"C",("W","C"):"C",("Y","C"):"C",("V","C"):"C",("G","Q"):"C",("I","Q"):"C",("L","Q"):"C",("M","Q"):"C",("F","Q"):"C",("P","Q"):"C",("S","Q"):"C",("T","Q"):"C",("W","Q"):"C",("Y","Q"):"C",("V","Q"):"C",("I","G"):"C",("L","G"):"C",("M","G"):"C",("F","G"):"C",("P","G"):"C",("S","G"):"C",("T","G"):"C",("W","G"):"C",("Y","G"):"C",("V","G"):"C",("L","I"):"C",("M","I"):"C",("F","I"):"C",("P","I"):"C",("S","I"):"C",("T","I"):"C",("W","I"):"C",("Y","I"):"C",("V","I"):"C",("M","L"):"C",("F","L"):"C",("P","L"):"C",("S","L"):"C",("T","L"):"C",("W","L"):"C",("Y","L"):"C",("V","L"):"C",("F","M"):"C",("P","M"):"C",("S","M"):"C",("T","M"):"C",("W","M"):"C",("Y","M"):"C",("V","M"):"C",("P","F"):"C",("S","F"):"C",("T","F"):"C",("W","F"):"C",("Y","F"):"C",("V","F"):"C",("S","P"):"C",("T","P"):"C",("W","P"):"C",("Y","P"):"C",("V","P"):"C",("T","S"):"C",("W","S"):"C",("Y","S"):"C",("V","S"):"C",("W","T"):"C",("Y","T"):"C",("V","T"):"C",("Y","W"):"C",("V","W"):"C",("V","Y"):"C",("R","*"):"R",("H","*"):"R",("K","*"):"R",("D","*"):"R",("E","*"):"R",("A","*"):"R",("N","*"):"R",("C","*"):"R",("Q","*"):"R",("G","*"):"R",("I","*"):"R",("L","*"):"R",("M","*"):"R",("F","*"):"R",("P","*"):"R",("S","*"):"R",("T","*"):"R",("W","*"):"R",("Y","*"):"R",("V","*"):"R",("*","R"):"R",("*","H"):"R",("*","K"):"R",("*","D"):"R",("*","E"):"R",("*","A"):"R",("*","N"):"R",("*","C"):"R",("*","Q"):"R",("*","G"):"R",("*","I"):"R",("*","L"):"R",("*","M"):"R",("*","F"):"R",("*","P"):"R",("*","S"):"R",("*","T"):"R",("*","W"):"R",("*","Y"):"R",("*","V"):"R"}
    aaSchemeDict2 = {("R","H"):"C",("R","K"):"C",("R","D"):"C",("R","E"):"C",("R","A"):"R",("R","N"):"C",("R","C"):"C",("R","Q"):"C",("R","G"):"C",("R","I"):"R",("R","L"):"R",("R","M"):"R",("R","F"):"R",("R","P"):"R",("R","S"):"C",("R","T"):"C",("R","W"):"R",("R","Y"):"C",("R","V"):"R",("H","K"):"C",("H","D"):"C",("H","E"):"C",("H","A"):"R",("H","N"):"C",("H","C"):"C",("H","Q"):"C",("H","G"):"C",("H","I"):"R",("H","L"):"R",("H","M"):"R",("H","F"):"R",("H","P"):"R",("H","S"):"C",("H","T"):"C",("H","W"):"R",("H","Y"):"C",("H","V"):"R",("K","D"):"C",("K","E"):"C",("K","A"):"R",("K","N"):"C",("K","C"):"C",("K","Q"):"C",("K","G"):"C",("K","I"):"R",("K","L"):"R",("K","M"):"R",("K","F"):"R",("K","P"):"R",("K","S"):"C",("K","T"):"C",("K","W"):"R",("K","Y"):"C",("K","V"):"R",("D","E"):"C",("D","A"):"R",("D","N"):"C",("D","C"):"C",("D","Q"):"C",("D","G"):"C",("D","I"):"R",("D","L"):"R",("D","M"):"R",("D","F"):"R",("D","P"):"R",("D","S"):"C",("D","T"):"C",("D","W"):"R",("D","Y"):"C",("D","V"):"R",("E","A"):"R",("E","N"):"C",("E","C"):"C",("E","Q"):"C",("E","G"):"C",("E","I"):"R",("E","L"):"R",("E","M"):"R",("E","F"):"R",("E","P"):"R",("E","S"):"C",("E","T"):"C",("E","W"):"R",("E","Y"):"C",("E","V"):"R",("A","N"):"R",("A","C"):"R",("A","Q"):"R",("A","G"):"R",("A","I"):"C",("A","L"):"C",("A","M"):"C",("A","F"):"C",("A","P"):"C",("A","S"):"R",("A","T"):"R",("A","W"):"C",("A","Y"):"R",("A","V"):"C",("N","C"):"C",("N","Q"):"C",("N","G"):"C",("N","I"):"R",("N","L"):"R",("N","M"):"R",("N","F"):"R",("N","P"):"R",("N","S"):"C",("N","T"):"C",("N","W"):"R",("N","Y"):"C",("N","V"):"R",("C","Q"):"C",("C","G"):"C",("C","I"):"R",("C","L"):"R",("C","M"):"R",("C","F"):"R",("C","P"):"R",("C","S"):"C",("C","T"):"C",("C","W"):"R",("C","Y"):"C",("C","V"):"R",("Q","G"):"C",("Q","I"):"R",("Q","L"):"R",("Q","M"):"R",("Q","F"):"R",("Q","P"):"R",("Q","S"):"C",("Q","T"):"C",("Q","W"):"R",("Q","Y"):"C",("Q","V"):"R",("G","I"):"R",("G","L"):"R",("G","M"):"R",("G","F"):"R",("G","P"):"R",("G","S"):"C",("G","T"):"C",("G","W"):"R",("G","Y"):"C",("G","V"):"R",("I","L"):"C",("I","M"):"C",("I","F"):"C",("I","P"):"C",("I","S"):"R",("I","T"):"R",("I","W"):"C",("I","Y"):"R",("I","V"):"C",("L","M"):"C",("L","F"):"C",("L","P"):"C",("L","S"):"R",("L","T"):"R",("L","W"):"C",("L","Y"):"R",("L","V"):"C",("M","F"):"C",("M","P"):"C",("M","S"):"R",("M","T"):"R",("M","W"):"C",("M","Y"):"R",("M","V"):"C",("F","P"):"C",("F","S"):"R",("F","T"):"R",("F","W"):"C",("F","Y"):"R",("F","V"):"C",("P","S"):"R",("P","T"):"R",("P","W"):"C",("P","Y"):"R",("P","V"):"C",("S","T"):"C",("S","W"):"R",("S","Y"):"C",("S","V"):"R",("T","W"):"R",("T","Y"):"C",("T","V"):"R",("W","Y"):"R",("W","V"):"C",("Y","V"):"R",("H","R"):"C",("K","R"):"C",("D","R"):"C",("E","R"):"C",("A","R"):"R",("N","R"):"C",("C","R"):"C",("Q","R"):"C",("G","R"):"C",("I","R"):"R",("L","R"):"R",("M","R"):"R",("F","R"):"R",("P","R"):"R",("S","R"):"C",("T","R"):"C",("W","R"):"R",("Y","R"):"C",("V","R"):"R",("K","H"):"C",("D","H"):"C",("E","H"):"C",("A","H"):"R",("N","H"):"C",("C","H"):"C",("Q","H"):"C",("G","H"):"C",("I","H"):"R",("L","H"):"R",("M","H"):"R",("F","H"):"R",("P","H"):"R",("S","H"):"C",("T","H"):"C",("W","H"):"R",("Y","H"):"C",("V","H"):"R",("D","K"):"C",("E","K"):"C",("A","K"):"R",("N","K"):"C",("C","K"):"C",("Q","K"):"C",("G","K"):"C",("I","K"):"R",("L","K"):"R",("M","K"):"R",("F","K"):"R",("P","K"):"R",("S","K"):"C",("T","K"):"C",("W","K"):"R",("Y","K"):"C",("V","K"):"R",("E","D"):"C",("A","D"):"R",("N","D"):"C",("C","D"):"C",("Q","D"):"C",("G","D"):"C",("I","D"):"R",("L","D"):"R",("M","D"):"R",("F","D"):"R",("P","D"):"R",("S","D"):"C",("T","D"):"C",("W","D"):"R",("Y","D"):"C",("V","D"):"R",("A","E"):"R",("N","E"):"C",("C","E"):"C",("Q","E"):"C",("G","E"):"C",("I","E"):"R",("L","E"):"R",("M","E"):"R",("F","E"):"R",("P","E"):"R",("S","E"):"C",("T","E"):"C",("W","E"):"R",("Y","E"):"C",("V","E"):"R",("N","A"):"R",("C","A"):"R",("Q","A"):"R",("G","A"):"R",("I","A"):"C",("L","A"):"C",("M","A"):"C",("F","A"):"C",("P","A"):"C",("S","A"):"R",("T","A"):"R",("W","A"):"C",("Y","A"):"R",("V","A"):"C",("C","N"):"C",("Q","N"):"C",("G","N"):"C",("I","N"):"R",("L","N"):"R",("M","N"):"R",("F","N"):"R",("P","N"):"R",("S","N"):"C",("T","N"):"C",("W","N"):"R",("Y","N"):"C",("V","N"):"R",("Q","C"):"C",("G","C"):"C",("I","C"):"R",("L","C"):"R",("M","C"):"R",("F","C"):"R",("P","C"):"R",("S","C"):"C",("T","C"):"C",("W","C"):"R",("Y","C"):"C",("V","C"):"R",("G","Q"):"C",("I","Q"):"R",("L","Q"):"R",("M","Q"):"R",("F","Q"):"R",("P","Q"):"R",("S","Q"):"C",("T","Q"):"C",("W","Q"):"R",("Y","Q"):"C",("V","Q"):"R",("I","G"):"R",("L","G"):"R",("M","G"):"R",("F","G"):"R",("P","G"):"R",("S","G"):"C",("T","G"):"C",("W","G"):"R",("Y","G"):"C",("V","G"):"R",("L","I"):"C",("M","I"):"C",("F","I"):"C",("P","I"):"C",("S","I"):"R",("T","I"):"R",("W","I"):"C",("Y","I"):"R",("V","I"):"C",("M","L"):"C",("F","L"):"C",("P","L"):"C",("S","L"):"R",("T","L"):"R",("W","L"):"C",("Y","L"):"R",("V","L"):"C",("F","M"):"C",("P","M"):"C",("S","M"):"R",("T","M"):"R",("W","M"):"C",("Y","M"):"R",("V","M"):"C",("P","F"):"C",("S","F"):"R",("T","F"):"R",("W","F"):"C",("Y","F"):"R",("V","F"):"C",("S","P"):"R",("T","P"):"R",("W","P"):"C",("Y","P"):"R",("V","P"):"C",("T","S"):"C",("W","S"):"R",("Y","S"):"C",("V","S"):"R",("W","T"):"R",("Y","T"):"C",("V","T"):"R",("Y","W"):"R",("V","W"):"C",("V","Y"):"R",("R","*"):"R",("H","*"):"R",("K","*"):"R",("D","*"):"R",("E","*"):"R",("A","*"):"R",("N","*"):"R",("C","*"):"R",("Q","*"):"R",("G","*"):"R",("I","*"):"R",("L","*"):"R",("M","*"):"R",("F","*"):"R",("P","*"):"R",("S","*"):"R",("T","*"):"R",("W","*"):"R",("Y","*"):"R",("V","*"):"R",("*","R"):"R",("*","H"):"R",("*","K"):"R",("*","D"):"R",("*","E"):"R",("*","A"):"R",("*","N"):"R",("*","C"):"R",("*","Q"):"R",("*","G"):"R",("*","I"):"R",("*","L"):"R",("*","M"):"R",("*","F"):"R",("*","P"):"R",("*","S"):"R",("*","T"):"R",("*","W"):"R",("*","Y"):"R",("*","V"):"R"}
    aaSchemeDict3 = {("R","H"):"C",("R","K"):"C",("R","D"):"R",("R","E"):"R",("R","A"):"R",("R","N"):"R",("R","C"):"R",("R","Q"):"R",("R","G"):"R",("R","I"):"R",("R","L"):"R",("R","M"):"R",("R","F"):"R",("R","P"):"R",("R","S"):"R",("R","T"):"R",("R","W"):"R",("R","Y"):"R",("R","V"):"R",("H","K"):"C",("H","D"):"R",("H","E"):"R",("H","A"):"R",("H","N"):"R",("H","C"):"R",("H","Q"):"R",("H","G"):"R",("H","I"):"R",("H","L"):"R",("H","M"):"R",("H","F"):"R",("H","P"):"R",("H","S"):"R",("H","T"):"R",("H","W"):"R",("H","Y"):"R",("H","V"):"R",("K","D"):"R",("K","E"):"R",("K","A"):"R",("K","N"):"R",("K","C"):"R",("K","Q"):"R",("K","G"):"R",("K","I"):"R",("K","L"):"R",("K","M"):"R",("K","F"):"R",("K","P"):"R",("K","S"):"R",("K","T"):"R",("K","W"):"R",("K","Y"):"R",("K","V"):"R",("D","E"):"C",("D","A"):"R",("D","N"):"C",("D","C"):"R",("D","Q"):"C",("D","G"):"R",("D","I"):"R",("D","L"):"R",("D","M"):"R",("D","F"):"R",("D","P"):"R",("D","S"):"R",("D","T"):"R",("D","W"):"R",("D","Y"):"R",("D","V"):"R",("E","A"):"R",("E","N"):"C",("E","C"):"R",("E","Q"):"C",("E","G"):"R",("E","I"):"R",("E","L"):"R",("E","M"):"R",("E","F"):"R",("E","P"):"R",("E","S"):"R",("E","T"):"R",("E","W"):"R",("E","Y"):"R",("E","V"):"R",("A","N"):"R",("A","C"):"R",("A","Q"):"R",("A","G"):"C",("A","I"):"R",("A","L"):"R",("A","M"):"R",("A","F"):"R",("A","P"):"C",("A","S"):"C",("A","T"):"C",("A","W"):"R",("A","Y"):"R",("A","V"):"R",("N","C"):"R",("N","Q"):"C",("N","G"):"R",("N","I"):"R",("N","L"):"R",("N","M"):"R",("N","F"):"R",("N","P"):"R",("N","S"):"R",("N","T"):"R",("N","W"):"R",("N","Y"):"R",("N","V"):"R",("C","Q"):"R",("C","G"):"R",("C","I"):"R",("C","L"):"R",("C","M"):"R",("C","F"):"R",("C","P"):"R",("C","S"):"R",("C","T"):"R",("C","W"):"R",("C","Y"):"R",("C","V"):"R",("Q","G"):"R",("Q","I"):"R",("Q","L"):"R",("Q","M"):"R",("Q","F"):"R",("Q","P"):"R",("Q","S"):"R",("Q","T"):"R",("Q","W"):"R",("Q","Y"):"R",("Q","V"):"R",("G","I"):"R",("G","L"):"R",("G","M"):"R",("G","F"):"R",("G","P"):"C",("G","S"):"C",("G","T"):"C",("G","W"):"R",("G","Y"):"R",("G","V"):"R",("I","L"):"C",("I","M"):"C",("I","F"):"R",("I","P"):"R",("I","S"):"R",("I","T"):"R",("I","W"):"R",("I","Y"):"R",("I","V"):"C",("L","M"):"C",("L","F"):"R",("L","P"):"R",("L","S"):"R",("L","T"):"R",("L","W"):"R",("L","Y"):"R",("L","V"):"C",("M","F"):"R",("M","P"):"R",("M","S"):"R",("M","T"):"R",("M","W"):"R",("M","Y"):"R",("M","V"):"C",("F","P"):"R",("F","S"):"C",("F","T"):"C",("F","W"):"R",("F","Y"):"R",("F","V"):"R",("P","S"):"C",("P","T"):"C",("P","W"):"R",("P","Y"):"R",("P","V"):"R",("S","T"):"C",("S","W"):"R",("S","Y"):"R",("S","V"):"R",("T","W"):"R",("T","Y"):"R",("T","V"):"R",("W","Y"):"C",("W","V"):"R",("Y","V"):"R",("H","R"):"C",("K","R"):"C",("D","R"):"R",("E","R"):"R",("A","R"):"R",("N","R"):"R",("C","R"):"R",("Q","R"):"R",("G","R"):"R",("I","R"):"R",("L","R"):"R",("M","R"):"R",("F","R"):"R",("P","R"):"R",("S","R"):"R",("T","R"):"R",("W","R"):"R",("Y","R"):"R",("V","R"):"R",("K","H"):"C",("D","H"):"R",("E","H"):"R",("A","H"):"R",("N","H"):"R",("C","H"):"R",("Q","H"):"R",("G","H"):"R",("I","H"):"R",("L","H"):"R",("M","H"):"R",("F","H"):"R",("P","H"):"R",("S","H"):"R",("T","H"):"R",("W","H"):"R",("Y","H"):"R",("V","H"):"R",("D","K"):"R",("E","K"):"R",("A","K"):"R",("N","K"):"R",("C","K"):"R",("Q","K"):"R",("G","K"):"R",("I","K"):"R",("L","K"):"R",("M","K"):"R",("F","K"):"R",("P","K"):"R",("S","K"):"R",("T","K"):"R",("W","K"):"R",("Y","K"):"R",("V","K"):"R",("E","D"):"C",("A","D"):"R",("N","D"):"C",("C","D"):"R",("Q","D"):"C",("G","D"):"R",("I","D"):"R",("L","D"):"R",("M","D"):"R",("F","D"):"R",("P","D"):"R",("S","D"):"R",("T","D"):"R",("W","D"):"R",("Y","D"):"R",("V","D"):"R",("A","E"):"R",("N","E"):"C",("C","E"):"R",("Q","E"):"C",("G","E"):"R",("I","E"):"R",("L","E"):"R",("M","E"):"R",("F","E"):"R",("P","E"):"R",("S","E"):"R",("T","E"):"R",("W","E"):"R",("Y","E"):"R",("V","E"):"R",("N","A"):"R",("C","A"):"R",("Q","A"):"R",("G","A"):"C",("I","A"):"R",("L","A"):"R",("M","A"):"R",("F","A"):"R",("P","A"):"C",("S","A"):"C",("T","A"):"C",("W","A"):"R",("Y","A"):"R",("V","A"):"R",("C","N"):"R",("Q","N"):"C",("G","N"):"R",("I","N"):"R",("L","N"):"R",("M","N"):"R",("F","N"):"R",("P","N"):"R",("S","N"):"R",("T","N"):"R",("W","N"):"R",("Y","N"):"R",("V","N"):"R",("Q","C"):"R",("G","C"):"R",("I","C"):"R",("L","C"):"R",("M","C"):"R",("F","C"):"R",("P","C"):"R",("S","C"):"R",("T","C"):"R",("W","C"):"R",("Y","C"):"R",("V","C"):"R",("G","Q"):"R",("I","Q"):"R",("L","Q"):"R",("M","Q"):"R",("F","Q"):"R",("P","Q"):"R",("S","Q"):"R",("T","Q"):"R",("W","Q"):"R",("Y","Q"):"R",("V","Q"):"R",("I","G"):"R",("L","G"):"R",("M","G"):"R",("F","G"):"R",("P","G"):"C",("S","G"):"C",("T","G"):"C",("W","G"):"R",("Y","G"):"R",("V","G"):"R",("L","I"):"C",("M","I"):"C",("F","I"):"R",("P","I"):"R",("S","I"):"R",("T","I"):"R",("W","I"):"R",("Y","I"):"R",("V","I"):"C",("M","L"):"C",("F","L"):"R",("P","L"):"R",("S","L"):"R",("T","L"):"R",("W","L"):"R",("Y","L"):"R",("V","L"):"C",("F","M"):"R",("P","M"):"R",("S","M"):"R",("T","M"):"R",("W","M"):"R",("Y","M"):"R",("V","M"):"C",("P","F"):"R",("S","F"):"C",("T","F"):"C",("W","F"):"R",("Y","F"):"R",("V","F"):"R",("S","P"):"C",("T","P"):"C",("W","P"):"R",("Y","P"):"R",("V","P"):"R",("T","S"):"C",("W","S"):"R",("Y","S"):"R",("V","S"):"R",("W","T"):"R",("Y","T"):"R",("V","T"):"R",("Y","W"):"C",("V","W"):"R",("V","Y"):"R",("R","*"):"R",("H","*"):"R",("K","*"):"R",("D","*"):"R",("E","*"):"R",("A","*"):"R",("N","*"):"R",("C","*"):"R",("Q","*"):"R",("G","*"):"R",("I","*"):"R",("L","*"):"R",("M","*"):"R",("F","*"):"R",("P","*"):"R",("S","*"):"R",("T","*"):"R",("W","*"):"R",("Y","*"):"R",("V","*"):"R",("*","R"):"R",("*","H"):"R",("*","K"):"R",("*","D"):"R",("*","E"):"R",("*","A"):"R",("*","N"):"R",("*","C"):"R",("*","Q"):"R",("*","G"):"R",("*","I"):"R",("*","L"):"R",("*","M"):"R",("*","F"):"R",("*","P"):"R",("*","S"):"R",("*","T"):"R",("*","W"):"R",("*","Y"):"R",("*","V"):"R"}
    aaSchemeDict4 = {("R","H"):"C",("R","K"):"C",("R","D"):"R",("R","E"):"R",("R","A"):"R",("R","N"):"R",("R","C"):"R",("R","Q"):"C",("R","G"):"R",("R","I"):"R",("R","L"):"R",("R","M"):"R",("R","F"):"C",("R","P"):"R",("R","S"):"R",("R","T"):"R",("R","W"):"C",("R","Y"):"C",("R","V"):"R",("H","K"):"C",("H","D"):"R",("H","E"):"R",("H","A"):"R",("H","N"):"R",("H","C"):"R",("H","Q"):"C",("H","G"):"R",("H","I"):"R",("H","L"):"R",("H","M"):"R",("H","F"):"C",("H","P"):"R",("H","S"):"R",("H","T"):"R",("H","W"):"C",("H","Y"):"C",("H","V"):"R",("K","D"):"R",("K","E"):"R",("K","A"):"R",("K","N"):"R",("K","C"):"R",("K","Q"):"C",("K","G"):"R",("K","I"):"R",("K","L"):"R",("K","M"):"R",("K","F"):"C",("K","P"):"R",("K","S"):"R",("K","T"):"R",("K","W"):"C",("K","Y"):"C",("K","V"):"R",("D","E"):"C",("D","A"):"R",("D","N"):"R",("D","C"):"R",("D","Q"):"R",("D","G"):"R",("D","I"):"R",("D","L"):"R",("D","M"):"R",("D","F"):"R",("D","P"):"R",("D","S"):"R",("D","T"):"R",("D","W"):"R",("D","Y"):"R",("D","V"):"R",("E","A"):"R",("E","N"):"R",("E","C"):"R",("E","Q"):"R",("E","G"):"R",("E","I"):"R",("E","L"):"R",("E","M"):"R",("E","F"):"R",("E","P"):"R",("E","S"):"R",("E","T"):"R",("E","W"):"R",("E","Y"):"R",("E","V"):"R",("A","N"):"C",("A","C"):"C",("A","Q"):"R",("A","G"):"C",("A","I"):"R",("A","L"):"R",("A","M"):"R",("A","F"):"R",("A","P"):"C",("A","S"):"C",("A","T"):"C",("A","W"):"R",("A","Y"):"R",("A","V"):"R",("N","C"):"C",("N","Q"):"R",("N","G"):"C",("N","I"):"R",("N","L"):"R",("N","M"):"R",("N","F"):"R",("N","P"):"C",("N","S"):"C",("N","T"):"C",("N","W"):"R",("N","Y"):"R",("N","V"):"R",("C","Q"):"R",("C","G"):"C",("C","I"):"R",("C","L"):"R",("C","M"):"R",("C","F"):"R",("C","P"):"C",("C","S"):"C",("C","T"):"C",("C","W"):"R",("C","Y"):"R",("C","V"):"R",("Q","G"):"R",("Q","I"):"R",("Q","L"):"R",("Q","M"):"R",("Q","F"):"C",("Q","P"):"R",("Q","S"):"R",("Q","T"):"R",("Q","W"):"C",("Q","Y"):"C",("Q","V"):"R",("G","I"):"R",("G","L"):"R",("G","M"):"R",("G","F"):"R",("G","P"):"C",("G","S"):"C",("G","T"):"C",("G","W"):"R",("G","Y"):"R",("G","V"):"R",("I","L"):"C",("I","M"):"C",("I","F"):"R",("I","P"):"R",("I","S"):"R",("I","T"):"R",("I","W"):"R",("I","Y"):"R",("I","V"):"C",("L","M"):"C",("L","F"):"R",("L","P"):"R",("L","S"):"R",("L","T"):"R",("L","W"):"R",("L","Y"):"R",("L","V"):"C",("M","F"):"R",("M","P"):"R",("M","S"):"R",("M","T"):"R",("M","W"):"R",("M","Y"):"R",("M","V"):"C",("F","P"):"R",("F","S"):"R",("F","T"):"R",("F","W"):"C",("F","Y"):"C",("F","V"):"R",("P","S"):"C",("P","T"):"C",("P","W"):"R",("P","Y"):"R",("P","V"):"R",("S","T"):"C",("S","W"):"R",("S","Y"):"R",("S","V"):"R",("T","W"):"R",("T","Y"):"R",("T","V"):"R",("W","Y"):"C",("W","V"):"R",("Y","V"):"R",("H","R"):"C",("K","R"):"C",("D","R"):"R",("E","R"):"R",("A","R"):"R",("N","R"):"R",("C","R"):"R",("Q","R"):"C",("G","R"):"R",("I","R"):"R",("L","R"):"R",("M","R"):"R",("F","R"):"C",("P","R"):"R",("S","R"):"R",("T","R"):"R",("W","R"):"C",("Y","R"):"C",("V","R"):"R",("K","H"):"C",("D","H"):"R",("E","H"):"R",("A","H"):"R",("N","H"):"R",("C","H"):"R",("Q","H"):"C",("G","H"):"R",("I","H"):"R",("L","H"):"R",("M","H"):"R",("F","H"):"C",("P","H"):"R",("S","H"):"R",("T","H"):"R",("W","H"):"C",("Y","H"):"C",("V","H"):"R",("D","K"):"R",("E","K"):"R",("A","K"):"R",("N","K"):"R",("C","K"):"R",("Q","K"):"C",("G","K"):"R",("I","K"):"R",("L","K"):"R",("M","K"):"R",("F","K"):"C",("P","K"):"R",("S","K"):"R",("T","K"):"R",("W","K"):"C",("Y","K"):"C",("V","K"):"R",("E","D"):"C",("A","D"):"R",("N","D"):"R",("C","D"):"R",("Q","D"):"R",("G","D"):"R",("I","D"):"R",("L","D"):"R",("M","D"):"R",("F","D"):"R",("P","D"):"R",("S","D"):"R",("T","D"):"R",("W","D"):"R",("Y","D"):"R",("V","D"):"R",("A","E"):"R",("N","E"):"R",("C","E"):"R",("Q","E"):"R",("G","E"):"R",("I","E"):"R",("L","E"):"R",("M","E"):"R",("F","E"):"R",("P","E"):"R",("S","E"):"R",("T","E"):"R",("W","E"):"R",("Y","E"):"R",("V","E"):"R",("N","A"):"C",("C","A"):"C",("Q","A"):"R",("G","A"):"C",("I","A"):"R",("L","A"):"R",("M","A"):"R",("F","A"):"R",("P","A"):"C",("S","A"):"C",("T","A"):"C",("W","A"):"R",("Y","A"):"R",("V","A"):"R",("C","N"):"C",("Q","N"):"R",("G","N"):"C",("I","N"):"R",("L","N"):"R",("M","N"):"R",("F","N"):"R",("P","N"):"C",("S","N"):"C",("T","N"):"C",("W","N"):"R",("Y","N"):"R",("V","N"):"R",("Q","C"):"R",("G","C"):"C",("I","C"):"R",("L","C"):"R",("M","C"):"R",("F","C"):"R",("P","C"):"C",("S","C"):"C",("T","C"):"C",("W","C"):"R",("Y","C"):"R",("V","C"):"R",("G","Q"):"R",("I","Q"):"R",("L","Q"):"R",("M","Q"):"R",("F","Q"):"C",("P","Q"):"R",("S","Q"):"R",("T","Q"):"R",("W","Q"):"C",("Y","Q"):"C",("V","Q"):"R",("I","G"):"R",("L","G"):"R",("M","G"):"R",("F","G"):"R",("P","G"):"C",("S","G"):"C",("T","G"):"C",("W","G"):"R",("Y","G"):"R",("V","G"):"R",("L","I"):"C",("M","I"):"C",("F","I"):"R",("P","I"):"R",("S","I"):"R",("T","I"):"R",("W","I"):"R",("Y","I"):"R",("V","I"):"C",("M","L"):"C",("F","L"):"R",("P","L"):"R",("S","L"):"R",("T","L"):"R",("W","L"):"R",("Y","L"):"R",("V","L"):"C",("F","M"):"R",("P","M"):"R",("S","M"):"R",("T","M"):"R",("W","M"):"R",("Y","M"):"R",("V","M"):"C",("P","F"):"R",("S","F"):"R",("T","F"):"R",("W","F"):"C",("Y","F"):"C",("V","F"):"R",("S","P"):"C",("T","P"):"C",("W","P"):"R",("Y","P"):"R",("V","P"):"R",("T","S"):"C",("W","S"):"R",("Y","S"):"R",("V","S"):"R",("W","T"):"R",("Y","T"):"R",("V","T"):"R",("Y","W"):"C",("V","W"):"R",("V","Y"):"R",("R","*"):"R",("H","*"):"R",("K","*"):"R",("D","*"):"R",("E","*"):"R",("A","*"):"R",("N","*"):"R",("C","*"):"R",("Q","*"):"R",("G","*"):"R",("I","*"):"R",("L","*"):"R",("M","*"):"R",("F","*"):"R",("P","*"):"R",("S","*"):"R",("T","*"):"R",("W","*"):"R",("Y","*"):"R",("V","*"):"R",("*","R"):"R",("*","H"):"R",("*","K"):"R",("*","D"):"R",("*","E"):"R",("*","A"):"R",("*","N"):"R",("*","C"):"R",("*","Q"):"R",("*","G"):"R",("*","I"):"R",("*","L"):"R",("*","M"):"R",("*","F"):"R",("*","P"):"R",("*","S"):"R",("*","T"):"R",("*","W"):"R",("*","Y"):"R",("*","V"):"R"}
    aaSchemeDict5 = {("R","H"):"C",("R","K"):"C",("R","D"):"R",("R","E"):"R",("R","A"):"R",("R","N"):"R",("R","C"):"R",("R","Q"):"R",("R","G"):"R",("R","I"):"R",("R","L"):"R",("R","M"):"R",("R","F"):"R",("R","P"):"R",("R","S"):"R",("R","T"):"R",("R","W"):"R",("R","Y"):"R",("R","V"):"R",("H","K"):"C",("H","D"):"R",("H","E"):"R",("H","A"):"R",("H","N"):"R",("H","C"):"R",("H","Q"):"R",("H","G"):"R",("H","I"):"R",("H","L"):"R",("H","M"):"R",("H","F"):"R",("H","P"):"R",("H","S"):"R",("H","T"):"R",("H","W"):"R",("H","Y"):"R",("H","V"):"R",("K","D"):"R",("K","E"):"R",("K","A"):"R",("K","N"):"R",("K","C"):"R",("K","Q"):"R",("K","G"):"R",("K","I"):"R",("K","L"):"R",("K","M"):"R",("K","F"):"R",("K","P"):"R",("K","S"):"R",("K","T"):"R",("K","W"):"R",("K","Y"):"R",("K","V"):"R",("D","E"):"C",("D","A"):"R",("D","N"):"R",("D","C"):"R",("D","Q"):"R",("D","G"):"R",("D","I"):"R",("D","L"):"R",("D","M"):"R",("D","F"):"R",("D","P"):"R",("D","S"):"R",("D","T"):"R",("D","W"):"R",("D","Y"):"R",("D","V"):"R",("E","A"):"R",("E","N"):"R",("E","C"):"R",("E","Q"):"R",("E","G"):"R",("E","I"):"R",("E","L"):"R",("E","M"):"R",("E","F"):"R",("E","P"):"R",("E","S"):"R",("E","T"):"R",("E","W"):"R",("E","Y"):"R",("E","V"):"R",("A","N"):"C",("A","C"):"C",("A","Q"):"C",("A","G"):"C",("A","I"):"C",("A","L"):"C",("A","M"):"C",("A","F"):"R",("A","P"):"C",("A","S"):"C",("A","T"):"C",("A","W"):"R",("A","Y"):"R",("A","V"):"C",("N","C"):"C",("N","Q"):"C",("N","G"):"C",("N","I"):"C",("N","L"):"C",("N","M"):"C",("N","F"):"R",("N","P"):"C",("N","S"):"C",("N","T"):"C",("N","W"):"R",("N","Y"):"R",("N","V"):"C",("C","Q"):"C",("C","G"):"C",("C","I"):"C",("C","L"):"C",("C","M"):"C",("C","F"):"R",("C","P"):"C",("C","S"):"C",("C","T"):"C",("C","W"):"R",("C","Y"):"R",("C","V"):"C",("Q","G"):"C",("Q","I"):"C",("Q","L"):"C",("Q","M"):"C",("Q","F"):"R",("Q","P"):"C",("Q","S"):"C",("Q","T"):"C",("Q","W"):"R",("Q","Y"):"R",("Q","V"):"C",("G","I"):"C",("G","L"):"C",("G","M"):"C",("G","F"):"R",("G","P"):"C",("G","S"):"C",("G","T"):"C",("G","W"):"R",("G","Y"):"R",("G","V"):"C",("I","L"):"C",("I","M"):"C",("I","F"):"R",("I","P"):"C",("I","S"):"C",("I","T"):"C",("I","W"):"R",("I","Y"):"R",("I","V"):"C",("L","M"):"C",("L","F"):"R",("L","P"):"C",("L","S"):"C",("L","T"):"C",("L","W"):"R",("L","Y"):"R",("L","V"):"C",("M","F"):"R",("M","P"):"C",("M","S"):"C",("M","T"):"C",("M","W"):"R",("M","Y"):"R",("M","V"):"C",("F","P"):"R",("F","S"):"R",("F","T"):"R",("F","W"):"C",("F","Y"):"C",("F","V"):"R",("P","S"):"C",("P","T"):"C",("P","W"):"R",("P","Y"):"R",("P","V"):"C",("S","T"):"C",("S","W"):"R",("S","Y"):"R",("S","V"):"R",("T","W"):"R",("T","Y"):"R",("T","V"):"C",("W","Y"):"C",("W","V"):"R",("Y","V"):"R",("H","R"):"C",("K","R"):"C",("D","R"):"R",("E","R"):"R",("A","R"):"R",("N","R"):"R",("C","R"):"R",("Q","R"):"R",("G","R"):"R",("I","R"):"R",("L","R"):"R",("M","R"):"R",("F","R"):"R",("P","R"):"R",("S","R"):"R",("T","R"):"R",("W","R"):"R",("Y","R"):"R",("V","R"):"R",("K","H"):"C",("D","H"):"R",("E","H"):"R",("A","H"):"R",("N","H"):"R",("C","H"):"R",("Q","H"):"R",("G","H"):"R",("I","H"):"R",("L","H"):"R",("M","H"):"R",("F","H"):"R",("P","H"):"R",("S","H"):"R",("T","H"):"R",("W","H"):"R",("Y","H"):"R",("V","H"):"R",("D","K"):"R",("E","K"):"R",("A","K"):"R",("N","K"):"R",("C","K"):"R",("Q","K"):"R",("G","K"):"R",("I","K"):"R",("L","K"):"R",("M","K"):"R",("F","K"):"R",("P","K"):"R",("S","K"):"R",("T","K"):"R",("W","K"):"R",("Y","K"):"R",("V","K"):"R",("E","D"):"C",("A","D"):"R",("N","D"):"R",("C","D"):"R",("Q","D"):"R",("G","D"):"R",("I","D"):"R",("L","D"):"R",("M","D"):"R",("F","D"):"R",("P","D"):"R",("S","D"):"R",("T","D"):"R",("W","D"):"R",("Y","D"):"R",("V","D"):"R",("A","E"):"R",("N","E"):"R",("C","E"):"R",("Q","E"):"R",("G","E"):"R",("I","E"):"R",("L","E"):"R",("M","E"):"R",("F","E"):"R",("P","E"):"R",("S","E"):"R",("T","E"):"R",("W","E"):"R",("Y","E"):"R",("V","E"):"R",("N","A"):"C",("C","A"):"C",("Q","A"):"C",("G","A"):"C",("I","A"):"C",("L","A"):"C",("M","A"):"C",("F","A"):"R",("P","A"):"C",("S","A"):"C",("T","A"):"C",("W","A"):"R",("Y","A"):"R",("V","A"):"C",("C","N"):"C",("Q","N"):"C",("G","N"):"C",("I","N"):"C",("L","N"):"C",("M","N"):"C",("F","N"):"R",("P","N"):"C",("S","N"):"C",("T","N"):"C",("W","N"):"R",("Y","N"):"R",("V","N"):"C",("Q","C"):"C",("G","C"):"C",("I","C"):"C",("L","C"):"C",("M","C"):"C",("F","C"):"R",("P","C"):"C",("S","C"):"C",("T","C"):"C",("W","C"):"R",("Y","C"):"R",("V","C"):"C",("G","Q"):"C",("I","Q"):"C",("L","Q"):"C",("M","Q"):"C",("F","Q"):"R",("P","Q"):"C",("S","Q"):"C",("T","Q"):"C",("W","Q"):"R",("Y","Q"):"R",("V","Q"):"C",("I","G"):"C",("L","G"):"C",("M","G"):"C",("F","G"):"R",("P","G"):"C",("S","G"):"C",("T","G"):"C",("W","G"):"R",("Y","G"):"R",("V","G"):"C",("L","I"):"C",("M","I"):"C",("F","I"):"R",("P","I"):"C",("S","I"):"C",("T","I"):"C",("W","I"):"R",("Y","I"):"R",("V","I"):"C",("M","L"):"C",("F","L"):"R",("P","L"):"C",("S","L"):"C",("T","L"):"C",("W","L"):"R",("Y","L"):"R",("V","L"):"C",("F","M"):"R",("P","M"):"C",("S","M"):"C",("T","M"):"C",("W","M"):"R",("Y","M"):"R",("V","M"):"C",("P","F"):"R",("S","F"):"R",("T","F"):"R",("W","F"):"C",("Y","F"):"C",("V","F"):"R",("S","P"):"C",("T","P"):"C",("W","P"):"R",("Y","P"):"R",("V","P"):"C",("T","S"):"C",("W","S"):"R",("Y","S"):"R",("V","S"):"R",("W","T"):"R",("Y","T"):"R",("V","T"):"C",("Y","W"):"C",("V","W"):"R",("V","Y"):"R",("R","*"):"R",("H","*"):"R",("K","*"):"R",("D","*"):"R",("E","*"):"R",("A","*"):"R",("N","*"):"R",("C","*"):"R",("Q","*"):"R",("G","*"):"R",("I","*"):"R",("L","*"):"R",("M","*"):"R",("F","*"):"R",("P","*"):"R",("S","*"):"R",("T","*"):"R",("W","*"):"R",("Y","*"):"R",("V","*"):"R",("*","R"):"R",("*","H"):"R",("*","K"):"R",("*","D"):"R",("*","E"):"R",("*","A"):"R",("*","N"):"R",("*","C"):"R",("*","Q"):"R",("*","G"):"R",("*","I"):"R",("*","L"):"R",("*","M"):"R",("*","F"):"R",("*","P"):"R",("*","S"):"R",("*","T"):"R",("*","W"):"R",("*","Y"):"R",("*","V"):"R"}
    aaSchemeDict6 = {("R","H"):"C",("R","K"):"C",("R","D"):"R",("R","E"):"R",("R","A"):"R",("R","N"):"R",("R","C"):"R",("R","Q"):"R",("R","G"):"R",("R","I"):"R",("R","L"):"R",("R","M"):"R",("R","F"):"R",("R","P"):"R",("R","S"):"R",("R","T"):"R",("R","W"):"R",("R","Y"):"R",("R","V"):"R",("H","K"):"C",("H","D"):"R",("H","E"):"R",("H","A"):"R",("H","N"):"R",("H","C"):"R",("H","Q"):"R",("H","G"):"R",("H","I"):"R",("H","L"):"R",("H","M"):"R",("H","F"):"R",("H","P"):"R",("H","S"):"R",("H","T"):"R",("H","W"):"R",("H","Y"):"R",("H","V"):"R",("K","D"):"R",("K","E"):"R",("K","A"):"R",("K","N"):"R",("K","C"):"R",("K","Q"):"R",("K","G"):"R",("K","I"):"R",("K","L"):"R",("K","M"):"R",("K","F"):"R",("K","P"):"R",("K","S"):"R",("K","T"):"R",("K","W"):"R",("K","Y"):"R",("K","V"):"R",("D","E"):"C",("D","A"):"R",("D","N"):"R",("D","C"):"R",("D","Q"):"R",("D","G"):"R",("D","I"):"R",("D","L"):"R",("D","M"):"R",("D","F"):"R",("D","P"):"R",("D","S"):"R",("D","T"):"R",("D","W"):"R",("D","Y"):"R",("D","V"):"R",("E","A"):"R",("E","N"):"R",("E","C"):"R",("E","Q"):"R",("E","G"):"R",("E","I"):"R",("E","L"):"R",("E","M"):"R",("E","F"):"R",("E","P"):"R",("E","S"):"R",("E","T"):"R",("E","W"):"R",("E","Y"):"R",("E","V"):"R",("A","N"):"R",("A","C"):"R",("A","Q"):"R",("A","G"):"C",("A","I"):"C",("A","L"):"C",("A","M"):"C",("A","F"):"C",("A","P"):"C",("A","S"):"R",("A","T"):"R",("A","W"):"C",("A","Y"):"R",("A","V"):"C",("N","C"):"C",("N","Q"):"C",("N","G"):"R",("N","I"):"R",("N","L"):"R",("N","M"):"R",("N","F"):"R",("N","P"):"R",("N","S"):"C",("N","T"):"C",("N","W"):"R",("N","Y"):"C",("N","V"):"R",("C","Q"):"C",("C","G"):"R",("C","I"):"R",("C","L"):"R",("C","M"):"R",("C","F"):"R",("C","P"):"R",("C","S"):"C",("C","T"):"C",("C","W"):"R",("C","Y"):"C",("C","V"):"R",("Q","G"):"R",("Q","I"):"R",("Q","L"):"R",("Q","M"):"R",("Q","F"):"R",("Q","P"):"R",("Q","S"):"C",("Q","T"):"C",("Q","W"):"R",("Q","Y"):"C",("Q","V"):"R",("G","I"):"C",("G","L"):"C",("G","M"):"C",("G","F"):"C",("G","P"):"C",("G","S"):"R",("G","T"):"R",("G","W"):"C",("G","Y"):"R",("G","V"):"C",("I","L"):"C",("I","M"):"C",("I","F"):"C",("I","P"):"C",("I","S"):"R",("I","T"):"R",("I","W"):"C",("I","Y"):"R",("I","V"):"C",("L","M"):"C",("L","F"):"C",("L","P"):"C",("L","S"):"R",("L","T"):"R",("L","W"):"C",("L","Y"):"R",("L","V"):"C",("M","F"):"C",("M","P"):"C",("M","S"):"R",("M","T"):"R",("M","W"):"C",("M","Y"):"R",("M","V"):"C",("F","P"):"C",("F","S"):"R",("F","T"):"R",("F","W"):"C",("F","Y"):"R",("F","V"):"C",("P","S"):"R",("P","T"):"R",("P","W"):"C",("P","Y"):"R",("P","V"):"C",("S","T"):"C",("S","W"):"R",("S","Y"):"C",("S","V"):"R",("T","W"):"R",("T","Y"):"C",("T","V"):"R",("W","Y"):"R",("W","V"):"C",("Y","V"):"R",("H","R"):"C",("K","R"):"C",("D","R"):"R",("E","R"):"R",("A","R"):"R",("N","R"):"R",("C","R"):"R",("Q","R"):"R",("G","R"):"R",("I","R"):"R",("L","R"):"R",("M","R"):"R",("F","R"):"R",("P","R"):"R",("S","R"):"R",("T","R"):"R",("W","R"):"R",("Y","R"):"R",("V","R"):"R",("K","H"):"C",("D","H"):"R",("E","H"):"R",("A","H"):"R",("N","H"):"R",("C","H"):"R",("Q","H"):"R",("G","H"):"R",("I","H"):"R",("L","H"):"R",("M","H"):"R",("F","H"):"R",("P","H"):"R",("S","H"):"R",("T","H"):"R",("W","H"):"R",("Y","H"):"R",("V","H"):"R",("D","K"):"R",("E","K"):"R",("A","K"):"R",("N","K"):"R",("C","K"):"R",("Q","K"):"R",("G","K"):"R",("I","K"):"R",("L","K"):"R",("M","K"):"R",("F","K"):"R",("P","K"):"R",("S","K"):"R",("T","K"):"R",("W","K"):"R",("Y","K"):"R",("V","K"):"R",("E","D"):"C",("A","D"):"R",("N","D"):"R",("C","D"):"R",("Q","D"):"R",("G","D"):"R",("I","D"):"R",("L","D"):"R",("M","D"):"R",("F","D"):"R",("P","D"):"R",("S","D"):"R",("T","D"):"R",("W","D"):"R",("Y","D"):"R",("V","D"):"R",("A","E"):"R",("N","E"):"R",("C","E"):"R",("Q","E"):"R",("G","E"):"R",("I","E"):"R",("L","E"):"R",("M","E"):"R",("F","E"):"R",("P","E"):"R",("S","E"):"R",("T","E"):"R",("W","E"):"R",("Y","E"):"R",("V","E"):"R",("N","A"):"R",("C","A"):"R",("Q","A"):"R",("G","A"):"C",("I","A"):"C",("L","A"):"C",("M","A"):"C",("F","A"):"C",("P","A"):"C",("S","A"):"R",("T","A"):"R",("W","A"):"C",("Y","A"):"R",("V","A"):"C",("C","N"):"C",("Q","N"):"C",("G","N"):"R",("I","N"):"R",("L","N"):"R",("M","N"):"R",("F","N"):"R",("P","N"):"R",("S","N"):"C",("T","N"):"C",("W","N"):"R",("Y","N"):"C",("V","N"):"R",("Q","C"):"C",("G","C"):"R",("I","C"):"R",("L","C"):"R",("M","C"):"R",("F","C"):"R",("P","C"):"R",("S","C"):"C",("T","C"):"C",("W","C"):"R",("Y","C"):"C",("V","C"):"R",("G","Q"):"R",("I","Q"):"R",("L","Q"):"R",("M","Q"):"R",("F","Q"):"R",("P","Q"):"R",("S","Q"):"C",("T","Q"):"C",("W","Q"):"R",("Y","Q"):"C",("V","Q"):"R",("I","G"):"C",("L","G"):"C",("M","G"):"C",("F","G"):"C",("P","G"):"C",("S","G"):"R",("T","G"):"R",("W","G"):"C",("Y","G"):"R",("V","G"):"C",("L","I"):"C",("M","I"):"C",("F","I"):"C",("P","I"):"C",("S","I"):"R",("T","I"):"R",("W","I"):"C",("Y","I"):"R",("V","I"):"C",("M","L"):"C",("F","L"):"C",("P","L"):"C",("S","L"):"R",("T","L"):"R",("W","L"):"C",("Y","L"):"R",("V","L"):"C",("F","M"):"C",("P","M"):"C",("S","M"):"R",("T","M"):"R",("W","M"):"C",("Y","M"):"R",("V","M"):"C",("P","F"):"C",("S","F"):"R",("T","F"):"R",("W","F"):"C",("Y","F"):"R",("V","F"):"C",("S","P"):"R",("T","P"):"R",("W","P"):"C",("Y","P"):"R",("V","P"):"C",("T","S"):"C",("W","S"):"R",("Y","S"):"C",("V","S"):"R",("W","T"):"R",("Y","T"):"C",("V","T"):"R",("Y","W"):"R",("V","W"):"C",("V","Y"):"R",("R","*"):"R",("H","*"):"R",("K","*"):"R",("D","*"):"R",("E","*"):"R",("A","*"):"R",("N","*"):"R",("C","*"):"R",("Q","*"):"R",("G","*"):"R",("I","*"):"R",("L","*"):"R",("M","*"):"R",("F","*"):"R",("P","*"):"R",("S","*"):"R",("T","*"):"R",("W","*"):"R",("Y","*"):"R",("V","*"):"R",("*","R"):"R",("*","H"):"R",("*","K"):"R",("*","D"):"R",("*","E"):"R",("*","A"):"R",("*","N"):"R",("*","C"):"R",("*","Q"):"R",("*","G"):"R",("*","I"):"R",("*","L"):"R",("*","M"):"R",("*","F"):"R",("*","P"):"R",("*","S"):"R",("*","T"):"R",("*","W"):"R",("*","Y"):"R",("*","V"):"R"}
    aaSchemeDict7 = {("R","H"):"C",("R","K"):"C",("R","D"):"R",("R","E"):"R",("R","A"):"R",("R","N"):"R",("R","C"):"R",("R","Q"):"R",("R","G"):"R",("R","I"):"R",("R","L"):"R",("R","M"):"R",("R","F"):"R",("R","P"):"R",("R","S"):"R",("R","T"):"R",("R","W"):"R",("R","Y"):"R",("R","V"):"R",("H","K"):"C",("H","D"):"R",("H","E"):"R",("H","A"):"R",("H","N"):"R",("H","C"):"R",("H","Q"):"R",("H","G"):"R",("H","I"):"R",("H","L"):"R",("H","M"):"R",("H","F"):"R",("H","P"):"R",("H","S"):"R",("H","T"):"R",("H","W"):"R",("H","Y"):"R",("H","V"):"R",("K","D"):"R",("K","E"):"R",("K","A"):"R",("K","N"):"R",("K","C"):"R",("K","Q"):"R",("K","G"):"R",("K","I"):"R",("K","L"):"R",("K","M"):"R",("K","F"):"R",("K","P"):"R",("K","S"):"R",("K","T"):"R",("K","W"):"R",("K","Y"):"R",("K","V"):"R",("D","E"):"C",("D","A"):"R",("D","N"):"R",("D","C"):"R",("D","Q"):"R",("D","G"):"R",("D","I"):"R",("D","L"):"R",("D","M"):"R",("D","F"):"R",("D","P"):"R",("D","S"):"R",("D","T"):"R",("D","W"):"R",("D","Y"):"R",("D","V"):"R",("E","A"):"R",("E","N"):"R",("E","C"):"R",("E","Q"):"R",("E","G"):"R",("E","I"):"R",("E","L"):"R",("E","M"):"R",("E","F"):"R",("E","P"):"R",("E","S"):"R",("E","T"):"R",("E","W"):"R",("E","Y"):"R",("E","V"):"R",("A","N"):"R",("A","C"):"R",("A","Q"):"R",("A","G"):"R",("A","I"):"C",("A","L"):"C",("A","M"):"C",("A","F"):"C",("A","P"):"C",("A","S"):"R",("A","T"):"R",("A","W"):"C",("A","Y"):"R",("A","V"):"C",("N","C"):"C",("N","Q"):"C",("N","G"):"C",("N","I"):"R",("N","L"):"R",("N","M"):"R",("N","F"):"R",("N","P"):"R",("N","S"):"C",("N","T"):"C",("N","W"):"R",("N","Y"):"C",("N","V"):"R",("C","Q"):"C",("C","G"):"C",("C","I"):"R",("C","L"):"R",("C","M"):"R",("C","F"):"R",("C","P"):"R",("C","S"):"C",("C","T"):"C",("C","W"):"R",("C","Y"):"C",("C","V"):"R",("Q","G"):"C",("Q","I"):"R",("Q","L"):"R",("Q","M"):"R",("Q","F"):"R",("Q","P"):"R",("Q","S"):"C",("Q","T"):"C",("Q","W"):"R",("Q","Y"):"C",("Q","V"):"R",("G","I"):"R",("G","L"):"R",("G","M"):"R",("G","F"):"R",("G","P"):"R",("G","S"):"C",("G","T"):"C",("G","W"):"R",("G","Y"):"C",("G","V"):"R",("I","L"):"C",("I","M"):"C",("I","F"):"C",("I","P"):"C",("I","S"):"R",("I","T"):"R",("I","W"):"C",("I","Y"):"R",("I","V"):"C",("L","M"):"C",("L","F"):"C",("L","P"):"C",("L","S"):"R",("L","T"):"R",("L","W"):"C",("L","Y"):"R",("L","V"):"C",("M","F"):"C",("M","P"):"C",("M","S"):"R",("M","T"):"R",("M","W"):"C",("M","Y"):"R",("M","V"):"C",("F","P"):"C",("F","S"):"R",("F","T"):"R",("F","W"):"C",("F","Y"):"R",("F","V"):"C",("P","S"):"R",("P","T"):"R",("P","W"):"C",("P","Y"):"R",("P","V"):"C",("S","T"):"C",("S","W"):"R",("S","Y"):"C",("S","V"):"R",("T","W"):"R",("T","Y"):"C",("T","V"):"R",("W","Y"):"R",("W","V"):"C",("Y","V"):"R",("H","R"):"C",("K","R"):"C",("D","R"):"R",("E","R"):"R",("A","R"):"R",("N","R"):"R",("C","R"):"R",("Q","R"):"R",("G","R"):"R",("I","R"):"R",("L","R"):"R",("M","R"):"R",("F","R"):"R",("P","R"):"R",("S","R"):"R",("T","R"):"R",("W","R"):"R",("Y","R"):"R",("V","R"):"R",("K","H"):"C",("D","H"):"R",("E","H"):"R",("A","H"):"R",("N","H"):"R",("C","H"):"R",("Q","H"):"R",("G","H"):"R",("I","H"):"R",("L","H"):"R",("M","H"):"R",("F","H"):"R",("P","H"):"R",("S","H"):"R",("T","H"):"R",("W","H"):"R",("Y","H"):"R",("V","H"):"R",("D","K"):"R",("E","K"):"R",("A","K"):"R",("N","K"):"R",("C","K"):"R",("Q","K"):"R",("G","K"):"R",("I","K"):"R",("L","K"):"R",("M","K"):"R",("F","K"):"R",("P","K"):"R",("S","K"):"R",("T","K"):"R",("W","K"):"R",("Y","K"):"R",("V","K"):"R",("E","D"):"C",("A","D"):"R",("N","D"):"R",("C","D"):"R",("Q","D"):"R",("G","D"):"R",("I","D"):"R",("L","D"):"R",("M","D"):"R",("F","D"):"R",("P","D"):"R",("S","D"):"R",("T","D"):"R",("W","D"):"R",("Y","D"):"R",("V","D"):"R",("A","E"):"R",("N","E"):"R",("C","E"):"R",("Q","E"):"R",("G","E"):"R",("I","E"):"R",("L","E"):"R",("M","E"):"R",("F","E"):"R",("P","E"):"R",("S","E"):"R",("T","E"):"R",("W","E"):"R",("Y","E"):"R",("V","E"):"R",("N","A"):"R",("C","A"):"R",("Q","A"):"R",("G","A"):"R",("I","A"):"C",("L","A"):"C",("M","A"):"C",("F","A"):"C",("P","A"):"C",("S","A"):"R",("T","A"):"R",("W","A"):"C",("Y","A"):"R",("V","A"):"C",("C","N"):"C",("Q","N"):"C",("G","N"):"C",("I","N"):"R",("L","N"):"R",("M","N"):"R",("F","N"):"R",("P","N"):"R",("S","N"):"C",("T","N"):"C",("W","N"):"R",("Y","N"):"C",("V","N"):"R",("Q","C"):"C",("G","C"):"C",("I","C"):"R",("L","C"):"R",("M","C"):"R",("F","C"):"R",("P","C"):"R",("S","C"):"C",("T","C"):"C",("W","C"):"R",("Y","C"):"C",("V","C"):"R",("G","Q"):"C",("I","Q"):"R",("L","Q"):"R",("M","Q"):"R",("F","Q"):"R",("P","Q"):"R",("S","Q"):"C",("T","Q"):"C",("W","Q"):"R",("Y","Q"):"C",("V","Q"):"R",("I","G"):"R",("L","G"):"R",("M","G"):"R",("F","G"):"R",("P","G"):"R",("S","G"):"C",("T","G"):"C",("W","G"):"R",("Y","G"):"C",("V","G"):"R",("L","I"):"C",("M","I"):"C",("F","I"):"C",("P","I"):"C",("S","I"):"R",("T","I"):"R",("W","I"):"C",("Y","I"):"R",("V","I"):"C",("M","L"):"C",("F","L"):"C",("P","L"):"C",("S","L"):"R",("T","L"):"R",("W","L"):"C",("Y","L"):"R",("V","L"):"C",("F","M"):"C",("P","M"):"C",("S","M"):"R",("T","M"):"R",("W","M"):"C",("Y","M"):"R",("V","M"):"C",("P","F"):"C",("S","F"):"R",("T","F"):"R",("W","F"):"C",("Y","F"):"R",("V","F"):"C",("S","P"):"R",("T","P"):"R",("W","P"):"C",("Y","P"):"R",("V","P"):"C",("T","S"):"C",("W","S"):"R",("Y","S"):"C",("V","S"):"R",("W","T"):"R",("Y","T"):"C",("V","T"):"R",("Y","W"):"R",("V","W"):"C",("V","Y"):"R",("R","*"):"R",("H","*"):"R",("K","*"):"R",("D","*"):"R",("E","*"):"R",("A","*"):"R",("N","*"):"R",("C","*"):"R",("Q","*"):"R",("G","*"):"R",("I","*"):"R",("L","*"):"R",("M","*"):"R",("F","*"):"R",("P","*"):"R",("S","*"):"R",("T","*"):"R",("W","*"):"R",("Y","*"):"R",("V","*"):"R",("*","R"):"R",("*","H"):"R",("*","K"):"R",("*","D"):"R",("*","E"):"R",("*","A"):"R",("*","N"):"R",("*","C"):"R",("*","Q"):"R",("*","G"):"R",("*","I"):"R",("*","L"):"R",("*","M"):"R",("*","F"):"R",("*","P"):"R",("*","S"):"R",("*","T"):"R",("*","W"):"R",("*","Y"):"R",("*","V"):"R"}
    criDict = {('*','A'):'R',('*','C'):'R',('*','D'):'R',('*','E'):'R',('*','F'):'R',('*','G'):'R',('*','H'):'R',('*','I'):'R',('*','K'):'R',('*','L'):'R',('*','M'):'R',('*','N'):'R',('*','P'):'R',('*','Q'):'R',('*','R'):'R',('*','S'):'R',('*','T'):'R',('*','V'):'R',('*','W'):'R',('*','Y'):'R',('A','*'):'R',('A','C'):'R',('A','D'):'R',('A','E'):'R',('A','F'):'C',('A','G'):'C',('A','H'):'R',('A','I'):'C',('A','K'):'R',('A','L'):'C',('A','M'):'C',('A','N'):'R',('A','P'):'C',('A','Q'):'R',('A','R'):'R',('A','S'):'C',('A','T'):'C',('A','V'):'C',('A','W'):'C',('A','Y'):'R',('C','*'):'R',('C','A'):'R',('C','D'):'R',('C','E'):'R',('C','F'):'R',('C','G'):'C',('C','H'):'R',('C','I'):'R',('C','K'):'R',('C','L'):'R',('C','M'):'R',('C','N'):'C',('C','P'):'R',('C','Q'):'C',('C','R'):'R',('C','S'):'C',('C','T'):'C',('C','V'):'R',('C','W'):'R',('C','Y'):'C',('D','*'):'R',('D','A'):'R',('D','C'):'R',('D','E'):'C',('D','F'):'R',('D','G'):'R',('D','H'):'R',('D','I'):'R',('D','K'):'R',('D','L'):'R',('D','M'):'R',('D','N'):'R',('D','P'):'R',('D','Q'):'R',('D','R'):'R',('D','S'):'R',('D','T'):'R',('D','V'):'R',('D','W'):'R',('D','Y'):'R',('E','*'):'R',('E','A'):'R',('E','C'):'R',('E','D'):'C',('E','F'):'R',('E','G'):'R',('E','H'):'R',('E','I'):'R',('E','K'):'R',('E','L'):'R',('E','M'):'R',('E','N'):'R',('E','P'):'R',('E','Q'):'R',('E','R'):'R',('E','S'):'R',('E','T'):'R',('E','V'):'R',('E','W'):'R',('E','Y'):'R',('F','*'):'R',('F','A'):'C',('F','C'):'R',('F','D'):'R',('F','E'):'R',('F','G'):'R',('F','H'):'R',('F','I'):'C',('F','K'):'R',('F','L'):'C',('F','M'):'C',('F','N'):'R',('F','P'):'C',('F','Q'):'R',('F','R'):'R',('F','S'):'R',('F','T'):'R',('F','V'):'C',('F','W'):'C',('F','Y'):'R',('G','*'):'R',('G','A'):'C',('G','C'):'C',('G','D'):'R',('G','E'):'R',('G','F'):'R',('G','H'):'R',('G','I'):'R',('G','K'):'R',('G','L'):'R',('G','M'):'R',('G','N'):'C',('G','P'):'C',('G','Q'):'C',('G','R'):'R',('G','S'):'C',('G','T'):'C',('G','V'):'R',('G','W'):'R',('G','Y'):'R',('H','*'):'R',('H','A'):'R',('H','C'):'R',('H','D'):'R',('H','E'):'R',('H','F'):'R',('H','G'):'R',('H','I'):'R',('H','K'):'C',('H','L'):'R',('H','M'):'R',('H','N'):'R',('H','P'):'R',('H','Q'):'R',('H','R'):'C',('H','S'):'R',('H','T'):'R',('H','V'):'R',('H','W'):'R',('H','Y'):'R',('I','*'):'R',('I','A'):'C',('I','C'):'R',('I','D'):'R',('I','E'):'R',('I','F'):'C',('I','G'):'R',('I','H'):'R',('I','K'):'R',('I','L'):'C',('I','M'):'C',('I','N'):'R',('I','P'):'C',('I','Q'):'R',('I','R'):'R',('I','S'):'R',('I','T'):'R',('I','V'):'C',('I','W'):'C',('I','Y'):'R',('K','*'):'R',('K','A'):'R',('K','C'):'R',('K','D'):'R',('K','E'):'R',('K','F'):'R',('K','G'):'R',('K','H'):'C',('K','I'):'R',('K','L'):'R',('K','M'):'R',('K','N'):'R',('K','P'):'R',('K','Q'):'R',('K','R'):'C',('K','S'):'R',('K','T'):'R',('K','V'):'R',('K','W'):'R',('K','Y'):'R',('L','*'):'R',('L','A'):'C',('L','C'):'R',('L','D'):'R',('L','E'):'R',('L','F'):'C',('L','G'):'R',('L','H'):'R',('L','I'):'C',('L','K'):'R',('L','M'):'C',('L','N'):'R',('L','P'):'C',('L','Q'):'R',('L','R'):'R',('L','S'):'R',('L','T'):'R',('L','V'):'C',('L','W'):'C',('L','Y'):'R',('M','*'):'R',('M','A'):'C',('M','C'):'R',('M','D'):'R',('M','E'):'R',('M','F'):'C',('M','G'):'R',('M','H'):'R',('M','I'):'C',('M','K'):'R',('M','L'):'C',('M','N'):'R',('M','P'):'C',('M','Q'):'R',('M','R'):'R',('M','S'):'R',('M','T'):'R',('M','V'):'C',('M','W'):'C',('M','Y'):'R',('N','*'):'R',('N','A'):'R',('N','C'):'C',('N','D'):'R',('N','E'):'R',('N','F'):'R',('N','G'):'C',('N','H'):'R',('N','I'):'R',('N','K'):'R',('N','L'):'R',('N','M'):'R',('N','P'):'R',('N','Q'):'C',('N','R'):'R',('N','S'):'C',('N','T'):'C',('N','V'):'R',('N','W'):'R',('N','Y'):'C',('P','*'):'R',('P','A'):'C',('P','C'):'R',('P','D'):'R',('P','E'):'R',('P','F'):'C',('P','G'):'C',('P','H'):'R',('P','I'):'C',('P','K'):'R',('P','L'):'C',('P','M'):'C',('P','N'):'R',('P','Q'):'R',('P','R'):'R',('P','S'):'C',('P','T'):'C',('P','V'):'C',('P','W'):'C',('P','Y'):'R',('Q','*'):'R',('Q','A'):'R',('Q','C'):'C',('Q','D'):'R',('Q','E'):'R',('Q','F'):'R',('Q','G'):'C',('Q','H'):'R',('Q','I'):'R',('Q','K'):'R',('Q','L'):'R',('Q','M'):'R',('Q','N'):'C',('Q','P'):'R',('Q','R'):'R',('Q','S'):'C',('Q','T'):'C',('Q','V'):'R',('Q','W'):'R',('Q','Y'):'C',('R','*'):'R',('R','A'):'R',('R','C'):'R',('R','D'):'R',('R','E'):'R',('R','F'):'R',('R','G'):'R',('R','H'):'C',('R','I'):'R',('R','K'):'C',('R','L'):'R',('R','M'):'R',('R','N'):'R',('R','P'):'R',('R','Q'):'R',('R','S'):'R',('R','T'):'R',('R','V'):'R',('R','W'):'R',('R','Y'):'R',('S','*'):'R',('S','A'):'C',('S','C'):'C',('S','D'):'R',('S','E'):'R',('S','F'):'R',('S','G'):'C',('S','H'):'R',('S','I'):'R',('S','K'):'R',('S','L'):'R',('S','M'):'R',('S','N'):'C',('S','P'):'C',('S','Q'):'C',('S','R'):'R',('S','T'):'C',('S','V'):'R',('S','W'):'R',('S','Y'):'C',('T','*'):'R',('T','A'):'C',('T','C'):'C',('T','D'):'R',('T','E'):'R',('T','F'):'R',('T','G'):'C',('T','H'):'R',('T','I'):'R',('T','K'):'R',('T','L'):'R',('T','M'):'R',('T','N'):'C',('T','P'):'C',('T','Q'):'C',('T','R'):'R',('T','S'):'C',('T','V'):'R',('T','W'):'R',('T','Y'):'C',('V','*'):'R',('V','A'):'C',('V','C'):'R',('V','D'):'R',('V','E'):'R',('V','F'):'C',('V','G'):'R',('V','H'):'R',('V','I'):'C',('V','K'):'R',('V','L'):'C',('V','M'):'C',('V','N'):'R',('V','P'):'C',('V','Q'):'R',('V','R'):'R',('V','S'):'R',('V','T'):'R',('V','W'):'C',('V','Y'):'R',('W','*'):'R',('W','A'):'C',('W','C'):'R',('W','D'):'R',('W','E'):'R',('W','F'):'C',('W','G'):'R',('W','H'):'R',('W','I'):'C',('W','K'):'R',('W','L'):'C',('W','M'):'C',('W','N'):'R',('W','P'):'C',('W','Q'):'R',('W','R'):'R',('W','S'):'R',('W','T'):'R',('W','V'):'C',('W','Y'):'C',('Y','*'):'R',('Y','A'):'R',('Y','C'):'C',('Y','D'):'R',('Y','E'):'R',('Y','F'):'R',('Y','G'):'R',('Y','H'):'R',('Y','I'):'R',('Y','K'):'R',('Y','L'):'R',('Y','M'):'R',('Y','N'):'C',('Y','P'):'R',('Y','Q'):'C',('Y','R'):'R',('Y','S'):'C',('Y','T'):'C',('Y','V'):'R',('Y','W'):'C'}
    totalSynSites = 0.0
    totalNonsynSites = 0.0
    totalC1Sites = 0.0
    totalR1Sites = 0.0
    totalC2Sites = 0.0
    totalR2Sites = 0.0
    totalC3Sites = 0.0
    totalR3Sites = 0.0
    totalC4Sites = 0.0
    totalR4Sites = 0.0
    totalC5Sites = 0.0
    totalR5Sites = 0.0
    totalC6Sites = 0.0
    totalR6Sites = 0.0
    totalC7Sites = 0.0
    totalR7Sites = 0.0
    totalMeanCSites = 0.0
    totalMeanRSites = 0.0
    codonNum = 0
    for codon in codonList:
        if 'N' in codon or '-' in codon:
            totalSynSites += 0.729166667
            totalNonsynSites += 2.270833333
            totalC1Sites += 1.395833333	
            totalR1Sites += 0.875	
            totalC2Sites += 1.270833333	
            totalR2Sites += 1	
            totalC3Sites += 0.708333333	
            totalR3Sites += 1.5625	
            totalC4Sites += 0.895833333	
            totalR4Sites += 1.375	
            totalC5Sites += 1.0625	
            totalR5Sites += 1.208333333	
            totalC6Sites += 0.854166667	
            totalR6Sites += 1.416666667	
            totalC7Sites += 0.8125	
            totalR7Sites += 1.458333333
            totalMeanCSites += 1.0 
            totalMeanRSites += 1.270833333
        else:
            currS = 0.0
            currN = 0.0
            currC1 = 0.0
            currC2 = 0.0
            currC3 = 0.0
            currC4 = 0.0
            currC5 = 0.0
            currC6 = 0.0
            currC7 = 0.0
            currR1 = 0.0
            currR2 = 0.0
            currR3 = 0.0
            currR4 = 0.0
            currR5 = 0.0
            currR6 = 0.0
            currR7 = 0.0
            currMeanC = 0.0
            currMeanR = 0.0
            site1 = codon[0]
            site2 = codon[1]
            site3 = codon[2]
            if site1 == 'A':
                mut1 = 'C' + site2 + site3
                mut2 = 'G' + site2 + site3
                mut3 = 'T' + site2 + site3
            elif site1 == 'C':
                mut1 = 'A' + site2 + site3
                mut2 = 'G' + site2 + site3
                mut3 = 'T' + site2 + site3
            elif site1 == 'G':
                mut1 = 'A' + site2 + site3
                mut2 = 'C' + site2 + site3
                mut3 = 'T' + site2 + site3
            elif site1 == 'T':
                mut1 = 'A' + site2 + site3
                mut2 = 'C' + site2 + site3
                mut3 = 'G' + site2 + site3
            if site2 == 'A':
                mut4 = site1 + 'C' + site3
                mut5 = site1 + 'G' + site3
                mut6 = site1 + 'T' + site3
            elif site2 == 'C':
                mut4 = site1 + 'A' + site3
                mut5 = site1 + 'G' + site3
                mut6 = site1 + 'T' + site3
            elif site2 == 'G':
                mut4 = site1 + 'A' + site3
                mut5 = site1 + 'C' + site3
                mut6 = site1 + 'T' + site3
            elif site2 == 'T':
                mut4 = site1 + 'A' + site3
                mut5 = site1 + 'C' + site3
                mut6 = site1 + 'G' + site3
            if site3 == 'A':
                mut7 = site1 + site2 + 'C'
                mut8 = site1 + site2 + 'G'
                mut9 = site1 + site2 + 'T'
            elif site3 == 'C':
                mut7 = site1 + site2 + 'A'
                mut8 = site1 + site2 + 'G'
                mut9 = site1 + site2 + 'T'
            elif site3 == 'G':
                mut7 = site1 + site2 + 'A'
                mut8 = site1 + site2 + 'C'
                mut9 = site1 + site2 + 'T'
            elif site3 == 'T':
                mut7 = site1 + site2 + 'A'
                mut8 = site1 + site2 + 'C'
                mut9 = site1 + site2 + 'G'
            if codonNum == 0:
                aaList = []
                if codon in startCodons:
                    currAA = 'M'
                else:
                    currAA = geneticCode[codon]
                if mut1 in startCodons:
                    aaList.append('M')
                else:
                    aaList.append(geneticCode[mut1])
                if mut2 in startCodons:
                    aaList.append('M')
                else:
                    aaList.append(geneticCode[mut2])
                if mut3 in startCodons:
                    aaList.append('M')
                else:
                    aaList.append(geneticCode[mut3])
                if mut4 in startCodons:
                    aaList.append('M')
                else:
                    aaList.append(geneticCode[mut4])
                if mut5 in startCodons:
                    aaList.append('M')
                else:
                    aaList.append(geneticCode[mut5])
                if mut6 in startCodons:
                    aaList.append('M')
                else:
                    aaList.append(geneticCode[mut6])
                if mut7 in startCodons:
                    aaList.append('M')
                else:
                    aaList.append(geneticCode[mut7])
                if mut8 in startCodons:
                    aaList.append('M')
                else:
                    aaList.append(geneticCode[mut8])
                if mut9 in startCodons:
                    aaList.append('M')
                else:
                    aaList.append(geneticCode[mut9])
            else:
                aaList = [geneticCode[mut1],geneticCode[mut2],geneticCode[mut3], geneticCode[mut4],geneticCode[mut5],geneticCode[mut6],geneticCode[mut7],geneticCode[mut8],geneticCode[mut9]]
                currAA = geneticCode[codon]
            for aa in aaList:
                if aa == currAA:
                    currS += 1.0
                else:
                    currN += 1.0
                    conRad1 = aaSchemeDict1[(currAA,aa)]
                    conRad2 = aaSchemeDict2[(currAA,aa)]
                    conRad3 = aaSchemeDict3[(currAA,aa)]
                    conRad4 = aaSchemeDict4[(currAA,aa)]
                    conRad5 = aaSchemeDict5[(currAA,aa)]
                    conRad6 = aaSchemeDict6[(currAA,aa)]
                    conRad7 = aaSchemeDict7[(currAA,aa)]
                    meanConRad = criDict[(currAA,aa)]
                    if conRad1 == 'R':
                        currR1 += 1
                    else:
                        currC1 += 1
                    if conRad2 == 'R':
                        currR2 += 1
                    else:
                        currC2 += 1
                    if conRad3 == 'R':
                        currR3 += 1
                    else:
                        currC3 += 1
                    if conRad4 == 'R':
                        currR4 += 1
                    else:
                        currC4 += 1
                    if conRad5 == 'R':
                        currR5 += 1
                    else:
                        currC5 += 1
                    if conRad6 == 'R':
                        currR6 += 1
                    else:
                        currC6 += 1
                    if conRad7 == 'R':
                        currR7 += 1
                    else:
                        currC7 += 1
                    if meanConRad == 'R':
                        currMeanR += 1
                    else:
                        currMeanC += 1
            currS /= 3.0
            currN /= 3.0
            currC1 /= 3.0
            currC2 /= 3.0
            currC3 /= 3.0
            currC4 /= 3.0
            currC5 /= 3.0
            currC6 /= 3.0
            currC7 /= 3.0
            currMeanC /= 3.0
            currR1 /= 3.0
            currR2 /= 3.0
            currR3 /= 3.0
            currR4 /= 3.0
            currR5 /= 3.0
            currR6 /= 3.0
            currR7 /= 3.0
            currMeanR /= 3.0
            totalSynSites += currS
            totalNonsynSites += currN
            totalC1Sites += currC1
            totalR1Sites += currR1
            totalC2Sites += currC2
            totalR2Sites += currR2
            totalC3Sites += currC3
            totalR3Sites += currR3
            totalC4Sites += currC4
            totalR4Sites += currR4
            totalC5Sites += currC5
            totalR5Sites += currR5
            totalC6Sites += currC6
            totalR6Sites += currR6
            totalC7Sites += currC7
            totalR7Sites += currR7
            totalMeanCSites += currMeanC
            totalMeanRSites += currMeanR
            codonNum += 1
    return [totalSynSites,totalNonsynSites,totalC1Sites,totalR1Sites,totalC2Sites,totalR2Sites,totalC3Sites,totalR3Sites,totalC4Sites,totalR4Sites,totalC5Sites,totalR5Sites,totalC6Sites,totalR6Sites,totalC7Sites,totalR7Sites,totalMeanCSites,totalMeanRSites]
    
polymorphismsSubstitutions(sys.argv[1],int(sys.argv[2]))
#uniquePolymorphismsSubstitutions(sys.argv[1],int(sys.argv[2]))