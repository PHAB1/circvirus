#!/usr/bin/env python

import os
import sys
import subprocess

def conv_B_M5(file,arch):
    blastOutFMT = 'sseqid sstart send slen qstart qend qlen evalue score length nident mismatch gaps sseq qseq qseqid'.split()
    blasrFMT = 'qName qLength qStart qEnd qStrand tName tLength tStart tEnd tStrand score numMatch numMismatch numIns numDel mapQV qAlignedSeq matchPattern tAlignedSeq'

    for line in file:
        fields = line.strip().split()
        
        nident = 0
        for snuc,qnuc in zip(fields[12],fields[13]):
            if snuc == qnuc:
                nident += 1
        
        fields = fields[0:10]+[str(nident)]+fields[10:]
        record = dict(zip(blastOutFMT, fields))
        output = [record['sseqid'], record['slen']] ##
        if int(record['sstart']) < int(record['send']):
            output += [ str(int(record['sstart'])-1), record['send'], "+ "]
        else:
            output += [ str(int(record['send'])-1), record['sstart'], "- "]
        output += [record['qseqid'], record['qlen'], str(int(record['qstart'])-1), record['qend'], '+']
        output += ['-3000']  ## a fake score
        output += [record['nident'], record['mismatch']]
        output += [str(record['qseq'].count('-'))]
        output += [str(record['sseq'].count('-'))]
        output += ['254'] ## fake mapQV
        output += [record['sseq']]
        aln = ''
        for i,j in zip(record['qseq'],record['sseq']):
            aln += '|' if i==j else '*'
        output += [aln]
        output += [record['qseq']]
        
        arch.write(' '.join(output)+'\n')

    arch.close()
