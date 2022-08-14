#!/usr/bin/python

import sys

def sequence_dic(file0):
        dic = {}
        with open(file0,'r') as f:
                file = f.read()
        lst = file.split('>')[1:]
        for i in range(len(lst)):
                transcript = lst[i].split('\n')[0].split('transcript:')[1].split(' ')[0].split('.')[0]
                seq = ''.join(lst[i].split('\n')[1:])
                dic[transcript] = seq
        return dic

def transcript_dic(transcriptPair):
        dic = {}
        with open(transcriptPair,'r') as f:
                file = f.read()
        lst = file.split('\n')[:-1]
        for i in range(len(lst)):
                gene,refT,queryT = lst[i].split('\t')
                dic[gene] = [refT,queryT]
        return dic

if __name__ == '__main__':
	transcriptPair = sys.argv[1] # gene \t refTransID \t queryTransID \n
        protRef = sys.argv[2] # human
        protQuery = sys.argv[3] # canine
        resultPath = sys.argv[4]
        # preparation
        protRefDic = sequence_dic(protRef)
        protQueryDic = sequence_dic(protQuery)
        transPairDic = transcript_dic(transcriptPair)
        # output
        geneLst = transPairDic.keys()
        for gene in geneLst:
                refT,queryT = transPairDic[gene]
                qTlst = queryT.split(';')
                for qT in qTlst:
                        if (refT in protRefDic.keys()) and (qT in protQueryDic.keys()): # check if gene sequences exist
                                out = open(resultPath + '/' + gene + '_' + refT + '_' + qT + '.fa','w')
                                refSeq = protRefDic[refT]
                                querySeq = protQueryDic[qT]
                                stringRef = '>' + gene + '_' + refT + '_ref' + '\n' + refSeq + '\n'
                                out.write(stringRef)
                                stringQuery = '>' + gene + '_' + qT + '_query' + '\n' + querySeq + '\n'
                                out.write(stringQuery)
                                out.close()


