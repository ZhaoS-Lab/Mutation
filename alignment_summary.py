#!/usr/bin/python

import sys
import glob

def sequence_alignment_analyzer(refSeq,querySeq):
        dic = {}
        idxR = 0
        idxQ = 0
        for i in range(len(refSeq)):
                if querySeq[i] != '-':
                        idxQ += 1
                if refSeq[i] != '-':
                        idxR += 1
                        dic[idxR] = [idxQ,refSeq[i],querySeq[i]]
        return dic
      

if __name__ == '__main__':
        clustalO_path = sys.argv[1]
        outputName = sys.argv[2]
        # analysis
        fileLst = glob.glob(clustalO_path + '/*.clustalO')
        out = open(outputName,'w')
        header = 'Gene\tRefTransID\tQueryTransID\tRefIdx\tQueryIdx\tRefAA\tQueryAA\n'
        out.write(header)
        for file0 in fileLst:
                with open(file0,'r') as f:
                        file = f.read()
                lst = file.split('>')[1:]
                refTrans = file0.split('/')[-1].split('_')[1]
                queryTrans = file0.split('/')[-1].split('_')[2].split('.')[0]
                refSeq = ''
                querySeq = ''
                gene = ''
                for i in range(len(lst)):
                        if 'ref' in lst[i]:
                                gene = lst[i].split('\n')[0].split('_')[0]
                                refSeq = ''.join(lst[i].split('\n')[1:])
                        elif 'query' in lst[i]:
                                querySeq = ''.join(lst[i].split('\n')[1:])
                        if (refSeq != '') and (querySeq != '') and (gene != ''):
                                alignDic = sequence_alignment_analyzer(refSeq,querySeq)
                                posLst = alignDic.keys()
                                posLst.sort()
                                for pos in posLst:
                                        idxQ,refAA,queryAA = alignDic[pos]
                                        string = gene + '\t' + refTrans + '\t' + queryTrans + '\t' + str(pos) + '\t' + str(idxQ) + '\t' + refAA + '\t' + queryAA + '\n'
                                        out.write(string)
        out.close()
