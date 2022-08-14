#!/usr/bin/python

import sys

def humanTrans_dic(humanTransLst):
        dic = {}
        with open(humanTransLst,'r') as f:
                file = f.read()
        lst = file.split('\n')[:-1]
        for i in range(len(lst)):
                gene,trans = lst[i].split('\t')
                dic[gene] = trans
        return dic

def canineTrans_dic(table0,humanTransDic):
        geneLst = humanTransDic.keys()
        with open(table0,'r') as f:
                file = f.read()
        dic = {}
        for gene in geneLst:
                dic[gene] = []
        lst = file.split('\n')[1:-1]
        for i in range(len(lst)):
                info = lst[i].split('\t')
                gene = info[5]
                trans = info[15]
                if gene in geneLst:
                        if trans not in dic[gene]:
                                dic[gene].append(trans)
        return dic

if __name__ == '__main__':
        table0 = sys.argv[1]
        humanTransLst = sys.argv[2]

        humanTransDic = humanTrans_dic(humanTransLst)
        dogTransDic = canineTrans_dic(table0,humanTransDic)

        out = open(humanTransLst.split('.')[0] + '_withCanineTrans.txt','w')

        geneLst = humanTransDic.keys()
        for gene in geneLst:
                human = humanTransDic[gene]
                dog = ';'.join(dogTransDic[gene])
                string = gene + '\t' + human + '\t' + dog + '\n'
                out.write(string)
        out.close()

        

