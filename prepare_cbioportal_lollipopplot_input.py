#!/usr/bin/python

import sys
import re

def align_dic(Human_canine_sequenceAlignment):
        with open(Human_canine_sequenceAlignment,'r') as f:
                file = f.read()
        lst = file.split('\n')[1:-1]
        dic = {}
        for i in range(len(lst)):
                gene,humanTrans,dogTrans,humanPos,dogPos,humanAA,dogAA = lst[i].split('\t')
                if gene not in dic.keys():
                        dic[gene] = {}
                if dogTrans not in dic[gene].keys():
                        dic[gene][dogTrans] = {}
                dic[gene][dogTrans][dogPos] = [humanPos,humanAA,dogAA]
        return dic

if __name__ == '__main__':
        file0 = sys.argv[1]
        Human_canine_sequenceAlignment = sys.argv[2]

        alignDic = align_dic(Human_canine_sequenceAlignment)

        with open(file0,'r') as f:
                file = f.read()
        lst = file.split('\n')[1:-1]

        out = open(file0.split('.')[0] + '_cBioPortalLollipopInput.txt','w')
        header = 'Hugo_Symbol\tProtein_Change\n'
        out.write(header)

        out2 = open(file0.split('.')[0] + '_cBioPortalLollipopPairs.txt','w')

        for i in range(len(lst)):
                if 'Somatic' in lst[i]:
                        info = lst[i].split('\t')
                        gene = info[5]
                        AAchange = info[11]
                        dogTrans = info[15]
                        if gene in alignDic.keys():
                                if len(re.findall(r'\d+', AAchange)) == 1:
                                        dogPos = re.findall(r'\d+', AAchange)[0]
                                        refAA = re.findall(r'\D+', AAchange)[0]
                                        altAA = re.findall(r'\D+', AAchange)[1]
                                        if dogPos in alignDic[gene][dogTrans].keys():
                                                humanPos,humanAA,dogAA = alignDic[gene][dogTrans][dogPos]
                                                adjAAchange = refAA + humanPos + altAA
                                                string = gene + '\t' + adjAAchange + '\n'
                                                out.write(string)
                                                out2.write(gene + '\t' + dogTrans + '\t' + AAchange + '\t' + adjAAchange + '\n')
                                        else:
                                                while (int(dogPos) > 1) and (dogPos not in alignDic[gene][dogTrans].keys()):
                                                        dogPos = str(int(dogPos) - 1)
                                                if dogPos in alignDic[gene][dogTrans].keys():
                                                        humanPos,humanAA,dogAA = alignDic[gene][dogTrans][dogPos]
                                                        adjAAchange = refAA + humanPos + altAA
                                                        string = gene + '\t' + adjAAchange + '\n'
                                                        out.write(string)
                                                        out2.write(gene + '\t' + dogTrans + '\t' + AAchange + '\t' + adjAAchange + '\n')
                                                else:
                                                        print ('## ' + gene + '\t' + dogTrans + '\t' + AAchange)
                                else: # indel
                                        dogPosLst = re.findall(r'\d+', AAchange)
                                        AALst = re.findall(r'\D+', AAchange)
                                        adjDogPosLst = []
                                        for dogPos in dogPosLst:
                                                if dogPos in alignDic[gene][dogTrans].keys():
                                                        humanPos,humanAA,dogAA = alignDic[gene][dogTrans][dogPos]
                                                        adjDogPosLst.append(humanPos)
                                                else:
                                                        while (int(dogPos) > 1) and (dogPos not in alignDic[gene][dogTrans].keys()):
                                                                dogPos = str(int(dogPos) - 1)
                                                        if dogPos in alignDic[gene][dogTrans].keys():
                                                                humanPos,humanAA,dogAA = alignDic[gene][dogTrans][dogPos]
                                                                adjDogPosLst.append(humanPos)
                                                        else:
                                                                print ('### ' + gene + '\t' + dogTrans + '\t' + AAchange)
                                        tmp = []
                                        if len(adjDogPosLst) == len(AALst) - 1:
                                                for i in range(len(adjDogPosLst)):
                                                        tmp.append(AALst[i])
                                                        tmp.append(adjDogPosLst[i])
                                                tmp.append(AALst[-1])
                                        elif len(adjDogPosLst) == len(AALst):
                                                for i in range(len(adjDogPosLst)):
                                                        tmp.append(adjDogPosLst[i])
                                                        tmp.append(AALst[i])
                                        if tmp != []:
                                                adjAAchange = ''.join(tmp)
                                                string = gene + '\t' + adjAAchange + '\n'
                                                out.write(string)
                                                out2.write(gene + '\t' + dogTrans + '\t' + AAchange + '\t' + adjAAchange + '\n')
                                        else:
                                                print adjDogPosLst
                                                print AALst
                                                print ('### ' + gene + '\t' + dogTrans + '\t' + AAchange)
        out.close()
        out2.close()



