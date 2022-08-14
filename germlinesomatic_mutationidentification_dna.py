#!/usr/bin/python

import sys
import re

def normal_mut_dic(normalMut):
        with open(normalMut,'r') as f:
                file = f.read()
        lst = file.split('\n')[1:-1]
        dic = {} # gene -> [mut(nuc),]
        for i in range(len(lst)):
                info = lst[i].split('\t')
                gene = info[5]
                transcriptID = info[15]
                AAchange = info[11]
                if gene not in dic.keys():
                        dic[gene] = {}
                if transcriptID not in dic[gene].keys():
                        dic[gene][transcriptID] = []
                if AAchange not in dic[gene][transcriptID]:
                        dic[gene][transcriptID].append(AAchange)
        return dic

def annovar_exonic_mut_dic(file0,geneList):
        with open(file0,'r') as f:
                file = f.read()
        lst = file.split('\n')[:-1]
        dic = {}
        for i in range(len(lst)):
                info = lst[i].split('\t')
                gene = info[2]
                if gene in geneList:
                        tagLst = info[3].split(',')[:-1]
                        for tag in tagLst:
                                if ('ENSCAFT' in tag) and ('p.' in tag):
                                        transcriptID = tag.split(':')[1]
                                        protChange = tag.split(':')[-1].split('p.')[-1]
                                        if gene not in dic.keys():
                                                dic[gene] = {}
                                        if transcriptID not in dic[gene].keys():
                                                dic[gene][transcriptID] = []
                                        if protChange not in dic[gene][transcriptID]:
                                                dic[gene][transcriptID].append(protChange)
        return dic

def pairIdx_dic(humanCaninePosPair):
        with open(humanCaninePosPair,'r') as f:
                file = f.read()
        lst = file.split('\n')[1:-1]
        dic = {}
        for i in range(len(lst)):
                gene,refPos,queryPos,refAA,queryAA = lst[i].split('\t')
                if gene not in dic.keys():
                        dic[gene] = {}
                dic[gene][refPos] = [queryPos,refAA,queryAA]
        return dic

def human2canine_prot_convert(prot,gene,pairIdxDic):
        pos2 = 'NA'
        if len(re.findall(r'\d+', prot)) == 1 and len(re.findall(r'\D+', prot)) == 2:
                pos = re.findall(r'\d+', prot)[0]
                nuc1,nuc2 = re.findall(r'\D+', prot)
                if gene in pairIdxDic.keys():
                        if pos in pairIdxDic[gene].keys():
                                pos2 = pairIdxDic[gene][pos][0]
        return pos2

def cosmic_somatic_dic(cosmic,geneList,pairIdxDic):
        with open(cosmic,'r') as f:
                file = f.read()
        lst = file.split('\n')[1:-1]
        dic = {}
        for i in range(len(lst)):
                gene = lst[i].split('\t')[0].split('_')[0]
                if gene in geneList:
                        prot = lst[i].split('\t')[20].split('p.')[1]
                        if gene not in dic.keys():
                                dic[gene] = []
                        if prot != "?":
                                dogPos = human2canine_prot_convert(prot,gene,pairIdxDic)
                                if dogPos not in dic[gene]:
                                        dic[gene].append(dogPos)
        return dic

def cBioPortal_somatic_dic(cBioPortal,geneList,pairIdxDic):
        with open(cBioPortal,'r') as f:
                file = f.read()
        headerLst = file.split('\n')[0].split('\t')[4:]
        lst = file.split('\n')[1:-1]
        # gene idx
        idxDic = {}
        maxIdx = -1
        for i in range(len(headerLst)):
                gene = headerLst[i]
                if ':' not in gene:
                        maxIdx = i
                        idxDic[i] = gene
        # analyze samples
        dic = {}
        for i in range(len(lst)):
                mutLst = lst[i].split('\t')[4:4+maxIdx+1]
                for j in range(len(mutLst)):
                        mut = mutLst[j]
                        gene = idxDic[j]
                        if mut != 'no alteration':
                                tmpLst = mut.split(',')
                                for tmp in tmpLst:
                                        if len(re.findall(r'\d+', tmp)) == 1 and len(re.findall(r'\D+', tmp)) == 2:
                                                dogPos = human2canine_prot_convert(tmp,gene,pairIdxDic)
                                                if gene not in dic.keys():
                                                        dic[gene] = []
                                                if dogPos not in dic[gene]:
                                                        dic[gene].append(dogPos)
        return dic

def mutation_frequency(matrix0,geneList):
        with open(matrix0,'r') as f:
                file = f.read()
        lst = file.split('\n')[1:-1]
        dic = {}
        for i in range(len(lst)):
                info = lst[i].split('\t')
                if info[16] != 'synonymous SNV':
                        sample = info[1]
                        gene = info[5]
                        if gene in geneList:
                                transcript = info[15]
                                AAchange = info[11]
                                VAF = float(info[14])
                                if gene not in dic.keys():
                                        dic[gene] = {}
                                if transcript not in dic[gene].keys():
                                        dic[gene][transcript] = {}
                                if AAchange not in dic[gene][transcript].keys():
                                        dic[gene][transcript][AAchange] = {}
                                if sample not in dic[gene][transcript][AAchange]:
                                        dic[gene][transcript][AAchange][sample] = VAF
        return dic

def check_annovar_snp(gene0,transcript0,AAchange0,dic,tag):
        result = '' # return tag if the AAchange is found in dbSNP, else ''
        if gene0 in dic.keys():
                if transcript0 in dic[gene0]:
                        if AAchange0 in dic[gene0][transcript0]:
                                result = tag
                else:
                        transcriptLst = dic[gene0].keys()
                        for transcript in transcriptLst:
                                if AAchange0 in dic[gene0][transcript]:
                                        result = tag
        return result
                
def distribution_density(dic,density):
        sampleLst = dic.keys()
        VAFlst = []
        for sample in sampleLst:
                VAFlst.append(dic[sample])
        satisfiedSampleNumber = 0
        for VAF in VAFlst:
                if (VAF >= 0.4 and VAF <= 0.6) or (VAF >= 0.9):
                        satisfiedSampleNumber += 1
        ratio = float(satisfiedSampleNumber) / len(VAFlst)
        flag = False
        if ratio >= density:
                flag = True # Germline
        return flag

if __name__ == '__main__':
        #matrix0 = '/scratch/yf94402/FidoCure/result/NewAnalysis_9-21-2021/FidoCure_NewAnalysis_MissingMetaAdded_2-11-2022_targetGenes_tumor.txt'
        matrix0 = '/scratch/yf94402/FidoCure/result/NewAnalysis_9-21-2021/FidoCure_NewAnalysis_MissingMetaAdded_2-11-2022_transcriptUpdate8-7-22_targetGenes_tumor.txt'
        # GeneList
        geneList = '/scratch/yf94402/FidoCure/data/TargetGeneList_69Genes.txt'
        # SNP
        #normalMut = '/scratch/yf94402/FidoCure/result/NewAnalysis_9-21-2021/FidoCure_NewAnalysis_MissingMetaAdded_2-11-2022_targetGenes_normal.txt'
        normalMut = '/scratch/yf94402/FidoCure/result/NewAnalysis_9-21-2021/FidoCure_NewAnalysis_MissingMetaAdded_2-11-2022_transcriptUpdate8-7-22_targetGenes_normal.txt'
        UCSC = '/scratch/yf94402/FidoCure/source/UCSC_PON/UCSC_Germline_PON_CanFam3.vcf-avinput.exonic_variant_function_WithGeneName'
        PanCancer_PON = '/scratch/yf94402/FidoCure/source/Mutect2_PON/pon.vcf-PASS-avinput.exonic_variant_function_WithGeneName'
        dbSNP = '/scratch/yf94402/FidoCure/source/lab_dbSNP/DbSNP_canFam3_version151-DogSD_Feb2020_V3.vcf-avinput.exonic_variant_function_WithGeneName'
        Literature1 = '/scratch/yf94402/FidoCure/source/LiteratureSNP/722g.990.SNP.INDEL.chrAll.vcf-PASS-avinput.exonic_variant_function_WithGeneName'
        Literature2 = '/scratch/yf94402/FidoCure/source/LiteratureSNP/dogs.590publicSamples.vcf-PASS-avinput.exonic_variant_function_WithGeneName'
        # Human somatic mutations
        cosmic = '/scratch/yf94402/FidoCure/source/Cosmic/CosmicMutantExport_FidoCure_gene.tsv'
        cBioPortal = '/scratch/yf94402/FidoCure/source/cBioPortal/alterations_across_samples.tsv'
        humanCaninePosPair = '/scratch/yf94402/FidoCure/result/ClustalO_allGenes/Human_canine_sequenceAlignment.txt'

        # analysis
        with open(geneList,'r') as f:
                file = f.read()
        geneList = file.split('\n')[:-1]

        normalMutDic = normal_mut_dic(normalMut)

        UCSCDic = annovar_exonic_mut_dic(UCSC,geneList)
        PONDic = annovar_exonic_mut_dic(PanCancer_PON,geneList)
        dbSNPDic = annovar_exonic_mut_dic(PanCancer_PON,geneList)
        Literature1Dic = annovar_exonic_mut_dic(Literature1,geneList)
        Literature2Dic = annovar_exonic_mut_dic(Literature2,geneList)

        pairIdxDic = pairIdx_dic(humanCaninePosPair)
        cosmicDic = cosmic_somatic_dic(cosmic,geneList,pairIdxDic)        
        cBioPortalDic = cBioPortal_somatic_dic(cBioPortal,geneList,pairIdxDic)

        mutFreqDic = mutation_frequency(matrix0,geneList)

        print cosmicDic
        print cBioPortalDic

        # germline/somatic labeling
        classificationDic = {}
        geneLst = mutFreqDic.keys()
        for gene in geneLst:
                transcriptLst = mutFreqDic[gene].keys()
                for transcript in transcriptLst:
                        AAchangeLst = mutFreqDic[gene][transcript]
                        for AAchange in AAchangeLst:
                                # germline labeling
                                germlineLst = []
                                tag = ''
                                if gene in normalMutDic.keys():
                                        if transcript in normalMutDic[gene].keys():
                                                if AAchangeLst in normalMutDic[gene][transcript]:
                                                        tag = 'Germline(NormalSample)'
                                if tag != '':
                                        germlineLst.append(tag)
                                tag = check_annovar_snp(gene,transcript,AAchange,UCSCDic,'Germline(UCSC)')
                                if tag != '':
                                        germlineLst.append(tag)
                                tag = check_annovar_snp(gene,transcript,AAchange,PONDic,'Germline(PON)')
                                if tag != '':
                                        germlineLst.append(tag)
                                tag = check_annovar_snp(gene,transcript,AAchange,dbSNPDic,'Germline(dbSNP)')
                                if tag != '':
                                        germlineLst.append(tag)
                                tag = check_annovar_snp(gene,transcript,AAchange,Literature1Dic,'Germline(Literature1)')
                                if tag != '':
                                        germlineLst.append(tag)
                                tag = check_annovar_snp(gene,transcript,AAchange,Literature2Dic,'Germline(Literature2)')
                                if tag != '':
                                        germlineLst.append(tag)
                                if germlineLst == []:
                                        germlineLst = []
                                        sampleLst = mutFreqDic[gene][transcript][AAchange].keys()
                                        if len(sampleLst) >= 5: # VAF distribution
                                                tmpDic = mutFreqDic[gene][transcript][AAchange] # tmpDic[sample] = VAF
                                                flag = distribution_density(tmpDic,0.8) # VAF density cutoff
                                                if flag:
                                                        germlineLst.append('Germline(VAFdist;n>=5)')
                                                else:
                                                        germlineLst.append('Somatic(VAFdist;n>=5)')
                                        else:
                                                pos = re.findall(r'\d+', AAchange)[0]
                                                if gene in cosmicDic.keys():
                                                        if pos in cosmicDic[gene]:
                                                                germlineLst.append('Somatic(Cosmic;n<5)')
                                                if gene in cBioPortalDic.keys():
                                                        if pos in cBioPortalDic[gene]:
                                                                germlineLst.append('Somatic(cBioPortal;n<5)')
                                        if germlineLst == []:
                                                germlineLst = ['Unclassifiable(n<5)']
                                # append classificationDic
                                if gene not in classificationDic.keys():
                                        classificationDic[gene] = {}
                                if transcript not in classificationDic[gene].keys():
                                        classificationDic[gene][transcript] = {}
                                classificationDic[gene][transcript][AAchange] = germlineLst

        # generate matrix
        with open(matrix0,'r') as f:
                file = f.read()
        header = file.split('\n')[0]
        newHeader = header + '\t' + 'Classification\n'
        lst = file.split('\n')[1:-1]

        out = open(matrix0.split('.')[0] + '_germlineSomaticClassification.txt','w')
        out.write(newHeader)

        for i in range(len(lst)):
                info = lst[i].split('\t')
                if info[16] != 'synonymous SNV':
                        sample = info[1]
                        gene = info[5]
                        if gene in geneList:
                                transcript = info[15]
                                AAchange = info[11]
                                classificationLst = classificationDic[gene][transcript][AAchange]
                                string = lst[i] + '\t' + ';'.join(classificationLst) + '\n'
                                out.write(string)
        out.close()


        

