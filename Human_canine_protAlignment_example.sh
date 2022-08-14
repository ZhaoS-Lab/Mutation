#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=clustalO
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem=10G

module load Clustal-Omega/1.2.4-GCC-8.3.0

data='/scratch/yf94402/FidoCure/source/NewAnalysis_7-24-2022' # the peptide files from the 2 species
result='/scratch/yf94402/FidoCure/result/ClustalO_human107_canine99' # output directory where you will see all the sequence alignments and result summary (from Step3)
script='/scratch/yf94402/FidoCure/script/NewAnalysis_9-21-2021/CrossSpecies_HomologousSiteIdentification'

mkdir $result
cd $result

### clustalO preparation
python $script/clustalo_preparaion.py \
$data/FidoCure_geneLst_cBioPortal_humanTransID_withCanineTrans.txt \
$data/Homo_sapiens.release107.GRCh38.pep.all.fa \
$data/Canis_familiaris.release99.CanFam3.1.pep.all.fa \
$result


# Reference/Query transcript pairs (gene + '\t' + refTranscript + '\t' + queryTranscript + '\n')
# reference protein file
# query protein file
# clustalO output path



### clustalO 
cd $result
for i in *.fa
do
clustalo --force -i $i -o ${i/.fa/.clustalO}
done



### alignment analyzer
python $script/alignment_summary.py \
$result \
Human_canine_sequenceAlignment.txt

# clustalO path
# output file name

