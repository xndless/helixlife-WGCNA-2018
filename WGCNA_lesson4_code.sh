# Installation
# IterativeWGCNA is reposited in the Python Package Index (PyPI) 
# and can be installed via pip or easy_install.
pip install iterativeWGCNA

# You can also use pip to install from the git master:
pip install git+git://github.com/cstoeckert/iterativeWGCNA.git

# To install the iterativeWGCNA package from the git master, 
# clone and then run the python setup.py script as folows:
git clone https://github.com/cstoeckert/iterativeWGCNA.git
cd iterativeWGCNA
python setup.py install
python setup.py install --user

###################################################################################################

mkdir helix-WGCNA-lesson4
cd helix-WGCNA-lesson4
wget https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-Data.zip
unzip FemaleLiver-Data.zip
ls
## ClinicalTraits.csv  FemaleLiver-Data.zip  GeneAnnotation.csv  LiverFemale3600.csv

cut -d , -f 1,9- LiverFemale3600.csv | \
perl -pe 's/,/\t/g' >inputFile_expr.txt

perl -i -pe 's/NA/0/g' inputFile_expr.txt

perl -pe 's/, / /' ClinicalTraits.csv | \
cut -d , -f 2,11-15,17-30 | \
perl -pe 's/,/\t/g;s/"//g' >inputFile_trait.txt

perl -i -pe 's/NA/0/g' inputFile_trait.txt

###################################################################################################

iterativeWGCNA --inputFile $PWD/inputFile_expr.txt \
--workingDir LiverFemale3600 \
--wgcnaParameters power=22,minModuleSize=30,saveTOMs=TRUE,minKMEtoStay=0.8,minCoreKME=0.8,\
networkType=signed,numericLabels=TRUE,maxBlockSize=30000 \
--enableWGCNAThreads \
--finalMergeCutHeight 0.05

tree
## ├── final-eigengenes.txt
## ├── final-kme_histogram.pdf
## ├── final-membership.txt
## ├── iterativeWGCNA.log
## ├── iterativeWGCNA-R.log
## ├── iterative-wgcna-run-summary.txt
## ├── merged-0.05-eigengenes.txt
## ├── merged-0.05-kme_histogram.pdf
## ├── merged-0.05-membership.txt
## ├── pass1
## │   ├── i1
## │   │   ├── eigengenes.txt
## │   │   ├── kme_histogram.pdf
## │   │   ├── membership.txt
## │   │   ├── P1_I1-TOM-block.1.RData
## │   │   ├── summary.txt
## │   │   ├── wgcna-blocks.RData
## │   │   ├── wgcna-kme_histogram.pdf
## │   │   └── wgcna-membership.txt
## │   │   ...
## │   ├── initial-pass-expression-set.txt
## │   ├── kme_histogram.pdf
## │   └── membership.txt

cd LiverFemale3600

ls pass*/i*/*-TOM-block.1.RData | while read file
do
folder=$(echo $file | perl -pe 's%(pass.+/i.+/).+-TOM-block.1.RData%$1%')
run=$(echo $file | perl -pe 's%pass.+/i.+/(.+)-TOM-block.1.RData%$1%')
cd $folder
Rscript ../../iterativeWGCNA_downstream_analysis_01.r $run
cd ../..
done

###################################################################################################

Rscript iterativeWGCNA_downstream_analysis_02.r \
merged-0.05-eigengenes.txt \
merged-0.05-membership.txt \
../inputFile_expr.txt \
01.expression_ME

###################################################################################################

Rscript iterativeWGCNA_downstream_analysis_03.r \
merged-0.05-eigengenes.txt \
../inputFile_expr.txt \
../inputFile_trait.txt \
02.relationship_analysis

###################################################################################################

mkdir 03.module_result
cut -f 2 merged-0.05-membership.txt | sed '1d' | sort -u | while read module
do
head -1 merged-0.05-membership.txt >03.module_result/${module}-module-gene.txt
grep -w $module merged-0.05-membership.txt >>03.module_result/${module}-module-gene.txt
done

cd 03.module_result
ls P*-module-gene.txt | while read MG
do
P=$(echo $MG | perl -pe 's/P(.+)_I.+_M.+-module-gene.txt/$1/')
I=$(echo $MG | perl -pe 's/P.+_I(.+)_M.+-module-gene.txt/$1/')
M=$(echo $MG | perl -pe 's/P.+_I.+_M(.+)-module-gene.txt/$1/')
cut -f 1 $MG >module_gene.list
python pick_targeted_nodes_edges.py \
../pass$P/i${I}/module_result/CytoscapeInput-edges-P${P}_I${I}.txt \
module_gene.list \
CytoscapeInput-edges-P${P}_I${I}_M${M}.txt

echo -e "Gene\tIntraModuleTO" >P${P}_I${I}_M${M}-MTO-sorting.txt
cat module_gene.list | while read gene
do
grep -w $gene CytoscapeInput-edges-P${P}_I${I}_M${M}.txt | \
awk 'BEGIN {OFS="\t";FS="\t"} {sum+=$3} END {print "'$gene'",sum}'
done | sort -k 2nr >>P${P}_I${I}_M${M}-MTO-sorting.txt
done
rm module_gene.list

ls *-MTO-sorting.txt | while read line
do
module=$(echo $line | perl -pe 's%(.+)-MTO-sorting.txt%$1%')
nlines=$(cat $line | wc -l)
nlines_P5=$(((nlines-1)/5))
head -n $[$nlines_P5+1] $module-MTO-sorting.txt >$module-5%hub.txt
done
mkdir hub_gene && mv *hub* hub_gene/
