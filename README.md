# Confirmatory Tools 
# WIsH 
Null directory is created from a diverse range of Alteromonas, Cellulophage, Cyanophage, Lactobacillus, Mycobacterium, Oenococcus, Pelagibacter, Prochlorococcus, Rhizobium, Synechococcus, and Thermus phage genomes that are known not to infect the bacterial model <br /> 
WIsH.sh 

# PHIRBO
Download taxdb.btd, taxdb.bti, taxdb.tar.gz from https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz <br />
Prepare blast results via blast_forphirbo.sh <br />
PHIRBO.sh 

# PHP
PHP.sh

# PHIST
PHIST.sh

# VirHostMatcher 
VHM.sh

# Exploratory Tools 
# vpf-class 

# vHULK

# RaFAH

# HostG

# VirHostMatcher-Net
apptainer build vhmn.sif vhmn.def <br /> 
apptainer shell vhmn.sif <br /> 
cd /opt/miniconda/bin  <br /> 
. activate <br /> 
cd <br /> 
python /opt/VirHostMatcher-Net/VirHostMatcher-Net.py -q phage_genomes -o phage_vhmn_out -i tmp -n 3 -t 8 <br /> 

# CHERRY
apptainer build cherry.sif cherry.def <br /> 
apptainer shell cherry.sif <br /> 
cd /opt/miniconda/bin  <br /> 
. activate <br /> 
cd <br />  
python run_Speed_up.py --contigs phage_genome.fasta --model pretrain --topk 1000

