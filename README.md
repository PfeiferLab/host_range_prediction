# Confirmatory Tools 
### WIsH 
Null directory is created from a diverse range of Alteromonas, Cellulophage, Cyanophage, Lactobacillus, Mycobacterium, Oenococcus, Pelagibacter, Prochlorococcus, Rhizobium, Synechococcus, and Thermus phage genomes that are known not to infect the bacterial model <br /> 
#create host model directory <br /> 
WIsH/WIsH -c build -g phage_hosts_genomes/ -m modelDir  <br /> 
#Run predict on null to get llikelihood.matrix in outputNullModelResultDir to feed into computeNullParameters.R  <br /> 
WIsH/WIsH -c predict -g null/ -m modelDir -r outputNullModelResultDir -b 1000   <br /> 
Rscript ../WIsH/computeNullParameters.R  <br /> 
WIsH/WIsH -c predict -g phage_genomes/ -m modelDir -r outputResultDir -b 1000 -n outputNullModelResultDir/nullParameters.tsv <br /> 
### PHIRBO
Download taxdb.btd, taxdb.bti, taxdb.tar.gz from https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz <br />
Prepare blast results via blast_forphirbo.sh <br />
python phirbo/phirbo.py phage_virusblast/ phage_hostsblast/ phage_phirbo/predictions.csv

### PHP
python3 PHP/countKmer.py -f phage_host_genomes -d phage_host_PHPkmer -n phage_PHPHostKmer -c -1 <br />
python3 PHP/PHP.py -v phage_genomes -o phage_PHPout -d phage_PHPkmer -n phage_PHPHostKmer <br />

### PHIST
python PHIST/phist.py phage_genome_dir phage_name phage_PHIRBO_outdir <br />
#results will be in out directory labelled as predictions.csv

### VirHostMatcher 
python /opt/VirHostMatcher-Net/VirHostMatcher-Net.py -q phage_genomes/ -o Phage_vhmn_output -n 1000 

# Exploratory Tools 
All exploratory tools require a docker image to be created before running the tools
apptainer build toolname.sif toolname.def <br /> 
apptainer shell toolname.sif <br /> 
cd /opt/miniconda/bin  <br /> 
. activate <br /> 

### vpf-class 
stack exec -- vpf-class --data-index data/index.yaml -i phage.fasta -o phage_test-classified

### vHULK
python vHULK.py -i phage_genomes/ -o phage_vhulk_outdir --all

### RaFAH
perl RaFAH.pl --predict --genomes_dir phage_genomes/ --extension .fasta --valid_ogs_file HP_Ranger_Model_3_Valid_Cols.txt --hmmer_db_file_name HP_Ranger_Model_3_Filtered_0.9_Valids.hmm --r_script_predict_file_name RaFAH_Predict_Host.R --r_model_file_name MMSeqs_Clusters_Ranger_Model_1+2+3_Clean.RData

### HostG
python run_Speed_up.py --contigs phage_genome.fasta --t 0 

### VirHostMatcher-Net
python /opt/VirHostMatcher-Net/VirHostMatcher-Net.py -q phage_genomes -o phage_vhmn_out -i tmp -n 3 -t 8 <br /> 

### CHERRY
python run_Speed_up.py --contigs phage_genome.fasta --model pretrain --topk 1000

