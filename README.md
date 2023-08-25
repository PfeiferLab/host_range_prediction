# Confirmatory Tools 

### Phirbo
Download taxdb.btd, taxdb.bti, taxdb.tar.gz from https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz <br />
Prepare rank-biased overlap blast results via blast_forphirbo.sh <br />
Run phribo by providing two input directories (for phages and bacteria) containing ranked lists from blast output, and an output file name 
python phirbo/phirbo.py phage_virusblast/ phage_hostsblast/ phage_phirbo/predictions.csv

### PHIST
PHIST takes as input two directories containing FASTA files (gzipped or not) with genomic sequences of viruses and candidate hosts  <br />
python PHIST/phist.py phage_genome_dir phage_name phage_PHIRBO_outdir <br />
results will be in out directory labelled as predictions.csv

### Prokaryotic virus Host Predictor (PHP)
First calculate the K-mer frequency of the host <br />
python3 PHP/countKmer.py -f phage_host_genomes -d phage_host_PHPkmer -n phage_PHPHostKmer -c -1 <br />
Then predict the infection relationship between the virus and the host
python3 PHP/PHP.py -v phage_genomes -o phage_PHPout -d phage_PHPkmer -n phage_PHPHostKmer <br />

### VirHostMatcher 
To run VHM create a folder containing virus fasta files and a folder containing host fasta files - no subfolders 
python /opt/VirHostMatcher-Net/VirHostMatcher-Net.py -q phage_genomes/ -o Phage_vhmn_output -n 1000 

### WIsH 
Null directory is created from a diverse range of Alteromonas, Cellulophage, Cyanophage, Lactobacillus, Mycobacterium, Oenococcus, Pelagibacter, Prochlorococcus, Rhizobium, Synechococcus, and Thermus phage genomes that are known not to infect the bacterial model listed in null.txt <br /> 
Create host model directory <br /> 
WIsH/WIsH -c build -g phage_hosts_genomes/ -m modelDir  <br /> 
Run predict on null to get llikelihood.matrix in outputNullModelResultDir to feed into computeNullParameters.R  <br /> 
WIsH/WIsH -c predict -g null/ -m modelDir -r outputNullModelResultDir -b 1000   <br /> 
Rscript ../WIsH/computeNullParameters.R  <br /> 
WIsH/WIsH -c predict -g phage_genomes/ -m modelDir -r outputResultDir -b 1000 -n outputNullModelResultDir/nullParameters.tsv <br /> 


# Exploratory Tools 
All exploratory tools require a docker image to be created before running the tools <br /> 
apptainer build toolname.sif toolname.def <br /> 
apptainer shell toolname.sif <br /> 
cd /opt/miniconda/bin  <br /> 
. activate <br /> 

### CHERRY
python run_Speed_up.py --contigs phage_genome.fasta --model pretrain --topk 1000

### HostG
python run_Speed_up.py --contigs phage_genome.fasta --t 0 

### Random Forest Assignment of Hosts (RaFAH)
perl RaFAH.pl --predict --genomes_dir phage_genomes/ --extension .fasta --valid_ogs_file HP_Ranger_Model_3_Valid_Cols.txt --hmmer_db_file_name HP_Ranger_Model_3_Filtered_0.9_Valids.hmm --r_script_predict_file_name RaFAH_Predict_Host.R --r_model_file_name MMSeqs_Clusters_Ranger_Model_1+2+3_Clean.RData

### viral Host UnveiLing Kit (vHULK)
python vHULK.py -i phage_genomes/ -o phage_vhulk_outdir --all

### VirHostMatcher-Net
python /opt/VirHostMatcher-Net/VirHostMatcher-Net.py -q phage_genomes -o phage_vhmn_out -i tmp -n 3 -t 8 <br /> 

### VPF-Class 
stack exec -- vpf-class --data-index data/index.yaml -i phage.fasta -o phage_test-classified


# Analysis and Plotting 

Figures 1A-D and Table 1 were created with dotplots_ani_combined_formanuscript.R <br />
Figure 2, 3A-B were created with upsetRforcomparitivepaper.R <br />

ANI data for heatmaps was generated by taxon_for_ANIs.R - which extracts the desired species names and their GCF numbers from the list of those available in the exploratory tools databases. The script taxons_forANIs.R produces an outfile that is directed into the command: <br /> 
bit-dl-ncbi-assemblies -w ncbi_accessions.txt -f fasta   <br /> 
The status of all genomes for this analysis was evaluted with: <br /> 
datasets summary genome accession $line --as-json-lines | \
dataformat tsv genome --fields accession,assminfo-status --elide-header 

The results fastas are then combined and used in the command:   <br /> 
anvi-script-compute-ani-for-fasta -f combined.fasta -o ani_output --method ANIb -T 10  <br />  

Supplementary Figure S1 was created with correlations.R <br />
Supplementary Figures 2, 3B were created with taxon_for_ANIs.R <br />
Supplementary Figure 3A was created with taxon_for_ANIs_localtemp.R <br /> 

Supplementary Figures 4A,B were created with ANIs.R <br /> 

KFS-EC3 is able to infect E. coli ATCC 10536 but not E. coli 15144, E. coli 2192, or E. coli 2196.  <br />
Genomes for these E. coli were filtered for the longest contig with long_script.py and then aligned to each other using blastn and annotated using prokka. <br />
synteny.R returns the regions present in 10536 and absent in all others and and synteny_2.R is returns the regions absent in 10536 but shared across all non-infectable bacteria
