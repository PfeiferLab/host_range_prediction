# Confirmatory Tools 

### Phirbo
https://github.com/aziele/phirbo <br />
Download taxdb.btd, taxdb.bti, taxdb.tar.gz from https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz <br />
Prepare rank-biased overlap blast results via `blast_forphirbo.sh` <br />
Run phribo by providing two input directories (i.e., for phages [`phage_virusblast/`] and bacteria [`phage_hostblast/`]) containing ranked lists from blast output, and an output file name (`phage_phirbo/predictions.csv`) <br />
`python phirbo/phirbo.py phage_virusblast/ phage_hostsblast/ phage_phirbo/predictions.csv`

### PHIST
https://github.com/refresh-bio/PHIST <br />
PHIST takes as input two directories containing FASTA files (gzipped or not) with genomic sequences of viruses (`phage_genome_dir`) and candidate hosts (`host_dir`) <br />
`python PHIST/phist.py phage_genome_dir host_dir phage_PHIST_outdir` <br />
Results will be in the output directory (`phage_PHIST_outdir`) labeled as `predictions.csv`

### Prokaryotic virus Host Predictor (PHP)
https://github.com/congyulu-bioinfo/PHP  <br />
First, calculate the K-mer frequency of the host <br />
`python3 PHP/countKmer.py -f phage_host_genomes -d phage_host_PHPkmer -n phage_PHPHostKmer -c -1` where: <br />
- `-f phage_host_genomes`: fasta file of prokaryotic genome sequence (_**one fasta per file**_)
- `-d phage_host_PHPkmer`: path of the prokaryotic k-mer file
- `-n phage_PHPHostKmer`: name of the prokaryotic k-mer file
- `-c -1`: the number of cores used; -1 represents the use of ALL cores

Then, predict the infection relationship between the virus and the host  <br />
`python3 PHP/PHP.py -v phage_genomes -o phage_PHPout -d phage_PHPkmer -n phage_PHPHostKmer` <br />
- `-v phage_genomes`: fasta file of query virus sequences (_**one virus genome per file**_)
- `-o phage_PHPout`: path of temp and result files
- `-d phage_PHPkmer`: path of the prokaryotic k-mer file
- `-n phage_PHPHostKmer`: name of the prokaryotic k-mer file

### VirHostMatcher 
https://github.com/jessieren/VirHostMatcher <br />
To run VHM create a folder containing virus fasta files and a folder containing host fasta files - no subfolders  <br />
`python VirHostMatcher/vhm.py -v phage_genomes/ -b phage_host_genomes/ -o Phage_vhm_output/` <br />
- `-v phage_genomes/`: path to virus folder
- `-b phage_host_genomes/`: path to host folder
- `-o Phage_vhm_output/`: path to output

### WIsH 
https://github.com/soedinglab/WIsH <br /> 
Null directory is created from a diverse range of Alteromonas, Cellulophage, Cyanophage, Lactobacillus, Mycobacterium, Oenococcus, Pelagibacter, Prochlorococcus, Rhizobium, Synechococcus, and Thermus phage genomes that are known not to infect the bacterial model listed in null.txt <br /> 
Create host model directory <br /> 
`WIsH/WIsH -c build -g phage_hosts_genomes/ -m modelDir`  <br /> 
- `-c build`: tell WIsH to create a model
- `-g phage_hosts_genomes/`: where bacterial genomes in FASTA format are stored
- `-m modelDir`: where model for each bacterial genome will be stored

Run predict on null to get llikelihood.matrix in outputNullModelResultDir to feed into computeNullParameters.R  <br /> 
```
WIsH/WIsH -c predict -g null/ -m modelDir -r outputNullModelResultDir -b 1000
Rscript ../WIsH/computeNullParameters.R
WIsH/WIsH -c predict -g phage_genomes/ -m modelDir -r outputResultDir -b 1000 -n outputNullModelResultDir/nullParameters.tsv
```
- `-c predict`: tell WIsH to run a prediction
- `-g null/` & `-g phage_genomes/`: where bacterial genomes in FASTA format are stored
- `-m modelDir`: where model for each bacterial genome will be stored
- `-r outputResultDir`: where results will be stored
- `-b 1000`: k option to output the k best prediction by likelihood
- `-n outputNullModelResultDir/nullParameters.tsv`: location of null model; to be used when p-values are desired

# Exploratory Tools 
All exploratory tools require a docker image to be created before running the tools <br /> 
```
apptainer build toolname.sif toolname.def
apptainer shell toolname.sif
cd /opt/miniconda/bin
. activate
```

### CHERRY
https://github.com/KennthShang/CHERRY <br /> 
Prepare the database by unpacking in install location via the command: `apptainer shell cherry.sif` <br /> 
Predict hosts for viruses using a fasta file containing the viral sequences as the input <br />  
`python run_Speed_up.py --contigs phage_genome.fasta --model pretrain --topk 1000` <br />
- `--contigs`: input fasta file
- `--model`: predicting host with pretrained (or retrained) parameters (default: pretrained)
- `--topk`: host prediction with topk score (default: 1)

### HostG
https://github.com/KennthShang/HostG  <br /> 
Prepare the database by unpacking in install location via the command: `apptainer shell hostg.sif` <br /> 
Predict hosts for viruses using a fasta file containing the viral sequences as the input <br />  
`python run_Speed_up.py --contigs phage_genome.fasta --t 0`
- `--contigs`: path to the contigs file
- `--t`: threshold value; predictions will only be outputted that have confidence (SoftMax) values above t (default: 0)

### Random Forest Assignment of Hosts (RaFAH)
https://github.com/felipehcoutinho/RaFAH <br />  
To setup RaFAH download and decompress all the files listed in the RaFAH github repositroy <br />  
Perform host predictions using a set of pre-computed model files <br />  
```
perl RaFAH.pl --predict --genomes_dir phage_genomes/ --extension .fasta --valid_ogs_file HP_Ranger_Model_3_Valid_Cols.txt --hmmer_db_file_name HP_Ranger_Model_3_Filtered_0.9_Valids.hmm --r_script_predict_file_name RaFAH_Predict_Host.R --r_model_file_name MMSeqs_Clusters_Ranger_Model_1+2+3_Clean.RData
```
- `--genomes_dir`: directory containing DNA sequences in FASTA format
- `--extension`: extension of the FASTA files to be analyzed (e.g., fasta, fna, fa; default: fasta)
- `--valid_ogs_file`: list of valid OGs used by a specific model file
- `--hmmer_db_file_name`: output file generated by querying Hmmer DB files against identified/provided protein sequences using hmmsearch
- `--r_script_predict_file_name`: R script used to perform predictions (.R extension)
- `--r_model_file_name`: R data file containing the pre-computed RF for a specific model (.RData extension)

### viral Host UnveiLing Kit (vHULK)
https://github.com/LaboratorioBioinformatica/vHULK <br />  
Run vHULK on target genera <br />  
`python vHULK.py -i phage_genomes/ -o phage_vhulk_outdir --all`
- `-i`: input directory; path to a folder containing metagenomic bins in .fa or .fasta format
- `-o`: output directory; location to store results in -- will be created if absent
- `--all`: write predictions for all input bins/genomes, even if they were skipped (i.e., size filtered or hmmscan failed)

### VirHostMatcher-Net
https://github.com/WeiliWw/VirHostMatcher-Net <br /> 
Install dependencies, build local modules, download required databases via the command: `apptainer shell vhmn.sif` <br /> 
Run VirHostMatcher-Net <br /> 
`python /opt/VirHostMatcher-Net/VirHostMatcher-Net.py -q phage_genomes -o phage_vhmn_out -i tmp -n 3 -t 8` <br /> 
- `-q`: directory containing query virus genomes with .fasta or .fa extension
- `-o`: output directory
- `-i`: directory storing intermediate results
- `-n`: topN; number of top predictions written to the output files (all predictions will be output if there is a tie in score)
- `-t`: number of threads to use

### VPF-Class 
https://github.com/biocom-uib/vpf-tools <br />  
Clone the repository and compile all targets using the command: `apptainer shell vhulk.sif` <br /> 
Run vHULK <br /> 
`stack exec -- vpf-class --data-index data/index.yaml -i phage.fasta -o phage_test-classified`
- `--data-index`: file that specifies classification levels
- `-i`: input file
- `-o`: output directory

# Analysis and Plotting 

Figures 1A-D and Table 1 were created with `dotplots_ani_combined_formanuscript.R` <br />
Figure 2, 3A-B were created with `upsetRforcomparitivepaper.R` <br />

ANI data for heatmaps was generated by `taxon_for_ANIs.R` - which extracts the desired species names and their GCF numbers from the list of those available in the exploratory tools databases. The script taxons_forANIs.R produces an outfile that is directed into the command: <br /> 
`bit-dl-ncbi-assemblies -w ncbi_accessions.txt -f fasta`   <br /> 
The status of all genomes for this analysis was evaluted with: <br /> 
```
datasets summary genome accession $line --as-json-lines | \
dataformat tsv genome --fields accession,assminfo-status --elide-header
```

The results fastas are then combined and used in the command:   <br /> 
`anvi-script-compute-ani-for-fasta -f combined.fasta -o ani_output --method ANIb -T 10`  <br />  

Supplementary Figure S1 was created with `correlations.R` <br />
Supplementary Figures 2, 3B were created with `taxon_for_ANIs.R` <br />
Supplementary Figure 3A was created with `taxon_for_ANIs_localtemp.R` <br /> 

Supplementary Figures 4A,B were created with `ANIs.R` <br /> 

KFS-EC3 is able to infect _E. coli_ ATCC 10536 but not _E. coli_ 15144, _E. coli_ 2192, or _E. coli_ 2196.  <br />
Genomes for these _E. coli_ were filtered for the longest contig with `long_script.py` and then aligned to each other using blastn and annotated using prokka. <br />
`synteny.R` returns the regions present in 10536 and absent in all others and `synteny_2.R` returns the regions absent in 10536 but shared across all non-infectable bacteria.
