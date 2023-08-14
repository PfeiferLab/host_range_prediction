#create host model directory 
WIsH/WIsH -c build -g phage_hosts_genomes/ -m modelDir 
#Run predict on null to get llikelihood.matrix in outputNullModelResultDir to feed into computeNullParameters.R 
WIsH/WIsH -c predict -g null/ -m modelDir -r outputNullModelResultDir -b 1000  
Rscript ../WIsH/computeNullParameters.R  
WIsH/WIsH -c predict -g phage_genomes/ -m modelDir -r outputResultDir -b 1000 -n outputNullModelResultDir/nullParameters.tsv 
