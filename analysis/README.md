# Analysis of ddRADSeq data for <i>Vittaria appalachiana</i>

This README details the methods and analysis of the ddRADSeq data for <i>Vittaria appalachiana</i> in Pelosi et al. (in prep). 

## Demulitplexing 

We used [stacks v2.64](https://catchenlab.life.illinois.edu/stacks/) to demultiplex reads. Samples were digested with nsiI and mspI at the University of Wisconsin - Madison genomics core. 
```
process_radtags -P -i fastq -b sample_barcodes.txt -o ./ -1 ../E-Sessa_Vittaria_S38_L003_R1_001.fastq \
-2 ../E-Sessa_Vittaria_S38_L003_R2_001.fastq --renz-1 nsiI --renz-2 mspI -c -q -r
```

Of the 458M reads, 444M were retained (96.9%). 

| File | Retained Reads | Low Quality | Barcode Not Found | RAD Cutsite Not Found | Total |
|------|----------------|-------------|-------------------|-----------------------|-------|
|Sessa_Vittaria_S38_L003_R1_001.fastq|444442231|64666|9517956|4505659|458530512|


## Plastome 
We used [RADOrgMiner v0.9](https://github.com/laczkol/RADOrgMiner) to map demulitplexed reads to the plastid genome of <i>Vittaria appalachiana</i>. The plastome assembly was downloaded from NCBI under accession NC_040219. 

```
RADOrgMiner.sh -r Vittaria.appalachiana_plastome.fa -np 12 -type PE -popmap PopMap.txt
RADOrgMiner.sh -r Vittaria.appalachiana_plastome.fa -np 24 -align no \
-call yes -minbc 15 -popmap PopMap.txt -type PE 
```

The resulting output files were used to analyze population genetics of the plastid in the R script `PlastidHaplotypes.R`. 

## Nuclear Analysis 

[iPyrad v0.9.53](https://ipyrad.readthedocs.io/en/master/) was used to de novo assemble RAD loci. 
```
ipyrad -p params-Subset_Clust95_dip.txt -c 24 -s 1234567 -f
```

To account for the peculularities of the biology of <i>V. appalachiania</i>, we mapped reads for each sample back to consensus sequences for each locus. We first converted the `.loci` output files from ipyrad to fasta files and then generated a consensus sequence for each locus. The python script [loci2fasta.py](https://github.com/btmartin721/file_converters/blob/master/loci2fasta.py) is from Bradley Martin. 
```
for file in *.loci; do python loci2fasta.py -L "$file";done
for file in *.fasta; do cons -sequence "$file" -outseq "$file".cons;done
cat *.cons > all_consensus.fasta
```
We then mapped reads for each sample to the consenus sequences. 
```
bwa index all_consensus.fasta
for file in $(cat samples.txt); do bwa mem "$file".R1.fq "$file".R2.fq > "$file".sam;done
for file in *.sam; do samtools sort "$file" -o "$file".sort;done
for file in *.sort; do gatk MarkDuplicates -I ${file} -M ${file}_deduplication_stats.txt -O ${file}_deduplicated.sam --REMOVE_DUPLICATES true;done
```
Generate mpileup and call variants with bcftools v1.1.9. 
```
bcftools mpileup -f all_consensus.fa -b samples.txt --threads 16 -a AD,DP > bcftools_allsites.mpileup
bcftools call bcftools_allsites.mpileup -m -f GQ -o allSites.vcf --ploidy 4 --threads 16
```
The resulting VCF was further filtered with [VCFtools v0.1.16](https://vcftools.github.io/index.html) to only retain SNPs that were present in at least 50% of samples with a minor allele frequency greater than 0.01 and minimum mean read depth of 4x. Two datasets were constructed: one that contains both variable and invariant sites and one that contains only variant sites. 
```
vcftools --vcf Subset_Clust95.vcf --out allSites.filtered.ma1.mm0.5.dp4 --recode --max-missing 0.5 --maf 0.01 --min-meanDP 4
```
Subsets of the original VCF were filtered with VCFtools to only include samples from within a colony to analyze the within-colony sequence data. 

The resulting VCFs were used to analyze nuclear population genetics in the R script `NuclearddRADAnalyses.R`. We also used pixy v1.2.10.beta2 to get estimates of pi, Fst, and Dxy for the variable + invariant data. 
```
pixy --vcf allSites.filtered.ma1.mm0.5.dp4.out.recode.vcf.gz --stats pi fst dxy --populations pops.pixy.txt --window_size 500 --n_cores 64 --output_prefix populations
pixy --vcf allSites.filtered.ma1.mm0.5.dp4.out.recode.vcf.gz --stats pi fst dxy --populations col.pixy.txt --window_size 500 --n_cores 64 --output_prefix colonies
```

