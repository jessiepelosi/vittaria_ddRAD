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

The resulting VCF was further filtered with [VCFtools v0.1.16](https://vcftools.github.io/index.html) to only retain bi-allelic SNPs that were present in at least 50% of samples with a minor allele frequency greater than 0.01:
```
vcftools --vcf Subset_Clust95.vcf --out Sub_Clust95_filt0.5_bi.recode.vcf --recode --min-alleles 1 --max-alleles 2 --max-missing 0.5 --maf 0.01
```
Subsets of the original VCF were filtered with VCFtools to only include samples from within a colony. The resulting VCFs were used to analyze nuclear population genetics in the R script `NuclearddRADAnalyses.R`. 

