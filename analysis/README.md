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

Samples that were collected from the same colony were merged. <-- check if this remains at end. 

## Plastome 
We used [RADOrgMiner v0.9](https://github.com/laczkol/RADOrgMiner) to map demulitplexed reads to the plastid genome of <i>Vittaria appalachiana</i>. The plastome assembly was downloaded from NCBI under accession NC_040219. 

```
RADOrgMiner.sh -r Vittaria.appalachiana_plastome.fa -np 12 -type PE -popmap PopMap.txt
RADOrgMiner.sh -r Vittaria.appalachiana_plastome.fa -np 24 -align no \
-call yes -minbc 15 -popmap PopMap.txt -type PE 
```
A second iteration was run with RADOrgMiner with a subset of the samples (XX/52). 

## Nuclear Analysis 

