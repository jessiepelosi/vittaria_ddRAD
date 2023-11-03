# Simulations

This directory contains information regarding two sets of simluations for <i>Vittaria appalachiana</i>: demographic history (changes in effective population size and splits) and genetic diversity. 

## Demography 

To generate a site frequency spectrum (SFS) from our VCF, we used [easySFS](https://github.com/isaacovercast/easySFS).  

```
python easySFS.py -i Subset_Clust95.vcf -p single_pop.txt --total-length 40398599 --preview
python easySFS.py -i Subset_Clust95.vcf -p single_pop.txt --total-length 40398599 --proj XX -o 1pop_full_out_1snp-perloc
```
Three populations (as defined in our DAPC clustering analysis) were also used for projection to an SFS. 
```
python easySFS.py -i Subset_Clust95.vcf -p 3pop.txt --total-length 40398599 --preview
python easySFS.py -i Subset_Clust95.vcf -p 3pop.txt --total-length 40398599 --proj XX,XX,XX -o 3pop_full_out_1snp-perloc
```


<b>MUTATION RATE</b>

Changes in effective population size without a specified model was assessed with [stairway plot 2](https://github.com/xiaoming-liu/stairway-plot-v2). 
```
java -cp stairway_plot_es Stairbuilder vittaria.blueprint
bash vittaria.blueprint.sh
```

Several models of popluation demography splits were tested with [dadi v2.1.1](https://dadi.readthedocs.io/en/latest/). 


## Genetic Diversity 
<B>ELISSA'S STUFF HERE</B>
