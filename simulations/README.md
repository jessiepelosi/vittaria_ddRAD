# Simulations

This directory contains information regarding two sets of simluations for <i>Vittaria appalachiana</i>: demographic history (changes in effective population size) and genetic diversity. 

## Demography 

To generate a site frequency spectrum (SFS) from our VCF, we used [easySFS](https://github.com/isaacovercast/easySFS).  

```
python easySFS.py -i Subset_Clust95.vcf -p single_pop.txt --total-length 40398599 --preview
python easySFS.py -i Subset_Clust95.vcf -p single_pop.txt --total-length 40398599 --proj XX -o 1pop_full_out_1snp-perloc
```

Changes in effective population size without a specified model was assessed with [stairway plot 2](https://github.com/xiaoming-liu/stairway-plot-v2). 
```
java -cp stairway_plot_es Stairbuilder vittaria.blueprint
bash vittaria.blueprint.sh
```

Six one-dimensional models of popluation demography splits were tested with [dadi v2.1.1](https://dadi.readthedocs.io/en/latest/). We used custom python and bash scripts to run each model with optimized parameters 100 times. 
```
for i in {1..100}; do bash RunDadi.sh -i clust1-50.sfs -r 5 -u 4.97e-9 -l 40398599 -o test -p $i -m [model];done
```
The models that were tested are: 
- Neutral (neutral; no change in population size)
- Two Epoch (two_epoch; one instantaneous change in population size)
- Exponential Growth (exponential_growth; population grows exponentially one time in past)
- Bottlegrowth (bottlegrowth: population undergoes bottleneck followed by exponential growth)
- Three Epoch (three_epoch; two instantaneous changes in population size)
- Three Epoch with Inbreeding (three_epoch_inbreeding; two instantaneous changes in population size with inbreeding) 

## Genetic Diversity 
<B>ELISSA'S STUFF HERE</B>
