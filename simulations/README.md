# Simulations

This directory contains information regarding two sets of simluations for <i>Vittaria appalachiana</i>: demographic history (changes in effective population size) and genetic diversity. 

## Demography 

### 1 Population 

To generate a site frequency spectrum (SFS) from our VCF, we used [easySFS](https://github.com/isaacovercast/easySFS).  

```
python easySFS.py -i Subset_Clust95.vcf -p single_pop.txt --total-length 40398599 --preview
python easySFS.py -i Subset_Clust95.vcf -p single_pop.txt --total-length 40398599 --proj 50 -o 1pop_full_out_1snp-perloc
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

### 3 Populations
To generate a site frequency spectrum (SFS) from our VCF, we used [easySFS](https://github.com/isaacovercast/easySFS).  

```
python easySFS.py -i Subset_Clust95.vcf -p 3pops.txt --total-length 40398599 --preview
python easySFS.py -i Subset_Clust95.vcf -p 3pops.txt --total-length 40398599 --proj 22,46,8 -o 3pop_full_out_1snp-perloc
```

XX three-dimensional models of popluation demography splits were tested with [dadi v2.1.1](https://dadi.readthedocs.io/en/latest/). We used scripts from the dadi_pipeline [three population pipeline](https://github.com/dportik/dadi_pipeline/tree/master/Three_Population_Pipeline) to optimize and run these models. 
```
clust2-clust1-clust3.sfs
```
The models that were tested are: 
- split_nomig (Split into three populations without migration)
- ancmig_adj_1 (Adjacent ancient migration, shortest isolation)
- ancmig_adj_2 (Adjacent ancient migration, shorter isolation)
- ancmig_adj_3 (Adjacent ancient migrations, longest isolation)
- sim_split_no_mig (Simultaneous split, no migration; from Barratt et al. 2018)
- sim_split_no_mig_size (Simulataneous split, no migration, size change; from Barratt et al. 2018)
- split_nomig_size (Divergence with no migration, size change; from Barratt et al. 2018)
- ancmig_2_size (Adjacent ancient migration, shorter isolation, size change; from Barratt et al. 2018) 

  
## Genetic Diversity 
<B>ELISSA'S STUFF HERE</B>
