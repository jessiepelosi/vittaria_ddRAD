# Simulations

This directory contains information regarding two sets of simluations for <i>Vittaria appalachiana</i>: demographic history (changes in effective population size) and genetic diversity. 

## Demography with 1 Population Models 

To generate a site frequency spectrum (SFS) from our VCF, we used [easySFS](https://github.com/isaacovercast/easySFS).  

```
python easySFS.py -i variantSites.filtered.ma1.mm0.5.dp4.thinned.recode.vcf -p single_pop.txt --total-length 1609347 --preview
python easySFS.py -i variantSites.filtered.ma1.mm0.5.dp4.thinned.recode.vcf -p single_pop.txt --total-length 1609347 --proj 146 -o variantSites_1pop_1snpperloc
```

Changes in effective population size without a specified model was assessed with [stairway plot 2](https://github.com/xiaoming-liu/stairway-plot-v2). 
```
java -cp stairway_plot_es Stairbuilder vittaria.blueprint
bash vittaria.blueprint.sh
```

Six one-dimensional models of popluation demography splits were tested with [dadi v2.1.1](https://dadi.readthedocs.io/en/latest/). We used a custom python scripts to run each model with optimized parameters (see `dadi_pipeline_cmds.py`).  

The models that were tested are: 
- Neutral (neutral; no change in population size)
- Two Epoch (two_epoch; one instantaneous change in population size)
- Exponential Growth (exponential_growth; population grows exponentially one time in past)
- Bottlegrowth (bottlegrowth: population undergoes bottleneck followed by exponential growth)
- Three Epoch (three_epoch; two instantaneous changes in population size)
- Three Epoch with Inbreeding (three_epoch_inbreeding; two instantaneous changes in population size with inbreeding)

