#!/bin/bash

# Run repeated trials of antic-herit-sim.r with different parameters to generate sim data for paper

echo -e "ngen\tlimit_n\tpyearmin\tpyearmax\tamode\trate\tamin\tamax\tantic\tanticp\therit\theritp\tyobslope\tyobslopep" > antic-herit-results.txt

# Table 1. 1000 trials of each scenario with limit_n = 100
for i in {1..1000}
do
    Rscript antic-herit-sim.r --seed $i --amode 1 --pyearmin 1700 --pyearmax 2000 --npairs 100000 -l 100 --amin 1989 --amax 2013 -q >> antic-herit-results.txt
    Rscript antic-herit-sim.r --seed $i --amode 1 --pyearmin 1700 --pyearmax 2000 --npairs 100000 -l 100 --amin 1950 --amax 2013 -q >> antic-herit-results.txt
    Rscript antic-herit-sim.r --seed $i --amode 1 --pyearmin 1700 --pyearmax 2000 --npairs 100000 -l 100 --amin 1880 --amax 2013 -q >> antic-herit-results.txt
    Rscript antic-herit-sim.r --seed $i --amode 4 --pyearmin 1700 --pyearmax 2000 --npairs 100000 -l 100 --amin 1989 --amax 2013 -q --rate .05 >> antic-herit-results.txt
    Rscript antic-herit-sim.r --seed $i --amode 4 --pyearmin 1700 --pyearmax 2000 --npairs 100000 -l 100 --amin 1950 --amax 2013 -q --rate .05 >> antic-herit-results.txt
    Rscript antic-herit-sim.r --seed $i --amode 4 --pyearmin 1700 --pyearmax 2000 --npairs 100000 -l 100 --amin 1880 --amax 2013 -q --rate .05 >> antic-herit-results.txt
    Rscript antic-herit-sim.r --seed $i --amode 4 --pyearmin 1700 --pyearmax 2000 --npairs 100000 -l 100 --amin 1989 --amax 2013 -q --rate .01 >> antic-herit-results.txt
    Rscript antic-herit-sim.r --seed $i --amode 4 --pyearmin 1700 --pyearmax 2000 --npairs 100000 -l 100 --amin 1950 --amax 2013 -q --rate .01 >> antic-herit-results.txt
    Rscript antic-herit-sim.r --seed $i --amode 4 --pyearmin 1700 --pyearmax 2000 --npairs 100000 -l 100 --amin 1880 --amax 2013 -q --rate .01 >> antic-herit-results.txt
    Rscript antic-herit-sim.r --seed $i --amode 2 --pyearmin 1700 --pyearmax 2000 --npairs 100000 -l 100 --amin 1989 --amax 2013 -q >> antic-herit-results.txt
done

# 1..2 at 1:38p, done 1:39p

# Table 2. Pocchiari stratifications with limit_n = 41 (their total n) or limit_n = 26 (what they had in Table 3 in the paper)
Rscript antic-herit-sim.r --seed 1 --amode 1 --pyearmin 1700 --pyearmax 2000 --npairs 100000 -l 26 --amin 1989 --amax 2013 -f > table2.txt

# need to re-write to spit out overall and stratified p val and antic years (9 pairs of values) for each simulation

# example of Simulation 1 for Results text
Rscript antic-herit-sim.r --seed 1 --amode 1 --pyearmin 1700 --pyearmax 2000 --npairs 100000 -l 100 --amin 1989 --amax 2013 -f | less -S

# another example of Simulation 1 with heritability being real - note that after ascertainment,
# even in a linear model with year of birth of included, parent age of onset is still a
# significant predictor of child age of onset
Rscript antic-herit-sim.r --seed 1 --amode 1 --pyearmin 1700 --pyearmax 2000 --npairs 100000 -l 100 --amin 1989 --amax 2013 -f --heritability_is_real | less -S


#!/bin/bash
cd /humgen/atgu1/fs03/eminikel/030e200k/analysis/2
echo -e "ngen\tlimit_n\tpyearmin\tpyearmax\tamode\trate\tamin\tamax\tantic\tanticp\therit\theritp\tyobslope\tyobslopep\therit_wyob\theritp_wyob" > antic-herit-results-h2.txt
# Table S1. like Table 1 but where age of onset is actually 100% heritable (50% from the affected parent)
for i in {1..1000}
do
    Rscript antic-herit-sim.r --seed $i --amode 1 --pyearmin 1700 --pyearmax 2000 --npairs 100000 -l 100 --amin 1989 --amax 2013 --heritability_is_real -q            >> antic-herit-results-h2.txt
    Rscript antic-herit-sim.r --seed $i --amode 1 --pyearmin 1700 --pyearmax 2000 --npairs 100000 -l 100 --amin 1950 --amax 2013 --heritability_is_real -q            >> antic-herit-results-h2.txt
    Rscript antic-herit-sim.r --seed $i --amode 1 --pyearmin 1700 --pyearmax 2000 --npairs 100000 -l 100 --amin 1880 --amax 2013 --heritability_is_real -q            >> antic-herit-results-h2.txt
    Rscript antic-herit-sim.r --seed $i --amode 4 --pyearmin 1700 --pyearmax 2000 --npairs 100000 -l 100 --amin 1989 --amax 2013 --heritability_is_real -q --rate .05 >> antic-herit-results-h2.txt
    Rscript antic-herit-sim.r --seed $i --amode 4 --pyearmin 1700 --pyearmax 2000 --npairs 100000 -l 100 --amin 1950 --amax 2013 --heritability_is_real -q --rate .05 >> antic-herit-results-h2.txt
    Rscript antic-herit-sim.r --seed $i --amode 4 --pyearmin 1700 --pyearmax 2000 --npairs 100000 -l 100 --amin 1880 --amax 2013 --heritability_is_real -q --rate .05 >> antic-herit-results-h2.txt
    Rscript antic-herit-sim.r --seed $i --amode 4 --pyearmin 1700 --pyearmax 2000 --npairs 100000 -l 100 --amin 1989 --amax 2013 --heritability_is_real -q --rate .01 >> antic-herit-results-h2.txt
    Rscript antic-herit-sim.r --seed $i --amode 4 --pyearmin 1700 --pyearmax 2000 --npairs 100000 -l 100 --amin 1950 --amax 2013 --heritability_is_real -q --rate .01 >> antic-herit-results-h2.txt
    Rscript antic-herit-sim.r --seed $i --amode 4 --pyearmin 1700 --pyearmax 2000 --npairs 100000 -l 100 --amin 1880 --amax 2013 --heritability_is_real -q --rate .01 >> antic-herit-results-h2.txt
    Rscript antic-herit-sim.r --seed $i --amode 2 --pyearmin 1700 --pyearmax 2000 --npairs 100000 -l 100 --amin 1989 --amax 2013 --heritability_is_real -q            >> antic-herit-results-h2.txt
done
# bsub -q priority -P $RANDOM -o tables1.out -e tables1.err -J tables1 ./tables1.bash

# re-generate images
Rscript antic-herit-sim.r --seed 1 --amode 1 --pyearmin 1700 --pyearmax 2000 --npairs 100000 -l 100 --amin 1989 --amax 2013 -f -i ~/d/sci/030e200k/analysis/2

