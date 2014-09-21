# aggregate-simulations.r
setwd("~/d/sci/030e200k/analysis/2")

require(sqldf)
options(stringsAsFactors=FALSE)

# Table 1
sim_results = read.table('antic-herit-results.txt',header=TRUE,sep='\t')
colnames(sim_results)[14] = "yobslopep" # temp patch

table1 = sqldf("
select   ngen, limit_n, pyearmin, pyearmax, amode, rate, amin, amax, -- that's 1-8
         avg(antic) mean_antic,
         cast(sum(case when anticp < .05 then 1 else 0 end) as float) / cast(count(*) as float) antic_sig,
         avg(herit) mean_herit,
         cast(sum(case when heritp < .05 then 1 else 0 end) as float) / cast(count(*) as float) herit_sig,
         avg(yobslope) mean_yobslope,
         cast(sum(case when yobslopep < .05 then 1 else 0 end) as float) / cast(count(*) as float) yobslope_sig
from     sim_results
group by 1,2,3,4,5,6,7,8
order by amode asc, rate desc, amin desc
;")

write.table(table1,"table1.txt",row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")

# Table S1
sim_results_s1 = read.table('antic-herit-results-h.txt',header=TRUE,sep='\t')

tables1 = sqldf("
select   ngen, limit_n, pyearmin, pyearmax, amode, rate, amin, amax, -- that's 1-8
         avg(antic) mean_antic,
         cast(sum(case when anticp < .05 then 1 else 0 end) as float) / cast(count(*) as float) antic_sig,
         avg(herit) mean_herit,
         cast(sum(case when heritp < .05 then 1 else 0 end) as float) / cast(count(*) as float) herit_sig,
         avg(yobslope) mean_yobslope,
         cast(sum(case when yobslopep < .05 then 1 else 0 end) as float) / cast(count(*) as float) yobslope_sig,
         avg(herit_wyob) mean_herit_wyob,
         cast(sum(case when heritp_wyob < .05 then 1 else 0 end) as float) / cast(count(*) as float) herit_wyob_sig
from     sim_results_s1
group by 1,2,3,4,5,6,7,8
order by amode asc, rate desc, amin desc
;")

write.table(tables1,"tables1.txt",row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")