setwd('~/d/sci/030e200k/analysis/2')
options(stringsAsFactors=FALSE)
require(sqldf)
require(survival)

# colors for plotting parents and children throughout paper
pcolor='#D18532' 
ccolor='#091090' 

# source datasets
mrc = read.table('mrc.txt',sep='\t',header=TRUE)
ucsf = read.table('ucsf.txt',sep='\t',header=TRUE)
goett = read.table('goett.txt',sep='\t',header=TRUE)
ancjdr = read.table('ancjdr.txt',sep='\t',header=TRUE)

# combine datasets
phen = rbind(mrc,ucsf,goett,ancjdr)

# schema of dataset
colnames(phen)
# [1] "study"              "fid"                "iid"                "affparent"          "father"            
# [6] "mother"             "sex"                "mutation_in_family" "gene_positive"      "gt129"             
# [11] "cis129"             "trans129"           "yob"                "age"                "censored"          
# [16] "dead"               "died_cjd"           "conf_non_cjd"       "duration_d"         "infcat"            
# [21] "notes"              "questionable"  

# check overall dataset size
dim(phen) # 513 individuals including unaffected relatives and individuals of unknown mutation status
sum(phen$gene_positive == 1) # 220 gene-positive individuals
sum(phen$gene_positive == 1 & phen$died_cjd) # 
sum((!is.na(phen$duration_d) | !is.na(phen$age) | !is.na(phen$gt129)) & phen$gene_positive == 1) 
# Of the 220, there are 212 about whom at least one key variable is known
# (the other 8 only contribute to genetic testing rate)
# size of dataset for each interesting variable:
sum((!is.na(phen$duration_d)) & phen$gene_positive == 1) 
# 60 gene-positive individuals for whom disease duration is known
sum((!is.na(phen$age)) & phen$gene_positive == 1) 
# 207 gene-positive individuals for whom age (either age at onset/death or censored age) is known
sum((!is.na(phen$gt129)) & phen$gene_positive == 1) 
# 108 gene-positive individuals for whom codon 129 genotype is known
sum((!is.na(phen$cis129)) & phen$gene_positive == 1) 
# 164 gene-positive individuals for whom the cis codon 129 is known
sum(( !is.na(phen$age) & !is.na(phen$gt129)) & phen$gene_positive == 1) 
# 104 gene-positive individuals for whom age and codon 129 genotype are known
sum((!is.na(phen$duration_d) & !is.na(phen$gt129)) & phen$gene_positive == 1) 
# 42 gene-positive individuals for whom duration and full codon 129 are known
sum((!is.na(phen$duration_d) & !is.na(phen$age) ) & phen$gene_positive == 1) 
# 59 gene-positive individuals for whom duration and age are known

# a few individuals were flagged as "questionable" if we were not absolutely certain
# about age at onset, diagnosis or paternity
sum(phen$questionable & phen$died_cjd) # of people who died of CJD
sum(phen$questionable & phen$gene_positive == 1) # of mutation-positive (including asymptomatics)
sum(phen$questionable) # of all individuals including unaffecteds

# uncomment this line to re-run without the "questionable" individuals
# phen = phen[!phen$questionable,]

# size and general characteristics of datasets for Table 3
sql_query = "
select   study,
         count(*) n_individuals,
         sum(gene_positive = 1) n_gene_positive,
         sum(case when gene_positive = 1 and infcat > 1 then 1 else 0 end) n_indirect,
         sum(case when died_cjd and age is not null then 1 else 0 end) n_ad,
         sum(case when died_cjd and age is not null and yob is not null then 1 else 0 end) n_ad_yob,
         sum(gene_positive = 1 and not died_cjd) n_asymptomatic
from     phen
group by study
order by study
;
"
descstats = sqldf(sql_query)
descstats$pct_indirect = descstats$n_indirect/descstats$n_gene_positive
descstats$pct_asymp = descstats$n_asymptomatic/descstats$n_gene_positive
descstats

# range of dates for direct ascertainment
sql_query = "
select   study, min(yob+age) fromyear, max(yob+age) toyear
from     phen
where    infcat =  1
and      died_cjd
group by study
order by study
;"
sqldf(sql_query)

# simple histogram of age of onset/death
hist(phen$age[phen$died_cjd],col='black',breaks=100)

# pub-quality histogram of age of onset/death (Fig 2a)
aa_hist_all = table(phen$age[phen$died_cjd])
aa_hist_both1 = table(phen$age[phen$died_cjd & phen$infcat == 1])
png('fig2a.aa.hist.e200k.png',res=300,pointsize=3,width=500,height=500)
par(lend=1)
plot(aa_hist_all,lwd=2,type='h',col='#999999',xaxt='n',xlab='Age of onset or death',ylab='Frequency',
     main = 'Age of onset/death distribution')
points(aa_hist_both1,col='black',type='h',lwd=2)
axis(side=1,at=c(30,40,50,60,70,80,90),labels=c(30,40,50,60,70,80,90))
abline(h=0,lwd=.5)
legend('topleft',c('direct','indirect'),lwd=2,col=c('black','#999999'))
mtext(side=3, line=-2, text="A", cex=2, adj=0, outer=TRUE)
dev.off()

# pub-quality histogram of year of onset/death (Fig 2b)
phen$yoa = phen$yob + phen$age
yoa_hist_all = table(phen$yoa[phen$died_cjd])
yoa_hist_both1 = table(phen$yoa[phen$died_cjd & phen$infcat == 1])
png('fig2b.yoa.hist.e200k.png',res=300,pointsize=3,width=500,height=500)
par(lend=1)
plot(yoa_hist_all,lwd=2,type='h',col='#999999',xaxt='n',xlab='Year of onset or death',ylab='Frequency',
     main = 'Year of onset/death distribution')
points(yoa_hist_both1,col='black',type='h',lwd=2)
axis(side=1,at=c(1940,1960,1980,2000,2013),labels=c(1940,1960,1980,2000,2013))
abline(h=0,lwd=.5)
legend('topleft',c('direct','indirect'),lwd=2,col=c('black','#999999'))
mtext(side=3, line=-2, text="B", cex=2, adj=0, outer=TRUE)
dev.off()

# check percentage of years of onset that are since the gene's discovery
sum(!is.na(phen$yoa) & phen$yoa > 1990)/sum(!is.na(phen$yoa)) # 92%
sum(!is.na(phen$yoa) & phen$yoa > 1988)/sum(!is.na(phen$yoa)) # 92%
# 92% are since 1991 inclusive, with no cases in 1989 or 1990

# is age of onset/death normally distributed?
shapiro.test(phen$age[phen$died_cjd])

# anova: does age of onset differ by study?
m = lm(age ~ study, data = subset(phen, died_cjd))
summary(m)
n = sum(phen$died_cjd & !is.na(phen$age) & !is.na(phen$study))
n 

# anova: does age of onset differ by sex?
m = lm(age ~ sex, data = subset(phen, died_cjd))
summary(m)
n = sum(phen$died_cjd & !is.na(phen$age) & !is.na(phen$sex))
n

# anova: does age of onset differ by mode of ascertainment?
phen$direct = phen$infcat==1
m = lm(age ~ direct, data = subset(phen, died_cjd))
summary(m)
n = sum(phen$died_cjd & !is.na(phen$age) & !is.na(phen$infcat))
n

# age of onset desc stats
mean(phen$age[phen$died_cjd],na.rm=TRUE)
sd(phen$age[phen$died_cjd],na.rm=TRUE)
sum(!is.na(phen$age[phen$died_cjd]))
median(phen$age[phen$died_cjd],na.rm=TRUE)

# overall survival curve for everyone
survivaldata = phen[phen$gene_positive == 1,c("age","died_cjd")]
mfit = survfit(Surv(age, died_cjd==1) ~ 1, data = survivaldata)
mfit # median 64, n = 208 observations, 159 events
# what does the survival curve look like overall?
plot(mfit,col='black',lwd=3,main='Survival of all E200K individuals',xlab='Age',ylab='Proportion surviving')
mfit$surv[mfit$time==80] # proportion surviving at age 80


# what are the sexes of the affected parents?
sql_query = "
select   parent.sex
from     phen parent, phen child
where    child.affparent = parent.iid
and      child.gene_positive = 1 
and      parent.gene_positive = 1
;
"
table(sqldf(sql_query))
# 47 F, 39 M

# anova: does age of onset differ by cis codon?
# linear model
m = lm(age ~ cis129, data = subset(phen, died_cjd))
summary(m)
# t test for same
t.test(age ~ cis129, data = subset(phen, died_cjd))
# n
sum(!is.na(phen$cis129) & !is.na(phen$age) & phen$died_cjd) 
sum(!is.na(phen$cis129) & !is.na(phen$age) & phen$died_cjd & phen$cis129 == 'M') 
sum(!is.na(phen$cis129) & !is.na(phen$age) & phen$died_cjd & phen$cis129 == 'V') 

# anova: does age of onset for cis 129M individuals differ by 129 genotype (i.e. trans codon)?
m = lm(age ~ gt129, data = subset(phen, died_cjd & cis129=='M'))
summary(m)
# t test for same
t.test(age ~ gt129, data = subset(phen, died_cjd & cis129=='M'))
# n
sum(phen$died_cjd & phen$cis129=='M' & phen$gt129=='MV' & !is.na(phen$gt129) & !is.na(phen$age),na.rm=TRUE)
sum(phen$died_cjd & phen$cis129=='M' & phen$gt129=='MM' & !is.na(phen$gt129) & !is.na(phen$age),na.rm=TRUE)

# cis 129M, MV vs MM survival analysis
# png('cx.129M.gt129.survival.analysis.png',width=600,height=400)
cism = subset(phen, cis129 == 'M' & gene_positive == 1 & !is.na(age))
survdata = data.frame(ages=cism$age, status=as.integer(!cism$censored), group=cism$gt129)
n = dim(survdata)[1]
msurv = with(survdata, Surv(ages, status==1))
diffobj = survdiff(Surv(ages,status==1)~group, data=survdata)
n = as.integer(sum(diffobj$n))
p = 1-pchisq(diffobj$chisq,df=1)
main = 'Survival of MM vs. MV genotypes with E200K cis 129M'
tobj = t.test(age ~ gt129, data = subset(phen, died_cjd & cis129=='M'))
p_t = tobj$p.value
subt = paste('n = ',n,' individuals, log-rank test p value = ',formatC(p,digits=2),
             ' t test p value = ',formatC(p_t,digits=2),sep='') 
plot(survfit(Surv(ages,status)~group,data=survdata),col=c('black','blue'),main=main,sub=subt,lwd=3)
legend('bottomleft',c('MM','MV'),col=c('black','blue'),lwd=3)
survfit(Surv(ages,status)~group,data=survdata)
# dev.off()
# note: there are 17 MV and 65 MM in the above analysis, which is slightly off from 
# expectation given the approx. 70%/30% allele ratio in Europeans
# however, there is bias here, as many individuals were not phased or not genotyped
# and for those I was only able to assign cis codons if at least one individual in the family
# was homozygous at codon 129. therefore 129 hets are underrepresented in this analysis.

# cis 129M vs cis 129V survival analysis
# png('cx.cis129.survival.analysis.png',width=600,height=400)
cisdata = subset(phen, gene_positive == 1 & !is.na(age) & !is.na(cis129))
survdata = data.frame(ages=cisdata$age, status=as.integer(!cisdata$censored), group=cisdata$cis129)
msurv = with(survdata, Surv(ages, status==1))
diffobj = survdiff(Surv(ages,status==1)~group, data=survdata)
n = as.integer(sum(diffobj$n))
p = 1-pchisq(diffobj$chisq,df=1)
tobj = t.test(age ~ cis129, data = subset(phen, died_cjd))
p_t = tobj$p.value
main = 'Survival of cis 129M vs. cis 129V haplotypes'
subt = paste('n = ',n,' individuals, log-rank test p value = ',formatC(p,digits=2),
             ' t test p value = ',formatC(p_t,digits=2),sep='') 
plot(survfit(Surv(ages,status)~group,data=survdata),col=c('black','blue'),main=main,sub=subt,lwd=3)
legend('bottomleft',c('cis 129M','cis 129V'),col=c('black','blue'),lwd=3)
survfit(Surv(ages,status)~group,data=survdata)
# dev.off()
# 145 M vs. 13 V

# to test above code, you can use this code above to add fake data to it and see how the MV curve changes
# testdata = data.frame(ages=as.numeric(c(30,29,29,27,26,25,24,23,22,21)), status=rep(1,10), group=rep("MV",10))
# survdata = rbind(survdata, testdata)

# pub-quality figure: yob/ao correlation (Fig 2c)
png('fig2c.yob.aa.png',res=300,pointsize=3,width=500,height=500)
m = lm(age ~ yob, data=subset(phen, died_cjd & infcat == 1))
summary(m)
slope1 = summary(m)$coefficients[2,1]
pval1 = summary(m)$coefficients[2,4]
m1 = m
m = lm(age ~ yob, data=subset(phen, died_cjd))
summary(m)
slope2 = summary(m)$coefficients[2,1]
pval2 = summary(m)$coefficients[2,4]
m2 = m
subtitle = ''
# uncomment this code to auto-generate a subtitle:
# subtitle=paste("slope = ",formatC(slope1,digits=2,format='f')," (p = ",
#                formatC(pval1,digits=2),") with only black points",", slope = ",
#                formatC(slope2,digits=2,format='f')," (p = ",
#                formatC(pval2,digits=2),") with all points",
#                sep='')
plot(phen$yob[phen$died_cjd], phen$age[phen$died_cjd], pch=19, col='#999999', xlab='Year of birth', ylab='Age of onset or death',
     main = 'Year of birth - age of onset/death correlation', sub = subtitle, ylim=c(0,90))
points(phen$yob[phen$died_cjd & phen$infcat == 1], phen$age[phen$died_cjd & phen$infcat == 1], pch=19, col='black')
abline(m1, col='black', lwd = .75)
abline(m2, col='black', lty='dashed', lwd= .75)
legend('bottomleft',c('direct','indirect'),col=c('black','gray'),pch=19)
mtext(side=3, line=-2, text="C", cex=2, adj=0, outer=TRUE)
dev.off()

slope1
pval1
slope2
pval2

# duration
sum(!is.na(phen$duration_d))
mean(phen$duration_d,na.rm=TRUE)
sd(phen$duration_d,na.rm=TRUE) 
summary(phen$duration_d)
median(phen$duration_d,na.rm=TRUE) 


# shape of duration distribution
hist(phen$duration_d)

# MM vs. MV duration for cis129M:
t.test(phen$duration_d[phen$gt129=='MM' & phen$cis129=='M'], phen$duration_d[phen$gt129=='MV' & phen$cis129=='M'])
# n
sum(!is.na(phen$duration_d[phen$gt129=='MM' & phen$cis129=='M'])) 
sum(!is.na(phen$duration_d[phen$gt129=='MV' & phen$cis129=='M'])) 

ks.test(phen$duration_d[phen$cis129=='M'], phen$duration_d[phen$cis129=='V']) 
wilcox.test(phen$duration_d[phen$cis129=='M'], phen$duration_d[phen$cis129=='V']) 
kruskal.test(phen$duration_d[!is.na(phen$cis129)], as.factor(phen$cis129[!is.na(phen$cis129)])) 

median(phen$duration_d[phen$cis129 == 'M' ],na.rm=TRUE) 
median(phen$duration_d[phen$cis129 == 'V' ],na.rm=TRUE)
sum(!is.na(phen$duration_d[phen$cis129=='M'])) 
sum(!is.na(phen$duration_d[phen$cis129=='V'])) 
sum(!is.na(phen$duration_d[phen$gt129=='VV'])) 

mean(phen$duration_d[phen$cis129 == 'M' & phen$gt129 == 'MM'],na.rm=TRUE) 
mean(phen$duration_d[phen$cis129 == 'M' & phen$gt129 == 'MV'],na.rm=TRUE) 
median(phen$duration_d[phen$cis129 == 'M' & phen$gt129 == 'MM'],na.rm=TRUE) 
sum(!is.na(phen$duration_d[phen$cis129 == 'M' & phen$gt129 == 'MM']))
median(phen$duration_d[phen$cis129 == 'M' & phen$gt129 == 'MV'],na.rm=TRUE) 
sum(!is.na(phen$duration_d[phen$cis129 == 'M' & phen$gt129 == 'MV']))
median(phen$duration_d[phen$cis129 == 'M' & is.na(phen$gt129)],na.rm=TRUE) 
sum(!is.na(phen$duration_d[phen$cis129 == 'M' & is.na(phen$gt129)]))

# Michael Geschwind has noted that in the above calculations, the median
# for cis M trans M is 136.98, while the median for all cis M is 137.1322
# and has questioned how the two medians can be so similar considering
# that the latter figure includes cis M trans V who have much longer
# duration. The answer is that the latter figure also includes cis M 
# trans unknown individuals (of whom there are 10), who have very short
# duration (median 87, see above). The median works out to be almost 
# the same. This is not an error.

ks.test(phen$duration_d[phen$cis129 == 'M' & phen$gt129 == 'MM'],phen$duration_d[phen$cis129 == 'M' & phen$gt129 == 'MV']) 
wilcox.test(phen$duration_d[phen$cis129 == 'M' & phen$gt129 == 'MM'],phen$duration_d[phen$cis129 == 'M' & phen$gt129 == 'MV'])
kruskal.test(phen$duration_d[phen$cis129 == 'M' & !is.na(phen$gt129)],as.factor(phen$gt129[phen$cis129 == 'M' & !is.na(phen$gt129)])) 

# is cis M vs. cis V difference due to haplotype or 129 zygosity?
sum(!is.na(phen$duration_d[phen$cis129=='V'])) 
sum(!is.na(phen$duration_d[phen$cis129=='V' & phen$gt129=='VV'])) 
sum(!is.na(phen$duration_d[phen$cis129=='V' & phen$gt129=='MV'])) 
# Of the cis 129V individuals with a known duration, 2 have an unknown 129 genotype,
# 2 are trans 129M and 2 are trans 129V. There is bias here, as cis 129 could
# in some cases only be phased if there was at least one 129 homozygote in the pedigree
# So that enriches for VV individuals among people with a known cis 129V haplotype
# On the other hand, because V is the minor allele, absent any bias we would expect
# a higher proportion of cis 129V individuals to be 129 heterozygous, which might
# predispose to longer disease duration. In sum it is difficult to conclude
# whether the longer disease duration among cis 129V individuals is due to 
# haplotpe or genotype.

# ANOVA: does duration differ by study?
m = lm(duration_d ~ study, data = phen)
summary(m)
sum(!is.na(phen$duration_d) & !is.na(phen$study))

# because non-normal, use Kruskal-Wallis test- does duration differ by study?
kruskal.test(phen$duration_d, as.factor(phen$study))

# does duration differ by mode of ascertainment?
kruskal.test(phen$duration_d, phen$direct)
ks.test(phen$duration_d[phen$direct], phen$duration_d[!phen$direct])
table(phen$direct[!is.na(phen$duration_d)])

# is duration normally distributed?
shapiro.test(phen$duration_d)

# year of birth desc stats
sum(!is.na(phen$yob)) 
sum(!is.na(phen$yob[phen$died_cjd])) 
# directly ascertained
sum(!is.na(phen$yob[phen$died_cjd & phen$infcat == 1])) 
sum(phen$died_cjd & phen$infcat == 1) 
sum(!is.na(phen$yob[phen$died_cjd & phen$infcat == 1]))/sum(phen$died_cjd & phen$infcat == 1) 
# indirectly ascertained
sum(!is.na(phen$yob[phen$died_cjd & phen$infcat > 1])) 
sum(phen$died_cjd & phen$infcat > 1) 
sum(!is.na(phen$yob[phen$died_cjd & phen$infcat > 1]))/sum(phen$died_cjd & phen$infcat > 1) 

# how many asymptomatic individuals are included?
sum(phen$gene_positive==1 & !phen$died_cjd) # 54
# double checking numbers:
# the year of death desc stats are based on 73 + 93 = 166
# the asymptomatics are 54 
# total = 220

# typical age of parent at birth of child in all individuals
sql_query = "
select   child.yob - parent.yob yobdiff
from     phen parent, phen child
where    parent.iid = child.mother or parent.iid = child.father
-- and      parent.died_cjd and child.died_cjd
;
"
yobdiff = sqldf(sql_query)
mean(yobdiff$yobdiff,na.rm=TRUE)
sd(yobdiff$yobdiff,na.rm=TRUE)

# what does the dist look like?
hist(yobdiff$yobdiff,breaks=100,col='black')

# check the outliers manually to look for any suspect data
sql_query = "
select   parent.iid, parent.sex, parent.yob, child.iid, child.sex, child.yob, 
         child.yob - parent.yob yobdiff
from     phen parent, phen child
where    (parent.iid = child.mother or parent.iid = child.father)
and      ((child.yob - parent.yob < 20) or (child.yob - parent.yob > 40))
and      (child.yob - parent.yob) is not null
;
"
# sqldf(sql_query)
# both of the parents older than 40 at child's birth were males
# all of the parents younger than 20 at child's birth were females

# test for anticipation among
# parent/child pairs where both are dead and neither is listed as dying of CJD
sql_query = "
select   parent.iid, child.iid, 
         parent.age parentage, child.age childage,
         parent.notes pnotes, child.notes cnotes
from     phen parent, phen  child
where    (parent.iid = child.mother or parent.iid = child.father)
and      (not parent.gene_positive = 1)
and      (not child.gene_positive = 1)
and      parent.dead and child.dead
--and      parent.conf_non_cjd and child.conf_non_cjd
--and      child.age > 40 and parent.age > 40
;
"
non_e200k_deaths = sqldf(sql_query)
t.test(non_e200k_deaths$parentage, non_e200k_deaths$childage, paired=TRUE, alternative='two.sided')

# further restrict to parent/child pairs whose cause of death is listed and is confidently
# not CJD
sql_query = "
select   parent.iid, child.iid, 
parent.age parentage, child.age childage,
parent.notes pnotes, child.notes cnotes
from     phen parent, phen  child
where    (parent.iid = child.mother or parent.iid = child.father)
and      (not parent.gene_positive = 1)
and      (not child.gene_positive = 1)
and      parent.dead and child.dead
and      parent.conf_non_cjd and child.conf_non_cjd
--and      child.age > 40 and parent.age > 40
;
"
non_e200k_deaths = sqldf(sql_query)
t.test(non_e200k_deaths$parentage, non_e200k_deaths$childage, paired=TRUE, alternative='two.sided')

# further restrict to individuals dying past the age of 40
sql_query = "
select   parent.iid, child.iid, 
         parent.age parentage, child.age childage,
parent.notes pnotes, child.notes cnotes
from     phen parent, phen  child
where    (parent.iid = child.mother or parent.iid = child.father)
and      (not parent.gene_positive = 1)
and      (not child.gene_positive = 1)
and      parent.dead and child.dead
and      parent.conf_non_cjd and child.conf_non_cjd
and      child.age > 40 and parent.age > 40
;
"
non_e200k_deaths = sqldf(sql_query)
t.test(non_e200k_deaths$parentage, non_e200k_deaths$childage, paired=TRUE, alternative='two.sided')

# get parent-child pairs where both have an age of onset/death
sql_query = "
select   parent.age parentage, parent.infcat parentinfcat,
         child.age childage, child.infcat childinfcat,
         parent.study study -- same for parent and child
from     phen parent, phen child
where    child.affparent = parent.iid
and      child.died_cjd
and      parent.died_cjd
and      parent.age is not null
and      child.age is not null
;
"
po = sqldf(sql_query)

# pairs where both individuals are directly observed
both1 = (po$parentinfcat == 1 & po$childinfcat == 1)
sum(both1)
# pairs where either individual was directly observed
one1  = (po$parentinfcat == 1 | po$childinfcat == 1)
sum(one1)
# all pairs with data available
allp  = rep(TRUE, dim(po)[1])
sum(allp)

# set xy limits for all versions of fig 2d
xylims = range(phen$age[phen$died_cjd],na.rm=TRUE)

# plot only the "both direct" pairs
# png('fig2d.ao.pc.plot1.png',width=500,height=500)
plot(NA,NA,xlim=xylims,ylim=xylims,xaxt='n',yaxt='n',xlab='',ylab='',
     main='Age of onset in parent/child pairs')
axis(side=1,col=pcolor,col.axis=pcolor)
mtext(side=1,col=pcolor,text=expression(bold('Parent age of onset/death')),padj=3)
axis(side=2,col=ccolor,col.axis=ccolor)
mtext(side=2,col=ccolor,text=expression(bold('Child age of onset/death')),padj=-3)
abline(a=0,b=1,col='black')
points(po$parentage[both1],po$childage[both1],pch=19,col='#000000')
legend('bottomleft',c('both direct'),
       pch=c(19),col=c('#000000'))
mtext(side=3, line=-2, text="D", cex=2, adj=0, outer=TRUE)
# dev.off()

# plot the "both direct" and "one direct" pairs
# png('fig2d.ao.pc.plot2.png',width=500,height=500)
plot(NA,NA,xlim=xylims,ylim=xylims,xaxt='n',yaxt='n',xlab='',ylab='',
     main='Age of onset in parent/child pairs')
axis(side=1,col=pcolor,col.axis=pcolor)
mtext(side=1,col=pcolor,text=expression(bold('Parent age of onset/death')),padj=3)
axis(side=2,col=ccolor,col.axis=ccolor)
mtext(side=2,col=ccolor,text=expression(bold('Child age of onset/death')),padj=-3)
abline(a=0,b=1,col='black')
points(po$parentage[one1],po$childage[one1],pch=19,col='#999999')
points(po$parentage[both1],po$childage[both1],pch=19,col='#000000')
legend('bottomleft',c('one direct','both direct'),
       pch=c(19,19),col=c('#999999','#000000'))
mtext(side=3, line=-2, text="D", cex=2, adj=0, outer=TRUE)
# dev.off()

# plot all pairs
# pub-quality Fig 2d
png('fig2d.ao.pc.plot3.png',res=300,pointsize=3,width=500,height=500)
plot(NA,NA,xlim=xylims,ylim=xylims,xaxt='n',yaxt='n',xlab='',ylab='',
     main='Age of onset in parent/child pairs')
axis(side=1,col=pcolor,col.axis=pcolor)
mtext(side=1,col=pcolor,text=expression(bold('Parent age of onset/death')),padj=3)
axis(side=2,col=ccolor,col.axis=ccolor)
mtext(side=2,col=ccolor,text=expression(bold('Child age of onset/death')),padj=-3)
abline(a=0,b=1,col='black')
points(po$parentage[allp],po$childage[allp],pch='O',col='#000000')
points(po$parentage[one1],po$childage[one1],pch=19,col='#999999')
points(po$parentage[both1],po$childage[both1],pch=19,col='#000000')
legend('bottomleft',c('all samples','one direct','both direct'),
       pch=c('O',19,19),col=c('#000000','#999999','#000000'))
legend('bottomleft',c('all samples','one direct','both direct'),
       pch=c(NA,19,19),col=c('#000000','#999999','#000000'),bty='n')
mtext(side=3, line=-2, text="D", cex=2, adj=0, outer=TRUE)
dev.off()


# test for anticipation in various subsets
t.test(po$parentage[both1], po$childage[both1], alternative='two.sided', paired=TRUE)
t.test(po$parentage[one1], po$childage[one1], alternative='two.sided', paired=TRUE)
t.test(po$parentage[allp], po$childage[allp], alternative='two.sided', paired=TRUE)

# how about unpaired?
t.test(po$parentage[allp], po$childage[allp], alternative='two.sided', paired=FALSE)
# surprisingly, unpaired t test gives you the same difference and a slightly smaller p value

# other ways of plotting anticipation
# bar plots of parent and child onset
# png('anticipation.barplot.both1.png',width=600,height=400)
tobj = t.test(po$parentage[both1], po$childage[both1], paired=TRUE)
p = tobj$p.value
n = tobj$parameter+1 # df + 1
diff_yrs = tobj$estimate
subt = paste(formatC(diff_yrs,digits=0,format='f'),' years anticipation at p = ',
             formatC(p,digits=2),' (two-sided paired t test), n = ',n,' pairs',sep='')
barplot(t(as.matrix(po[both1,c("parentage","childage")])),beside=TRUE,
        col=c(pcolor,ccolor),border=NA,xaxt='n',
        xlab='Parent-child pairs, both observed directly',
        ylab='Age of onset or death',
        main='Anticipation in directly observed parent-child pairs',
        sub=subt)
# dev.off()

# png('anticipation.barplot.one1.png',width=600,height=400)
tobj = t.test(po$parentage[one1], po$childage[one1], paired=TRUE)
p = tobj$p.value
n = tobj$parameter+1 # df + 1
diff_yrs = tobj$estimate
subt = paste(formatC(diff_yrs,digits=0,format='f'),' years anticipation at p = ',
             formatC(p,digits=2),' (two-sided paired t test), n = ',n,' pairs',sep='')
barplot(t(as.matrix(po[one1,c("parentage","childage")])),beside=TRUE,
        col=c(pcolor,ccolor),border=NA,xaxt='n',
        xlab='Parent-child pairs, one observed directly',
        ylab='Age of onset or death',
        main='Anticipation in parent-child pairs with one direct observation',
        sub=subt)
# dev.off()

# png('anticipation.barplot.allp.png',width=600,height=400)
tobj = t.test(po$parentage[allp], po$childage[allp], paired=TRUE)
p = tobj$p.value
n = tobj$parameter+1 # df + 1
diff_yrs = tobj$estimate
subt = paste(formatC(diff_yrs,digits=0,format='f'),' years anticipation at p = ',
             formatC(p,digits=2),' (two-sided paired t test), n = ',n,' pairs',sep='')
barplot(t(as.matrix(po[allp,c("parentage","childage")])),beside=TRUE,
        col=c(pcolor,ccolor),border=NA,xaxt='n',
        xlab='Parent-child pairs, all',
        ylab='Age of onset or death',
        main='Anticipation in all parent-child pairs',
        sub=subt)
# dev.off()

# anticipation in subsets of data
t.test(po$parentage[allp & po$study=='ancjdr'], po$childage[allp & po$study=='ancjdr'], paired=TRUE)
t.test(po$parentage[allp & po$study=='goett'], po$childage[allp & po$study=='goett'], paired=TRUE)
t.test(po$parentage[allp & po$study=='mrc'], po$childage[allp & po$study=='mrc'], paired=TRUE)
t.test(po$parentage[allp & po$study=='ucsf'], po$childage[allp & po$study=='ucsf'], paired=TRUE)
# and in all:
t.test(po$parentage[allp], po$childage[allp], paired=TRUE)

# yob-ao correlation in subsets of data
summary(lm(age ~ yob, data = subset(phen, study=='ancjdr' & died_cjd)))
summary(lm(age ~ yob, data = subset(phen, study=='goett' & died_cjd)))
summary(lm(age ~ yob, data = subset(phen, study=='mrc' & died_cjd)))
summary(lm(age ~ yob, data = subset(phen, study=='ucsf' & died_cjd)))

# anticipation when stratifying by 129 genotype. get pairs
# where parent and child both have a known 129 genotype and 
# those genotypes are the same.

# among MM:
sql_query = "
select   parent.age parentage, parent.infcat parentinfcat,
child.age childage, child.infcat childinfcat,
parent.study study -- same for parent and child
from     phen parent, phen child
where    child.affparent = parent.iid
and      child.died_cjd
and      parent.died_cjd
and      parent.age is not null
and      child.age is not null
and      parent.gt129 = 'MM'
and      child.gt129 = 'MM'
;
"
po129mm = sqldf(sql_query)
t.test(po129mm$parentage, po129mm$childage, paired=TRUE) 

# among MV:
sql_query = "
select   parent.age parentage, parent.infcat parentinfcat,
child.age childage, child.infcat childinfcat,
parent.study study -- same for parent and child
from     phen parent, phen child
where    child.affparent = parent.iid
and      child.died_cjd
and      parent.died_cjd
and      parent.age is not null
and      child.age is not null
and      parent.gt129 = 'MV' and parent.cis129 ='M'
and      child.gt129 = 'MV' and child.cis129 ='M'
;
"
po129mv = sqldf(sql_query)
t.test(po129mv$parentage, po129mv$childage, paired=TRUE) # not enough observations

# test for heritability based on parent/offspring regression
# png('cx.paroff.regression.png',width=600,height=400)
m = lm(childage ~ parentage, data=po)
p = summary(m)$coefficients[2,4]
sum(!is.na(po$childage) & !is.na(po$parentage))

subt = paste('p = ',formatC(p,digits=2),sep='')
plot(po$parentage, po$childage, xaxt='n', yaxt='n', pch=19,
     main = 'No correlation between parent and child age of onset/death',
     xlab = '', ylab = '', sub=subt)
axis(side=1,col=pcolor,col.axis=pcolor)
mtext(side=1,col=pcolor,text=expression(bold('Parent')),padj=3)
axis(side=2,col=ccolor,col.axis=ccolor)
mtext(side=2,col=ccolor,text=expression(bold('Child')),padj=-3)
# dev.off()

# test for heritability among different subsets of data
cor.test(po$parentage, po$childage)
cor.test(po$parentage[one1], po$childage[one1])
cor.test(po$parentage[both1], po$childage[both1])

# get all parent-child pairs
sql_query = "
select   parent.study study, 
         parent.gene_positive p_gene_positive,
         parent.died_cjd p_died_cjd,
         parent.age p_age,
         parent.yob p_yob,
         child.gene_positive c_gene_positive,
         child.died_cjd c_died_cjd,
         child.age c_age,
         child.yob c_yob
from     phen parent, phen child
where    child.affparent = parent.iid
;
"
pc = sqldf(sql_query)

# number of gene-positive parent-child pairs
sum(pc$p_gene_positive == 1 & pc$c_gene_positive == 1)

# descriptive stats for parent-child pairs
sql_query = "
select   study,
         count(*),
         sum(p_gene_positive = 1 and c_gene_positive = 1) both_gp
from     pc
group by study
order by study
"
pc_descstats = sqldf(sql_query)
pc_descstats

descstats$gp_pairs = pc_descstats$both_gp

# Create Table 3
summary = c('all',sum(descstats$n_gene_positive),sum(descstats$gp_pairs),
            sum(descstats$pct_indirect*descstats$n_gene_positive)/sum(descstats$n_gene_positive),
            sum(descstats$pct_asymp*descstats$n_gene_positive)/sum(descstats$n_gene_positive))
descstats_for_tbl3 = rbind(descstats[,c('study','n_gene_positive','gp_pairs','pct_indirect','pct_asymp')],
                          summary)
descstats_for_tbl3

write.table(descstats_for_tbl3,'descstats.txt',sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE)

# anticipation survival analysis
# parent-child pairs where both are gene positive
sql_query = "
select   parent.age parentage, parent.censored parentcensored, parent.infcat parentinfcat,
child.age childage, child.censored childcensored, child.infcat childinfcat,
parent.study study, -- same for parent and child
parent.iid piid, child.iid ciid
from     phen parent, phen child
where    child.affparent = parent.iid
and      child.gene_positive = 1
and      parent.gene_positive = 1
and      parent.age is not null
and      child.age is not null
;
"
po_all = sqldf(sql_query)
dim(po_all)
table(po_all$study)
sum(is.na(phen$age)) 
sum(!is.na(phen$age) & phen$affparent %in% (phen$iid)) 
sum(!is.na(phen$age) & phen$affparent %in% (phen$iid) & 
      !is.na(phen$age[match(phen$affparent,phen$iid)]))
table(phen[,c("gene_positive","study")])

# create a variable for difference in parent and child age per R1's suggestion
po_all$agediff = po_all$parentage - po_all$childage

# examine number of children of affected parents whose mutation status is not known, in the German data
goett_parents = phen$iid[phen$study=='goett' & phen$died_cjd]
table(phen$gene_positive[!is.na(phen$age) & phen$affparent %in% goett_parents & !phen$died_cjd]) 

# rearrange data to perform overall survival analysis
po_all_parents = po_all[,1:3]
po_all_children = po_all[,4:6]
po_all_parents$grp = 'parent'
po_all_children$grp = 'child'
colnames(po_all_parents) = c('age','censored','infcat','grp')
colnames(po_all_children) = c('age','censored','infcat','grp')
po_all_long = rbind(po_all_parents, po_all_children)

# parent vs. child survival analysis, all affected pairs
# note that in some cases, an individual may be a parent in one pair
# and a child in another pair
# png('cx.po.survival.analysis.png',width=600,height=400)
survdata = data.frame(ages=po_all_long$age, status=as.integer(!po_all_long$censored), group=po_all_long$grp)
msurv = with(survdata, Surv(ages, status==1))
diffobj = survdiff(Surv(ages,status==1)~group, data=survdata)
n = as.integer(sum(diffobj$n))
diffobj
n
p = 1-pchisq(diffobj$chisq,df=1)
main = 'Survival of parents and children with the E200K mutation'
subt = paste('n = ',n,' individuals, log-rank test p value = ',formatC(p,digits=2),sep='') 
plot(survfit(Surv(ages,status)~group,data=survdata),col=c(ccolor,pcolor),main=main,sub=subt,lwd=3)
legend('bottomleft',c('children','parents'),col=c(ccolor,pcolor),lwd=3)
sf = survfit(Surv(ages,status)~group,data=survdata)
sf # summary of survival curves
medians = summary(sf)$table[,"median"] # extract medians from summary
as.integer(medians[2] - medians[1]) # take difference of medians of survival curves
# dev.off()

# Ben Neale's suggestion - parent-child survival analysis w/o double-counting
# randomly select one pair from each pedigree for survival analysis
po_pairs = sqldf("
select   parent.age parentage, parent.censored parentcensored, parent.infcat parentinfcat,
child.age childage, child.censored childcensored, child.infcat childinfcat,
parent.study study, parent.fid fid, -- same for parent and child
parent.iid piid, child.iid ciid
from     phen parent, phen child
where    child.affparent = parent.iid
and      child.gene_positive = 1
and      parent.gene_positive = 1
and      parent.age is not null
and      child.age is not null
;")
dim(po_pairs)

# first try the randomization once
set.seed(1)
po_pairs$rand = runif(0,1,n=dim(po_pairs)[1]) # assign a random number to each pair
min_rands = sqldf("
                  select   p.fid, min(p.rand) minrand
                  from     po_pairs p
                  group by 1
                  order by 1
                  ;")
# only use the pair in each family with the lowest random number
po_use = sqldf("
               select   po.*
               from     po_pairs po, min_rands mr
               where    po.fid = mr.fid
               and      po.rand = mr.minrand
               ;")
dim(po_use)

# rearrange data to perform survival analysis
po_mr_parents = po_use[,1:3]
po_mr_children = po_use[,4:6]
po_mr_parents$grp = 'parent'
po_mr_children$grp = 'child'
colnames(po_mr_parents) = c('age','censored','infcat','grp')
colnames(po_mr_children) = c('age','censored','infcat','grp')
po_mr_long = rbind(po_mr_parents, po_mr_children)
survdata = data.frame(ages=po_mr_long$age, status=as.integer(!po_mr_long$censored), group=po_mr_long$grp)
msurv = with(survdata, Surv(ages, status==1))
diffobj = survdiff(Surv(ages,status==1)~group, data=survdata)
n = as.integer(sum(diffobj$n))
diffobj
n
p = 1-pchisq(diffobj$chisq,df=1)
main = 'Survival of parents and children with the E200K mutation'
subt = paste('n = ',n,' individuals, log-rank test p value = ',formatC(p,digits=2),sep='') 
plot(survfit(Surv(ages,status)~group,data=survdata),col=c(ccolor,pcolor),main=main,sub=subt,lwd=3)
legend('bottomleft',c('children','parents'),col=c(ccolor,pcolor),lwd=3)
sf = survfit(Surv(ages,status)~group,data=survdata)
sf # summary of survival curves
medians = summary(sf)$table[,"median"] # extract medians from summary
as.integer(medians[2] - medians[1]) # take difference of medians of survival curves

# that worked, now repeat 1000 times to report aggregate figures
n_iter = 1000
surv_results = data.frame(seed=1:n_iter,p=numeric(n_iter),meddiff=integer(n_iter))
for (seed in 1:n_iter) {
  set.seed(seed)
  po_pairs$rand = runif(0,1,n=dim(po_pairs)[1])
  min_rands = sqldf("
                    select   p.fid, min(p.rand) minrand
                    from     po_pairs p
                    group by 1
                    order by 1
                    ;")
  po_use = sqldf("
                 select   po.*
                 from     po_pairs po, min_rands mr
                 where    po.fid = mr.fid
                 and      po.rand = mr.minrand
                 ;")
  # rearrange data to perform survival analysis
  po_mr_parents = po_use[,1:3]
  po_mr_children = po_use[,4:6]
  po_mr_parents$grp = 'parent'
  po_mr_children$grp = 'child'
  colnames(po_mr_parents) = c('age','censored','infcat','grp')
  colnames(po_mr_children) = c('age','censored','infcat','grp')
  po_mr_long = rbind(po_mr_parents, po_mr_children)
  survdata = data.frame(ages=po_mr_long$age, status=as.integer(!po_mr_long$censored), group=po_mr_long$grp)
  msurv = with(survdata, Surv(ages, status==1))
  diffobj = survdiff(Surv(ages,status==1)~group, data=survdata)
  n = as.integer(sum(diffobj$n))
  p = 1-pchisq(diffobj$chisq,df=1)
  sf = survfit(Surv(ages,status)~group,data=survdata)
  medians = summary(sf)$table[,"median"] # extract medians from summary
  surv_results$p[seed] = p
  surv_results$meddiff[seed] = as.integer(medians[2] - medians[1]) # take difference of medians of survival curves
}
mean(surv_results$meddiff)
sum(surv_results$p < .05)

# exact same code as above but for UCSF and MRC only
po_pairs = sqldf("
                    select   parent.age parentage, parent.censored parentcensored, parent.infcat parentinfcat,
                    child.age childage, child.censored childcensored, child.infcat childinfcat,
                    parent.study study, parent.fid fid, -- same for parent and child
                    parent.iid piid, child.iid ciid
                    from     phen parent, phen child
                    where    child.affparent = parent.iid
                    and      child.gene_positive = 1
                    and      parent.gene_positive = 1
                    and      parent.age is not null
                    and      child.age is not null
                    and      parent.study in ('ucsf','mrc')
                    ;
                    ")
dim(po_pairs)
set.seed(1)
po_pairs$rand = runif(0,1,n=dim(po_pairs)[1])
min_rands = sqldf("
                  select   p.fid, min(p.rand) minrand
                  from     po_pairs p
                  group by 1
                  order by 1
                  ;")
po_use = sqldf("
               select   po.*
               from     po_pairs po, min_rands mr
               where    po.fid = mr.fid
               and      po.rand = mr.minrand
               ;")
dim(po_use)

# rearrange data to perform survival analysis
po_mr_parents = po_use[,1:3]
po_mr_children = po_use[,4:6]
po_mr_parents$grp = 'parent'
po_mr_children$grp = 'child'
colnames(po_mr_parents) = c('age','censored','infcat','grp')
colnames(po_mr_children) = c('age','censored','infcat','grp')
po_mr_long = rbind(po_mr_parents, po_mr_children)
survdata = data.frame(ages=po_mr_long$age, status=as.integer(!po_mr_long$censored), group=po_mr_long$grp)
msurv = with(survdata, Surv(ages, status==1))
diffobj = survdiff(Surv(ages,status==1)~group, data=survdata)
n = as.integer(sum(diffobj$n))
diffobj
n
p = 1-pchisq(diffobj$chisq,df=1)
main = 'Survival of parents and children with the E200K mutation'
subt = paste('n = ',n,' individuals, log-rank test p value = ',formatC(p,digits=2),sep='') 
plot(survfit(Surv(ages,status)~group,data=survdata),col=c(ccolor,pcolor),main=main,sub=subt,lwd=3)
legend('bottomleft',c('children','parents'),col=c(ccolor,pcolor),lwd=3)
sf = survfit(Surv(ages,status)~group,data=survdata)
sf # summary of survival curves
medians = summary(sf)$table[,"median"] # extract medians from summary
as.integer(medians[2] - medians[1]) # take difference of medians of survival curves

n_iter = 1000
surv_results = data.frame(seed=1:n_iter,p=numeric(n_iter),meddiff=integer(n_iter))
for (seed in 1:n_iter) {
  set.seed(seed)
  po_pairs$rand = runif(0,1,n=dim(po_pairs)[1])
  min_rands = sqldf("
                    select   p.fid, min(p.rand) minrand
                    from     po_pairs p
                    group by 1
                    order by 1
                    ;")
  po_use = sqldf("
                 select   po.*
                 from     po_pairs po, min_rands mr
                 where    po.fid = mr.fid
                 and      po.rand = mr.minrand
                 ;")
  # rearrange data to perform survival analysis
  po_mr_parents = po_use[,1:3]
  po_mr_children = po_use[,4:6]
  po_mr_parents$grp = 'parent'
  po_mr_children$grp = 'child'
  colnames(po_mr_parents) = c('age','censored','infcat','grp')
  colnames(po_mr_children) = c('age','censored','infcat','grp')
  po_mr_long = rbind(po_mr_parents, po_mr_children)
  survdata = data.frame(ages=po_mr_long$age, status=as.integer(!po_mr_long$censored), group=po_mr_long$grp)
  msurv = with(survdata, Surv(ages, status==1))
  diffobj = survdiff(Surv(ages,status==1)~group, data=survdata)
  n = as.integer(sum(diffobj$n))
  p = 1-pchisq(diffobj$chisq,df=1)
  sf = survfit(Surv(ages,status)~group,data=survdata)
  medians = summary(sf)$table[,"median"] # extract medians from summary
  surv_results$p[seed] = p
  surv_results$meddiff[seed] = as.integer(medians[2] - medians[1]) # take difference of medians of survival curves
}
mean(surv_results$meddiff)
sum(surv_results$p < .05)


#  sibling pair regression for heritability?
sql_query = "
select   sib1.age sib1age, sib2.age sib2age
from     phen sib1, phen sib2
where    sib1.affparent = sib2.affparent
and      sib1.iid > sib2.iid
and      sib1.died_cjd and sib2.died_cjd
;
"
sibs = sqldf(sql_query)

cor.test(sibs$sib1age, sibs$sib2age)
m = lm(sib1age ~ sib2age, data=sibs)
summary(m)
sum(!is.na(sibs$sib1age) & !is.na(sibs$sib2age))

# what fraction of at-risk individuals have undergone predictive testing?
gene_pos_iids = phen$iid[phen$gene_positive == 1]
# at risk are those with a gene pos parent and self not yet dead of CJD
phen$at_risk = (phen$mother %in% gene_pos_iids | phen$father %in% gene_pos_iids) & !phen$died_cjd
# also at risk are those with an at risk parent who did not test negative
# so iterate for a few generations here (manually checked to see when no additional
# individuals were being added)
at_risk_iids = phen$iid[phen$at_risk & (phen$gene_positive != 0)]
phen$at_risk = phen$at_risk | (phen$mother %in% at_risk_iids | phen$father %in% at_risk_iids)
at_risk_iids = phen$iid[phen$at_risk & (phen$gene_positive != 0)]
phen$at_risk = phen$at_risk | (phen$mother %in% at_risk_iids | phen$father %in% at_risk_iids)
at_risk_iids = phen$iid[phen$at_risk & (phen$gene_positive != 0)]
phen$at_risk = phen$at_risk | (phen$mother %in% at_risk_iids | phen$father %in% at_risk_iids)

# summary of at-risk people by mutation status (-1 = unknown, 0 = negative, 1 = positive)
table(phen$gene_positive[phen$at_risk])

# add in the manual counts of at-risk individuals who were omitted from the phen
# spreadsheet, from the UCSF pedigrees
ucsf_at_risk = read.table('ucsf-at-risk-counts.txt',header=TRUE)
# N.B. in the manual UCSF counting, as above, at risk is rigorously defined as a direct 
# descendent of someone who tests positive, without any negative tests in intervening generations.
# Sibs of single affecteds w/o other family history are NOT considered at risk.
ucsf_at_risk_n = sum(ucsf_at_risk$n_atrisk) 
ucsf_at_risk_n

tested = sum(phen$gene_positive[phen$at_risk] %in% c(0,1))
untested = sum(phen$gene_positive[phen$at_risk] == -1) + ucsf_at_risk_n

proportion_tested = tested / (tested+untested)
proportion_tested

# % of people ascertained directly who have no family history
sum(is.na(phen$affparent[phen$infcat==1 & phen$died_cjd]))/sum(phen$infcat==1 & phen$died_cjd)

# note that it goes down to 42% if you also include asymptomatics seen directly,  i.e. 
sum(is.na(phen$affparent[phen$infcat==1]))/sum(phen$infcat==1) # = 42%

