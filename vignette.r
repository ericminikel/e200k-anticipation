setwd('~/d/j/cureffi/media/2014/10/')

# colors for plotting parents and children
pcolor='#D18532' 
ccolor='#091090' 

# example of two parent-child pairs
ex1 = data.frame(
  family = c(1,1,2,2),
  plot_height = c(4,3,2,1),
  year_of_birth = c(1940,1970,1940,1970),
  lifespan = c(60,40,40,60),
  color = c(pcolor,ccolor,pcolor,ccolor)
  )

png('example-of-two-families.png',width=750,height=500,res=100)
par(mar=c(2,6,1,1))
plot(NA,NA,xlim=c(1930,2030),ylim=c(.5,4.5),axes=FALSE,ylab='',xlab='')
# draw a rectangle for the period of ascertainment, and label it
rect(xleft=1989,ybottom=.5,xright=2013,ytop=4.5,col='#C6C3B5',border=NA)
text(x=(1989+2013)/2,y=.5,pos=3,label='1989-2013',col='#565345')
# draw lines from parent to child like on a pedigree chart
points(c(1940,1940),c(4,3.5),type='l',lwd=1)
points(c(1940,1970),c(3.5,3.5),type='l',lwd=1)
points(c(1970,1970),c(3.5,3),type='l',lwd=1)
points(c(1940,1940),c(2,1.5),type='l',lwd=1)
points(c(1940,1970),c(1.5,1.5),type='l',lwd=1)
points(c(1970,1970),c(1.5,1),type='l',lwd=1)
# draw each individual's birth and lifespan
for (i in 1:4) {
  points(ex1$year_of_birth[i],ex1$plot_height[i],pch=19,cex=5,col=ex1$color[i])
  points(c(ex1$year_of_birth[i],ex1$year_of_birth[i]+ex1$lifespan[i]),rep(ex1$plot_height[i],2),type='l',lwd=10,col=ex1$color[i])
  text(ex1$year_of_birth[i]+ex1$lifespan[i]/2,ex1$plot_height[i],pos=3,col=ex1$color[i],label=paste(ex1$lifespan[i],'years'))
}
# add axes
axis(side=2,at=c(.7,2.3),labels=c('',''),lwd.ticks=0)
axis(side=2,at=c(2.7,4.3),labels=c('',''),lwd.ticks=0)
axis(side=2,at=c(3.5,1.5),labels=c('Family 1','Family 2'),las=1,lwd=0,lwd.ticks=0)
axis(side=1,at=c(1940,1970,2000,2030),lwd=0,lwd.ticks=1)
dev.off()


# simulation
set.seed(1)
parent_yob = runif(n=100000,min=1700,max=2000) # generate parents born from 1700 through 2000
child_yob = parent_yob + rnorm(n=100000,m=28,s=6) # child born when parent is 28+-6 years old
parent_onset_age = round(rnorm(n=100000, m=64, s=10)) # parent's age of onset is 64+-10
child_onset_age = round(rnorm(n=100000, m=64, s=10)) # child's age of onset is 64+-10
parent_onset_year = parent_yob + parent_onset_age # figure out what year parent has onset
child_onset_year = child_yob + child_onset_age # figure out what year child has onset

# if we can magically ascertain _everyone_ regardless of year of onset, is there anticipation?
t.test(parent_onset_age,child_onset_age,paired=TRUE,alternative='two.sided')
# no significant difference in age of onset

# what if we can only ascertain pairs where both get sick in 1989-2013?
parent_ascertainable = parent_onset_year >= 1989 & parent_onset_year <= 2013
child_ascertainable = child_onset_year >= 1989 & child_onset_year <= 2013
pair_ascertainable = parent_ascertainable & child_ascertainable
t.test(parent_onset_age[pair_ascertainable],child_onset_age[pair_ascertainable],paired=TRUE,alternative='two.sided')
# highly significant 17-year difference in age of onset

png('anticipation-02-anticipation-visualized.png',width=750,height=500,res=100)
par(mar=c(4,4,1,1))
plot(parent_onset_age[pair_ascertainable],child_onset_age[pair_ascertainable],pch=19,
     xlim=c(0,100),ylim=c(0,100),axes=FALSE,xlab='',ylab='')
axis(side=1,at=(0:10)*10,col=pcolor,col.axis=pcolor)
mtext(side=1,col=pcolor,text=expression(bold('Parent age of onset/death')),padj=3)
axis(side=2,at=(0:10)*10,col=ccolor,col.axis=ccolor,las=1)
mtext(side=2,col=ccolor,text=expression(bold('Child age of onset/death')),padj=-3)
abline(a=0,b=1,col='black')
dev.off()

# FamiLinx analysis
# SQL query that generated the flat table
# select   r.Child_id, r.Parent_id, (yc.Dyear - yc.Byear) child_ad, (yp.Dyear - yp.Byear) parent_ad
# from     relationship r, years yc, years yp
# where    r.Child_id = yc.Id
# and      r.Parent_id = yp.Id
# and      yc.Byear > -1 and yc.Dyear > -1
# and      yp.Byear > -1 and yp.Dyear > -1
# into outfile '~/d/sci/039famil/familinx/adpairs.txt'
# ;

familinx_pairs = read.table('~/d/sci/039famil/familinx/adpairs.txt',header=FALSE)
colnames(familinx_pairs) = c('child_id','parent_id','child_age_death','parent_age_death')
is_valid = familinx_pairs$child_age_death >= 0 & familinx_pairs$child_age_death <= 120 &
  familinx_pairs$parent_age_death >= 0 & familinx_pairs$parent_age_death <= 120
t.test(familinx_pairs$parent_age_death[is_valid],familinx_pairs$child_age_death[is_valid],paired=TRUE,alternative='two.sided')
# highly significant 13.6 year difference in age of death
over_40 = is_valid & familinx_pairs$child_age_death >= 40 & familinx_pairs$parent_age_death >= 40
t.test(familinx_pairs$parent_age_death[over_40],familinx_pairs$child_age_death[over_40],paired=TRUE,alternative='two.sided')
# 0.85 years anticipation at p < 2.2e-16

# plot(familinx_pairs$parent_age_death[is_valid],familinx_pairs$child_age_death[is_valid],pch='.')
# plot(familinx_pairs$parent_age_death[over_40],familinx_pairs$child_age_death[over_40],pch='.')
# over_50 = is_valid & familinx_pairs$child_age_death >= 50 & familinx_pairs$parent_age_death >= 50
# t.test(familinx_pairs$parent_age_death[over_50],familinx_pairs$child_age_death[over_50],paired=TRUE,alternative='two.sided')


set.seed(1)
parent_yob = runif(n=100000,min=1700,max=2000) # generate parents born from 1700 through 2000
child_yob = parent_yob + rnorm(n=100000,m=28,s=6) # child born when parent is 28+-6 years old
parent_hypo_onset_age = round(rnorm(n=100000, m=64, s=10)) # parent's age of onset is 64+-10
child_hypo_onset_age = round(rnorm(n=100000, m=64, s=10)) # child's age of onset is 64+-10
parent_hypo_onset_year = parent_yob + parent_hypo_onset_age # figure out what year parent has onset
child_hypo_onset_year = child_yob + child_hypo_onset_age # figure out what year child has onset
# read in the actuarial life tables
life = read.table('~/d/sci/src/e200k-anticipation/ssa.actuarial.table.2009.txt',sep='\t',header=TRUE)
life$age = as.integer(life$age)
life$msur = as.numeric(life$msur)
life$fsur = as.numeric(life$fsur)
# calculate percent surviving.  this is 1-CDF.
life$mpctsur = life$msur/100000 
life$fpctsur = life$fsur/100000
life$allpctsur = (life$msur + life$fsur) / 200000
# calculate probability, from age 0:119 of dying at each possible age.  this is the PDF.
life$pdf = life$allpctsur - c(life$allpctsur[-1],0)
# simulate when the people would have died of other causes
parent_intercurrent_age = sample(0:119, size=100000, prob=life$pdf, replace=TRUE)
child_intercurrent_age = sample(0:119, size=100000, prob=life$pdf, replace=TRUE)

# if we can magically ascertain _everyone_ regardless of year of onset, is there anticipation?
t.test(parent_hypo_onset_age,child_hypo_onset_age) # no.

# what if we can only ascertain pairs where both get sick in 1989-2013?
parent_ascertainable = parent_hypo_onset_age < parent_intercurrent_age & parent_hypo_onset_year > 1989 & parent_hypo_onset_year < 2013
child_ascertainable = child_hypo_onset_age < child_intercurrent_age & child_hypo_onset_year > 1989 & child_hypo_onset_year < 2013
pair_ascertainable = parent_ascertainable & child_ascertainable
t.test()





t.test(parent_hypo_onset_age,child_hypo_onset_age,paired=TRUE,alternative='two.sided')
t.test(parent_hypo_onset_age[pair_ascertainable],child_hypo_onset_age[pair_ascertainable],paired=TRUE,alternative='two.sided')



t.test(parent_hypo_onset_age,child_hypo_onset_age,paired=TRUE,alternative='two.sided')
t.test(parent_hypo_onset_age[pair_ascertainable],child_hypo_onset_age[pair_ascertainable],paired=TRUE,alternative='two.sided')

cor.test(c(parent_yob[pair_ascertainable],child_yob[pair_ascertainable]),c(parent_hypo_onset_age[pair_ascertainable],child_hypo_onset_age[pair_ascertainable]))