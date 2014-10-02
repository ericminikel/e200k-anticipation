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

year_of_birth = c(parent_yob[pair_ascertainable],child_yob[pair_ascertainable])
age_of_onset = c(parent_onset_age[pair_ascertainable],child_onset_age[pair_ascertainable])
cor.test(year_of_birth,age_of_onset)

png('~/d/j/cureffi/media/2014/10/yob-ao-correlation.png',width=750,height=500,res=100)
plot(year_of_birth,age_of_onset,pch=20,axes=FALSE,xlab='Year of birth',
     ylab='Age of onset',main='Year of birth and age of onset are correlated\nin simulated ascertained pairs')
abline(a=1989,b=-1,col='gray',lwd=2)
abline(a=2013,b=-1,col='gray',lwd=2)
axis(side=1,at=c(190:198)*10,lwd=0,lwd.ticks=1,cex.axis=.8)
axis(side=2,at=c(2:10)*10,lwd=0,lwd.ticks=1,las=1,cex.axis=.8)
points(c(1951,1951),c(55,80),col='red',type='l',lwd=2)
text(x=1951,y=80,pos=4,label='55 years old',col='red')
points(c(1980,1980),c(25,60),col='red',type='l',lwd=2)
text(x=1980,y=60,pos=2,label='25 years old',col='red')
dev.off()

mean(age_of_onset[round(year_of_birth)==1951])
mean(age_of_onset[round(year_of_birth) > 1979])

png('~/d/j/cureffi/media/2014/10/parent-child-correlation.png',width=750,height=500,res=100)
m = lm(child_onset_age[pair_ascertainable] ~ parent_onset_age[pair_ascertainable])
plot(parent_onset_age[pair_ascertainable],child_onset_age[pair_ascertainable],pch=20,axes=FALSE,
     xlab='',ylab='',main='Ascertainment bias causes parent and child age of onset to be correlated',
     ylim=c(0,100),xlim=c(0,100))
axis(side=1,at=(0:10)*10,col=pcolor,col.axis=pcolor)
mtext(side=1,col=pcolor,text=expression(bold('Parent age of onset/death')),padj=3)
axis(side=2,at=(0:10)*10,col=ccolor,col.axis=ccolor,las=1)
mtext(side=2,col=ccolor,text=expression(bold('Child age of onset/death')),padj=-3)
abline(a=0,b=1,col='black')
abline(m,col='red')
summary(m)
dev.off()
