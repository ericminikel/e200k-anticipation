#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/r_3.0.2/bin/Rscript

# Eric Vallabh Minikel
# Demonstrates that "windowing" of year of onset is sufficient to create robust 
# and highly significant false positive signals of anticipation and heritability

# This simulation will generate n parent/child pairs and then ascertain them
# according to user-specified criteria. The criteria I changed most often
# can be changed via command-line parameters; others are still hard-coded.
# Once the pairs are ascertained, a variety of tests are performed and 
# plots generated to show how the ascertainment method creates bias.
# In "fullmode" this script will print the output of many tests to stdout
# which is useful for getting a sense of how the simulation works.
# If "imgdir" is set, it will save plots to the specified directory,
# otherwise it creates no plots.
# By default it runs in sparse mode where all it does is output the 
# figures for anticipation (years), heritability (%) and year of birth
# age of onset correlation (slope units). This mode is useful for 
# running large batches with different parameters.

suppressPackageStartupMessages(require(survival))
options(stringsAsFactors=FALSE)

suppressPackageStartupMessages(require(optparse)) # http://cran.r-project.org/web/packages/optparse/optparse.pdf

option_list = list(
  make_option(c("-n", "--npairs"), action="store", default=40000,
              type='integer', help="Number of parent/child pairs"),
  make_option(c("--pyearmin"), action="store", default=1800,
              type='integer', help="Minimum parent birth year"),
  make_option(c("--pyearmax"), action="store", default=1980,
              type='integer', help="Maximum parent birth year"),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print verbose output [default %default]"),
  make_option(c("-q", "--quiet"), action="store_false", dest="verbose",
              help="Do not print verbose output (this is the default)"),
  make_option(c("-i", "--imgdir"), action="store", default="",
              type='character', help="Directory to save images to (default is no images)"),
  make_option(c("-f", "--fullmode"), action="store_true", default=FALSE,
              help="Run in full output mode"),
  make_option(c("-a", "--amode"), action="store", default=1,
              type="integer", help="Ascertainment mode [default %default]"),
  make_option(c("-r", "--rate"), action="store", default=.05,
              type="numeric", help="Rate of decline in ascertainment [default %default]"),
  make_option(c("--amin"), action="store", default=1989,
              type="integer", help="Minimum ascertainable year [default %default]"),
  make_option(c("--amax"), action="store", default=2013,
              type="integer", help="Maximum ascertainable year [default %default]"),
  make_option(c("-s","--seed"), action="store", default=2222,
              type="integer", help="Seed for randomization [default %default]"),
  make_option(c("-l","--limit_n"), action="store", default=NA,
              type="integer", help="Limit to this many ascertained pairs"),
  make_option(c("--heritability_is_real"), action="store_true", default=FALSE,
              help="Simulate heritability of age of onset being real")
)
opt = parse_args(OptionParser(option_list=option_list))

# uncomment to run in interactive mode instead of from command line
# opt = list()
# opt[["n"]] = 100000
# opt[["verbose"]] = FALSE
# opt[["pyearmin"]] = 1700
# opt[["pyearmax"]] = 2000
# opt[["imgdir"]] = ""
# opt[["fullmode"]] = TRUE
# opt[["amode"]] = 1
# opt[["rate"]] = .05
# opt[["amin"]] = 1989
# opt[["amax"]] = 2013
# opt[["seed"]] = 1
# opt[["limit_n"]] = NA
# opt[["heritability_is_real"]] = FALSE


# set random seed
seed = opt$seed
set.seed(seed) # 2222

# bother to create images??
imgs = opt$imgdir != ""
imgdir = opt$imgdir



######
# graphical parameters
pcolor='#D18532' #'orange'
ccolor='#091090' #'violet'

######
# load & process U.S. actuarial life tables for life expectancy censoring
# http://www.ssa.gov/oact/STATS/table4c6.html
# life = read.table('~/d/sci/src/e200k-anticipation/ssa.actuarial.table.2009.txt',sep='\t',header=TRUE)
life = read.table('ssa.actuarial.table.2009.txt',sep='\t',header=TRUE)
life$age = as.integer(life$age)
life$msur = as.numeric(life$msur)
life$fsur = as.numeric(life$fsur)
# calculate percent surviving.  this is 1-CDF.
life$mpctsur = life$msur/100000 
life$fpctsur = life$fsur/100000
life$allpctsur = (life$msur + life$fsur) / 200000
# calculate probability, from age 0:119 of dying at each possible age.  this is the PDF.
life$pdf = life$allpctsur - c(life$allpctsur[-1],0)
# note: here is how to sample one person's age at death from the life expectancy distribution:
# sample(0:119, size=1, prob=life$pdf, replace=TRUE)

#######
# model parameters you can tweak
nparents = opt$n # 40000 # 40000
onset_mean = 64 # orig 65
onset_sd = 10 # orig 15
parent_yob_min = opt$pyearmin # 1800 # originally 1870
parent_yob_max = opt$pyearmax # 2000 # originally 1980
fertility_age_mean = 28 # orig 26
fertility_age_sd = 6 # orig 3
ascertainable_year_min = opt$amin # 1989 # try 1880 - results still significant
ascertainable_year_max = opt$amax # 2013
ascertainable_year_min_mode3 = 1969 # never made this adjustable via command line & not used in paper
ascertainable_dec_rate_mode4 = opt$rate # .01
ascertain_mode = opt$amode # 1
# 1 = parent and child must BOTH be in ascertainable range for pair to be asceratined (Simulations 1-3 in paper)
# 2 = only one must be in range; if so, the other can also be ascertained regardless of how early year of onset (Simulation 10 in paper)
# 3 = like mode 1, but the second individual must still be within a more generous range (e.g. 1969-2013) instead of the strict range (1989-2013)
#     we didn't end up using mode 3 in the paper
# 4 = like mode 1, but the probability of ascertaining the second individual declines by X% for each year before the min date. (Simulation 4-9 in paper)

if (opt$verbose) {
  print(opt,file=stderr())
}

# accepts a list of vectors of identical length and returns one vector with the first non-NA value
coalesce = function(...) {
  # convert input arguments into a list of vectors
  input_list = list(...)
  # check that all input vectors are of same length
  vectorlength = length(input_list[[1]])
  for (j in 1:length(input_list)) {
    if(length(input_list[[j]]) != vectorlength) {
      stop(paste("Not all vectors are of same length. First vector length: ",vectorlength,". Vector #",j,"'s length: ",length(input_list[[j]]),sep=""))
    }
  }
  # create a result vector to fill with first non-NA values
  result = rep(NA,vectorlength)
  # fill with first non-NA value
  for (i in 1:length(result)) {
    for (j in 1:length(input_list)) {
      if(!is.na(input_list[[j]][i])) {
        result[i] = input_list[[j]][i]
        break
      }
    }
  }
  return(result)
}


#######
# simulation
parent_yob = round(runif(n=nparents, min=parent_yob_min, max=parent_yob_max))
child_yob = round(rnorm(n=nparents, m=fertility_age_mean, s=fertility_age_sd)) + parent_yob
parent_hypo_onset_age = round(rnorm(n=nparents, m=onset_mean, s=onset_sd)) # parent's hypothetical onset age, if they lived long enough to have onset
parent_hypo_onset_year = parent_yob + parent_hypo_onset_age
child_hypo_onset_age = round(rnorm(n=nparents, m=onset_mean, s=onset_sd))
child_hypo_onset_year = child_yob + child_hypo_onset_age
parent_intercurrent_age = sample(0:119, size=nparents, prob=life$pdf, replace=TRUE)
child_intercurrent_age = sample(0:119, size=nparents, prob=life$pdf, replace=TRUE)

# test results of having actual heritability 
if (opt$heritability_is_real){
  child_hypo_onset_age = rowMeans(cbind(round(rnorm(n=nparents, m=onset_mean, s=onset_sd)),parent_hypo_onset_age))
  child_hypo_onset_year = child_yob + child_hypo_onset_age
}


# determine "case" status
parent_case = parent_hypo_onset_age < parent_intercurrent_age # is this person a "case" i.e. do they have an age of onset prior to death of intercurrent illness?
child_case = child_hypo_onset_age < child_intercurrent_age # is this person a "case" i.e. do they have an age of onset prior to death of intercurrent illness?

# "hypothetical" scenario: include all individuals who had onset in their lifetimes, i.e. did not die of intercurrent disease, even if their death is after 2013
parent_onset_age = parent_hypo_onset_age
parent_onset_age[!parent_case] = NA
parent_onset_year = parent_hypo_onset_year
parent_onset_year[!parent_case] = NA
child_onset_age = child_hypo_onset_age
child_onset_age[!child_case] = NA
child_onset_year = child_hypo_onset_year
child_onset_year[!child_case] = NA

# evaluate how the "competing risks" due to life expectancy change the estimated age of onset.
mean(c(parent_hypo_onset_age,child_hypo_onset_age))
mean(c(parent_onset_age,child_onset_age),na.rm=TRUE)
# do a survival analysis on all individuals, without ascertainment, but accounting for intercurrent death censoring
all_ages = c(pmin(parent_hypo_onset_age, parent_intercurrent_age), pmin(child_hypo_onset_age,  child_intercurrent_age))
all_status = c(parent_hypo_onset_age < parent_intercurrent_age, child_hypo_onset_age < child_intercurrent_age) # TRUE = event / died of E200K. FALSE = no event / died intercurrently
survdata = data.frame(ages=all_ages,status=all_status)
mfit = survfit(Surv(ages,status==1) ~ 1, data = survdata)
mfit


if (imgs) {
    # Fig 1A without lines
    png(file.path(imgdir,'fig1a.withoutlines.png'),res=300,pointsize=3,width=500,height=500)
    plot(NA,NA,xlim=c(1800,2070),ylim=c(1,50),main='Lifespan of parent/child pairs',xlab='Year',ylab='Parent/child pair id')
    for (i in 1:50) {
        points(c(parent_yob[i],parent_hypo_onset_year[i]),c(i,i),col=pcolor,type='l',lwd=3)
        points(parent_yob[i],i,col=pcolor,pch=15,cex=1.2)
        points(parent_hypo_onset_year[i],i,col=pcolor,pch=4,cex=1.2)
        points(c(child_yob[i],child_hypo_onset_year[i]),c(i,i),col=ccolor,type='l',lwd=2)
        points(child_yob[i],i,col=ccolor,pch=15,cex=1)
        points(child_hypo_onset_year[i],i,col=ccolor,pch=4,cex=1)
    }
    dev.off()
}

if (imgs) {
    # visualization of simulation - first 30 pairs
    # Fig 1A is this with nparents=40000, ascertain_mode = 1, yob range = 1800-2000 and asc range = 1989-2013
    png(file.path(imgdir,'fig1a.png'),res=300,pointsize=3,width=500,height=500)
    plot(NA,NA,xlim=c(1870,2070),ylim=c(1,50),main='Lifespan of parent/child pairs',xlab='Year',ylab='Parent/child pair id')
    
    abline(v=1989,col='red')
    abline(v=2013,col='red')

    for (i in 1:50) {
        points(c(parent_yob[i],parent_hypo_onset_year[i]),c(i,i),col=pcolor,type='l',lwd=1.2)
        points(parent_yob[i],i,col=pcolor,pch=15,cex=1.3)
        points(parent_hypo_onset_year[i],i,col=pcolor,pch=4,cex=1.3)
        points(c(child_yob[i],child_hypo_onset_year[i]),c(i,i),col=ccolor,type='l',lwd=1)
        points(child_yob[i],i,col=ccolor,pch=15,cex=1.3)
        points(child_hypo_onset_year[i],i,col=ccolor,pch=4,cex=1.3)
    }
    mtext(side=3, line=-2, text="A", cex=2, adj=0, outer=TRUE)
    
    # legend('bottomleft',c('parents','child','birth','onset'),pch=c(16,16,15,4),col=c(pcolor,ccolor,'black','black'),cex=1)

    dev.off()
}

# distribution of year of onset in hypothetical data
if (imgs) {
    png(file.path(imgdir,'hypo.hist.yoo.png'),width=500,height=500)
    hist(c(parent_onset_year, child_onset_year), breaks=100, col='black', xlim=c(parent_yob_min, parent_yob_max+120))
    dev.off()
}

# no correlation betwen yob and ao in hypothetical data
if (imgs) {
    png(file.path(imgdir,'hypo.hist.yob.png'),width=500,height=500)
    plot(c(parent_yob, child_yob), c(parent_onset_age, child_onset_age), ylim=c(0,120))
    dev.off()
}
if (opt$fullmode) {
    all_yob = c(parent_yob, child_yob)
    all_onset = c(parent_onset_age, child_onset_age)
    m = lm(all_onset ~ all_yob)
    summary(m)
}

# no trend (no heritability) in hypothetical data
if (imgs) {
    png(file.path(imgdir,'hypo.ao.corr.png'),width=500,height=500)
    plot(parent_onset_age, child_onset_age)
    dev.off()
}
if (opt$fullmode) {
    m = lm(child_onset_age ~ parent_onset_age)
    summary(m)
}

# no difference (no anticipation) in hypothetical data
if (imgs) {
    png(file.path(imgdir,'hypo.ao.antic.barplot.png'),width=500,height=500)
    barplot(rbind(parent_onset_age, child_onset_age), beside=TRUE, col=c(pcolor,ccolor), border=NA)
    dev.off()
}
if (opt$fullmode) {
  cat("t test of all simulated individuals\n")
  t.test(parent_onset_age, child_onset_age, alternative="two.sided", paired=TRUE)
}

# now actually only ascertain the ascertainable cases
parent_onset_age[parent_onset_year > ascertainable_year_max] = NA
parent_onset_year[parent_onset_year > ascertainable_year_max] = NA
child_onset_age[child_onset_year > ascertainable_year_max] = NA
child_onset_year[child_onset_year > ascertainable_year_max] = NA

if (opt$fullmode) {
  cat("t test after right truncation of year of onset only\n")
  t.test(parent_onset_age, child_onset_age, alternative="two.sided", paired=TRUE)
  # note that this only introduces a small difference - the magnitude depends on how 
  # far back in history your simulated pairs go. when parent_yob is 1700 to 2000,
  # you get .35 years anticipation, when it's 1500 to 2000, you get .17 years anticipation,
  # when it's 1000 to 2000, you get .10 years, and so on.
}

if (ascertain_mode == 1) {
    # mode 1: you can ONLY ascertain people with onset in the given range.
    # therefore parent and child must each be in ascertainable range to be ascertained
    parent_asc = parent_onset_year > ascertainable_year_min & parent_onset_year < ascertainable_year_max # will this person be ascertained in the study
    child_asc = child_onset_year > ascertainable_year_min & child_onset_year < ascertainable_year_max
    both_asc = parent_asc & child_asc # will this parent-child _pair_ be ascertained in the study
} else if (ascertain_mode == 2) { 
    # mode 2: you FIRST ascertain people within range, then look to ascertain their family members
    # therefore parent can be ascertained if child is in range and vice versa
    parent_asc = parent_onset_year > ascertainable_year_min & parent_onset_year < ascertainable_year_max # will this person be ascertained in the study
    child_asc = child_onset_year > ascertainable_year_min & child_onset_year < ascertainable_year_max
    both_asc = parent_asc = child_asc = parent_asc | child_asc # if one is ascertained, both can be ascertained
} else if (ascertain_mode == 3) {
    # mode 3: you FIRST ascertain people within range, then look to ascertain their family members within SECOND, MORE GENEROUS range
    # therefore parent can be ascertained if child is in range and vice versa, but cases who died very long ago will not be ascertained
    parent_asc = parent_onset_year > ascertainable_year_min & parent_onset_year < ascertainable_year_max # will this person be ascertained in the study
    child_asc = child_onset_year > ascertainable_year_min & child_onset_year < ascertainable_year_max
    either_asc = parent_asc | child_asc # if one is ascertained, both can be ascertained
    parent_asc = either_asc & parent_onset_year > ascertainable_year_min_mode3
    child_asc = either_asc & child_onset_year > ascertainable_year_min_mode3
    both_asc = parent_asc & child_asc
} else if (ascertain_mode == 4) {
    # mode 4: FIRST ascertain people within range, then ascertain relatives from earlier but with a declining probability of ascertainment
    # depending on how long ago.
    parent_in_range = parent_onset_year > ascertainable_year_min & parent_onset_year < ascertainable_year_max
    parent_in_range[is.na(parent_in_range)] = FALSE
    child_in_range = child_onset_year > ascertainable_year_min & child_onset_year < ascertainable_year_max
    child_in_range[is.na(child_in_range)] = FALSE
    either_asc = parent_in_range | child_in_range
    parent_p = pmin(pmax(1.0*(parent_in_range), 
                   1.0-ascertainable_dec_rate_mode4*(ascertainable_year_min - parent_onset_year),
                   0.0,na.rm=TRUE),1.0)
    parent_asc = either_asc & (runif(n=nparents,min=0,max=1) < parent_p)
    child_p = pmin(pmax(1.0*(child_in_range),
                  1.0-ascertainable_dec_rate_mode4*(ascertainable_year_min - child_onset_year), 
                  0.0,na.rm=TRUE),1.0)
    child_asc = either_asc & (runif(n=nparents,min=0,max=1) < child_p)
    both_asc = parent_asc & child_asc
}

# convert NA to false
parent_asc[is.na(parent_asc)] = FALSE
child_asc[is.na(child_asc)] = FALSE
both_asc[is.na(both_asc)] = FALSE

if (opt$fullmode) {
  cat("t test after ascertainment, without limiting of n\n")
  t.test(parent_onset_age[both_asc], child_onset_age[both_asc], alternative="two.sided", paired=TRUE)
}




### survival analysis of ascertained data
# cannot combine this option with limit_n as currently written.
# therefore this code is run before the limit_n clause

# function to make a data frame for survival tests
# age0, age1 = vectors of ages for two groups 0 and 1
# stat0, stat1 = vectors of TRUE = died / had an event, FALSE = censored, matching age0 and age1
# groupnames = names of the two groups being compared
getsurvdf = function(age1,age2,stat1,stat2,groupnames=c(0,1)) { 
    ages = c(age1,age2) # make a single vector with times for both
    group = c(rep(groupnames[1],length(age1)),rep(groupnames[2],length(age2))) 
    status = c(stat1,stat2) 
    survdf = data.frame(ages,status,group) # convert to data frame
    return (survdf) # return the data frame
}

if (opt$fullmode) {
    # # choose pairs to include
    child_is_aw = child_intercurrent_age + child_yob > ascertainable_year_max & child_hypo_onset_age + child_yob > ascertainable_year_max
    parent_is_aw = parent_intercurrent_age + parent_yob > ascertainable_year_max & parent_hypo_onset_age + parent_yob > ascertainable_year_max
    
    parent_age = pmin(parent_hypo_onset_age, parent_intercurrent_age, ascertainable_year_max - parent_yob)
    child_age =  pmin(child_hypo_onset_age,  child_intercurrent_age,  ascertainable_year_max - child_yob)
    parent_event = parent_hypo_onset_age <= ascertainable_year_max - parent_yob & parent_intercurrent_age > parent_onset_age
    child_event  = child_hypo_onset_age  <= ascertainable_year_max - child_yob  & child_intercurrent_age  > child_onset_age
    
    
    # now try when only half of alive-and-wells are included
    # p_known = 0 # probability of an alive-and-well individual being known as such. 
    surv_results = data.frame(p_known = seq(0,1,.1), parent_median=rep(0,11), child_median=rep(0,11))
    for (p_known in seq(0,1,.1)) {
        child_known_aw = child_is_aw & sample(c(TRUE,FALSE),prob=c(p_known,1-p_known),size=length(child_is_aw),replace=TRUE)
        parent_known_aw = parent_is_aw & sample(c(TRUE,FALSE),prob=c(p_known,1-p_known),size=length(parent_is_aw),replace=TRUE)
        surv_include_half = (parent_asc & child_asc) | (parent_asc & child_known_aw) | (child_asc & parent_known_aw)
        
        survdf = getsurvdf(parent_age[surv_include_half],child_age[surv_include_half],parent_event[surv_include_half],child_event[surv_include_half])
        msurv = with(survdf,Surv(ages, status==1))
        mfit = survfit(Surv(ages,status==1) ~ group, data = survdf)
        diffobj = survdiff(Surv(ages,status==1) ~ group, data=survdf)
        p = 1-pchisq(diffobj$chisq,df=1)
        subt = '' # paste('log-rank test p value = ',formatC(p,digits=2),sep='')
        medians = summary(mfit)$table[,"median"] # extract medians from summary
        med_diff = as.integer(medians[2] - medians[1]) # take difference of medians of survival curves
      
        surv_results[surv_results$p_known==p_known,c("parent_median","child_median")] = medians
    }
    
    png(file.path(imgdir,'fig1e.new.png'),res=300,pointsize=3,width=500,height=500)
    plot(NA,NA,xlim=c(0,1),ylim=c(54,72),
        xlab='Ascertainment rate of alive and well individuals',ylab='Median age of onset',
        main='Survival of parents vs. children',yaxt='n',xaxt='n')
    axis(side=1,at=(0:10)/10,labels=paste((0:10)*10,"%",sep=""),cex.axis=.9)
    axis(side=2,at=seq(54,72,2),labels=seq(54,72,2))
    points(surv_results$p_known,surv_results$parent_median,type='b',col=pcolor,pch=18,lwd=1.2,cex=1.3)
    points(surv_results$p_known,surv_results$child_median, type='b',col=ccolor,pch=18,lwd=1.2,cex=1.3)
    legend('bottomright',c('parents','children'),col=c(pcolor,ccolor),lwd=1.2,pch=18)
    mtext(side=3, line=-2, text="E", cex=2, adj=0, outer=TRUE)
    dev.off()
}


# simple survival analysis with all available pairs.
# this is equivalent to running the above loop with p_known = 1


# include pairs with one member alive and well in 2013. don't include pairs with an intercurrent death for this
# particular analysis
# surv_include = (parent_asc & child_asc) | (parent_asc & child_is_aw) | (child_asc & parent_is_aw)
# survdf = getsurvdf(parent_age[surv_include],child_age[surv_include],parent_event[surv_include],child_event[surv_include])
# msurv = with(survdf,Surv(ages, status==1))
# mfit = survfit(Surv(ages,status==1) ~ group, data = survdf)
# diffobj = survdiff(Surv(ages,status==1) ~ group, data=survdf)
# p = 1-pchisq(diffobj$chisq,df=1)
# subt = '' # paste('log-rank test p value = ',formatC(p,digits=2),sep='')
# medians = summary(mfit)$table[,"median"] # extract medians from summary
# med_diff = as.integer(medians[2] - medians[1]) # take difference of medians of survival curves





# if desired by user, limit the final n post-ascertainment
# this feature is only fully implemented in batch mode, results unpredictable in interactive mode
if (!is.na(opt$limit_n)) {
  both_asc_indices = which(both_asc) # convert bool to indices
  both_asc_indices_limited = both_asc_indices[1:opt$limit_n] # limit number of indices
  both_asc = 1:nparents %in% both_asc_indices_limited # convert back to boolean vector
  # now handle the parent_asc and child_asc vectors. assume you get lone parents and
  # children up to the last both_asc index
  max_asc_index = max(which(both_asc))
  parent_asc[max_asc_index:nparents] = FALSE
  child_asc[max_asc_index:nparents] = FALSE
}



# now re-assess anticipation, heritability etc now that we've ascertained

if (imgs) {
    png(file.path(imgdir,'yoa.all.asc.indiv.hist.png'),width=500,height=500)
    # distribution of year of onset in all ascertained individuals
    hist(c(parent_onset_year[parent_asc], child_onset_year[child_asc]), breaks=10, col='black', xlim=c(parent_yob_min, parent_yob_max+120))
    dev.off()
    png(file.path(imgdir,'yoa.both.asc.pairs.png'),width=500,height=500)
    # or you can plot only for _pairs_ that are ascertained
    hist(c(parent_onset_year[both_asc], child_onset_year[both_asc]), breaks=10, col='black', xlim=c(parent_yob_min, parent_yob_max+120))
    dev.off()
}

# very strong correlation betwen yob and ao in ascertained data
all_yob_cens = c(parent_yob[parent_asc], child_yob[child_asc])
all_onset_cens = c(parent_onset_age[parent_asc], child_onset_age[child_asc])
m = lm(all_onset_cens ~ all_yob_cens)
if (opt$fullmode) {
    summary(m)
}

if (imgs) {
    slope = summary(m)$coefficients[2,1]
    pval  = summary(m)$coefficients[2,4]
    subtitle = '' # paste('slope of ',formatC(slope,digits=2,format='f'),' at p = ',formatC(pval,digits=2),sep='')
    png(file.path(imgdir,'fig1c.png'),res=300,pointsize=3,width=500,height=500)
    plot(c(parent_yob[parent_asc], child_yob[child_asc]), c(parent_onset_age[parent_asc], child_onset_age[child_asc]), pch=19, ylim=c(0,120),
         xlab='Year of birth', ylab='Age of onset', main='Artifactual year of birth - age of onset correlation', sub=subtitle)
    abline(m, col='red', lwd=2)
    mtext(side=3, line=-2, text="C", cex=2, adj=0, outer=TRUE)
    dev.off()
}



# show heritability in ascertained data via parent-offspring regression
# there is no correlation here in mode 2.
m = lm(child_onset_age[both_asc] ~ parent_onset_age[both_asc])
if(opt$fullmode) {
    summary(m) 
}

if (imgs) {
    png(file.path(imgdir,'fig1d.png'),res=300,pointsize=3,width=500,height=500)
    slope = summary(m)$coefficients[2,1]
    pval  = summary(m)$coefficients[2,4]
    subtitle = '' # paste('slope of ',formatC(slope,digits=2,format='f'),' at p = ',formatC(pval,digits=2),sep='')
    plot(parent_onset_age[both_asc], child_onset_age[both_asc], pch=19, xaxt='n', yaxt='n', xlab='', ylab='', main='Age of onset of ascertained pairs', sub=subtitle)
    axis(side=1,col=pcolor,col.axis=pcolor)
    mtext(side=1,col=pcolor,text=expression(bold('Parent age of onset')),padj=3)
    axis(side=2,col=ccolor,col.axis=ccolor)
    mtext(side=2,col=ccolor,text=expression(bold('Child age of onset')),padj=-3)
    abline(m,col='red')
    mtext(side=3, line=-2, text="D", cex=2, adj=0, outer=TRUE)
    dev.off()
}

if (opt$fullmode) {
    # including year of birth in model abolishes heritability
    m = lm(child_onset_age[both_asc] ~ parent_onset_age[both_asc])
    summary(m) # child_yob significant, parent_onset_age not significant

    m = lm(child_onset_age[both_asc] ~ parent_onset_age[both_asc] + child_yob[both_asc])
    summary(m) # child_yob significant, parent_onset_age not significant  
}



# large difference (anticipation) in ascertained data
if (opt$fullmode) {
    cat("t test after ascertainment, with n limits\n")
    t.test(parent_onset_age[both_asc], child_onset_age[both_asc], alternative="two.sided", paired=TRUE)  # 21 years, p = 1e-9
}
antic = t.test(parent_onset_age[both_asc], child_onset_age[both_asc], alternative="two.sided", paired=TRUE)$estimate
pval = t.test(parent_onset_age[both_asc], child_onset_age[both_asc], alternative="two.sided", paired=TRUE)$p.value
subtitle = '' # paste(formatC(antic,digits=1,format='f'),' years difference at p = ',formatC(pval,digits=2),sep='')
temp1 = rbind(parent_onset_age[both_asc], child_onset_age[both_asc])
temp1 = temp1[,!is.na(temp1[1,]) & !is.na(temp1[2,])]
if (imgs) {
    png(file.path(imgdir,'fig1b.png'),res=300,pointsize=3,width=500,height=500)
    tobj = t.test(parent_onset_age[both_asc], child_onset_age[both_asc], alternative="two.sided", paired=TRUE)  # 21 years, p = 1e-9
    barplot(temp1[,1:25], beside=TRUE, col=c(pcolor,ccolor), border=NA, xlab='Parent/child pairs', ylab='Age of onset', 
            main='Age of onset of ascertained pairs', sub=subtitle)
    legend('bottomleft',c('parents','children'),col=c(pcolor,ccolor),pch=15,bg='white')
    mtext(side=3, line=-2, text="B", cex=2, adj=0, outer=TRUE)
    dev.off()
}

if (opt$fullmode) {
    # test the robustness of observed anticipation to the stratifications 
    # applied by Pocchiari 2013
    early_birth_cohort = child_yob < 1939 & both_asc
    late_birth_cohort = child_yob >= 1939 & both_asc
    early_death_year_cohort = child_onset_year < 2000 & both_asc
    late_death_year_cohort = child_onset_year >= 2000 & both_asc
    early_death_age_cohort = child_onset_age < 61 & both_asc
    late_death_age_cohort = child_onset_age >= 61 & both_asc
    parent_early_death_age_cohort = parent_onset_age < 70 & both_asc
    parent_late_death_age_cohort = parent_onset_age >= 70 & both_asc
    cohorts = data.frame(early_birth_cohort,late_birth_cohort,early_death_year_cohort,
        late_death_year_cohort,early_death_age_cohort,late_death_age_cohort,
        parent_early_death_age_cohort,parent_late_death_age_cohort)
    for (i in 1:(dim(cohorts)[2])) {
        cohort = cohorts[,i]
        cohortname = names(cohorts)[i]
        print(paste("COHORT: ",cohortname,sep=""))
        if(length(child_onset_age[cohort]) > 1) {
            print(t.test(child_onset_age[cohort], parent_onset_age[cohort], paired=TRUE, alternative="less"))
        } else {
            print(paste("Not enough data points",asc_child_death_age))
        }
    }
}


if (opt$fullmode) {
    print(summary(mfit,times=c(40:60)),file=stdout()) # display survival results
    print(diffobj,file=stdout())
}

if (imgs) {
    png(file.path(imgdir,'fig1e.png'),res=300,pointsize=3,width=500,height=500)
    par(lend=1)
    plot(survfit(Surv(ages,status==1)~group,data=survivaldata),lwd=c(3,1.5),col=c(pcolor,ccolor),xlab='Age',ylab='Percent without onset',
         main='Survival analysis of parents vs. children', sub=subt, yaxt='n')
    axis(side=2, at=seq(0,1,.2), labels=paste(100*seq(0,1,.2),"%",sep=""))
    legend('bottomleft',c('parents','children'),col=c(pcolor,ccolor),lwd=c(3,1.5),pch=15,bg='white')
    mtext(side=3, line=-2, text="E", cex=2, adj=0, outer=TRUE)
    dev.off()
}

# prepare summary stats to print out
ttest = t.test(parent_onset_age[both_asc],child_onset_age[both_asc],paired=TRUE,alternative="two.sided")
antic = ttest$estimate
anticp = ttest$p.value
# if(ttest$p.value < .001) {
#   antic = ttest$estimate
# } else {
#   antic = 'ns'
# }
m = lm(child_onset_age[both_asc] ~ parent_onset_age[both_asc])
herit = 2*summary(m)$coefficients[2,1]
heritp = summary(m)$coefficients[2,4]
# if(summary(m)$coefficients[2,4] < .001) {
#   herit = 2*summary(m)$coefficients[2,1]
# } else {
#   herit = 'ns'
# }
m = lm(all_onset_cens ~ all_yob_cens)
yobslope = summary(m)$coefficients[2,1]
yobslopep = summary(m)$coefficients[2,4]
# if(summary(m)$coefficients[2,4] < .001) {
#   yobslope = summary(m)$coefficients[2,1]
# } else {
#   yobslope = 'ns'
# }
m = lm(child_onset_age[both_asc] ~ parent_onset_age[both_asc] + child_yob[both_asc])
herit_wyob = 2*summary(m)$coefficients[2,1]
heritp_wyob = summary(m)$coefficients[2,4]

output = data.frame(opt$npairs, opt$limit_n, opt$pyearmin, opt$pyearmax, opt$amode, opt$rate, opt$amin, opt$amax, antic, anticp, herit, heritp, yobslope, yobslopep, herit_wyob, heritp_wyob)
write.table(output, file=stdout(), sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)

# antic
# herit
# yobslope
