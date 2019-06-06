powerFunc2Gr<- function(m1,sd1,m2,sd2,n){
  ## sample data according to inputs
  group1 <- rnorm(n/2,m1,sd1)
  group2 <- rnorm(n/2,m2,sd2)
  ## test differences between samples: is it significant?
  ts <- t.test(group1,group2)
  ps <- ts[3]  ## the p-value of the test
  ## the function will output whatever is on the last line-- here, output the pvalue
  ps  
}

## try running the function! what happens if you run it several times with the same values?
powerFunc2Gr(1,1,2,1,10)

## because we're sampling, we get different results!
## so, embed in function 'replicate' to run repeatedly.
out <- replicate(1000,powerFunc2Gr(m1=1,sd1=1,m2=2,sd2=1,n=10))
power = mean(out<0.05)

## can also get more elaborate: loop over multiple sample sizes!
## look at increasing samples by 10 over 10 steps
power <- rbind(seq(1:10)*10,rep(0,10))
rownames(power)<- c("totalN","power")

## here, we loop over s for subjects. we could loop over any of the input variables.
for(s in 1:10){
  out <- replicate(1000,powerFunc2Gr(1,1,2,1,s*10))  ##the default is that the inputs will be assigned to variables in the same order as defined when you wrote the function
  power[2,s] = mean(out<0.05)
}



### adding more conditions: same logic, different test
powerFunc3Gr<- function(m1,sd1,m2,sd2,m3,sd3,n){
  ## draw samples
  group1 <- rnorm(n/3,m1,sd1)
  group2 <- rnorm(n/3,m2,sd2)
  group3 <- rnorm(n/3,m3,sd3)
  ## format simulated data into the form the statistical model wants
  gs <- c(rep("g1",n/3), rep("g2",n/3), rep("g3",n/3))
  rs <- c(group1,group2,group3)
  ds <- as.data.frame(cbind(gs,rs))
  ds$rs <- as.numeric(as.character(ds$rs))
  ## run the test, and output only the pvalue
  lm1 <- lm(rs~gs,ds)
  ps <- unlist(anova(lm1)[5])[1]  ## this little unlist trick extracts the pvalue from a model table
  ps  
}

out <- replicate(1000,powerFunc3Gr(1,1,2,1,3,1,12))
power= mean(out<0.05)


### fully-crossed repeated measures between design with continuous outcome
## here, we'll hard-code the means and sds per group
powerFuncLMB<- function(s,i){
  ## here s=n subjects and i= n items
  # !! make sure that s and i are divisible by number of groups in the design !! 
  ## each person sees each item once
  ## observations per group, per person = i/groups
  sg <- i/4
  ## observations per group, per item = s/groups
  ig <- s/4
  ## observations per group, combining all people = (s*i)/groups
  g <- s*i/4
    
  ## set up dependant measure per group
  ##in these data, I put in 2 main effects
  ## concretely: the difference between 1 and 2 is the same as 3 and 4, and the diference between 1 and 3 is the same as 2 and 4
  group1 <- rnorm(g,m=1,sd=1)
  group2 <- rnorm(g,2,1)
  group3 <- rnorm(g,3,1)
  group4 <- rnorm(g,4,1)
  rs <- c(group1,group2,group3,group4)
  
  ## stack together the group identifiers, as a 2x2
  ## we are going to code the levels as .5 and -.5, corresponding to effects coding
  ## groups 1 and 2 belong to v1, level .5, groups 3 and 4 belong to v1, level -.5
  ## groups 1 and 3 belong to v2, level .5, groups 2 and 4 belong to v2, level -.5
  v1 <- c(rep(.5,g*2), rep(-.5,g*2))
  v2 <- rep(c(rep(.5,g), rep(-.5,g)),2)
  
  ## assign participants to observations
  ## repeat sequence of subjects sg times, and then repeat that for all 4 groups
  ## and append 's' to the numeric for ease of tracking
  ss <- paste0('s',rep(rep(seq(1:s),sg),4))
  
  ## now paste together and sort by subject, by condition
  ds <- as.data.frame(cbind(rs,v1,v2,ss))
  ds <- ds[order(ss),]
  
  ##assign items to observations
  ## here, we want to make sure we rotate through as in a latin square-- what we are setting up here is like having 4 lists in your experiment.
  ## so, participant 1 has items assigned to conditions that start at number 1
  ## participant 2 starts condition 1 at ig+1, with 1:ig assigned to condition 4
  ## pariticipant 3 stars condition 1 at 2ig+1, participant 4 starts condtion 1 at 3ig+1
  ## repeat this ig times to fill out across participants and append 'i' to beginning
  ii <- paste0('i',rep(c(1:i, (ig+1):i, 1:ig, (2*ig+1):i, 1:(2*ig), (3*ig+1):i, 1:(3*ig)),ig))
  
  ##stick it on
  ds$ii <- ii
  
  ##update types
  ds$rs <- as.numeric(as.character(ds$rs))
  ds$v1 <- as.numeric(as.character(ds$v1))
  ds$v2 <- as.numeric(as.character(ds$v2))
  
  
  lmer1 <- lmer(rs~v1*v2 + (1|ii) + (1|ss),ds,REML=F)
  
  ## to test p-values from a mixed effect model, compare to a model without the effect of interest
  ## examine that interaction: it should have power = .05 (our alpha/significance level)
  lmer2 <- lmer(rs~v1+v2 + (1|ii) + (1|ss),ds,REML=F)
  ps <- anova(lmer1,lmer2)[2,8]
  ps  
}

out <- replicate(1000,powerFuncLMB(s=12,i=12))
power= mean(out<0.05)   ###power is at just about 0.05-- what we expect!

out <- replicate(10,powerFuncLMB(s=4,i=4))
mean(out<0.05) 
## with few subjects, few items, and few replicates, the power level can be pretty divergent from what we expect
## this is why underpowered designs can reflect false positives
## make sure to run enough replicates that you're estimating the central tendency of the effect, not the extremes


## fully crossed repeated measures design with binomial outcome
## start by copying the above and make a few changes...
powerFuncBMB<- function(s,i){
  ## here s=n subjects and i= n items
  # !! make sure that s and i are divisible by number of groups in the design !! 
  ## each person sees each item once
  ## each person sees all items
  ## observations per group, per person = i/groups
  sg <- i/4
  ## observations per group, per item = s/groups
  ig <- s/4
  ## observations per group, combining all people = (s*i)/groups
  g <- s*i/4
  
  ## set up dependant measure per group: here, draw from binomial distribution
  ## get g samples from one bernoulli trial each, with overall probability = prob
  ## this will leave us with a string of zeroes and ones, with the average probability of 1 = prob
  ## here, I put in an interaction effect with no main effects.
  group1 <- rbinom(g,size=1,prob=.1)
  group2 <- rbinom(g,1,.1)
  group3 <- rbinom(g,1,.1)
  group4 <- rbinom(g,1,.5)
  rs <- c(group1,group2,group3,group4)
  
  ## stack together the group identifiers, as a 2x2
  ## we are going to code the levels as .5 and -.5, corresponding to effects coding
  ## groups 1 and 2 belong to v1, level .5, groups 3 and 4 belong to v1, level -.5
  ## groups 1 and 3 belong to v2, level .5, groups 2 and 4 belong to v2, level -.5
  v1 <- c(rep(.5,g*2), rep(-.5,g*2))
  v2 <- rep(c(rep(.5,g), rep(-.5,g)),2)
  
  ## assign participants to observations
  ## repeat sequence of subjects sg times, and then repeat that for all 4 groups
  ## and append 's' to the numeric for ease of tracking
  ss <- paste0('s',rep(rep(seq(1:s),sg),4))
  
  ## now paste together and sort by subject, by condition
  ds <- as.data.frame(cbind(rs,v1,v2,ss))
  ds <- ds[order(ss),]
  
  ##assign items to observations
  ## here, we want to make sure we rotate through as in a latin square-- what we are setting up here is like having 4 lists in your experiment.
  ## so, participant 1 has items assigned to conditions that start at number 1
  ## participant 2 starts condition 1 at ig+1, with 1:ig assigned to condition 4
  ## pariticipant 3 stars condition 1 at 2ig+1, participant 4 starts condtion 1 at 3ig+1
  ## repeat this ig times to fill out across participants and append 'i' to beginning
  ii <- paste0('i',rep(c(1:i, (ig+1):i, 1:ig, (2*ig+1):i, 1:(2*ig), (3*ig+1):i, 1:(3*ig)),ig))
  
  ##stick it on
  ds$ii <- ii
  
  ##update types
  ds$rs <- as.numeric(as.character(ds$rs))
  ds$v1 <- as.numeric(as.character(ds$v1))
  ds$v2 <- as.numeric(as.character(ds$v2))
  
  glmer1 <- glmer(rs~v1*v2 + (1|ii) + (1|ss),ds,family='binomial')
  
  ## to test p-values from a mixed effect model, compare to a model without the effect of interest
  ## examine that interaction: we want to know how many observations we need to see the effect
  glmer2 <- glmer(rs~v1+v2 + (1|ii) + (1|ss),ds,family='binomial')
  ps <- anova(glmer1,glmer2)[2,8]
  ps  
}

out <- replicate(10,powerFuncBMB(s=144,i=144)) ## this model is slow, so I am running few replicates for demo purposes.
power= mean(out<0.05)   ## power increases as s and i increase, and increases the most when you increase them both together.
