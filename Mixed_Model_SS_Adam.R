## load the data and have a look at it

setwd("~/GitHub/CC-Linear-mixed-models")
# https://ourcodingclub.github.io/tutorials/mixed-models/
library(tidyverse)
library(lme4)
library(emmeans)
library(MuMIn)
library(patchwork)
ss_data <- read_csv("stepping_stones_data.csv")


head(ss_data)

## Let's say we want to know how the body length affects test scores.

## Have a look at the data distribution:

hist(ss_data$Met_Pwr)  # somewhat normal, but not bimodal

## It is good practice to  standardise your explanatory variables before proceeding - you can use scale() to do that:
## this means that is is standardized toa  percentage of the distribution...

ss_data$COV_s <- scale(ss_data$COV)

## Back to our question: is Met_Int affected by COV?

###---- Fit all data in one analysis -----###

## One way to analyse this data would be to try fitting a linear model to all our data, ignoring the sites and the mountain ranges for now.

library(lme4)
library(dplyr)

basic.lm <- lm(Met_Pwr ~ COV, data = ss_data)

summary(basic.lm)

## Let's plot the data with ggplot2

library(ggplot2)

ggplot(ss_data, aes(x = Subject, y = Met_Pwr)) +
  geom_point()+
  geom_smooth(method = "lm")


### Assumptions?

## Plot the residuals - the red line should be close to being flat, like the dashed grey line

plot(basic.lm, which = 1)  # not perfect, but look alright

## Have a quick look at the  qqplot too - point should ideally fall onto the diagonal dashed line

plot(basic.lm, which = 2)  # a bit off at the extremes, but that's often the case; again doesn't look too bad


## However, what about observation independence? Are our data independent?
## We collected multiple samples from 18 participants...
## It's perfectly plausible that the data from within each participant are more similar to each other than the data from different participants - they are correlated. Pseudoreplication isn't our friend.

## Have a look at the data to see if above is true
boxplot(Met_Pwr ~ Subject, data = ss_data)  # certainly looks like something is going on here

# ## We could also plot it colouring points by subject
# ggplot(ss_data, aes(x = COV2, y = Met_Int, colour = Subject))+
#   geom_point(size = 2)+
#   theme_classic()+
#   theme(legend.position = "none")
# ### WANT MORE COLOR ### 
# ## From the above plots it looks like our met cost vary both in the COV and in their subject characteristics. This confirms that our observations from within each of the met cost aren't independent. We can't ignore that.

## So what do we do?

###----- Run multiple analyses -----###


## We could run many separate analyses and fit a regression for each of the subjects.

## Lets have a quick look at the data split by subjects
## We use the facet_wrap to do that

ggplot(aes(COV, Met_Pwr), data = ss_data) + geom_point() +
  facet_wrap(~ Subject) +
  xlab("COV") + ylab("Metabolic power")

##----- Modify the model -----###
## We want to use all the data, but account for the data coming from different subjects

## let's add subjects as a fixed effect to our basic.lm

#### HELP WITH INTERP: https://www.r-bloggers.com/2009/11/r-tutorial-series-simple-linear-regression/

Subjects.lm <- lm(Met_Pwr ~ COV + Subject, data = ss_data)
summary(Subjects.lm)

#### COV IS ACTUALLY VERY SIGNIFICANT, but we are uncertain whether the variability in the subject characteristics
## are influencing the values that we are getting. We want to account for that in our interpretation of the data. 
# thus we will try to interpret that variability 

###----- Mixed effects models -----###

# A mixed model is a good choice here: it will allow us to use all the data we 
# have (higher sample size) and account for the correlations between data coming 
# from the subjects. We will also estimate fewer parameters and 
# avoid problems with multiple comparisons that we would encounter while using 
# separate regressions. We are going to work in lme4, so load the package 
# (or use install.packages if you don’t have lme4 on your computer
# install.packages("sjstats")
library(lme4)
library(sjstats)

mixed.lmer <- lmer(Met_Pwr ~ COV + (1 | Subject), data = ss_data)
summary(mixed.lmer)




### this is Dan's stuff that should be ignored for now
a <- summary(mixed.lmer2) # from DAN
r_sq <- r.squaredGLMM(mixed.lmer2)
#r_sq value tells me that the correlation for estimated metabolic power for each condition of variability is 
# 0.68 for no pooling  


## will output ICC, adjusted ICC and unadjusted ICC
## MATCHES WHAT DAN WROTE ABOUT IN HIS EMAIL

# take variance for =Subjects and divide by total variance
0.03823/(0.03823+0.02212) #63%

#we can see that Subject was clearly important: it explains a lot of the variation
#we know that because we can take the variance for Subjects and divide it by the total variance


#So the difference between Subjects explains ~63% of the variance that is "left over"
# after the variance explained by our fixed effects... which is the relationship between metabolic cost 
# and the coefficient of variablity

plot(mixed.lmer)  # looks alright, no patterns evident

qqnorm(resid(mixed.lmer))
qqline(resid(mixed.lmer))  # points fall moderately nicely onto the line - good!

#points are fitting better onto the line... but not perfectly...
# https://data.library.virginia.edu/understanding-q-q-plots/ 
# but I think we are OK...


(mm_plot <- ggplot(ss_data, aes(x = COV, y = Met_Pwr, colour = Subject)) +
    facet_wrap(~Subject, nrow=3) +   # a panel for each mountain range
    geom_point(alpha = 0.5) +
    theme_classic() +
    geom_line(data = cbind(ss_data, pred = predict(mixed.lmer)), aes(y = pred), linewidth = 1) +  # adding predicted line from mixed model 
    theme(legend.position = "none",
          panel.spacing = unit(2, "lines"))  # adding space between panels
)

# data presented for each subject with an adjusted intercept but a constant slope

# https://mfviz.com/hierarchical-models/


mixed.ranslope <- lmer(Met_Pwr ~ COV + (1 + COV|Subject), data = ss_data) 
summary(mixed.ranslope)

### plot
(mm_plot <- ggplot(ss_data, aes(x = COV, y = Met_Pwr, colour = Subject)) +
    facet_wrap(~Subject, nrow=3) +   # a panel for each mountain range
    geom_point(alpha = 0.5) +
    theme_classic() +
    geom_line(data = cbind(ss_data, pred = predict(mixed.ranslope)), aes(y = pred), linewidth = 1) +  # adding predicted line from mixed model #and applies geometric but applies it individually
    theme(legend.position = "none",
          panel.spacing = unit(2, "lines"))  # adding space between panels
)



install.packages("ggeffects")
library(ggeffects)  # install the package first if you haven't already, then load it

# Extract the prediction data frame
# Plot the predictions 

#extract predictions from dataframe
pred.mm <- ggpredict(mixed.lmer, terms = c("COV"))  # this gives overall predictions for the model
# c() returns as a vector list

# final plotted slope across each subject
(ggplot(pred.mm) + 
    geom_line(aes(x = x, y = predicted)) +          # slope
    geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
                fill = "lightgrey", alpha = 0.5) +  # error band
    geom_point(data = ss_data,                      # adding the raw data (scaled values)
               aes(x = COV, y = Met_Pwr, colour = Subject)) + 
    labs(x = "COV", y = "Met_Pwr", 
         title = "COV does affect Met_Int in Subjects") + 
    theme_minimal()
)

         

#plot each of the predicted slope wrapped to each subject independently
ggplot(ss_data, aes(x = COV, y = Met_Pwr, colour = Subject)) +
  facet_wrap(~Subject, nrow=3) +   # a panel for each subject
  geom_point(alpha = 0.5) +
  theme_classic() +
  geom_line(data = cbind(ss_data, pred = predict(mixed.ranslope)), aes(y = pred), size = 1) +  # adding predicted line from mixed model #and applies geometric but applies it individually
  theme(legend.position = "none",
        panel.spacing = unit(2, "lines"))  # adding space between panels



# combining and incorperating unique slope
predslm = predict(mixed.ranslope, interval = "confidence")
head(predslm)

require(ggplot2)
install.packages("ggpubr")

install.packages("rempsyc")

pub_plot <- ggplot(ss_data, aes(x=COV, y=Met_Pwr, colour=Subject))+
  geom_point(alpha = 0.25)+
  geom_line(data = cbind(ss_data, pred = predict(mixed.ranslope)), aes(y = pred), size = 1)
  # geom_ribbon( aes(ymin = lwr, ymax = upr, fill = Subject, color = NULL), alpha = .15) 
  
pub_plot + theme_bw() +  theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


coef(mixed.ranslope)
predict(mixed.ranslope)












# #cbind combines two rows of data together
# #basically we add predicted to the end of ss_data, temporarily that gets fed in. 
# 
# temp <- cbind(ss_data, pred = predict(mixed.ranslope))
# pred = predict(mixed.ranslope) ## is our 18 participants and the four data points for predicted... 
# view(pred)
# 
# ggplot(ss_data, aes(x = COV, y = Met_Pwr, colour = Subject)) + 
#   geom_line(data = cbind(ss_data, pred = predict(mixed.ranslope)), size = 1)  # adding predicted line from mixed model #and applies geometric but applies it individually
# 
# ggplot(ss_data, aes(x = COV, y = Met_Pwr, colour = Subject)) + 
#   geom_line(data = cbind(ss_data, pred = predict(mixed.ranslope)), aes(y = pred), size = 1)
# #this is what the predicted slopes would look like, mapped for the mixed.ranslope LMER. 
#   
# 
# 
# 
# 
