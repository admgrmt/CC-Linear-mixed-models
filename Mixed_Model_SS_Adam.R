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

hist(ss_data$Met_Int)  # somewhat normal, but not bimodal

## It is good practice to  standardise your explanatory variables before proceeding - you can use scale() to do that:
## this means that is is standardized toa  percentage of the distribution...

ss_data$COV2 <- scale(ss_data$COV)

## Back to our question: is Met_Int affected by COV?

###---- Fit all data in one analysis -----###

## One way to analyse this data would be to try fitting a linear model to all our data, ignoring the sites and the mountain ranges for now.

library(lme4)
library(dplyr)

basic.lm <- lm(Met_Int ~ COV2, data = ss_data)

summary(basic.lm)

## Let's plot the data with ggplot2

library(ggplot2)

ggplot(ss_data, aes(x = COV2, y = Met_Int)) +
  geom_point()+
  geom_smooth(method = "lm")


### Assumptions?

## Plot the residuals - the red line should be close to being flat, like the dashed grey line

plot(basic.lm, which = 1)  # not perfect, but look alright

## Have a quick look at the  qqplot too - point should ideally fall onto the diagonal dashed line

plot(basic.lm, which = 2)  # a bit off at the extremes, but that's often the case; again doesn't look too bad


## However, what about observation independence? Are our data independent?
## We collected multiple samples from 18 participants...
## It's perfectly plausible that the data from within each participant are more similar to each other than the data from different particiapnts - they are correlated. Pseudoreplication isn't our friend.

## Have a look at the data to see if above is true
boxplot(Met_Int ~ Subject, data = ss_data)  # certainly looks like something is going on here

## We could also plot it colouring points by subject
ggplot(ss_data, aes(x = COV2, y = Met_Int, colour = Subject))+
  geom_point(size = 2)+
  theme_classic()+
  theme(legend.position = "none")
### WANT MORE COLOR ### 
## From the above plots it looks like our mountain ranges vary both in the dragon body length and in their test scores. This confirms that our observations from within each of the ranges aren't independent. We can't ignore that.

## So what do we do?

###----- Run multiple analyses -----###


## We could run many separate analyses and fit a regression for each of the subjects.

## Lets have a quick look at the data split by subjects
## We use the facet_wrap to do that

ggplot(aes(COV2, Met_Int), data = ss_data) + geom_point() +
  facet_wrap(~ Subject) +
  xlab("COV") + ylab("Metabolic power")

##----- Modify the model -----###
## We want to use all the data, but account for the data coming from different subjects

## let's add subjects as a fixed effect to our basic.lm

Subjects.lm <- lm(Met_Int ~ COV2 + Subject, data = ss_data)
summary(Subjects.lm)

#### COV2 IS ACTUALLY VERY SIGNIFICANT
#### But some others are not significant. But let’s think about what we are doing here for a second. 
# The above model is estimating the difference in Met_int between the subjects - 
# we can see all of them in the model output returned by summary(). But we are not 
# interested in quantifying Met_Int for each specific subject: we just want
# to know whether COV2 affects Met_Int and we want to simply control for the 
# variation coming from the subjects. # This is what we refer to as 
# “random factors” and so we arrive at mixed effects models. Ta-daa!


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

## Fixed and random effects

# Fixed are variables that we expect to have an impact on the dependent variable
# (metabolic cost). They are explanatory (in SLR). Want to make
# conclusions about COV2 and Metabolic cost. COV is fixed and Met_Int is dependent

# random effects are grouping factors that we try to control. Always categorical
# we cant force R to make these effects cotniuous. We ar enot intersted in impact per
# say, but we know they might be influencing

# data for this sample could just be SAMPLE OF ALL POSSIBILITIES. E.g., with
# unlimited time and funding we could have sampled each subject across the world, but
# we are needing to generalize. We don't care that Subject 03 and 04 are differnet,
# but we know that there might be something (i.e. technique, vision, rx time) about them that is making things diferent
# and we would like to know how much variation is attributed to THEM when we predict 
# metabolic cost for people "off the street"

# Our data is just a sample of these 18 subjects. We hope that are results are genrealizable
# to people "off the street, however, we know that metabolic cost within a subject might
# be correlated, so we want to account for that. 

# IF WE CHOSE THESE SUBJECTS PRIORI, and were interested in these
# SPECIFIC PEOPLE and wanted to make predictions ABOUT THEM, it would be FITTED AS A
# FIXED effect.

# should this be fixed or random? Data is for 4 DP for each participant...
# fits each subject with their own intercept...clear

mixed.lmer2 <- lmer(Met_Int ~ COV2 + (1 + COV2|Subject), data = ss_data)
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


#So the difference between Subjects explains ~60% of the variance that is "left over"
# after the variance explained by our fixed effects...

plot(mixed.lmer)  # looks alright, no patterns evident

qqnorm(resid(mixed.lmer))
qqline(resid(mixed.lmer))  # points fall nicely onto the line - good!

#points are fitting better onto the line... but not perfectly...
# https://data.library.virginia.edu/understanding-q-q-plots/ 
# but I think we are OK...

(mm_plot <- ggplot(ss_data, aes(x = COV, y = Met_Int, colour = Subject)) +
    facet_wrap(~Subject, nrow=3) +   # a panel for each mountain range
    geom_point(alpha = 0.5) +
    theme_classic() +
    geom_line(data = cbind(ss_data, pred = predict(mixed.lmer2)), aes(y = pred), size = 1) +  # adding predicted line from mixed model 
    theme(legend.position = "none",
          panel.spacing = unit(2, "lines"))  # adding space between panels
)

# Types of random effects
# Crossed random effects
## it does not seem we have either fully or crossed effects given that each of the COV conditions was
## independent to each Subject's step length and was not

# I don't think we are also expecting any nested random effects..

# I don't know if I need to do this... but this is similar to what Dan plotted the second time...
mixed.ranslope <- lmer(Met_Int ~ COV2 + (1 + COV2|Subject), data = ss_data) 
summary(mixed.ranslope)

### plot
(mm_plot <- ggplot(ss_data, aes(x = COV, y = Met_Int, colour = Subject)) +
    facet_wrap(~Subject, nrow=3) +   # a panel for each mountain range
    geom_point(alpha = 0.5) +
    theme_classic() +
    geom_line(data = cbind(ss_data, pred = predict(mixed.ranslope)), aes(y = pred), size = 1) +  # adding predicted line from mixed model #and applies geometric but applies it individually
    theme(legend.position = "none",
          panel.spacing = unit(2, "lines"))  # adding space between panels
)



# continuing...

#Presenting your model results
# Once you get your model, you have to present it in a nicer form.
# Plotting model predictions
install.packages("ggeffects")
library(ggeffects)  # install the package first if you haven't already, then load it

# Extract the prediction data frame
# Plot the predictions 

#extract predictions from dataframe
pred.mm <- ggpredict(mixed.lmer2, terms = c("COV2"))  # this gives overall predictions for the model
# c() returns as a vector list

(ggplot(pred.mm) + 
    geom_line(aes(x = x, y = predicted)) +          # slope
    geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
                fill = "lightgrey", alpha = 0.5) +  # error band
    geom_point(data = ss_data,                      # adding the raw data (scaled values)
               aes(x = COV2, y = Met_Int, colour = Subject)) + 
    labs(x = "COV2 (indexed)", y = "Met_Int", 
         title = "COV does affect Met_Int in Subjects") + 
    theme_minimal()
)

#this is basically the slope generated along the path

#plot each of the predicted slopes 
#want to array it down the line

ggplot(ss_data, aes(x = COV, y = Met_Int, colour = Subject)) +
  facet_wrap(~Subject, nrow=3) +   # a panel for each subject
  geom_point(alpha = 0.5) +
  theme_classic() +
  geom_line(data = cbind(ss_data, pred = predict(mixed.ranslope)), aes(y = pred), size = 1) +  # adding predicted line from mixed model #and applies geometric but applies it individually
  theme(legend.position = "none",
        panel.spacing = unit(2, "lines"))  # adding space between panels

#cbind combines two rows of data together
#basically we add predicted to the end of ss_data, temporarily that gets fed in. 

temp <- cbind(ss_data, pred = predict(mixed.ranslope))
pred = predict(mixed.ranslope) ## is our 18 participants and the four data points for predicted... 
view(pred)

ggplot(ss_data, aes(x = COV, y = Met_Int, colour = Subject)) + 
  geom_line(data = cbind(ss_data, pred = predict(mixed.ranslope)), size = 1)  # adding predicted line from mixed model #and applies geometric but applies it individually

ggplot(ss_data, aes(x = COV, y = Met_Int, colour = Subject)) + 
  geom_line(data = cbind(ss_data, pred = predict(mixed.ranslope)), aes(y = pred), size = 1)  # adding predicted line from mixed model #and applies geometric but applies it individually
#this is what the predicted slopes would look like, mapped for the mixed.ranslope LMER. 
  



