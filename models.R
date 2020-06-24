library(lme4)
library(picante)
library(MuMIn)
library(tidyverse)


birds.full <- read_csv("./BEFTA birds 2013-2019 final 12.06.2020.csv")
birds <- birds.full %>%
  # remove flyovers and observations > 70 m [to only include records in plot]
  subset(distance < 30) %>%
  # remove records of uncertain ID [doubtful records, those marked 'Unidentified', Sunda frogmouth etc had value of 1 in 'uncertain' column]
  subset(uncertain == 0) %>%
  # remove red junglefowl / domestic chicken
  subset(common != "Red Junglefowl") %>%
  #remove Eurasian tree sparrow
  subset(common != "Eurasian Tree Sparrow") %>%
  # merge the bulbuls together into one species; rename Cinereous to Japanese Tit (following recent paper)
  mutate(common = str_replace_all(common, c("Olive-winged Bulbul" = "Bulbul",
                                            "Yellow-vented Bulbul" = "Bulbul",
                                            "Sooty-headed Bulbul" = "Bulbul",
                                            "Cinereous Tit" = "Japanese Tit"))) %>%
  # adjust Latin names too
  mutate(latin = str_replace_all(latin, c("Pycnonotus plumosus" = "Pycnonotus sp.",
                                          "Pycnonotus goiavier" = "Pycnonotus sp.",
                                          "Pycnonotus aurigaster" = "Pycnonotus sp.",
                                          "Parus cinereus" = "Parus minor"))) 

#define functions

qqfunc <- function(model){
  N <- length(resid(model))
  sigma <- summary(model)$sigma
  par(mfrow=c(3,3))
  rnum <- sample(1:9,1)
  for(i in 1:(rnum-1)){
    x<- rnorm(N, 0, sigma)
    qqnorm(x,main=i)
    qqline(x)
  }
  qqnorm(resid(model), main=rnum)
  qqline(resid(model))
  for(i in (rnum+1):9){
    x<- rnorm(N, 0, sigma)
    qqnorm(x, main=i)
    qqline(x)
  }
  return(rnum)
}

panel.cor <- function(x, y, digits=1, prefix="", cex.cor = 6)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r1=cor(x,y,use="pairwise.complete.obs")
  r <- abs(cor(x, y,use="pairwise.complete.obs"))
  txt <- format(c(r1, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) { cex <- 0.9/strwidth(txt) } else {
    cex = cex.cor}
  text(0.5, 0.5, txt, cex = cex * r)
}

panel.smooth2=function (x, y, col = par("col"), bg = NA, pch = par("pch"),
                        cex = 1, col.smooth = "black", span = 2/3, iter = 3, ...)
{
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok))
    lines(stats::lowess(x[ok], y[ok], f = span, iter = iter),
          col = 1, ...)
}

panel.lines2=function (x, y, col = par("col"), bg = NA, pch = par("pch"),
                       cex = 1, ...)
{
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)){
    tmp=lm(y[ok]~x[ok])
    abline(tmp)}
}

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="white", ...)
}

### MODELS ###

model.data <- birds %>%
  mutate(stage = fct_relevel(stage, "Before", "After")) %>%
  group_by(round, plot) %>%
  mutate(n = length(common)) %>%
  mutate(sr = n_distinct(common)) %>%
  mutate(triplet = str_replace_all(plot, c("C10" = "T1",
                                           "C11" = "T1",
                                           "C12" = "T1",
                                           "C17" = "T2",
                                           "C18" = "T2",
                                           "C19" = "T2",
                                           "D28" = "T3",
                                           "D29" = "T3",
                                           "D30" = "T3",
                                           "F04" = "T4",
                                           "F05" = "T4",
                                           "F06" = "T4",
                                           "G07" = "T5",
                                           "G08" = "T5",
                                           "G09" = "T5",
                                           "G14" = "T6",
                                           "G15" = "T6",
                                           "G16" = "T6")))
view(model.data)

sr.model.0 <- glmer(sr ~ (1|triplet/plot), data=model.data, family = poisson)
summary(sr.model.0)
qqfunc(sr.model.0)
plot(sr.model.0)
AIC(sr.model.0)

sr.model.1 <- glmer(sr ~ stage*treatment + (1|triplet/plot), data=model.data, family=poisson)
summary(sr.model.1)
qqfunc(sr.model.1)
plot(sr.model.1)
AIC(sr.model.1)

n.model.0 <- glmer(n ~ (1|triplet/plot), data=model.data, family = poisson())
summary(n.model.0)
qqfunc(n.model.0)
plot(n.model.0)
AIC(n.model.0)

n.model.1 <- glmer(n ~ stage*treatment + (1|triplet/plot), data=model.data, family = poisson)
summary(n.model.1)
qqfunc(n.model.1)
plot(n.model.1)
AIC(n.model.1)

# hmm qqfunc plots all look good, but seems like variance in residuals increases with values - i.e. over-dispersion, if that's the right word? 
# check model diagnostics and Zuur et al. 2016
# try glmer.nb()
#  glmer.nb(y~x1+x2+(1/randomvariable),family=quasipoisson,data=data)
# Bolker TREE paper a while back...

# Individual species models

# I am tired - check that this is doing the right thing. 
bwp.model.data <- model.data %>%
  subset(latin == "Prinia familiaris") %>%
  mutate(count = length(latin))

view(bwp.model.data)
