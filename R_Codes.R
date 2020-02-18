#### Importing Dataset
wd <- read.csv("Dataset.csv", header = TRUE)
summary(wd)
attach(wd)
Order2 = subset (Order, Order != "PERISSODACTYLA" & Order != "PRIMATES" & Order != "XENARTHRA")
a= which(Order=="PERISSODACTYLA")
b= which(Order=="PRIMATES")
c= which(Order=="XENARTHRA")
removed=c(a, b, c)
wd2=wd[-removed,]
attach(wd2)

#########################################################################
############ Generalized Linear Model and Contrast Analysis #############
#########################################################################
library(RT4Bio)

# Ectoparasite richness
m.richness = glm.nb(Ectoparasite_richness~Order)
anova(m.richness,test="Chisq")
rdiagnostic(m.richness)
summary(m.richness)
coms(qvar="Order",mma=m.richness)
new = ifelse(wd2$Order=="RODENTIA","RODENTIA","OTHERS")
wd2$Order= new
wd2$Order
m.richness2 = glm.nb(Ectoparasite_richness~wd2$Order)
anova(m.richness2,test="Chisq")
rdiagnostic(m.richness2)
summary(m.richness2)

# Closeness centrality
m.centrality = glm(Closeness_centrality~Order,family = gaussian)
anova(m.centrality,test="F")
rdiagnostic(m.centrality)
summary(m.centrality)
coms(qvar="Order",mma=m.centrality)
m.centrality2 = glm(Closeness_centrality~wd2$Order,family = gaussian)
anova(m.centrality2,test="F")
rdiagnostic(m.centrality2)
summary(m.centrality2)


#########################################################################
############ Generalized Additive Mixed Model  ##########################
############ (Including all mammal orders) ##############################
#########################################################################

#### Importing Dataset
library(mgcv)
DATOS <- read.csv("Dataset.csv", header = TRUE)
DATOS$lweight <- log(DATOS$Weight+1)
dim(DATOS)
names(DATOS)
DATOS$Rrange <- log(1/DATOS$Proportional_range)
DATOS$Rdiversity <- log((1/DATOS$Diversity_field)+1)
DATOS$ID <- 1:nrow(DATOS)
attach(DATOS)
DATOS2 <- na.omit(DATOS)
attach(DATOS2)

# Ectoparasite richness
THETA <- (1/(var(DATOS2$Ectoparasite_richness)-mean(DATOS2$Ectoparasite_richness)))*mean(DATOS2$Ectoparasite_richness)^2
M <- gamm(Ectoparasite_richness ~ s(lweight, k=10) + s(Rrange, k=20, bs="cr") + s(Rdiversity, k=10, bs="cc"), random = list(Order= ~1+ID), family = negbin(theta=THETA), data=DATOS2)
DEV_RESIDUAL <- 2*sum(DATOS2$Ectoparasite_richness*log(DATOS2$Ectoparasite_richness/exp(predict(M$gam))))
DEV_NULL <- 2*sum(DATOS2$Ectoparasite_richness*log(DATOS2$Ectoparasite_richness/exp(mean(predict(M$gam)))))
1 - DEV_RESIDUAL/DEV_NULL #Pseudo R2
anova(M$gam)

# Closeness centrality
DATOS3 <- subset(DATOS2, DATOS2$Closeness_centrality>0)
M <- gamm(1/Closeness_centrality ~ s(lweight, k=4) + s(Rrange, k=3) + s(Rdiversity, k=4), random = list(Order= ~1+ID), family = Gamma("sqrt"), data=DATOS3)
DEV_NULL <-deviance(lm(sqrt(1/DATOS3$Closeness_centrality) ~1))
DEV_RESIDUAL <-deviance(lm(as.vector(M$gam$residuals) ~1))
1 - DEV_RESIDUAL/DEV_NULL #Pseudo R2
anova(M$gam)


#########################################################################
############ Generalized Additive Mixed Model  ##########################
############ (Only for Chiroptera and Rodentia together) ################
#########################################################################

#### Importing Dataset
library(mgcv)
DATOS <- read.csv("Dataset.csv", header = TRUE)
DATOS$lweight <- log(DATOS$Weight+1)
dim(DATOS)
names(DATOS)
DATOS$Rrange <- log(1/DATOS$Proportional_range)
DATOS$Rdiversity <- log((1/DATOS$Diversity_field)+1)
DATOS$ID <- 1:nrow(DATOS)
attach(DATOS)
DATOS2 <- na.omit(DATOS)
attach(DATOS2)
DATOS2 = subset (DATOS2, DATOS2$Order == "RODENTIA" | DATOS2$Order == "CHIROPTERA")
attach(DATOS2)

# Ectoparasite richness
THETA <- (1/(var(DATOS2$Ectoparasite_richness)-mean(DATOS2$Ectoparasite_richness)))*mean(DATOS2$Ectoparasite_richness)^2
M <- gamm(Ectoparasite_richness ~ s(lweight, k=10) + s(Rrange, k=20, bs="cr") + s(Rdiversity, k=10, bs="cc"), random = list(Order= ~1+ID), family = negbin(theta=THETA), data=DATOS2)
DEV_RESIDUAL <- 2*sum(DATOS2$Ectoparasite_richness*log(DATOS2$Ectoparasite_richness/exp(predict(M$gam))))
DEV_NULL <- 2*sum(DATOS2$Ectoparasite_richness*log(DATOS2$Ectoparasite_richness/exp(mean(predict(M$gam)))))
anova(M$gam)
1 - DEV_RESIDUAL/DEV_NULL #Pseudo R2

# Closeness centrality
DATOS3 <- subset(DATOS2, DATOS2$Closeness_centrality>0)
M <- gamm(1/Closeness_centrality ~ s(lweight, k=4) + s(Rrange, k=3) + s(Rdiversity, k=4), random = list(Order= ~1+ID), family = Gamma("sqrt"), data=DATOS3)
anova(M$gam)
DEV_NULL <-deviance(lm(sqrt(1/DATOS3$Closeness_centrality) ~1))
DEV_RESIDUAL <-deviance(lm(as.vector(M$gam$residuals) ~1))
1 - DEV_RESIDUAL/DEV_NULL #Pseudo R2

#########################################################################
############ Generalized Additive Mixed Model  ##########################
############ (Only for Rodentia) ########################################
#########################################################################

#### Importing Dataset
library(mgcv)
DATOS <- read.csv("Dataset.csv", header = TRUE)
DATOS$lweight <- log(DATOS$Weight+1)
dim(DATOS)
names(DATOS)
DATOS$Rrange <- log(1/DATOS$Proportional_range)
DATOS$Rdiversity <- log((1/DATOS$Diversity_field)+1)
DATOS$ID <- 1:nrow(DATOS)
attach(DATOS)
DATOS2 <- na.omit(DATOS)
attach(DATOS2)
DATOS2 = subset (DATOS2, DATOS2$Order == "RODENTIA")
attach(DATOS2)

# Ectoparasite richness
THETA <- (1/(var(DATOS2$Ectoparasite_richness)-mean(DATOS2$Ectoparasite_richness)))*mean(DATOS2$Ectoparasite_richness)^2
M <- gamm(Ectoparasite_richness ~ s(lweight, k=10) + s(Rrange, k=20, bs="cr") + s(Rdiversity, k=10, bs="cc"), family = negbin(theta=THETA), data=DATOS2)
DEV_RESIDUAL <- 2*sum(DATOS2$Ectoparasite_richness*log(DATOS2$Ectoparasite_richness/exp(predict(M$gam))))
DEV_NULL <- 2*sum(DATOS2$Ectoparasite_richness*log(DATOS2$Ectoparasite_richness/exp(mean(predict(M$gam)))))
anova(M$gam)
1 - DEV_RESIDUAL/DEV_NULL #Pseudo R2

# Closeness centrality
DATOS3 <- subset(DATOS2, DATOS2$Closeness_centrality>0)
M <- gamm(1/Closeness_centrality ~ s(lweight, k=4) + s(Rrange, k=3) + s(Rdiversity, k=4), random = list(Order= ~1+ID), family = Gamma("sqrt"), data=DATOS3)
anova(M$gam)
DEV_NULL <-deviance(lm(sqrt(1/DATOS3$Closeness_centrality) ~1))
DEV_RESIDUAL <-deviance(lm(as.vector(M$gam$residuals) ~1))
1 - DEV_RESIDUAL/DEV_NULL #Pseudo R2


#########################################################################
############ Generalized Additive Mixed Model  ##########################
############ (Only for Chiroptera) ######################################
#########################################################################

#### Importing Dataset
library(mgcv)
DATOS <- read.csv("Dataset.csv", header = TRUE)
DATOS$lweight <- log(DATOS$Weight+1)
dim(DATOS)
names(DATOS)
DATOS$Rrange <- log(1/DATOS$Proportional_range)
DATOS$Rdiversity <- log((1/DATOS$Diversity_field)+1)
DATOS$ID <- 1:nrow(DATOS)
attach(DATOS)
DATOS2 <- na.omit(DATOS)
attach(DATOS2)
DATOS2 = subset (DATOS2, DATOS2$Order == "CHIROPTERA")
attach(DATOS2)

# Ectoparasite richness
THETA <- (1/(var(DATOS2$Ectoparasite_richness)-mean(DATOS2$Ectoparasite_richness)))*mean(DATOS2$Ectoparasite_richness)^2
M <- gamm(Ectoparasite_richness ~ s(lweight, k=10) + s(Rrange, k=20, bs="cr") + s(Rdiversity, k=10, bs="cc"), family = negbin(theta=THETA), data=DATOS2)
DEV_RESIDUAL <- 2*sum(DATOS2$Ectoparasite_richness*log(DATOS2$Ectoparasite_richness/exp(predict(M$gam))))
DEV_NULL <- 2*sum(DATOS2$Ectoparasite_richness*log(DATOS2$Ectoparasite_richness/exp(mean(predict(M$gam)))))
anova(M$gam)
1 - DEV_RESIDUAL/DEV_NULL #Pseudo R2

# Closeness centrality
DATOS3 <- subset(DATOS2, DATOS2$Closeness_centrality>0)
M <- gamm(1/Closeness_centrality ~ s(lweight, k=4) + s(Rrange, k=3) + s(Rdiversity, k=4), random = list(Order= ~1+ID), family = Gamma("sqrt"), data=DATOS3)
anova(M$gam)
DEV_NULL <-deviance(lm(sqrt(1/DATOS3$Closeness_centrality) ~1))
DEV_RESIDUAL <-deviance(lm(as.vector(M$gam$residuals) ~1))
1 - DEV_RESIDUAL/DEV_NULL #Pseudo R2