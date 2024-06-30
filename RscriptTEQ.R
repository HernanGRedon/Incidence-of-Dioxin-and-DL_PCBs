###############################################################
# Mothers
###############################################################
library(deSolve)
library(future.apply)
library(future)

plan(multisession(workers = 7))

minconc <- 0
maxconc <- 32
step <- 0.01

adiposeconc <- future_sapply(X = seq(minconc,maxconc,step),FUN = function(i) {
  dioxinConcCountry <- 10
  
  # Define initial parameters
  params <- with(list(
    dose_pg_per_kgBW = i,
    f_min = 0.01,
    f_max = 0.7,
    ke_adult = 0.05,
    ke_infant = 0.05,
    K_half = 100,
    ka = 0.0025,
    wa = 0.25,
    wh = 0.03,
    absorb = 0.97,
    durBreastFed = 12,
    milkIntake = 800,
    fatContent = 0.035,
    dioxinConc = dioxinConcCountry), {
      cb_init = dioxinConcCountry / 10
      list(dose_pg_per_kgBW=dose_pg_per_kgBW, f_min=f_min, f_max=f_max, ke_adult=ke_adult, 
           ke_infant=ke_infant, K_half=K_half, ka=ka, wa=wa, wh=wh, absorb=absorb, 
           durBreastFed=durBreastFed, milkIntake=milkIntake, fatContent=fatContent, 
           dioxinConc=dioxinConc, cb_init=cb_init)
    })
  
  
  # Time vector
  times <- seq(0.00002, 420, 0.001)  # From 0 to 108 months
  
  # Initial states
  state <- c(Cbody = params$cb_init, Ch = ((params$cb_init/params$wh)*params$f_min + ((params$f_max-params$f_min)*params$cb_init)/ (params$K_half + params$cb_init)),
             Ca = ((params$cb_init/params$wa)* (1-params$f_min + ((params$f_max-params$f_min)*params$cb_init)/ (params$K_half + params$cb_init))))
  
  # Model function
  modelpbpk <- function(times, state, params) {
    with(as.list(c(state, params)), {
      
      years = round(times / 12,digits = 6)
      
      BW = 0.0006*(years + 0.083333)^3 - 0.0912*(years + 0.083333)^2 + 4.32*(years + 0.083333) + 3.6520
      BW_old = 0.0006*(years)^3 - 0.0912*(years)^2 + 4.32*years+ 3.6520
      
      intake = ifelse(times < durBreastFed,
                      yes = ((milkIntake * fatContent * dioxinConc * 30) / 1000) / BW,
                      no =  ((dose_pg_per_kgBW*BW/1000) * 30 * absorb / BW))
      
      ke = ifelse(times < 12,yes =  ke_infant,no =  ke_adult)
      
      hfactor = f_min + (((f_max - f_min) * Cbody)/ (K_half + Cbody))
      hepatic_change = ke * Cbody * hfactor
      adipose_change = ka * Ca * wa
      
      BW_change = (Cbody/BW)*(BW - BW_old)
      
      
      dCbody = intake - hepatic_change - adipose_change - BW_change
      dCh = (dCbody/wh) * hfactor
      dCa = (dCbody/wa) * (1-hfactor)
      
      #print(Ca)
      return(list(c(dCbody, dCh,dCa)))
      
      ########################################################
      
    })
  }
  # Solve the ODE
  output <- ode(y = state, times = times, func = modelpbpk, parms = params, atol = 1e-6, method = "euler")
  output <- as.data.frame(output)
  
  solution <- output["420000",]
  solution <- solution$Ca
  return(solution)
})

plan(sequential)  

adiposeconc <- as.data.frame(adiposeconc)

adiposeconc$value <- seq(minconc,maxconc,step)

writexl::write_xlsx(adiposeconc, "Adiposeconc.xlsx")




##############################################################################################################

# We used the N values as input for the TOXICR library. 
Nvalue <-33

#############################################################################################################

library(ToxicR)
# the health outcome selected is sperm concentration
healthoutcome<- matrix(0,nrow=4,ncol=2)
colnames(healthoutcome) <- c("Lower","Upper")
healthoutcome[,1] <- c(45,42.4,28.2,25) # lower
healthoutcome[,2] <- c(72.1,63.3,52.9,47.7) # higher

# standard error based on the 95%CI
# SE = (Upper[2]-Lower[1])/2*Z) Z value for 95CI% is 1.96.


STDError <- c()
for (i in 1:nrow(healthoutcome)) {
  result <- (as.numeric(healthoutcome[i,2])-as.numeric(healthoutcome[i,1]))/(2*1.96)
STDError <- c(STDError,result)
print(healthoutcome[i,2])
print(healthoutcome[i,1])
}

#Standard deviation
# STDEV = SE /sqrt(N)

STDDevfun <- function(STDError, X) {
  return(STDError/sqrt(X))
}

STDDev <- STDDevfun(STDError,Nvalue)


# Dose values are average of the reported values
M2<- matrix(0,nrow=4,ncol=4)

colnames(M2) <- c("Dose","Resp","N","StDev")
M2[,1] <- c(1,2.11,3.2,5.1) # Exposure to dioxins
M2[,2] <- c(57,51.8,38.6,34.5) # Response: Sperm concentration
M2[,3] <- as.numeric(Nvalue) #N values
M2[,4] <- as.numeric(STDDev) #standard deviation

M2

#mcmc fitting procedure
#model_fitting = ma_continuous_fit(D =  M2[,1],Y = M2[,2],fit_type = "mcmc",BMR_TYPE = "sd",BMR = 1.1,samples = 100000,threads = 7,burnin = 10000)

#print(summary(model_fitting))

#plot(model_fitting)

#the model with the highest weight is Exponential-Aerts Distribution: Log-Normal (0.336)

#Model selection for smoking, alcohol and abstinence days below 2 days, obese

#MCMC
model <- single_continuous_fit(M2[,1,drop=FALSE], M2[,2:4], BMR_TYPE="sd", BMR=1.1, ewald = TRUE,
                               distribution = "lognormal",fit_type="mcmc",model_type = "exp-aerts",samples = 1000000,threads = 7,burnin = 100000)

#plot(model)
print(summary(model))

#plot(model)

rm(list = setdiff(ls(), c("resultstable", "model")))

######################################################################################################################################################

#match the breast milk data with the adipose concentration data

#Change back to adiposeconc.xlsx
adiposeconc <- readxl::read_xlsx("Adiposeconc30years.xlsx")


###############################################################
# load concentration data
concentrationdata <-  readxl::read_xlsx("ferg2-CTTF-dioxin-TEQ-Sciensanoformat.xlsx",sheet = 1)

#selection of columns
library(dplyr)
library(tidyr)
concentrationdata <- concentrationdata %>% dplyr::select(REF_YEAR_END,REF_LOCATION,SUBSTANCE,VALUE_MEAN_TEQ,VALUE_P000_TEQ,VALUE_P100_TEQ)

#total_TEQ
concentrationdata_bycountry <- concentrationdata %>% 
  group_by(concentrationdata$REF_LOCATION,concentrationdata$REF_YEAR_END) %>%
  summarise(exposure = sum(VALUE_MEAN_TEQ))


# backcalculating the dietary exposure
concentrationdata_bycountry$estimatedexp <- c()

for (i in 1:nrow(concentrationdata_bycountry)) {
  
  concentration <- concentrationdata_bycountry[,3]
  concentration <- round(as.numeric(concentration[i,]),digits = 2)
  
  
  a <- adiposeconc %>%
    mutate(diff = abs(adiposeconc - concentration)) %>%
    filter(diff == min(diff)) %>%
    select(value)
  
  print(a)
  concentrationdata_bycountry$estimatedexp[i] <- as.numeric(a)
  
  
}

concentrationdata_bycountry$estimatedexp <- ifelse(concentrationdata_bycountry$estimatedexp == 0, 0.01, concentrationdata_bycountry$estimatedexp)

colnames(concentrationdata_bycountry) <- c("Location","Year","exposure","estimatedexp")


#plotting by region and country

WHO_regions <- readxl::read_xlsx("Copy of FERG2_Final sub-regional clusters.xlsx",sheet = 1)

concentrationdata_bycountry <-merge.data.frame(concentrationdata_bycountry,WHO_regions[,c("Country","WHO_short","Sub regional clusters FERG2")],by.x = "Location",by.y = "Country")

colnames(concentrationdata_bycountry) <- c("Location","Year","exposure","estimatedexp","ClusterA","ClusterB")


###############################################################
#Kids PBPK
###############################################################
library(deSolve)
differenceofspermconcen <- c()
serumconc <- c()

for( i in 1:nrow(concentrationdata_bycountry)) {
  
  dioxinConcCountry <- as.numeric(concentrationdata_bycountry[i,3])
  doseintake <- as.numeric(concentrationdata_bycountry[i,4])
  # Define initial parameters
  params <- with(list(
    dose_pg_per_kgBW = doseintake,
    f_min = 0.01,
    f_max = 0.7,
    ke_adult = 0.05,
    ke_infant = 0.05,
    K_half = 100,
    ka = 0.0025,
    wa = 0.25,
    wh = 0.03,
    absorb = 0.97,
    durBreastFed = 12,
    milkIntake = 800,
    fatContent = 0.035,
    dioxinConc = dioxinConcCountry), {
      cb_init = dioxinConcCountry / 10
      list(dose_pg_per_kgBW=dose_pg_per_kgBW, f_min=f_min, f_max=f_max, ke_adult=ke_adult, 
           ke_infant=ke_infant, K_half=K_half, ka=ka, wa=wa, wh=wh, absorb=absorb, 
           durBreastFed=durBreastFed, milkIntake=milkIntake, fatContent=fatContent, 
           dioxinConc=dioxinConc, cb_init=cb_init)
    })
  
  # Time vector
  times <- seq(0, 108, 0.01)  # From 0 to 108 months
  
  state <- c(Cbody = params$cb_init, Ch = ((params$cb_init/params$wh)*params$f_min + ((params$f_max-params$f_min)*params$cb_init)/ (params$K_half + params$cb_init)),
             Ca = ((params$cb_init/params$wa)* (1-params$f_min + ((params$f_max-params$f_min)*params$cb_init)/ (params$K_half + params$cb_init))))
  
  # Model function
  modelpbpk <- function(time, state, params) {
    with(as.list(c(state, params)), {
      years = time / 12
      BW = 0.00058*(years + 0.083333)^3 - 0.0948*(years + 0.083333)^2 + 4.8434*(years + 0.083333) + 2.2785
      BW_old = 0.00058*(years)^3 - 0.0948*(years)^2 + 4.8434*years+ 2.2785
      
      intake = ifelse(time < durBreastFed,
                      yes = milkIntake * fatContent * dioxinConc * 30 / 1000 / BW,
                      no =  ((dose_pg_per_kgBW*BW/1000) * 30 * absorb / BW))
      
      ke = ifelse(time < 12,yes =  ke_infant,no =  ke_adult)
      
      hfactor = f_min + (f_max - f_min) * Cbody / (K_half + Cbody)
      hepatic_change = ke * Cbody * hfactor
      adipose_change = ka * Ca * wa
      BW_change = (Cbody/BW)*(BW - BW_old)
      
      
      dCbody = intake - hepatic_change - adipose_change - BW_change
      dCh = (dCbody/wh) * (f_min + (((f_max - f_min)*Cbody)/(K_half+Cbody)))  
      dCa = (dCbody/wa) * (1-(f_min + (((f_max - f_min)*Cbody)/(K_half+Cbody))))  
      
      return(list(c(dCbody, dCh,dCa)))
    })
  }
  
  # Solve the ODE
  output <- ode(y = state, times = times, func = modelpbpk, parms = params, rtol = 1e-6, method = "euler")
  output <- as.data.frame(output)
  
  serumat8 <- output[output$time== 96,]
  serumat9 <- output[output$time== 108,]
  serumaverage8to9 <- as.numeric(serumat8[4] + serumat9[4])/2
  
  # Define the parameters
  a <- as.data.frame(model$mcmc_result)
  b <- as.data.frame(model$mcmc_result)
  c <- as.data.frame(model$mcmc_result)
  d <- as.data.frame(model$mcmc_result)
  
  a <- a[2]
  b <- b[3]
  c <- c[4]
  d <- d[5]
  
  a_2.5 <- as.numeric(exp(a$PARM_samples.1))
  b_2.5 <- as.numeric(b$PARM_samples.2)
  c_2.5 <- as.numeric(c$PARM_samples.3)
  d_2.5 <- as.numeric(d$PARM_samples.4)
  
  
  model_plotting <- function(a_2.5, b_2.5, c_2.5, d_2.5, x) {
    yvalue <- a_2.5 * (1 + ((c_2.5 - 1) * (1 - exp(-b_2.5 * (x^d_2.5)))))
    return(yvalue)
  }
  
  ci_fuction <-
    function(x) {
      c(quantile(x, probs = c(0.025, 0.975)))
    }
  
  # Calculate y values using vectorized function
  y_values_2.5 <- model_plotting(a_2.5, b_2.5, c_2.5, d_2.5, serumaverage8to9)
  
  CI <- ci_fuction(y_values_2.5-a_2.5)
  
  
  #print(y_values_97.5-a_97.5)
  
  differenceofspermconcen <- rbind(differenceofspermconcen,c(mean(y_values_2.5-a_2.5),as.numeric(CI[1]),as.numeric(CI[2])))
  serumconc <- rbind(serumconc,serumaverage8to9)
  print(mean(y_values_2.5-a_2.5))
}

differenceofspermconcen <- cbind.data.frame(concentrationdata_bycountry[,c(1:2)],differenceofspermconcen)

differenceofspermconcenpercountry <- as.vector(differenceofspermconcen$differenceofspermconcen)

colnames(differenceofspermconcen) <- c("Location","Year","Mean","Low","High")

serumconc <- cbind.data.frame(concentrationdata_bycountry[,c(1:2)],serumconc)

colnames(serumconc) <- c("Location","Year","Mean")

#plotting by region and country

WHO_regions <- readxl::read_xlsx("Copy of FERG2_Final sub-regional clusters.xlsx",sheet = 1)

differenceofspermconcen <-merge.data.frame(differenceofspermconcen,WHO_regions[,c("Country","WHO_short","Sub regional clusters FERG2")],by.x = "Location",by.y = "Country")

colnames(differenceofspermconcen) <- c("Location","Year","Mean","Low","High","clusterA","ClusterB")

serumconc <-merge.data.frame(serumconc,WHO_regions[,c("Country","WHO_short","Sub regional clusters FERG2")],by.x = "Location",by.y = "Country")

colnames(serumconc) <- c("Location","Year","Mean","clusterA","ClusterB")


#########################################################################################################
# packages
library(tidyverse)
library(mc2d)
library(ggplot2)
# settings
set.seed(123)

mean_median_ci <-
  function(x) {
    c(mean = mean(x),
      median = median(x),
      quantile(x, probs = c(0.025, 0.975)))
  }

# fit a distribution through the percentiles of common sperm concentration

library(rriskDistributions)

#fitting values for a full Male population
fittingval <- get.lnorm.par(p = c(0.025,0.05,0.1,0.24,0.5,0.75,0.9,0.95,0.975),
                            q = c(4,9,17,36,64,100,192,192,237))


differenceofspermconcen$population <- 100000

# FOR ALL population
Incidenceval <- differenceofspermconcen[,c("Location","Year","clusterA","ClusterB")]
Incidenceval$incidenceper100kmean <- 0
Incidenceval$incidenceper100k2.5 <- 0
Incidenceval$incidenceper100k97.5 <- 0

for(i in 1:nrow(differenceofspermconcen)){
  
  pnormal<- plnorm(15, meanlog = fittingval[1], sd = fittingval[2], lower.tail = T)- plnorm(15+differenceofspermconcen$Mean[i], meanlog = fittingval[1], sd = fittingval[2], lower.tail = T) 
  
  
  print(pnormal)
  incidence <- differenceofspermconcen$population[i] * pnormal
  
  Incidenceval[i,5] <- incidence
  #2.5%
  pnormal2.5<- plnorm(15, meanlog = fittingval[1], sd = fittingval[2], lower.tail = T) - plnorm(15+differenceofspermconcen$Low[i], meanlog = fittingval[1], sd = fittingval[2], lower.tail = T)  
  
  
  pnormal2.5
  Incidenceval[i,6] <- differenceofspermconcen$population[i] * pnormal2.5
  
  #97.5%
  pnormal97.5<- plnorm(15, meanlog = fittingval[1], sd = fittingval[2], lower.tail = T) - plnorm(15+differenceofspermconcen$High[i], meanlog = fittingval[1], sd = fittingval[2], lower.tail = T) 
  
  pnormal97.5
  Incidenceval[i,7] <- differenceofspermconcen$population[i] * pnormal97.5
}

Incidenceval$usualincidence <- differenceofspermconcen$population * plnorm(15, meanlog = fittingval[1], sd = fittingval[2], lower.tail = T) 


Incidenceval


#plotting by region and country

WHO_regions <- readxl::read_xlsx("Copy of FERG2_Final sub-regional clusters.xlsx",sheet = 1)

# adjust the year by +18 years accounting for the boys.

#Countries
ggplot(Incidenceval, aes(x = as.factor(Year), y = incidenceper100kmean, fill = Location)) + 
  geom_col(position = position_dodge(width = 0.9), alpha = 0.9) + 
  labs(title = "Mean by Region over Time (per 100,000 18 year old men)",
       x = "Year",
       y = "Mean Incidence") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +  # Remove the legend
  facet_wrap(~ clusterA, scales = "free_y") +
  scale_x_discrete(drop = FALSE)  # Ensure all years are shown

##################################################################################################
##################################################################################################
# Plotting Incidence

Incidenceval$Year  <- Incidenceval$Year + 18
#Africa
Incidenceval2 <- Incidenceval 

AFR <- ggplot(Incidenceval2, aes(x = as.factor(Year), y = incidenceper100kmean, fill = ClusterB)) + 
  geom_col(position = position_dodge2(preserve = "single", width = 0.9), alpha = 0.9) + 
  geom_errorbar(aes(ymin = incidenceper100k2.5, ymax = incidenceper100k97.5), 
                position = position_dodge2(preserve = "single", width = 0.9), 
                alpha = 0.9, width = 0.9) + 
  labs(title = "Mean by Region over Time (per 100,000 18 year old men)",
       x = "Year",
       y = "Mean Incidence") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  facet_wrap(~ clusterA, scales = "free_y") +
  scale_x_discrete(drop = FALSE)

AFR

ggsave(filename = "AFR_TEQ_incidence.jpeg", AFR,dpi = 600,width = 17,height = 8)

##################################################################################################
# Plotting Concentration

#Africa
concentrationdata_bycountry2 <- concentrationdata_bycountry 

AFR <-  ggplot(concentrationdata_bycountry2, aes(x = as.factor(Year), y = exposure, fill = ClusterB)) + 
  geom_col(position = position_dodge2(preserve = "single", width = 0.9), alpha = 0.9) +  
  labs(title = "Mean by Region over Time",
       x = "Year",
       y = "Mean concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  facet_wrap(~ ClusterA, scales = "free_y") +
  scale_x_discrete(drop = FALSE)


AFR

ggsave(filename = "AFR_TEQ_concen.jpeg", AFR,dpi = 600,width = 17,height = 8)

###############################################################################################
##################################################################################################
# Plotting Dietary Concentration

#Africa
concentrationdata_bycountry2 <- concentrationdata_bycountry 

AFR <-  ggplot(concentrationdata_bycountry2, aes(x = as.factor(Year), y = estimatedexp, fill = ClusterB)) + 
  geom_col(position = position_dodge2(preserve = "single", width = 0.9), alpha = 0.9) +  
  labs(title = "Mean by Region over Time",
       x = "Year",
       y = "Mean concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  facet_wrap(~ ClusterA, scales = "free_y") +
  scale_x_discrete(drop = FALSE)


AFR

ggsave(filename = "AFR_TEQ_diet_concen.jpeg", AFR,dpi = 600,width = 17,height = 8)
 
###############################################################################################
##################################################################################################
# Plotting Sperm difference

#Africa
differenceofspermconcen2 <- differenceofspermconcen

AFR <-  ggplot(differenceofspermconcen2, aes(x = as.factor(Year), y = estimatedexp, fill = ClusterB)) + 
  geom_col(position = position_dodge2(preserve = "single", width = 0.9), alpha = 0.9) +  
  labs(title = "Mean by Region over Time",
       x = "Year",
       y = "Mean concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  facet_wrap(~ clusterA, scales = "free_y") +
  scale_x_discrete(drop = FALSE)


AFR


ggsave(filename = "AFR_TEQ_diet_concen.jpeg", AFR,dpi = 600,width = 17,height = 8)


##################################################################################################

# BY Subregion


##################################################################################################
# Plotting Incidence

Incidenceval$Year  <- Incidenceval$Year + 18
#Africa
Incidenceval2 <- Incidenceval %>% dplyr::filter(clusterA == "AFR") 

AFR <- ggplot(Incidenceval2, aes(x = as.factor(Year), y = incidenceper100kmean, fill = Location)) + 
  geom_col(position = position_dodge2(preserve = "single", width = 0.9), alpha = 0.9) + 
  geom_errorbar(aes(ymin = incidenceper100k2.5, ymax = incidenceper100k97.5), 
                position = position_dodge2(preserve = "single", width = 0.9), 
                alpha = 0.9, width = 0.9) + 
  labs(title = "Mean by Region over Time (per 100,000 18 year old men)",
       x = "Year",
       y = "Mean Incidence") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  facet_wrap(~ ClusterB, scales = "free_y") +
  scale_x_discrete(drop = FALSE)

AFR

ggsave(filename = "AFR_TEQ_incidence.jpeg", AFR,dpi = 600,width = 17,height = 8)

#America
Incidenceval2 <- Incidenceval %>% dplyr::filter(clusterA == "AMR") 

AMR <- ggplot(Incidenceval2, aes(x = as.factor(Year), y = incidenceper100kmean, fill = Location)) + 
  geom_col(position = position_dodge2(preserve = "single", width = 0.9), alpha = 0.9) + 
  geom_errorbar(aes(ymin = incidenceper100k2.5, ymax = incidenceper100k97.5), 
                position = position_dodge2(preserve = "single", width = 0.9), 
                alpha = 0.9, width = 0.9) + 
  labs(title = "Mean by Region over Time (per 100,000 18 year old men)",
       x = "Year",
       y = "Mean Incidence") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  facet_wrap(~ ClusterB, scales = "free_y") +
  scale_x_discrete(drop = FALSE)

AMR

ggsave(filename = "AMR_TEQ_incidence.jpeg", AMR,dpi = 600,width = 17,height = 8)


#EMR
Incidenceval2 <- Incidenceval %>% dplyr::filter(clusterA == "EMR") 

EMR <- ggplot(Incidenceval2, aes(x = as.factor(Year), y = incidenceper100kmean, fill = Location)) + 
  geom_col(position = position_dodge2(preserve = "single", width = 0.9), alpha = 0.9) + 
  geom_errorbar(aes(ymin = incidenceper100k2.5, ymax = incidenceper100k97.5), 
                position = position_dodge2(preserve = "single", width = 0.9), 
                alpha = 0.9, width = 0.9) + 
  labs(title = "Mean by Region over Time (per 100,000 18 year old men)",
       x = "Year",
       y = "Mean Incidence") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  facet_wrap(~ ClusterB, scales = "free_y") +
  scale_x_discrete(drop = FALSE)
EMR

ggsave(filename = "EMR_TEQ_incidence.jpeg", EMR,dpi = 600,width = 17,height = 8)

#Europe
Incidenceval2 <- Incidenceval %>% dplyr::filter(clusterA == "EUR") 

EUR <- ggplot(Incidenceval2, aes(x = as.factor(Year), y = incidenceper100kmean, fill = Location)) + 
  geom_col(position = position_dodge2(preserve = "single", width = 0.9), alpha = 0.9) + 
  geom_errorbar(aes(ymin = incidenceper100k2.5, ymax = incidenceper100k97.5), 
                position = position_dodge2(preserve = "single", width = 0.9), 
                alpha = 0.9, width = 0.9) + 
  labs(title = "Mean by Region over Time (per 100,000 18 year old men)",
       x = "Year",
       y = "Mean Incidence") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  facet_wrap(~ ClusterB, scales = "free_y") +
  scale_x_discrete(drop = FALSE)

EUR

ggsave(filename = "EUR_TEQ_incidence.jpeg", EUR,dpi = 600,width = 17,height = 8)

#SEAR
Incidenceval2 <- Incidenceval %>% dplyr::filter(clusterA == "SEAR") 

SEAR <- ggplot(Incidenceval2, aes(x = as.factor(Year), y = incidenceper100kmean, fill = Location)) + 
  geom_col(position = position_dodge2(preserve = "single", width = 0.9), alpha = 0.9) + 
  geom_errorbar(aes(ymin = incidenceper100k2.5, ymax = incidenceper100k97.5), 
                position = position_dodge2(preserve = "single", width = 0.9), 
                alpha = 0.9, width = 0.9) + 
  labs(title = "Mean by Region over Time (per 100,000 18 year old men)",
       x = "Year",
       y = "Mean Incidence") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  facet_wrap(~ ClusterB, scales = "free_y") +
  scale_x_discrete(drop = FALSE)

SEAR

ggsave(filename = "SEAR_TEQ_incidence.jpeg", SEAR,dpi = 600,width = 17,height = 8)

#WPR
Incidenceval2 <- Incidenceval %>% dplyr::filter(clusterA == "WPR") 

WPR <- ggplot(Incidenceval2, aes(x = as.factor(Year), y = incidenceper100kmean, fill = Location)) + 
  geom_col(position = position_dodge2(preserve = "single", width = 0.9), alpha = 0.9) + 
  geom_errorbar(aes(ymin = incidenceper100k2.5, ymax = incidenceper100k97.5), 
                position = position_dodge2(preserve = "single", width = 0.9), 
                alpha = 0.9, width = 0.9) + 
  labs(title = "Mean by Region over Time (per 100,000 18 year old men)",
       x = "Year",
       y = "Mean Incidence") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  facet_wrap(~ ClusterB, scales = "free_y") +
  scale_x_discrete(drop = FALSE)

WPR

ggsave(filename = "WPR_TEQ_incidence.jpeg", WPR,dpi = 600,width = 17,height = 8)

##################################################################################################
# Plotting Concentration

#Africa
concentrationdata_bycountry2 <- concentrationdata_bycountry %>% dplyr::filter(ClusterA == "AFR") 

AFR <-  ggplot(concentrationdata_bycountry2, aes(x = as.factor(Year), y = exposure, fill = Location)) + 
  geom_col(position = position_dodge2(preserve = "single", width = 0.9), alpha = 0.9) +  
  labs(title = "Mean by Region over Time",
       x = "Year",
       y = "Mean concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  facet_wrap(~ ClusterB, scales = "free_y") +
  scale_x_discrete(drop = FALSE)


AFR

ggsave(filename = "AFR_TEQ_concen.jpeg", AFR,dpi = 600,width = 17,height = 8)

#America
concentrationdata_bycountry2 <- concentrationdata_bycountry %>% dplyr::filter(ClusterA == "AMR") 

AMR <-  ggplot(concentrationdata_bycountry2, aes(x = as.factor(Year), y = exposure, fill = Location)) + 
  geom_col(position = position_dodge2(preserve = "single", width = 0.9), alpha = 0.9) +  
  labs(title = "Mean by Region over Time",
       x = "Year",
       y = "Mean concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  facet_wrap(~ ClusterB, scales = "free_y") +
  scale_x_discrete(drop = FALSE)


AMR

ggsave(filename = "AMR_TEQ_concen.jpeg", AMR,dpi = 600,width = 17,height = 8)


#EMR
concentrationdata_bycountry2 <- concentrationdata_bycountry %>% dplyr::filter(ClusterA == "EMR") 

EMR <-  ggplot(concentrationdata_bycountry2, aes(x = as.factor(Year), y = exposure, fill = Location)) + 
  geom_col(position = position_dodge2(preserve = "single", width = 0.9), alpha = 0.9) +  
  labs(title =  "Mean by Region over Time",
       x = "Year",
       y = "Mean concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  facet_wrap(~ ClusterB, scales = "free_y") +
  scale_x_discrete(drop = FALSE)

EMR

ggsave(filename = "EMR_TEQ_concen.jpeg", EMR,dpi = 600,width = 17,height = 8)

#Europe
concentrationdata_bycountry2 <- concentrationdata_bycountry %>% dplyr::filter(ClusterA == "EUR") 

EUR <-  ggplot(concentrationdata_bycountry2, aes(x = as.factor(Year), y = exposure, fill = Location)) + 
  geom_col(position = position_dodge2(preserve = "single", width = 0.9), alpha = 0.9) +  
  labs(title =  "Mean by Region over Time",
       x = "Year",
       y = "Mean concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  facet_wrap(~ ClusterB, scales = "free_y") +
  scale_x_discrete(drop = FALSE)


EUR

ggsave(filename = "EUR_TEQ_concen.jpeg", EUR,dpi = 600,width = 17,height = 8)

#SEAR
concentrationdata_bycountry2 <- concentrationdata_bycountry %>% dplyr::filter(ClusterA == "SEAR") 

SEAR <-  ggplot(concentrationdata_bycountry2, aes(x = as.factor(Year), y = exposure, fill = Location)) + 
  geom_col(position = position_dodge2(preserve = "single", width = 0.9), alpha = 0.9) +  
  labs(title = "Mean by Region over Time",
       x = "Year",
       y = "Mean concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  facet_wrap(~ ClusterB, scales = "free_y") +
  scale_x_discrete(drop = FALSE)


SEAR

ggsave(filename = "SEAR_TEQ_concen.jpeg", SEAR,dpi = 600,width = 17,height = 8)

#WPR
concentrationdata_bycountry2 <- concentrationdata_bycountry %>% dplyr::filter(ClusterA == "WPR") 

WPR <-  ggplot(concentrationdata_bycountry2, aes(x = as.factor(Year), y = exposure, fill = Location)) + 
  geom_col(position = position_dodge2(preserve = "single", width = 0.9), alpha = 0.9) +  
  labs(title =  "Mean by Region over Time",
       x = "Year",
       y = "Mean concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  facet_wrap(~ ClusterB, scales = "free_y") +
  scale_x_discrete(drop = FALSE)


WPR

ggsave(filename = "WPR_TEQ_concen.jpeg", WPR,dpi = 600,width = 17,height = 8)

###############################################################################################
##################################################################################################
# Plotting Dietary Concentration

#Africa
concentrationdata_bycountry2 <- concentrationdata_bycountry %>% dplyr::filter(ClusterA == "AFR") 

AFR <-  ggplot(concentrationdata_bycountry2, aes(x = as.factor(Year), y = estimatedexp, fill = Location)) + 
  geom_col(position = position_dodge2(preserve = "single", width = 0.9), alpha = 0.9) +  
  labs(title = "Mean by Region over Time",
       x = "Year",
       y = "Mean concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  facet_wrap(~ ClusterB, scales = "free_y") +
  scale_x_discrete(drop = FALSE)


AFR

ggsave(filename = "AFR_TEQ_diet_concen.jpeg", AFR,dpi = 600,width = 17,height = 8)

#America
concentrationdata_bycountry2 <- concentrationdata_bycountry %>% dplyr::filter(ClusterA == "AMR") 

AMR <-   ggplot(concentrationdata_bycountry2, aes(x = as.factor(Year), y = estimatedexp, fill = Location)) + 
  geom_col(position = position_dodge2(preserve = "single", width = 0.9), alpha = 0.9) +  
  labs(title = "Mean by Region over Time",
       x = "Year",
       y = "Mean concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  facet_wrap(~ ClusterB, scales = "free_y") +
  scale_x_discrete(drop = FALSE)


AMR

ggsave(filename = "AMR_TEQ_diet_concen.jpeg", AMR,dpi = 600,width = 17,height = 8)


#EMR
concentrationdata_bycountry2 <- concentrationdata_bycountry %>% dplyr::filter(ClusterA == "EMR") 

EMR <-   ggplot(concentrationdata_bycountry2, aes(x = as.factor(Year), y = estimatedexp, fill = Location)) + 
  geom_col(position = position_dodge2(preserve = "single", width = 0.9), alpha = 0.9) +  
  labs(title = "Mean by Region over Time",
       x = "Year",
       y = "Mean concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  facet_wrap(~ ClusterB, scales = "free_y") +
  scale_x_discrete(drop = FALSE)

EMR

ggsave(filename = "EMR_TEQ_diet_concen.jpeg", EMR,dpi = 600,width = 17,height = 8)

#Europe
concentrationdata_bycountry2 <- concentrationdata_bycountry %>% dplyr::filter(ClusterA == "EUR") 

EUR <-   ggplot(concentrationdata_bycountry2, aes(x = as.factor(Year), y = estimatedexp, fill = Location)) + 
  geom_col(position = position_dodge2(preserve = "single", width = 0.9), alpha = 0.9) +  
  labs(title = "Mean by Region over Time",
       x = "Year",
       y = "Mean concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  facet_wrap(~ ClusterB, scales = "free_y") +
  scale_x_discrete(drop = FALSE)


EUR

ggsave(filename = "EUR_TEQ_diet_concen.jpeg", EUR,dpi = 600,width = 17,height = 8)

#SEAR
concentrationdata_bycountry2 <- concentrationdata_bycountry %>% dplyr::filter(ClusterA == "SEAR") 

SEAR <-   ggplot(concentrationdata_bycountry2, aes(x = as.factor(Year), y = estimatedexp, fill = Location)) + 
  geom_col(position = position_dodge2(preserve = "single", width = 0.9), alpha = 0.9) +  
  labs(title = "Mean by Region over Time",
       x = "Year",
       y = "Mean concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  facet_wrap(~ ClusterB, scales = "free_y") +
  scale_x_discrete(drop = FALSE)


SEAR

ggsave(filename = "SEAR_TEQ_diet_concen.jpeg", SEAR,dpi = 600,width = 17,height = 8)

#WPR
concentrationdata_bycountry2 <- concentrationdata_bycountry %>% dplyr::filter(ClusterA == "WPR") 

WPR <-   ggplot(concentrationdata_bycountry2, aes(x = as.factor(Year), y = estimatedexp, fill = Location)) + 
  geom_col(position = position_dodge2(preserve = "single", width = 0.9), alpha = 0.9) +  
  labs(title = "Mean by Region over Time",
       x = "Year",
       y = "Mean concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  facet_wrap(~ ClusterB, scales = "free_y") +
  scale_x_discrete(drop = FALSE)


WPR

ggsave(filename = "WPR_TEQ_diet_concen.jpeg", WPR,dpi = 600,width = 17,height = 8)

###############################################################################################
##################################################################################################
# Plotting Sperm difference

#Africa
differenceofspermconcen2 <- differenceofspermconcen %>% dplyr::filter(clusterA == "AFR") 

AFR <-  ggplot(differenceofspermconcen2, aes(x = as.factor(Year), y = estimatedexp, fill = Location)) + 
  geom_col(position = position_dodge2(preserve = "single", width = 0.9), alpha = 0.9) +  
  labs(title = "Mean by Region over Time",
       x = "Year",
       y = "Mean concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  facet_wrap(~ ClusterB, scales = "free_y") +
  scale_x_discrete(drop = FALSE)


AFR


ggsave(filename = "AFR_TEQ_diet_concen.jpeg", AFR,dpi = 600,width = 17,height = 8)

#America
concentrationdata_bycountry2 <- concentrationdata_bycountry %>% dplyr::filter(ClusterA == "AMR") 

AMR <-   ggplot(concentrationdata_bycountry2, aes(x = as.factor(Year), y = estimatedexp, fill = Location)) + 
  geom_col(position = position_dodge2(preserve = "single", width = 0.9), alpha = 0.9) +  
  labs(title = "Mean by Region over Time",
       x = "Year",
       y = "Mean concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  facet_wrap(~ ClusterB, scales = "free_y") +
  scale_x_discrete(drop = FALSE)


AMR

ggsave(filename = "AMR_TEQ_diet_concen.jpeg", AMR,dpi = 600,width = 17,height = 8)


#EMR
concentrationdata_bycountry2 <- concentrationdata_bycountry %>% dplyr::filter(ClusterA == "EMR") 

EMR <-   ggplot(concentrationdata_bycountry2, aes(x = as.factor(Year), y = estimatedexp, fill = Location)) + 
  geom_col(position = position_dodge2(preserve = "single", width = 0.9), alpha = 0.9) +  
  labs(title = "Mean by Region over Time",
       x = "Year",
       y = "Mean concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  facet_wrap(~ ClusterB, scales = "free_y") +
  scale_x_discrete(drop = FALSE)

EMR

ggsave(filename = "EMR_TEQ_diet_concen.jpeg", EMR,dpi = 600,width = 17,height = 8)

#Europe
concentrationdata_bycountry2 <- concentrationdata_bycountry %>% dplyr::filter(ClusterA == "EUR") 

EUR <-   ggplot(concentrationdata_bycountry2, aes(x = as.factor(Year), y = estimatedexp, fill = Location)) + 
  geom_col(position = position_dodge2(preserve = "single", width = 0.9), alpha = 0.9) +  
  labs(title = "Mean by Region over Time",
       x = "Year",
       y = "Mean concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  facet_wrap(~ ClusterB, scales = "free_y") +
  scale_x_discrete(drop = FALSE)


EUR

ggsave(filename = "EUR_TEQ_diet_concen.jpeg", EUR,dpi = 600,width = 17,height = 8)

#SEAR
concentrationdata_bycountry2 <- concentrationdata_bycountry %>% dplyr::filter(ClusterA == "SEAR") 

SEAR <-   ggplot(concentrationdata_bycountry2, aes(x = as.factor(Year), y = estimatedexp, fill = Location)) + 
  geom_col(position = position_dodge2(preserve = "single", width = 0.9), alpha = 0.9) +  
  labs(title = "Mean by Region over Time",
       x = "Year",
       y = "Mean concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  facet_wrap(~ ClusterB, scales = "free_y") +
  scale_x_discrete(drop = FALSE)


SEAR

ggsave(filename = "SEAR_TEQ_diet_concen.jpeg", SEAR,dpi = 600,width = 17,height = 8)

#WPR
concentrationdata_bycountry2 <- concentrationdata_bycountry %>% dplyr::filter(ClusterA == "WPR") 

WPR <-   ggplot(concentrationdata_bycountry2, aes(x = as.factor(Year), y = estimatedexp, fill = Location)) + 
  geom_col(position = position_dodge2(preserve = "single", width = 0.9), alpha = 0.9) +  
  labs(title = "Mean by Region over Time",
       x = "Year",
       y = "Mean concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  facet_wrap(~ ClusterB, scales = "free_y") +
  scale_x_discrete(drop = FALSE)


WPR

ggsave(filename = "WPR_TEQ_diet_concen.jpeg", WPR,dpi = 600,width = 17,height = 8)

###############################################################################################


###############################################################################################
# Load the necessary library
library(ggplot2)

# Set the parameters for the log-normal distributions
meanlog1 <- 4
meanlog2 <- 3.8
sdlog <- 0.5

# Create a data frame with values for both distributions
x <- seq(0, 50, length.out = 1000)
y1 <- dlnorm(x, meanlog1, sdlog)
y2 <- dlnorm(x, meanlog2, sdlog)

data <- data.frame(
  x = rep(x, 2),
  y = c(y1, y2),
  group = factor(rep(c("Meanlog 4", "Meanlog 3.8"), each = 1000))
)

# Calculate the mean values of the log-normal distributions
mean1 <- exp(meanlog1 + (sdlog^2) / 2)
mean2 <- exp(meanlog2 + (sdlog^2) / 2)

# Plot the distributions using ggplot2
ggplot(data, aes(x = x, y = y, color = group)) +
  geom_line(size = 1) +
  geom_area(data = subset(data, x <= 15 & group == "Meanlog 4"), aes(fill = group), alpha = 0.6) +
  geom_area(data = subset(data, x <= 15 & group == "Meanlog 3.8"), aes(fill = group), alpha = 0.2) +
  labs(title = "Sperm concentration distributions with shifted mean",
       x = "Value",
       y = "Density") +
  theme_minimal() +
  scale_color_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("blue", "red")) +
  ylim(0, max(y1, y2) * 1.2)


###############################################################################################
# Add time frame
Incidenceval$timeframe <- 0

output <- character (nrow(Incidenceval))
condition <- (Incidenceval$Year) < 2000  # condition check outside the loop
condition2 <- (Incidenceval$Year) >= 2000 & (Incidenceval$Year) < 2010
condition3 <- (Incidenceval$Year) >= 2010 

for (i in 1:nrow(Incidenceval)) {
  if (condition[i]) {
    output[i] <- "Before 2000"
  } else if (condition2[i]) {
    output[i] <- "Between 2000 and 2010"
  } else if (condition3[i]) {
    output[i] <- "After 2010"
  }
}

Incidenceval$timeframe <- output

Incidenceval <- Incidenceval %>% 
  group_by(Incidenceval$Location,Incidenceval$timeframe,clusterA,ClusterB) %>%
  summarise(incidenceper100kmean  = mean(incidenceper100kmean),
            incidenceper100k2.5  = mean(incidenceper100k2.5),
            incidenceper100k97.5  = mean(incidenceper100k97.5))

colnames(Incidenceval) <- c("Location","timeframe","clusterA","ClusterB","incidenceper100kmean","incidenceper100k97.5","incidenceper100k2.5")


library(ggplot2)
Incidenceval_byregion_clusterA <- Incidenceval %>% 
  group_by(Incidenceval$clusterA,Incidenceval$timeframe) %>%
  summarise(mean_exposure = mean(incidenceper100kmean),
            p_2.5th  = mean(incidenceper100k2.5),
            p_97.5th = mean(incidenceper100k97.5))

colnames(Incidenceval_byregion_clusterA) <- c("Region","timeframe","mean_exposure","p_2.5th","p_97.5th")

Incidenceval_byregion_clusterB <- Incidenceval %>% 
  group_by(Incidenceval$ClusterB,Incidenceval$timeframe) %>%
  summarise(mean_exposure = mean(incidenceper100kmean),
            p_2.5th  = mean(incidenceper100k2.5),
            p_97.5th = mean(incidenceper100k97.5))

colnames(Incidenceval_byregion_clusterB) <- c("Region","timeframe","mean_exposure","p_2.5th","p_97.5th")


Incidenceval_byregion_clusterA$year <- 0
Incidenceval_byregion_clusterB$year <- 0

##############################################################################
#absolute values
#Region
ggplot(Incidenceval_byregion_clusterA, aes(x = as.factor(year), y = mean_exposure, fill = Region)) + 
  geom_col(position = position_dodge(width = 0.9), alpha = 0.9) + 
  labs(title = "Mean by Region over Time ( per 100,000 18 year old men)",
       x = "",
       y = "Mean Incidence",
       fill = "WHO Region") +
  theme_minimal() +theme(axis.text.x = element_blank())+
  facet_grid(~ Region, scales = "free_y") + facet_wrap(~ timeframe) +
  scale_x_discrete(drop = FALSE)  # Ensure all years are shown

#Subregion
ggplot(Incidenceval_byregion_clusterB, aes(x = as.factor(year), y = mean_exposure, fill = Region)) + 
  geom_col(position = position_dodge(width = 0.9), alpha = 0.9) + 
  labs(title = "Mean by SubRegion over Time( per 100,000 18 year old men)",
       x = "",
       y = "Mean Incidence",
       fill = "WHO Sub-Region") +
  theme_minimal() +theme(axis.text.x = element_blank())+
  facet_grid(~ Region, scales = "free_y") + facet_wrap(~ timeframe) +
  scale_x_discrete(drop = FALSE)  # Ensure all years are shown


##############################################################################
# Serum concentration


# Add time frame
serumconc$timeframe <- 0

output <- character (nrow(serumconc))
condition <- (serumconc$Year) < 2000  # condition check outside the loop
condition2 <- (serumconc$Year) >= 2000 & (serumconc$Year) < 2010
condition3 <- (serumconc$Year) >= 2010 

for (i in 1:nrow(serumconc)) {
  if (condition[i]) {
    output[i] <- "Before 2000"
  } else if (condition2[i]) {
    output[i] <- "Between 2000 and 2010"
  } else if (condition3[i]) {
    output[i] <- "After 2010"
  }
}

serumconc$timeframe <- output


serumconc <- serumconc %>% 
  group_by(serumconc$Location,serumconc$timeframe,clusterA,ClusterB) %>%
  summarise(serumconc = mean(Mean))

colnames(serumconc) <- c("Location","timeframe","ClusterA","ClusterB","Mean")

serumconc$year <- 0

#############################################################################
#Absolute values
# Region
ggplot(serumconc, aes(x = as.factor(year), y = Mean, fill = ClusterA)) + 
  geom_col(position = position_dodge(width = 0.9), alpha = 0.9) + 
  labs(title = "Mean by Country over Time( per 100,000 18 year old men)",
       x = "",
       y = "Serum concentration",
       fill = "WHO Sub-Region") +
  theme_minimal() +theme(axis.text.x = element_blank())+
  facet_grid(~ Region, scales = "free_y") + facet_wrap(~ timeframe) +
  scale_x_discrete(drop = FALSE)  # Ensure all years are shown


# SubRegion
ggplot(serumconc, aes(x = as.factor(year), y = Mean, fill = ClusterB)) + 
  geom_col(position = position_dodge(width = 0.9), alpha = 0.9) + 
  labs(title = "Mean by Country over Time( per 100,000 18 year old men)",
       x = "",
       y = "Serum concentration",
       fill = "WHO Sub-Region") +
  theme_minimal() +
  facet_grid(~ Region, scales = "free_y") + facet_wrap(~ timeframe) +
  scale_x_discrete(drop = FALSE)  # Ensure all years are shown


