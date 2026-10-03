library(nhanesA)
library(plyr)
library(dplyr)
library(pps)
library(survey)
library(base)
library(methods)
library(datasets)
library(utils)
library(grDevices)
library(graphics)
library(stats)
library(grid)
library(Matrix)
library(survival)
library(sampling)
library(samplingbook)
library(lattice)
library(vcd)
library(simPop)
library(caret)

## estimating probabilities

predict.nnet.1 <- function(nnet.obj, train.data, new.data, wght=nnet.obj$wts ) {
  len<- length(attr(nnet.obj$terms,"variables"))
  x.var <-  as.character(attr(nnet.obj$terms,"variables")[3:len])
  old.x <- as.data.frame(train.data[,x.var])
  new.x <- new.data[,x.var]
  new.len<- nrow(new.x)
  model.x<- model.matrix( as.formula(paste( c( "~",  nnet.obj$terms[[3]]),collapse = "" )  ),
                          data= rbind(new.x,old.x) )[1: new.len, ] 
  
  #n.layer<- length(nnet.obj$n)
  hidden.n<- (prod(nnet.obj$n[1:2])+ nnet.obj$n[2])
  hidden.x <-
    model.x %*% matrix( wght[1:hidden.n], ncol= nnet.obj$n[2])
  hidden.x <- 1/(1+exp(- hidden.x))
  out.x <-  cbind(1,  hidden.x)  %*%  matrix( wght[- (1:hidden.n)], ncol= nnet.obj$n[3])
  
  1/(1+exp(-out.x))
}  
## d f(x_i)/ d w_j 
predict.nnet.grad <- function(nnet.obj, train.data, new.data, wght=nnet.obj$wts ) {
  len<- length(attr(nnet.obj$terms,"variables"))
  x.var <-  as.character(attr(nnet.obj$terms,"variables")[3:len])
  old.x <- as.data.frame(train.data[,x.var])
  new.x <- new.data[,x.var]
  new.len<- nrow(new.x)
  model.x<- matrix( model.matrix( as.formula(paste( c( "~",  nnet.obj$terms[[3]]),collapse = "" )  ),
                                  data= rbind(new.x,old.x) )[1: new.len, ], nrow = new.len )
  
  #n.layer<- length(nnet.obj$n)
  hidden.n<- (prod(nnet.obj$n[1:2])+ nnet.obj$n[2])
  hidden.x1 <-   model.x %*% matrix( wght[1:hidden.n], ncol= nnet.obj$n[2])
  hidden.x <- 1/(1+exp(- hidden.x1))
  out.x <-  cbind(1,  hidden.x)  %*%  matrix( wght[- (1:hidden.n)], ncol= nnet.obj$n[3])
  
  #1/(1+exp(-out.x))
  
  temp.hidden  <- t( t( exp(- hidden.x1)/(1+exp(- hidden.x1))^2) *  wght[- (1:hidden.n)][-1] )
  temp.hidden1 <-  model.x*temp.hidden[,1]
  i<-2
  while( i <=  nnet.obj$n[2]) {
    temp.hidden1 <- cbind(temp.hidden1,
                          model.x*temp.hidden[,i])
    i<-i+1
  } 
  
  list( grad =   c(exp(-out.x)/(1+exp(-out.x))^2) * cbind(  
    temp.hidden1, 1,  hidden.x), 
    raw  =  1/(1+exp(-out.x))
  ) 
  
}
### confidence interval for ANN estimates 
Prev.CI.nnet<-function( nnet.obj, train.data, new.data, resp.var="Disease1", coverage= c(0.025,0.975)){
  
  alpha<- nnet.obj$decay
  temp.J1<- predict.nnet.grad(nnet.obj,train.data=train.data, new.data= train.data)
  
  J1<-   temp.J1$grad
  
  ## (21) 
  #H1<-    J1%*%solve(t(J1)%*%J1+alpha*diag(length(nnet.obj$wts)))%*%t(J1)
  
  temp.JJ <- t(J1)%*%J1
  temp.JJ.1 <- solve(temp.JJ+alpha*diag(length(nnet.obj$wts)))
  
  # RSS <-  sum( (train.data[,c(resp.var)] - predict(nnet.obj, type="raw"))^2)
  RSS <-  sum( (train.data[,c(resp.var)] - temp.J1$raw)^2)
  ## s^2   (23)
  train.n<-dim(train.data)[1]  ## number of obs in training 
  s.2 <- RSS/( train.n - 2*sum( diag( temp.JJ.1%*%temp.JJ ) ) +
                 sum( diag( temp.JJ.1%*%temp.JJ%*%temp.JJ.1%*%temp.JJ ) )  )
  
  ## var(y|x) (20)
  temp.g<- predict.nnet.grad(nnet.obj,train.data=train.data,  new.data= new.data, wght=nnet.obj$wts )
  g0<- temp.g$grad
  
  nn <- dim(new.data)[1]
  ## var(\bar y|x) = sum (var(y|x))/ nn^2 (2) 
  mean.g<-  colMeans(g0[1:nn,])
  var.avg.y<- s.2 * ( 1/nn  + mean.g%*%temp.JJ.1%*%temp.JJ%*%temp.JJ.1%*%(mean.g))
  sd.avg.y<- sqrt(var.avg.y)
  
  # Prev.nnet.noweight <- predict(nnet.obj, newdata = new.data, type="raw")
  
  Prev.nnet.noweight <- temp.g$raw
  
  mean(Prev.nnet.noweight) + qnorm(coverage)* c(sd.avg.y)
  
}




# 
# setwd("D:/GoogleDriveOakland/MLSurvey")
# 
# load("G:/My Drive/MLSurvey/CLSAPpl1.Rdata")
# 
# rm(list=setdiff(ls(), "Ppl.Data"))



### how to generate the population used 
# ### work on a smaller population dataset
# set.seed(112300)
# Ppl.Data2<-  Ppl.Data  %>%   
#   group_by(Sampling_Strata) %>%
#   slice_sample(prop=0.1)
# 
# # save(Ppl.Data2, file="USPpl2.Rdata")
# ## remove Ppl.Data and  collection garbage to save space 
# rm( "Ppl.Data")
# gc()


############################################## start from here #####################################################
setwd("G:/My Drive/MLSurvey")

load( file="USPpl2.Rdata")
#rm(list=setdiff(ls(), "Ppl.Data2"))



# Ppl.Data2$Ppl_Size_State <-Ppl.Data2$Ppl_Size_Prov
# #Ppl.Data2$Ppl_Size_Prov <-NULL
# Ppl.Data2$Sample.size.state<-Ppl.Data2$Sample.size.Prov
# #Ppl.Data2$Sample.size.Prov<-NULL
# Ppl.Data2$WGHTS_state_TRM<-Ppl.Data2$ WGHTS_PROV_TRM
# #Ppl.Data2$ WGHTS_PROV_TRM<-NULL


##### sampling information 

Prov.vec   <- c( "AB", "BC", "MB","NB","NL","NS","ON","PE","QC","SK")
Prov.DCS   <- c( "AB", "BC", "MB",     "NL","NS","ON",     "QC") 
Prov.NonDCS<- c( "NB", "PE", "SK")      
# state.prop  <- c(71000,79000,32000,14000,11000,16000,82000,2800,2e+05,12000)

Age.grp    <- c("45-54","55-64","65-74","75-85")
Sex.grp    <- c("female","male")



Sampling.mat<-expand.grid(Prov=Prov.vec, DCS.vec  =c("DCS","Non_DCS")
                          ,SEX =c("F","M"), Age.grp = Age.grp )  

Sampling.mat<-Sampling.mat[-which(Sampling.mat$Prov%in% Prov.NonDCS 
                                  & Sampling.mat$DCS.vec=="DCS"),]
# rownames(Sampling.mat)<-1:136
# state.Size   <- c(  IL= 90 ,IN=100,IA=90,MD=60,MA=60,MI=60,NJ=85,NM=55,NC=200,OH=50)
# DCS.Size    <- c( DCS= 420, NON.DCS = 280)/700
# #Sex.Size <- c(   F= 510,       M = 350)/860
# Sex.Size <- c(   F= 500,       M = 490)/990

rownames(Sampling.mat)<-1:136
Prov.Size   <- c(  AB= 90 , BC=100,MB=90,NB=60,NL=60,NS=60,ON=85,PE=55,QC=200,SK=50)
DCS.Size    <- c( DCS= 420, NON.DCS = 280)/700
Sex.Size <- c(   F= 500,       M = 490)/990

### sample size in each sub strata  
Sampling.mat$Sample_size<- round( with(Sampling.mat,Prov.Size[Prov]*
                                         Sex.Size[SEX] * DCS.Size[DCS.vec]^(1-Prov%in% Prov.NonDCS)  ) )
Sampling.mat$Sampling_Strata<-  with(Sampling.mat, 
                                     paste(Prov, DCS.vec, SEX, Age.grp,  sep="_") )
temp.vec <- table(Ppl.Data2$Sampling_Strata)   ### calculate the population size strata  here 
Sampling.mat<-merge(Sampling.mat
                    , data.frame(Sampling_Strata=names(temp.vec ), N_h=c(temp.vec) )
                    , by="Sampling_Strata" ) 



library(dplyr)



## sample sizes 
Sampling.mat$n_h <- Sampling.mat$Sample_size
## matching varible 
Sampling.mat$GEOSTRAT_TRM <- with(Sampling.mat, paste( Prov, DCS.vec,sep="_") ) 

## basic design weight, initial weight
Sampling.mat$basic_weight <- Sampling.mat$N_h /Sampling.mat$n_h 
head(Sampling.mat)





library(dplyr)

Temp<- Ppl.Data2 %>% dplyr::group_by( Sampling_Strata,WGHTS_PROV_TRM ) %>% dplyr::summarise( sample.size = max(Sample_size)) %>% 
  dplyr::group_by( WGHTS_PROV_TRM ) %>%  dplyr::summarize( Sample.size.state = sum(sample.size)*1)

Ppl.Data2<-  dplyr::left_join(Ppl.Data2,  Temp, by = "WGHTS_PROV_TRM")

Ppl.Data2<- Ppl.Data2 %>%    dplyr::group_by(WGHTS_PROV_TRM) %>%  dplyr::mutate(Ppl_Size_PROV= dplyr::n() ) %>%  dplyr::ungroup()

###########################Sampling Designed ended; Now work on Population disease ##########################################################
Ppl.Data2$WGHTS_state_TRM <- Ppl.Data2$WGHTS_PROV_TRM 
### modify the Ppl first before the sampling: 
### Creat a scenario to predict something:  

set.seed(23471)
Ppl.Data2$Gen1 <- rbinom( dim(Ppl.Data2)[1], 1, prob= 0.1 + 0.2*Ppl.Data2$WGHTS_state_TRM %in% c("BC","ON","QC") )
## smoking @ 2017: https://uwaterloo.ca/tobacco-use-canada/adult-tobacco-use/smoking-provinces
smoking.prov.vec <- c( "BC", "AB", "SK","MB", "ON","QC", "NB", "NS","PE","NL" )
smoking.prevalence<-c(15.6,18.9,17.8,14.5,12.9,15.7,13.7,18.5,11.8,20.1)/100
names(smoking.prevalence)<-smoking.prov.vec

## order it to alphabetical order
smoking.prevalence<-smoking.prevalence[order(smoking.prov.vec)]
smoking.prov.vec  <-smoking.prov.vec[order(smoking.prov.vec)]

Ppl.Data2$Smoking <- rbinom( dim(Ppl.Data2)[1], 1, prob=smoking.prevalence[Ppl.Data2$WGHTS_state_TRM])
expit<- function(x){ 1/(1+exp(-x))}

Ppl.Data2$Disease1 <- rbinom( dim(Ppl.Data2)[1], 1, prob= expit(-2 + Ppl.Data2$Gen1
                                                                +  (Ppl.Data2$SEX_ASK_TRM=="M") +  Ppl.Data2$Smoking  )  )
## prevalence of disease one (our target)
mean(Ppl.Data2$Disease1)



###########################Population disease done; Simulation begin  ##########################################################
library(survey)
logit<-function(p){ log(p/(1-p))}
expit<-function(x) {1/(1+exp(-x))}
dexpit<-function(x) {exp(-x)/(1+exp(-x))^2}
library(caret)

for(Sim.index in 1: 1000){
  print( paste("Now working on Sim.index=",Sim.index, sep=""))
  set.seed(Sim.index)
  strat_sample <- Ppl.Data2 %>% dplyr::group_by(WGHTS_state_TRM) %>%dplyr:: sample_n(Sample.size.state)
  
  strat_sample$ID2<- strat_sample$entity_id
  #slice_sample(n=Sample_size)
  # Change  the sampling vector
  strat_sample$Sampling_Strata <-  strat_sample$WGHTS_state_TRM 
  
  #Sample.Data <- Ppl.Data %>% group_by(Sampling_Strata) %>% sample_n(Sample_size)
  
  ## verify if they are the same (difference should be between -1,1)
  
  table(strat_sample$Sampling_Strata)- matrix( as.numeric( 
    strat_sample %>%    group_by(Sampling_Strata) %>% 
      summarise( mean( Sample.size.state ) ) %>% as_tibble %>% unlist), ncol=2  )[,2] 
  
  
  
  ## inflation weight
  strat_sample$inflation.weight<-  strat_sample$basic_weight
  ## analytic weight
  # strat_sample <- strat_sample %>%  left_join(
  #   strat_sample %>% group_by(Sampling_Strata) %>% summarise(sum.analytic.weight = sum(inflation.weight))
  #   , by= "Sampling_Strata" )%>%  
  #   mutate( analytic.weight = inflation.weight/ sum.analytic.weight* Sample_size )  
  
  
  
  # Logistics GREG estimators ( need to modify it)
  #############################################################################################################
  #http://www.asasrms.org/Proceedings/y2010/Files/308721_61470.pdf
  
  
  
  ## prevalence (our target)
  # mean(Ppl.Data2$Disease1)
  
  prev2<-mean(Ppl.Data2$Disease1);  
  # log( prev2/(1-prev2) )
  
  
  # strat_sample$Disease1<- as.numeric( strat_sample$Disease1)
  
  svy.degn  <-  svydesign( ids= ~entity_id  , strata = ~ Sampling_Strata, 
                           weights = ~inflation.weight, data =  strat_sample )
  
  svy.logit.model <-svyglm( Disease1 ~ Smoking + SEX_ASK_TRM   + WGHTS_state_TRM + AGE_NMBR_TRM+ DCS.groups + ED_UDR11_TRM+Gen1,  
                            family= quasibinomial, design=svy.degn)
  
  
  logit.model <- glm( Disease1 ~ Smoking + SEX_ASK_TRM   + WGHTS_state_TRM + AGE_NMBR_TRM+ DCS.groups + ED_UDR11_TRM+Gen1, 
                      family=binomial(link = "logit"), data=strat_sample)
  
  
  ## the SE estimate of the linear predictor average
  
  # temp<-cbind( 1, Ppl.Data2[-strat_sample$ID2, c("Smoking" , "SEX_ASK_TRM", "WGHTS_state_TRM", "AGE_NMBR_TRM", "DCS.groups", "ED_UDR11_TRM")] )
  # 
  temp <- model.matrix(~ Smoking + SEX_ASK_TRM   + WGHTS_state_TRM + AGE_NMBR_TRM+ DCS.groups + ED_UDR11_TRM+Gen1,  
                       Ppl.Data2[-strat_sample$ID2, ])
  
  # temp$SEX_ASK_TRM<- as.numeric(temp$SEX_ASK_TRM == "F")
  # 
  # temp<- as.matrix(temp)
  ## the point estimate by averaging linear predictor
  (est1.pt <- mean( expit( temp%*%coef(logit.model) )  ))
  
  #sqrt(diag(temp%*%vcov(logit.model)%*%t(temp))) ## std error of the linear predictor 
  
  (est1.se <-c(sqrt( colMeans( c(dexpit(temp%*%coef(logit.model))) * temp ) %*%vcov(logit.model)  %*%  colMeans( c( dexpit(temp%*%coef(logit.model))) * temp ) )))
  
  #Not this one (est1.se <- as.numeric( sqrt( colSums(temp)%*%vcov(logit.model)%*%colSums(temp))) ) 
  
  ## the Wald CI of the prevalence
  (est1.CI <-  (est1.pt+ qnorm(p=c(0.025,0.975))*est1.se ) )
  
  
  
  ## 5.1.1 Model-Unbiased Prediction Estimator for y
  
  svy.logit.pd<-predict( svy.logit.model, newdata = Ppl.Data2[-strat_sample$ID2 ,]) 
  
  # sqrt(diag(temp[1:2,]%*%vcov(svy.logit.model)%*%t(temp[1:2,]))) ## std error of the linear predictor 
  # SE(svy.logit.pd)[1:2]
  
  ## the point estimate by averaging linear predictor
  (svy.est1.pt <- mean( expit(svy.logit.pd)))
  
  ## the SE estimate of the linear predictor average  (by  deltas method )
  (svy.est1.se <-c(sqrt( colMeans( c(dexpit(temp%*%coef(logit.model))) * temp ) %*%vcov(svy.logit.model)  %*%  colMeans( c( dexpit(temp%*%coef(logit.model))) * temp ) )))
  
  #Do not run this one: (svy.est1.se <-  sqrt( colSums(temp)%*%vcov(svy.logit.model)%*%colSums(temp))  )
  
  ## the Wald CI of the prevalence
  (svy.est1.CI <- svy.est1.pt + qnorm(p=c(0.025,0.975))*svy.est1.se ) 
  
  ## model assisted to estimate the population total ( mean ) (7)
  # strat_sample$residue<- strat_sample$Disease1 - expit(predict( svy.logit.model, newdata = strat_sample))
  
  # (DIFF <- (sum( expit(predict( svy.logit.model, newdata = Ppl.Data2))) +
  #             sum( strat_sample$inflation.weight* strat_sample$residue))/nrow(Ppl.Data2))
  # 
  ## model assisted to estimate the population total ( mean ) variance (8)
  ## re define the design again  
  
  # svy.degn  <-  svydesign( ids= ~entity_id  , strata = ~ Sampling_Strata, 
  #                          weights = ~inflation.weight, data =  strat_sample )
  # 
  # SE.DIFF<- SE( svymean( ~residue, design = svy.degn)  )
  
  ## the Wald CI of the prevalence
  # (svy.DIFF.CI <- DIFF + qnorm(p=c(0.025,0.975))*c(SE.DIFF) ) 
  
  ## svy ratio (ignoring logistic regression) not good 
  ##svymean( ~Disease1, design = svy.degn) + qnorm(p=c(0.025,0.975))*c( SE(svymean( ~Disease1, design = svy.degn)))
  
  ## prevalence
  mean(Ppl.Data2$Disease1)
  ######################################### add the ann  stuffs here  ################################################
  
  if(Sim.index ==1  ){
    strat_sample$Disease1C<- paste( "Disease",strat_sample$Disease1,sep="")
    
    fitControl <- trainControl(method = "repeatedcv", 
                               number = 10, #10, 
                               repeats = 5, #5, 
                               classProbs = TRUE, 
                               summaryFunction = twoClassSummary)
    
    nnetGrid <-  expand.grid(size = seq(from = 1, to = 10, by =5),
                             decay = seq(from = 0.1, to = 0.5, by = 0.3))
    
    
    
    nnetFit <- train(Disease1C ~ Smoking + SEX_ASK_TRM   + WGHTS_state_TRM + AGE_NMBR_TRM+ DCS.groups + ED_UDR11_TRM+Gen1, 
                     data = strat_sample,
                     method = "nnet", trace = FALSE, 
                     metric = "ROC",
                     trControl = fitControl,
                     tuneGrid = nnetGrid,
                     verbose = FALSE)
    save(nnetFit, file="nnetFit.RData")
    
  }
  
  load( file="nnetFit.RData")
  strat_sample$Disease1<-  as.numeric(sub("Disease","",as.character(strat_sample$Disease1)) )
  
  nnetFit.final <- nnet::nnet( (Disease1) ~  Smoking + SEX_ASK_TRM   + WGHTS_state_TRM + AGE_NMBR_TRM+ DCS.groups + ED_UDR11_TRM+Gen1,
                               data = strat_sample, 
                               Hess= TRUE,
                               trace = FALSE,  ##   false is faster  but print out  
                               size = nnetFit$bestTune$size,  decay = nnetFit$bestTune$decay, maxit = 1000 )
  
  ANN.CI<- Prev.CI.nnet(nnetFit.final,train.data=strat_sample, 
                        new.data= Ppl.Data2 
                        ,resp.var="Disease1")
  
  ANN.prev <- mean( ANN.CI)
  ANN.se  <-  diff( ANN.CI)/2/qnorm(0.975)
  
  
  
  ######################################### add the svy.ann  stuffs here  ################################################ 
  
  if(Sim.index ==1  ){
    strat_sample$Disease1C<- paste( "Disease",strat_sample$Disease1,sep="")
    
    fitControl <- trainControl(method = "repeatedcv", 
                               number = 10, #10, 
                               repeats = 5, #5, 
                               classProbs = TRUE, 
                               summaryFunction = twoClassSummary)
    
    nnetGrid <-  expand.grid(size = seq(from = 1, to = 10, by =5),
                             decay = seq(from = 0.1, to = 0.5, by = 0.3))
    
    
    nnetFit.wt <- train(Disease1C ~ Smoking + SEX_ASK_TRM   + WGHTS_state_TRM + AGE_NMBR_TRM+ DCS.groups + ED_UDR11_TRM+Gen1,
                        data = strat_sample,
                        method = "nnet",trace = FALSE, 
                        metric = "ROC",
                        weights =  strat_sample$inflation.weight,
                        trControl = fitControl,
                        tuneGrid = nnetGrid,
                        verbose = FALSE)
    
    
    save(nnetFit.wt, file="nnetFitWt.RData")
  }
  
  load( file="nnetFitWt.RData")
  b.size <- as.numeric(nnetFit.wt$bestTune[1])
  b.decay<- as.numeric(nnetFit.wt$bestTune[2])
  strat_sample$Disease1<-  as.numeric(sub("Disease","",as.character(strat_sample$Disease1)) )
  nnetFit.wt1<- nnet::nnet( (Disease1) ~ Smoking + SEX_ASK_TRM   + WGHTS_state_TRM + AGE_NMBR_TRM+ DCS.groups + ED_UDR11_TRM+Gen1,
                            data=strat_sample,   
                            weights =  strat_sample$inflation.weight,
                            trace=FALSE,
                            size=b.size, decay=b.decay, maxit=1000)
  
  strat_sample$residue<-as.numeric(strat_sample$Disease1)- 
    predict(nnetFit.wt1, newdata =  strat_sample,type="raw")
  
  (DIFF.nnet <- (sum( predict(nnetFit.wt1, newdata =  Ppl.Data2,type="raw")) +
                   sum( strat_sample$inflation.weight* strat_sample$residue))/nrow(Ppl.Data2))
  
  ## model assisted to estimate the population total ( mean ) variance (8)
  ## re define the design again  
  
  svy.degn  <-  svydesign( ids= ~entity_id  , strata = ~ Sampling_Strata, 
                           weights = ~inflation.weight, data =  strat_sample )
  
  SE.DIFF.nnet<- SE( svymean( ~residue, design = svy.degn)  )
  
  ## the Wald CI of the prevalence
  (svy.DIFF.nnet.CI <- DIFF.nnet + qnorm(p=c(0.025,0.975))*c(SE.DIFF.nnet) )   #0.2444929 0.2615878
  
  
  #################################################
  svy.degn  <-  svydesign( ids= ~entity_id  , strata = ~ Sampling_Strata, 
                           weights = ~inflation.weight, data =  strat_sample )
  ## svy ratio (ignoring logistic regression) not good 
  svy.mean.mth  <- svymean( ~Disease1, design = svy.degn)
  svy.mean.prev <- c( svy.mean.mth[1])
  svy.mean.se   <- c(SE(svy.mean.mth))
  svy.mean.CI <- svy.mean.prev + qnorm(p=c(0.025,0.975))*c(svy.mean.se)
  
  ################ Section for saving the results   ####################################
  # Sim.index<- 1
  if (Sim.index ==1 ){
    Name.vec<- c("Sim.index", "logistic.prev" ,"logistic.se"    , "logist.CI.L" , "logist.CI.U", 
                 "svy.logit.prev","svy.logistic.se",  "svy.logist.CI.L" , "svy.logist.CI.U",
                 
                 "ann.prev" ,"ann.se"    , "ann.CI.L" , "ann.CI.U",
                 "svy.ann.prev" ,"svy.ann.se"    , "svy.ann.CI.L" , "svy.ann.CI.U",
                 "svy.mean.prev" ,"svy.mean.se"    , "svy.mean.CI.L" , "svy.mean.CI.U")
    write.table( t(Name.vec), file= "result.csv",  append = FALSE, sep=",",col.names = FALSE)
  }
  
  value.vec<- c(Sim.index, est1.pt, est1.se , est1.CI,
                svy.est1.pt, svy.est1.se,svy.est1.CI,
                ANN.prev, ANN.se,ANN.CI,
                DIFF.nnet , SE.DIFF.nnet, svy.DIFF.nnet.CI,
                svy.mean.prev, svy.mean.se, svy.mean.CI )
  ### only run for the first time
  
  write.table( t(value.vec), file= "SimulationofDisease1withGEN1.csv",  append = TRUE, sep=",",col.names = FALSE)
} ## for loop  Sim.index 




