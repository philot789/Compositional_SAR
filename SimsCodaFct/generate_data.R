## date: 10.06.2024
## description: generate compositional data with 6 components for given locations 
##              for t=15 time points without dependency and with temporal, 
##              spatial and temporal-spatial depencies


# load RData Dats
load("")#insert path
dta<-data.frame(x=rep(Dats[,1],each=15),y=rep(Dats[,2],each=15),t=rep(c(1:15),length(Dats[,1])))

library(dplyr)
library(ggplot2)
library(manipulate)

# no temporal or spatial dependence:
set.seed(101)
## sample from uniform distr
comps<-matrix(runif(nrow(Dats)*length(t)*6),nrow=nrow(Dats)*length(t), ncol=6)
comps<-comps*100/rowSums(comps)
comps<-data.frame(comps)
colnames(comps)<-c("c1","c2","c3","c4","c5","c6")
dta_undep<-cbind(dta,comps)
## save data
saveRDS(dta_undep,"dta_undependent.rds")


# temporal dependence:
set.seed(101)
dta_temp<-data.frame(x=rep(Dats[,1],each=15),y=rep(Dats[,2],each=15),pd=rep(Dats[,3],each=15),t=rep(c(1:15),length(Dats[,1])),c1=NA,c2=NA,c3=NA,c4=NA,c5=NA,c6=NA)
dta_temp<-dta_temp[order(dta_temp$pd),]
ind<-c(1:nrow(Dats))
t<-c(1:15)
## define temporal structure
temporal_structure<- sin(seq(0,  2*pi, length.out = 15))
plot(t,temporal_structure)
temporal_structure<-diff(temporal_structure)
## plot temporal structure
for (i in ind){
  times_ind<-c((1+(i-1)*length(t)):(length(t)+(i-1)*length(t)))
  start_comp<-runif(6,min = 0, max=100)
  ## for an observation with pd=mean(pd)
  start_comp<-start_comp*100/sum(start_comp)
  max_distance<-max(cumsum(temporal_structure[temporal_structure>0]))
  min_distance<- min(cumsum(temporal_structure[temporal_structure<0]))
  factor_min<-abs((-start_comp[6])/min_distance)
  factor_max<-abs((100-start_comp[6])/max_distance)
  factor<-min(factor_min,factor_max)
  dta_temp[times_ind[1],c(5:10)]<-start_comp
  for (tp in times_ind[-1]){
    time<-(tp-1)%%length(t)
    dta_temp[tp,c(5:10)]<-dta_temp[tp-1,c(5:10)]
    dta_temp[tp,10]<-dta_temp[tp-1,10]+temporal_structure[time]*factor
    dta_temp[tp,c(5:10)]<-dta_temp[tp,c(5:10)]*100/sum(dta_temp[tp,c(5:10)])
  }
}
## save data 
saveRDS(dta_temp,"data_temporal.rds")
## plot data for different locations in dependence of time
manipulate(
  plot((dta_temp%>%filter(x==unique(dta_temp$x)[i])%>%select(t))[,1],
       (dta_temp%>%filter(x==unique(dta_temp$x)[i])%>%select(c6))[,1],
       ylab="C6",xlab="t", main=paste0("x=",unique(dta_temp$x)[i])), 
  i = slider( min = 1, max = length(unique(dta_temp$pd)) ) ) 


# spatial dependence:
set.seed(101)
dta_spatial<-data.frame(x=rep(Dats[,1],each=15),y=rep(Dats[,2],each=15),pd=rep(Dats[,3],each=15),t=rep(c(1:15),length(Dats[,1])),c1=NA,c2=NA,c3=NA,c4=NA,c5=NA,c6=NA)
dta_spatial<-dta_spatial[order(dta_spatial$pd),]
dta_spatial<-dta_spatial[order(dta_spatial$t),]
## plot used spatial structure (for t=1)
ggplot(data=dta_spatial%>%filter(t==1),aes(x=x,y=y,color=pd-mean(pd)))+
  geom_point()
t<-c(1:15)
for (tp in t){
  ind<-c((1+(tp-1)*nrow(Dats)):(nrow(Dats)+(tp-1)*nrow(Dats)))
  start_comp<-runif(6,min = 0, max=100)
  ## for an observation with pd=mean(pd)
  start_comp<-start_comp*100/sum(start_comp)
  max_distance<-max(dta_spatial$pd-mean(dta_spatial$pd))
  min_distance<- min(dta_spatial$pd-mean(dta_spatial$pd))
  factor_min<-(-start_comp[6])/min_distance
  factor_max<-(100-start_comp[6])/max_distance
  factor<-min(factor_min,factor_max)
  for (i in ind){
    dta_spatial[i,c(5:10)]<-start_comp
    dta_spatial[i,10]<-dta_spatial[i,10]+0.5*factor*(dta_spatial[i,3]-mean(dta_spatial$pd))
    dta_spatial[i,c(5:10)]<-dta_spatial[i,c(5:10)]*100/sum(dta_spatial[i,c(5:10)])
  }
}
## save data 
saveRDS(dta_spatial,"data_spatial.rds")
## plot data for different time points in dependence of pd
manipulate(
  plot((dta_spatial%>%filter(t==unique(dta_spatial$t)[i])%>%select(pd))[,1],
       (dta_spatial%>%filter(t==unique(dta_spatial$t)[i])%>%select(c6))[,1],
       ylab="C6",xlab="pd", main=paste0("t=",unique(dta_spatial$t)[i])), 
  i = slider( min = 1, max = length(unique(dta_spatial$t)) ) ) 


# temporal and spatial
set.seed(101)
dta_spatial_temp<-data.frame(x=rep(Dats[,1],each=15),y=rep(Dats[,2],each=15),pd=rep(Dats[,3],each=15),t=rep(c(1:15),length(Dats[,1])),c1=NA,c2=NA,c3=NA,c4=NA,c5=NA,c6=NA)
dta_spatial_temp<-dta_spatial_temp[order(dta_spatial_temp$t),]
dta_spatial_temp<-dta_spatial_temp[order(dta_spatial_temp$pd),]
t<-c(1:15)
ind<-c(1:nrow(Dats))
pds<-unique(dta_spatial_temp$pd)
temp_spatial_structure<-t((mean(pds)-pds)%*%t(mean(t)-t))
## plot temporal-spatial structure
structure_dat<-data.frame(pd=dta_spatial_temp$pd,t=dta_spatial_temp$t, diff=as.vector(temp_spatial_structure))
ggplot(data=structure_dat,aes(x=pd,y=t,z=diff))+
  geom_contour_filled()
## sample random start composition
start_comp<-runif(6,min = 0, max=100)
start_comp<-start_comp*100/sum(start_comp)
start_comp[6]
max(temp_spatial_structure)
min(temp_spatial_structure)
factor_min<-abs((-start_comp[6])/min(temp_spatial_structure))
factor_max<-abs((100-start_comp[6])/max(temp_spatial_structure))
factor<-min(factor_min,factor_max)
for (i in ind){
  times_ind<-c((1+(i-1)*length(t)):(length(t)+(i-1)*length(t)))
  for (tp in t){
    dta_spatial_temp[ times_ind[tp],c(5:10)]<-start_comp
    dta_spatial_temp[ times_ind[tp],10]<-dta_spatial_temp[times_ind[tp],10]+temp_spatial_structure[tp,i]*factor
    dta_spatial_temp[times_ind[tp],c(5:10)]<-dta_spatial_temp[times_ind[tp],c(5:10)]*100/sum(dta_spatial_temp[times_ind[tp],c(5:10)])
  }
}
## save data 
saveRDS(dta_spatial_temp,"data_spatial_temp.rds")
## plot data for different time points in dependence of pd
manipulate(
  plot((dta_spatial_temp%>%filter(t==unique(dta_spatial_temp$t)[i])%>%select(pd))[,1],
       (dta_spatial_temp%>%filter(t==unique(dta_spatial_temp$t)[i])%>%select(c6))[,1],
       ylab="C6",xlab="pd", main=paste0("t=",unique(dta_spatial_temp$t)[i])), 
  i = slider( min = 1, max = length(unique(dta_spatial_temp$t)) ) ) 
## plot data for different locations in dependence of time
manipulate(
  plot((dta_spatial_temp%>%filter(pd==unique(dta_spatial_temp$pd)[i])%>%select(t))[,1],
       (dta_spatial_temp%>%filter(pd==unique(dta_spatial_temp$pd)[i])%>%select(c6))[,1],
       ylab="C6",xlab="t", main=paste0("pd=",unique(dta_spatial_temp$pd)[i])), 
  i = slider( min = 1, max = length(unique(dta_spatial_temp$pd)) ) ) 



