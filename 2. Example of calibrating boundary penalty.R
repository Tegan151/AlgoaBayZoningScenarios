#######################################################
#This code accompanies the manuscript:
# Scenario planning from the bottom up supports equitable and ecosystem-based approaches to marine spatial planning 
# Tegan Carpenter-Kling*, Hannah Truter, Anne Lemahieu, Bernadette Snow, Mia Strand, Nina Rivers, James Blignaut, Rozanne Bester, Lea Nupnau, Amanda T. Lombard
##Code created by Tegan Carpenter-Kling
##Institute for Coastal and Marine Research, Nelson Mandela University, South Africa
##Please direct queries to tegan.carpenterkling@gmail.com
#########################################################
##This code provides an example of how the manuscript identified the optimal boundary penalty for zoning scenarios using prioritizR
##For more information, please see methods section in manuscript
#########################################################
# This code does not explain prioritizR's inner workings and function. For more detail about prioritizR please see the developer's website: https://prioritizr.net/
###################################

library(prioritizr)
rm(list=ls())

##Step 1. Bring in PrioritizR input data and set up data for the algorithm 
setwd('~/PrioritizR_Input_BoE/')
list.files()
polygon<- readRDS('polygon.rds')

p1 <- readRDS('p1.rds')
p2 <- readRDS('p2.rds')
p3 <- readRDS('p3.rds')
p4 <- readRDS('p4.rds')
p5 <- readRDS('p5.rds')
p6 <- readRDS('p6.rds')

targets <- readRDS('targets.rds')
feature <- gsub('_ConsRes','',names(p1))

zone_features <- zones(names(p1),names(p2),names(p3),names(p4),names(p5),names(p6), zone_names = c("ConsRes", 'CmmNnEx','CmmExtr','CmmFshr','Trnsprt','Mining'),feature_names = feature)

#Set up clumping scheme for boundary penalty (for more info see: ?? add_boundary_penalties)
bm_ply <- boundary_matrix(polygon)
bm_ply <- rescale_matrix(bm_ply)
zm6 <- diag(6)
print(zm6)
colnames(zm6) <- c("ConsRes", 'CmmNnEx','CmmExtr','CmmFshr','Trnsprt','Mining')
rownames(zm6) <- colnames(zm6)

#Select range of boundary penalties to test. 
#Initially, we tested 10-fold increases in the boundary penalty between 0 and 0.1 however these values were inadequate to identify each scenario’s optimal boundary length value and therefore additional iterations with values between 0.13 and 0.01 at 0.015 intervals were run. 
pen <- c(0.1,0.01,0.001,0.0001,0,0.130, 0.115, 0.085, 0.070, 0.055 ,0.040 ,0.025)
label <- c('0-1','0-01','0-001','0-0001','0','0-130', '0-115', '0-085', '0-070', '0-055' ,'0-040' ,'0-025')


fin <- as.data.frame(matrix(ncol=3,nrow=length(pen)))
colnames(fin) <- c('id','cost','boundary')
fin[,'id'] <- pen

#These loops will solve zoning scenarios with each specified boundary penalty, with (loop1) and without (loop2) mandatory allocation of planning units and save the results to your working directory. 
#CAUTION this take a long time to run, to reduce running time, try increasing the gap
setwd('~/CalibrationExampleResults/')

#LOOP1 - problems with mandatory allocation of planning units
for(i in 1:length(pen)){
  p <-
    problem(
      polygon,
      zones(names(p1),names(p2), names(p3),names(p4),names(p5),names(p6),
             zone_names = c("ConsRes",'CmmNnEx','CmmExtr','CmmFshr','Trnsprt','Mining')),
      cost_column = c("ConsRes",'CmmNnEx','CmmExtr','CmmFshr','Trnsprt','Mining') ) %>%
    add_locked_in_constraints(c("lock_in_ConsRes",'lock_in_CmmNnEx','lock_in_CmmExtr','lock_in_CmmFshr','lock_in_Trnsprt','lock_in_Mining')) %>%
    add_min_set_objective()%>%
    add_absolute_targets(targets) %>%
    add_boundary_penalties(penalty = pen[i],  zones=zm6,edge_factor= rep(0.5,6),data=bm_ply) %>%
    add_mandatory_allocation_constraints()%>%
    add_gurobi_solver(numeric_focus = T,threads = parallel::detectCores(TRUE)-1,gap=0.05)
  s <- solve(p)
  cost <- eval_cost_summary(p, s[, c("solution_1_ConsRes","solution_1_CmmNnEx", "solution_1_CmmExtr","solution_1_CmmFshr",'solution_1_Trnsprt','solution_1_Mining')])
  boundary <- eval_boundary_summary(p, s[, c("solution_1_ConsRes","solution_1_CmmNnEx", "solution_1_CmmExtr","solution_1_CmmFshr",'solution_1_Trnsprt','solution_1_Mining')])
  fin$cost[fin$id == pen[i]] <- as.numeric(cost[1,2])
  fin$boundary[fin$id == pen[i]] <- as.numeric(boundary[1,2])
  saveRDS(s,paste0('EmptyBay_Man_',label[i],'.rds'))
  rm(s,p)
}
write.csv(fin, 'EmptyBay_Man_CostBoundary.csv',row.names = F)

#LOOP2 - problems without mandatory allocation of planning units
rm(fin)
fin <- as.data.frame(matrix(ncol=3,nrow=length(pen)))
colnames(fin) <- c('id','cost','boundary')
fin[,'id'] <- pen
for(i in 1:length(pen)){
  p <-
    problem(
      polygon,
      zones(names(p1),names(p2), names(p3),names(p4),names(p5),names(p6),
            zone_names = c("ConsRes",'CmmNnEx','CmmExtr','CmmFshr','Trnsprt','Mining')),
      cost_column = c("ConsRes",'CmmNnEx','CmmExtr','CmmFshr','Trnsprt','Mining') ) %>%
    add_locked_in_constraints(c("lock_in_ConsRes",'lock_in_CmmNnEx','lock_in_CmmExtr','lock_in_CmmFshr','lock_in_Trnsprt','lock_in_Mining')) %>%
    add_min_set_objective()%>%
    add_absolute_targets(targets) %>%
    add_boundary_penalties(penalty = pen[i],  zones=zm6,edge_factor= rep(0.5,6),data=bm_ply) %>%
    # add_mandatory_allocation_constraints()%>%
    add_gurobi_solver(numeric_focus = T,threads = parallel::detectCores(TRUE)-1,gap=0.05)
  s <- solve(p)
  cost <- eval_cost_summary(p, s[, c("solution_1_ConsRes","solution_1_CmmNnEx", "solution_1_CmmExtr","solution_1_CmmFshr",'solution_1_Trnsprt','solution_1_Mining')])
  boundary <- eval_boundary_summary(p, s[, c("solution_1_ConsRes","solution_1_CmmNnEx", "solution_1_CmmExtr","solution_1_CmmFshr",'solution_1_Trnsprt','solution_1_Mining')])
  fin$cost[fin$id == pen[i]] <- as.numeric(cost[1,2])
  fin$boundary[fin$id == pen[i]] <- as.numeric(boundary[1,2])
  saveRDS(s,paste0('EmptyBay_NoMan_',label[i],'.rds'))
  rm(s,p)
}
write.csv(fin, 'EmptyBay_NoMan_CostBoundary.csv',row.names = F)



# Step 2 Identify optimal bounary penalty for each zoning problem######
# The optimal boundary penalty will now be identified as the inflection point in the trade-off curve between the total cost and sum of perimeters of the solutions' zones, following the methods described by Ardron et al. 2010 (https://pacmara.org/wp-content/uploads/2018/05/Marxan-Good-Practices-Handbook-v2-2013.pdf. 
##################
library(ggplot2)
library(topsis)
library(ggpubr)
rm(list=ls())
###Empty bay####
df <- read.csv('EmptyBay_NoMan_CostBoundary.csv')
df$scen <- 'EB_NoMan'
df1 <- read.csv('EmptyBay_Man_CostBoundary.csv')
df1$scen <- 'EB_Man'
df <- rbind(df,df1)
df$id <- as.character(df$id)

id <- unique(df$scen)
fin <- as.data.frame(matrix(ncol=4))
colnames(fin) <- colnames(df)
for(i in 1:length(id)){
  temp <- df[df$scen == id[i],]
  d <- as.matrix(temp[,c(2,3)])
  topsis_results <- topsis(
    decision =d,
    weights = c(1, 1),
    impacts = c("-", "-")
  )
  r <- topsis_results$alt.row[topsis_results$rank == 1]
  fin <- rbind(fin,temp[r,])
}
fin

fin <- fin[!is.na(fin$id),]
a <- fin[fin$scen %in% c(c("EB_Man")),]
b <- fin[fin$scen %in% c("EB_NoMan"),]
a$scen <- factor(a$scen,
                 levels=c("EB_Man"),
                 labels = c('Empty Bay' ))

b$scen <- factor(b$scen,
                 levels=c("EB_NoMan"),
                 labels = c( 'Empty Bay' ))

c <- df[df$scen %in% c(c("EB_Man")),]
d <- df[df$scen %in% c("EB_NoMan"),]
c$scen <- factor(c$scen,
                 levels=c("EB_Man"),
                 labels = c('Empty Bay' ))

d$scen <- factor(d$scen,
                 levels=c("EB_NoMan"),
                 labels = c( 'Empty Bay' ))


aa <- ggplot(c,aes(cost,boundary))+
  geom_point()+
  geom_line()+
  ggtitle('Mandatory allocations set')+
  ylab('Sum of perimeters of the solution zones')+
  xlab('Total cost of solution')+
  geom_point(data=a,aes(cost,boundary),col='red')+
  geom_text(data=a,aes(cost,boundary,label=id),hjust=1,vjust=-0.5,col='red',size=5)+
  facet_wrap(~scen,scale='free',ncol=1)+theme_bw()
aa
bb <- ggplot(d,aes(cost,boundary))+
  geom_point()+
  geom_line()+
  ggtitle('Mandatory allocations not set')+
  ylab('Sum of perimeters of the solution zones')+
  xlab('Total cost of solution')+
  geom_point(data=b,aes(cost,boundary),col='red')+
  geom_text(data=b,aes(cost,boundary,label=id),hjust=1,vjust=-0.5,col='red',size=5)+
  facet_wrap(~scen,scale='free',ncol=1)+theme_bw()
bb
fig <- ggarrange(aa, bb, 
                 ncol = 2)


fig
