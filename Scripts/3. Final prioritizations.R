#######################################################
#This code accompanies the manuscript:
# Scenario planning from the bottom up supports equitable and ecosystem-based approaches to marine spatial planning 
# Tegan Carpenter-Kling*, Hannah Truter, Anne Lemahieu, Bernadette Snow, Mia Strand, Nina Rivers, James Blignaut, Rozanne Bester, Lea Nupnau, Amanda T. Lombard
##Code created by Tegan Carpenter-Kling
##Institute for Coastal and Marine Research, Nelson Mandela University, South Africa
##Please direct queries to tegan.carpenterkling@gmail.com
#########################################################
##This code provides the final code of the zoning scenarios prioiritizations and mapping presented in the manuscript 
##The data necessary to run this code can be found on github
##For more information, please see methods section in manuscript
#########################################################
# This code does not explain prioritizR's inner workings and function. For more detail about prioritizR please see the developer's website: https://prioritizr.net/
###################################

library(sf)
library(prioritizr)
library(ggplot2)
rm(list=ls())

#Bay of plenty######
setwd('~')

polygon<- readRDS('PrioritizR_Input_BoP/polygon.rds')

p1 <- readRDS('PrioritizR_Input_BoP/p1.rds')
p2 <- readRDS('PrioritizR_Input_BoP/p2.rds')
p3 <- readRDS('PrioritizR_Input_BoP/p3.rds')
p4 <- readRDS('PrioritizR_Input_BoP/p4.rds')
p5 <- readRDS('PrioritizR_Input_BoP/p5.rds')


targets <- readRDS('PrioritizR_Input_BoP/targets.rds')
feature <- gsub('_ConsRes','',names(p1))

zone_features <- zones(names(p1),names(p2),names(p3),names(p4),names(p5), zone_names = c("ConsRes", 'CmmNnEx','CmmExtr','CmmFshr','Trnsprt'),feature_names = feature)


#Set up penalty factor
bm_ply <- boundary_matrix(polygon)
bm_ply <- rescale_matrix(bm_ply)
zm6 <- diag(5)
colnames(zm6) <- c("ConsRes", 'CmmNnEx','CmmExtr','CmmFshr','Trnsprt')
rownames(zm6) <- colnames(zm6)


p <-
    problem(
      polygon,
      zones(names(p1),names(p2), names(p3),names(p4),names(p5),
            zone_names = c("ConsRes", 'CmmNnEx','CmmExtr','CmmFshr','Trnsprt')),
      cost_column = c("ConsRes", 'CmmNnEx','CmmExtr','CmmFshr','Trnsprt') ) %>%
    add_locked_in_constraints(c("lock_in_ConsRes",'lock_in_CmmNnEx','lock_in_CmmExtr','lock_in_CmmFshr','lock_in_Trnsprt')) %>%
    add_min_set_objective()%>%
    add_absolute_targets(targets) %>%
    add_boundary_penalties(penalty = 0.055,  zones=zm6,edge_factor= rep(0.5,5),data=bm_ply) %>%
    add_mandatory_allocation_constraints()%>%
    add_gurobi_solver(numeric_focus = T,threads = parallel::detectCores(TRUE)-1,gap=0.05)
  s <- solve(p)
  
  saveRDS(s,'Final_prioritizations/BayOfPlenty_Man.rds')
  
  p <-
    problem(
      polygon,
      zones(names(p1),names(p2), names(p3),names(p4),names(p5),
            zone_names = c("ConsRes", 'CmmNnEx','CmmExtr','CmmFshr','Trnsprt')),
      cost_column = c("ConsRes", 'CmmNnEx','CmmExtr','CmmFshr','Trnsprt') ) %>%
    add_locked_in_constraints(c("lock_in_ConsRes",'lock_in_CmmNnEx','lock_in_CmmExtr','lock_in_CmmFshr','lock_in_Trnsprt')) %>%
    add_min_set_objective()%>%
    add_absolute_targets(targets) %>%
    add_boundary_penalties(penalty = 0.1,  zones=zm6,edge_factor= rep(0.5,5),data=bm_ply) %>%
    # add_mandatory_allocation_constraints()%>%
    add_gurobi_solver(numeric_focus = T,threads = parallel::detectCores(TRUE)-1,gap=0.05)
  s <- solve(p)
  
  saveRDS(s,'Final_prioritizations/BayOfPlenty_NoMan.rds')
  
  
  
#Business as usual#####  
  rm(list=ls())
  
  polygon<- readRDS('PrioritizR_Input_BAU/polygon.rds')
  
  p1 <- readRDS('PrioritizR_Input_BAU/p1.rds')
  p2 <- readRDS('PrioritizR_Input_BAU/p2.rds')
  p3 <- readRDS('PrioritizR_Input_BAU/p3.rds')
  p4 <- readRDS('PrioritizR_Input_BAU/p4.rds')
  p5 <- readRDS('PrioritizR_Input_BAU/p5.rds')
  p6 <- readRDS('PrioritizR_Input_BAU/p6.rds')
  
  
  targets <- readRDS('PrioritizR_Input_BAU/targets.rds')
  feature <- gsub('_ConsRes','',names(p1))
  
  zone_features <- zones(names(p1),names(p2),names(p3),names(p4),names(p5),names(p6), zone_names = c("ConsRes", 'CmmNnEx','CmmExtr','CmmFshr','Trnsprt','Mining'),feature_names = feature)
  
  #Set up penalty factor
  bm_ply <- boundary_matrix(polygon)
  bm_ply <- rescale_matrix(bm_ply)
  zm6 <- diag(6)
  print(zm6)
  colnames(zm6) <- c("ConsRes", 'CmmNnEx','CmmExtr','CmmFshr','Trnsprt','Mining')
  rownames(zm6) <- colnames(zm6)
  
p <-
    problem(
      polygon,
      zones(names(p1),names(p2), names(p3),names(p4),names(p5),names(p6),
            zone_names = c("ConsRes", 'CmmNnEx','CmmExtr','CmmFshr','Trnsprt','Mining')),
      cost_column = c("ConsRes", 'CmmNnEx','CmmExtr','CmmFshr','Trnsprt','Mining') ) %>%
    add_locked_in_constraints(c("lock_in_ConsRes",'lock_in_CmmNnEx','lock_in_CmmExtr','lock_in_CmmFshr','lock_in_Trnsprt','lock_in_Mining')) %>%
    add_min_set_objective()%>%
    add_absolute_targets(targets) %>%
    add_boundary_penalties(penalty = 0.025,  zones=zm6,edge_factor= rep(0.5,6),data=bm_ply) %>%
    add_mandatory_allocation_constraints()%>%
    add_gurobi_solver(numeric_focus = T,threads = parallel::detectCores(TRUE)-1,gap=0.05)
  s <- solve(p)
  
  saveRDS(s,'Final_prioritizations/BAU_Man.rds')
  
  p <-
      problem(
        polygon,
        zones(names(p1),names(p2), names(p3),names(p4),names(p5),names(p6),
              zone_names = c("ConsRes", 'CmmNnEx','CmmExtr','CmmFshr','Trnsprt','Mining')),
        cost_column = c("ConsRes", 'CmmNnEx','CmmExtr','CmmFshr','Trnsprt','Mining') ) %>%
        add_locked_in_constraints(c("lock_in_ConsRes",'lock_in_CmmNnEx','lock_in_CmmExtr','lock_in_CmmFshr','lock_in_Trnsprt','lock_in_Mining')) %>%
    add_min_set_objective()%>%
    add_absolute_targets(targets) %>%
    add_boundary_penalties(penalty = 0.025,  zones=zm6,edge_factor= rep(0.5,6),data=bm_ply) %>%
    # add_mandatory_allocation_constraints()%>%
    add_gurobi_solver(numeric_focus = T,threads = parallel::detectCores(TRUE)-1,gap=0.05)
  s <- solve(p)
  
  saveRDS(s,'Final_prioritizations/BAU_NoMan.rds')
  
  

  
#Empty bay###### 
  rm(list=ls())

  polygon<- readRDS('PrioritizR_Input_BoE/polygon.rds')
  
  p1 <- readRDS('PrioritizR_Input_BoE/p1.rds')
  p2 <- readRDS('PrioritizR_Input_BoE/p2.rds')
  p3 <- readRDS('PrioritizR_Input_BoE/p3.rds')
  p4 <- readRDS('PrioritizR_Input_BoE/p4.rds')
  p5 <- readRDS('PrioritizR_Input_BoE/p5.rds')
  p6 <- readRDS('PrioritizR_Input_BoE/p6.rds')
  
  targets <- readRDS('PrioritizR_Input_BoE/targets.rds')
  feature <- gsub('_ConsRes','',names(p1))
  
  zone_features <- zones(names(p1),names(p2),names(p3),names(p4),names(p5),names(p6), zone_names = c("ConsRes", 'CmmNnEx','CmmExtr','CmmFshr','Trnsprt','Mining'),feature_names = feature)
  
  
  #Set up penalty factor
  bm_ply <- boundary_matrix(polygon)
  bm_ply <- rescale_matrix(bm_ply)
  zm6 <- diag(6)
  print(zm6)
  colnames(zm6) <- c("ConsRes", 'CmmNnEx','CmmExtr','CmmFshr','Trnsprt','Mining')
  rownames(zm6) <- colnames(zm6)
  
  
  p <-
    problem(
      polygon,
      zones(names(p1),names(p2), names(p3),names(p4),names(p5),names(p6),
            zone_names = c("ConsRes", 'CmmNnEx','CmmExtr','CmmFshr','Trnsprt','Mining')),
      cost_column = c("ConsRes", 'CmmNnEx','CmmExtr','CmmFshr','Trnsprt','Mining') ) %>%
    add_locked_in_constraints(c("lock_in_ConsRes",'lock_in_CmmNnEx','lock_in_CmmExtr','lock_in_CmmFshr','lock_in_Trnsprt','lock_in_Mining')) %>%
    add_min_set_objective()%>%
    add_absolute_targets(targets) %>%
    add_boundary_penalties(penalty = 0.13,  zones=zm6,edge_factor= rep(0.5,6),data=bm_ply) %>%
    add_mandatory_allocation_constraints()%>%
    add_gurobi_solver(numeric_focus = T,threads = parallel::detectCores(TRUE)-1,gap=0.05)
  s <- solve(p)
  
  saveRDS(s,'Final_prioritizations/EmptyBay_Man.rds')
  
  p <-
    problem(
      polygon,
      zones(names(p1),names(p2), names(p3),names(p4),names(p5),names(p6),
            zone_names = c("ConsRes", 'CmmNnEx','CmmExtr','CmmFshr','Trnsprt','Mining')),
      cost_column = c("ConsRes", 'CmmNnEx','CmmExtr','CmmFshr','Trnsprt','Mining') ) %>%
    add_locked_in_constraints(c("lock_in_ConsRes",'lock_in_CmmNnEx','lock_in_CmmExtr','lock_in_CmmFshr','lock_in_Trnsprt','lock_in_Mining')) %>%
    add_min_set_objective()%>%
    add_absolute_targets(targets) %>%
    add_boundary_penalties(penalty = 0.04,  zones=zm6,edge_factor= rep(0.5,6),data=bm_ply) %>%
    # add_mandatory_allocation_constraints()%>%
    add_gurobi_solver(numeric_focus = T,threads = parallel::detectCores(TRUE)-1,gap=0.05)
  s <- solve(p)
  
  saveRDS(s,'Final_prioritizations/EmptyBay_NoMan.rds')
  
  
 
  
#Untouched paradise#############
  rm(list=ls())
  
  polygon<- readRDS('PrioritizR_Input_UTP/polygon.rds')
  
  p1 <- readRDS('PrioritizR_Input_UTP/p1.rds')
  p2 <- readRDS('PrioritizR_Input_UTP/p2.rds')
  p3 <- readRDS('PrioritizR_Input_UTP/p3.rds')
  p4 <- readRDS('PrioritizR_Input_UTP/p4.rds')
  p5 <- readRDS('PrioritizR_Input_UTP/p5.rds')
  
  targets <- readRDS('PrioritizR_Input_UTP/targets.rds')
  feature <- gsub('_ConsRes','',names(p1))
  
  zone_features <- zones(names(p1),names(p2),names(p3),names(p4),names(p5), zone_names = c("ConsRes", 'CmmNnEx','CmmExtr','CmmFshr','Trnsprt'),feature_names = feature)
  
  #Set up penalty factor
  bm_ply <- boundary_matrix(polygon)
  bm_ply <- rescale_matrix(bm_ply)
  zm6 <- diag(5)
  print(zm6)
  colnames(zm6) <- c("ConsRes", 'CmmNnEx','CmmExtr','CmmFshr','Trnsprt')
  rownames(zm6) <- colnames(zm6)
  
p <-
    problem(
      polygon,
      zones(names(p1),names(p2), names(p3),names(p4),names(p5),
            zone_names = c("ConsRes", 'CmmNnEx','CmmExtr','CmmFshr','Trnsprt')),
      cost_column = c("ConsRes", 'CmmNnEx','CmmExtr','CmmFshr','Trnsprt') ) %>%
    add_locked_in_constraints(c("lock_in_ConsRes",'lock_in_CmmNnEx','lock_in_CmmExtr','lock_in_CmmFshr','lock_in_Trnsprt')) %>%
    add_min_set_objective()%>%
    add_absolute_targets(targets) %>%
    add_boundary_penalties(penalty = 0.07,  zones=zm6,edge_factor= rep(0.5,5),data=bm_ply) %>%
    add_mandatory_allocation_constraints()%>%
    add_gurobi_solver(numeric_focus = T,threads = parallel::detectCores(TRUE)-1,gap=0.05)
  s <- solve(p)
  
  saveRDS(s,'Final_prioritizations/Untouched_Man.rds')
  
  p <-
    problem(
      polygon,
      zones(names(p1),names(p2), names(p3),names(p4),names(p5),
            zone_names = c("ConsRes", 'CmmNnEx','CmmExtr','CmmFshr','Trnsprt')),
      cost_column = c("ConsRes", 'CmmNnEx','CmmExtr','CmmFshr','Trnsprt') ) %>%
    add_locked_in_constraints(c("lock_in_ConsRes",'lock_in_CmmNnEx','lock_in_CmmExtr','lock_in_CmmFshr','lock_in_Trnsprt')) %>%
    add_min_set_objective()%>%
    add_absolute_targets(targets) %>%
    add_boundary_penalties(penalty = 0.07,  zones=zm6,edge_factor= rep(0.5,5),data=bm_ply) %>%
    # add_mandatory_allocation_constraints()%>%
    add_gurobi_solver(numeric_focus = T,threads = parallel::detectCores(TRUE)-1,gap=0.05)
  s <- solve(p)
  
  saveRDS(s,'Final_prioritizations/Untouched_NoMan.rds')
  
#Plot prioritization######
  library(ggplot2)
  rm(list=ls())
  setwd('~')
  #Bay of plenty - mandatory####
  ras <- readRDS('Final_prioritizations/BayOfPlenty_Man.rds')
  ras$ANCHOR <- 0
  ras$ANCHOR[ras$ANCHOR_ConsRes>0] <- 1
  ras$MARICUL <- 0
  ras$MARICUL[ras$MARICUL_ConsRes>0] <- 1
  
  shp <- read_sf('SetUpExampleData/MPAs.shp')
  ras$MPA <- shp$MPA1
  ras$concon <- 0
  ras$conres <- 0
  ras$concon[ras$MPA == "Controlled" ] <- 1
  ras$conres[ras$MPA %in% c("Restricted", "No-take") ] <- 1
  ras$solution_1_ConsRes[ras$concon == 1] <- 0
  ras$solution_1_ConsRes[ras$conres == 1] <- 0
  ras$solution <- category_vector(
    ras[, c("solution_1_ConsRes","conres",'concon',"solution_1_CmmNnEx", "solution_1_CmmExtr","solution_1_CmmFshr",'ANCHOR','MARICUL')]
  )
  
  unique(ras$solution)
  r <- ras[, c("solution")]
  # plot solution
  r$solution[r$solution == 0] <- 'Not allocated'
  r$solution[r$solution == 1] <- 'solution_1_ConsRes'
  r$solution[r$solution == 2] <- 'concon'
  r$solution[r$solution == 3] <- 'conres'
  r$solution[r$solution == 4] <- 'solution_1_CmmNnEx'
  r$solution[r$solution == 5] <- 'solution_1_CmmExtr'
  r$solution[r$solution == 6] <- 'solution_1_CmmFshr'
  r$solution[r$solution == 7] <- 'ANCHOR'
  r$solution[r$solution == 8] <- 'MARICUL'
  # r$solution[r$solution == 8] <- 'mining'
  r$solution <- factor(r$solution,levels= c("solution_1_ConsRes",'conres','concon',"solution_1_CmmNnEx", "solution_1_CmmExtr","solution_1_CmmFshr",'ANCHOR','MARICUL'),labels=c('Bio',"MPARestricted","MPAControlled","CommunityNonExtractive" , "CommunityExtractive","Commercialfisheries" ,"Anchorage","Mariculture"))
  
  land <- read_sf('SetUpExampleData/coast.shp')
  
  land <- st_transform(land, crs = st_crs(ras))
  
  ggplot()+
    geom_sf(data=st_union(r[r$solution == 'Bio',]),fill='#42C713',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'MPARestricted',]),fill='#D2F3AE',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'MPAControlled',]),fill='#BBD999',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'CommunityNonExtractive',]),fill='#3249DC',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'CommunityExtractive',]),fill='#4BCAC2',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'Commercialfisheries',]),fill='#8F4DEF',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'Anchorage',]),fill='#F4D03F',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'Mariculture',]),fill='#E88E3E',col='transparent')+
    geom_sf(data=r,fill='transparent',col='grey50',linewidth=0.05)+
    geom_sf(data=land,col='black',fill='grey90')+
    coord_sf(xlim = c(46137 ,158130 ),ylim=c(-3869394 ,-3810394 ))+
    
    theme_bw() +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          strip.background=element_blank()) +
    xlab('Longitude')+
    ylab('Latitude')
  ggsave('Final_scenario_maps/BoP_Man.tiff',dpi=300)
  
  
  #Bay of plenty - No mandatory####
  rm(list=ls())
  
  ras <- readRDS('Final_prioritizations/BayOfPlenty_NoMan.rds')
  ras$ANCHOR <- 0
  ras$ANCHOR[ras$ANCHOR_ConsRes>0] <- 1
  ras$MARICUL <- 0
  ras$MARICUL[ras$MARICUL_ConsRes>0] <- 1
  
  shp <- read_sf('SetUpExampleData/MPAs.shp')
  ras$MPA <- shp$MPA1
  ras$concon <- 0
  ras$conres <- 0
  ras$concon[ras$MPA == "Controlled" ] <- 1
  ras$conres[ras$MPA %in% c("Restricted", "No-take") ] <- 1
  ras$solution_1_ConsRes[ras$concon == 1] <- 0
  ras$solution_1_ConsRes[ras$conres == 1] <- 0
  
  ras$solution <- category_vector(
    ras[, c("solution_1_ConsRes","conres",'concon',"solution_1_CmmNnEx", "solution_1_CmmExtr","solution_1_CmmFshr",'ANCHOR','MARICUL')]
  )
  
  r <- ras[, c("solution")]
  # plot solution
  r$solution[r$solution == 0] <- 'Not allocated'
  r$solution[r$solution == 1] <- 'solution_1_ConsRes'
  r$solution[r$solution == 2] <- 'concon'
  r$solution[r$solution == 3] <- 'conres'
  r$solution[r$solution == 4] <- 'solution_1_CmmNnEx'
  r$solution[r$solution == 5] <- 'solution_1_CmmExtr'
  r$solution[r$solution == 6] <- 'solution_1_CmmFshr'
  r$solution[r$solution == 7] <- 'ANCHOR'
  r$solution[r$solution == 8] <- 'MARICUL'
  # r$solution[r$solution == 8] <- 'mining'
  r$solution <- factor(r$solution,levels= c("solution_1_ConsRes",'conres','concon',"solution_1_CmmNnEx", "solution_1_CmmExtr","solution_1_CmmFshr",'ANCHOR','MARICUL'),labels=c('Bio',"MPARestricted","MPAControlled","CommunityNonExtractive" , "CommunityExtractive","Commercialfisheries" ,"Anchorage","Mariculture"))
  
  land <- read_sf('SetUpExampleData/coast.shp')
  land <- st_transform(land, crs = st_crs(ras))
  
  ggplot()+
    geom_sf(data=st_union(r[r$solution == 'Bio',]),fill='#42C713',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'MPARestricted',]),fill='#D2F3AE',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'MPAControlled',]),fill='#BBD999',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'CommunityNonExtractive',]),fill='#3249DC',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'CommunityExtractive',]),fill='#4BCAC2',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'Commercialfisheries',]),fill='#8F4DEF',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'Anchorage',]),fill='#F4D03F',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'Mariculture',]),fill='#E88E3E',col='transparent')+
    geom_sf(data=r,fill='transparent',col='grey50',linewidth=0.05)+
    geom_sf(data=land,col='black',fill='grey90')+
    coord_sf(xlim = c(46137 ,158130 ),ylim=c(-3869394 ,-3810394 ))+
    
    theme_bw() +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          strip.background=element_blank()) +
    xlab('Longitude')+
    ylab('Latitude')
  ggsave('Final_scenario_maps/BoP_NoMan.tiff',dpi=300)
  
  #Business as usual - No mandatory####
  rm(list=ls())
  
  ras <- readRDS('Final_prioritizations/BAU_NoMan.rds')
  ras$ANCHOR <- 0
  ras$ANCHOR[ras$ANCHOR_ConsRes>0] <- 1
  ras$MARICUL <- 0
  ras$MARICUL[ras$MARICUL_ConsRes>0] <- 1
  
  shp <- read_sf('SetUpExampleData/MPAs.shp')
  ras$MPA <- shp$MPA1
  ras$concon <- 0
  ras$conres <- 0
  ras$concon[ras$MPA == "Controlled" ] <- 1
  ras$conres[ras$MPA %in% c("Restricted", "No-take") ] <- 1
  ras$solution_1_ConsRes[ras$concon == 1] <- 0
  ras$solution_1_ConsRes[ras$conres == 1] <- 0
  
  ras$solution <- category_vector(
    ras[, c("solution_1_ConsRes","conres",'concon',"solution_1_CmmNnEx", "solution_1_CmmExtr","solution_1_CmmFshr",'ANCHOR','MARICUL',"solution_1_Mining")]
  )
  
  unique(ras$solution)
  r <- ras[, c("solution")]
  # plot solution
  r$solution[r$solution == 0] <- 'Not allocated'
  r$solution[r$solution == 1] <- 'solution_1_ConsRes'
  r$solution[r$solution == 2] <- 'concon'
  r$solution[r$solution == 3] <- 'conres'
  r$solution[r$solution == 4] <- 'solution_1_CmmNnEx'
  r$solution[r$solution == 5] <- 'solution_1_CmmExtr'
  r$solution[r$solution == 6] <- 'solution_1_CmmFshr'
  r$solution[r$solution == 7] <- 'ANCHOR'
  r$solution[r$solution == 8] <- 'MARICUL'
  r$solution[r$solution == 9] <- 'mining'
  r$solution <- factor(r$solution,levels= c("solution_1_ConsRes",'conres','concon',"solution_1_CmmNnEx", "solution_1_CmmExtr","solution_1_CmmFshr",'ANCHOR','MARICUL',"mining"),labels=c('Bio',"MPARestricted","MPAControlled","CommunityNonExtractive" , "CommunityExtractive","Commercialfisheries" ,"Anchorage","Mariculture",'Mining'))
  
  land <- read_sf('SetUpExampleData/coast.shp')
  land <- st_transform(land, crs = st_crs(ras))
  
  
  
  
  
  ggplot()+
    geom_sf(data=st_union(r[r$solution == 'Bio',]),fill='#42C713',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'MPARestricted',]),fill='#D2F3AE',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'MPAControlled',]),fill='#BBD999',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'CommunityNonExtractive',]),fill='#3249DC',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'CommunityExtractive',]),fill='#4BCAC2',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'Commercialfisheries',]),fill='#8F4DEF',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'Anchorage',]),fill='#F4D03F',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'Mariculture',]),fill='#E88E3E',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'Mining',]),fill='black',col='transparent')+
    geom_sf(data=r,fill='transparent',col='grey50',linewidth=0.05)+
    geom_sf(data=land,col='black',fill='grey90')+
    coord_sf(xlim = c(46137 ,158130 ),ylim=c(-3869394 ,-3810394 ))+
    
    theme_bw() +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          strip.background=element_blank()) +
    xlab('Longitude')+
    ylab('Latitude')
  ggsave('Final_scenario_maps/BAU_NoMan.tiff',dpi=300)
  
  #Business as usual - mandatory####
  rm(list=ls())
  ras <- readRDS('Final_prioritizations/BAU_Man.rds')
  ras$ANCHOR <- 0
  ras$ANCHOR[ras$ANCHOR_ConsRes>0] <- 1
  ras$MARICUL <- 0
  ras$MARICUL[ras$MARICUL_ConsRes>0] <- 1
  
  shp <- read_sf('SetUpExampleData/MPAs.shp')
  ras$MPA <- shp$MPA1
  ras$concon <- 0
  ras$conres <- 0
  ras$concon[ras$MPA == "Controlled" ] <- 1
  ras$conres[ras$MPA %in% c("Restricted", "No-take") ] <- 1
  ras$solution_1_ConsRes[ras$concon == 1] <- 0
  ras$solution_1_ConsRes[ras$conres == 1] <- 0
  ras$solution <- category_vector(
    ras[, c("solution_1_ConsRes","conres",'concon',"solution_1_CmmNnEx", "solution_1_CmmExtr","solution_1_CmmFshr",'ANCHOR','MARICUL',"solution_1_Mining")]
  )
  
  unique(ras$solution)
  r <- ras[, c("solution")]
  # plot solution
  r$solution[r$solution == 0] <- 'Not allocated'
  r$solution[r$solution == 1] <- 'solution_1_ConsRes'
  r$solution[r$solution == 2] <- 'concon'
  r$solution[r$solution == 3] <- 'conres'
  r$solution[r$solution == 4] <- 'solution_1_CmmNnEx'
  r$solution[r$solution == 5] <- 'solution_1_CmmExtr'
  r$solution[r$solution == 6] <- 'solution_1_CmmFshr'
  r$solution[r$solution == 7] <- 'ANCHOR'
  r$solution[r$solution == 8] <- 'MARICUL'
  r$solution[r$solution == 9] <- 'mining'
  r$solution <- factor(r$solution,levels= c("solution_1_ConsRes",'conres','concon',"solution_1_CmmNnEx", "solution_1_CmmExtr","solution_1_CmmFshr",'ANCHOR','MARICUL',"mining"),labels=c('Bio',"MPARestricted","MPAControlled","CommunityNonExtractive" , "CommunityExtractive","Commercialfisheries" ,"Anchorage","Mariculture",'Mining'))
  
  land <- read_sf('SetUpExampleData/coast.shp')
  land <- st_transform(land, crs = st_crs(ras))
ggplot()+
    geom_sf(data=st_union(r[r$solution == 'Bio',]),fill='#42C713',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'MPARestricted',]),fill='#D2F3AE',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'MPAControlled',]),fill='#BBD999',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'CommunityNonExtractive',]),fill='#3249DC',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'CommunityExtractive',]),fill='#4BCAC2',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'Commercialfisheries',]),fill='#8F4DEF',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'Anchorage',]),fill='#F4D03F',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'Mariculture',]),fill='#E88E3E',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'Mining',]),fill='black',col='transparent')+
    geom_sf(data=r,fill='transparent',col='grey50',linewidth=0.05)+
    geom_sf(data=land,col='black',fill='grey90')+
    coord_sf(xlim = c(46137 ,158130 ),ylim=c(-3869394 ,-3810394 ))+
    
    theme_bw() +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          strip.background=element_blank()) +
    xlab('Longitude')+
    ylab('Latitude')
  ggsave('Final_scenario_maps/BAU_Man.tiff',dpi=300)
  
  #Empty - mandatory####
  rm(list=ls())
  
  ras <- readRDS('Final_prioritizations/EmptyBay_Man.rds')
  ras$ANCHOR <- 0
  ras$ANCHOR[ras$ANCHOR_ConsRes>0] <- 1
  ras$MARICUL <- 0
  ras$MARICUL[ras$MARICUL_ConsRes>0] <- 1
  
  ras$solution <- category_vector(
    ras[, c("solution_1_ConsRes","solution_1_CmmNnEx", "solution_1_CmmExtr","solution_1_CmmFshr",'ANCHOR','MARICUL',"solution_1_Mining")]
  )
  
  unique(ras$solution)
  r <- ras[, c("solution")]
  # plot solution
  r$solution[r$solution == 0] <- 'Not allocated'
  r$solution[r$solution == 1] <- 'solution_1_ConsRes'
  r$solution[r$solution == 2] <- 'solution_1_CmmNnEx'
  r$solution[r$solution == 3] <- 'solution_1_CmmExtr'
  r$solution[r$solution == 4] <- 'solution_1_CmmFshr'
  r$solution[r$solution == 5] <- 'ANCHOR'
  r$solution[r$solution == 6] <- 'MARICUL'
  r$solution[r$solution == 7] <- 'mining'
  r$solution <- factor(r$solution,levels= c("solution_1_ConsRes","solution_1_CmmNnEx", "solution_1_CmmExtr","solution_1_CmmFshr",'ANCHOR','MARICUL',"mining"),labels=c('Bio',"CommunityNonExtractive" , "CommunityExtractive","Commercialfisheries" ,"Anchorage","Mariculture",'Mining'))
  
  land <- read_sf('SetUpExampleData/coast.shp')
  land <- st_transform(land, crs = st_crs(ras))
  ggplot()+
    geom_sf(data=st_union(r[r$solution == 'Bio',]),fill='#42C713',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'CommunityNonExtractive',]),fill='#3249DC',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'CommunityExtractive',]),fill='#4BCAC2',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'Commercialfisheries',]),fill='#8F4DEF',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'Anchorage',]),fill='#F4D03F',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'Mariculture',]),fill='#E88E3E',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'Mining',]),fill='black',col='transparent')+
    geom_sf(data=r,fill='transparent',col='grey50',linewidth=0.05)+
    geom_sf(data=land,col='black',fill='grey90')+
    coord_sf(xlim = c(46137 ,158130 ),ylim=c(-3869394 ,-3810394 ))+
    
    theme_bw() +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          strip.background=element_blank()) +
    xlab('Longitude')+
    ylab('Latitude')
  ggsave('Final_scenario_maps/EB_Man.tiff',dpi=300)
  
  #Empty - mandatory####
  rm(list=ls())
  ras <- readRDS('Final_prioritizations/EmptyBay_NoMan.rds')
  ras$ANCHOR <- 0
  ras$ANCHOR[ras$ANCHOR_ConsRes>0] <- 1
  ras$MARICUL <- 0
  ras$MARICUL[ras$MARICUL_ConsRes>0] <- 1
  ras$solution <- category_vector(
    ras[, c("solution_1_ConsRes","solution_1_CmmNnEx", "solution_1_CmmExtr","solution_1_CmmFshr",'ANCHOR','MARICUL',"solution_1_Mining")]
  )
  ggplot(ras,aes(fill=factor(solution)))+geom_sf()
  unique(ras$solution)
  r <- ras[, c("solution")]
  # plot solution
  r$solution[r$solution == 0] <- 'Not allocated'
  r$solution[r$solution == 1] <- 'solution_1_ConsRes'
  r$solution[r$solution == 2] <- 'solution_1_CmmNnEx'
  r$solution[r$solution == 3] <- 'solution_1_CmmExtr'
  r$solution[r$solution == 4] <- 'solution_1_CmmFshr'
  r$solution[r$solution == 5] <- 'ANCHOR'
  r$solution[r$solution == 6] <- 'MARICUL'
  r$solution[r$solution == 7] <- 'mining'
  r$solution <- factor(r$solution,levels= c("solution_1_ConsRes","solution_1_CmmNnEx", "solution_1_CmmExtr","solution_1_CmmFshr",'ANCHOR','MARICUL',"mining"),labels=c('Bio',"CommunityNonExtractive" , "CommunityExtractive","Commercialfisheries" ,"Anchorage","Mariculture",'Mining'))
  
  land <- read_sf('SetUpExampleData/coast.shp')
  land <- st_transform(land, crs = st_crs(ras))
  ggplot()+
    geom_sf(data=st_union(r[r$solution == 'Bio',]),fill='#42C713',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'CommunityNonExtractive',]),fill='#3249DC',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'CommunityExtractive',]),fill='#4BCAC2',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'Commercialfisheries',]),fill='#8F4DEF',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'Anchorage',]),fill='#F4D03F',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'Mariculture',]),fill='#E88E3E',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'Mining',]),fill='black',col='transparent')+
    geom_sf(data=r,fill='transparent',col='grey50',linewidth=0.05)+
    geom_sf(data=land,col='black',fill='grey90')+
    coord_sf(xlim = c(46137 ,158130 ),ylim=c(-3869394 ,-3810394 ))+
    
    theme_bw() +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          strip.background=element_blank()) +
    xlab('Longitude')+
    ylab('Latitude')
  ggsave('Final_scenario_maps/EB_NoMan.tiff',dpi=300)
  
  
  #Untouched - mandatory####
  rm(list=ls())
  ras <- readRDS('Final_prioritizations/Untouched_Man.rds')
  ras$ANCHOR <- 0
  ras$ANCHOR[ras$ANCHOR_ConsRes>0] <- 1
  ras$MARICUL <- 0
  ras$MARICUL[ras$MARICUL_ConsRes>0] <- 1
  
  shp <- read_sf('SetUpExampleData/MPAs.shp')
  ras$MPA <- shp$MPA1
  ras$concon <- 0
  ras$conres <- 0
  ras$concon[ras$MPA == "Controlled" ] <- 1
  ras$conres[ras$MPA %in% c("Restricted", "No-take") ] <- 1
  ras$solution_1_ConsRes[ras$concon == 1] <- 0
  ras$solution_1_ConsRes[ras$conres == 1] <- 0
  ras$solution <- category_vector(
    ras[, c("solution_1_ConsRes","conres",'concon',"solution_1_CmmNnEx", "solution_1_CmmExtr","solution_1_CmmFshr",'ANCHOR','MARICUL')]
  )
  
  unique(ras$solution)
  r <- ras[, c("solution")]
  # plot solution
  r$solution[r$solution == 0] <- 'Not allocated'
  r$solution[r$solution == 1] <- 'solution_1_ConsRes'
  r$solution[r$solution == 2] <- 'concon'
  r$solution[r$solution == 3] <- 'conres'
  r$solution[r$solution == 4] <- 'solution_1_CmmNnEx'
  r$solution[r$solution == 5] <- 'solution_1_CmmExtr'
  r$solution[r$solution == 6] <- 'solution_1_CmmFshr'
  r$solution[r$solution == 7] <- 'ANCHOR'
  r$solution[r$solution == 8] <- 'MARICUL'
  # r$solution[r$solution == 8] <- 'mining'
  r$solution <- factor(r$solution,levels= c("solution_1_ConsRes",'conres','concon',"solution_1_CmmNnEx", "solution_1_CmmExtr","solution_1_CmmFshr",'ANCHOR','MARICUL'),labels=c('Bio',"MPARestricted","MPAControlled","CommunityNonExtractive" , "CommunityExtractive","Commercialfisheries" ,"Anchorage","Mariculture"))
  
  land <- read_sf('SetUpExampleData/coast.shp')
  land <- st_transform(land, crs = st_crs(ras))
  
  ggplot()+
    geom_sf(data=st_union(r[r$solution == 'Bio',]),fill='#42C713',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'MPARestricted',]),fill='#D2F3AE',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'MPAControlled',]),fill='#BBD999',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'CommunityNonExtractive',]),fill='#3249DC',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'CommunityExtractive',]),fill='#4BCAC2',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'Commercialfisheries',]),fill='#8F4DEF',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'Anchorage',]),fill='#F4D03F',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'Mariculture',]),fill='#E88E3E',col='transparent')+
    geom_sf(data=r,fill='transparent',col='grey50',linewidth=0.05)+
    geom_sf(data=land,col='black',fill='grey90')+
    coord_sf(xlim = c(46137 ,158130 ),ylim=c(-3869394 ,-3810394 ))+
    
    theme_bw() +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          strip.background=element_blank()) +
    xlab('Longitude')+
    ylab('Latitude')
  ggsave('Final_scenario_maps/UnT_Man.tiff',dpi=300)
  
  #Untouched - mandatory####
  ras <- readRDS('Final_prioritizations/Untouched_NoMan.rds')
  ras$ANCHOR <- 0
  ras$ANCHOR[ras$ANCHOR_ConsRes>0] <- 1
  ras$MARICUL <- 0
  ras$MARICUL[ras$MARICUL_ConsRes>0] <- 1
  
  shp <- read_sf('SetUpExampleData/MPAs.shp')
  ras$MPA <- shp$MPA1
  ras$concon <- 0
  ras$conres <- 0
  ras$concon[ras$MPA == "Controlled" ] <- 1
  ras$conres[ras$MPA %in% c("Restricted", "No-take") ] <- 1
  ras$solution_1_ConsRes[ras$concon == 1] <- 0
  ras$solution_1_ConsRes[ras$conres == 1] <- 0
  # ggplot()+
  #   geom_sf(data=ras[ras$solution_1_Trnsprt > 0,],fill='blue')+
  #   geom_sf(data=ras[ras$ANCHOR > 0,],fill='transparent',col='hotpink')+
  #   geom_sf(data=ras[ras$MARICUL > 0,],fill='transparent',col='green')
  
  ras$solution <- category_vector(
    ras[, c("solution_1_ConsRes","conres",'concon',"solution_1_CmmNnEx", "solution_1_CmmExtr","solution_1_CmmFshr",'ANCHOR','MARICUL')]
  )
  
  unique(ras$solution)
  r <- ras[, c("solution")]
  # plot solution
  r$solution[r$solution == 0] <- 'Not allocated'
  r$solution[r$solution == 1] <- 'solution_1_ConsRes'
  r$solution[r$solution == 2] <- 'concon'
  r$solution[r$solution == 3] <- 'conres'
  r$solution[r$solution == 4] <- 'solution_1_CmmNnEx'
  r$solution[r$solution == 5] <- 'solution_1_CmmExtr'
  r$solution[r$solution == 6] <- 'solution_1_CmmFshr'
  r$solution[r$solution == 7] <- 'ANCHOR'
  r$solution[r$solution == 8] <- 'MARICUL'
  # r$solution[r$solution == 8] <- 'mining'
  r$solution <- factor(r$solution,levels= c("solution_1_ConsRes",'conres','concon',"solution_1_CmmNnEx", "solution_1_CmmExtr","solution_1_CmmFshr",'ANCHOR','MARICUL'),labels=c('Bio',"MPARestricted","MPAControlled","CommunityNonExtractive" , "CommunityExtractive","Commercialfisheries" ,"Anchorage","Mariculture"))
  
  land <- read_sf('SetUpExampleData/coast.shp')
  land <- st_transform(land, crs = st_crs(ras))
  
  ggplot()+
    geom_sf(data=st_union(r[r$solution == 'Bio',]),fill='#42C713',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'MPARestricted',]),fill='#D2F3AE',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'MPAControlled',]),fill='#BBD999',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'CommunityNonExtractive',]),fill='#3249DC',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'CommunityExtractive',]),fill='#4BCAC2',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'Commercialfisheries',]),fill='#8F4DEF',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'Anchorage',]),fill='#F4D03F',col='transparent')+
    geom_sf(data=st_union(r[r$solution == 'Mariculture',]),fill='#E88E3E',col='transparent')+
    geom_sf(data=r,fill='transparent',col='grey50',linewidth=0.05)+
    geom_sf(data=land,col='black',fill='grey90')+
    coord_sf(xlim = c(46137 ,158130 ),ylim=c(-3869394 ,-3810394 ))+
    
    theme_bw() +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          strip.background=element_blank()) +
    xlab('Longitude')+
    ylab('Latitude')
  ggsave('Final_scenario_maps/UnT_NoMan.tiff',dpi=300)
  
  
