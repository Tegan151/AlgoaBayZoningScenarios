#######################################################
#This code accompanies the manuscript:
# Scenario planning from the bottom up supports equitable and ecosystem-based approaches to marine spatial planning 
# Tegan Carpenter-Kling*, Hannah Truter, Anne Lemahieu, Bernadette Snow, Mia Strand, Nina Rivers, James Blignaut, Rozanne Bester, Lea Nupnau, Amanda T. Lombard
##Code created by Tegan Carpenter-Kling
##Institute for Coastal and Marine Research, Nelson Mandela University, South Africa
##Please direct queries to tegan.carpenterkling@gmail.com
#########################################################
##The code provides examples of how to:
#   -prepare spatial data in the R environment for RioritizR with zones
#   -extrapolate intensity data of fishing with an area to outside
#   -create cost layers for zones
#   -create a layer of 'locked in' planning units
#   -calculate area based targets as described by Holness et al. 2022 (doi: 10.1016/j.biocon.2022.109574)
#   -create and save inputs for prioritizR with zones
#########################################################

library(sf)
library(ggplot2)
library(terra)

setwd('~/SetUpExampleData/')

############################################################################################################
##1. Examples of how create feature layers
############################################################################################################

########################################
####1.1. Create planning unit layer#####
########################################

rm(list=ls())
#Read in polygon of planning domain
pu <- read_sf("PlanningDomain.shp")
#rasterize into regular sized planning units
r <- rast(pu,res=1000)#The unit of 'res' will depend on the CRS of your dataframe. It is suggested to work with project data (i.e. units = m)
v <- as.polygons(r)
writeVector(v, "pulayer.shp",overwrite=TRUE)

v <- read_sf("pulayer.shp")
colnames(v)
p <- st_intersection(v,pu)
ggplot(p)+geom_sf()
write_sf(p,'pulayer.shp',delete_dsn=T)


####################################################################
####1.2 Example of assigning a polygon feature to planning units####
####################################################################

rm(list=ls())
shp <- read_sf('pulayer.shp')#Read in planningunits shape file
temp <- read_sf('AB_Mariculture_legal_WGS84.shp')#Read in shapefile of polygon feature e.g. mariculture 
temp <- st_transform(temp,st_crs(shp))#project polygon
t <- st_intersection(temp,shp)#intersect polygon with planning units
ggplot(t)+geom_sf()
t2 <- shp[shp$FID %in% t$FID,] #subset planning units that intersected with polygon

t <- t %>%
  dplyr::mutate(area_meters = st_area(t), #workout area of feature in each intersected planning unit
                area_meters2 = st_area(t2), #workout each planning units area
                ProposedMariculture = as.numeric((area_meters/area_meters2)*100))#work out the percentage coverage 

t <- t[,c('ProposedMariculture')]

ggplot()+
  geom_sf(data=t,aes(fill=ProposedMariculture))

shp <- st_join(shp,t,left=T,largest=T) #join to original file
shp <- shp %>% st_sf(sf_column_name = 'geometry')
shp$ProposedMariculture[is.na(shp$ProposedMariculture)] <- 0
shp$ProposedMariculture[shp$ProposedMariculture<50] <- 0 #make all values less than 50%, 0 (i.e. absent)
shp$ProposedMariculture[shp$ProposedMariculture>0] <- 100# make all values greater than 50%, 100 (i.e. present)
ggplot()+
  geom_sf(data=shp,aes(fill=ProposedMariculture)) #plot

#########################################################################################################
####1.3 Example of removing fishing from MPA and extrapolate values to planning units outside of MPA#####
#########################################################################################################
rm(list=ls())

shp <- read_sf('small_pelagic_fishery.shp')
a <- read_sf('MPAs.shp')
ggplot()+
  geom_sf(data=shp,aes(fill = SMALLPE))+
  geom_sf(data=a[a$MPA1 == 'Restricted',],fill='transparent',col='hotpink')

shp1 <- st_join(shp,a,left=T,largest=T)#join MPA and fishing feature layers

# Remove fishing activity from restricted and no take MPAs using formula in manuscript
shp1$sa <- shp1$SMALLPE
shp1$sa[shp1$MPA1 %in%c('Restricted','No-take')] <- 0#Make all values of fishing intensity within MPA 0
shp1$sa2 <- shp1$sa*(sum(shp1$SMALLPE)/sum(shp1$sa))#Extrapolate these values outside of MPA
shp1$SMALLPE <- shp1$sa2
shp1 <- shp1[,!(names(shp1) %in% c('sa2','sa'))]
shp1 <- shp1[,!(colnames(shp1) %in% c('MPA1','MPA2','MPA'))]
ggplot()+
  geom_sf(data=shp1,aes(fill =SMALLPE))

write_sf(shp1,'small_pelagic_fishery_corrected.shp')
######################################################
###2. Create cost layers #############################
######################################################

#Priority conservation cost layer
rm(list=ls())

shp <- read_sf('features.shp')#Read in multipolygon which contains all of your feature layers, including those which are only being used to create cost layers

#create dataframe with feature names and assign each feature to a category
df <- as.data.frame(matrix(ncol=2,nrow=ncol(shp)-1))
colnames(df) <- c('MSP_id','category')
df$MSP_id <- colnames(shp)[!(colnames(shp) %in% c('geometry'))]

#Add categories
df$category[df$MSP_id %in% c('Llthgnt','Aaqudns','Aargyrzn','AbalnRf','AgISRS','AglDI','AglES','AglhBE','AglhERS','AglhsIS','AglhsMS','AglhSOS','AglhSRS','AglISS','AglMSRC','AglSMS','AgOSRC','AgSIS1','AgSMS1','Ainodrs','Ajapncs','AntTrnR','APBreed','Baitbll','BllShrk','BrnzWhl','BttlnDl','BWhale','Ccristcp','CGBreed','Cgibbcps','Clatcps','Cnasuts','Cnufar','ComDol','CrGnntA','CrPngnA','DamTrnF','Dhttntt','DmTrnNs','Dsargus','DTBreed','Emachnt','Emargnt','Erays','FishLrv','GanntUs','Gcrvdns','GllyShr','HammerJ','Hghprfl','Hmpbcks',"HBdlphn",'Hplblph','Lamia','Minke','Orca','Osharks','Pcmmrsn','PengnUs','Porodrm','Pprerbt','Prupstrs','Psaltrx','Pundlss','Raggies','Reef','Rholubi','RsTrnRs','SABAPbn','SABAPrc','Sdrbnns','Slaland','SqdSpwnn','SqudSpwn','Srays','SRghtWh','STBreed','Turtles','Urobnsn','WhtShrk','WrmTmAE','OthrTrR','SCP',"SandSharks"  )] <- 'Biodiversity'

df$category[df$MSP_id %in%c("ANCHOR") ] <- 'Anchorage'

df$category[df$MSP_id %in%c("DUMP") ] <- 'Dredge'

df$category[df$MSP_id %in% c("SHIPPIN","SHIPINT")] <- 'Shipping'

df$category[df$MSP_id %in% c("CUMULAT")] <- 'CUMULAT'

df$category[df$MSP_id %in% c("defence")] <- 'defence'

df$category[df$MSP_id %in% c("OffshrL")] <- 'mining'

df$category[df$MSP_id %in% c("MARICUL",'PrpsdMr')] <- 'ADZ'

df$category[df$MSP_id %in% c("Blstwtr", "wastwtr")] <- 'waste'

df$category[df$MSP_id %in%c("Brdwtch","CstWhlW","DiveRec","KiteBrd","LargBts","LifeSvC","OWSwim","RcrtnlB","CR_BH","CR_CC","CR_CR","CR_WS","ShrkCgD","SmallBt","SUP","Surfing","SurfSki") ] <- 'RecreationNonExtract'

df$category[df$MSP_id %in% c("SQUID", "TRAWL", "SHARKFI", "LINEFIS","SMALLPE",'SSF')] <- 'CommFishery'

df$category[df$MSP_id %in%c("FshTrps","CH_CC","CH_CR","CH_WS","CL_BH","CL_CC","CL_CR","CL_WS" ,"MddnsnP","OthrHrS","CS_BH","CS_WS","Wrecks" ) ] <- 'Heritage'

df$category[df$MSP_id %in%c('RcrtnlK','RcrtnlSh','RcrtnlSk','RcrtnlSp','RcrtnSC') ] <- 'RecreationExtract'

df$category[df$MSP_id %in%c("LS_BH","LS_CC","LS_CR","LS_WS","SbsstnF") ] <- 'SubFish'
write.csv(df,'featuresBYcategory.csv',row.names = F)#Save for later

#Create costs layer for each zone by accumulating waited layers (weighted by 0, 0.5 or 1 according to table 4 in manuscript) and then scale between 0 and 100

#Conservation priority zone
df$ConsPriWaiting <- NA
df$ConsPriWaiting[df$category %in% c('RecreationExtract','CUMULAT','Dredge','Anchorage','CommFishery','ADZ','SubFish','waste',"Shipping","mining")] <- 1
df$ConsPriWaiting[df$category %in% c('RecreationNonExtract','Heritage')] <- 0.5
df$ConsPriWaiting[df$category %in% c('Biodiversity',"defence")] <- 0
summary(is.na(df$ConsPriWaiting))
r <- shp[,c(df$MSP_id[df$ConsPriWaiting ==1])]
r1 <- shp[,c(df$MSP_id[df$ConsPriWaiting ==0.5])]
r1[,colnames(r1)[!(colnames(r1)  %in% c('geometry'))]] <-st_drop_geometry(r1)*0.5
r <- st_join(r,r1,left=T,largest=T)
r$ConsPri <- rowSums(st_drop_geometry(r),na.rm = T)
head(r)
r <- r[,c('ConsPri')]
ggplot(r)+
  geom_sf(aes(fill=ConsPri))+
  scale_fill_viridis_c(option = "H")+
  ggtitle('Conservation priority zone cost layer')
costlayers <- r ##save as new layer


#Community non-extraction priority zone

df$CommNonExtrWaiting <- NA
df$CommNonExtrWaiting[df$category %in% c('Anchorage','ADZ','RecreationExtract','CommFishery','SubFish','waste','CUMULAT')] <- 1
df$CommNonExtrWaiting[df$category %in% c('Dredge','RecreationNonExtract','Heritage','Biodiversity')] <- 0
df$CommNonExtrWaiting[df$category %in% c('mining','Shipping','defence')] <- 0.5
summary(is.na(df))
r <- shp[,c(df$MSP_id[df$CommNonExtrWaiting ==1])]
r1 <- shp[,c(df$MSP_id[df$CommNonExtrWaiting ==0.5])]
r1[,colnames(r1)[!(colnames(r1)  %in% c('geometry'))]] <- st_drop_geometry(r1)*0.5
r <- st_join(r,r1,left=T,largest=T)
r$CommNonExtr <- rowSums(st_drop_geometry(r),na.rm = T)
head(r)
r <- r[,c('CommNonExtr')]
ggplot(r)+
  geom_sf(aes(fill=CommNonExtr))+
  scale_fill_viridis_c(option = "H")+
  ggtitle('Community non-extractive priority zone cost layer')

costlayers <- st_join(costlayers, r,left=T,largest=T )


#Community extractive priority zone

df$CommExtrWaiting <- NA
df$CommExtrWaiting[df$category %in% c('Dredge','Anchorage','ADZ')] <- 1
df$CommExtrWaiting[df$category %in% c('CUMULAT','RecreationNonExtract','Heritage','RecreationExtract','SubFish')] <- 0
df$CommExtrWaiting[df$category %in% c('CommFishery','Biodiversity','mining','waste','Shipping','defence')] <- 0.5

r <- shp[,c(df$MSP_id[df$CommExtrWaiting ==1])]
r1 <- shp[,c(df$MSP_id[df$CommExtrWaiting ==0.5])]
r1[,colnames(r1)[!(colnames(r1)  %in% c('geometry'))]] <- st_drop_geometry(r1)*0.5
r <- st_join(r,r1,left=T,largest=T)
r$CommExtr <- rowSums(st_drop_geometry(r),na.rm = T)
head(r)
r <- r[,c('CommExtr')]
ggplot(r)+
  geom_sf(aes(fill=CommExtr))+
  scale_fill_viridis_c(option = "H")+
  ggtitle('Community extractive priority zone cost layer')

costlayers <- st_join(costlayers, r,left=T,largest=T )



#Commercial fishery priority zone###
df$CommFishery <- NA
df$CommFishery[df$category %in% c('ADZ','Dredge','Anchorage')] <- 1
df$CommFishery[df$category %in% c('CUMULAT','CommFishery')] <- 0
df$CommFishery[df$category %in% c('RecreationNonExtract','Heritage','Biodiversity','RecreationExtract','SubFish','mining','Shipping','waste','defence')] <- 0.5

r <- shp[,c(df$MSP_id[df$CommFishery ==1])]
r1 <- shp[,c(df$MSP_id[df$CommFishery ==0.5])]
r1[,colnames(r1)[!(colnames(r1)  %in% c('geometry'))]] <- st_drop_geometry(r1)*0.5
r <- st_join(r,r1,left=T,largest=T)
r$CommFishery <- rowSums(st_drop_geometry(r),na.rm = T)
head(r)
r <- r[,c('CommFishery')]

ggplot(r)+
  geom_sf(aes(fill=CommFishery))+
  scale_fill_viridis_c(option = "H")+
  ggtitle('Commercial fisheries priority zone cost layer')

costlayers <- st_join(costlayers, r,left=T,largest=T )



#cost layer for mariculutre and ship to ship bunkering zone###
ras1 <- shp
ras1$MS <- 100
ras1$MS <- ras1$MS - ras1$ANCHOR
ras1$MS <- ras1$MS - ras1$MARICUL

r <- ras1[,c('MS')]

ggplot(r)+
  geom_sf(aes(fill=MS))+
  ggtitle('Ship to ship bunkering and mariculutre zone cost layer')

costlayers <- st_join(costlayers, r,left=T,largest=T )

#Scale cost layers between 0 - 100
min <- min(c(costlayers$ConsPri,costlayers$CommNonExtr ,costlayers$CommExtr,costlayers$CommFishery))
max <- max(c(costlayers$ConsPri,costlayers$CommNonExtr ,costlayers$CommExtr,costlayers$CommFishery))
costlayers$ConsPri <- ((costlayers$ConsPri - min) / (max - min)) *100
costlayers$CommNonExtr  <- ((costlayers$CommNonExtr  - min )/ (max - min)) *100
costlayers$CommExtr <- ((costlayers$CommExtr - min) / (max - min)) *100
costlayers$CommFishery <- ((costlayers$CommFishery - min) / (max - min)) *100

# Add 1 to all values so there are no zeros in cost layers
costlayers$ConsPri <-  costlayers$ConsPri+1
costlayers$CommExtr <- costlayers$CommExtr+1
costlayers$CommNonExtr <- costlayers$CommNonExtr+1
costlayers$CommFishery <- costlayers$CommFishery+1
costlayers$MS <- costlayers$MS+1

write_sf(costlayers,'CostPolygon.shp', delete_dsn = T)
######################################################################
###3. Create pu layer with features locked in or out for each zone#####
######################################################################
rm(list=ls())
polygon <- read_sf('features.shp')
t <- read_sf('MPAs.shp')
head(t)
# Create a layer for each zone with no planning units locked in (i.e. FALSE)
t$lock_in_ConsPri <- FALSE
t$lock_in_CmmExtr <- FALSE
t$lock_in_CmmNnEx <- FALSE
t$lock_in_CmmFshr <- FALSE
t$lock_in_MS <- FALSE
# Lock in MPAs into conservation priority zone
t$lock_in_ConsPri[t$MPA1 %in% c('Restricted','No-take')]  <- TRUE
t$lock_in_ConsPri[t$MPA1 %in% c('Controlled')]<- TRUE
# Lock in ship to ship bunkering and mariculutre into MS zone and Fishtraps into non-extractive priority zone
polygon <- polygon[,c('ANCHOR','FshTrps','MARICUL')]
t <- st_join(polygon,t,left=T,largest=T)
head(t)
t$lock_in_MS[t$ANCHOR>0] <- T
t$lock_in_MS[t$MARICUL>0] <- T
t$lock_in_CmmNnEx[t$FshTrps>0] <- T

polygon <- t[,c( 'lock_in_ConsPri','lock_in_CmmNnEx','lock_in_CmmExtr','lock_in_CmmFshr','lock_in_MS', "geometry"    )]

polygon$lock_in_ConsPri[is.na(polygon$lock_in_ConsPri)] <- FALSE
polygon$lock_in_CmmNnEx[is.na(polygon$lock_in_CmmNnEx)] <- FALSE
polygon$lock_in_CmmExtr[is.na(polygon$lock_in_CmmExtr)] <- FALSE
polygon$lock_in_CmmFshr[is.na(polygon$lock_in_CmmFshr)] <- FALSE
polygon$lock_in_MS[is.na(polygon$lock_in_MS)] <- FALSE



polygon$lock_in_ConsPri <- as.logical(polygon$lock_in_ConsPri)
polygon$lock_in_CmmNnEx <- as.logical(polygon$lock_in_CmmNnEx)
polygon$lock_in_CmmExtr <- as.logical(polygon$lock_in_CmmExtr)
polygon$lock_in_CmmFshr <- as.logical(polygon$lock_in_CmmFshr)
polygon$lock_in_MS <- as.logical(polygon$lock_in_MS)


# Join to cost layer
p<- read_sf('CostPolygon.shp')
polygon <- st_join(p,polygon,left=T,largest=T)
head(polygon)
####make sure that areas arent locked in and out ie a bit of the anchorage overlaps with the mpa
for(i in 1:nrow(polygon)){
  if(polygon$lock_in_ConsPri[i]  && polygon$lock_in_MS[i] == T){
    polygon$lock_in_ConsPri[i] <- F
    
  }
}

write_sf(polygon,'CostPolygon_and_locks.shp',delete_dsn=T)

############################################################################
###4. Calculating targets using area-based target setting detailed in Holness et al. 2022 (doi: 10.1016/j.biocon.2022.109574)
############################################################################
rm(list=ls())

shp <- read_sf('features.shp')#Read in multipolygon which contains all of your feature layers
shp <- shp[,!(colnames(shp) %in% c("CUMULAT", "DUMP"  ,  "SHIPINT", "OffshrL", "defence", "Blstwtr", "wastwtr", "SHIPPIN",'PrpsdMr'))]# remove features that were only used to create cost layers are not features in the prioritization 

#create dataframe to save new targets to
df <- read.csv('featuresBYcategory.csv')
df$Total_PU <- NA
df$Target_PU <- NA
df <- df[!(df$MSP_id %in% c("CUMULAT", "DUMP"  ,  "SHIPINT", "OffshrL", "defence", "Blstwtr", "wastwtr", "SHIPPIN",'PrpsdMr')),  ]


######Calculate targets using the area-based formula from Holness et al. 2022#
id <- unique(df$MSP_id)
for(i in 1:length(id)){
  r <- shp[,id[i]]
  colnames(r) <- c('feature','geometry')
  if(class(r$feature)=="character"){r$feature <- as.numeric(r$feature)}
  r$feature[is.na(r$feature)] <- 0
  r <- r %>%
    dplyr::mutate(area_meters = st_area(r),
                  area = as.numeric((area_meters/1000000)))
  r$feature <- r$feature*r$area
  a <- sum(r$feature[r$feature>0])
  df$Total_PU[df$MSP_id == id[i]] <- a
  if(a <= 1000){
    df$Target_PU[df$MSP_id == id[i]] <- a*0.6
  }
  if(a > 1000 & a <= 5000 ){
    df$Target_PU[df$MSP_id == id[i]] <- 1000*0.6 + (a-1000)*0.5
  }
  if(a > 5000 & a <= 10000 ){
    df$Target_PU[df$MSP_id == id[i]] <- 1000*0.6 + 4000*0.5 + (a-5000)*0.4
  }
  if(a > 10000 & a <= 50000 ){
    df$Target_PU[df$MSP_id == id[i]] <- 1000*0.6 + 4000*0.5 + 5000*0.4 + (a-10000)*0.3
  }
  
  if(a > 50000 ){
    df$Target_PU[df$MSP_id == id[i]] <- 1000*0.6 + 4000*0.5 + 5000*0.4 + 10000*0.3+ (a-20000)*0.2
  }
}

df$target <- (df$Target_PU/df$Total_PU)*100

#Pre-set targets for mariculture and ship-to-ship bunkering##
df$target[df$MSP_id %in% c('ANCHOR','MARICUL')] <- 100

#add zone targets#
df$Biodiversity <- NA
df$Biodiversity[df$category =='Biodiversity'] <- 1
df$Biodiversity[is.na(df$Biodiversity)] <- 0

df$CommmunityUseExtractive <- NA  
df$CommmunityUseExtractive[df$category %in%  c('RecreationExtract','SubFish')] <- 1
df$CommmunityUseExtractive[is.na(df$CommmunityUseExtractive)] <- 0

df$CommunityUseNonExtractive <- NA
df$CommunityUseNonExtractive[df$category ==  'Heritage'] <- 1
df$CommunityUseNonExtractive[df$category ==  'RecreationNonExtract'] <- 1
df$CommunityUseNonExtractive[is.na(df$CommunityUseNonExtractive)] <- 0

df$CommercialFisheries <- NA
df$CommercialFisheries[df$category ==  'CommFishery'] <- 1
df$CommercialFisheries[is.na(df$CommercialFisheries)] <- 0

df$MS <- NA
df$MS[df$MSP_id %in% c('ANCHOR','MARICUL')] <- 1
df$MS[is.na(df$MS)] <- 0
summary(is.na(df))


##Recalculate targets for activities that are not zoned for the conservation priority zone but can occur in the current Addo and Sardinia MPA (see methods)


polygon <- read_sf('CostPolygon_and_locks.shp')

summary(is.na(polygon))
colnames(polygon) <- c( "ConsPri","CmmNnEx","CmmExtr","CmmFshr",'Trnsprt','lock_in_ConsPri','lock_in_CmmNnEx','lock_in_CmmExtr','lock_in_CmmFshr','lock_in_MS', "geometry"   )

polygon$lock_in_ConsPri <- as.logical(polygon$lock_in_ConsPri)
polygon$lock_in_CmmNnEx <- as.logical(polygon$lock_in_CmmNnEx)
polygon$lock_in_CmmExtr <- as.logical(polygon$lock_in_CmmExtr)
polygon$lock_in_CmmFshr <- as.logical(polygon$lock_in_CmmFshr)
polygon$lock_in_MS <- as.logical(polygon$lock_in_MS)

polygon <- st_join(polygon,shp,left=T,largest=T)


polygon <-polygon[,c( 'lock_in_ConsPri','lock_in_CmmNnEx', names(shp))]

#Crop polygon to MPA restricted and controlled
res <- polygon[polygon$lock_in_ConsPri == T,]
fish <- polygon[polygon$lock_in_CmmNnEx == T,]

df <- df[df$MSP_id %in% names(polygon),]

id <- names(polygon)[!(names(polygon) %in% c("lock_in_ConsPri",'lock_in_CmmNnEx',"geometry"))]
df$conres_area <- NA
df$fish_area <- NA

for(i in 1:length(id)){
  r <- res[,id[i]]
  
  colnames(r) <- c('feature','geometry')
  
  df$conres_area[df$MSP_id == id[i]] <- sum(r$feature[r$feature>0])
  
  f <- fish[,id[i]]
  colnames(f) <- c('feature','geometry')
  df$fish_area[df$MSP_id == id[i]] <- sum(f$feature[f$feature>0])
  rm(c,r,f)
}
head(df)

#Change targeted PU for those targets that were preset
df$Target_PU[df$MSP_id %in% c('MARICUL','ANCHOR')] <- df$Total_PU[df$MSP_id %in% c('MARICUL','ANCHOR')]

#Now for all activities that are allowed in the MPA restricted and controlled zone, remove these PU from target_PU

#MPA Restricted
temp <- df[df$category %in% c('Heritage','RecreationNonExtract'),]
df <-  df[!(df$category %in% c('Heritage','RecreationNonExtract')),]
temp$Target_PU <- temp$Target_PU - temp$conres_area 
df <- rbind(temp,df)
a <- df$MSP_idf[df$target<0]

#Remove features were more than targeted PUs are already met by the restricted MPA
df <- df[df$Target_PU>0,]


#Remove features were more than targeted PUs are already met by the controlled MPA

temp <- df[df$category %in% c('CommFishery','SubFish','RecreationExtract'),]
df <-  df[!(df$category %in% c('CommFishery','SubFish','RecreationExtract')),]
temp$Target_PU <- temp$Target_PU - temp$conres_area 
temp[temp$Target_PU<0,]
#Remove features were more than targeted PUs are already met by the controlled MPA
temp <- temp[temp$Target_PU>0,]
df <- rbind(temp,df)

#Check for conservation layers were entire target is met in MPA restricted and controlled existing zones
temp <- df[df$category %in% c('Conservation'),]
temp$Target_PU <- temp$Target_PU - temp$conres_area
temp <- temp[temp$Target_PU<0,]
df <- df[!(df$MSP_id %in% temp$MSP_id),]

#Remove conservationand non extractive recreational features were more than targeted PUs are already met by the fish trips

temp <- df[df$category %in% c('Conservation','RecreationNonExtract','Heritage'),]
temp$Target_PU <- temp$Target_PU - temp$fish_area 
temp <- temp[temp$fish_area>0,]
temp <- temp[temp$MSP_id != 'FshTrps',]
df <- df[!(df$MSP_id %in% temp$MSP_id),]
df <- rbind(temp,df)
df <- df[!(df$MSP_id %in% c('AntTrnR','RsTrnRs','OthrTrR')),]
temp <- temp[temp$Target_PU<0,]
df <- df[!(df$MSP_id %in% temp$MSP_id),]
write.csv(df,'targets.csv',row.names = F)
#############################################
##5. Create and save input for prioritizR####
#############################################
rm(list = ls())
df <- read.csv('targets.csv')
ras <- read_sf('features.shp')

#Absolute targets
df$zone_1 <- round((df$Biodiversity*df$Target_PU ),0)
df$zone_2 <- round((df$CommunityUseNonExtractive*df$Target_PU ),0)
df$zone_3 <- round((df$CommmunityUseExtractive*df$Target_PU ),0)
df$zone_4 <- round((df$CommercialFisheries*df$Target_PU ),0)
df$zone_5 <- round((df$MS*df$Target_PU ),0)
head(df)

ras <- ras[,names(ras) %in% df$MSP_id]



#join feature data to polygon

p1 <- st_drop_geometry(ras)
p2 <- st_drop_geometry(ras)
p3 <- st_drop_geometry(ras)
p4 <- st_drop_geometry(ras)
p5 <- st_drop_geometry(ras)



names(p1) <- paste(names(ras),'ConsPri',sep='_')
names(p2) <- paste(names(ras),'CmmNnEx',sep='_')
names(p3) <- paste(names(ras),'CmmExtr',sep='_')
names(p4) <- paste(names(ras),'CmmFshr',sep='_')
names(p5) <- paste(names(ras),'MS',sep='_')

polygon <- read_sf('CostPolygon_and_locks.shp')
colnames(polygon) <- c( "ConsPri","CmmNnEx","CmmExtr","CmmFshr",'MS','lock_in_ConsPri','lock_in_CmmNnEx','lock_in_CmmExtr','lock_in_CmmFshr','lock_in_MS', "geometry"   )

polygon$lock_in_ConsPri[is.na(polygon$lock_in_ConsPri)] <- 0
polygon$lock_in_CmmNnEx[is.na(polygon$lock_in_CmmNnEx)] <- 0
polygon$lock_in_CmmExtr[is.na(polygon$lock_in_CmmExtr)] <- 0
polygon$lock_in_CmmFshr[is.na(polygon$lock_in_CmmFshr)] <- 0
polygon$lock_in_MS[is.na(polygon$lock_in_MS)] <- 0


polygon$lock_in_ConsPri <- as.logical(polygon$lock_in_ConsPri)
polygon$lock_in_CmmNnEx <- as.logical(polygon$lock_in_CmmNnEx)
polygon$lock_in_CmmExtr <- as.logical(polygon$lock_in_CmmExtr)
polygon$lock_in_CmmFshr <- as.logical(polygon$lock_in_CmmFshr)
polygon$lock_in_MS <- as.logical(polygon$lock_in_MS)
polygon <- cbind(polygon, p1,p2,p3,p4,p5)



saveRDS(polygon,'polygon.rds')
saveRDS(p1,'p1.rds')
saveRDS(p2,'p2.rds')
saveRDS(p3,'p3.rds')
saveRDS(p4,'p4.rds')
saveRDS(p5,'p5.rds')


#Re-order targets and save as matrix
temp <- as.data.frame(matrix(ncol=ncol(df),nrow=nrow(df)))
colnames(temp) <- colnames(df)
feature <- names(st_drop_geometry(ras))

for(i in 1:length(feature)){
  temp[i,] <- df[df$MSP_id == feature[i],]
}

d <- temp


targets <- matrix(
  cbind(d$zone_1,d$zone_2,d$zone_3,d$zone_4,d$zone_5),
  nrow=nrow(d),
  ncol=5,
  dimnames = list(
    d$MSP_id,c( "ConsPri","CmmNnEx","CmmExtr","CmmFshr",'MS')
  )
)


saveRDS(targets,'targets.rds')




