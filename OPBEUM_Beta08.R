list.of.packages <- c("data.table","DBI","doParallel","e1071",
                      "FNN","foreach","foreign","geosphere","ggmap","ggplot2","ggrepel","gmapsdistance",
                      "graphics","Grid2Polygons","gridExtra","Imap","kernlab","maps","maptools","nnet",
                      "osmar","osrm","parallel","plm","plyr","progress","rangeMapper","raster","RDSTK","rgdal",
                      "rgeos","RPostgreSQL","RSQLite","snow","sp","spatstat","stringr","tm","utils","wordcloud")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(nnet)
library(graphics)
library(e1071)
library(plm)
library(kernlab)
library(data.table)
library(stringr)
library(parallel)
library(doParallel)
library(snow)
library(sp)
library(spatstat)
library(maptools)
library(rgeos)
library(ggplot2)
library(plyr)
library(RSQLite)
library(DBI)
library(foreign)
library(rgdal)
library(maps)
library(geosphere)
library(gridExtra)
library(ggmap)
library(wordcloud)
library(tm)
library(ggrepel)
library(gmapsdistance)
library(Imap)
library(utils)
library(osmar)
library(RPostgreSQL)
library(raster)
library(progress)
library(spatstat)
library(RDSTK)
library(osrm)
library(Grid2Polygons)
library(FNN)
library(rangeMapper)
library(foreach)
library(dplyr)

require(compiler)
enableJIT(3)

basedir<-"R://Data/Projects/Purple_Line"

setwd(basedir)

##Using cell centroids for distance or something else? (default is cell cenntroid, if using other shapefile enter name below)
cell_centroid<-"Name"   #<--name shapefile here

##name of the shapefile for the study area
studyarea<-"purpleline_1mi_buffer"   #<--name shapefile here

##What size grid do you want (in meters)
gridsize<-50

######################
##   BE Variables   ##
######################

#connect to databases
drv <- dbDriver("PostgreSQL")
con <- dbConnect(drv, user = "postgres", dbname = "osm", host = "localhost")
con_rt <- dbConnect(drv, user = "postgres", dbname = "opbeumDB", host = "localhost")


#----> CREATE GRID CELLS<----

  ## Set-up Inputs/OUTPUTS

# ##ODs from real station data
#   GRID_BIKESHARE_REALODs = fread(paste0("Input/","cabi_OD_01.csv"), sep=",",header = TRUE)
#     names(GRID_BIKESHARE_REALODs)[names(GRID_BIKESHARE_REALODs)=="Origin"] <- "TERMINAL_N_O"
#     names(GRID_BIKESHARE_REALODs)[names(GRID_BIKESHARE_REALODs)=="Destination"] <- "TERMINAL_N_D"


   #study area shapefile
  place.spr.original<-readOGR(dsn = "Input", layer = paste0(studyarea))
  place.spr.original.proj<-spTransform( place.spr.original, CRS( "+init=epsg:3347" ) ) 
  #place.spr <- spTransform(place.spr.original,CRS("+proj=longlat +datum=WGS84"))
  
  # #Bike-share station locations
  # bs.stations<-readOGR(dsn = "Input", layer = "Capital_Bike_Share_Locations") 

# place.spr<-place.sp #use for specific area
# place.spr<-place.sp[z,] #use if shape loop

### define SpatialGrid object
place.bb<-bbox(place.spr.original.proj)
place.cs <- c(3.28084, 3.28084)*gridsize  # cell size 300m x 300m
# 1 ft = 3.28084 m
place.cc <- place.bb[, 1] + (place.cs /2)  # cell offset
place.cd <- ceiling(diff(t(place.bb))/place.cs)  # number of cells per direction
place.grd <- GridTopology(cellcentre.offset=place.cc, cellsize=place.cs, cells.dim=place.cd) #create the grid

##make grid topology into SP dataframe
place_sp_grd <- SpatialGridDataFrame(place.grd,
                                     data=data.frame(id=1:prod(place.cd)),
                                     proj4string=CRS(proj4string(place.spr.original.proj)))

##convert SP dataframe to polygons
#all_cells<-Grid2Polygons(place_sp_grd) #depricated
all_cells<-as(place_sp_grd, "SpatialPolygonsDataFrame") 


##clip cells to polygon
all_cells <- gIntersection(place.spr.original.proj, all_cells, byid = TRUE, drop_lower_td = TRUE)

######SWAP PROJECTED WITH WGS 
#transform to WGS
all_cells.proj <- all_cells
all_cells<- spTransform(all_cells,CRS("+proj=longlat +datum=WGS84"))

all_cells_id <- as.data.frame(matrix(0, ncol = 1, nrow = length(all_cells)))
all_cells_id$V1<- 1:nrow(all_cells_id)
colnames(all_cells_id)[1] <- "PageNumber"


# Extract polygon ID's
all_cells_id <- sapply(slot(all_cells, "polygons"), function(x) slot(x, "ID"))

# Create dataframe with correct rownames
all_cells.df <- data.frame( ID=1:length(all_cells), row.names = all_cells_id)

all_cells <- SpatialPolygonsDataFrame(all_cells,all_cells.df)

all_cells <-spTransform(all_cells, CRS("+proj=longlat +datum=WGS84"))
colnames(all_cells@data)[1] <- "PageNumber"

if(cell_centroid=="Name"){
#get the centroid of the cells for routing
origins<- gCentroid(all_cells,byid=TRUE)
origins <- SpatialPointsDataFrame(origins,all_cells.df)
  colnames(origins@data)[1] <- "PageNumber"

}else{
## Set-up Inputs/OUTPUTS
origins<-readOGR(dsn = "Input", layer = cell_centroid) #input origins shapefile
origins <-spTransform(origins, CRS("+proj=longlat +datum=WGS84"))
}


#writeOGR(obj=all_cells, dsn=".", layer="all_cells", driver="ESRI Shapefile") 

##
#Build Itersection Query
##
INT_FINAL<- data.frame(ref_count=integer(),
                       lat=numeric(), 
                       lon=numeric(),
                       PageNumber=integer(),
                       stringsAsFactors=FALSE) 

time1<-system.time({
  
  #Chunk grid cells 
  Splitrows=1000
  splitnumb_raw<-(nrow(all_cells)/Splitrows)
  if(nrow(all_cells)<Splitrows){
    Splitrows<-nrow(all_cells)
    splitnumb_l<-1
    leftover<-0
  }else{
    splitnumb_l<-floor(nrow(all_cells)/Splitrows)
    leftover<-splitnumb_raw-floor(splitnumb_l)
  }
  
  #Loops through cell chunks
  for(i in 1:splitnumb_l) {
    if(i==splitnumb_l){
      X1=(i-1)*Splitrows
      X2=Splitrows*i+Splitrows*leftover
    }else{i
      X1=(i-1)*Splitrows
      X2=Splitrows*i
    }
    
    #ith grid cell (row) bounding box and corners
    all_cells_cut<-all_cells[X1:X2,]
    all_cells_cut <-spTransform(all_cells, CRS("+proj=longlat +datum=WGS84"))
    
    bb<-bbox(all_cells_cut)
    
    ##Find intersections & Cul da sacks
    select<-"SELECT 
    ref_count,
    ST_Y(ST_Transform(geom_vertex,4326)) AS lat, 
    ST_X(ST_Transform(geom_vertex,4326)) AS lon 
    FROM at_2po_vertex
    WHERE  
    geom_vertex && 
    ST_MakeEnvelope("
    
    bb2<-paste(bb[1,1],",",bb[2,1],",",bb[1,2],",",bb[2,2],",","4283")
    q_tsig <- paste(select,bb2,")")
    
    #Pull SQL result
    tsig_result <- dbGetQuery(con_rt, q_tsig)
    
    #Convert node XY to spatial data frame
    node.spdf <- SpatialPointsDataFrame(coords = tsig_result[,c(3,2)], data = tsig_result,
                                        proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
    
    #make point coordinate system same as cells
    proj4string(node.spdf) <- CRS("+proj=longlat +ellps=WGS84") 
    proj4string(all_cells_cut) <- CRS("+proj=longlat +ellps=WGS84") 
    
    #Use just the PageNumber column
    all_cells_cut_red<-all_cells_cut[,1]
    
    #spatial join nodes with cells
    PageNumber = over(node.spdf,all_cells_cut_red)
    tsig_result2<-cbind(tsig_result,PageNumber)
    INT_FINAL=rbind(INT_FINAL,tsig_result2)
    
  }
})

print(time1)

#Add 1 to node directions
INT_FINAL$ref_count=INT_FINAL$ref_count+1

#Split intersections from cul de sacs
INTERS <- plyr::count(INT_FINAL[ which(INT_FINAL$ref_count>1), ], c('PageNumber'))
colnames(INTERS)[2] <- "Intersections"
INTERS3_4 <- plyr::count(INT_FINAL[ which(INT_FINAL$ref_count>2 & INT_FINAL$ref_count<5), ], c('PageNumber'))
colnames(INTERS3_4)[2] <- "Intersections3_4way"
INTERS5 <- plyr::count(INT_FINAL[ which(INT_FINAL$ref_count>4), ], c('PageNumber'))
colnames(INTERS5)[2] <- "Intersections_5plus"
INTCDS <- plyr::count(INT_FINAL[ which(INT_FINAL$ref_count==1), ], c('PageNumber'))
colnames(INTCDS)[2] <- "CulDeSacs"
LINKS <-aggregate(cbind(INT_FINAL$ref_count)~INT_FINAL$PageNumber, data=INT_FINAL, sum, na.rm=TRUE)
names(LINKS) <- c("PageNumber", "LINKS")

#Merge with full cell set
all_cells.df<-as.data.frame(all_cells[,1])
int.temp<-merge(all_cells.df,INTERS,all=T)
int.temp<-merge(int.temp,INTERS3_4,all=T)
int.temp<-merge(int.temp,INTERS5,all=T)
int.temp<-merge(int.temp,INTCDS,all=T)
int.temp<-merge(int.temp,LINKS,all=T)

#Convert NAs to 0s
int.temp[c("Intersections", "CulDeSacs","LINKS")][is.na(int.temp[c("Intersections", "CulDeSacs","LINKS")])] <- 0

int.temp$INTDns<-int.temp$Intersections/2.223948
int.temp$CDSDns<-int.temp$CulDeSacs/2.223948
int.temp$LINKDns<-int.temp$LINKS/2.223948
int.temp$CNR<-int.temp$LINKS/(int.temp$LINKS+int.temp$CulDeSacs)
int.temp$ALPHA<-abs((int.temp$LINKS-(int.temp$LINKS+int.temp$CulDeSacs)+1)/(2*(int.temp$LINKS+int.temp$CulDeSacs)-5))
int.temp$BETA<-int.temp$LINKS/(int.temp$Intersections+int.temp$CulDeSacs)
int.temp$GAMMA<-int.temp$LINKS/(3*((int.temp$Intersections+int.temp$CulDeSacs)-2))
int.temp$CYCL<-int.temp$LINKS-(int.temp$Intersections+int.temp$CulDeSacs)+2

###Get LEHD Data
LEHD_FINAL<- data.frame(tot_jobs=numeric(),
                        n11=numeric(),
                        n21=numeric(),
                        n22=numeric(),
                        n23=numeric(),
                        n31_33=numeric(),
                        n42=numeric(),
                        n44_45=numeric(),
                        n48_49=numeric(),
                        n51=numeric(),
                        n52=numeric(),
                        n53 =numeric(),
                        n54=numeric(),
                        n55=numeric(),
                        n56=numeric(),
                        n61=numeric(),
                        n62=numeric(),
                        n71=numeric(),
                        n72=numeric(),
                        n81=numeric(),
                        n92=numeric(),
                        PageNumber=integer(),
                        BlockAcres=numeric(),
                        Blockinter=numeric(),
                        stringsAsFactors=FALSE) 
time_LEHD<-system.time({
  
  #Chunk grid cells 
  Splitrows=100
  splitnumb_raw<-(nrow(all_cells)/Splitrows)
  if(nrow(all_cells)<Splitrows){
    Splitrows<-nrow(all_cells)
    splitnumb_l<-1
    leftover<-0
  }else{
    splitnumb_l<-floor(nrow(all_cells)/Splitrows)
    leftover<-splitnumb_raw-floor(splitnumb_l)
  }
  
  #Loops through cell chunks
  for(i in 1:splitnumb_l) {
    if(i==splitnumb_l){
      X1=(i-1)*Splitrows
      X2=Splitrows*i+Splitrows*leftover
    }else{i
      X1=(i-1)*Splitrows
      X2=Splitrows*i
    }
    
    #ith grid cell (row) bounding box and corners
    #all_cells <-spTransform(all_cells, CRS("+proj=longlat +datum=WGS84"))
    all_cells_cut<-all_cells[X1:X2,]
    
    bb<-bbox(all_cells_cut)
    
    ##Find employment
    select<-"SELECT 
    C000_y  AS  tot_jobs,
    CNS01_y AS	n11,
    CNS02_y AS	n21,
    CNS03_y AS	n22,
    CNS04_y AS	n23,
    CNS05_y AS	n31_33,
    CNS06_y AS	n42,
    CNS07_y AS	n44_45,
    CNS08_y AS	n48_49,
    CNS09_y AS	n51,
    CNS10_y AS	n52,
    CNS11_y AS	n53, 
    CNS12_y AS	n54,
    CNS13_y AS	n55,
    CNS14_y AS	n56,
    CNS15_y AS	n61,
    CNS16_y AS	n62,
    CNS17_y AS	n71,
    CNS18_y AS	n72,	
    CNS19_y AS	n81,
    CNS20_y AS	n92,
    ST_AsText(geom) AS geom
    
    FROM blocks_lehd_2011_us
    WHERE  
    geom && 
    ST_MakeEnvelope("
    
    bb2<-paste(bb[1,1],",",bb[2,1],",",bb[1,2],",",bb[2,2],",","4283")
    q_tsig <- paste(select,bb2,")")
    
    #Pull SQL result
    tsig_result <- dbGetQuery(con, q_tsig)
    tsig_result$ID<-row.names.data.frame(tsig_result)
    
    #tsig_result1<-tsig_result
    poly.spdf<-WKT2SpatialPolygonsDataFrame(tsig_result, 'geom','ID')
    #plot(poly.spdf)
    
    # #Convert node XY to spatial data frame
    # node.spdf <- SpatialPointsDataFrame(coords = tsig_result[,c(10,9)], data = tsig_result,
    #                                     proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
    # 
    
    #make point coordinate system same as cells
    proj4string(poly.spdf) <- CRS("+proj=longlat +ellps=WGS84") 
    proj4string(all_cells_cut) <- CRS("+proj=longlat +ellps=WGS84") 
    #all_cells_cut.proj <- spTransform(all_cells_cut, CRS( "+init=epsg:3347" ) ) #for area
    poly.spdf.proj <- spTransform( poly.spdf, CRS( "+init=epsg:3347" ) ) #for area
  
    #Use just the PageNumber column
    
    all_cells_cut_red<-all_cells_cut[,1]
    
    #spatial join nodes with cells
    #PageNumber = over(poly.spdf,all_cells_cut_red)
    
    #X1<-gIntersection(poly.spdf, all_cells_cut_red, byid=T,id=all_cells_cut_red$PageNumber)
    #PageNumber<-sapply(X1@polygons, function(x) x@PageNumber)
    
    X1<-raster::intersect(poly.spdf, all_cells_cut_red)
    PageNumber<-as.data.frame(X1$PageNumber)
    ID<-as.data.frame(X1$ID)
    Res<-cbind(PageNumber,ID)
    colnames(Res)<-c("PageNumber","ID")
    Res<-unique(Res)
    
    all_block_area<-data.frame(round(gArea(poly.spdf.proj, byid = TRUE)*0.000247105,2))
    ABA_ID<-as.data.frame(poly.spdf.proj$ID)
    all_block_area<-cbind(all_block_area,ABA_ID)
    colnames(all_block_area) <- c("BlockAcres","ID")
    
    X1.proj <- spTransform( X1, CRS( "+init=epsg:3347" ) ) #for area
    
    all_block_intersect<-as.data.frame(round(gArea(X1.proj, byid = TRUE)*0.000247105,2))
    all_block_intersect<-cbind(all_block_intersect,ID)
    colnames(all_block_intersect) <- c("Blockinter","ID")
    
    blockprop<-merge(all_block_intersect,all_block_area,"ID",all.x=T)
    
    tsig_result2<-merge(blockprop,tsig_result,"ID",all.x=T)
    tsig_result3<-cbind(tsig_result2[c("tot_jobs",	"n11",	"n21",	"n22", "n23",	"n31_33",	"n42",	"n44_45",	"n48_49",	"n51",	"n52",	"n53",	"n54",	"n55",	"n56",	"n61",	"n62",	"n71",	"n72",	"n81",	"n92","Blockinter", "BlockAcres")],Res[1])
    LEHD_FINAL=rbind(LEHD_FINAL,tsig_result3)
    
  }
})

print(time_LEHD)

#get cell area
all_cells_area<-data.frame(round(gArea(all_cells.proj, byid = TRUE)*0.000247105,2))
all_cells_area$PageNumber<- 1:nrow(all_cells_area)
colnames(all_cells_area)[1] <- "CellAcres"

LEHD_FINAL<-merge(all_cells_area,LEHD_FINAL,"PageNumber",all.x=T)

#Convert NAs to 0s
LEHD_FINAL[is.na(LEHD_FINAL)] <- 0

#proportion totals
LEHD_FINAL$acrepercent<-LEHD_FINAL$Blockinter/LEHD_FINAL$BlockAcres
#LEHD_FINAL$acrepercent[LEHD_FINAL$acrepercent > 1] <- 1
LEHD_FINAL$tot_jobs<-round(LEHD_FINAL$tot_jobs*LEHD_FINAL$acrepercent)
LEHD_FINAL$n11<-round(LEHD_FINAL$n11*LEHD_FINAL$acrepercent)
LEHD_FINAL$n21<-round(LEHD_FINAL$n21*LEHD_FINAL$acrepercent)
LEHD_FINAL$n22<-round(LEHD_FINAL$n22*LEHD_FINAL$acrepercent)
LEHD_FINAL$n23<-round(LEHD_FINAL$n23*LEHD_FINAL$acrepercent)
LEHD_FINAL$n31_33<-round(LEHD_FINAL$n31_33*LEHD_FINAL$acrepercent)
LEHD_FINAL$n42<-round(LEHD_FINAL$n42*LEHD_FINAL$acrepercent)
LEHD_FINAL$n44_45<-round(LEHD_FINAL$n44_45*LEHD_FINAL$acrepercent)
LEHD_FINAL$n48_49<-round(LEHD_FINAL$n48_49*LEHD_FINAL$acrepercent)
LEHD_FINAL$n51<-round(LEHD_FINAL$n51*LEHD_FINAL$acrepercent)
LEHD_FINAL$n52<-round(LEHD_FINAL$n52*LEHD_FINAL$acrepercent)
LEHD_FINAL$n53<-round(LEHD_FINAL$n53*LEHD_FINAL$acrepercent)
LEHD_FINAL$n54<-round(LEHD_FINAL$n54*LEHD_FINAL$acrepercent)
LEHD_FINAL$n55<-round(LEHD_FINAL$n55*LEHD_FINAL$acrepercent)
LEHD_FINAL$n56<-round(LEHD_FINAL$n56*LEHD_FINAL$acrepercent)
LEHD_FINAL$n61<-round(LEHD_FINAL$n61*LEHD_FINAL$acrepercent)
LEHD_FINAL$n62<-round(LEHD_FINAL$n62*LEHD_FINAL$acrepercent)
LEHD_FINAL$n71<-round(LEHD_FINAL$n71*LEHD_FINAL$acrepercent)
LEHD_FINAL$n72<-round(LEHD_FINAL$n71*LEHD_FINAL$acrepercent)
LEHD_FINAL$n81<-round(LEHD_FINAL$n81*LEHD_FINAL$acrepercent)
LEHD_FINAL$n92<-round(LEHD_FINAL$n92*LEHD_FINAL$acrepercent)


#Agg varaibles
LEHD_FINAL.agg<-ddply(LEHD_FINAL,~PageNumber,summarise, 
                      tot_jobs=sum(tot_jobs),
                      n11=sum(n11),
                      n21=sum(n21),
                      n22=sum(n22),
                      n23=sum(n23),
                      n31_33=sum(n31_33),
                      n42=sum(n42),
                      n44_45=sum(n44_45),
                      n48_49=sum(n48_49),
                      n51=sum(n51),
                      n52=sum(n52),
                      n53=sum(n53),
                      n54=sum(n54),
                      n55=sum(n55),
                      n56=sum(n56),
                      n61=sum(n61),
                      n62=sum(n62),
                      n71=sum(n71),
                      n72=sum(n72),
                      n81=sum(n81),
                      n92=sum(n92),
                      CellAcres=mean(CellAcres))


##Classify jobs
office_3<-as.data.frame(LEHD_FINAL.agg$n51+LEHD_FINAL.agg$n52+LEHD_FINAL.agg$n53+LEHD_FINAL.agg$n55+LEHD_FINAL.agg$n92)
colnames(office_3)<-"office_3"
retail_3<-as.data.frame(LEHD_FINAL.agg$n44_45+LEHD_FINAL.agg$n72)
colnames(retail_3)<-"retail_3"
service_3<-as.data.frame(LEHD_FINAL.agg$n54+LEHD_FINAL.agg$n56+LEHD_FINAL.agg$n61+LEHD_FINAL.agg$n62+LEHD_FINAL.agg$n71+LEHD_FINAL.agg$n81)
colnames(service_3)<-"service_3"

office_5<-as.data.frame(LEHD_FINAL.agg$n51+LEHD_FINAL.agg$n52+LEHD_FINAL.agg$n53+LEHD_FINAL.agg$n92)
colnames(office_5)<-"office_5"
retail_5<-as.data.frame(LEHD_FINAL.agg$n44_45)
colnames(retail_5)<-"retail_5"
service_5<-as.data.frame(LEHD_FINAL.agg$n54+LEHD_FINAL.agg$n56+LEHD_FINAL.agg$n61+LEHD_FINAL.agg$n62+LEHD_FINAL.agg$n81)
colnames(service_5)<-"service_5"
industrial_5<-as.data.frame(LEHD_FINAL.agg$n11+LEHD_FINAL.agg$n21+LEHD_FINAL.agg$n22+LEHD_FINAL.agg$n23+LEHD_FINAL.agg$n31_33+LEHD_FINAL.agg$n42+LEHD_FINAL.agg$n48_49)
colnames(industrial_5)<-"industrial_5"
entertainment_5<-as.data.frame(LEHD_FINAL.agg$n71+LEHD_FINAL.agg$n72)
colnames(entertainment_5)<-"entertainment_5"

office_8<-as.data.frame(LEHD_FINAL.agg$n51+LEHD_FINAL.agg$n52+LEHD_FINAL.agg$n53+LEHD_FINAL.agg$n55)
colnames(office_8)<-"office_8"
retail_8<-as.data.frame(LEHD_FINAL.agg$n44_45)
colnames(retail_8)<-"retail_8"
service_8<-as.data.frame(LEHD_FINAL.agg$n54+LEHD_FINAL.agg$n56+LEHD_FINAL.agg$n81)
colnames(service_8)<-"service_8"
industrial_8<-as.data.frame(LEHD_FINAL.agg$n11+LEHD_FINAL.agg$n21+LEHD_FINAL.agg$n22+LEHD_FINAL.agg$n23+LEHD_FINAL.agg$n31_33+LEHD_FINAL.agg$n42+LEHD_FINAL.agg$n48_49)
colnames(industrial_8)<-"industrial_8"
entertainment_8<-as.data.frame(LEHD_FINAL.agg$n71+LEHD_FINAL.agg$n72)
colnames(entertainment_8)<-"entertainment_8"
education_8<-as.data.frame(LEHD_FINAL.agg$n61)
colnames(education_8)<-"education_8"
healthcare_8<-as.data.frame(LEHD_FINAL.agg$n62)
colnames(healthcare_8)<-"healthcare_8"
pubadmin_8<-as.data.frame(LEHD_FINAL.agg$n92)
colnames(pubadmin_8)<-"pubadmin_8"

##prep entopy
office_3e<-as.data.frame(office_3/sum(office_3)*log(office_3/sum(office_3)))
office_3e[is.na(office_3e)] <- 0
retail_3e<-as.data.frame(retail_3/sum(retail_3)*log(retail_3/sum(retail_3)))
retail_3e[is.na(retail_3e)] <- 0
service_3e<-as.data.frame(service_3/sum(service_3)*log(service_3/sum(service_3)))
service_3e[is.na(service_3e)] <- 0

office_5e<-as.data.frame(office_5/sum(office_5)*log(office_5/sum(office_5)))
office_5e[is.na(office_5e)] <- 0
retail_5e<-as.data.frame(retail_5/sum(retail_5)*log(retail_5/sum(retail_5)))
retail_5e[is.na(retail_5e)] <- 0
service_5e<-as.data.frame(service_5/sum(service_5)*log(service_5/sum(service_5)))
service_5e[is.na(service_5e)] <- 0
industrial_5e<-as.data.frame(industrial_5/sum(industrial_5)*log(industrial_5/sum(industrial_5)))
industrial_5e[is.na(industrial_5e)] <- 0  
entertainment_5e<-as.data.frame(entertainment_5/sum(entertainment_5)*log(entertainment_5/sum(entertainment_5)))
entertainment_5e[is.na(entertainment_5e)] <- 0 

office_8e<-as.data.frame(office_8/sum(office_8)*log(office_8/sum(office_8)))
office_8e[is.na(office_8e)] <- 0
retail_8e<-as.data.frame(retail_8/sum(retail_8)*log(retail_8/sum(retail_8)))
retail_8e[is.na(retail_8e)] <- 0
service_8e<-as.data.frame(service_8/sum(service_8)*log(service_8/sum(service_8)))
service_8e[is.na(service_8e)] <- 0
industrial_8e<-as.data.frame(industrial_8/sum(industrial_8)*log(industrial_8/sum(industrial_8)))
industrial_8e[is.na(industrial_8e)] <- 0  
entertainment_8e<-as.data.frame(entertainment_8/sum(entertainment_8)*log(entertainment_8/sum(entertainment_8)))
entertainment_8e[is.na(entertainment_8e)] <- 0 
education_8e<-as.data.frame(education_8/sum(education_8)*log(education_8/sum(education_8)))
education_8e[is.na(education_8e)] <- 0 
healthcare_8e<-as.data.frame(healthcare_8/sum(healthcare_8)*log(healthcare_8/sum(healthcare_8)))
healthcare_8e[is.na(healthcare_8e)] <- 0 
pubadmin_8e<-as.data.frame(pubadmin_8/sum(pubadmin_8)*log(pubadmin_8/sum(pubadmin_8)))
pubadmin_8e[is.na(pubadmin_8e)] <- 0 

#Calc Employment variables
emp.temp<-LEHD_FINAL.agg[1]
emp.temp$EmpDens<-LEHD_FINAL.agg$tot_jobs/LEHD_FINAL.agg$CellAcres
emp.temp$EmpRetDens<-retail_5$retail_5/LEHD_FINAL.agg$CellAcres
emp.temp$Emp<-LEHD_FINAL.agg$tot_jobs
emp.temp$EmpRet<-retail_5$retail_5

#Retail, office and serice jobs entopy
#emp.temp$Entr_3Classes<-(-1)*(LEHD_FINAL.agg$retail_jobs/sum(LEHD_FINAL.agg$retail_jobs)*log(LEHD_FINAL.agg$retail_jobs/sum(LEHD_FINAL.agg$retail_jobs))+
#                              LEHD_FINAL.agg$office_jobs/sum(LEHD_FINAL.agg$office_jobs)*log(LEHD_FINAL.agg$office_jobs/sum(LEHD_FINAL.agg$office_jobs))+
#                              LEHD_FINAL.agg$service_jobs/sum(LEHD_FINAL.agg$office_jobs)*log(LEHD_FINAL.agg$service_jobs/sum(LEHD_FINAL.agg$office_jobs))/log(3))


# #Retail, office and serice jobs entopy
# ent.temp<-cbind(ret,off,srv,agr,min,ind)
# colnames(ent.temp) <- c("ret","off","srv","agr","min","ind")  

#3-class entropy
emp.temp$Entr_3C<-(-1)*((retail_3e$retail_3+
                           office_3e$office_3+
                           service_3e$service_3)/log(3))
#5-class entropy
emp.temp$Entr_5C<-(-1)*((retail_5e$retail_5+
                           office_5e$office_5+
                           service_5e$service_5+
                           industrial_5e$industrial_5+
                           entertainment_5e$entertainment_5)/log(5))
#8-class entropy
emp.temp$Entr_8C<-(-1)*((retail_8e$retail_8+
                           office_8e$office_8+
                           service_8e$service_8+
                           industrial_8e$industrial_8+
                           entertainment_8e$entertainment_8+
                           education_8e$education_8+
                           healthcare_8e$healthcare_8+
                           pubadmin_8e$pubadmin_8)/log(8))



###Get census Data
CENSUS_FINAL<- data.frame(pop=numeric(),
                          housing=numeric(),
                          lat=numeric(), 
                          lon=numeric(),
                          PageNumber=integer(),
                          BlockAcres=numeric(),
                          Blockinter=numeric(),
                          stringsAsFactors=FALSE) 
time_census<-system.time({
  #Chunk grid cells 
  Splitrows=100
  splitnumb_raw<-(nrow(all_cells)/Splitrows)
  if(nrow(all_cells)<Splitrows){
    Splitrows<-nrow(all_cells)
    splitnumb_l<-1
    leftover<-0
  }else{
    splitnumb_l<-floor(nrow(all_cells)/Splitrows)
    leftover<-splitnumb_raw-floor(splitnumb_l)
  }
  
  #Loops through cell chunks
  for(i in 1:splitnumb_l) {
    if(i==splitnumb_l){
      X1=(i-1)*Splitrows
      X2=Splitrows*i+Splitrows*leftover
    }else{i
      X1=(i-1)*Splitrows
      X2=Splitrows*i
    }
    
    #ith grid cell (row) bounding box and corners
    all_cells <-spTransform(all_cells, CRS("+proj=longlat +datum=WGS84"))
    all_cells_cut<-all_cells[X1:X2,]
    bb<-bbox(all_cells_cut)
    
    ##Find intersections & Cul da sacks
    select<-"SELECT 
    pop10                   AS  pop,
    housing10 							AS	housing,
    ST_AsText(geom) AS geom
    FROM blocks_census2010_us
    WHERE  
    geom && 
    ST_MakeEnvelope("
    
    bb2<-paste(bb[1,1],",",bb[2,1],",",bb[1,2],",",bb[2,2],",","4283")
    q_tsig <- paste(select,bb2,")")
    
    #Pull SQL result
    tsig_result <- dbGetQuery(con, q_tsig)
    tsig_result$ID<-row.names.data.frame(tsig_result)
    
    
    #tsig_result1<-tsig_result
    poly.spdf<-WKT2SpatialPolygonsDataFrame(tsig_result, 'geom','ID')
    #plot(test.spdf)
    
    # #Convert node XY to spatial data frame
    # node.spdf <- SpatialPointsDataFrame(coords = tsig_result[,c(4,3)], data = tsig_result,
    #                                     proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
    # 
    #make point coordinate system same as cells
    proj4string(poly.spdf) <- CRS("+proj=longlat +ellps=WGS84") 
    proj4string(all_cells_cut) <- CRS("+proj=longlat +ellps=WGS84") 
    poly.spdf.proj <- spTransform( poly.spdf, CRS( "+init=epsg:3347" ) ) #for area
    
    #Use just the PageNumber column
    all_cells_cut_red<-all_cells_cut[,1]
    
    #spatial join nodes with cells
    #PageNumber = over(poly.spdf,all_cells_cut_red)
    
    X1<-raster::intersect(poly.spdf, all_cells_cut_red)
    PageNumber<-as.data.frame(X1$PageNumber)
    ID<-as.data.frame(X1$ID)
    Res<-cbind(PageNumber,ID)
    colnames(Res)<-c("PageNumber","ID")
    Res<-unique(Res)
    
    all_block_area<-data.frame(round(gArea(poly.spdf.proj, byid = TRUE)*0.000247105,2))
    ABA_ID<-as.data.frame(poly.spdf.proj$ID)
    all_block_area<-cbind(all_block_area,ABA_ID)
    colnames(all_block_area) <- c("BlockAcres","ID")
    
    X1.proj <- spTransform( X1, CRS( "+init=epsg:3347" ) ) #for area
    
    all_block_intersect<-as.data.frame(round(gArea(X1.proj, byid = TRUE)*0.000247105,2))
    all_block_intersect<-cbind(all_block_intersect,ID)
    colnames(all_block_intersect) <- c("Blockinter","ID")
    
    blockprop<-merge(all_block_intersect,all_block_area,"ID",all.x=T)
    
    tsig_result2<-merge(blockprop,tsig_result,"ID",all.x=T)
    tsig_result3<-cbind(tsig_result2[c("pop","housing","Blockinter", "BlockAcres")],Res[1])
    CENSUS_FINAL=rbind(CENSUS_FINAL,tsig_result3)
    
  }
})

print(time_census)


#get cell area
all_cells_area<-data.frame(round(gArea(all_cells.proj, byid = TRUE)*0.000247105,2))
all_cells_area$PageNumber<- 1:nrow(all_cells_area)
colnames(all_cells_area)[1] <- "CellAcres"

CENSUS_FINAL<-merge(CENSUS_FINAL,all_cells_area,"PageNumber")

#Convert NAs to 0s
CENSUS_FINAL[is.na(CENSUS_FINAL)] <- 0

#proportion totals
CENSUS_FINAL$acrepercent<-CENSUS_FINAL$Blockinter/CENSUS_FINAL$BlockAcres
#CENSUS_FINAL$acrepercent[CENSUS_FINAL$acrepercent > 1] <- 1
CENSUS_FINAL$pop<-round(CENSUS_FINAL$pop*CENSUS_FINAL$acrepercent)
CENSUS_FINAL$housing<-round(CENSUS_FINAL$housing*CENSUS_FINAL$acrepercent)
#Agg varaibles
CENSUS_FINAL.agg<-ddply(CENSUS_FINAL,~PageNumber,summarise, 
                        pop=sum(pop),
                        housing=sum(housing),
                        CellAcres=mean(CellAcres))

blocks_count<-as.data.frame(table(CENSUS_FINAL$PageNumber))
colnames(blocks_count)<-c("PageNumber","Blocks")

#Calc census variables
cen.temp<-CENSUS_FINAL.agg[1]
cen.temp$PopDns<-CENSUS_FINAL.agg$pop/CENSUS_FINAL.agg$CellAcres
cen.temp$HUDns<-CENSUS_FINAL.agg$housing/CENSUS_FINAL.agg$CellAcres
cen.temp$Pop<-CENSUS_FINAL.agg$pop
cen.temp$HU<-CENSUS_FINAL.agg$housing
cen.temp$Blocks<-blocks_count$Blocks
cen.temp$AreaAcre<-CENSUS_FINAL.agg$CellAcres

##MERGE final data sets
BASEcells<-as.data.frame(all_cells.df[1])
BASEcells<-merge(BASEcells,cen.temp,"PageNumber",all.x=T)
BASEcells<-merge(BASEcells,emp.temp,"PageNumber",all.x=T)
BASEcells<-merge(BASEcells,int.temp,"PageNumber",all.x=T)


##clean 2nd to lasttime
BASEcells[is.na(BASEcells)] <- 0

#Calc last BE variables
BASEcells$ActvDns<-(BASEcells$Pop+BASEcells$Emp)/BASEcells$AreaAcre
RSE<-as.data.frame(retail_5$retail_5+service_5$service_5+entertainment_5$entertainment_5)
colnames(RSE)<-"RSE"
BASEcells$RSEDns<-(RSE$RSE)/BASEcells$AreaAcre
pop1<-as.data.frame(BASEcells$Pop)
colnames(pop1)<-"pop1"
pop1$pop1[pop1$pop1==0]<-1
BASEcells$EPB1<-BASEcells$Emp/pop1$pop1 
alpha<-sum(BASEcells$Emp)/sum(pop1$pop1)
BASEcells$EPB1_n<-(BASEcells$Emp/pop1$pop1)/alpha
BASEcells$EPB2<-retail_5$retail_5/pop1$pop1
beta<-sum(retail_5)/sum(pop1$pop1)
BASEcells$EPB2_n<-(retail_5$retail_5/pop1$pop1)/beta
BASEcells$EPB3<-RSE$RSE/pop1$pop1
gamma<-beta<-sum(RSE)/sum(pop1$pop1)
BASEcells$EPB3_n<-(RSE$RSE/pop1$pop1)/(RSE$RSE)

##last clean
BASEcells[is.na(BASEcells)] <- 0
BASEcells<-as.data.frame(sapply(BASEcells,as.numeric))
#BASEcells<-do.call(data.frame,sapply(BASEcells, function(x) replace(x, is.infinite(x),NA)))

#Join BE data to grid cells 
#all_cells@data = data.frame(all_cells@data, BASEcells[match(all_cells@data[,'PageNumber'], BASEcells[,'PageNumber']),])

#} #end county loop if looping

##-------------Calc netwotk distance to transit----------------##
if (file.exists(paste0("temp/",studyarea,".csv"))) file.remove(paste0("temp/",studyarea,".csv"))

TRANSIT_FINAL<- data.frame(
  originID=numeric(), 
  o_x=numeric(), 
  o_y=numeric(), 
  RailDIS_n=numeric(),
  rail_x=numeric(), 
  rail_y=numeric(),
  BusDIS_n=numeric(),
  bus_x=numeric(), 
  bus_y=numeric(),
  stringsAsFactors=FALSE)

#Chunk grid cells 
Splitrows=1000
splitnumb_raw<-(nrow(origins)/Splitrows)
if(nrow(origins)<Splitrows){
  Splitrows<-nrow(origins)
  splitnumb_l<-1
  leftover<-0
}else{
  splitnumb_l<-floor(nrow(origins)/Splitrows)
  leftover<-splitnumb_raw-floor(splitnumb_l)
}

#Loops through cell chunks
#for(j in 1:3) {
  for(j in 1:splitnumb_l) {
    
  if(j==splitnumb_l){
    X1=(j-1)*Splitrows
    X2=Splitrows*j+Splitrows*leftover
  }else{j
    X1=(j-1)*Splitrows
    X2=Splitrows*j
  }  
  
  num.cores <- detectCores() ## detect how many cores on your machine
  cl <- makeCluster(num.cores, type = "SOCK")
  registerDoParallel(cl)
  getDoParWorkers() ## check if all cores are working
  clusterEvalQ(cl,library(RPostgreSQL))
  clusterEvalQ(cl,library(DBI))
  clusterEvalQ(cl,library(RSQLite))
  clusterEvalQ(cl,library(sp))
  clusterEvalQ(cl,library(rgeos))
  
  origins_cut<-origins[X1:X2,]
  
  #Read in bike station locations
  timeTR<-system.time({
    OO<-nrow(origins_cut)

    TRANSIT_FINAL=foreach(i=1:OO,.combine = rbind,.errorhandling='remove') %dopar% {

      drv <- dbDriver("PostgreSQL")
      con <- dbConnect(drv, user = "postgres", dbname = "osm", host = "localhost")
      con_rt <- dbConnect(drv, user = "postgres", dbname = "opbeumDB", host = "localhost")
      
      origins_bb<-bbox(origins_cut[i,])
      #get origin ID and YX
      originID<-origins_cut$PageNumber[i] #change this to the ID column for the loaded shape'
      bbl<-origins_bb*.9995  
      bbh<-origins_bb*1.0005  
      O_xy<-coordinates(origins_cut[i,])
      
      #transform to planar coordinates
      origins.proj <- spTransform( origins_cut[i,], CRS( "+init=epsg:3347" ) ) 
      
      ##Find destinations in area
      select<-"SELECT 
      route_type,
      ST_Y(geom) AS lat, 
      ST_X(geom) AS lon 
      FROM gtfs_stops_us
      WHERE  
      geom && 
      ST_MakeEnvelope("
      
      bb2<-paste(bbl[1,1],",",bbl[2,1],",",bbh[1,2],",",bbh[2,2],",","4283")
      q_tsig <- paste(select,bb2,")")
      
      #Pull SQL result
      tsig_result <- dbGetQuery(con, q_tsig)
      
      #get only bus stations
      tsig_result_bus <- tsig_result[ which(tsig_result$route_type==3), ]
      
      #get only rail stations
      tsig_result_rail <- tsig_result[ which(tsig_result$route_type<3), ]
      
      ##
      ##get ORIGIN OSM NODE
      ##
      
      bbl<-origins_bb*.9995  
      bbh<-origins_bb*1.0005 
      
      ##Find closest OSM node or the Origin
      select<-"SELECT 
      osm_id,
      ST_Y(ST_Transform(geom_vertex,4326)) AS lat, 
      ST_X(ST_Transform(geom_vertex,4326)) AS lon 
      FROM at_2po_vertex
      WHERE  
      geom_vertex && 
      ST_MakeEnvelope("
      
      bb2<-paste(bbl[1],",",bbl[2],",",bbh[1],",",bbh[2],",","4283")
      q_tsig <- paste(select,bb2,")")
      
      #Pull SQL result
      tsig_result <- dbGetQuery(con_rt, q_tsig)
      
      #Convert node XY to spatial data frame
      node.spdf <- SpatialPointsDataFrame(coords = tsig_result[,c(3,2)], data = tsig_result,
                                          proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
      #transform to planor coordinates
      node.spdf.proj <- spTransform( node.spdf, CRS( "+init=epsg:3347" ) ) 
      
      #get closest osm node
      near_node<-gDistance(origins.proj, spgeom2=node.spdf.proj, byid=T, hausdorff=FALSE, densifyFrac = NULL)
      node_dist<-min(near_node)
      node_id<-node.spdf$osm_id[which(near_node==min(near_node))]
      
      #get the 'source' node for routing
      select<-"SELECT source
      FROM at_2po_4pgr
      WHERE osm_source_id =" 
      q_node_origin <- paste(select,node_id)
      
      qnode_result_origin <- dbGetQuery(con_rt, q_node_origin)[1,]
      
      #********************
      #Calc RAIL diatance**
      #********************
      if(nrow(tsig_result_rail)>0){
        
        transit.spdf_rail <- SpatialPointsDataFrame(coords = tsig_result_rail[,c(3,2)], data = tsig_result_rail,
                                                    proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
        
        #transform to planar coordinates
        transit.spdf_rail.proj <- spTransform( transit.spdf_rail, CRS( "+init=epsg:3347" ) ) 
        
        #get closest rail location
        near_rail<-gDistance(origins.proj, spgeom2=transit.spdf_rail.proj, byid=T, hausdorff=FALSE, densifyFrac = NULL)
        rail_dist<-min(near_rail)
        rail_xy<-unique(coordinates(transit.spdf_rail[which(near_rail==min(near_rail)),]))
        rail_row<-which(near_rail==min(near_rail))[1]
        
        ##
        ##get RAIL OSM NODE
        ##
        
        bbl<-rail_xy*.9995  
        bbh<-rail_xy*1.0005  
        
        ##Find closest OSM node
        select<-"SELECT 
        osm_id,
        ST_Y(ST_Transform(geom_vertex,4326)) AS lat, 
        ST_X(ST_Transform(geom_vertex,4326)) AS lon 
        FROM at_2po_vertex
        WHERE  
        geom_vertex && 
        ST_MakeEnvelope("
        
        bb2<-paste(bbl[1],",",bbl[2],",",bbh[1],",",bbh[2],",","4283")
        q_tsig <- paste(select,bb2,")")
        
        
        #Pull SQL result
        tsig_result <- dbGetQuery(con_rt, q_tsig)
        
        if(nrow(tsig_result)>0){ 
          
          #Convert node XY to spatial data frame
          node.spdf <- SpatialPointsDataFrame(coords = tsig_result[,c(3,2)], data = tsig_result,
                                              proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
          #transform to planor coordinates
          node.spdf.proj <- spTransform( node.spdf, CRS( "+init=epsg:3347" ) ) 
          
          #get closest osm node
          near_node_rail<-gDistance(transit.spdf_rail.proj[rail_row,], spgeom2=node.spdf.proj, byid=T, hausdorff=FALSE, densifyFrac = NULL)
          node_rail_dist<-min(near_node_rail)
          node_rail_id<-node.spdf$osm_id[which(near_node_rail==min(near_node_rail))]
          
          #get the 'source' node for routing
          select<-"SELECT target
          FROM at_2po_4pgr
          WHERE osm_target_id =" 
          q_node_rail <- paste(select,node_rail_id)
          
          qnode_result_rail <- dbGetQuery(con_rt, q_node_rail)[1,]
          
          
          #----------------------
          
          ##get OD dist to rail
          
          q1<-   "SELECT seq, id1 AS node, id2 AS edge,km AS cost
          FROM pgr_astar('
          SELECT id AS id,
          source,
          target,
          cost,
          x1, y1, x2, y2
          FROM at_2po_4pgr as r,
          (SELECT ST_Expand(ST_Extent(geom_way),0.01) as box  FROM at_2po_4pgr as l1  
          WHERE l1.source ="
          o_node<-qnode_result_origin
          q2<-"OR l1.target =" 
          d_node<-qnode_result_rail
          q3<-") as box
          WHERE r.geom_way && box.box',"
          q4<-","
          q5<-", false, false)as r INNER JOIN at_2po_4pgr as g ON r.id2 = g.id ;"
          
          q_railD <- paste0(q1,o_node,q2,d_node,q3,o_node,q4,d_node,q5)
          
          d_result_rail <- dbGetQuery(con_rt, q_railD)
          
          RailDIS_n=sum(d_result_rail$cost)*0.621371
          
          #----------------------------  
        } else{
          RailDIS_n<-9999
          rail_xy<-matrix(0, 1, 2)
        }
        
      } else{
        RailDIS_n<-9999
        rail_xy<-matrix(0, 1, 2)
      }
      
      #********************
      #Calc BUS diatance**
      #********************
      if(nrow(tsig_result_bus)>0){
        
        #Convert node XY to spatial data frame
        transit.spdf_bus <- SpatialPointsDataFrame(coords = tsig_result_bus[,c(3,2)], data = tsig_result_bus,
                                                   proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))    
        #transform to planar coordinates
        transit.spdf_bus.proj <- spTransform( transit.spdf_bus, CRS( "+init=epsg:3347" ) ) 
        
        #get closest bus location
        near_bus<-gDistance(origins.proj, spgeom2=transit.spdf_bus.proj, byid=T, hausdorff=FALSE, densifyFrac = NULL)
        bus_dist<-min(near_bus)
        bus_xy<-unique(coordinates(transit.spdf_bus[which(near_bus==min(near_bus)),]))
        bus_row<-which(near_bus==min(near_bus))[1]
        
        ##
        ##get ORIGIN OSM NODE
        ##
        
        ##Find closest OSM node or the Origin
        select<-"SELECT 
        osm_id,
        ST_Y(ST_Transform(geom_vertex,4326)) AS lat, 
        ST_X(ST_Transform(geom_vertex,4326)) AS lon 
        FROM at_2po_vertex
        WHERE  
        geom_vertex && 
        ST_MakeEnvelope("
        
        bb2<-paste(bbl[1],",",bbl[2],",",bbh[1],",",bbh[2],",","4283")
        q_tsig <- paste(select,bb2,")")
        
        #Pull SQL result
        tsig_result_bus <- dbGetQuery(con_rt, q_tsig)
        
        if(nrow(tsig_result_bus)>0){  
          
          #Convert node XY to spatial data frame
          node.spdf <- SpatialPointsDataFrame(coords = tsig_result_bus[,c(3,2)], data = tsig_result_bus,
                                              proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
          #transform to planor coordinates
          node.spdf.proj <- spTransform( node.spdf, CRS( "+init=epsg:3347" ) ) 
          
          #get closest osm node
          near_node<-gDistance(origins.proj, spgeom2=node.spdf.proj, byid=T, hausdorff=FALSE, densifyFrac = NULL)
          node_dist<-min(near_node)
          node_id<-node.spdf$osm_id[which(near_node==min(near_node))]
          
          
          #----------------------
          
          ##
          ##get bus OSM NODE
          ##
          
          bbl<-bus_xy*.9995  
          bbh<-bus_xy*1.0005  
          
          ##Find closest OSM node
          select<-"SELECT 
          osm_id,
          ST_Y(ST_Transform(geom_vertex,4326)) AS lat, 
          ST_X(ST_Transform(geom_vertex,4326)) AS lon 
          FROM at_2po_vertex
          WHERE  
          geom_vertex && 
          ST_MakeEnvelope("
          
          bb2<-paste(bbl[1],",",bbl[2],",",bbh[1],",",bbh[2],",","4283")
          q_tsig <- paste(select,bb2,")")
          
          #Pull SQL result
          tsig_result <- dbGetQuery(con_rt, q_tsig)
          
          #Convert node XY to spatial data frame
          node.spdf <- SpatialPointsDataFrame(coords = tsig_result[,c(3,2)], data = tsig_result,
                                              proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
          #transform to planor coordinates
          node.spdf.proj <- spTransform( node.spdf, CRS( "+init=epsg:3347" ) ) 
          
          #get closest osm node
          near_node_bus<-gDistance(transit.spdf_bus.proj[bus_row,], spgeom2=node.spdf.proj, byid=T, hausdorff=FALSE, densifyFrac = NULL)
          node_bus_dist<-min(near_node_bus)
          node_bus_id<-tsig_result$osm_id[which(near_node_bus==min(near_node_bus))]
          
          #get the 'source' node for routing
          select<-"SELECT target
          FROM at_2po_4pgr
          WHERE osm_source_id =" 
          q_node_bus <- paste(select,node_bus_id)
          
          qnode_result_bus <- dbGetQuery(con_rt, q_node_bus)[1,]
          
          #----------------------------  
          
          ##get OD dist to Bus
          
          q1<-   "SELECT seq, id1 AS node, id2 AS edge,km AS cost
          FROM pgr_astar('
          SELECT id AS id,
          source,
          target,
          cost,
          x1, y1, x2, y2
          FROM at_2po_4pgr as r,
          (SELECT ST_Expand(ST_Extent(geom_way),0.01) as box  FROM at_2po_4pgr as l1  
          WHERE l1.source ="
          o_node<-qnode_result_origin
          q2<-"OR l1.target =" 
          d_node<-qnode_result_bus
          q3<-") as box
          WHERE r.geom_way && box.box',"
          q4<-","
          q5<-", false, false)as r INNER JOIN at_2po_4pgr as g ON r.id2 = g.id ;"
          
          q_busD <- paste0(q1,o_node,q2,d_node,q3,o_node,q4,d_node,q5)
          
          d_result_bus <- dbGetQuery(con_rt, q_busD)
          
          BusDIS_n=sum(d_result_bus$cost)*0.621371
          
          #----------------------------  
          
        } else{
          BusDIS_n<-9999
          bus_xy<-matrix(0, 1, 2)
        }
        
      } else{
        BusDIS_n<-9999
        bus_xy<-matrix(0, 1, 2)
      }
      
      #end distance calcs  
      
      dbDisconnect(con)
      dbDisconnect(con_rt)
      
      #prep column names
      colnames(O_xy) <- c("o_x","o_y")
      colnames(rail_xy) <- c("rail_x","rail_y")
      colnames(bus_xy) <- c("bus_x","bus_y")
      
      #Merge all location results
      cbind.data.frame(originID,O_xy,RailDIS_n,rail_xy,BusDIS_n,bus_xy,row.names = NULL)
      
    }  
    
    write.table(TRANSIT_FINAL, paste0("temp/",studyarea,".csv"), row.names = F, col.names = F, append = T, sep=",",quote=F)
    
    stopCluster(cl)
  })
  print(timeTR)
  print(splitnumb_l-j)
}

TRANSIT_FINAL2 = fread(paste0("temp/",studyarea,".csv"), sep=",",header = F)
colnames(TRANSIT_FINAL2)<-c("PageNumber","o_x","o_y","RailDIS_n","rail_x","rail_y","BusDIS_n","bus_x","bus_y")
#write.table(TRANSIT_FINAL2, paste0(studyarea,"2",".csv"), row.names = F, col.names = T, append = F, sep=",",quote=F)

########################
## End Transit router ##
########################






#Merge transit distances with final dataset
BASEcells<-merge(BASEcells,TRANSIT_FINAL2,all.x=T)

#####CREATE STEVE'S NEW VARIABLE!
BASEcells$ACTIVITY<-BASEcells$Pop +BASEcells$Emp


#Join BE data to grid cells 
all_cells@data = data.frame(all_cells@data, BASEcells[match(all_cells@data[,'PageNumber'], BASEcells[,'PageNumber']),])

#Write outputs (cSV and Shapefile)
writeOGR(obj=all_cells, dsn="Output", layer="OPBEUM_RESULT_GRID", driver="ESRI Shapefile",overwrite_layer=T,check_exists=T) 
write.table(BASEcells, paste0("Output/","OPBEUM_RESULT_GRID.csv"), row.names = F, col.names = T, append = F, sep=",",quote=F)


###multimodal access add-on----

###------------------------###
###------Bing Geocoder-----###
###------------------------###

#list needed packages
list.of.packages <- c("rgdal","RCurl","RJSONIO","rgeos","maptools","broom","ggplot2","dplyr","rjson","chron","sp","leaflet","reshape","KernSmooth","htmlwidgets","data.table")

#load packages
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

#load packages
lapply(list.of.packages, require, character.only = TRUE)


#set Bing maps API key (https://msdn.microsoft.com/en-us/library/ff428642.aspx)
BingMapsKey<-"4SSYi1cN6ZqhHBFBGpxF~irqp0i-N366dGuqrTU6IWw~AvL7kNX1wh9dsO0GZa3QpkR1SDQKbwjc__cHl3sfL7vQ7Sfn_kJAp6cxSUmd0-XD"





#Walking Routes (this sets up the API query for walking travel)
walking <- function(origin,destination, BingMapsKey){
  require(RCurl)
  require(RJSONIO)
  u <- URLencode(paste0("http://dev.virtualearth.net/REST/V1/Routes/Walking?wp.0=",origin,"&wp.1=",destination,"&key=",BingMapsKey))
  d <- getURL(u)
  j <- RJSONIO::fromJSON(d,simplify = FALSE) 
  if (j$statusCode == 200) {
    j<-j
  }
  else {    
    j <- NA
  }
  j
}

#query Bing distance API and save results
relay_pairs_sub<-slice(relay_pairs,1:100)

mlat.df<-data.frame(matrix(0, nrow = nrow(relay_pairs_sub), ncol = 100,
                           dimnames = list(NULL, paste0("mlat", 1:100))))
mlon.df<-data.frame(matrix(0, nrow = nrow(relay_pairs_sub), ncol = 100,
                           dimnames = list(NULL, paste0("mlon", 1:100))))

for (a in 1:nrow(relay_pairs_sub)) {
  tryCatch({
    orig<-paste(relay_pairs_sub[a,2],relay_pairs_sub[a,3],sep=",")
    dest<-paste(relay_pairs_sub[a,5],relay_pairs_sub[a,6],sep=",")
    w<-walking(orig,dest,BingMapsKey)
    ml<-length(w$resourceSets[[1]]$resources[[1]]$routeLegs[[1]]$itineraryItems)
    mlat.df[a,1]<-relay_pairs_sub[a,2]
    mlon.df[a,1]<-relay_pairs_sub[a,3]
    mlat.df[a,ml+2]<-relay_pairs_sub[a,5]
    mlon.df[a,ml+2]<-relay_pairs_sub[a,6]
    for (b in 1:length(w$resourceSets[[1]]$resources[[1]]$routeLegs[[1]]$itineraryItems)) {
      mlat.df[a,b+1]<-w$resourceSets[[1]]$resources[[1]]$routeLegs[[1]]$itineraryItems[[b]]$maneuverPoint$coordinates[1]
      mlon.df[a,b+1]<-w$resourceSets[[1]]$resources[[1]]$routeLegs[[1]]$itineraryItems[[b]]$maneuverPoint$coordinates[2]
      
    }
    
  }, error=function(e){})
}


################
## END OPBEUM ##
################

##Add leaflet map

#load ppl stations
ppl_stations<-readOGR(dsn = "Input", layer = "purple_line_stops")
ppl_stations<-spTransform( ppl_stations,CRS("+proj=longlat +datum=WGS84")) 
#load  lines
ppl_line<-readOGR(dsn = "Input", layer = "purple_line")
ppl_line<-spTransform( ppl_line,CRS("+proj=longlat +datum=WGS84")) 

metro_line<-readOGR(dsn = "Input", layer = "Metro__Lines")
metro_line<-spTransform( metro_line,CRS("+proj=longlat +datum=WGS84")) 

metro_bus_line<-readOGR(dsn = "Input", layer = "Metro_Bus_Lines")
metro_bus_line<-spTransform( metro_bus_line,CRS("+proj=longlat +datum=WGS84")) 


#load grids
all_cells<-readOGR(dsn = "Output", layer = "OPBEUM_RESULT_GRID")


# required libraries
library(leaflet, quietly = T, warn.conflicts = F)

# start basemap
map <- leaflet() %>% 
  
  # add ocean basemap
  addProviderTiles(providers$Esri.OceanBasemap) %>%
  
  # add another layer with place names
  addProviderTiles(providers$Hydda.RoadsAndLabels, group = 'Place names') %>%
  
  
  # focus map in a certain area / zoom level
  setView(lng = -77.02450, lat = 38.99478, zoom = 13) %>%
  
  # add layers control
  addLayersControl(overlayGroups = c('Place names',
                                     'Purple Line Stations',
                                     'Purple Line',
                                     'WMATA Rail Lines',
                                     'WMATA Bus Lines',
                                     'Rail Distance',
                                     'Employment Density'),
                   options = layersControlOptions(collapsed = FALSE),
                   position = 'topright') %>%
  
  # list groups to hide on startup
  hideGroup(c('WMATA Bus Lines','WMATA Rail Lines','Employment Density'))

# add points
map <- map %>%
  addCircleMarkers(data = as.data.frame(ppl_stations@coords), ~coords.x1, ~coords.x2,
                   weight = 0.5,
                   col = 'black', 
                   fillColor = 'darkslategrey',
                   radius = 4, 
                   fillOpacity = 0.9, 
                   stroke = T, 
                   label = ~paste0('Point at: ', 
                                   as.character(round(coords.x1,3)), ', ', 
                                   as.character(round(coords.x2,3))), 
                   group = 'Purple Line Stations')
# add lines

map <- map %>%
  addPolylines(data = metro_line,
               weight = 3,
               color = 'green',
               opacity = .50,
               popup = 'WMATA Rail', 
               smoothFactor = 3,
               group = 'WMATA Rail Lines') 

map <- map %>%
  addPolylines(data = metro_bus_line,
               weight = 2,
               color = 'blue',
               opacity = .50,
               popup = 'WMATA Bus', 
               smoothFactor = 3,
               group = 'WMATA Bus Lines') 

map <- map %>%
addPolylines(data = ppl_line,
             weight = 3,
             color = 'purple',
             popup = 'Purple Line', 
             smoothFactor = 3,
             group = 'Purple Line') 
# add polygons
bins <- c(0, .25,.5, .75, Inf)
pal <- colorBin("YlOrRd", domain = all_cells$RlDIS_n, bins = bins)


# add polygons
map <- map %>%
  addPolygons(data=all_cells,
              weight = 1, 
              color = 'grey', 
              fillColor = ~pal(RlDIS_n),
              fill = T, 
              fillOpacity = 0.25, 
              stroke = T, 
              dashArray = c(5,5), 
              smoothFactor = 3,
              group = 'Rail Distance')

bins2 <- c(0, 1,2, 5, Inf)
pal2 <- colorBin("Blues", domain = all_cells$EmpDens, bins = bins2)

map <- map %>%
  addPolygons(data=all_cells,
              weight = 1,
              color = 'grey',
              fillColor = ~pal2(EmpDens),
              fill = T,
              fillOpacity = 0.25,
              stroke = T,
              dashArray = c(5,5),
              smoothFactor = 3,
              group = 'Employment Density')

# add more features
map <- map %>% 
  
  # add a map scalebar
  addScaleBar(position = 'topright') %>%
  
  # add measurement tool
  addMeasure(
    primaryLengthUnit = "kilometers",
    secondaryLengthUnit = 'miles', 
    primaryAreaUnit = "hectares",
    secondaryAreaUnit="acres", 
    position = 'topleft')

# show map    
map    

# # # save a stand-alone, interactive map as an html file
 library(htmlwidgets)
saveWidget(widget = map, file = 'map.html', selfcontained = T)

# # # save a snapshot as a png file
library(mapview)
mapshot(map, file = 'map.png')
map
