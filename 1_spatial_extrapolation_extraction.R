 ### script to calculate multiple overlap metrics between species and gaps and fishing hours

# load libraries ####
datapath="/Users/EcoCast/Dropbox"

library(raster)
library(glue)
library(sf)
library(tidyverse)
library(scales)
library(maps)
library(lubridate)
source(glue("{datapath}/OLE/github/OLE_Projects/utilities/load_libraries.R"))
to360 <- function(x) {x %% 360}
round_any = function(x, accuracy, f=round){f(x/ accuracy) * accuracy}

# define directories and paths ####
outdir="UserDefined";dir.create(outdir)
template=raster(glue("{datapath}/OLE/spatial_data/template.grd")) ## this is a user defined raster to assign fishing effort / gaps to
spdir="UserDefined"
template_1deg=template
res(template_1deg)=1
nums=seq(1,ncell(template_1deg))
template_1deg=setValues(template_1deg,nums)
template_1deg_coords=rasterToPoints(template_1deg) %>% 
  as.data.frame() %>% rename(cell_lat=y,cell_lon=x,cells=Sea.level.anomaly)

data=maps::map("world2",fill=T)
IDs <- sapply(strsplit(data$names, ":"), function(x) x[1])
wrld_simpl <- map2SpatialPolygons(data, IDs=IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))
wrld=SpatialPolygons(wrld_simpl@polygons,proj4string=wrld_simpl@proj4string) %>% 
  gBuffer(., byid=TRUE, width=0)

# load / locate species gaps fishing
## gaps ####
rawgaps=read.csv(glue("{datapath}/ole_gap_events_v20230313.csv")) ## these are already gaps

gaps2=rawgaps %>% mutate(date=as.Date(gap_start)) %>% ## restrict gaps to area of interest (AOI)
  mutate(off_lat_flag=case_when(off_lat>=10 & off_lat<60~"AOI",
                                T~"False")) %>% 
  mutate(on_lat_flag=case_when(on_lat>=10 &on_lat<60~"AOI",
                               T~"False")) %>% 
  mutate(off_lon_flag=case_when(off_lon>=(-179.9) & off_lon<(-100.1)~"AOI",
                                T~"False")) %>% 
  mutate(on_lon_flag=case_when(on_lon>=(-179.9) & on_lon<(-100.1)~"AOI",
                               T~"False")) %>% 
  filter((off_lat_flag=="AOI" & off_lon_flag=="AOI") |(on_lat_flag=="AOI" & on_lon_flag=="AOI") ) %>% 
  mutate(lat=off_lat,lon=off_lon) %>%
  mutate(lat_end=on_lat,lon_end=on_lon) %>%
  mutate(lon=to360(lon)) %>% 
  mutate(lon_end=to360(lon_end)) %>% 
  mutate(date=as.Date(gap_start)) %>% 
  dplyr::select(gap_start,gap_end,gap_hours,lon,lat,lon_end,lat_end,date,vessel_class,flag,gap_id,disabling_event)

gap_dates=unique(gaps2$date)

## Spatially extrapolate gaps ####
outdirsub=glue("{datapath}/UserDefined");dir.create(outdirsub)
outdirgaps=glue("{datapath}/UserDefined");dir.create(outdirgaps)
#24*14 = 336 hours in two weeks for two week filter

glimpse(gaps2)
gaps3=gaps2 %>% mutate(gap_start=as.POSIXct(gap_start)) %>% 
  mutate(gap_end=as.POSIXct(gap_end)) %>% dplyr::select(-date) 

ids=unique(gaps3$gap_id)

library(foreach)
library(doParallel, quietly = TRUE)
registerDoParallel(20)

    foreach(i=1:length(ids),
                    .export = c("gaps3","template_1deg","template_1deg_coords","ids","outdirgaps"),
    .combine=rbind,.packages = c("glue","tidyverse","raster","maps","sp"),.verbose=T) %dopar% {
  
     if(!file.exists(glue("{outdirgaps}/{ids[i]}.csv"))){
      
      dat=gaps3 %>% filter(gap_id==ids[i])
  ends=dat %>% dplyr::select(-c(lon,lat)) %>% rename(lon=lon_end,lat=lat_end) %>% 
    mutate(type="end")
  starts=dat %>% dplyr::select(-c(lon_end,lat_end)) %>% 
    mutate(type="start")
  dat2=rbind(starts,ends)
  
  xy <- dat2[,c(4,5)]
  
  spdf=xy %>% as.matrix() %>%  st_linestring() %>% st_sfc() %>% st_as_sf() %>% mutate(junk=1) %>% 
   as_Spatial()

  test=raster::extract(template_1deg,spdf) %>% unlist()
  ncels=length(test)
  dat3=dat[rep(seq_len(nrow(dat)), each = ncels), ] %>% 
    mutate(cells=test) #%>%  ## data frame with all template cells
  
  ## now dividing up time
  days=seq(as.Date(dat$gap_start),as.Date(dat$gap_end),by=1)
  
  percell=length(days)/ncels ## number of days in each cell
  mess=rep(percell,ncels) ## perpixel chunks of time
  smart.round <- function(x, digits = 0) {
    up <- 10 ^ digits
    x <- x * up
    y <- floor(x)
    indices <- tail(order(x-y), round(sum(x)) - sum(y))
    y[indices] <- y[indices] + 1
    y / up
  } ## handling odd numbers of days percell
  perpixelround=smart.round(mess)

  dat_df=as.data.frame(days)
  
  empty=list()
  for(ii in 1:length(test)){
    pixel=test[ii] 
    a=rep(pixel,perpixelround[ii])
    empty[[length(empty)+1]]=a
  }
  
  sequence=unlist(empty)
  dat_df2=dat_df %>% mutate(cells=sequence) %>% 
    left_join(.,template_1deg_coords,by="cells")%>% 
    left_join(dat3,by="cells") %>% as.data.frame()
  
   write.csv(dat_df2,glue("{outdirgaps}/{ids[i]}.csv"))
     }else{NULL}
  
            }


## appending individual gaps into table ####
    dir=outdirgaps
    outdir="UserDefined";dir.create(outdir)

    files=list.files(dir,full.names = T)
    
    trial=files %>%
      lapply(read_csv) %>%
      bind_rows
    
    write.csv(trial,glue("{outdir}/gaps_interpolated.csv"))
    
## extracting species data to interpolated gaps ####
  inter_gaps=read.csv(glue("{outdir}/gaps_interpolated.csv"))
gap_dates=unique(inter_gaps$days)
outdir_inter="UserDefined";dir.create(outdir_inter)
buffer_rad=(111.1*1000)/2 ## because gaps are interpolated at 1 degree, need a 1 degree mean of sdms, this is a radius

library(foreach)
library(doParallel, quietly = TRUE)
registerDoParallel(20)
system.time(print(
   foreach(i=1:length(gap_dates),.export = c("inter_gaps","spdir","outdir_inter","buffer_rad"),
          .combine=rbind,.packages = c("glue","tidyverse","raster","maps"),.verbose=T) %dopar% {
            print(gap_dates[i])
            get_date=gap_dates[i]
            if(!file.exists(glue("{outdir_inter}/gaps_species_overlap_{get_date}.csv"))){
              
              fullist=c(albacoretuna_TOPP= NA_real_,
                        `black-footedAlbatross_Dallas`= NA_real_,
                        blueShark_TOPP= NA_real_,
                        californiaSeaLion_TOPP= NA_real_,
                        blueWhale_TOPP= NA_real_,
                        elephantSeal_TOPP= NA_real_,
                        laysanAlbatross_Dallas= NA_real_,
                        leatherbackTurtle_TOPP= NA_real_,
                        makoShark_TOPP= NA_real_,
                        pacificBluefinTuna_TOPP= NA_real_,
                        salmonShark_TOPP= NA_real_,
                        sootyShearwater_TOPP= NA_real_
                        ,whiteShark_TOPP= NA_real_,
                        yellowfinTuna_TOPP= NA_real_)
              
              newmean=function(x){
                mean(x,na.rm=T)}
              
              ras1=list.files(spdir,pattern=".grd",full.names = T)  %>% 
                grep(get_date,.,value=T)
              if(length(ras1)==0){NULL}
              else {
                ras=ras1 %>% stack()
                names=list.files(spdir,pattern=".grd",full.names = F) %>% 
                  grep(get_date,.,value=T) %>% 
                  gsub(glue("_{get_date}.grd"),"",.)
                names(ras)=names
                
                extraction=raster::extract(ras,inter_gaps %>% filter(days==get_date) %>% dplyr::select(cell_lon,cell_lat),fun=newmean,na.rm=T,buffer=buffer_rad) %>% as.data.frame() # extract points from raster stack
                extraction_full=add_column(extraction, !!!fullist[setdiff(names(fullist),names(extraction))]) # add in any missing columns (in case a variable w.n. available for a specific day)
                
                master=cbind((inter_gaps %>% filter(days==get_date)),extraction_full)%>%
                  dplyr::select(c(flag,vessel_class,cell_lat,cell_lon,days,
                                  cells,gap_start,gap_end,gap_hours,lon,lat,lon_end,lat_end,
                                  gap_id,disabling_event,
                                  albacoretuna_TOPP,black.footedAlbatross_Dallas,blueShark_TOPP,blueWhale_TOPP,
                                  californiaSeaLion_TOPP,elephantSeal_TOPP,laysanAlbatross_Dallas,leatherbackTurtle_TOPP,
                                  makoShark_TOPP,pacificBluefinTuna_TOPP,salmonShark_TOPP,sootyShearwater_TOPP,
                                  whiteShark_TOPP,yellowfinTuna_TOPP)) # heather's code line
                rownames(master)=seq(1:nrow(master))
                write.csv(master,glue("{outdir_inter}/gaps_interpolated_species_overlap_{get_date}.csv")) 
              }
            }
          }
))

overlap_interpolated=list.files(outdir_inter,full.names = T) %>% 
  lapply(read_csv) %>%
  bind_rows
write.csv(overlap_interpolated,"UserDefined/species_gaps_interpolated.csv")

## observed fishery activity ####
dir=glue("{datapath}/UserDefined")
fish=list.files(dir,full.names = T) %>% 
  lapply(read_csv) %>%
  bind_rows

### new species/fishery overlap: extracto ####
fish_cleaned=fish %>% 
  mutate(lon_bin=to360(lon_bin)) %>% 
  filter(date<=as.Date("2022-12-31"))

fish_dates=unique(fish_cleaned$date)
outdir_inter="UserDefined";dir.create(outdir_inter)
# buffer_rad=(111.1*1000)/2 ## because gaps are interpolated at 1 degree, need a 1 degree mean of sdms, this is a radius

library(foreach)
library(doParallel, quietly = TRUE)
registerDoParallel(20)
foreach(i=1:length(fish_dates),.export = c("fish_cleaned","spdir","outdir_inter","buffer_rad"),
                  .combine=rbind,.packages = c("glue","tidyverse","raster","maps"),.verbose=T) %dopar% {
            print(fish_dates[i])
            get_date=fish_dates[i]
            if(!file.exists(glue("{outdir_inter}/fish_extracto_species_overlap_{get_date}.csv"))){
              
              fullist=c(albacoretuna_TOPP= NA_real_,
                        `black-footedAlbatross_Dallas`= NA_real_,
                        blueShark_TOPP= NA_real_,
                        californiaSeaLion_TOPP= NA_real_,
                        blueWhale_TOPP= NA_real_,
                        elephantSeal_TOPP= NA_real_,
                        laysanAlbatross_Dallas= NA_real_,
                        leatherbackTurtle_TOPP= NA_real_,
                        makoShark_TOPP= NA_real_,
                        pacificBluefinTuna_TOPP= NA_real_,
                        salmonShark_TOPP= NA_real_,
                        sootyShearwater_TOPP= NA_real_
                        ,whiteShark_TOPP= NA_real_,
                        yellowfinTuna_TOPP= NA_real_)
              
              newmean=function(x){
                mean(x,na.rm=T)}
              
              ras1=list.files(spdir,pattern=".grd",full.names = T)  %>% 
                grep(get_date,.,value=T)
              if(length(ras1)==0){NULL}
              else {
                ras=ras1 %>% raster::stack()
                names=list.files(spdir,pattern=".grd",full.names = F) %>% 
                  grep(get_date,.,value=T) %>% 
                  gsub(glue("_{get_date}.grd"),"",.)
                names(ras)=names
                
                extraction=raster::extract(ras,fish_cleaned %>% filter(date==get_date) %>% dplyr::select(lon_bin,lat_bin)) %>% as.data.frame() # extract points from raster stack
                extraction_full=add_column(extraction, !!!fullist[setdiff(names(fullist),names(extraction))]) # add in any missing columns (in case a variable w.n. available for a specific day)
                
                master=cbind((fish_cleaned %>% filter(date==get_date)),extraction_full)%>%
                  dplyr::select(c(flag,vessel_class,lat_bin,lon_bin,date,
                                  hours,fishing_hours,hours_in_gaps_over_12,hours_in_gaps_under_12,fishing_hours_in_gaps_over_12,fishing_hours_in_gaps_under_12,
                                   albacoretuna_TOPP,black.footedAlbatross_Dallas,blueShark_TOPP,blueWhale_TOPP,
                                  californiaSeaLion_TOPP,elephantSeal_TOPP,laysanAlbatross_Dallas,leatherbackTurtle_TOPP,
                                  makoShark_TOPP,pacificBluefinTuna_TOPP,salmonShark_TOPP,sootyShearwater_TOPP,
                                  whiteShark_TOPP,yellowfinTuna_TOPP)) # heather's code line
                rownames(master)=seq(1:nrow(master))
                write.csv(master,glue("{outdir_inter}/fish_extracto_species_overlap_{get_date}.csv")) 
              }
              
            }
          }

overlap_fish_sp=list.files(outdir_inter,full.names = T) %>% 
  lapply(read_csv) %>%
  bind_rows
write.csv(overlap_fish_sp,"/Users/EcoCast/Dropbox/OLE/overlap_OLE/overlap_03_14_23/overlap_aggregated_dfs/species_fishery_extracto.csv")

