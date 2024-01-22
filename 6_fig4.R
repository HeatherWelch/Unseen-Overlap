## RFMO analysis
# for all of these, only species and gears available in the RFMO databases

datapath="UserDefined"
sf::sf_use_s2(FALSE)
library(glue)
library(sf)
source(glue("{datapath}/Dropbox/OLE/github/OLE_Projects/utilities/load_libraries.R"))
outdir="UserDefined";dir.create(outdir)

# templates and extents ####
lat_lon=data.frame(label=c("40°N  ","20°N","220°","240°"),
                   x=c(208,208,220,240),
                   y=c(40,20,9,9)
) ## for drawing parallels and meridians

df=data.frame(xmin = 180,
              xmax = 260,
              ymin = 10,
              ymax = 62)

df210=data.frame(xmin = 210,
              xmax = 260,
              ymin = 10,
              ymax = 62)

e <- as(extent(to360(-180), to360(-100), 10, 62), 'SpatialPolygons')
e210 <- as(extent(210, to360(-100), 10, 62), 'SpatialPolygons')

data=maps::map("world2",fill=T)
IDs <- sapply(strsplit(data$names, ":"), function(x) x[1])
wrld_simpl <- map2SpatialPolygons(data, IDs=IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))
wrld=SpatialPolygons(wrld_simpl@polygons,proj4string=wrld_simpl@proj4string) %>% 
  gBuffer(., byid=TRUE, width=0)
wrld2=st_as_sf(wrld) %>% st_crop( xmin=180, xmax=260, ymax=62, ymin=10) %>% as_Spatial()
wrld210=st_as_sf(wrld) %>% st_crop( xmin=210, xmax=260, ymax=62, ymin=10) %>% as_Spatial()

eez=st_read("/Users/heatherwelch/Dropbox/woody/World_EEZ_v11_20191118_HR_0_360\ 2/eez_v11_0_360.shp") 

eez2=eez%>% 
  st_simplify(preserveTopology=TRUE, dTolerance = .2)

a=as_Spatial(eez)
b=crop(a,e)
b210=crop(a,e210)

## 5 degree RFMO
template=raster()
res(template)=5
nums=seq(1,ncell(template))
template=setValues(template,nums)

template360=raster::shift(raster::rotate(raster::shift(template,dx=180)),dx=180)
coords360=rasterToPoints(template360) %>% 
  as.data.frame() %>% rename(cell_lat=y,cell_lon=x,cells=layer) %>% 
  mutate(cells=as.factor(cells))
template360_SA=crop(template360,e)
template360_SA_shp=rasterToPolygons(template360_SA) 

## 1 degree gaps
template_overlap=raster("/Users/HeatherWelch/Dropbox/OLE/spatial_data/template.grd")
template_overlap_1deg=template_overlap
res(template_overlap_1deg)=1
nums=seq(1,ncell(template_overlap_1deg))
template_overlap_1deg=setValues(template_overlap_1deg,nums)
template_overlap_1deg_coords=rasterToPoints(template_overlap_1deg) %>% 
  as.data.frame() %>% rename(cell_lat=y,cell_lon=x,cells=Sea.level.anomaly) %>% 
  mutate(cells=as.factor(cells))

# overlap gaps data ####
datdir="UserDefined"
master_join=read.csv(glue("{datdir}/cleaned_metrics_03-23-23.csv"))
glimpse(master_join)
masterAIS=read.csv(glue("{datdir}/AIS_metrics_03-23-23.csv"))

# RFMO bycatch data ####
### 2. IATTC batch data for ll and ps, no USA ####
# source: https://www.iattc.org/en-US/Data/Public-domain
# a.	Contain flag
# c.	Lats and lons are cell mid points
# mako, blue, no salmon or white sharks

iattc_bycatch_ll=read.csv("/Users/heatherwelch/Dropbox/OLE/data_OLE/iattc_effort_data/PublicLLShark/PublicLLSharkNum.csv")%>% 
  filter(Flag!="USA") %>% 
  filter(Year>=2017&Year<=2021) %>% 
  mutate(lat=LatC5,lon360=to360(LonC5)) %>% 
  dplyr::select(c(lat,lon360,BSHn,MAKn))%>% 
  gather(species,n,-c(lat,lon360)) %>% 
  group_by(lat,lon360) %>% 
  summarise(number=sum(n,na.rm=T)) %>% 
  mutate(source="iattc_bycatch_ll")


iattc_bycatch_ps=read.csv("/Users/heatherwelch/Dropbox/OLE/data_OLE/iattc_effort_data/PublicPSShark/PublicPSSharkFlag.csv") %>% 
  filter(Flag!="USA") %>% 
  filter(Year>=2017&Year<=2021) %>% 
  filter(Flag!="USA") %>% 
  mutate(lat=LatC1,lon360=to360(LonC1)) %>% 
  dplyr::select(c(lat,lon360,BSHn)) %>% 
  gather(species,n,-c(lat,lon360)) %>% 
  group_by(lat,lon360) %>% 
  summarise(number=sum(n,na.rm=T)) %>% 
  mutate(source="iattc_bycatch_ps")

master_rfmo_sharks=do.call("rbind",list(iattc_bycatch_ll,iattc_bycatch_ps)) %>% 
  ungroup() %>% 
  mutate(cells=as.factor(raster::extract(template360,.[,2:1]))) %>% 
  group_by(cells) %>% 
  summarise(number=sum(number,na.rm=T)) %>% 
  ungroup() %>% 
  left_join(.,coords360)


# overlaps-shark unseen data disabling, no USA ####
sharks_overlap=master_join %>% 
  filter(!is.na(species)) %>%
  filter(flag!="USA") %>% 
  mutate(year=year(as.Date(date))) %>% 
  filter(year<=2021) %>% 
  filter(species=="blueShark_TOPP"|
           species=="makoShark_TOPP") %>% 
  filter(vessel_class=="tuna_purse_seines"|
           vessel_class=="drifting_longlines") %>% 
  
  mutate(overlap_hours_in_gaps_over_12=hours_in_gaps_over_12*hab.suit,
         overlap_hours_in_gaps_under_12=hours_in_gaps_under_12*hab.suit,
         overlap_total_gap_2w=total_gap_2w*hab.suit,
         overlap_total_gap=total_gap*hab.suit,
         overlap_total_gap_disabling=total_gap_disabling*hab.suit,
         overlap_total_gap_disabling_2w=total_gap_disabling_2w*hab.suit) %>%
  
  mutate(cells=as.factor(raster::extract(template360,.[,c("cell_lon","cell_lat")]))) %>% 
  dplyr::select(-c(cell_lon,cell_lat,date,vessel_class,flag,species,hab.suit,flag,year,X)) %>% 
  group_by(cells) %>%
  summarise(overlap_hours_in_gaps_over_12=sum(overlap_hours_in_gaps_over_12,na.rm=T),
            overlap_hours_in_gaps_under_12=sum(overlap_hours_in_gaps_under_12,na.rm=T),
            overlap_total_gap_disabling_2w=sum(overlap_total_gap_disabling_2w,na.rm=T),
            overlap_total_gap_disabling=sum(overlap_total_gap_disabling,na.rm=T),
            overlap_total_gap_2w=sum(overlap_total_gap_2w,na.rm=T),
            overlap_total_gap=sum(overlap_total_gap,na.rm=T)
            ) %>% 
    
  mutate(overlap_observed=overlap_hours_in_gaps_over_12+overlap_hours_in_gaps_under_12) %>%  ## observed ais
  mutate(overlap_total_2w_activity=overlap_total_gap_2w+overlap_hours_in_gaps_under_12) %>% # 2w observed + gaps
  mutate(overlap_total_activity=overlap_total_gap+overlap_hours_in_gaps_under_12) %>% #  observed + gaps
  left_join(.,coords360) 


# bycatch figure ####
# 1, clip RFMO to my data
master_bycatch=sharks_overlap %>% 
  left_join(.,master_rfmo_sharks) %>% 
  filter(cell_lon>=210&cell_lon<=260) %>%  ## iattc only, can't deal w US vessels wcpfc
  filter(cell_lat>=10&cell_lat<=62)  %>% .[complete.cases(.),]

master_bycatch_noout=master_bycatch %>% 
  filter(cells!=1021 & cells!=950 & cells!=1022 & cells!=1023 & cells!=1094 & cells!=1092  & cells!=1093 ) ## 7 anomalous areas
mod_noout=lm(overlap_total_gap_disabling_2w~number,data=master_bycatch_noout)
coefficients=coef(mod_noout)
mod=lm(overlap_total_gap_disabling_2w~number,data=master_bycatch)
master_bycatch_preds=master_bycatch %>% 
  mutate(pred_high=predict(mod,.),
         pred_low=predict(mod_noout,.))

master_bycatch_out=master_bycatch %>% 
  filter(cells==1021 | cells==950 | cells==1022 | cells==1023 | cells==1094 | cells==1092  | cells==1093 ) 

max_y=max(master_bycatch$overlap_total_gap_disabling_2w)
scatter_sharks=ggplot()+
  geom_point(data=master_bycatch_out,aes(x=number,y=overlap_total_gap_disabling_2w),color="red",size=3.5,shape=0,stroke=2)+
  geom_point(data=master_bycatch_noout,aes(x=number,y=overlap_total_gap_disabling_2w),size=3)+
  geom_smooth(data=master_bycatch,aes(x=number,y=overlap_total_gap_disabling_2w),se=F,method=lm,color="blue")+
  geom_smooth(data=master_bycatch_noout,aes(x=number,y=overlap_total_gap_disabling_2w),se=F,method=lm,color="black")+
   geom_segment(aes(x = 0, xend = 8390.0, y = coefficients[1] + 3.556e-05, yend =  coefficients[1] + 3.556e-05*8390.0),linetype="dashed",color="black") +
  scale_y_continuous("Unseen overlap with sharks",breaks=c(0,max_y),
                     labels=c("Lower","Higher"))+
  theme_classic()+
  xlab("Reported number of sharks caught (IATTC)")

scatter_sharks

png(glue("{outdir}/scatter_sharks.png"),width=14,height=10,units='cm',res=400,type = "cairo")
par(ps=10)
par(mar=c(4,4,1,1))
par(cex=1)
print({scatter_sharks})
# gg_hm
dev.off()


### maps ####
lat_lon_full=data.frame(label=c("40°N    ","20°N   ","220°","240°"),
                        x=c(208,208,220,240),
                        y=c(40,20,9,9)
) ## for drawing parallels and meridians

template360_SA_shp_clip=template360_SA_shp %>% 
  st_as_sf() %>% 
  filter(layer=='1023'|
           layer=="950"|
           layer=="1094"|
           layer=="1021"|
           layer=="1022"|
           layer=="1093"|
           layer=="1092") %>% 
  as_Spatial()

rfmo_shakrs=ggplot() +
  geom_rect(data=df210,color = "black",fill="white",
            aes(xmin = xmin,
                xmax = xmax,
                ymin = ymin,
                ymax = ymax),alpha=1)+
  geom_tile(data=master_bycatch, aes(x=cell_lon, y=cell_lat,fill=number)) +
  scale_fill_gradientn("Number",colours = pals::parula(100),na.value=NA,
                       breaks=c(1,4000,8000),
                       labels=c("0","4000","8000"))+
  geom_segment(aes(x=210,xend=260,y=60,yend=60),color="grey",linetype="dashed")+
  geom_segment(aes(x=210,xend=260,y=40,yend=40),color="grey",linetype="dashed")+
  geom_segment(aes(x=210,xend=260,y=20,yend=20),color="grey",linetype="dashed")+
  geom_segment(aes(x=220,xend=220,y=10,yend=62),color="grey",linetype="dashed")+
  geom_segment(aes(x=240,xend=240,y=10,yend=62),color="grey",linetype="dashed")+
  geom_text(data=lat_lon_full,aes(x=x,y=y,label=label),color="darkgrey",size=3)+
  geom_path(data=fortify(b210),aes(x=long, y = lat, group=group),color="black",fill=NA,size=.5)+
  geom_polygon(data=fortify(template360_SA_shp_clip),aes(x=long, y = lat, group=group),color="red",fill=NA,size=1)+
  coord_map("conic", lat0 = 30,xlim = c(210, 260), ylim = c(10,62))+
  geom_polygon(data=fortify(wrld210),aes(x=long, y = lat, group=group),color="#688f74",fill="grey",alpha=1,size=.2)+
  theme_void() +
  theme(legend.position = "bottom")

# rfmo_shakrs

png(glue("{outdir}/rfmo_sharks1.png"),width=12,height=10.5,units='cm',res=400,type = "cairo")
par(ps=10)
par(mar=c(4,4,1,1))
par(cex=1)
print({rfmo_shakrs})
# gg_hm
dev.off()

overlap_shakrs=ggplot() +
  geom_rect(data=df210,color = "black",fill="white",
            aes(xmin = xmin,
                xmax = xmax,
                ymin = ymin,
                ymax = ymax),alpha=1)+
  geom_tile(data=master_bycatch, aes(x=cell_lon, y=cell_lat,fill=overlap_total_gap_disabling_2w)) +
  scale_fill_gradientn("Overlap",colours = pals::parula(100),na.value=NA,
                       breaks=c(0,.45),
                       labels=c("Lower","Higher"))+
  geom_segment(aes(x=210,xend=260,y=60,yend=60),color="grey",linetype="dashed")+
  geom_segment(aes(x=210,xend=260,y=40,yend=40),color="grey",linetype="dashed")+
  geom_segment(aes(x=210,xend=260,y=20,yend=20),color="grey",linetype="dashed")+
  geom_segment(aes(x=220,xend=220,y=10,yend=62),color="grey",linetype="dashed")+
  geom_segment(aes(x=240,xend=240,y=10,yend=62),color="grey",linetype="dashed")+
  geom_text(data=lat_lon_full,aes(x=x,y=y,label=label),color="darkgrey",size=3)+
  geom_path(data=fortify(b210),aes(x=long, y = lat, group=group),color="black",fill=NA,size=.5)+
  geom_polygon(data=fortify(template360_SA_shp_clip),aes(x=long, y = lat, group=group),color="red",fill=NA,size=1)+
  coord_map("conic", lat0 = 30,xlim = c(210, 260), ylim = c(10,62))+
  geom_polygon(data=fortify(wrld210),aes(x=long, y = lat, group=group),color="#688f74",fill="grey",alpha=1,size=.2)+
  theme_void() +
  theme(legend.position = "bottom")

png(glue("{outdir}/overlap_sharks1.png"),width=12,height=10.5,units='cm',res=400,type = "cairo")
par(ps=10)
par(mar=c(4,4,1,1))
par(cex=1)
print({overlap_shakrs})
# gg_hm
dev.off()

## correlations ####
cordat_full=master_bycatch %>% 
  dplyr::select(number,number_overlap) %>% 
  cor()

cordat_noout=master_bycatch_noout %>% 
  dplyr::select(number,number_overlap) %>% 
  cor()


## Unseen catch based on lines of best fit ####
# filter(cells!=1021 & cells!=950 & cells!=1022 & cells!=1023 & cells!=1094 & cells!=1092  & cells!=1093 ) ## 7 anomalous areas
anomalies=master_bycatch %>% 
  filter(cells=='1021'|
           cells=="950"|
           cells=="1022"|
           cells=="1093"|
           cells=="1094"|
           cells=="1092"|
           cells=="1023")

lm_full=lm(overlap_total_gap_disabling_2w~number,data=master_bycatch)
lm_no_anom=lm(overlap_total_gap_disabling_2w~number,data=master_bycatch_noout)

intercept_full=coefficients(lm_full)[1]
intercept_no_anom=coefficients(lm_no_anom)[1]
slope_full=coefficients(lm_full)[2]
slope_no_anom=coefficients(lm_no_anom)[2]

empty=list()
for(i in 1:nrow(anomalies)){
  dat=anomalies[i,]
  number=dat %>% pull(number)
  overlap_total_gap_disabling_2w=dat %>% pull(overlap_total_gap_disabling_2w)
  
  ## full line
  x_full=(overlap_total_gap_disabling_2w-intercept_full)/slope_full
  diff_full=x_full-number
  
  ## partial line
  x_no_anom=(overlap_total_gap_disabling_2w-intercept_no_anom)/slope_no_anom
  diff_no_anom=x_no_anom-number
  
  df=data.frame(
    cell=dat %>% pull(cells),
    number1=number,
    number_overlap=overlap_total_gap_disabling_2w,
    absolute_full=x_full,
    absolute_no_anom=x_no_anom,
    diff_full=diff_full,
    diff_no_anom=diff_no_anom
  )

  empty[[length(empty)+1]]=df
}

stats_master=do.call("rbind",empty) 
stats_master2=stats_master%>% 
  summarise(diff_full=sum(diff_full),
            diff_no_anom=sum(diff_no_anom),
            number=sum(number1)) 

# which flags are in those four cells? ####
# filter(cells!=1021 & cells!=950 & cells!=1022 & cells!=1023 & cells!=1094 & cells!=1092  & cells!=1093 )
my_data=master_join %>% 
  filter(!is.na(species)) %>%
  filter(flag!="USA") %>% 
  mutate(year=year(as.Date(date))) %>% 
  filter(year<=2021) %>% 
  filter(species=="blueShark_TOPP"|
           species=="makoShark_TOPP") %>% 
  filter(vessel_class=="tuna_purse_seines"|
           vessel_class=="drifting_longlines") %>% 
  
  mutate(overlap_hours_in_gaps_over_12=hours_in_gaps_over_12*hab.suit,
         overlap_hours_in_gaps_under_12=hours_in_gaps_under_12*hab.suit,
         overlap_total_gap_2w=total_gap_2w*hab.suit,
         overlap_total_gap=total_gap*hab.suit,
         overlap_total_gap_disabling=total_gap_disabling*hab.suit,
         overlap_total_gap_disabling_2w=total_gap_disabling_2w*hab.suit) %>%
  
  mutate(cells=as.factor(raster::extract(template360,.[,c("cell_lon","cell_lat")]))) %>% 
  filter(cells=='1021'|
           cells=="950"|
           cells=="1022"|
           cells=="1093"|
           cells=="1094"|
           cells=="1092"|
           cells=="1023") %>%
  dplyr::select(-c(cell_lon,cell_lat,date,vessel_class,species,hab.suit,cells,year,X)) %>% 
  group_by(flag) %>%
  summarise(overlap_hours_in_gaps_over_12=sum(overlap_hours_in_gaps_over_12,na.rm=T),
            overlap_hours_in_gaps_under_12=sum(overlap_hours_in_gaps_under_12,na.rm=T),
            overlap_total_gap_disabling_2w=sum(overlap_total_gap_disabling_2w,na.rm=T),
            overlap_total_gap_disabling=sum(overlap_total_gap_disabling,na.rm=T),
            overlap_total_gap_2w=sum(overlap_total_gap_2w,na.rm=T),
            overlap_total_gap=sum(overlap_total_gap,na.rm=T)
  ) %>% 
  
  mutate(overlap_observed=overlap_hours_in_gaps_over_12+overlap_hours_in_gaps_under_12) %>%  ## observed ais
  mutate(overlap_total_2w_activity=overlap_total_gap_2w+overlap_hours_in_gaps_under_12) %>% # 2w observed + gaps
  mutate(overlap_total_activity=overlap_total_gap+overlap_hours_in_gaps_under_12)  #  observed + gaps
  

iattc_ll=read.csv("/Users/heatherwelch/Dropbox/OLE/data_OLE/iattc_effort_data/PublicLLShark/PublicLLSharkNum.csv")%>% 
   filter(Flag!="USA") %>% 
  filter(Year>=2017&Year<=2021) %>% 
  mutate(lat=LatC5,lon360=to360(LonC5)) %>% 
  dplyr::select(c(lat,lon360,BSHn,MAKn,Flag))%>% 
  gather(species,n,-c(lat,lon360,Flag)) %>% 
  group_by(lat,lon360,Flag) %>% 
  summarise(number=sum(n,na.rm=T)) %>% 
  mutate(source="iattc_bycatch_ll")


iattc_ps=read.csv("/Users/heatherwelch/Dropbox/OLE/data_OLE/iattc_effort_data/PublicPSShark/PublicPSSharkFlag.csv") %>% 
  filter(Year>=2017&Year<=2021) %>% 
   filter(Flag!="USA") %>% 
  mutate(lat=LatC1,lon360=to360(LonC1)) %>% 
  dplyr::select(c(lat,lon360,BSHn,Flag)) %>% 
  gather(species,n,-c(lat,lon360,Flag)) %>% 
  group_by(lat,lon360,Flag) %>% 
  summarise(number=sum(n,na.rm=T)) %>% 
  mutate(source="iattc_bycatch_ps")

master_rfmo_sharks_flag=do.call("rbind",list(iattc_ll,iattc_ps)) %>% 
  ungroup() %>% 
  mutate(cells=as.factor(raster::extract(template360,.[,2:1]))) %>% 
  filter(cells=='1021'|
           cells=="950"|
           cells=="1022"|
           cells=="1093"|
           cells=="1094"|
           cells=="1092"|
           cells=="1023") %>%
  group_by(Flag) %>% 
  summarise(number=sum(number,na.rm=T)) 

test=full_join(my_data,master_rfmo_sharks_flag,by=c("flag"="Flag"))%>%
  mutate(overlap_total_gap_disabling_2w_total=sum(overlap_total_gap_disabling_2w),
         number_total=sum(number,na.rm=T)) %>%
  mutate(overlap_total_gap_disabling_2w_perc=overlap_total_gap_disabling_2w/overlap_total_gap_disabling_2w_total*100,
         number_perc=number/number_total*100) %>% 
  dplyr::select(c(overlap_total_gap_disabling_2w_perc,number_perc)) %>% arrange(desc(number_perc))

# which flags are in those four cells BY CELL? ####
# filter(cells!=1021 & cells!=950 & cells!=1022 & cells!=1023 & cells!=1094 & cells!=1092  & cells!=1093 )
my_data=master_join %>% 
  filter(!is.na(species)) %>%
  filter(flag!="USA") %>% 
  mutate(year=year(as.Date(date))) %>% 
  filter(year<=2021) %>% 
  filter(species=="blueShark_TOPP"|
           species=="makoShark_TOPP") %>% 
  filter(vessel_class=="tuna_purse_seines"|
           vessel_class=="drifting_longlines") %>% 
  
  mutate(overlap_hours_in_gaps_over_12=hours_in_gaps_over_12*hab.suit,
         overlap_hours_in_gaps_under_12=hours_in_gaps_under_12*hab.suit,
         overlap_total_gap_2w=total_gap_2w*hab.suit,
         overlap_total_gap=total_gap*hab.suit,
         overlap_total_gap_disabling=total_gap_disabling*hab.suit,
         overlap_total_gap_disabling_2w=total_gap_disabling_2w*hab.suit) %>%
  
  mutate(cells=as.factor(raster::extract(template360,.[,c("cell_lon","cell_lat")]))) %>% 
  filter(cells=='1021'|
           cells=="950"|
           cells=="1022"|
           cells=="1093"|
           cells=="1094"|
           cells=="1092"|
           cells=="1023") %>%
  dplyr::select(-c(cell_lon,cell_lat,date,vessel_class,species,hab.suit,year,X)) %>% 
  group_by(flag,cells) %>%
  summarise(overlap_hours_in_gaps_over_12=sum(overlap_hours_in_gaps_over_12,na.rm=T),
            overlap_hours_in_gaps_under_12=sum(overlap_hours_in_gaps_under_12,na.rm=T),
            overlap_total_gap_disabling_2w=sum(overlap_total_gap_disabling_2w,na.rm=T),
            overlap_total_gap_disabling=sum(overlap_total_gap_disabling,na.rm=T),
            overlap_total_gap_2w=sum(overlap_total_gap_2w,na.rm=T),
            overlap_total_gap=sum(overlap_total_gap,na.rm=T)
  ) %>% 
  
  mutate(overlap_observed=overlap_hours_in_gaps_over_12+overlap_hours_in_gaps_under_12) %>%  ## observed ais
  mutate(overlap_total_2w_activity=overlap_total_gap_2w+overlap_hours_in_gaps_under_12) %>% # 2w observed + gaps
  mutate(overlap_total_activity=overlap_total_gap+overlap_hours_in_gaps_under_12)  #  observed + gaps


iattc_ll=read.csv("/Users/heatherwelch/Dropbox/OLE/data_OLE/iattc_effort_data/PublicLLShark/PublicLLSharkNum.csv")%>% 
  filter(Flag!="USA") %>% 
  filter(Year>=2017&Year<=2021) %>% 
  mutate(lat=LatC5,lon360=to360(LonC5)) %>% 
  dplyr::select(c(lat,lon360,BSHn,MAKn,Flag))%>% 
  gather(species,n,-c(lat,lon360,Flag)) %>% 
  group_by(lat,lon360,Flag) %>% 
  summarise(number=sum(n,na.rm=T)) %>% 
  mutate(source="iattc_bycatch_ll")


iattc_ps=read.csv("/Users/heatherwelch/Dropbox/OLE/data_OLE/iattc_effort_data/PublicPSShark/PublicPSSharkFlag.csv") %>% 
  filter(Year>=2017&Year<=2021) %>% 
  filter(Flag!="USA") %>% 
  mutate(lat=LatC1,lon360=to360(LonC1)) %>% 
  dplyr::select(c(lat,lon360,BSHn,Flag)) %>% 
  gather(species,n,-c(lat,lon360,Flag)) %>% 
  group_by(lat,lon360,Flag) %>% 
  summarise(number=sum(n,na.rm=T)) %>% 
  mutate(source="iattc_bycatch_ps")

master_rfmo_sharks_flag=do.call("rbind",list(iattc_ll,iattc_ps)) %>% 
  ungroup() %>% 
  mutate(cells=as.factor(raster::extract(template360,.[,2:1]))) %>% 
  filter(cells=='1021'|
           cells=="950"|
           cells=="1022"|
           cells=="1093"|
           cells=="1094"|
           cells=="1092"|
           cells=="1023") %>%
  rename(flag=Flag) %>% 
  group_by(flag,cells) %>% 
  summarise(number=sum(number,na.rm=T)) 

b=master_rfmo_sharks_flag %>% 
  filter(flag=="CHN"|flag=="MEX")
  
b_mine=my_data %>% 
  dplyr::select(c(flag,cells,overlap_total_gap_disabling_2w)) %>% 
  filter(flag=="CHN"|flag=="MEX")

b_master=full_join(b,b_mine) %>% 
  mutate(number=replace_na(number,0))

together=left_join(b_master,stats_master,by=c("cells"="cell")) %>% 
  mutate(percent_overlap=overlap_total_gap_disabling_2w/number_overlap) %>%  ## how much of disabling in each cell is each country responsible for
  mutate(iso_diff_full=diff_full*percent_overlap,
         iso_diff_no_anom=diff_no_anom*percent_overlap
         )

country_stats= together %>% 
  group_by(flag) %>% 
  summarise(number=sum(number),
            iso_diff_full=sum(iso_diff_full),
            iso_diff_no_anom=sum(iso_diff_no_anom)) %>% 
  ungroup() %>% 
  mutate(lower_total=sum(iso_diff_full),
         upper_total=sum(iso_diff_no_anom)) %>% 
  mutate(iso_lower=iso_diff_full/lower_total,
         iso_upper=iso_diff_no_anom/upper_total)
  
# Mexico target catch ####
iattc_catch_ps=read.csv("/Users/heatherwelch/Dropbox/OLE/data_OLE/iattc_effort_data/PublicPSTuna_1/PublicPSTunaFlag.csv") %>% 
  filter(Year>=2017&Year<=2022) %>% 
  filter(Flag=="MEX") %>% 
  mutate(lat=LatC1,lon360=to360(LonC1)) %>% 
  dplyr::select(-c(Year, Month, Flag, LatC1, LonC1, NumSets)) %>% 
  gather(species,number,-c(lat,lon360)) %>% 
  mutate(cells=as.factor(raster::extract(template360,.[,2:1]))) %>% 
  group_by(cells) %>% 
  summarise(number=sum(number,na.rm=T))%>% 
  ungroup() %>% 
  left_join(.,coords360)%>% 
  filter(cell_lat<62 & cell_lat>10)%>% 
  filter(cell_lon<260 & cell_lon>180)

## no mexican longline catch
iattc_catch_ll=read.csv("/Users/heatherwelch/Dropbox/OLE/data_OLE/iattc_effort_data/PublicLLTunaBillfish/PublicLLTunaBillfishNum.csv") %>% 
  filter(Year>=2017&Year<=2022) %>% 
  filter(Flag=="MEX") %>% 
  mutate(lat=LatC5,lon360=to360(LonC5)) %>% 
  dplyr::select(-c(Year,Month,Flag,LatC5,LonC5,Hooks,
                   ALBmt,BETmt,PBFmt,SKJmt,TUNmt,YFTmt,
                   BLMmt,BILmt,BUMmt,MLSmt,SFAmt,SSPmt,
                   SWOmt
                   )) %>% 
  gather(species,n,-c(lat,lon360)) %>% 
  group_by(lat,lon360) %>% 
  summarise(n=sum(n,na.rm=T))%>% 
  filter(lat<62 & lat>10)%>% 
  filter(lon360<260 & lon360>180)

  
rfmo_mex_tuna=ggplot() +
  geom_rect(data=df210,color = "black",fill="white",
            aes(xmin = xmin,
                xmax = xmax,
                ymin = ymin,
                ymax = ymax),alpha=1)+
  geom_tile(data=iattc_catch_ps, aes(x=cell_lon, y=cell_lat,fill=number)) +
  scale_fill_gradientn("Number",colours = pals::parula(100),na.value=NA,
                       breaks=c(10000,20000,30000),
                       labels=c("10k","20k","30k"))+
  geom_segment(aes(x=210,xend=260,y=60,yend=60),color="grey",linetype="dashed")+
   geom_segment(aes(x=210,xend=260,y=40,yend=40),color="grey",linetype="dashed")+
   geom_segment(aes(x=210,xend=260,y=20,yend=20),color="grey",linetype="dashed")+
   geom_segment(aes(x=220,xend=220,y=10,yend=62),color="grey",linetype="dashed")+
  geom_segment(aes(x=240,xend=240,y=10,yend=62),color="grey",linetype="dashed")+
    geom_text(data=lat_lon,aes(x=x,y=y,label=label),color="darkgrey",size=2)+
  geom_polygon(data=fortify(b210),aes(x=long, y = lat, group=group),color="black",fill=NA)+
  geom_polygon(data=fortify(template360_SA_shp_clip),aes(x=long, y = lat, group=group),color="red",fill=NA)+
  coord_map("conic", lat0 = 30,xlim = c(210, 260), ylim = c(10,62))+
   geom_polygon(data=fortify(wrld210),aes(x=long, y = lat, group=group),color="black",fill="grey",alpha=1)+
  theme_void() +
  theme(legend.position = "bottom")

png(glue("{outdir}/rfmo_mex_tuna.png"),width=12,height=10.5,units='cm',res=400,type = "cairo")
par(ps=10)
par(mar=c(4,4,1,1))
par(cex=1)
print({rfmo_mex_tuna})
# gg_hm
dev.off()

# China target catch ####
#no china ps catch
iattc_catch_ps=read.csv("/Users/heatherwelch/Dropbox/OLE/data_OLE/iattc_effort_data/PublicPSTuna_1/PublicPSTunaFlag.csv") %>% 
  filter(Year>=2017&Year<=2022) %>% 
  filter(Flag=="CHN") %>% 
  mutate(lat=LatC1,lon360=to360(LonC1)) %>% 
  dplyr::select(-c(Year, Month, Flag, LatC1, LonC1, NumSets)) %>% 
  gather(species,number,-c(lat,lon360)) %>% 
  mutate(cells=as.factor(raster::extract(template360,.[,2:1]))) %>% 
  group_by(cells) %>% 
  summarise(number=sum(number,na.rm=T))%>% 
  ungroup() %>% 
  left_join(.,coords360)%>% 
  filter(cell_lat<62 & cell_lat>10)%>% 
  filter(cell_lon<260 & cell_lon>180)

iattc_catch_ll=read.csv("/Users/heatherwelch/Dropbox/OLE/data_OLE/iattc_effort_data/PublicLLTunaBillfish/PublicLLTunaBillfishNum.csv") %>% 
  filter(Year>=2017&Year<=2022) %>% 
  filter(Flag=="CHN") %>% 
  mutate(lat=LatC5,lon360=to360(LonC5)) %>% 
  dplyr::select(-c(Year,Month,Flag,LatC5,LonC5,Hooks,
                   ALBmt,BETmt,PBFmt,SKJmt,TUNmt,YFTmt,
                   BLMmt,BILmt,BUMmt,MLSmt,SFAmt,SSPmt,
                   SWOmt
  )) %>% 
  gather(species,number,-c(lat,lon360)) %>% 
  mutate(cells=as.factor(raster::extract(template360,.[,2:1]))) %>% 
  group_by(cells) %>% 
  summarise(number=sum(number,na.rm=T))%>% 
  ungroup() %>% 
  left_join(.,coords360)%>% 
  filter(cell_lat<62 & cell_lat>10)%>% 
  filter(cell_lon<260 & cell_lon>180)

rfmo_chn_tuna=ggplot() +
  geom_rect(data=df210,color = "black",fill="white",
            aes(xmin = xmin,
                xmax = xmax,
                ymin = ymin,
                ymax = ymax),alpha=1)+
  geom_tile(data=iattc_catch_ll, aes(x=cell_lon, y=cell_lat,fill=number)) +
  scale_fill_gradientn("Number",colours = pals::parula(100),na.value=NA,
                       breaks=c(10000,20000,30000),
                       labels=c("10k","20k","30k"))+
  geom_segment(aes(x=210,xend=260,y=60,yend=60),color="grey",linetype="dashed")+
  geom_segment(aes(x=210,xend=260,y=40,yend=40),color="grey",linetype="dashed")+
  geom_segment(aes(x=210,xend=260,y=20,yend=20),color="grey",linetype="dashed")+
  geom_segment(aes(x=220,xend=220,y=10,yend=62),color="grey",linetype="dashed")+
  geom_segment(aes(x=240,xend=240,y=10,yend=62),color="grey",linetype="dashed")+
  geom_text(data=lat_lon,aes(x=x,y=y,label=label),color="darkgrey",size=2)+
  geom_polygon(data=fortify(b210),aes(x=long, y = lat, group=group),color="black",fill=NA)+
  geom_polygon(data=fortify(template360_SA_shp_clip),aes(x=long, y = lat, group=group),color="red",fill=NA)+
  coord_map("conic", lat0 = 30,xlim = c(210, 260), ylim = c(10,62))+
  geom_polygon(data=fortify(wrld210),aes(x=long, y = lat, group=group),color="black",fill="grey",alpha=1)+
  theme_void() +
  theme(legend.position = "bottom")

png(glue("{outdir}/rfmo_chn_tuna.png"),width=12,height=10.5,units='cm',res=400,type = "cairo")
par(ps=10)
par(mar=c(4,4,1,1))
par(cex=1)
print({rfmo_chn_tuna})
# gg_hm
dev.off()

# Mexico sharks during study period ####
iattc_ll=read.csv("/Users/heatherwelch/Dropbox/OLE/data_OLE/iattc_effort_data/PublicLLShark/PublicLLSharkNum.csv")%>% 
  filter(Flag=="MEX") %>% 
  filter(Year>=2017&Year<=2021) %>% 
  mutate(lat=LatC5,lon360=to360(LonC5)) %>% 
  dplyr::select(c(lat,lon360,BSHn,MAKn,Flag))%>% 
  gather(species,n,-c(lat,lon360,Flag)) %>% 
  group_by(lat,lon360,Flag) %>% 
  summarise(number=sum(n,na.rm=T)) %>% 
  mutate(source="iattc_bycatch_ll")


iattc_ps=read.csv("/Users/heatherwelch/Dropbox/OLE/data_OLE/iattc_effort_data/PublicPSShark/PublicPSSharkFlag.csv") %>% 
  filter(Year>=2017&Year<=2021) %>% 
  filter(Flag=="MEX") %>% 
  mutate(lat=LatC1,lon360=to360(LonC1)) %>% 
  dplyr::select(c(lat,lon360,BSHn,Flag)) %>% 
  gather(species,n,-c(lat,lon360,Flag)) %>% 
  group_by(lat,lon360,Flag) %>% 
  summarise(number=sum(n,na.rm=T)) %>% 
  mutate(source="iattc_bycatch_ps")

master_rfmo_sharks_flag=do.call("rbind",list(iattc_ll,iattc_ps)) %>% 
  ungroup() %>% 
  mutate(cells=as.factor(raster::extract(template360,.[,2:1]))) %>% 
  group_by(cells) %>% 
  summarise(number=sum(number,na.rm=T)) %>% 
  ungroup() %>% 
  left_join(.,coords360) %>% 
  filter(cell_lat<62 & cell_lat>10)%>% 
  filter(cell_lon<260 & cell_lon>210)
  
mex_sharks_studyperiod=ggplot() +
  geom_rect(data=df210,color = "black",fill="white",
            aes(xmin = xmin,
                xmax = xmax,
                ymin = ymin,
                ymax = ymax),alpha=1)+
  geom_tile(data=master_rfmo_sharks_flag, aes(x=cell_lon, y=cell_lat,fill=number)) +
  scale_fill_gradientn("Number",colours = pals::parula(100),na.value=NA)+
  geom_segment(aes(x=210,xend=260,y=60,yend=60),color="grey",linetype="dashed")+
  geom_segment(aes(x=210,xend=260,y=40,yend=40),color="grey",linetype="dashed")+
  geom_segment(aes(x=210,xend=260,y=20,yend=20),color="grey",linetype="dashed")+
  geom_segment(aes(x=220,xend=220,y=10,yend=62),color="grey",linetype="dashed")+
  geom_segment(aes(x=240,xend=240,y=10,yend=62),color="grey",linetype="dashed")+
  geom_text(data=lat_lon,aes(x=x,y=y,label=label),color="darkgrey",size=2)+
  geom_polygon(data=fortify(b210),aes(x=long, y = lat, group=group),color="black",fill=NA)+
  geom_polygon(data=fortify(template360_SA_shp_clip),aes(x=long, y = lat, group=group),color="red",fill=NA)+
  coord_map("conic", lat0 = 30,xlim = c(210, 260), ylim = c(10,62))+
  geom_polygon(data=fortify(wrld210),aes(x=long, y = lat, group=group),color="black",fill="grey",alpha=1)+
  theme_void() +
  theme(legend.position = "bottom")

png(glue("{outdir}/mex_sharks_studyperiod.png"),width=12,height=10.5,units='cm',res=400,type = "cairo")
par(ps=10)
par(mar=c(4,4,1,1))
par(cex=1)
print({mex_sharks_studyperiod})
# gg_hm
dev.off()

# China sharks during study period ####
iattc_ll=read.csv("/Users/heatherwelch/Dropbox/OLE/data_OLE/iattc_effort_data/PublicLLShark/PublicLLSharkNum.csv")%>% 
  filter(Flag=="CHN") %>% 
  filter(Year>=2017&Year<=2021) %>% 
  mutate(lat=LatC5,lon360=to360(LonC5)) %>% 
  dplyr::select(c(lat,lon360,BSHn,MAKn,Flag))%>% 
  gather(species,n,-c(lat,lon360,Flag)) %>% 
  group_by(lat,lon360,Flag) %>% 
  summarise(number=sum(n,na.rm=T)) %>% 
  mutate(source="iattc_bycatch_ll")


iattc_ps=read.csv("/Users/heatherwelch/Dropbox/OLE/data_OLE/iattc_effort_data/PublicPSShark/PublicPSSharkFlag.csv") %>% 
  filter(Year>=2017&Year<=2021) %>% 
  filter(Flag=="CHN") %>% 
  mutate(lat=LatC1,lon360=to360(LonC1)) %>% 
  dplyr::select(c(lat,lon360,BSHn,Flag)) %>% 
  gather(species,n,-c(lat,lon360,Flag)) %>% 
  group_by(lat,lon360,Flag) %>% 
  summarise(number=sum(n,na.rm=T)) %>% 
  mutate(source="iattc_bycatch_ps")

master_rfmo_sharks_flag=do.call("rbind",list(iattc_ll,iattc_ps)) %>% 
  ungroup() %>% 
  mutate(cells=as.factor(raster::extract(template360,.[,2:1]))) %>% 
  group_by(cells) %>% 
  summarise(number=sum(number,na.rm=T)) %>% 
  ungroup() %>% 
  left_join(.,coords360) %>% 
  filter(cell_lat<62 & cell_lat>10)%>% 
  filter(cell_lon<260 & cell_lon>210)

chn_sharks_studyperiod=ggplot() +
  geom_rect(data=df210,color = "black",fill="white",
            aes(xmin = xmin,
                xmax = xmax,
                ymin = ymin,
                ymax = ymax),alpha=1)+
  geom_tile(data=master_rfmo_sharks_flag, aes(x=cell_lon, y=cell_lat,fill=number)) +
  scale_fill_gradientn("Number",colours = pals::parula(100),na.value=NA,
                         breaks=c(0,4000,8000),
                         labels=c("0","4k","8k"))+
  geom_segment(aes(x=210,xend=260,y=60,yend=60),color="grey",linetype="dashed")+
  geom_segment(aes(x=210,xend=260,y=40,yend=40),color="grey",linetype="dashed")+
  geom_segment(aes(x=210,xend=260,y=20,yend=20),color="grey",linetype="dashed")+
  geom_segment(aes(x=220,xend=220,y=10,yend=62),color="grey",linetype="dashed")+
  geom_segment(aes(x=240,xend=240,y=10,yend=62),color="grey",linetype="dashed")+
  geom_text(data=lat_lon,aes(x=x,y=y,label=label),color="darkgrey",size=2)+
  geom_polygon(data=fortify(b210),aes(x=long, y = lat, group=group),color="black",fill=NA)+
  geom_polygon(data=fortify(template360_SA_shp_clip),aes(x=long, y = lat, group=group),color="red",fill=NA)+
  coord_map("conic", lat0 = 30,xlim = c(210, 260), ylim = c(10,62))+
  geom_polygon(data=fortify(wrld210),aes(x=long, y = lat, group=group),color="black",fill="grey",alpha=1)+
  theme_void() +
  theme(legend.position = "bottom")

png(glue("{outdir}/chn_sharks_studyperiod.png"),width=12,height=10.5,units='cm',res=400,type = "cairo")
par(ps=10)
par(mar=c(4,4,1,1))
par(cex=1)
print({chn_sharks_studyperiod})
# gg_hm
dev.off()

# Mexico sharks historical ####
iattc_ll=read.csv("/Users/heatherwelch/Dropbox/OLE/data_OLE/iattc_effort_data/PublicLLShark/PublicLLSharkNum.csv")%>% 
  filter(Flag=="MEX") %>% 
  # filter(Year>=2017&Year<=2021) %>% 
  mutate(lat=LatC5,lon360=to360(LonC5)) %>% 
  dplyr::select(c(lat,lon360,BSHn,MAKn,Flag))%>% 
  gather(species,n,-c(lat,lon360,Flag)) %>% 
  group_by(lat,lon360,Flag) %>% 
  summarise(number=sum(n,na.rm=T)) %>% 
  mutate(source="iattc_bycatch_ll")


iattc_ps=read.csv("/Users/heatherwelch/Dropbox/OLE/data_OLE/iattc_effort_data/PublicPSShark/PublicPSSharkFlag.csv") %>% 
  # filter(Year>=2017&Year<=2021) %>% 
  filter(Flag=="MEX") %>% 
  mutate(lat=LatC1,lon360=to360(LonC1)) %>% 
  dplyr::select(c(lat,lon360,BSHn,Flag)) %>% 
  gather(species,n,-c(lat,lon360,Flag)) %>% 
  group_by(lat,lon360,Flag) %>% 
  summarise(number=sum(n,na.rm=T)) %>% 
  mutate(source="iattc_bycatch_ps")

master_rfmo_sharks_flag=do.call("rbind",list(iattc_ll,iattc_ps)) %>% 
  ungroup() %>% 
  mutate(cells=as.factor(raster::extract(template360,.[,2:1]))) %>% 
  group_by(cells) %>% 
  summarise(number=sum(number,na.rm=T)) %>% 
  ungroup() %>% 
  left_join(.,coords360) %>% 
  filter(cell_lat<62 & cell_lat>10)%>% 
  filter(cell_lon<260 & cell_lon>210)

mex_sharks_historical=ggplot() +
  geom_rect(data=df210,color = "black",fill="white",
            aes(xmin = xmin,
                xmax = xmax,
                ymin = ymin,
                ymax = ymax),alpha=1)+
  geom_tile(data=master_rfmo_sharks_flag, aes(x=cell_lon, y=cell_lat,fill=number)) +
  scale_fill_gradientn("Number",colours = pals::parula(100),na.value=NA,
                         breaks=c(500,1000,1500),
                         labels=c(".5k","1k","1.5k"))+
  geom_segment(aes(x=210,xend=260,y=60,yend=60),color="grey",linetype="dashed")+
  geom_segment(aes(x=210,xend=260,y=40,yend=40),color="grey",linetype="dashed")+
  geom_segment(aes(x=210,xend=260,y=20,yend=20),color="grey",linetype="dashed")+
  geom_segment(aes(x=220,xend=220,y=10,yend=62),color="grey",linetype="dashed")+
  geom_segment(aes(x=240,xend=240,y=10,yend=62),color="grey",linetype="dashed")+
  geom_text(data=lat_lon,aes(x=x,y=y,label=label),color="darkgrey",size=2)+
  geom_polygon(data=fortify(b210),aes(x=long, y = lat, group=group),color="black",fill=NA)+
  geom_polygon(data=fortify(template360_SA_shp_clip),aes(x=long, y = lat, group=group),color="red",fill=NA)+
  coord_map("conic", lat0 = 30,xlim = c(210, 260), ylim = c(10,62))+
  geom_polygon(data=fortify(wrld210),aes(x=long, y = lat, group=group),color="black",fill="grey",alpha=1)+
  theme_void() +
  theme(legend.position = "bottom")

png(glue("{outdir}/mex_sharks_historical.png"),width=12,height=10.5,units='cm',res=400,type = "cairo")
par(ps=10)
par(mar=c(4,4,1,1))
par(cex=1)
print({mex_sharks_historical})
# gg_hm
dev.off()

# China sharks historical ####
iattc_ll=read.csv("/Users/heatherwelch/Dropbox/OLE/data_OLE/iattc_effort_data/PublicLLShark/PublicLLSharkNum.csv")%>% 
  filter(Flag=="CHN") %>% 
  # filter(Year>=2017&Year<=2021) %>% 
  mutate(lat=LatC5,lon360=to360(LonC5)) %>% 
  dplyr::select(c(lat,lon360,BSHn,MAKn,Flag))%>% 
  gather(species,n,-c(lat,lon360,Flag)) %>% 
  group_by(lat,lon360,Flag) %>% 
  summarise(number=sum(n,na.rm=T)) %>% 
  mutate(source="iattc_bycatch_ll")


iattc_ps=read.csv("/Users/heatherwelch/Dropbox/OLE/data_OLE/iattc_effort_data/PublicPSShark/PublicPSSharkFlag.csv") %>% 
  # filter(Year>=2017&Year<=2021) %>% 
  filter(Flag=="CHN") %>% 
  mutate(lat=LatC1,lon360=to360(LonC1)) %>% 
  dplyr::select(c(lat,lon360,BSHn,Flag)) %>% 
  gather(species,n,-c(lat,lon360,Flag)) %>% 
  group_by(lat,lon360,Flag) %>% 
  summarise(number=sum(n,na.rm=T)) %>% 
  mutate(source="iattc_bycatch_ps")

master_rfmo_sharks_flag=do.call("rbind",list(iattc_ll,iattc_ps)) %>% 
  ungroup() %>% 
  mutate(cells=as.factor(raster::extract(template360,.[,2:1]))) %>% 
  group_by(cells) %>% 
  summarise(number=sum(number,na.rm=T)) %>% 
  ungroup() %>% 
  left_join(.,coords360) %>% 
  filter(cell_lat<62 & cell_lat>10)%>% 
  filter(cell_lon<260 & cell_lon>210)

chn_sharks_historical=ggplot() +
  geom_rect(data=df210,color = "black",fill="white",
            aes(xmin = xmin,
                xmax = xmax,
                ymin = ymin,
                ymax = ymax),alpha=1)+
  geom_tile(data=master_rfmo_sharks_flag, aes(x=cell_lon, y=cell_lat,fill=number)) +
  scale_fill_gradientn("Number",colours = pals::parula(100),na.value=NA,
                       breaks=c(5000,10000),
                       labels=c("5k","10k"))+
  geom_segment(aes(x=210,xend=260,y=60,yend=60),color="grey",linetype="dashed")+
  geom_segment(aes(x=210,xend=260,y=40,yend=40),color="grey",linetype="dashed")+
  geom_segment(aes(x=210,xend=260,y=20,yend=20),color="grey",linetype="dashed")+
  geom_segment(aes(x=220,xend=220,y=10,yend=62),color="grey",linetype="dashed")+
  geom_segment(aes(x=240,xend=240,y=10,yend=62),color="grey",linetype="dashed")+
  geom_text(data=lat_lon,aes(x=x,y=y,label=label),color="darkgrey",size=2)+
  geom_polygon(data=fortify(b210),aes(x=long, y = lat, group=group),color="black",fill=NA)+
  geom_polygon(data=fortify(template360_SA_shp_clip),aes(x=long, y = lat, group=group),color="red",fill=NA)+
  coord_map("conic", lat0 = 30,xlim = c(210, 260), ylim = c(10,62))+
  geom_polygon(data=fortify(wrld210),aes(x=long, y = lat, group=group),color="black",fill="grey",alpha=1)+
  theme_void() +
  theme(legend.position = "bottom")

png(glue("{outdir}/chn_sharks_historical.png"),width=12,height=10.5,units='cm',res=400,type = "cairo")
par(ps=10)
par(mar=c(4,4,1,1))
par(cex=1)
print({chn_sharks_historical})
# gg_hm
dev.off()

## mexico by the numbers ####

iattc_bycatch_ll=read.csv("/Users/heatherwelch/Dropbox/OLE/data_OLE/iattc_effort_data/PublicLLShark/PublicLLSharkNum.csv")%>% 
  filter(Flag=="MEX") %>% 
  # filter(Year>=2017&Year<=2021) %>% 
  mutate(lat=LatC5,lon360=to360(LonC5)) %>% 
  dplyr::select(c(lat,lon360,BSHn,MAKn,Year))%>% 
  gather(species,n,-c(lat,lon360,Year)) %>% 
  group_by(lat,lon360,species,Year) %>% 
  summarise(number=sum(n,na.rm=T)) %>% 
  mutate(source="iattc_bycatch_ll")


iattc_bycatch_ps=read.csv("/Users/heatherwelch/Dropbox/OLE/data_OLE/iattc_effort_data/PublicPSShark/PublicPSSharkFlag.csv") %>% 
  filter(Flag=="MEX") %>%   
  mutate(lat=LatC1,lon360=to360(LonC1)) %>% 
  dplyr::select(c(lat,lon360,BSHn,Year)) %>% 
  gather(species,n,-c(lat,lon360,Year)) %>% 
  group_by(species,Year,lat,lon360) %>% 
  summarise(number=sum(n,na.rm=T)) %>% 
  mutate(source="iattc_bycatch_ps")

master_rfmo_sharks=do.call("rbind",list(iattc_bycatch_ll,iattc_bycatch_ps)) %>% 
  ungroup() %>% 
  mutate(cells=as.factor(raster::extract(template360,.[,c("lon360","lat")]))) %>% 
  left_join(.,coords360) %>% 
  filter(cell_lon>=210&cell_lon<=260) %>%  ## iattc only, can't deal w US vessels wcpfc
  filter(cell_lat>=10&cell_lat<=62) %>% 
  group_by(cells,Year,source,species) %>% 
  summarise(number=sum(number,na.rm=T)) %>% 
  ungroup() %>% 
  left_join(.,coords360) 

m2=master_rfmo_sharks %>% 
  group_by(source,Year) %>% 
  summarise(n=sum(number,na.rm=T))

## china by the numbers ####

iattc_bycatch_ll=read.csv("/Users/heatherwelch/Dropbox/OLE/data_OLE/iattc_effort_data/PublicLLShark/PublicLLSharkNum.csv")%>% 
  filter(Flag=="CHN") %>% 
  # filter(Year>=2017&Year<=2021) %>% 
  mutate(lat=LatC5,lon360=to360(LonC5)) %>% 
  dplyr::select(c(lat,lon360,BSHn,MAKn,Year))%>% 
  gather(species,n,-c(lat,lon360,Year)) %>% 
  group_by(lat,lon360,species,Year) %>% 
  summarise(number=sum(n,na.rm=T)) %>% 
  mutate(source="iattc_bycatch_ll")


iattc_bycatch_ps=read.csv("/Users/heatherwelch/Dropbox/OLE/data_OLE/iattc_effort_data/PublicPSShark/PublicPSSharkFlag.csv") %>% 
  filter(Flag=="CHN") %>%   
  mutate(lat=LatC1,lon360=to360(LonC1)) %>% 
  dplyr::select(c(lat,lon360,BSHn,Year)) %>% 
  gather(species,n,-c(lat,lon360,Year)) %>% 
  group_by(species,Year,lat,lon360) %>% 
  summarise(number=sum(n,na.rm=T)) %>% 
  mutate(source="iattc_bycatch_ps")

master_rfmo_sharks=do.call("rbind",list(iattc_bycatch_ll,iattc_bycatch_ps)) %>% 
  ungroup() %>% 
  mutate(cells=as.factor(raster::extract(template360,.[,c("lon360","lat")]))) %>% 
  left_join(.,coords360) %>% 
  filter(cell_lon>=210&cell_lon<=260) %>%  ## iattc only, can't deal w US vessels wcpfc
  filter(cell_lat>=10&cell_lat<=62) %>% 
  group_by(cells,Year,source,species) %>% 
  summarise(number=sum(number,na.rm=T)) %>% 
  ungroup() %>% 
  left_join(.,coords360) 

m2_chn=master_rfmo_sharks %>% 
  group_by(source,Year) %>% 
  summarise(n=sum(number,na.rm=T))