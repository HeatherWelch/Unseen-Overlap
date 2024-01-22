### two case-studies: leatherback and laysan

library(patchwork)
library(nngeo)
library(sf)
library(scales)
library(paletteer)
library(circlize)

source("/Users/HeatherWelch/Dropbox/OLE/github/OLE_Projects/utilities/load_libraries.R")
boundingboxdir="UserDefined" # minimum bounding polygon to constrain species predictions
lbst_bb=glue("{boundingboxdir}/leatherbackTurtle_TOPP.shp") %>% st_read() %>% as_Spatial()
layalb_bb=glue("{boundingboxdir}/laysanAlbatross_TOPP.shp") %>% st_read() %>% as_Spatial()

## load in data ####
datdir="UserDefined/cleaned_overlap_aggregated_dfs"
master_join=read.csv(glue("{datdir}/cleaned_metrics_03-23-23.csv"))
glimpse(master_join)
masterAIS=read.csv(glue("{datdir}/AIS_metrics_03-23-23.csv"))

outdir="UserDefined";dir.create(outdir)

## load in spatial data ####
template=raster("/Users/HeatherWelch/Dropbox/OLE/spatial_data/template.grd") ## this is a user defined raster to assign fishing effort / gaps to
template_1deg=template
res(template_1deg)=1
nums=seq(1,ncell(template_1deg))
template_1deg=setValues(template_1deg,nums)
template_1deg_coords=rasterToPoints(template_1deg) %>% 
  as.data.frame() %>% rename(cell_lat=y,cell_lon=x,cells=Sea.level.anomaly)

sf::sf_use_s2(FALSE)
e <- as(extent(to360(-180), to360(-100), 10, 62), 'SpatialPolygons')
var="sst"
data=maps::map("world2",fill=T)
IDs <- sapply(strsplit(data$names, ":"), function(x) x[1])
wrld_simpl <- map2SpatialPolygons(data, IDs=IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))
wrld=SpatialPolygons(wrld_simpl@polygons,proj4string=wrld_simpl@proj4string) %>% 
  gBuffer(., byid=TRUE, width=0)
wrld2=st_as_sf(wrld) %>% st_crop( xmin=180, xmax=260, ymax=62, ymin=10) %>% as_Spatial()

df=data.frame(xmin = 180,
              xmax = 260,
              ymin = 10,
              ymax = 62)

eez=st_read("/Users/heatherwelch/Dropbox/woody/World_EEZ_v11_20191118_HR_0_360\ 2/eez_v11_0_360.shp") 
eez2=eez%>% 
  st_simplify(preserveTopology=TRUE, dTolerance = .2)

a=as_Spatial(eez)
b=crop(a,e)


## case_study panels####
case_study=master_join %>% 
  filter(!is.na(species)) %>%
  mutate(overlap_hours_in_gaps_over_12=hours_in_gaps_over_12*hab.suit,
         overlap_hours_in_gaps_under_12=hours_in_gaps_under_12*hab.suit,
         overlap_total_gap_2w=total_gap_2w*hab.suit) 


### lbst all flags ####
lbst_all=case_study %>% filter(species=="leatherbackTurtle_TOPP") %>% 
  dplyr::select(c(cell_lon,cell_lat,overlap_hours_in_gaps_over_12,overlap_hours_in_gaps_under_12,
                  overlap_total_gap_2w,overlap_hours_in_gaps_under_12)) %>% 
  group_by(cell_lon,cell_lat) %>% 
  summarise_all(sum,na.rm=T) %>% 
  ungroup() %>% 
  
  ## overlap metrics
  mutate(overlap_observed=overlap_hours_in_gaps_over_12+overlap_hours_in_gaps_under_12) %>%  ## observed ais
  mutate(overlap_total_2w_activity=overlap_total_gap_2w+overlap_hours_in_gaps_under_12) %>% # 2w observed + gaps
  mutate(overlap_perc_increase_2w_gaps=(overlap_total_2w_activity-overlap_observed)/overlap_total_2w_activity) %>% 
  
  ## stuff for plotting
  mutate(overlap_perc_increase_2w_gaps_rec=case_when(overlap_perc_increase_2w_gaps<0~0,
                                                     overlap_perc_increase_2w_gaps>.20~.20,
                                                     T~overlap_perc_increase_2w_gaps))%>%  
  mutate(overlap_perc_increase_2w_gaps_rec_round=as.character(round(overlap_perc_increase_2w_gaps_rec,3))) %>% 
  
  mutate(overlap_total_2w_activity_scale=rescale(overlap_total_2w_activity,c(.00001,1)))%>% 
  mutate(overlap_total_2w_activity_rec=case_when(overlap_total_2w_activity_scale>.002~.002,
                                                 T~overlap_total_2w_activity_scale))%>%  
  mutate(overlap_total_2w_activity_rec_round=as.character(round(overlap_total_2w_activity_rec,4))) %>% 
  .[complete.cases(.),]

d <- expand.grid(overlap_perc_increase_2w_gaps_rec_round = seq(0, .20, 0.001), overlap_total_2w_activity_rec_round = seq(0, .002, 0.0001)) %>%
  mutate(fill_val = overlap_perc_increase_2w_gaps_rec_round,
         transparency =  overlap_total_2w_activity_rec_round)

p=ggplot(d, aes(x=overlap_perc_increase_2w_gaps_rec_round, y=overlap_total_2w_activity_rec_round, fill = fill_val,alpha=transparency)) +
  geom_tile() +
  scale_fill_gradientn(colours = pals::parula(nrow(d)),na.value="black")+
  theme_classic() +
  theme(legend.position = "none")

pal_d <- ggplot_build(p)$data[[1]] %>%
  dplyr::select(x, y, fill,alpha) %>%
  mutate(overlap_perc_increase_2w_gaps_rec_round = as.character(x),
         overlap_total_2w_activity_rec_round = as.character(y)) %>% 
  dplyr::select(-c(x,y))

lbst_all2=lbst_all %>% 
  dplyr::select(c(cell_lon,cell_lat,overlap_perc_increase_2w_gaps_rec_round,overlap_total_2w_activity_rec_round)) %>% 
  left_join(.,pal_d)

lat_lon_full=data.frame(label=c("60°N      ","40°N     ","20°N    ","200°","220°","240°"),
                        x=c(178,178,178,200,220,240),
                        y=c(60,40,20,8.5,8.5,8.5))

lbst_plot=ggplot() +
  geom_tile(data=lbst_all2, aes(x=cell_lon, y=cell_lat,fill=fill,alpha=alpha)) +
  scale_fill_identity() +
  geom_segment(aes(x=180,xend=260,y=60,yend=60),color="grey",linetype="dashed")+
  geom_segment(aes(x=180,xend=260,y=40,yend=40),color="grey",linetype="dashed")+
  geom_segment(aes(x=180,xend=260,y=20,yend=20),color="grey",linetype="dashed")+
  geom_segment(aes(x=200,xend=200,y=10,yend=62),color="grey",linetype="dashed")+
  geom_segment(aes(x=220,xend=220,y=10,yend=62),color="grey",linetype="dashed")+
  geom_segment(aes(x=240,xend=240,y=10,yend=62),color="grey",linetype="dashed")+
  geom_text(data=lat_lon_full,aes(x=x,y=y,label=label),color="darkgrey",size=3)+
  # geom_polygon(data=fortify(lbst_bb),aes(x=long, y = lat, group=group),color="darkgrey",fill=NA)+
  geom_path(data=fortify(b),aes(x=long, y = lat, group=group),color="black",fill=NA,size=.5)+
  coord_map("conic", lat0 = 30,xlim = c(180, 260), ylim = c(10,62))+
  geom_polygon(data=fortify(wrld2),aes(x=long, y = lat, group=group),color="#688f74",fill="grey",alpha=1,size=.2)+
  geom_rect(data=df,color = "black",fill=NA,
            aes(xmin = xmin,
                xmax = xmax,
                ymin = ymin,
                ymax = ymax),alpha=1)+
  theme_void() +
  theme(legend.position = "none")

png(glue("{outdir}/Fig3_lbst.png"),width=12,height=10.5,units='cm',res=400,type = "cairo")
par(ps=10)
par(mar=c(4,4,1,1))
par(cex=1)
print({lbst_plot})
# gg_hm
dev.off()

### laysan all flags ####
laysab_all=case_study %>% filter(species=="laysanAlbatross_Dallas") %>% 
  dplyr::select(c(cell_lon,cell_lat,overlap_hours_in_gaps_over_12,overlap_hours_in_gaps_under_12,
                  overlap_total_gap_2w,overlap_hours_in_gaps_under_12)) %>% 
  group_by(cell_lon,cell_lat) %>% 
  summarise_all(sum,na.rm=T) %>% 
  ungroup() %>% 
  
  ## overlap metrics
  mutate(overlap_observed=overlap_hours_in_gaps_over_12+overlap_hours_in_gaps_under_12) %>%  ## observed ais
  mutate(overlap_total_2w_activity=overlap_total_gap_2w+overlap_hours_in_gaps_under_12) %>% # 2w observed + gaps
  mutate(overlap_perc_increase_2w_gaps=(overlap_total_2w_activity-overlap_observed)/overlap_total_2w_activity) %>% 
  
  ## stuff for plotting
  mutate(overlap_perc_increase_2w_gaps_rec=case_when(overlap_perc_increase_2w_gaps<0~0,
                                                     overlap_perc_increase_2w_gaps>.20~.20,
                                                     T~overlap_perc_increase_2w_gaps))%>%  
  mutate(overlap_perc_increase_2w_gaps_rec_round=as.character(round(overlap_perc_increase_2w_gaps_rec,3))) %>% 
  
  mutate(overlap_total_2w_activity_scale=rescale(overlap_total_2w_activity,c(.00001,1)))%>% 
  mutate(overlap_total_2w_activity_rec=case_when(overlap_total_2w_activity_scale>.002~.002,
                                                 T~overlap_total_2w_activity_scale))%>%  
  mutate(overlap_total_2w_activity_rec_round=as.character(round(overlap_total_2w_activity_rec,4))) %>% 
  .[complete.cases(.),]

d <- expand.grid(overlap_perc_increase_2w_gaps_rec_round = seq(0, .20, 0.001), overlap_total_2w_activity_rec_round = seq(0, .002, 0.0001)) %>%
  mutate(fill_val = overlap_perc_increase_2w_gaps_rec_round,
         transparency =  overlap_total_2w_activity_rec_round)

p=ggplot(d, aes(x=overlap_perc_increase_2w_gaps_rec_round, y=overlap_total_2w_activity_rec_round, fill = fill_val,alpha=transparency)) +
  geom_tile() +
  scale_fill_gradientn(colours = pals::parula(nrow(d)),na.value="black")+
  theme_classic() +
  theme(legend.position = "none")

pal_d <- ggplot_build(p)$data[[1]] %>%
  dplyr::select(x, y, fill,alpha) %>%
  mutate(overlap_perc_increase_2w_gaps_rec_round = as.character(x),
         overlap_total_2w_activity_rec_round = as.character(y)) %>% 
  dplyr::select(-c(x,y))

laysab_all2=laysab_all %>% 
  dplyr::select(c(cell_lon,cell_lat,overlap_perc_increase_2w_gaps_rec_round,overlap_total_2w_activity_rec_round)) %>% 
  left_join(.,pal_d)


laysab_plot=ggplot() +
  geom_tile(data=laysab_all2, aes(x=cell_lon, y=cell_lat,fill=fill,alpha=alpha)) +
  scale_fill_identity() +
  geom_segment(aes(x=180,xend=260,y=60,yend=60),color="grey",linetype="dashed")+
  geom_segment(aes(x=180,xend=260,y=40,yend=40),color="grey",linetype="dashed")+
  geom_segment(aes(x=180,xend=260,y=20,yend=20),color="grey",linetype="dashed")+
  geom_segment(aes(x=200,xend=200,y=10,yend=62),color="grey",linetype="dashed")+
  geom_segment(aes(x=220,xend=220,y=10,yend=62),color="grey",linetype="dashed")+
  geom_segment(aes(x=240,xend=240,y=10,yend=62),color="grey",linetype="dashed")+
  geom_text(data=lat_lon_full,aes(x=x,y=y,label=label),color="darkgrey",size=3)+
  # geom_polygon(data=fortify(lbst_bb),aes(x=long, y = lat, group=group),color="darkgrey",fill=NA)+
  geom_path(data=fortify(b),aes(x=long, y = lat, group=group),color="black",fill=NA,size=.5)+
  coord_map("conic", lat0 = 30,xlim = c(180, 260), ylim = c(10,62))+
  geom_polygon(data=fortify(wrld2),aes(x=long, y = lat, group=group),color="#688f74",fill="grey",alpha=1,size=.2)+
  geom_rect(data=df,color = "black",fill=NA,
            aes(xmin = xmin,
                xmax = xmax,
                ymin = ymin,
                ymax = ymax),alpha=1)+
  theme_void() +
  theme(legend.position = "none")

png(glue("{outdir}/Fig3_laysab.png"),width=12,height=10.5,units='cm',res=400,type = "cairo")
par(ps=10)
par(mar=c(4,4,1,1))
par(cex=1)
print({laysab_plot})
# gg_hm
dev.off()

### legend  ####
p2=ggplot(d, aes(x=overlap_perc_increase_2w_gaps_rec_round, 
                 y=rescale(overlap_total_2w_activity_rec_round,c(0,1)), fill = fill_val,alpha=transparency)) +
  geom_tile() +
  scale_fill_gradientn(colours = pals::parula(100),na.value="black")+
  theme_classic() +
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0,0),breaks=c(0,.05,.10,.15,.20),labels=c("0%","5%","10%","15%",">20%"))+
  scale_y_continuous(expand = c(0,0),breaks=c(0,1),labels = c("Lower","Higher"))+
  xlab(NULL)+
  ylab(NULL)+
  theme(legend.position = "none",
        plot.margin = margin(.5,.5,.5,.5, "cm")
  )

png(glue("{outdir}/Fig3_legend.png"),width=5.5,height=4.5,units='cm',res=400,type = "cairo")
par(ps=10)
par(mar=c(4,4,4,4))
par(cex=1)
print({p2})
# gg_hm
dev.off()

## percents by location/flag ####
eez2=eez %>% 
  filter(TERRITORY1=="Alaska"|
           TERRITORY1=="Hawaii"|
           TERRITORY1=="Mexico"|
           TERRITORY1=="Canada"|
           TERRITORY1=="United States")
eez3=eez2 %>% st_remove_holes(.) %>% 
  st_buffer(dist=.2)

derived_map_sub=master_join %>% 
  filter(!is.na(species)) %>%
 filter(species=="leatherbackTurtle_TOPP"|species=="laysanAlbatross_Dallas") %>% 
  mutate(overlap_hours_in_gaps_over_12=hours_in_gaps_over_12*hab.suit,
         overlap_hours_in_gaps_under_12=hours_in_gaps_under_12*hab.suit,
         overlap_total_gap_2w=total_gap_2w*hab.suit,
         overlap_total_gap_disabling_2w=total_gap_disabling_2w*hab.suit) %>%
  dplyr::select(-c(cells,date,vessel_class,hab.suit)) %>% 
  group_by(cell_lat,cell_lon,flag,species) %>% 
  summarise_all(sum,na.rm=T) %>% 
  ungroup() %>% 
  
  ## overlap metrics
  mutate(overlap_observed=overlap_hours_in_gaps_over_12+overlap_hours_in_gaps_under_12) %>%  ## observed ais
  mutate(overlap_total_2w_activity=overlap_total_gap_2w+overlap_hours_in_gaps_under_12) # 2w observed + gaps

derived_map_sub_eez=derived_map_sub %>% 
  mutate(y=cell_lat,x=cell_lon)
coordinates(derived_map_sub_eez)=~x+y
crs(derived_map_sub_eez)=crs(eez)
derived_map_sub_eez2=st_as_sf(derived_map_sub_eez)

derived_map_sub_eez2_join <- st_join(derived_map_sub_eez2, eez3) %>% as.data.frame() %>% 
  dplyr::select(-geometry) %>% 
  mutate(TERRITORY1=case_when(is.na(TERRITORY1)~"High seas",
                              T~TERRITORY1)) %>% 
  dplyr::select(species,flag,TERRITORY1,overlap_observed,overlap_total_2w_activity,overlap_total_gap_2w)

full_domain_all_flags=derived_map_sub_eez2_join %>% 
  group_by(species,flag) %>% 
  summarise(overlap_total_gap_2w=sum(overlap_total_gap_2w,na.rm = T)) %>% 
  mutate(total=sum(overlap_total_gap_2w,na.rm=T)) %>% 
  mutate(perc=percent(overlap_total_gap_2w/total,1)) %>% 
  filter(perc!="0%")

highseas_domain_all_flags=derived_map_sub_eez2_join %>% 
  filter(TERRITORY1=="High seas") %>% 
  group_by(species,flag) %>% 
  summarise(overlap_total_gap_2w=sum(overlap_total_gap_2w,na.rm = T)) %>% 
  mutate(total=sum(overlap_total_gap_2w,na.rm=T)) %>% 
  mutate(perc=percent(overlap_total_gap_2w/total,1)) %>% 
  filter(perc!="0%")

not_highseas_domain_all_flags=derived_map_sub_eez2_join %>% 
  filter(TERRITORY1!="High seas") %>% 
  group_by(species,flag) %>% 
  summarise(overlap_total_gap_2w=sum(overlap_total_gap_2w,na.rm = T)) %>% 
  mutate(total=sum(overlap_total_gap_2w,na.rm=T)) %>% 
  mutate(perc=percent(overlap_total_gap_2w/total,1)) %>% 
  filter(perc!="0%")

domain_no_flag=derived_map_sub_eez2_join %>% 
  # filter(TERRITORY1!="High seas") %>% 
  group_by(species,TERRITORY1) %>% 
  summarise(overlap_total_gap_2w=sum(overlap_total_gap_2w,na.rm = T)) %>% 
  mutate(total=sum(overlap_total_gap_2w,na.rm=T)) %>% 
  mutate(perc=percent(overlap_total_gap_2w/total,1)) %>% 
  filter(perc!="0%")

## stats to report in paper ####
stats=derived_map_sub_eez2_join %>% 
  mutate(TERRITORY1=case_when(TERRITORY1=="Hawaii"~"US EEZ",
                              TERRITORY1=="Alaska"~"US EEZ",
                              TERRITORY1=="United States"~"US EEZ",
                              TERRITORY1=="Mexico"~"MEX EEZ",
                              TERRITORY1=="High seas"~"High seas",
                              TERRITORY1=="Canada"~"CAN EEZ",
                              T~TERRITORY1)) %>% 
  filter(species=="leatherbackTurtle_TOPP"|species=="laysanAlbatross_Dallas") %>% 
  group_by(species,TERRITORY1,flag) %>% 
  summarise(overlap_total_gap_2w=sum(overlap_total_gap_2w,na.rm = T)) %>% 
  ungroup() %>%
  mutate(total=sum(overlap_total_gap_2w,na.rm=T)) %>% 
  mutate(perc=percent(overlap_total_gap_2w/total,1)) %>% 
  filter(perc!="0%")

lbst_stats=stats %>% 
  mutate(TERRITORY1=case_when(TERRITORY1=="Hawaii"~"US EEZ",
                              TERRITORY1=="Alaska"~"US EEZ",
                              TERRITORY1=="United States"~"US EEZ",
                              TERRITORY1=="Mexico"~"MEX EEZ",
                              TERRITORY1=="High seas"~"High seas",
                              TERRITORY1=="Canada"~"CAN EEZ",
                              T~TERRITORY1)) %>% 
  filter(species=="leatherbackTurtle_TOPP") %>% 
  group_by(species,TERRITORY1,flag) %>% 
  summarise(overlap_total_gap_2w=sum(overlap_total_gap_2w,na.rm = T)) %>% 
  ungroup() %>%
  mutate(total=sum(overlap_total_gap_2w,na.rm=T)) %>% 
  mutate(perc=percent(overlap_total_gap_2w/total,1)) %>% 
  filter(perc!="0%")

lbst_stats_US_EEZ=stats %>% 
  mutate(TERRITORY1=case_when(TERRITORY1=="Hawaii"~"US EEZ",
                              TERRITORY1=="Alaska"~"US EEZ",
                              TERRITORY1=="United States"~"US EEZ",
                              TERRITORY1=="Mexico"~"MEX EEZ",
                              TERRITORY1=="High seas"~"High seas",
                              TERRITORY1=="Canada"~"CAN EEZ",
                              T~TERRITORY1)) %>% 
  filter(species=="leatherbackTurtle_TOPP") %>% 
  filter(TERRITORY1=="US EEZ") %>%
  group_by(species,TERRITORY1,flag) %>% 
  summarise(overlap_total_gap_2w=sum(overlap_total_gap_2w,na.rm = T)) %>% 
  ungroup() %>%
  mutate(total=sum(overlap_total_gap_2w,na.rm=T)) %>% 
  mutate(perc=percent(overlap_total_gap_2w/total,1)) %>% 
  filter(perc!="0%")

lbst_stats_highseas=stats %>% 
  mutate(TERRITORY1=case_when(TERRITORY1=="Hawaii"~"US EEZ",
                              TERRITORY1=="Alaska"~"US EEZ",
                              TERRITORY1=="United States"~"US EEZ",
                              TERRITORY1=="Mexico"~"MEX EEZ",
                              TERRITORY1=="High seas"~"High seas",
                              TERRITORY1=="Canada"~"CAN EEZ",
                              T~TERRITORY1)) %>% 
  filter(species=="leatherbackTurtle_TOPP") %>% 
  filter(TERRITORY1=="High seas") %>%
  group_by(species,TERRITORY1,flag) %>% 
  summarise(overlap_total_gap_2w=sum(overlap_total_gap_2w,na.rm = T)) %>% 
  ungroup() %>%
  mutate(total=sum(overlap_total_gap_2w,na.rm=T)) %>% 
  mutate(perc=percent(overlap_total_gap_2w/total,1)) %>% 
  filter(perc!="0%")

layalb_stats=stats %>% 
  mutate(TERRITORY1=case_when(TERRITORY1=="Hawaii"~"US EEZ",
                              TERRITORY1=="Alaska"~"US EEZ",
                              TERRITORY1=="United States"~"US EEZ",
                              TERRITORY1=="Mexico"~"MEX EEZ",
                              TERRITORY1=="High seas"~"High seas",
                              TERRITORY1=="Canada"~"CAN EEZ",
                              T~TERRITORY1)) %>% 
  filter(species=="laysanAlbatross_Dallas") %>% 
  group_by(species,TERRITORY1,flag) %>% 
  summarise(overlap_total_gap_2w=sum(overlap_total_gap_2w,na.rm = T)) %>% 
  ungroup() %>%
  mutate(total=sum(overlap_total_gap_2w,na.rm=T)) %>% 
  mutate(perc=percent(overlap_total_gap_2w/total,1)) %>% 
  filter(perc!="0%")

layalb_stats_highseas=stats %>% 
  mutate(TERRITORY1=case_when(TERRITORY1=="Hawaii"~"US EEZ",
                              TERRITORY1=="Alaska"~"US EEZ",
                              TERRITORY1=="United States"~"US EEZ",
                              TERRITORY1=="Mexico"~"MEX EEZ",
                              TERRITORY1=="High seas"~"High seas",
                              TERRITORY1=="Canada"~"CAN EEZ",
                              T~TERRITORY1)) %>% 
  filter(species=="laysanAlbatross_Dallas") %>% 
  filter(TERRITORY1=="High seas") %>%
  group_by(species,TERRITORY1,flag) %>% 
  summarise(overlap_total_gap_2w=sum(overlap_total_gap_2w,na.rm = T)) %>% 
  ungroup() %>%
  mutate(total=sum(overlap_total_gap_2w,na.rm=T)) %>% 
  mutate(perc=percent(overlap_total_gap_2w/total,1)) %>% 
  filter(perc!="0%")

### leatherback ####
test_chord=derived_map_sub_eez2_join %>% 
  filter(species=="leatherbackTurtle_TOPP") %>%
  mutate(TERRITORY1=case_when(TERRITORY1=="Hawaii"~"US EEZ",
                              TERRITORY1=="United States"~"US EEZ",
                              TERRITORY1=="Mexico"~"MEX EEZ",
                              TERRITORY1=="High seas"~"High seas",
                              TERRITORY1=="Canada"~"CAN EEZ",
                              T~TERRITORY1)) %>% 
  mutate(flag=case_when(flag=="USA"~"US",
                        flag=="MEX"~"MEX",
                        flag=="CAN"~"CAN",
                       
                        T~"DWF")) %>% 
  group_by(species,TERRITORY1,flag) %>% 
  summarise(overlap_total_gap_2w=sum(overlap_total_gap_2w,na.rm = T)) %>% 
  ungroup() %>% 
  filter(overlap_total_gap_2w>3) %>%
  dplyr::select(flag,TERRITORY1,overlap_total_gap_2w) 
 
sectors_chord=list(unique(test_chord$flag),unique(test_chord$TERRITORY1)) %>%
  unlist()  %>% 
  as.data.frame() %>%
  .[complete.cases(.),] %>% 
  as.data.frame() %>%
 rename("Geo"=".") %>% 
  mutate(colors=c(rep("#EFC000FF",length(unique(test_chord$flag))),rep("#0073C2FF",length(unique(test_chord$TERRITORY1)))))

arrow_colors=test_chord %>% 
  mutate(color=case_when(flag=="DWF"~"#A73030FF",
                         flag=="CHN"~"#3B3B3BFF",
                         # flag=="DWF"~"#EFC000FF",
                         flag=="US"~"#868686FF",
                         flag=="MEX"~"#458553",
                         flag=="CAN"~"#4A6990FF")) %>% 
  pull(color)

alpha=test_chord %>% 
  mutate(alpha=case_when(flag=="CAN"~.5,
                         flag=="MEX"~.1,
                         flag=="DWF"~.1,
                         flag=="VUT"~.1,
                         flag=="US"~.8,
                         flag=="TWN"~.1)) %>% 
  pull(alpha)

 mycolors_chord=sectors_chord$colors

png(glue("{outdir}/chrod_lbst.png"),width=10,height=10,units='cm',res=400)
par(ps=10)

par(oma=c(.2,.2,.2,.2))
par(cex=1)
circos.clear()
 
a=chordDiagram(
  x = test_chord, 
  grid.col = mycolors_chord,
  col=arrow_colors,
  transparency = alpha,
  directional = 1,
  direction.type = c("arrows"), 
  annotationTrack = "grid",
  annotationTrackHeight = c(0.05, 0.1),
  link.arr.type = "big.arrow", 
  diffHeight = -uh(5, "mm"),
  link.sort = TRUE, 
  link.largest.ontop = TRUE,
  link.border = "black",
  big.gap = 12.0,
  small.gap = 2,
  preAllocateTracks = list(list(track.height = 0.1),list(bg.border = "white"))
)

circos.trackPlotRegion(
  track.index = 1, 
  bg.border = NA, 
  panel.fun = function(x, y) {
    
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.index = get.cell.meta.data("sector.index") %>% gsub("_T","",.)%>% gsub("_F","",.) 
    
    # Add names to the sector. 
    circos.text(
      x = mean(xlim), 
      y=ylim[1] - cm_h(9), 
      labels = sector.index, 
      facing = "clockwise", 
      niceFacing = TRUE,
      adj = c(0, .5),
      cex = .7
    )
    
    # Add graduation on axis
  }
)

dev.off()

### laysan ####
test_chord=derived_map_sub_eez2_join %>% 
  filter(species=="laysanAlbatross_Dallas") %>% 
  mutate(TERRITORY1=case_when(TERRITORY1=="Hawaii"~"US EEZ",
                              TERRITORY1=="Alaska"~"US EEZ",
                              TERRITORY1=="United States"~"USA EEZ",
                              TERRITORY1=="Mexico"~"MEX EEZ",
                              TERRITORY1=="High seas"~"High seas",
                              TERRITORY1=="Canada"~"CAN EEZ",
                              T~TERRITORY1)) %>% 
  mutate(flag=case_when(flag=="USA"~"US",
                        flag=="MEX"~"MEX",
                        flag=="CAN"~"CAN",
                        T~"DWF")) %>% 
  group_by(species,TERRITORY1,flag) %>% 
  summarise(overlap_total_gap_2w=sum(overlap_total_gap_2w,na.rm = T)) %>% 
  ungroup() %>% 
  filter(overlap_total_gap_2w>3) %>%
  dplyr::select(flag,TERRITORY1,overlap_total_gap_2w) 
  
sectors_chord=list(unique(test_chord$flag),unique(test_chord$TERRITORY1)) %>%
  unlist()  %>% 
  as.data.frame() %>%
  .[complete.cases(.),] %>% 
  as.data.frame() %>%
  rename("Geo"=".") %>% 
  mutate(colors=c(rep("#EFC000FF",length(unique(test_chord$flag))),rep("#0073C2FF",length(unique(test_chord$TERRITORY1)))))

arrow_colors=test_chord %>% 
  mutate(color=case_when(flag=="DWF"~"#A73030FF",
                         flag=="CHN"~"#3B3B3BFF",
                         # flag=="DWF"~"#EFC000FF",
                         flag=="US"~"#868686FF",
                         flag=="MEX"~"#458553",
                         flag=="CAN"~"#4A6990FF")) %>% 
  pull(color)

alpha=test_chord %>% 
  mutate(alpha=case_when(flag=="CAN"~.5,
                         flag=="MEX"~.1,
                         flag=="DWF"~.1,
                         flag=="VUT"~.1,
                         flag=="US"~.8,
                         flag=="TWN"~.1)) %>% 
  pull(alpha)

labels=list(unique(test_chord$flag),unique(test_chord$TERRITORY1)) %>%
  unlist()  %>% 
  as.data.frame() %>% 
  rename("Geo"=".") %>%
  mutate(length=nchar(Geo)) %>% 
  mutate(lenght2=rescale(length,c(-1.5,0))) %>% 
  pull(lenght2)

mycolors_chord=sectors_chord$colors

png(glue("{outdir}/chrod_layalb.png"),width=10,height=10,units='cm',res=400)
par(ps=10)
par(oma=c(.2,.2,.2,.2))
par(cex=1)
circos.clear()

a=chordDiagram(
  x = test_chord, 
   grid.col = mycolors_chord,
   col=arrow_colors,
  transparency = alpha,
  directional = 1,
  direction.type = c("arrows"), 
  annotationTrack = "grid",
  annotationTrackHeight = c(0.05, 0.1),
  link.arr.type = "big.arrow", 
  diffHeight = -uh(5, "mm"),
  link.sort = TRUE, 
  link.largest.ontop = TRUE,
  link.border = "black",
  big.gap = 12.0,
  preAllocateTracks = list(list(track.height = 0.1),list(bg.border = "white"))
)

circos.trackPlotRegion(
  track.index = 1, 
  bg.border = NA, 
  panel.fun = function(x, y) {
    
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.index = get.cell.meta.data("sector.index") %>% gsub("_T","",.)%>% gsub("_F","",.) 
    
    # Add names to the sector. 
    circos.text(
      x = mean(xlim), 
      y=ylim[1] - cm_h(9), 
      labels = sector.index, 
      facing = "clockwise", 
      niceFacing = TRUE,
      adj = c(0, .5),
      cex = .7
    )
    
    # Add graduation on axis
  }
)


dev.off()

