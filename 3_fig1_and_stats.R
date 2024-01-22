## calculating stats

library(patchwork)
library(nngeo)
library(sf)
source("/Users/HeatherWelch/Dropbox/OLE/github/OLE_Projects/utilities/load_libraries.R")

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

## total fishing vessel activity metrics ####
FA_totals=masterAIS %>% 
  dplyr::select(-c(cells,date,vessel_class,flag,X)) %>% 
  summarise_all(sum,na.rm=T) %>% 
  mutate(observed=hours_in_gaps_over_12+hours_in_gaps_under_12) %>%  ## observed ais
  mutate(total_activity=total_gap+hours_in_gaps_under_12) %>%  ## full observed + gaps
  mutate(total_2w_activity=total_gap_2w+hours_in_gaps_under_12) %>% # 2w observed + gaps
  mutate(perc_gaps=observed/total_activity) %>%  ## what percent of total activity is observed
  mutate(perc_gaps_2w=observed/total_2w_activity) %>%  ## what percent of total activity is observed 2w
  mutate(perc_unseen=total_gap_disabling/total_activity, ## full what percent of total activity is disabling
         perc_unseen_2w=total_gap_disabling_2w/total_2w_activity)  ## full what percent of total activity is disabling 2w

FA_totals_flag=masterAIS %>% 
  dplyr::select(-c(cells,date,vessel_class,X)) %>% 
  group_by(flag) %>% 
  summarise_all(sum,na.rm=T) %>% 
  mutate(observed=hours_in_gaps_over_12+hours_in_gaps_under_12) %>%  ## observed ais
  mutate(total_activity=total_gap+hours_in_gaps_under_12) %>%  ## full observed + gaps
  mutate(total_2w_activity=total_gap_2w+hours_in_gaps_under_12)  # 2w observed + gaps

arrange(FA_totals_flag,desc(total_activity)) %>% pull(flag)

FA_totals_vessel_class=masterAIS %>% 
  mutate(vessel_class=case_when(vessel_class=="other_purse_seines" | vessel_class=="other_seines" ~ "non_tuna_purse_seines",
                                T~vessel_class)) %>% 
  dplyr::select(-c(cells,date,flag,X)) %>% 
  group_by(vessel_class) %>% 
  summarise_all(sum,na.rm=T) %>% 
  mutate(observed=hours_in_gaps_over_12+hours_in_gaps_under_12) %>%  ## observed ais
  mutate(total_activity=total_gap+hours_in_gaps_under_12) %>%  ## full observed + gaps
  mutate(total_2w_activity=total_gap_2w+hours_in_gaps_under_12)  # 2w observed + gaps

arrange(FA_totals_vessel_class,desc(total_activity)) %>% pull(vessel_class)

## total overlap metrics ####
master_join_derived=master_join %>% 
  filter(!is.na(species)) %>% 
  mutate(overlap_hours_in_gaps_over_12=hours_in_gaps_over_12*hab.suit,
         overlap_hours_in_gaps_under_12=hours_in_gaps_under_12*hab.suit,
         overlap_total_gap=total_gap*hab.suit,
         overlap_total_gap_2w=total_gap_2w*hab.suit,
         overlap_total_gap_disabling=total_gap_disabling*hab.suit,
         overlap_total_gap_disabling_2w=total_gap_disabling_2w*hab.suit) %>% 
  dplyr::select(c(overlap_hours_in_gaps_over_12,overlap_hours_in_gaps_under_12,overlap_total_gap,overlap_total_gap_2w,overlap_total_gap_disabling,overlap_total_gap_disabling_2w)) %>% 
  # group_by(cell_lat,cell_lon,flag,vessel_class,species) %>% 
  summarise_all(.funs=sum,na.rm=T) %>% 
  mutate(overlap_observed=overlap_hours_in_gaps_over_12+overlap_hours_in_gaps_under_12) %>%  ## observed ais
  
  mutate(overlap_total_activity=overlap_total_gap+overlap_hours_in_gaps_under_12) %>%  ## full observed + gaps
  mutate(overlap_total_2w_activity=overlap_total_gap_2w+overlap_hours_in_gaps_under_12) %>% # 2w observed + gaps
  
  mutate(overlap_perc_gaps=overlap_observed/overlap_total_activity) %>%  ## what percent of total activity is observed
  mutate(overlap_perc_2w_gaps=overlap_observed/overlap_total_2w_activity) %>%  ## what percent of total activity is observed 2w
  
  mutate(overlap_perc_increase_gaps=(overlap_total_activity-overlap_observed)/overlap_total_activity) %>%  ## % increase
  mutate(overlap_perc_increase_2w_gaps=(overlap_total_2w_activity-overlap_observed)/overlap_total_2w_activity) %>%  ##  % increase 2w
  
  mutate(overlap_perc_unseen=overlap_total_gap_disabling/overlap_total_activity, ## full what percent of total activity is disabling
         overlap_perc_unseen_2w=overlap_total_gap_disabling_2w/overlap_total_2w_activity)  ## full what percent of total activity is disabling

## figure 1. oberved overlap vs total overlap 2w ####
derived_map=master_join %>% 
  filter(!is.na(species)) %>%
  mutate(overlap_hours_in_gaps_over_12=hours_in_gaps_over_12*hab.suit,
         overlap_hours_in_gaps_under_12=hours_in_gaps_under_12*hab.suit,
         overlap_total_gap_2w=total_gap_2w*hab.suit,
         overlap_total_gap_disabling_2w=total_gap_disabling_2w*hab.suit) %>%
  dplyr::select(-c(cells,date,vessel_class,flag,species,hab.suit)) %>% 
  group_by(cell_lat,cell_lon) %>% 
    summarise_all(sum,na.rm=T) %>% 
  ungroup() %>% 
    
  ## overlap metrics
  mutate(overlap_observed=overlap_hours_in_gaps_over_12+overlap_hours_in_gaps_under_12) %>%  ## observed ais
  mutate(overlap_total_2w_activity=overlap_total_gap_2w+overlap_hours_in_gaps_under_12) %>% # 2w observed + gaps
  mutate(overlap_perc_2w_gaps=overlap_observed/overlap_total_2w_activity) %>%  ## what percent of total activity is observed 2w
  mutate(overlap_perc_increase_2w_gaps=(overlap_total_2w_activity-overlap_observed)/overlap_total_2w_activity) %>%  ## what percent of total activity is observed 2w
  mutate(overlap_perc_unseen_2w=overlap_total_gap_disabling_2w/overlap_total_2w_activity) #%>%   ## full what percent of total activity is disabling
  
unseenVseen=derived_map %>% 
  dplyr::select(c(overlap_total_2w_activity,overlap_perc_increase_2w_gaps,cell_lat,cell_lon)) %>% 
                  
  mutate(overlap_perc_increase_2w_gaps_rec=case_when(overlap_perc_increase_2w_gaps<0~0,
                                                    overlap_perc_increase_2w_gaps>.25~.25,
                                                     T~overlap_perc_increase_2w_gaps))%>%  
  mutate(overlap_perc_increase_2w_gaps_rec_round=as.character(round(overlap_perc_increase_2w_gaps_rec,3))) %>% 
  
  mutate(overlap_total_2w_activity_rec=case_when(overlap_total_2w_activity>20~20,
                                                     T~overlap_total_2w_activity))%>%  
  mutate(overlap_total_2w_activity_rec_round=as.character(round(overlap_total_2w_activity_rec,1))) %>% 
  .[complete.cases(.),]
  
d <- expand.grid(overlap_perc_increase_2w_gaps_rec_round = seq(0, .25, 0.001), overlap_total_2w_activity_rec_round = seq(0, 20, .1)) %>%
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

unseenVseen2=unseenVseen %>% 
  dplyr::select(c(cell_lon,cell_lat,overlap_perc_increase_2w_gaps_rec_round,overlap_total_2w_activity_rec_round)) %>% 
  left_join(.,pal_d)

lat_lon_full=data.frame(label=c("60°N      ","40°N     ","20°N    ","200°","220°","240°"),
                        x=c(178,178,178,200,220,240),
                        y=c(60,40,20,8.5,8.5,8.5))


fig1=ggplot() +
  geom_tile(data=unseenVseen2, aes(x=cell_lon, y=cell_lat,fill=fill,alpha=alpha)) +
  scale_fill_identity() +
  geom_segment(aes(x=180,xend=260,y=60,yend=60),color="grey",linetype="dashed")+
  geom_segment(aes(x=180,xend=260,y=40,yend=40),color="grey",linetype="dashed")+
  geom_segment(aes(x=180,xend=260,y=20,yend=20),color="grey",linetype="dashed")+
  geom_segment(aes(x=200,xend=200,y=10,yend=62),color="grey",linetype="dashed")+
  geom_segment(aes(x=220,xend=220,y=10,yend=62),color="grey",linetype="dashed")+
  geom_segment(aes(x=240,xend=240,y=10,yend=62),color="grey",linetype="dashed")+
  geom_text(data=lat_lon_full,aes(x=x,y=y,label=label),color="darkgrey",size=3)+
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


png(glue("{outdir}/fig1_conical.png"),width=12,height=10.5,units='cm',res=400,type = "cairo")
par(ps=10)
par(mar=c(4,4,1,1))
par(cex=1)
print({fig1})
# gg_hm
dev.off()

p2=ggplot(d, aes(x=overlap_perc_increase_2w_gaps_rec_round, 
                 y=rescale(overlap_total_2w_activity_rec_round,c(0,1)), fill = fill_val,alpha=transparency)) +
  geom_tile() +
  scale_fill_gradientn(colours = pals::parula(100),na.value="black")+
  theme_classic() +
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0,0),breaks=c(0,.08,.16,.25),labels=c("0%","8%","16%",">25%"))+
  scale_y_continuous(breaks=c(0,1),labels = c("Lower","Higher"),expand = c(0,0))+
  xlab(NULL)+
  ylab(NULL)+
  theme(legend.position = "none",
        plot.margin = margin(.5,.5,.5,.5, "cm")
  )

png(glue("{outdir}/fig1_legend.png"),width=5.5,height=4.5,units='cm',res=400,type = "cairo")
par(ps=10)
par(mar=c(4,4,4,4))
par(cex=1)
print({p2})
# gg_hm
dev.off()

