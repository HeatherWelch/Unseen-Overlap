## combining spatially extrapolated gaps, fishing activity, and species distributions

source("/Users/HeatherWelch/Dropbox/OLE/github/OLE_Projects/utilities/load_libraries.R")

outdir="UserDefined";dir.create(outdir)

## random template to link fishing and gaps ## this is a user defined raster to assign fishing effort / gaps to
template=raster("/Users/HeatherWelch/Dropbox/OLE/spatial_data/template.grd")
template_1deg=template
res(template_1deg)=1
nums=seq(1,ncell(template_1deg))
template_1deg=setValues(template_1deg,nums)
template_1deg_coords=rasterToPoints(template_1deg) %>% 
  as.data.frame() %>% rename(cell_lat=y,cell_lon=x,cells=Sea.level.anomaly) %>% 
  mutate(cells=as.factor(cells))

## read in overlaps data ####
fishing=read.csv("UserDefined/species_fishery_extracto.csv")

fishing_clip=fishing %>% 
  filter(lon_bin>=to360(-179.9) & lon_bin<=to360(-100.1)) %>% 
  filter(lat_bin>=(10) & lat_bin<=(60))


gaps=read.csv("UserDefined/species_gaps_interpolated.csv") %>% 
  group_by(gap_id) %>% 
  mutate(n_cells=n()) %>% 
  mutate(gap_hours_per_cell=gap_hours/n_cells) %>% 
  mutate(gaps_2w=case_when(gap_hours>(24*14)~"no", ## gaps longer than 2 weeks (upper bound)
                           T~"yes")) %>%   ## gaps shorter than 2 weeks (lower bound)
  filter(cell_lat>=10 & cell_lat<60) %>%
  filter(cell_lon>=to360(-179.9) & cell_lon<=to360(-100.1)) %>% 
  rename(date=days)

# normalizing gap and fishing metrics together ####

fishing_clip_n=fishing_clip %>% 
  dplyr::select(c(lon_bin,lat_bin,hours_in_gaps_under_12,hours_in_gaps_over_12,date,vessel_class,flag)) %>% 
  mutate(cells=raster::extract(template_1deg,.[,1:2]))  %>% 
  dplyr::select(-c(lon_bin,lat_bin)) %>% 
  group_by(cells,date,vessel_class,flag) %>% 
  summarise(hours_in_gaps_under_12=sum(hours_in_gaps_under_12,na.rm=T),
            hours_in_gaps_over_12=sum(hours_in_gaps_over_12,na.rm=T)) %>% 
  ungroup() %>% 
  mutate(hours_in_gaps_under_12=na_if(hours_in_gaps_under_12,0))%>% 
  mutate(hours_in_gaps_over_12=na_if(hours_in_gaps_over_12,0))

gaps2=gaps %>% 
  ungroup() %>% 
  dplyr::select(c(cell_lon,cell_lat,gap_hours_per_cell,date,vessel_class,flag)) %>% 
  mutate(cells=raster::extract(template_1deg,.[,1:2]))  %>% 
  dplyr::select(-c(cell_lon,cell_lat)) %>% 
  group_by(cells,date,vessel_class,flag) %>% 
  summarise(total_gap=sum(gap_hours_per_cell,na.rm=T))%>% 
  ungroup()

gaps_disabling=gaps %>% 
  ungroup() %>% 
  filter(disabling_event=="TRUE") %>%
  dplyr::select(c(cell_lon,cell_lat,gap_hours_per_cell,date,vessel_class,flag)) %>% 
  mutate(cells=raster::extract(template_1deg,.[,1:2]))  %>% 
  dplyr::select(-c(cell_lon,cell_lat)) %>% 
  group_by(cells,date,vessel_class,flag) %>% 
  summarise(total_gap_disabling=sum(gap_hours_per_cell,na.rm=T))%>% 
  ungroup()

gaps_2w=gaps %>% 
  ungroup() %>% 
  filter(gaps_2w=="yes") %>% 
  dplyr::select(c(cell_lon,cell_lat,gap_hours_per_cell,date,vessel_class,flag)) %>% 
  mutate(cells=raster::extract(template_1deg,.[,1:2]))  %>% 
  dplyr::select(-c(cell_lon,cell_lat)) %>% 
  group_by(cells,date,vessel_class,flag) %>% 
  summarise(total_gap_2w=sum(gap_hours_per_cell,na.rm=T))%>% 
  ungroup()

gaps_2w_disabling=gaps %>% 
  ungroup() %>% 
  filter(gaps_2w=="yes") %>% 
  filter(disabling_event=="TRUE") %>% 
  dplyr::select(c(cell_lon,cell_lat,gap_hours_per_cell,date,vessel_class,flag)) %>% 
  mutate(cells=raster::extract(template_1deg,.[,1:2]))  %>% 
  dplyr::select(-c(cell_lon,cell_lat)) %>% 
  group_by(cells,date,vessel_class,flag) %>% 
  summarise(total_gap_disabling_2w=sum(gap_hours_per_cell,na.rm=T))%>% 
  ungroup()


master=full_join(fishing_clip_n,gaps2) %>% 
  full_join(.,gaps_disabling) %>% 
  full_join(.,gaps_2w) %>% 
  full_join(.,gaps_2w_disabling) %>% 
  mutate(cells=as.factor(cells)) %>% 
  gather(key, values,-c(cells,date,vessel_class,flag)) %>% 
  mutate(n_values=rescale(values,c(.00001,1))) %>% 
  dplyr::select(-values) %>% 
  spread(key,n_values)

master_test=master %>% 
  summarise(hours_in_gaps_under_12=sum(hours_in_gaps_under_12,na.rm=T),
            hours_in_gaps_over_12=sum(hours_in_gaps_over_12,na.rm=T),
            total_gap=sum(total_gap,na.rm=T),
            total_gap_disabling=sum(total_gap_disabling,na.rm=T),
            total_gap_2w=sum(total_gap_2w,na.rm=T),
            total_gap_disabling_2w=sum(total_gap_disabling_2w,na.rm=T)
  ) %>% 
  mutate(observed=hours_in_gaps_over_12+hours_in_gaps_under_12) %>%  ## observed ais
  mutate(total_activity=total_gap+hours_in_gaps_under_12) %>%  ## full observed + gaps
  mutate(total_2w_activity=total_gap_2w+hours_in_gaps_under_12) %>% # 2w observed + gaps
  mutate(perc_gaps=observed/total_activity) %>%  ## what percent of total activity is observed
  
  mutate(perc_unseen=total_gap_disabling/total_activity,
         perc_unseen_2w=total_gap_disabling_2w/total_2w_activity)

## now going to overlap ####

## only need to group by day and location, same species conditions regardless of who encounters it, commented out code demonstrates
# test=fishing_clip %>% 
#   filter(date=="2017-01-01") %>% 
#   filter(lat_bin== 57.875,
#          lon_bin==207.625
#          ) %>% 
#   dplyr::select(-c(X,...1,hours_in_gaps_over_12,hours,hours_in_gaps_under_12,fishing_hours,
#                    fishing_hours_in_gaps_over_12,fishing_hours_in_gaps_under_12,vessel_class,flag)) %>% 
#   distinct()
         

fishing_clip_species=fishing_clip %>% 
  dplyr::select(-c(X,...1,hours_in_gaps_over_12,hours,hours_in_gaps_under_12,fishing_hours,
                   fishing_hours_in_gaps_over_12,fishing_hours_in_gaps_under_12,vessel_class,flag)) %>% 
  gather(species,hab.suit,-c(date,lat_bin,lon_bin)) %>% 
  mutate(hab.suit=round(hab.suit,3)) %>% 
  .[complete.cases(.),] %>% ## get rid of pixels missing a given species
  distinct() %>%  ## get rid of duplicates, e.g. usa and can trawlers each encountered albacore on the same day in teh same pixel
  mutate(cells=raster::extract(template_1deg,.[,2:1])) %>% 
  group_by(cells,date,species) %>% 
  summarise(hab.suit_fishing=mean(hab.suit,na.rm=T)) %>% 
  ungroup()
            
gaps_clip_species=gaps %>% 
  ungroup() %>%
  dplyr::select(-c(X,...1,
                   cells,gap_start,gap_end,lon,lat,lon_end,lat_end,gap_id,
                   gap_hours,n_cells,gaps_2w,disabling_event,gap_hours_per_cell,vessel_class,flag)) %>% 
  gather(species,hab.suit,-c(date,cell_lat,cell_lon)) %>% 
  mutate(hab.suit=round(hab.suit,3)) %>% 
  .[complete.cases(.),] %>% ## get rid of pixels missing a given species
  distinct() %>%  ## get rid of duplicates, e.g. usa and can trawlers each encountered albacore on the same day in teh same pixel
  mutate(cells=raster::extract(template_1deg,.[,2:1])) %>% 
  group_by(cells,date,species) %>% 
  summarise(hab.suit_gaps=mean(hab.suit,na.rm=T))%>% 
  ungroup()

gaps_fish_species=
  full_join(gaps_clip_species,fishing_clip_species) %>% 
  mutate(hab.suit=case_when(is.na(hab.suit_fishing) & !is.na(hab.suit_gaps) ~ hab.suit_gaps,## if no fishing value, give it the gaps value
         !is.na(hab.suit_fishing) & is.na(hab.suit_gaps) ~ hab.suit_fishing, ## if no gaps value, give it the fishing value
         !is.na(hab.suit_fishing) & !is.na(hab.suit_gaps) ~ (hab.suit_fishing+hab.suit_gaps)/2)) %>% ## both have values, average
  mutate(hab.suit_n=scales::rescale(hab.suit,c(.00001,1))) %>% 
  dplyr::select(-c(hab.suit_gaps,hab.suit_fishing,hab.suit)) %>% 
  rename(hab.suit=hab.suit_n) %>% 
  mutate(cells=as.factor(cells))
 
## join AIS data + species data back together ####
## master_join:
# hours_in_gaps_over_12: normalized vessel activity hours in gaps over 12 hours
# hours_in_gaps_under_12: normalized vessel activity hours in gaps under 12 hours
# total_gap: normalized interpolated gap hours, all lengths
# total_gap_2w: normalized interpolated gap hours, under two weeks
# total_gap_disabling: normalized interpolated disabling hours, all lengths
# total_gap_disabling_2w: normalized interpolated disabling hours, under two weeks
# hab.suit: normalized species habitat suitability

master_join=full_join(master,gaps_fish_species) %>% 
  left_join(template_1deg_coords)

write.csv(master,glue("{outdir}/AIS_metrics_03-23-23.csv"))
write.csv(gaps_fish_species,glue("{outdir}/species_metrics_03-23-23.csv"))
write.csv(master_join,glue("{outdir}/cleaned_metrics_03-23-23.csv"))


