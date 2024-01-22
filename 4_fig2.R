### baselines (observed), upper and lower estimates for species habitat and fisheries hours and unseen overlap
library(patchwork)
library(nngeo)
library(sf)
library(patchwork)
source("/Users/HeatherWelch/Dropbox/OLE/github/OLE_Projects/utilities/load_libraries.R")

## load in data ####
datdir="UserDefined/cleaned_overlap_aggregated_dfs"
master_join=read.csv(glue("{datdir}/cleaned_metrics_03-23-23.csv"))
glimpse(master_join)
masterAIS=read.csv(glue("{datdir}/AIS_metrics_03-23-23.csv"))

outdir="UserDefined";dir.create(outdir)

## species baselines
total_habitat=list.files("UserDefined/stats_08_03_23b",full.names = T) %>% 
  lapply(read_csv) %>% 
  bind_rows %>% 
  group_by(species) %>% 
  summarise(total.habitat=sum(habitat_sum))

## vessel baselines
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

## calculating oberved vs observed+unseen species habitat covered by vessels ####
observed_actvity_species=master_join %>% 
  mutate(species=case_when(species=="black.footedAlbatross_Dallas"~"black-footedAlbatross_Dallas",
                           T~species)) %>% 
  filter(!is.na(hours_in_gaps_over_12) | !is.na(hours_in_gaps_under_12) ) %>% ## pixels with observed activity only
  dplyr::select(c(species,hab.suit,date,cell_lon,cell_lat)) %>%
  distinct() %>% 
  group_by(species) %>% 
  summarise(observed_overlap=sum(hab.suit,na.rm=T))%>% 
  left_join(.,total_habitat) %>% 
  mutate(observed_hab_suit_perc=observed_overlap/total.habitat)%>% 
  dplyr::select(species,observed_hab_suit_perc,observed_overlap) %>% 
  .[complete.cases(.),]

total_activity_species=master_join %>% 
  mutate(species=case_when(species=="black.footedAlbatross_Dallas"~"black-footedAlbatross_Dallas",
                           T~species)) %>% 
  filter(!is.na(total_gap) | !is.na(hours_in_gaps_under_12)) %>% 
  dplyr::select(c(species,hab.suit,date,cell_lon,cell_lat)) %>%
  distinct() %>% 
  group_by(species) %>% 
  summarise(total_activity_overlap=sum(hab.suit,na.rm=T))%>% 
  left_join(.,total_habitat) %>% 
  mutate(total_activity_hab_suit_perc=total_activity_overlap/total.habitat) %>% 
  dplyr::select(species,total_activity_hab_suit_perc,total_activity_overlap)%>% 
  .[complete.cases(.),]

total_activity_2w_species=master_join %>% 
  mutate(species=case_when(species=="black.footedAlbatross_Dallas"~"black-footedAlbatross_Dallas",
                           T~species)) %>% 
  filter(!is.na(total_gap_2w) | !is.na(hours_in_gaps_under_12)) %>% 
  dplyr::select(c(species,hab.suit,date,cell_lon,cell_lat)) %>%
  distinct() %>% 
  group_by(species) %>% 
  summarise(total_activity_2w_overlap=sum(hab.suit,na.rm=T))%>% 
  left_join(.,total_habitat) %>% 
  mutate(total_activity_2w_hab_suit_perc=total_activity_2w_overlap/total.habitat) %>% 
  dplyr::select(species,total_activity_2w_hab_suit_perc,total_activity_2w_overlap)%>% 
  .[complete.cases(.),]

master_species=left_join(total_activity_2w_species,total_activity_species) %>% 
  left_join(.,observed_actvity_species) %>% 
  dplyr::select(c(species,total_activity_2w_hab_suit_perc,total_activity_hab_suit_perc,observed_hab_suit_perc)) %>% 
  mutate(total_activity_hab_diff=total_activity_hab_suit_perc-observed_hab_suit_perc,
         total_activity_2w_diff=total_activity_2w_hab_suit_perc-observed_hab_suit_perc) %>% 
  dplyr::select(c(species,observed_hab_suit_perc,total_activity_2w_diff,total_activity_hab_diff)) %>% 
  gather(metric,value,-c(species)) %>% 
  mutate(metric=factor(metric,levels=rev(c("observed_hab_suit_perc","total_activity_2w_diff","total_activity_hab_diff")))) %>% 
  mutate(nice_names=gsub("_TOPP","",species)) %>% 
  mutate(nice_names=gsub("_Dallas","",nice_names)) %>%
  mutate(nice_names=replace(nice_names,nice_names=="sootyShearwater","Sooty shearwater")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="pacificBluefinTuna","Bluefin tuna")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="makoShark","Mako shark")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="leatherbackTurtle","Leatherback turtle")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="elephantSeal","Elephant seal")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="black-footedAlbatross","Black-footed albatross")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="blueShark","Blue shark")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="blueWhale","Blue whale")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="californiaSeaLion","California sea lion")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="albacoretuna","Albacore tuna")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="salmonShark","Salmon shark")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="whiteShark","White shark")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="laysanAlbatross","Laysan albatross")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="yellowfinTuna","Yellowfin tuna"))

order=master_species %>% 
  filter(metric=="observed_hab_suit_perc") %>% 
  arrange(desc(value)) %>% 
  pull(nice_names)

master_species=master_species %>% 
  mutate(nice_names=factor(nice_names,levels=order))

## stats for abstract, how much of species habitat contains fishing activity
all_habitat=total_habitat %>% pull(total.habitat) %>% sum()
observed_actvity_species_total=master_join %>% 
  mutate(species=case_when(species=="black.footedAlbatross_Dallas"~"black-footedAlbatross_Dallas",
                           T~species)) %>% 
  filter(!is.na(hours_in_gaps_over_12) | !is.na(hours_in_gaps_under_12) ) %>% ## pixels with observed activity only
  dplyr::select(c(species,hab.suit,date,cell_lon,cell_lat)) %>%
  distinct() %>% 
  summarise(observed_overlap=sum(hab.suit,na.rm=T))%>% 
  mutate(total.habitat=all_habitat) %>% 
  mutate(observed_hab_suit_perc=observed_overlap/total.habitat)%>% 
  dplyr::select(observed_hab_suit_perc,observed_overlap) %>% 
  .[complete.cases(.),]

## calculating oberved vs observed+unseen vessel activity covered by species habitat ####
total_vessel=master_join %>% 
  dplyr::select(-c(cells,date,vessel_class,flag,cell_lon,cell_lat,X)) %>% 
  group_by(species) %>% 
  summarise_all(sum,na.rm=T) %>% 
  ungroup() ## this just simplifies the dataset

vessels=total_vessel %>% 
  mutate(observed=hours_in_gaps_over_12+hours_in_gaps_under_12) %>% 
  mutate(total_activity=total_gap+hours_in_gaps_under_12) %>%  ## full observed + gaps
  mutate(total_2w_activity=total_gap_2w+hours_in_gaps_under_12) %>%  # 2w observed + gaps
  dplyr::select(c(species,observed,total_activity,total_2w_activity,total_gap,total_gap_2w)) %>% 
  mutate(species=case_when(species=="black.footedAlbatross_Dallas"~"black-footedAlbatross_Dallas",
                           T~species)) %>% 
  mutate(observed_vessels_perc=observed/20835.39)%>% ## from FA_totals, total observed vessel activity
  mutate(total_activity_vessels_perc=total_activity/25771.4) %>% ## from FA_totals, total observed+unseen vessel activity
  mutate(total_activity_2w_vessels_perc=total_2w_activity/22173.18) %>% ## from FA_totals, total observed+unseen_2w vessel activity
  dplyr::select(species,observed_vessels_perc,total_activity_vessels_perc,total_activity_2w_vessels_perc)%>% 
  .[complete.cases(.),]

vessels_common_denom=total_vessel %>% 
  mutate(observed=hours_in_gaps_over_12+hours_in_gaps_under_12) %>% 
  mutate(total_activity=total_gap+hours_in_gaps_under_12) %>%  ## full observed + gaps
  mutate(total_2w_activity=total_gap_2w+hours_in_gaps_under_12) %>%  # 2w observed + gaps
  dplyr::select(c(species,observed,total_activity,total_2w_activity,total_gap,total_gap_2w)) %>% 
  mutate(species=case_when(species=="black.footedAlbatross_Dallas"~"black-footedAlbatross_Dallas",
                           T~species)) %>% 
  mutate(observed_vessels_perc=observed/20835.39)%>% ## from FA_totals, total observed vessel activity
  mutate(total_activity_vessels_perc=total_activity/20835.39) %>% ## from FA_totals, total observed+unseen vessel activity
  mutate(total_activity_2w_vessels_perc=total_2w_activity/20835.39) %>% ## from FA_totals, total observed+unseen_2w vessel activity
  dplyr::select(species,observed_vessels_perc,total_activity_vessels_perc,total_activity_2w_vessels_perc)%>% 
  .[complete.cases(.),]

stats_for_fig=vessels_common_denom %>% 
  mutate(perc=percent(observed_vessels_perc,1)) %>% 
  mutate(nice_names=gsub("_TOPP","",species)) %>% 
  mutate(nice_names=gsub("_Dallas","",nice_names)) %>%
  mutate(nice_names=replace(nice_names,nice_names=="sootyShearwater","Sooty shearwater")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="pacificBluefinTuna","Bluefin tuna")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="makoShark","Mako shark")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="leatherbackTurtle","Leatherback turtle")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="elephantSeal","Elephant seal")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="black-footedAlbatross","Black-footed albatross")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="blueShark","Blue shark")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="blueWhale","Blue whale")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="californiaSeaLion","California sea lion")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="albacoretuna","Albacore tuna")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="salmonShark","Salmon shark")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="whiteShark","White shark")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="laysanAlbatross","Laysan albatross")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="yellowfinTuna","Yellowfin tuna")) %>% 
  mutate(nice_names=factor(nice_names,levels=order)) %>% 
  dplyr::select(c(nice_names,perc)) %>% 
  arrange(nice_names)

vessels_hours=total_vessel %>% 
  mutate(observed=hours_in_gaps_over_12+hours_in_gaps_under_12) %>% 
  mutate(total_activity=total_gap+hours_in_gaps_under_12) %>%  ## full observed + gaps
  mutate(total_2w_activity=total_gap_2w+hours_in_gaps_under_12) %>%  # 2w observed + gaps
  dplyr::select(c(species,observed,total_activity,total_2w_activity,total_gap,total_gap_2w)) %>% 
  mutate(species=case_when(species=="black.footedAlbatross_Dallas"~"black-footedAlbatross_Dallas",
                           T~species)) %>% 
  .[complete.cases(.),] %>% 
  dplyr::select(-c(total_2w_activity,total_activity)) %>% 
  gather(metric,value,-c(species)) %>% 
  mutate(metric=factor(metric,levels=rev(c("observed","total_gap_2w","total_gap")))) %>% 
  mutate(nice_names=gsub("_TOPP","",species)) %>% 
  mutate(nice_names=gsub("_Dallas","",nice_names)) %>%
  mutate(nice_names=replace(nice_names,nice_names=="sootyShearwater","Sooty shearwater")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="pacificBluefinTuna","Bluefin tuna")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="makoShark","Mako shark")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="leatherbackTurtle","Leatherback turtle")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="elephantSeal","Elephant seal")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="black-footedAlbatross","Black-footed albatross")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="blueShark","Blue shark")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="blueWhale","Blue whale")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="californiaSeaLion","California sea lion")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="albacoretuna","Albacore tuna")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="salmonShark","Salmon shark")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="whiteShark","White shark")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="laysanAlbatross","Laysan albatross")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="yellowfinTuna","Yellowfin tuna")) %>% 
  mutate(nice_names=factor(nice_names,levels=order))


## plots  ####
vessels_plot=ggplot(vessels_hours)+
  geom_bar(aes(x=value,y=nice_names,fill=metric),stat="identity",color="black")+
  theme_classic()+
  scale_fill_manual("",values=c("observed"="#0073C2FF",
                                   "total_gap_2w"="#EFC000FF",
                                   "total_gap"="#A73030FF"),
                    labels=c("observed"="Observed",
                              "total_gap_2w"="Lower estimate",
                              "total_gap"="Upper estimate"))+
  scale_x_continuous(expand=c(0,0),
                     breaks=c(0,5000,10000,15000,20000),
                     labels=c("0","5k","10k","15k","20k"))+
  xlab("Hours")+
  ylab(NULL)+
  theme(legend.position = "none")

species=ggplot(master_species)+
  geom_bar(aes(x=value,y=nice_names,fill=metric),stat="identity",color="black")+
  theme_classic()+
  scale_fill_manual("",values=c("observed_hab_suit_perc"="#0073C2FF",
                                "total_activity_2w_diff"="#EFC000FF",
                                "total_activity_hab_diff"="#A73030FF"),
                    labels=c("observed_hab_suit_perc"="Observed",
                             "total_activity_2w_diff"="Lower estimate",
                             "total_activity_hab_diff"="Upper estimate"))+
  scale_x_continuous(labels = percent,expand=c(0,0))+
  # xlab("% of total species habitat")+
  xlab("Percent")+
  ylab(NULL)+
  theme(legend.position = "none")

legend=ggplot(vessels_hours)+
  geom_bar(aes(x=value,y=nice_names,fill=metric),stat="identity",color="black")+
  theme_classic()+
  scale_fill_manual("",values=c("observed"="#0073C2FF",
                                "total_gap_2w"="#EFC000FF",
                                "total_gap"="#A73030FF"),
                    labels=c("observed"="Observed",
                             "total_gap_2w"="Lower estimate",
                             "total_gap"="Upper estimate"))+
  scale_x_continuous(expand=c(0,0),
                     breaks=c(0,5000,10000,15000,20000),
                     labels=c("0","5k","10k","15k","20k"))+
  xlab("Hours")+
  ylab(NULL)

## % increase overlap ####
fig2_species=master_join %>% 
  filter(!is.na(species)) %>%
  
  mutate(overlap_hours_in_gaps_over_12=hours_in_gaps_over_12*hab.suit,
         overlap_hours_in_gaps_under_12=hours_in_gaps_under_12*hab.suit,
         overlap_total_gap_2w=total_gap_2w*hab.suit,
         overlap_total_gap=total_gap*hab.suit,
         overlap_total_gap_disabling_2w=total_gap_disabling_2w*hab.suit,
         overlap_total_gap_disabling=total_gap_disabling*hab.suit) %>%
  dplyr::select(-c(cells,date,vessel_class,flag,hab.suit)) %>% 
  group_by(species) %>% 
  summarise_all(sum,na.rm=T) %>% 
  ungroup() %>% 
  
  ## overlap metrics
  mutate(overlap_observed=overlap_hours_in_gaps_over_12+overlap_hours_in_gaps_under_12) %>%  ## observed ais
  mutate(overlap_total_2w_activity=overlap_total_gap_2w+overlap_hours_in_gaps_under_12) %>% # 2w observed + gaps
  mutate(overlap_total_activity=overlap_total_gap+overlap_hours_in_gaps_under_12) %>%  # observed + gaps

  mutate(overlap_perc_2w_gaps=overlap_observed/overlap_total_2w_activity) %>%  ## what percent of total activity is observed 2w
  mutate(overlap_perc_increase_2w_gaps=(overlap_total_2w_activity-overlap_observed)/overlap_total_2w_activity) %>%  ## what percent of total activity is observed 2w
  mutate(overlap_perc_unseen_2w=overlap_total_gap_disabling_2w/overlap_total_2w_activity) %>%   ## full what percent of total activity is disabling 2w
  
  mutate(overlap_perc_gaps=overlap_observed/overlap_total_activity) %>%  ## what percent of total activity is observed 
  mutate(overlap_perc_increase_gaps=(overlap_total_activity-overlap_observed)/overlap_total_activity) %>%  ## what percent of total activity is observed 
  mutate(overlap_perc_unseen=overlap_total_gap_disabling/overlap_total_activity) %>%   ## full what percent of total activity is disabling
  
  mutate(nice_names=gsub("_TOPP","",species)) %>% 
  mutate(nice_names=gsub("_Dallas","",nice_names)) %>%
  mutate(nice_names=replace(nice_names,nice_names=="sootyShearwater","Sooty shearwater")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="pacificBluefinTuna","Bluefin tuna")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="makoShark","Mako shark")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="leatherbackTurtle","Leatherback turtle")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="elephantSeal","Elephant seal")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="black.footedAlbatross","Black-footed albatross")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="blueShark","Blue shark")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="blueWhale","Blue whale")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="californiaSeaLion","California sea lion")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="albacoretuna","Albacore tuna")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="salmonShark","Salmon shark")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="whiteShark","White shark")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="laysanAlbatross","Laysan albatross")) %>%
  mutate(nice_names=replace(nice_names,nice_names=="yellowfinTuna","Yellowfin tuna")) %>% 
  mutate(nice_names=factor(nice_names,levels=c(order))) 

fig_2_overlap=ggplot(fig2_species)+
  geom_segment(aes(y=nice_names,yend=nice_names,x=overlap_perc_increase_2w_gaps,xend=overlap_perc_increase_gaps))+
  geom_point(aes(y=nice_names,x=overlap_perc_increase_2w_gaps),size=3,color="#EFC000FF")+
  geom_point(aes(y=nice_names,x=overlap_perc_increase_gaps),size=3,color="#A73030FF")+
  theme_classic()+
  ylab(NULL)+
  xlab("Percent")+
  scale_x_continuous(labels=percent)+
  theme(legend.position = "none")

png(glue("{outdir}/overlap_fig2.png"),width=8,height=13.6,units='cm',res=400,type = "cairo")
par(ps=10)
par(mar=c(4,4,1,1))
par(cex=1)
print({fig_2_overlap})
# gg_hm
dev.off()

png(glue("{outdir}/vessels_fig2.png"),width=8,height=13.6,units='cm',res=400,type = "cairo")
par(ps=10)
par(mar=c(4,4,1,1))
par(cex=1)
print({vessels_plot})
# gg_hm
dev.off()

png(glue("{outdir}/habitat_fig2.png"),width=8,height=13.6,units='cm',res=400,type = "cairo")
par(ps=10)
par(mar=c(4,4,1,1))
par(cex=1)
print({species})
# gg_hm
dev.off()

png(glue("{outdir}/legend_fig2.png"),width=8,height=13.6,units='cm',res=400,type = "cairo")
par(ps=10)
par(mar=c(4,4,1,1))
par(cex=1)
print({legend})
# gg_hm
dev.off()

