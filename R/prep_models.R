add_rarity = function(){ 

melted_rares = melted_cov %>%
  filter(Presence > 0) %>% 
  group_by(Species) %>% 
  summarise(n_occu = n()) %>% 
  mutate(decile_rank = ntile(n_occu, 10)) %>% 
  # mutate(Rarity = ifelse(decile_rank > 5, "Common", "Rare"))
  mutate(Rarity = ifelse(n_occu > 100, "Common", "Rare"))

ggplot(melted_rares, aes(x = n_occu, colour = Rarity)) + 
  geom_histogram(binwidth = 5) + 
  xlim(0, 500)

return(melted_rares)

}

prep_data_occu = function() {
  
  traits$Species = traits$CURRENT_SPECIES_NAME
  
  df_occu_traits = OS_occu %>% 
    mutate(log_full = log(RR_full),
           log_high = log(RR_high),
           log_light = log(RR_light)) %>% 
    left_join(traits, by = "Species") %>% 
    left_join(melted_rares %>% dplyr::select(Species, Rarity, decile_rank), by = "Species")  %>% 
    mutate(log_length = log(MaxLength)) 
  
  df_occu_full = df_occu_traits %>% 
    drop_na(log_full, Trophic.Level, log_length)
  df_occu_high = df_occu_traits %>% 
    drop_na(log_high, Trophic.Level, log_length)
  df_occu_light = df_occu_traits %>% 
    drop_na(log_light, Trophic.Level, log_length)
  
  list_occu = list(df_occu_full, df_occu_high, df_occu_light)
}

prep_data_abun = function(){
  
  traits$Species = traits$CURRENT_SPECIES_NAME
  
  df_ab_traits = df_ab %>%
    mutate(log_full = log(IRR_full),
           log_high = log(IRR_high),
           log_light = log(IRR_light)) %>% 
    left_join(traits, by = "Species") %>% 
    left_join(melted_rares %>% dplyr::select(Species, Rarity, decile_rank), by = "Species")  %>% 
    mutate(log_length = log(MaxLength))
  
  df_ab_full = df_ab_traits %>% 
    drop_na(log_full, Trophic.Level, log_length)
  df_ab_high = df_ab_traits %>% 
    drop_na(log_high, Trophic.Level, log_length)
  df_ab_light = df_ab_traits %>% 
    drop_na(log_light, Trophic.Level, log_length) 
  
  list_abun = list(df_ab_full, df_ab_high, df_ab_light)
  
}

prep_data_bm = function() {
  
  traits$Species = traits$CURRENT_SPECIES_NAME
  
  df_bm_traits = df_bm_product %>%
    mutate(log_full = log(IRR_product_full),
           log_high = log(IRR_product_high),
           log_light = log(IRR_product_light)) %>% 
    left_join(traits, by = "Species") %>% 
    left_join(melted_rares %>% dplyr::select(Species, Rarity, decile_rank), by = "Species")  %>% 
    mutate(log_length = log(MaxLength))
  
  df_bm_full = df_bm_traits %>% 
    drop_na(log_full, Trophic.Level, log_length)  
  df_bm_high = df_bm_traits %>% 
    drop_na(log_high, Trophic.Level, log_length)
  df_bm_light = df_bm_traits %>% 
    drop_na(log_light, Trophic.Level, log_length)
  
  list_bm = list(df_bm_full, df_bm_high, df_bm_light)
  
}

