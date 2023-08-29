tidy_sites = function() {
  
  sites_correct$Effectiveness_cor = "out"  #Set to outside all sites that were unprotected when surveyed
  
  covariates_corrected = dplyr::left_join(covariates, 
                                          sites_correct[,c(1,21)], 
                                          by = "SurveyID")
  
  covariates_corrected$Effectiveness = replace(covariates_corrected$Effectiveness, 
                                               covariates_corrected$Effectiveness_cor == "out",
                                               "out") # Correct it on the original covariates data
  
  covar_jb = read.csv(here::here("data", "raw_data", "Site_Env_Socio_MPA_metadata_AllenAtlas_JB.csv"))
  covar_graham = read.csv(here::here("data", "raw_data", "NEOLI_Graham2014.csv"))
  covar_jb_graham = covar_jb |> 
    left_join(covar_graham |> dplyr::select(SiteCode, area..km2.)) |>
    rename(Area_km2 = "area..km2.")
  
  covar_jb_graham$Area_km2 = as.numeric(covar_jb_graham$Area_km2)
  
  covar_jb_select = covar_jb_graham %>% select(SurveyID, DiffYear, Area_km2)
  covar_jb_sites = covar_jb_select[covar_jb_select$SurveyID %in% covariates_corrected$SurveyID,]
  
  covariates_corrected = covariates_corrected %>% 
    left_join(covar_jb_sites, by = "SurveyID")
  covariates_corrected = covariates_corrected[-which(duplicated(covariates_corrected$SurveyID)),]
  
  return(covariates_corrected)
}

make_ssmatrix_melted = function(n) {
  
  tab <- table(fish_data$valid_name_FishBase)
  boxplot(tab, outline = F)
  tab <- tab[tab>n]
  fish_filt  <- fish_data[fish_data$valid_name_FishBase %in% names(tab),]
  fish_filt$Species = as.factor(fish_filt$valid_name_FishBase)
  
  ssmatrix = t(fossil::create.matrix(fish_filt, 
                                     tax.name = "Species", 
                                     locality = "SurveyID",
                                     abund = F))
  
  melted = ssmatrix |>
    reshape2::melt() |>
    rename(SurveyID = "Var1", 
           Species = "Var2", 
           Presence = "value") |>
    left_join(covariates_corrected |> 
                dplyr::select(SurveyID, 
                              Effectiveness, 
                              DiffYear,
                              Area_km2), by = "SurveyID") |>
    left_join(sites_info[,c(1:4,10)], by = "SurveyID") |>
    drop_na(Effectiveness)

  melted$Effectiveness = fct_collapse(melted$Effectiveness, 
                                      "Highly Protected" = c("High No take", 
                                                             "High No take multizoned",
                                                             "Medium No take", 
                                                             "Medium No take multizoned"),
                                      "Lightly Protected" = c("High Restricted take", 
                                                              "High Restricted take multizoned", 
                                                              "Medium Restricted take",
                                                              "Medium Restricted take multizoned",
                                                              "Low No take", 
                                                              "Low No take multizoned",
                                                              "Low Restricted take",
                                                              "Low Restricted take multizoned"),
                                      "Unprotected" = c("Medium Fishing", "Low Fishing", "out"))
  
  melted$Effectiveness = as.character(melted$Effectiveness)
  melted$Effectiveness[melted$Effectiveness == "Highly Protected" & melted$DiffYear >= 10 & melted$Area_km2 >= 100] = "Fully Protected"
  melted$Effectiveness = as.factor(melted$Effectiveness)
  melted$Effectiveness[melted$Effectiveness != "Unprotected" & melted$DiffYear < 1] = "Unprotected"
  
  return(melted)
  
}

rid_rare_outside = function() { # Get rid of species with 0-4 occurrences in the reference level
  
  # Some values were dropped, make sure we still have at least N occurrences
  
  melted_sum = melted |>
    group_by(Species) |>
    summarise(Presence = sum(Presence)) |>
    dplyr::filter(Presence >= 30)
  
  spec_vec = as.character(melted_sum$Species)
  
  melted_filt = melted[melted$Species %in% spec_vec,]
  
  # Get rid of rares outside
  
  melted_fsum = melted_filt |>
    filter(Effectiveness == "Unprotected") |>
    group_by(Species) |>
    summarise(Presence = sum(Presence)) |>
    filter(Presence > 4) 
  
  spec_vec = as.character(melted_fsum$Species)
  
  melted_filt = melted_filt[melted_filt$Species %in% spec_vec,]
  melted_filt$Species = as.factor(as.character(melted_filt$Species))

}

add_covariates = function(){
  
  all_cov = hab_filt |>
    left_join(socio_filt, by = "SurveyID") |>
    left_join(env_filt, by = "SurveyID") |>
    left_join(fine_habitat[,-c(2:4)], by = "SurveyID")
  
  melted_cov = melted_filt |> 
    left_join(all_cov, by = "SurveyID")
  
  return(melted_cov)
  
}

compute_ES_occu = function() {
  
  melted_cov$Effectiveness = relevel(as.factor(melted_cov$Effectiveness), ref = "Unprotected")
  
  all_sp_glm_qual = lapply(1:length(unique(melted_cov$Species)), function(i) {
  # all_sp_glm_qual = lapply(101:150, function(i) {
    
    subdf = filter(melted_cov, Species == unique(melted_cov$Species)[i])
    
    pres_sp1 = subdf %>%
      filter(Presence > 0)
    
    pres_sp1$Ecoregion = as.factor(as.character(pres_sp1$Ecoregion))
    
    eco_table = table(pres_sp1$Ecoregion)
    eco_tab = names(eco_table[eco_table > 4])
    
    abs_sp1 = subdf %>%
      filter(Presence < 1)
    
    abs_sp1 = abs_sp1[abs_sp1$Ecoregion %in% eco_tab,]
    pres_sp1 = pres_sp1[pres_sp1$Ecoregion %in% eco_tab,]
    
    presabs_sp1 = rbind(pres_sp1, abs_sp1)
    
    eco_weights = presabs_sp1 |>
      group_by(Ecoregion) |> 
      mutate(weights = (1/n())) |> 
      ungroup()
    
    pca_prep = all_cov[all_cov$SurveyID %in% unique(presabs_sp1$SurveyID),]
    pca_prep = pca_prep |> 
      left_join(eco_weights %>% dplyr::select(SurveyID, weights), by = "SurveyID")
    
    # ACP = FactoMineR::PCA(pca_prep[, -c(1, 97)],
    #                       ncp = 25,
    #                       scale.unit = T,
    #                       row.w = pca_prep$weights,
    #                       graph = F)

    pca = ade4::dudi.pca(pca_prep[, -c(1, 97)],
                         center = T, scale = T, 
                         row.w = pca_prep$weights,
                         nf = 25, 
                         scannf = F)
    
    # cum_var = which(factoextra::get_eig(ACP)[,3] > 70)[1]
    # PC_coords1 = as.data.frame(ACP$ind$coord)
    
    PC_coords = pca$li
    
    eigs = factoextra::get_eig(pca)
    cum_var = which(eigs[,3] > 70)[1]
    
    pca_prep_and_coords = cbind(pca_prep, PC_coords)
    presabs_sp1 = presabs_sp1[, 1:10] |> 
      left_join(pca_prep_and_coords[, c(1, 98:122)], by = "SurveyID") 
                        
    istherefull = ifelse(sum(presabs_sp1$Presence[which(presabs_sp1$Effectiveness == "Fully Protected")]) > 4, T, F)
    istherehigh = ifelse(sum(presabs_sp1$Presence[which(presabs_sp1$Effectiveness == "Highly Protected")]) > 4, T, F)
    istherelight = ifelse(sum(presabs_sp1$Presence[which(presabs_sp1$Effectiveness == "Lightly Protected")]) > 4, T, F)
    isthereunp = ifelse(sum(presabs_sp1$Presence[which(presabs_sp1$Effectiveness == "Unprotected")]) > 4, T, F)
    
    presabs_sp1_nofull = presabs_sp1 %>%
      filter(Effectiveness != "Fully Protected")
    presabs_sp1_nohigh = presabs_sp1 %>%
      filter(Effectiveness != "Highly Protected")
    presabs_sp1_nolight = presabs_sp1 %>% 
      filter(Effectiveness != "Lightly Protected")
    presabs_sp1_nofullhigh = presabs_sp1 %>% 
      filter(Effectiveness != "Highly Protected" & Effectiveness != "Fully Protected")
    presabs_sp1_nofulllight = presabs_sp1 %>% 
      filter(Effectiveness != "Lightly Protected" & Effectiveness != "Fully Protected")
    presabs_sp1_nohighlight = presabs_sp1 %>% 
      filter(Effectiveness != "Lightly Protected" & Effectiveness != "Highly Protected")
    
    if(istherefull == F){presabs_sp1 = presabs_sp1_nofull}
    if(istherehigh == F){presabs_sp1 = presabs_sp1_nohigh}
    if(istherelight == F){presabs_sp1 = presabs_sp1_nolight}
    if(istherelight == F & istherefull == F){presabs_sp1 = presabs_sp1_nofulllight}
    if(istherehigh == F & istherefull == F){presabs_sp1 = presabs_sp1_nofullhigh}
    if(istherelight == F & istherehigh == F){presabs_sp1 = presabs_sp1_nohighlight}
    
    if(isthereunp == F){return(NULL)} 
    if(istherefull == F & istherehigh == F & istherelight == F){return(NULL)}
    
    names_var = c("Effectiveness", names(as.data.frame(PC_coords[, c(1:cum_var)])))
    formula_var = paste0(names_var, collapse = "+")
    full_formula = paste("Presence ~ ", formula_var, sep = "")
    
    model = glm(formula = full_formula,
                family = "binomial",
                data = presabs_sp1)
    
    if(model$converged == F){return(NULL)}
    
    sum_isthere = sum(istherefull, istherehigh, istherelight)
    if(sum_isthere < 2){vif = sqrt((regclass::VIF(model)))} else{vif = regclass::VIF(model)[, 3]}
    
    while(max(vif) > sqrt(5)){ 
      
      u = which.max(vif)
      names_var = names_var[-u]
      formula_var = paste0(names_var, collapse = "+")
      full_formula = paste("Presence ~ ", formula_var, sep = "")
      
      model = glm(formula = full_formula,
                  family = "binomial",
                  data = presabs_sp1)            
      
      if(sum(str_detect(names_var, "Effectiveness")) == 0 | length(names_var) < 2){return(NULL)}
      if(sum_isthere < 2){vif = sqrt((regclass::VIF(model)))} else{vif = regclass::VIF(model)[, 3]}
      
      full_formula = full_formula
      
    }
    
    # full_formula_sp = as.formula(paste(full_formula, "+ Matern(1 | SiteLongitude + SiteLatitude)", sep = ""))
    # separ = spaMM::is_separated.formula(formula = full_formula_sp,
                                        # data = presabs_sp1)
    
    # if(separ == T){return(NULL)} 
    # 
    # sp_model = spaMM::fitme(Presence ~ Effectiveness + Matern(1 | SiteLongitude + SiteLatitude),
    #                         family = "binomial",
    #                         method = "PQL/L",
    #                         data = presabs_sp1)
    
    sp_model = spaMM::fitme(full_formula_sp,
                            family = "binomial",
                            method = "PQL/L",
                            data = presabs_sp1)
    
    # model = glm(full_formula, 
    #             family = "binomial", 
    #             data = presabs_sp1)
    
    # mod_coefs = as.data.frame(sp_model$fixef)
    mod_coefs = as.data.frame(model$coefficients)
    
    intercept = mod_coefs["(Intercept)",]
    
    estimate_outside = ifelse(sum(stringr::str_detect(row.names(mod_coefs),"Intercept")) == 1,
                              boot::inv.logit(intercept),
                              NA)
    
    estimate_full = ifelse(sum(stringr::str_detect(row.names(mod_coefs),"EffectivenessFully Protected")) == 1,
                           boot::inv.logit(intercept+mod_coefs["EffectivenessFully Protected",]),
                           NA)
    
    estimate_high = ifelse(sum(stringr::str_detect(row.names(mod_coefs),"EffectivenessHighly Protected")) == 1,
                           boot::inv.logit(intercept+mod_coefs["EffectivenessHighly Protected",]),
                           NA)
    
    estimate_light = ifelse(sum(stringr::str_detect(row.names(mod_coefs),"EffectivenessLightly Protected")) == 1,
                            boot::inv.logit(intercept+mod_coefs["EffectivenessLightly Protected",]),
                            NA)
    
    # OR_full = ifelse(sum(stringr::str_detect(row.names(mod_coefs),"EffectivenessFully Protected")) == 1,
    #                  exp(mod_coefs["EffectivenessFully Protected",]),
    #                  NA)
    # 
    # OR_high = ifelse(sum(stringr::str_detect(row.names(mod_coefs),"EffectivenessHighly Protected")) == 1,
    #                  exp(mod_coefs["EffectivenessHighly Protected",]),
    #                  NA)
    # 
    # OR_light = ifelse(sum(stringr::str_detect(row.names(mod_coefs),"EffectivenessLightly Protected")) == 1,
    #                   exp(mod_coefs["EffectivenessLightly Protected",]),
    #                   NA)
    
    RR_full = estimate_full/estimate_outside
    RR_high = estimate_high/estimate_outside
    RR_light = estimate_light/estimate_outside
    
    p = spaMM::predict.HLfit(sp_model, type = "response")
    # p = predict(model, type = "response")
    roc.mod = pROC::roc(presabs_sp1$Presence, as.numeric(p))
    AUC = pROC::auc(roc.mod)
    
    Species = unique(melted_cov$Species)[i]
    
    coef = data.frame(Species,
                      estimate_outside,
                      estimate_full, 
                      estimate_high,
                      estimate_light,
                      RR_full,
                      RR_high,
                      RR_light,
                      AUC)
    
    cat(i, "\n")
    return(coef)
    
  }) 
  
  OS_occu = data.table::rbindlist(all_sp_glm_qual)
  save(OS_occu, file = here::here("outputs", "OS_occu.RData"))
  return(OS_occu)
  
}

make_ssmatrix_melted_bm = function(n) {
  
  tab <- table(fish_data$valid_name_FishBase)
  boxplot(tab, outline = F)
  tab <- tab[tab>n]
  fish_filt  <- fish_data[fish_data$valid_name_FishBase %in% names(tab),]
  fish_filt$Species = as.factor(fish_filt$valid_name_FishBase)
  
  fish_filt = fish_filt |>
    mutate(trunc_bm = trunc(Biomass))
  
  fish_filt1 = fish_filt |>
    group_by(Species) |>
    mutate(percentile_99 = ifelse(max(Num) > 499, quantile(Num, 0.99), quantile(Num, 1))) |>
    filter(Num < percentile_99) |>
    ungroup()
  
  ssmatrix_bm = t(fossil::create.matrix(as.data.frame(fish_filt1), 
                                          tax.name = "Species", 
                                          locality = "SurveyID",
                                          abund = T,
                                          abund.col = "trunc_bm"))
  
  melted_bm = ssmatrix_bm |>
    reshape2::melt() |>
    rename(SurveyID = "Var1", 
           Species = "Var2", 
           trunc_bm = "value") |>
    left_join(covariates_corrected |> 
                dplyr::select(SurveyID, 
                              Effectiveness, 
                              DiffYear,
                              Area_km2), by = "SurveyID") |>
    left_join(sites_info[,c(1:4,10)], by = "SurveyID") |>
    drop_na(Effectiveness)
  
  melted_bm$Effectiveness = fct_collapse(melted_bm$Effectiveness, 
                                      "Highly Protected" = c("High No take", 
                                                             "High No take multizoned",
                                                             "Medium No take", 
                                                             "Medium No take multizoned"),
                                      "Lightly Protected" = c("High Restricted take", 
                                                              "High Restricted take multizoned", 
                                                              "Medium Restricted take",
                                                              "Medium Restricted take multizoned",
                                                              "Low No take", 
                                                              "Low No take multizoned",
                                                              "Low Restricted take",
                                                              "Low Restricted take multizoned"),
                                      "Unprotected" = c("Medium Fishing", "Low Fishing", "out"))
  
  melted_bm$Effectiveness = as.character(melted_bm$Effectiveness)
  melted_bm$Effectiveness[melted_bm$Effectiveness == "Highly Protected" & melted_bm$DiffYear >= 10 & melted_bm$Area_km2 >= 100] = "Fully Protected"
  melted_bm$Effectiveness = as.factor(melted_bm$Effectiveness)
  melted_bm$Effectiveness[melted_bm$Effectiveness != "Unprotected" & melted_bm$DiffYear < 1] = "Unprotected"
  
  return(melted_bm)
  
}

rid_rare_outside_bm = function() { # Get rid of species with 0-4 occurrences in the reference level
  
  # Some values were dropped, make sure we still have at least N occurrences
  
  melted_sum = melted_bm |>
    filter(trunc_bm > 0) |>
    group_by(Species) |>
    summarise(trunc_bm = n()) |>
    filter(trunc_bm >= 30)
  
  spec_vec = as.character(melted_sum$Species)
  
  melted_filt_bm = melted_bm[melted_bm$Species %in% spec_vec,]
  
  # Get rid of rares outside
  
  melted_fsum = melted_filt_bm |>
    filter(Effectiveness == "Unprotected" & trunc_bm > 0) |>
    group_by(Species) |>
    summarise(trunc_bm = n()) |>
    filter(trunc_bm > 4) 
  
  spec_vec = as.character(melted_fsum$Species)
  
  melted_filt_bm = melted_filt_bm[melted_filt_bm$Species %in% spec_vec,]
  melted_filt_bm$Species = as.factor(as.character(melted_filt_bm$Species))
  
  save(melted_filt_bm, file = here::here("data", "data", "melted_filt_bm.RData"))
  
  return(melted_filt_bm)
}

add_covariates_bm = function(){
  
  all_cov = hab_filt |>
    left_join(socio_filt, by = "SurveyID") |>
    left_join(env_filt, by = "SurveyID") |>
    left_join(fine_habitat[,-c(2:4)], by = "SurveyID")
  
  melted_cov = melted_filt_bm |> 
    left_join(all_cov, by = "SurveyID")
  
  return(melted_cov)
  
}

compute_ES_bm = function() {
  
  melted_cov$Effectiveness = relevel(as.factor(melted_cov$Effectiveness), ref = "Unprotected")
  
  all_sp_glm_bm = lapply(1:length(unique(melted_cov$Species)), function(i) {
  # all_sp_glm_bm = lapply(101:150, function(i) {
    
    subdf = filter(melted_cov, Species == unique(melted_cov$Species)[i])
    
    pres_sp1 = subdf %>%
      filter(trunc_bm > 0)
    
    eco_weights = pres_sp1 |>
      group_by(Ecoregion) |> 
      mutate(weights = (1/n())) |> 
      ungroup()
    
    pca_prep = all_cov[all_cov$SurveyID %in% unique(pres_sp1$SurveyID),]
    pca_prep = pca_prep |> 
      left_join(eco_weights %>% dplyr::select(SurveyID, weights), by = "SurveyID")
    
    pca = ade4::dudi.pca(pca_prep[, -c(1, 97)],
                         center = T, scale = T, 
                         row.w = pca_prep$weights,
                         nf = 25, 
                         scannf = F)
    
    PC_coords = pca$li
    eigs = factoextra::get_eig(pca)
    cum_var = which(eigs[,3] > 70)[1]
    
    pca_prep_and_coords = cbind(pca_prep, PC_coords)
    pres_sp1 = pres_sp1[, 1:10] |> 
      left_join(pca_prep_and_coords[, c(1, 98:122)], by = "SurveyID") 
    
    istherefull = ifelse(length(pres_sp1$trunc_bm[which(pres_sp1$Effectiveness == "Fully Protected")]) > 9, T, F)
    istherehigh = ifelse(length(pres_sp1$trunc_bm[which(pres_sp1$Effectiveness == "Highly Protected")]) > 9, T, F)
    istherelight = ifelse(length(pres_sp1$trunc_bm[which(pres_sp1$Effectiveness == "Lightly Protected")]) > 9, T, F)
    isthereunp = ifelse(length(pres_sp1$trunc_bm[which(pres_sp1$Effectiveness == "Unprotected")]) > 9, T, F)
    
    pres_sp1_nofull = pres_sp1 %>%
      filter(Effectiveness != "Fully Protected")
    pres_sp1_nohigh = pres_sp1 %>%
      filter(Effectiveness != "Highly Protected")
    pres_sp1_nolight = pres_sp1 %>% 
      filter(Effectiveness != "Lightly Protected")
    pres_sp1_nofullhigh = pres_sp1 %>% 
      filter(Effectiveness != "Highly Protected" & Effectiveness != "Fully Protected")
    pres_sp1_nofulllight = pres_sp1 %>% 
      filter(Effectiveness != "Lightly Protected" & Effectiveness != "Fully Protected")
    pres_sp1_nohighlight = pres_sp1 %>% 
      filter(Effectiveness != "Lightly Protected" & Effectiveness != "Highly Protected")
    
    if(istherefull == F){pres_sp1 = pres_sp1_nofull}
    if(istherehigh == F){pres_sp1 = pres_sp1_nohigh}
    if(istherelight == F){pres_sp1 = pres_sp1_nolight}
    if(istherelight == F & istherefull == F){pres_sp1 = pres_sp1_nofulllight}
    if(istherehigh == F & istherefull == F){pres_sp1 = pres_sp1_nofullhigh}
    if(istherelight == F & istherehigh == F){pres_sp1 = pres_sp1_nohighlight}
    
    if(isthereunp == F){return(NULL)} 
    if(istherefull == F & istherehigh == F & istherelight == F){return(NULL)}
    
    names_var = names(pres_sp1[, c(4, 11:(10+cum_var))])
    formula_var = paste0(names_var, collapse = "+")
    full_formula = paste("log10(trunc_bm) ~ ", formula_var, sep = "")
    
    tryCatch({model = glm(formula = full_formula,
                          data = pres_sp1)}, error = function(e) model <<- NULL)
    
    if(is.null(model) == T){return(NULL)}
    sum_isthere = sum(istherefull, istherehigh, istherelight)
    if(sum_isthere < 2){vif = sqrt((regclass::VIF(model)))} else{vif = regclass::VIF(model)[, 3]}
    
    tryCatch({
      while(max(vif) > 2.23){ 
        
        u = which.max(vif)
        names_var = names_var[-u]
        formula_var = paste0(names_var, collapse = "+")
        full_formula = paste("log10(trunc_bm) ~ ", formula_var, sep = "")
        
        model = glm(formula = full_formula,
                    data = pres_sp1)
        
        if(sum_isthere < 2){vif = sqrt((regclass::VIF(model)))} else{vif = regclass::VIF(model)[, 3]}
        
        full_formula = full_formula
      }}, error = function(e) return(NULL))
    
    full_formula_sp = as.formula(paste(full_formula, "+ Matern(1 | SiteLongitude + SiteLatitude)", sep = ""))
    
    tryCatch({sp_model = spaMM::fitme(full_formula_sp,
                                      family = "gaussian",
                                      method = "PQL/L",
                                      data = pres_sp1)}, error = function(e) sp_model <<- NULL)

    model = glm(full_formula, data = pres_sp1)
    
    if(is.null(sp_model) == F){mod_coefs = as.data.frame(sp_model$fixef)}
    if(is.null(sp_model) == F){intercept = mod_coefs["(Intercept)",]}
    
    # mod_coefs = as.data.frame(model$coefficients)
    # intercept = mod_coefs["(Intercept)",]
    
    estimate_outside = ifelse(sum(stringr::str_detect(row.names(mod_coefs),"Intercept")) == 1,
                              10^(intercept),
                              NA)
    
    estimate_full = ifelse(sum(stringr::str_detect(row.names(mod_coefs),"EffectivenessFully Protected")) == 1,
                           10^(intercept+mod_coefs["EffectivenessFully Protected",]),
                           NA)
    
    estimate_high = ifelse(sum(stringr::str_detect(row.names(mod_coefs),"EffectivenessHighly Protected")) == 1,
                           10^(intercept+mod_coefs["EffectivenessHighly Protected",]),
                           NA)
    
    estimate_light = ifelse(sum(stringr::str_detect(row.names(mod_coefs),"EffectivenessLightly Protected")) == 1,
                            10^(intercept+mod_coefs["EffectivenessLightly Protected",]),
                            NA)
    
    IRR_full = estimate_full/estimate_outside
    IRR_high = estimate_high/estimate_outside
    IRR_light = estimate_light/estimate_outside
    
    tryCatch({p = spaMM::predict.HLfit(sp_model)}, error = function(e) p <<- NA)
    # p = predict(model)
    tryCatch({RSQ = (cor(10^(p), pres_sp1$trunc_bm))^2}, error = function(e) RSQ <<- NA)
    
    Species = unique(melted_cov$Species)[i]
    
    tryCatch({coef = data.frame(Species,
                                estimate_outside,
                                estimate_full, 
                                estimate_high,
                                estimate_light,
                                IRR_full,
                                IRR_high,
                                IRR_light,
                                RSQ)}, error = function(e) coef <<- NULL)
    
    cat(i, "\n")
    return(coef)
    
  }) 
  
  df_bm = data.table::rbindlist(all_sp_glm_bm)
  save(df_bm, file = here::here("outputs", "df_bm.RData"))
  df_product = df_bm %>% 
    left_join(OS_occu %>% 
                dplyr::select(Species, 
                              estimate_outside, 
                              estimate_full, 
                              estimate_high,
                              estimate_light, 
                              AUC), 
              by = "Species") %>% 
    rowwise() %>% 
    mutate(product_full = estimate_full.x*estimate_full.y,
           product_high = estimate_high.x*estimate_high.y,
           product_light = estimate_light.x*estimate_light.y,
           product_outside = estimate_outside.x*estimate_outside.y) %>% 
    mutate(IRR_product_full = product_full/product_outside,
           IRR_product_high = product_high/product_outside,
           IRR_product_light = product_light/product_outside) %>% 
    drop_na(estimate_outside.y)
  
  return(df_product)
}

make_ssmatrix_melted_ab = function(n) {
  
  tab <- table(fish_data$valid_name_FishBase)
  boxplot(tab, outline = F)
  tab <- tab[tab>n]
  fish_filt  <- fish_data[fish_data$valid_name_FishBase %in% names(tab),]
  fish_filt$Species = as.factor(fish_filt$valid_name_FishBase)
  
  fish_filt1 = fish_filt |>
    group_by(Species) |>
    mutate(percentile_99 = ifelse(max(Num) > 499, quantile(Num, 0.99), quantile(Num, 1))) |>
    filter(Num < percentile_99) |>
    ungroup()
  
  ssmatrix_ab = t(fossil::create.matrix(as.data.frame(fish_filt1), 
                                        tax.name = "Species", 
                                        locality = "SurveyID",
                                        abund = T,
                                        abund.col = "Num"))
  
  melted_ab = ssmatrix_ab |>
    reshape2::melt() |>
    rename(SurveyID = "Var1", 
           Species = "Var2", 
           Abundance = "value") |>
    left_join(covariates_corrected |> 
                dplyr::select(SurveyID, 
                              Effectiveness, 
                              DiffYear,
                              Area_km2), by = "SurveyID") |>
    left_join(sites_info[,c(1:4,10)], by = "SurveyID") |>
    drop_na(Effectiveness)
  
  melted_ab$Effectiveness = fct_collapse(melted_ab$Effectiveness, 
                                         "Highly Protected" = c("High No take", 
                                                                "High No take multizoned",
                                                                "Medium No take", 
                                                                "Medium No take multizoned"),
                                         "Lightly Protected" = c("High Restricted take", 
                                                                 "High Restricted take multizoned", 
                                                                 "Medium Restricted take",
                                                                 "Medium Restricted take multizoned",
                                                                 "Low No take", 
                                                                 "Low No take multizoned",
                                                                 "Low Restricted take",
                                                                 "Low Restricted take multizoned"),
                                         "Unprotected" = c("Medium Fishing", "Low Fishing", "out"))
  
  melted_ab$Effectiveness = as.character(melted_ab$Effectiveness)
  melted_ab$Effectiveness[melted_ab$Effectiveness == "Highly Protected" & melted_ab$DiffYear >= 10 & melted_ab$Area_km2 >= 100] = "Fully Protected"
  melted_ab$Effectiveness = as.factor(melted_ab$Effectiveness)
  melted_ab$Effectiveness[melted_ab$Effectiveness != "Unprotected" & melted_ab$DiffYear < 1] = "Unprotected"
  
  return(melted_ab)
  
}

rid_rare_outside_ab = function() { # Get rid of species with 0-4 occurrences in the reference level
  
  # Some values were dropped, make sure we still have at least N occurrences
  
  melted_sum = melted_ab |>
    filter(Abundance > 0) |>
    group_by(Species) |>
    summarise(Abundance = n()) |>
    filter(Abundance >= 30)
  
  spec_vec = as.character(melted_sum$Species)
  
  melted_filt_ab = melted_ab[melted_ab$Species %in% spec_vec,]
  
  # Get rid of rares outside
  
  melted_fsum = melted_filt_ab |>
    filter(Effectiveness == "Unprotected" & Abundance > 0) |>
    group_by(Species) |>
    summarise(Abundance = n()) |>
    filter(Abundance > 4) 
  
  spec_vec = as.character(melted_fsum$Species)
  
  melted_filt_ab = melted_filt_ab[melted_filt_ab$Species %in% spec_vec,]
  melted_filt_ab$Species = as.factor(as.character(melted_filt_ab$Species))
  
  save(melted_filt_ab, file = here::here("data", "data", "melted_filt_ab.RData"))
  
  return(melted_filt_ab)
}

add_covariates_ab = function(){
  
  all_cov = hab_filt |>
    left_join(socio_filt, by = "SurveyID") |>
    left_join(env_filt, by = "SurveyID") |>
    left_join(fine_habitat[,-c(2:4)], by = "SurveyID")
  
  melted_cov = melted_filt_ab |> 
    left_join(all_cov, by = "SurveyID")
  
  return(melted_cov)
  
}

compute_ES_ab = function() {
  
  melted_cov$Effectiveness = relevel(as.factor(melted_cov$Effectiveness), ref = "Unprotected")
  
  all_sp_glm_ab = lapply(1:length(unique(melted_cov$Species)), function(i) {
    # all_sp_glm_ab = lapply(101:150, function(i) {
    
    subdf = filter(melted_cov, Species == unique(melted_cov$Species)[i])
    
    pres_sp1 = subdf %>%
      filter(Abundance > 0)
    
    eco_weights = pres_sp1 |>
      group_by(Ecoregion) |> 
      mutate(weights = (1/n())) |> 
      ungroup()
    
    pca_prep = all_cov[all_cov$SurveyID %in% unique(pres_sp1$SurveyID),]
    pca_prep = pca_prep |> 
      left_join(eco_weights %>% dplyr::select(SurveyID, weights), by = "SurveyID")
    
    pca = ade4::dudi.pca(pca_prep[, -c(1, 97)],
                         center = T, scale = T, 
                         row.w = pca_prep$weights,
                         nf = 25, 
                         scannf = F)
    
    PC_coords = pca$li
    eigs = factoextra::get_eig(pca)
    cum_var = which(eigs[,3] > 70)[1]
    
    pca_prep_and_coords = cbind(pca_prep, PC_coords)
    pres_sp1 = pres_sp1[, 1:10] |> 
      left_join(pca_prep_and_coords[, c(1, 98:122)], by = "SurveyID") 
    
    istherefull = ifelse(length(pres_sp1$Abundance[which(pres_sp1$Effectiveness == "Fully Protected")]) > 9, T, F)
    istherehigh = ifelse(length(pres_sp1$Abundance[which(pres_sp1$Effectiveness == "Highly Protected")]) > 9, T, F)
    istherelight = ifelse(length(pres_sp1$Abundance[which(pres_sp1$Effectiveness == "Lightly Protected")]) > 9, T, F)
    isthereunp = ifelse(length(pres_sp1$Abundance[which(pres_sp1$Effectiveness == "Unprotected")]) > 9, T, F)
    
    pres_sp1_nofull = pres_sp1 %>%
      filter(Effectiveness != "Fully Protected")
    pres_sp1_nohigh = pres_sp1 %>%
      filter(Effectiveness != "Highly Protected")
    pres_sp1_nolight = pres_sp1 %>% 
      filter(Effectiveness != "Lightly Protected")
    pres_sp1_nofullhigh = pres_sp1 %>% 
      filter(Effectiveness != "Highly Protected" & Effectiveness != "Fully Protected")
    pres_sp1_nofulllight = pres_sp1 %>% 
      filter(Effectiveness != "Lightly Protected" & Effectiveness != "Fully Protected")
    pres_sp1_nohighlight = pres_sp1 %>% 
      filter(Effectiveness != "Lightly Protected" & Effectiveness != "Highly Protected")
    
    if(istherefull == F){pres_sp1 = pres_sp1_nofull}
    if(istherehigh == F){pres_sp1 = pres_sp1_nohigh}
    if(istherelight == F){pres_sp1 = pres_sp1_nolight}
    if(istherelight == F & istherefull == F){pres_sp1 = pres_sp1_nofulllight}
    if(istherehigh == F & istherefull == F){pres_sp1 = pres_sp1_nofullhigh}
    if(istherelight == F & istherehigh == F){pres_sp1 = pres_sp1_nohighlight}
    
    if(isthereunp == F){return(NULL)} 
    if(istherefull == F & istherehigh == F & istherelight == F){return(NULL)}
    
    names_var = names(pres_sp1[, c(4, 11:(10+cum_var))])
    formula_var = paste0(names_var, collapse = "+")
    full_formula = paste("Abundance ~ ", formula_var, sep = "")
    
    tryCatch({model = MASS::glm.nb(formula = full_formula,
                                   data = pres_sp1)}, error = function(e) model <<- NULL)
    
    if(is.null(model) == T){return(NULL)}
    sum_isthere = sum(istherefull, istherehigh, istherelight)
    if(sum_isthere < 2){vif = sqrt((regclass::VIF(model)))} else{vif = regclass::VIF(model)[, 3]}
    
    tryCatch({
      while(max(vif) > 2.23){ 
        
        u = which.max(vif)
        names_var = names_var[-u]
        formula_var = paste0(names_var, collapse = "+")
        full_formula = paste("Abundance ~ ", formula_var, sep = "")
        
        model = glm(formula = full_formula,
                    data = pres_sp1)
        
        if(sum_isthere < 2){vif = sqrt((regclass::VIF(model)))} else{vif = regclass::VIF(model)[, 3]}
        
        full_formula = full_formula
      }}, error = function(e) return(NULL))
    
    full_formula_sp = as.formula(paste(full_formula, "+ Matern(1 | SiteLongitude + SiteLatitude)", sep = ""))
    
    tryCatch({sp_model = spaMM::fitme(full_formula_sp,
                                      family = "negbin",
                                      method = "PQL/L",
                                      data = pres_sp1)}, error = function(e) sp_model <<- NULL)
    
    if(is.null(sp_model) == F){mod_coefs = as.data.frame(sp_model$fixef)}
    if(is.null(sp_model) == F){intercept = mod_coefs["(Intercept)",]}
    
    estimate_outside = ifelse(sum(stringr::str_detect(row.names(mod_coefs),"Intercept")) == 1,
                              exp(intercept),
                              NA)
    
    estimate_full = ifelse(sum(stringr::str_detect(row.names(mod_coefs),"EffectivenessFully Protected")) == 1,
                           exp(intercept+mod_coefs["EffectivenessFully Protected",]),
                           NA)
    
    estimate_high = ifelse(sum(stringr::str_detect(row.names(mod_coefs),"EffectivenessHighly Protected")) == 1,
                           exp(intercept+mod_coefs["EffectivenessHighly Protected",]),
                           NA)
    
    estimate_light = ifelse(sum(stringr::str_detect(row.names(mod_coefs),"EffectivenessLightly Protected")) == 1,
                            exp(intercept+mod_coefs["EffectivenessLightly Protected",]),
                            NA)
    
    IRR_full = estimate_full/estimate_outside
    IRR_high = estimate_high/estimate_outside
    IRR_light = estimate_light/estimate_outside
    
    tryCatch({p = spaMM::predict.HLfit(sp_model, type = "response")}, error = function(e) p <<- NA)
    tryCatch({RSQ = cor(p, pres_sp1$Abundance)^2}, error = function(e) RSQ <<- NA)
    
    Species = unique(melted_cov$Species)[i]
    
    tryCatch({coef = data.frame(Species,
                                estimate_outside,
                                estimate_full, 
                                estimate_high,
                                estimate_light,
                                IRR_full,
                                IRR_high,
                                IRR_light,
                                RSQ)}, error = function(e) coef <<- NULL)
    
    cat(i, "\n")
    
    return(coef)
    
  }) 
  
  
  df_ab = data.table::rbindlist(all_sp_glm_ab)
  save(df_ab, file = here::here("outputs", "df_ab.RData"))
  df_product_ab = df_ab %>% 
    left_join(OS_occu %>% 
                dplyr::select(Species, 
                              estimate_outside, 
                              estimate_full, 
                              estimate_high,
                              estimate_light, 
                              AUC), 
              by = "Species") %>% 
    rowwise() %>% 
    mutate(product_full = estimate_full.x*estimate_full.y,
           product_high = estimate_high.x*estimate_high.y,
           product_light = estimate_light.x*estimate_light.y,
           product_outside = estimate_outside.x*estimate_outside.y) %>% 
    mutate(IRR_product_full = product_full/product_outside,
           IRR_product_high = product_high/product_outside,
           IRR_product_light = product_light/product_outside) %>% 
    drop_na(estimate_outside.y)
  
  return(df_product_ab)
  
}

compute_ES_ab_wzero = function() {
  
  melted_cov$Effectiveness = relevel(as.factor(melted_cov$Effectiveness), ref = "Unprotected")
  
  all_sp_glm_ab = lapply(1:length(unique(melted_cov$Species)), function(i) {
    # all_sp_glm_ab = lapply(101:150, function(i) {
    
    subdf = filter(melted_cov, Species == unique(melted_cov$Species)[i])
    
    pres_sp1 = subdf %>%
      filter(Abundance > 0)
    
    pres_sp1$Ecoregion = as.factor(as.character(pres_sp1$Ecoregion))
    
    eco_table = table(pres_sp1$Ecoregion)
    eco_tab = names(eco_table[eco_table > 4])
    
    abs_sp1 = subdf %>%
      filter(Abundance < 1)
    
    abs_sp1 = abs_sp1[abs_sp1$Ecoregion %in% eco_tab,]
    pres_sp1 = pres_sp1[pres_sp1$Ecoregion %in% eco_tab,]
    
    presabs_sp1 = rbind(pres_sp1, abs_sp1)
    presabs_sp1 = pres_sp1
    
    eco_weights = presabs_sp1 |>
      group_by(Ecoregion) |> 
      mutate(weights = (1/n())) |> 
      ungroup()
    
    pca_prep = all_cov[all_cov$SurveyID %in% unique(presabs_sp1$SurveyID),]
    pca_prep = pca_prep |> 
      left_join(eco_weights %>% dplyr::select(SurveyID, weights), by = "SurveyID")
    
    pca = ade4::dudi.pca(pca_prep[, -c(1, 97)],
                         center = T, scale = T, 
                         row.w = pca_prep$weights,
                         nf = 25, 
                         scannf = F)
    
    PC_coords = pca$li
    eigs = factoextra::get_eig(pca)
    cum_var = which(eigs[,3] > 70)[1]
    
    pca_prep_and_coords = cbind(pca_prep, PC_coords)
    pres_sp1 = presabs_sp1[, 1:10] |> 
      left_join(pca_prep_and_coords[, c(1, 98:122)], by = "SurveyID") 
    
    istherefull = ifelse(length(pres_sp1$Abundance[which(pres_sp1$Effectiveness == "Fully Protected")]) > 9, T, F)
    istherehigh = ifelse(length(pres_sp1$Abundance[which(pres_sp1$Effectiveness == "Highly Protected")]) > 9, T, F)
    istherelight = ifelse(length(pres_sp1$Abundance[which(pres_sp1$Effectiveness == "Lightly Protected")]) > 9, T, F)
    isthereunp = ifelse(length(pres_sp1$Abundance[which(pres_sp1$Effectiveness == "Unprotected")]) > 9, T, F)
    
    pres_sp1_nofull = pres_sp1 %>%
      filter(Effectiveness != "Fully Protected")
    pres_sp1_nohigh = pres_sp1 %>%
      filter(Effectiveness != "Highly Protected")
    pres_sp1_nolight = pres_sp1 %>% 
      filter(Effectiveness != "Lightly Protected")
    pres_sp1_nofullhigh = pres_sp1 %>% 
      filter(Effectiveness != "Highly Protected" & Effectiveness != "Fully Protected")
    pres_sp1_nofulllight = pres_sp1 %>% 
      filter(Effectiveness != "Lightly Protected" & Effectiveness != "Fully Protected")
    pres_sp1_nohighlight = pres_sp1 %>% 
      filter(Effectiveness != "Lightly Protected" & Effectiveness != "Highly Protected")
    
    if(istherefull == F){pres_sp1 = pres_sp1_nofull}
    if(istherehigh == F){pres_sp1 = pres_sp1_nohigh}
    if(istherelight == F){pres_sp1 = pres_sp1_nolight}
    if(istherelight == F & istherefull == F){pres_sp1 = pres_sp1_nofulllight}
    if(istherehigh == F & istherefull == F){pres_sp1 = pres_sp1_nofullhigh}
    if(istherelight == F & istherehigh == F){pres_sp1 = pres_sp1_nohighlight}
    
    if(isthereunp == F){return(NULL)} 
    if(istherefull == F & istherehigh == F & istherelight == F){return(NULL)}
    
    names_var = names(pres_sp1[, c(4, 11:(10+cum_var))])
    formula_var = paste0(names_var, collapse = "+")
    full_formula = paste("Abundance ~ ", formula_var, sep = "")
    
    tryCatch({model = MASS::glm.nb(formula = full_formula,
                                   data = pres_sp1)}, error = function(e) model <<- NULL)
    
    if(is.null(model) == T){return(NULL)}
    sum_isthere = sum(istherefull, istherehigh, istherelight)
    if(sum_isthere < 2){vif = sqrt((regclass::VIF(model)))} else{vif = regclass::VIF(model)[, 3]}
    
    tryCatch({
      while(max(vif) > 2.23){ 
        
        u = which.max(vif)
        names_var = names_var[-u]
        formula_var = paste0(names_var, collapse = "+")
        full_formula = paste("Abundance ~ ", formula_var, sep = "")
        
        model = MASS::glm.nb(formula = full_formula,
                    data = pres_sp1)
        
        if(sum_isthere < 2){vif = sqrt((regclass::VIF(model)))} else{vif = regclass::VIF(model)[, 3]}
        
        full_formula = full_formula
      }}, error = function(e) return(NULL))
    
    full_formula_sp = as.formula(paste(full_formula, "+ Matern(1 | SiteLongitude + SiteLatitude)", sep = ""))
    
    tryCatch({sp_model = spaMM::fitme(full_formula_sp,
                                      family = "negbin",
                                      method = "PQL/L",
                                      data = pres_sp1)}, error = function(e) sp_model <<- NULL)
    
    if(is.null(sp_model) == F){mod_coefs = as.data.frame(sp_model$fixef)}
    if(is.null(sp_model) == F){intercept = mod_coefs["(Intercept)",]}
    
    estimate_outside = ifelse(sum(stringr::str_detect(row.names(mod_coefs),"Intercept")) == 1,
                              exp(intercept),
                              NA)
    
    estimate_full = ifelse(sum(stringr::str_detect(row.names(mod_coefs),"EffectivenessFully Protected")) == 1,
                           exp(intercept+mod_coefs["EffectivenessFully Protected",]),
                           NA)
    
    estimate_high = ifelse(sum(stringr::str_detect(row.names(mod_coefs),"EffectivenessHighly Protected")) == 1,
                           exp(intercept+mod_coefs["EffectivenessHighly Protected",]),
                           NA)
    
    estimate_light = ifelse(sum(stringr::str_detect(row.names(mod_coefs),"EffectivenessLightly Protected")) == 1,
                            exp(intercept+mod_coefs["EffectivenessLightly Protected",]),
                            NA)
    
    IRR_full = estimate_full/estimate_outside
    IRR_high = estimate_high/estimate_outside
    IRR_light = estimate_light/estimate_outside
    
    tryCatch({p = spaMM::predict.HLfit(sp_model, type = "response")}, error = function(e) p <<- NA)
    tryCatch({p = predict(model, type = "response")}, error = function(e) p <<- NA)
    tryCatch({RSQ = cor(p, pres_sp1$Abundance)^2}, error = function(e) RSQ <<- NA)
    
    Species = unique(melted_cov$Species)[i]
    
    tryCatch({coef = data.frame(Species,
                                estimate_outside,
                                estimate_full, 
                                estimate_high,
                                estimate_light,
                                IRR_full,
                                IRR_high,
                                IRR_light,
                                RSQ)}, error = function(e) coef <<- NULL)
    
    cat(i, "\n")
    
    return(coef)
    
  }) 
  
  
  df_ab = data.table::rbindlist(all_sp_glm_ab)
  save(df_ab, file = here::here("outputs", "df_ab.RData"))
  df_product_ab = df_ab %>% 
    left_join(OS_occu %>% 
                dplyr::select(Species, 
                              estimate_outside, 
                              estimate_full, 
                              estimate_high,
                              estimate_light, 
                              AUC), 
              by = "Species") %>% 
    rowwise() %>% 
    mutate(product_full = estimate_full.x*estimate_full.y,
           product_high = estimate_high.x*estimate_high.y,
           product_light = estimate_light.x*estimate_light.y,
           product_outside = estimate_outside.x*estimate_outside.y) %>% 
    mutate(IRR_product_full = product_full/product_outside,
           IRR_product_high = product_high/product_outside,
           IRR_product_light = product_light/product_outside) %>% 
    drop_na(estimate_outside.y)
  
  return(df_product_ab)
  
}

#### Histograms ####

histo_occu = OS_occu %>%
  rename(Full = RR_full,
         High = RR_high, 
         Light = RR_light) %>% 
  dplyr::select(Species, Full, High, Light) %>% 
  reshape2::melt() %>% 
  mutate(log_RR = log(value)) %>% 
  drop_na(value)

ggplot(histo_occu, aes(x = log_RR, y = variable, fill = stat(x))) + 
  ggridges:: geom_density_ridges_gradient(scale = 2, size = 1) +
  scale_fill_viridis_c(name = "log_RR", option = "C") +
  theme_minimal()

histo_abun = df_ab %>% 
  dplyr::select(Species, IRR_full, IRR_high, IRR_light) %>%
  rename(Full = IRR_full,
         High = IRR_high, 
         Light = IRR_light) %>%   
  reshape2::melt() %>% 
  mutate(log_RR = log(value)) %>% 
  drop_na(value)

ggplot(histo_abun, aes(x = log_RR, y = variable, fill = stat(x))) + 
  ggridges:: geom_density_ridges_gradient(scale = 2, size = 1) +
  scale_fill_viridis_c(name = "log_RR", option = "C") + 
  theme_minimal()

histo_biom = df_product %>% 
  dplyr::select(Species, IRR_product_full, IRR_product_high, IRR_product_light) %>%
  rename(Full = IRR_product_full,
         High = IRR_product_high, 
         Light = IRR_product_light) %>%   
  reshape2::melt() %>% 
  mutate(log_RR = log(value)) %>% 
  drop_na(value)

ggplot(histo_biom, aes(x = log_RR, y = variable, fill = stat(x))) + 
  ggridges:: geom_density_ridges_gradient(scale = 2, size = 1) +
  scale_fill_viridis_c(name = "log_RR", option = "C") + 
  theme_minimal()

histo_occu$Type = "Occurrence"
histo_abun$Type = "Abundance"
histo_biom$Type = "Biomass"

all_hist = rbind(histo_occu, histo_abun, histo_biom) 
all_hist = all_hist %>% 
  mutate(Type = fct_relevel(Type, "Occurrence", "Abundance", "Biomass")) %>% 
  mutate(variable = fct_relevel(variable, "Full", "High", "Light"))

ggplot(all_hist, aes(x = log_RR, y = variable, fill = stat(x))) +
  ggridges::geom_density_ridges_gradient(scale = 1.5, size = 0.5, alpha = 0.8) +
  scale_fill_viridis_c(name = "log_RR", option = "C") + 
  facet_wrap(~Type) +
  xlim(-3, 3) +
  geom_vline(xintercept = 0) +
  theme_minimal()

ggplot(all_hist, aes(x = log_RR, y = Type)) +
  ggridges::geom_density_ridges_gradient(scale = 1.5, 
                                         size = 0.5, 
                                         alpha = 0.8,
                                         show.legend = F) +
  facet_wrap(~variable) +
  xlim(-3, 3) +
  geom_vline(xintercept = 0) +
  xlab("Log(Effect size)") +
  theme_minimal()

ggplot(all_hist, aes(y = log_RR, x = variable, fill = variable)) + 
  geom_violin(alpha = 0.75) +
  theme(legend.position = 'top', 
        legend.spacing.x = unit(1.0, 'cm')) +
  scale_fill_manual(values = c("#092147", "#1A488E", "#A0BACC")) 
  


  
  

all_hist$category <- as.factor(paste(all_hist$variable, all_hist$Type, sep="_"))
levels(all_hist$category)

# GRAPH NICO

F_H = ggplot(all_hist %>% filter(variable != "Light"), aes(y = log_RR, x = Type, fill = variable)) + 
  geom_boxplot(outlier.shape = NA, 
               show.legend = F,
               width = 0.5,
               alpha = 0.75) + 
  facefuns::geom_split_violin(trim = F, 
                              show.legend = F, 
                              alpha = 0.5) +
  theme_minimal() + 
  xlab("Full vs. High") +
  ylab("Effect size of protection") +
  ylim(-2, 2) +
  theme(axis.title.x.top = element_text(vjust = +1),
        axis.text.x.top = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_fill_manual(values = c("#092147", "#1A488E")) +
  scale_x_discrete(position = "top")


F_L = ggplot(all_hist %>% filter(variable != "High"), aes(y = log_RR, x = Type, fill = variable)) + 
  geom_boxplot(outlier.shape = NA, 
               alpha = 0.75,
               width = 0.5,
               show.legend = F) + 
  facefuns::geom_split_violin(trim = F, 
                              show.legend = F, 
                              alpha = 0.5) +
  xlab("Full vs. Light") +
  theme_minimal() + 
  ylim(-2, 2) +
  theme(axis.title.y = element_blank(),
        axis.title.x.top = element_text(vjust = +1),
        axis.text.y = element_blank(),
        axis.text.x.top = element_text(size = 10),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_fill_manual(values = c("#092147", "#A0BACC")) +
  scale_x_discrete(position = "top")
  
H_L = ggplot(all_hist %>% filter(variable != "Full"), aes(y = log_RR, x = Type, fill = variable)) + 
  geom_boxplot(outlier.shape = NA, 
               alpha = 0.75,
               width = 0.5,
               show.legend = F) + 
  facefuns::geom_split_violin(trim = F, 
                              show.legend = F, 
                              alpha = 0.5) +
  xlab("High vs. Light") +
  theme_minimal() + 
  ylim(-2, 2) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x.top = element_text(vjust = +1),
        axis.text.x.top = element_text(size = 10),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_fill_manual(values = c("#1A488E", "#A0BACC")) +
  scale_x_discrete(position = "top")

all_facets = gridExtra::grid.arrange(F_H, 
                                     F_L, 
                                     H_L, 
                                     ncol = 3)

all_facets = cowplot::plot_grid(F_H, F_L, H_L, ncol = 3)

#### Models ####

melted_rares = melted_cov %>%
  filter(Presence > 0) %>% 
  group_by(Species) %>% 
  summarise(n_occu = n()) %>% 
  mutate(decile_rank = ntile(n_occu, 10)) %>% 
  # mutate(Rarity = ifelse(decile_rank > 5, "Common", "Rare"))
  mutate(Rarity = ifelse(n_occu > 100, "Common", "Rare"))

ggplot(melted_rares, aes(x = n_occu, colour = "black")) + 
  geom_histogram(binwidth = 5) + 
  xlim(0, 500)

# melted_rares_pres = melted_cov %>%
#   filter(Presence > 0) %>% 
#   group_by(Species) %>% 
#   summarise(n_occu = n())  
# 
# melted_rares_abs = melted_cov %>%
#   group_by(Species, Ecoregion) %>%
#   filter(Presence > 0) %>% 
#   summarise(n_occu = n()) %>% 
#   filter(n_occu > 0) 
# 
# melted_rares_all = melted_cov %>% 
#   group_by(Species, Ecoregion) %>% 
#   summarise(n_all = n())
# 
# melted_rares_presabs = melted_rares_abs %>% 
#   left_join(melted_rares_all, by = c("Species", "Ecoregion")) %>% 
#   group_by(Species) %>% 
#   summarise(n_occu_sum = sum(n_occu),
#             n_all_sum = sum(n_all)) %>% 
#   mutate(Rate = n_occu_sum/n_all_sum) %>% 
#   mutate(decile_rank = ntile(Rate, 10)) %>% 
#   mutate(Rarity = ifelse(decile_rank > 5, "Common", "Rare"))

# IUCN_status_detailled$Species = gsub("_", " ", IUCN_status_detailled$species)
# IUCN_status_detailled$IUCN_status[is.na(IUCN_status_detailled$IUCN_status)] = "DD"
# IUCN_status_detailled$IUCN_status = relevel(as.factor(IUCN_status_detailled$IUCN_status), ref = "NT")
# IUCN_status_detailled$IUCN_status1 = ifelse(IUCN_status_detailled$IUCN_status == "CR" | 
#                                              IUCN_status_detailled$IUCN_status == "EN" | 
#                                              IUCN_status_detailled$IUCN_status == "VU", "Thr", "DD/NThr")
# IUCN_status_detailled$IUCN_status1 = relevel(as.factor(IUCN_status_detailled$IUCN_status1), ref = "DD/NThr")

traits$Species = traits$CURRENT_SPECIES_NAME
# traits = traits %>% left_join(IUCN_status_detailled, by = "Species")

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

# df_ab_traits = df_product_ab %>%
#   mutate(log_full = log(IRR_product_full),
#          log_high = log(IRR_product_high),
#          log_light = log(IRR_product_light)) %>%
#   left_join(traits, by = "Species") %>%
#   left_join(melted_rares_presabs %>% dplyr::select(Species, Rarity, decile_rank), by = "Species")  %>%
#   mutate(log_length = log(MaxLength))

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

df_bm_traits = df_product %>%
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

#### OCCU ####

mod_occu_full = glm(log_full ~ (Trophic.Level * log_length) * Rarity,
                    data = df_occu_full)
mod_occu_high = glm(log_high ~ (Trophic.Level * log_length) * Rarity,
                    data = df_occu_high)
mod_occu_light = glm(log_light ~ (Trophic.Level * log_length) * Rarity,
                     data = df_occu_light)

mod_occu_full_AIC = MASS::stepAIC(mod_occu_full)
mod_occu_high_AIC = MASS::stepAIC(mod_occu_high)
mod_occu_light_AIC = MASS::stepAIC(mod_occu_light)

summary(mod_occu_full_AIC)
summary(mod_occu_high_AIC)
summary(mod_occu_light_AIC)

res = DHARMa::simulateResiduals(mod_occu_full_AIC)
plot(res)
p = predict(mod_occu_full_AIC, type = "response")
RSQ = cor(p, df_occu_full$log_full)

res = DHARMa::simulateResiduals(mod_occu_high_AIC)
plot(res)
p = predict(mod_occu_high_AIC, type = "response")
RSQ = cor(p, df_occu_high$log_high)

res = DHARMa::simulateResiduals(mod_occu_light_AIC)
plot(res)
p = predict(mod_occu_light_AIC, type = "response")
RSQ = cor(p, df_occu_light$log_light)

interaction_occu_full = visreg::visreg(mod_occu_full_AIC, 
                                  type = "conditional",
                                  "Trophic.Level",
                                  by = "Rarity",
                                  trans = exp, 
                                  gg = F,
                                  nn = 2000)$fit

interaction_occu_full$group = "Occurrence"
interaction_occu_full$Protection = "Full"

ggplot(interaction_occu_full, aes(x = Trophic.Level, y = visregFit)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.5) +
  facet_wrap(~Rarity) + 
  ylab("Effect size on occurrence of full protection") +
  xlab("Trophic level") +
  theme_minimal() +
  theme(strip.text.x = element_text(size = 10, face = "bold"))

#### ABUNDANCE ####

mod_ab_full = glm(log_full ~ (Trophic.Level * log_length) * Rarity,
                    data = df_ab_full)
mod_ab_high = glm(log_high ~ (Trophic.Level * log_length) * Rarity,
                    data = df_ab_high)
mod_ab_light = glm(log_light ~ (Trophic.Level * log_length) * Rarity,
                     data = df_ab_light)

mod_ab_full_AIC = MASS::stepAIC(mod_ab_full)
mod_ab_high_AIC = MASS::stepAIC(mod_ab_high)
mod_ab_light_AIC = MASS::stepAIC(mod_ab_light)

summary(mod_ab_full_AIC) # RarityRare is higher because Rare species here = high trophic.levels !
hist(df_ab_full$Trophic.Level[df_ab_full$Rarity == "Rare"]) # With a low number of species modeled
summary(mod_ab_high_AIC)
summary(mod_ab_light_AIC)

res = DHARMa::simulateResiduals(mod_ab_full_AIC)
plot(res)
p = predict(mod_ab_full_AIC, type = "response")
RSQ = cor(p, df_ab_full$log_full)

res = DHARMa::simulateResiduals(mod_ab_high_AIC)
plot(res)
p = predict(mod_ab_high_AIC, type = "response")
RSQ = cor(p, df_ab_high$log_high)

res = DHARMa::simulateResiduals(mod_ab_light_AIC)
plot(res)
p = predict(mod_ab_light_AIC, type = "response")
RSQ = cor(p, df_ab_light$log_light)

# interaction_ab_full = visreg::visreg(mod_ab_full_AIC, 
#                                      type = "conditional",
#                                      "Trophic.Level",
#                                      by = "Rarity",
#                                      trans = exp, 
#                                      gg = F,
#                                      nn = 500)$fit
# 
# ggplot(interaction_ab_full, aes(x = Trophic.Level, y = visregFit)) +
#   geom_line(linewidth = 1) +
#   geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.5) +
#   facet_wrap(~Rarity) + 
#   ylab("Effect Size") +
#   theme_minimal()

interaction_ab_high = visreg::visreg(mod_ab_high_AIC, 
                                     type = "conditional",
                                     "Trophic.Level",
                                     by = "Rarity",
                                     trans = exp, 
                                     gg = F,
                                     nn = 2000)$fit

interaction_ab_high$group = "Abundance"
interaction_ab_high$Protection = "High"

ggplot(interaction_ab_high, aes(x = Trophic.Level, y = visregFit)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.5) +
  facet_wrap(~Rarity) + 
  ylab("Effect size on abundance of high protection") +
  xlab("Trophic level") +
  theme_minimal() +
  theme(strip.text.x = element_text(size = 10, face = "bold"))

#### BIOMASS ####

mod_bm_full = glm(log_full ~ (Trophic.Level * log_length) * Rarity,
                  data = df_bm_full)
mod_bm_high = glm(log_high ~ (Trophic.Level * log_length) * Rarity,
                  data = df_bm_high)
mod_bm_light = glm(log_light ~ (Trophic.Level * log_length) * Rarity,
                   data = df_bm_light)

mod_bm_full_AIC = MASS::stepAIC(mod_bm_full)
mod_bm_high_AIC = MASS::stepAIC(mod_bm_high)
mod_bm_light_AIC = MASS::stepAIC(mod_bm_light)

summary(mod_bm_full_AIC)
summary(mod_bm_high_AIC)
summary(mod_bm_light_AIC)

res = DHARMa::simulateResiduals(mod_bm_full)
plot(res)
p = predict(mod_bm_full, type = "response")
RSQ = cor(p, df_bm_full$log_full)

res = DHARMa::simulateResiduals(mod_bm_high)
plot(res)
p = predict(mod_bm_high, type = "response")
RSQ = cor(p, df_bm_high$log_high)

res = DHARMa::simulateResiduals(mod_bm_light)
plot(res)
p = predict(mod_bm_light, type = "response")
RSQ = cor(p, df_bm_light$log_light)


interaction_bm_full = visreg::visreg(mod_bm_full_AIC, 
                                     type = "conditional",
                                     "log_length",
                                     by = "Trophic.Level",
                                     trans = exp, 
                                     gg = F,
                                     nn = 2000)$fit

supp.labs <- c("2" = "Low","3.3" = "Intermediate", "4" = "High")
interaction_bm_full$Trophic.Level = as.factor(interaction_bm_full$Trophic.Level)

ggplot(interaction_bm_full, aes(x = log_length, y = visregFit)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), 
              alpha = 0.5, 
              fill = "#0d0887") +
  ylab("Effect size on biomass of full protection") +
  xlab("Maximum length (log)") +
  facet_grid(~Trophic.Level, labeller = labeller(Trophic.Level = supp.labs)) +
  theme_minimal() +
  theme(strip.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(vjust = +3),
        axis.title.x = element_text(vjust = -0.75))


interaction_bm_full2 = visreg::visreg(mod_bm_full_AIC, 
                                  type = "conditional",
                                  "Trophic.Level",
                                  by = "Rarity",
                                  trans = exp, 
                                  gg = F,
                                  nn = 2000)$fit

interaction_bm_full2$group = "Biomass"
interaction_bm_full2$Protection = "Full"

ggplot(interaction_bm_full2, aes(x = Trophic.Level, y = visregFit)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.5) +
  ylab("Effect size on biomass of high protection") +
  xlab("Trophic level") +
  facet_grid(~Rarity) +
  theme_minimal() +
  theme(strip.text.x = element_text(size = 10, face = "bold"))

interaction_bm_high_rare = visreg::visreg(mod_bm_high_AIC, 
                                          type = "conditional",
                                          "log_length",
                                          by = "Trophic.Level",
                                          cond = list(Rarity = "Rare"),
                                          trans = exp, 
                                          gg = F,
                                          nn = 2000)$fit
interaction_bm_high_rare$Rarity = "Rare"

interaction_bm_high_common = visreg::visreg(mod_bm_high_AIC, 
                                            type = "conditional",
                                            "log_length",
                                            by = "Trophic.Level",
                                            cond = list(Rarity = "Common"),
                                            trans = exp, 
                                            gg = F,
                                            nn = 2000)$fit
interaction_bm_high_common$Rarity = "Common"

both_bm_high = rbind(interaction_bm_high_common, interaction_bm_high_rare)

supp.labs <- c("2" = "Low trophic level","3.3" = "Intermediate trophic level", "4" = "High trophic level")

ggplot(both_bm_high, aes(x = log_length, y = visregFit)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.5) +
  ylab("Effect size on biomass of high protection") +
  xlab("Maximum length (log)") +
  facet_grid(Rarity ~ Trophic.Level, scales = 'free', labeller = labeller(Trophic.Level = supp.labs)) +
  theme_minimal() +
  theme(strip.text = element_text(size = 10, face = "bold"),
        panel.spacing = unit(0.5, "lines"),
        axis.title.y = element_text(vjust = +3),
        axis.title.x = element_text(vjust = -0.75))

# Group all Rarity:Trophic.Level

all_rarity_trophic = rbind(interaction_occu_full %>% dplyr::select(-c(log_length, log_full)),
                           interaction_ab_high %>% dplyr::select(-c(log_high)),
                           interaction_bm_full2 %>% dplyr::select(-c(log_length, log_full)))

all_rarity_trophic = all_rarity_trophic %>% 
  mutate(group = fct_relevel(group, "Occurrence", "Abundance", "Biomass"))

# ggplot(all_rarity_trophic, aes(x = Trophic.Level, y = visregFit, group = group)) +
#   geom_line(aes(linetype = Protection, color = group), linewidth = 1, alpha = 1) +
#   ylab("Effect size of protection") +
#   xlab("Trophic level") +
#   facet_grid(~Rarity) +
#   theme_minimal() +
#   theme(strip.text = element_text(size = 10, face = "bold"),
#         panel.spacing = unit(0.5, "lines"),
#         axis.title.y = element_text(vjust = +3),
#         axis.title.x = element_text(vjust = -0.75))

ggplot(all_rarity_trophic, aes(x = Trophic.Level, y = visregFit, group = group)) +
  geom_line(linewidth = 1,
            show.legend = F) +
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr, fill = group), 
              alpha = 0.75,
              show.legend = F) +
  ylab("Effect size of protection") +
  xlab("Trophic level") +
  facet_grid(group~Rarity, scales = "fixed") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12, face = "bold"),
        panel.spacing = unit(0.5, "lines"),
        axis.title.y = element_text(vjust = +2.65, 
                                    size = 15),
        axis.title.x = element_text(vjust = -0.75,
                                    size = 15),
        axis.text = element_text(size = 12)) +
  scale_fill_manual(values = c("#092147", "#1A488E", "#092147")) +
  labs(fill = "Type")

#### correlations ####

OF_OH = cor(OS_occu$RR_full, OS_occu$RR_high, use = "pairwise.complete.obs")
OF_OL = cor(OS_occu$RR_full, OS_occu$RR_light, use = "pairwise.complete.obs")
OH_OL = cor(OS_occu$RR_light, OS_occu$RR_high, use = "pairwise.complete.obs")

AF_AH = cor(df_ab$IRR_full, df_ab$IRR_high, use = "pairwise.complete.obs")
AF_AL = cor(df_ab$IRR_full, df_ab$IRR_light, use = "pairwise.complete.obs")
AH_AL = cor(df_ab$IRR_light, df_ab$IRR_high, use = "pairwise.complete.obs")

BF_BH = cor(df_product$IRR_product_full, df_product$IRR_product_high, use = "pairwise.complete.obs")
BF_BL = cor(df_product$IRR_product_full, df_product$IRR_product_light, use = "pairwise.complete.obs")
BL_BH = cor(df_product$IRR_product_light, df_product$IRR_product_high, use = "pairwise.complete.obs")

occu_ab = OS_occu %>% 
  left_join(df_ab, by = "Species")

OF_AF = cor(occu_ab$RR_full, occu_ab$IRR_full, use = "pairwise.complete.obs")
OF_AH = cor(occu_ab$RR_full, occu_ab$IRR_high, use = "pairwise.complete.obs")
OF_AL = cor(occu_ab$RR_full, occu_ab$IRR_light, use = "pairwise.complete.obs")

OH_AF = cor(occu_ab$RR_high, occu_ab$IRR_full, use = "pairwise.complete.obs")
OH_AH = cor(occu_ab$RR_high, occu_ab$IRR_high, use = "pairwise.complete.obs")
OH_AL = cor(occu_ab$RR_high, occu_ab$IRR_light, use = "pairwise.complete.obs")

OL_AF = cor(occu_ab$RR_light, occu_ab$IRR_full, use = "pairwise.complete.obs")
OL_AH = cor(occu_ab$RR_light, occu_ab$IRR_high, use = "pairwise.complete.obs")
OL_AL = cor(occu_ab$RR_light, occu_ab$IRR_light, use = "pairwise.complete.obs")

occu_bm = OS_occu %>% 
  left_join(df_product, by = "Species")

OF_BF = cor(occu_bm$RR_full, occu_bm$IRR_product_full, use = "pairwise.complete.obs")
OF_BH = cor(occu_bm$RR_full, occu_bm$IRR_product_high, use = "pairwise.complete.obs")
OF_BL = cor(occu_bm$RR_full, occu_bm$IRR_product_light, use = "pairwise.complete.obs")

OH_BF = cor(occu_bm$RR_high, occu_bm$IRR_product_full, use = "pairwise.complete.obs")
OH_BH = cor(occu_bm$RR_high, occu_bm$IRR_product_high, use = "pairwise.complete.obs")
OH_BL = cor(occu_bm$RR_high, occu_bm$IRR_product_light, use = "pairwise.complete.obs")

OL_BF = cor(occu_bm$RR_light, occu_bm$IRR_product_full, use = "pairwise.complete.obs")
OL_BH = cor(occu_bm$RR_light, occu_bm$IRR_product_high, use = "pairwise.complete.obs")
OL_BL = cor(occu_bm$RR_light, occu_bm$IRR_product_light, use = "pairwise.complete.obs")

ab_bm = df_product %>% 
  left_join(df_ab, by = "Species")

AF_BF = cor(ab_bm$IRR_full.y, ab_bm$IRR_product_full, use = "pairwise.complete.obs")
AF_BH = cor(ab_bm$IRR_full.y, ab_bm$IRR_product_high, use = "pairwise.complete.obs")
AF_BL = cor(ab_bm$IRR_full.y, ab_bm$IRR_product_light, use = "pairwise.complete.obs")

AH_BF = cor(ab_bm$IRR_high.y, ab_bm$IRR_product_full, use = "pairwise.complete.obs")
AH_BH = cor(ab_bm$IRR_high.y, ab_bm$IRR_product_high, use = "pairwise.complete.obs")
AH_BL = cor(ab_bm$IRR_high.y, ab_bm$IRR_product_light, use = "pairwise.complete.obs")

AL_BF = cor(ab_bm$IRR_light.y, ab_bm$IRR_product_full, use = "pairwise.complete.obs")
AL_BH = cor(ab_bm$IRR_light.y, ab_bm$IRR_product_high, use = "pairwise.complete.obs")
AL_BL = cor(ab_bm$IRR_light.y, ab_bm$IRR_product_light, use = "pairwise.complete.obs")

intra_type = c(AF_AH, AF_AL, AH_AL, OF_OH, OF_OL, OH_OL, BF_BH, BF_BL, BL_BH)
intra_prot = c(OF_AF, OF_BF, AF_BF, OH_AH, OH_BH, AH_BH, OL_AL, OL_BL, AL_BL)

cors = rbind(AF_AH, AF_AL, AH_AL, OF_OH, OF_OL, OH_OL, BF_BH, BF_BL, BL_BH,
             OF_AF, OF_BF, AF_BF, OH_AH, OH_BH, AH_BH, OL_AL, OL_BL, AL_BL,
             OF_AH, OF_AL, OH_AF, OH_AL, AF_BH, AF_BL, OF_BH, OF_BL, OH_BF, 
             OH_BL, OL_AF, OL_AH, OL_BF, OL_BH, AL_BF, AL_BH, AH_BL, AH_BF)

x_ocabm = c("Abundance", "Abundance", "Abundance", "Occurrence", "Occurrence", "Occurrence", "Biomass", "Biomass", "Biomass",
            "Occurrence", "Occurrence", "Abundance", "Occurrence", "Occurrence", "Abundance", "Occurrence", "Occurrence", "Abundance",
            "Occurrence", "Occurrence", "Occurrence", "Occurrence", "Abundance", "Abundance", "Occurrence", "Occurrence", "Occurrence",
            "Occurrence", "Occurrence", "Occurrence", "Occurrence", "Occurrence", "Abundance", "Abundance", "Abundance", "Abundance")
y_ocabm = c("Abundance", "Abundance", "Abundance", "Occurrence", "Occurrence", "Occurrence", "Biomass", "Biomass", "Biomass",
            "Abundance", "Biomass", "Biomass", "Abundance", "Biomass", "Biomass", "Abundance", "Biomass", "Biomass",
            "Abundance", "Abundance", "Abundance", "Abundance", "Biomass", "Biomass", "Biomass", "Biomass", "Biomass",
            "Biomass", "Abundance", "Abundance", "Biomass", "Biomass", "Biomass", "Biomass", "Biomass", "Biomass")

x_prot = c("Full", "Full", "High", "Full", "Full", "High", "Full", "Full", "Light", 
           "Full", "Full", "Full", "High", "High", "High", "Light", "Light", "Light", 
           "Full", "Full", "High", "High", "Full", "Full", "Full", "Full", "High", 
           "High", "Light", "Light", "Light", "Light", "Light", "Light", "High", "High")

y_prot = c("High", "Light", "Light", "High", "Light", "Light", "High", "Light", "High", 
           "Full", "Full", "Full", "High", "High", "High", "Light", "Light", "Light",
           "High", "Light", "Full", "Light", "High", "Light", "High", "Light", "Full", 
           "Light", "Full", "High", "Full", "High", "Full", "High", "Light", "Full")
cors_xy = as.data.frame(cbind(cors, x_ocabm, y_ocabm, x_prot, y_prot))
cors_xy$V1 = as.numeric(cors_xy$V1)

table(x_prot)
table(y_prot)

mean(OF_AF, OH_AH, OL_AL)
mean(OF_BF, OH_BH, OL_BL)
mean(AF_BF, AH_BH, AL_BL)

mean(OF_OH, AF_AH, BF_BH)
mean(OF_OL, AF_AL, BF_BL)
mean(OH_OL, AH_AL, BL_BH)
