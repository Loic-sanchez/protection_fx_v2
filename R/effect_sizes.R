compute_ES_occu = function(x) {
  
  melted_cov$Effectiveness = relevel(as.factor(melted_cov$Effectiveness), ref = "Unprotected")
  
  all_sp_glm_qual = mclapply(1:length(unique(melted_cov$Species)), mc.cores = x, function(i) {

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
    
    pca = ade4::dudi.pca(pca_prep[, -c(1, 97)],
                         center = T, scale = T, 
                         row.w = pca_prep$weights,
                         nf = 25, 
                         scannf = F)
    
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
    
    full_formula_sp = as.formula(paste(full_formula, "+ Matern(1 | SiteLongitude + SiteLatitude)", sep = ""))
    separ = spaMM::is_separated.formula(formula = full_formula_sp,
                                        data = presabs_sp1)
    
    if(separ == T){return(NULL)}
    
    sp_model = spaMM::fitme(full_formula_sp,
                            family = "binomial",
                            method = "PQL/L",
                            data = presabs_sp1)
    
    mod_coefs = as.data.frame(sp_model$fixef)
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
    
    RR_full = estimate_full/estimate_outside
    RR_high = estimate_high/estimate_outside
    RR_light = estimate_light/estimate_outside
    
    p = spaMM::predict.HLfit(sp_model, type = "response")
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

compute_ES_abun = function(x) {
  
melted_cov_ab$Effectiveness = relevel(as.factor(melted_cov_ab$Effectiveness), ref = "Unprotected")

all_sp_glm_ab = parallel::mclapply(1:length(unique(melted_cov_ab$Species)), mc.cores = x, function(i) {

  subdf = filter(melted_cov_ab, Species == unique(melted_cov_ab$Species)[i])
  
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
  
  eco_weights = presabs_sp1 %>%
    group_by(Ecoregion) %>% 
    mutate(weights = (1/n())) %>% 
    ungroup()
  
  pca_prep = all_cov[all_cov$SurveyID %in% unique(presabs_sp1$SurveyID),]
  pca_prep = pca_prep %>% 
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
  pres_sp1 = presabs_sp1[, 1:10] %>% 
    left_join(pca_prep_and_coords[, c(1, 98:122)], by = "SurveyID") 
  
  istherefull = ifelse(length(pres_sp1$Abundance[which(pres_sp1$Effectiveness == "Fully Protected" & pres_sp1$Abundance > 0)]) > 9, T, F)
  istherehigh = ifelse(length(pres_sp1$Abundance[which(pres_sp1$Effectiveness == "Highly Protected" & pres_sp1$Abundance > 0)]) > 9, T, F)
  istherelight = ifelse(length(pres_sp1$Abundance[which(pres_sp1$Effectiveness == "Lightly Protected" & pres_sp1$Abundance > 0)]) > 9, T, F)
  isthereunp = ifelse(length(pres_sp1$Abundance[which(pres_sp1$Effectiveness == "Unprotected" & pres_sp1$Abundance > 0)]) > 9, T, F)
  
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
  tryCatch({RSQ = cor(p, pres_sp1$Abundance)^2}, error = function(e) RSQ <<- NA)
  
  Species = unique(melted_cov_ab$Species)[i]
  
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

return(df_ab)

}

compute_ES_bm = function(x) {
  
  melted_cov_bm$Effectiveness = relevel(as.factor(melted_cov_bm$Effectiveness), ref = "Unprotected")
  
  all_sp_glm_bm = mclapply(1:length(unique(melted_cov_bm$Species)), mc.cores = x, function(i) {

    subdf = filter(melted_cov_bm, Species == unique(melted_cov_bm$Species)[i])
    
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
    
    if(is.null(sp_model) == F){mod_coefs = as.data.frame(sp_model$fixef)}
    if(is.null(sp_model) == F){intercept = mod_coefs["(Intercept)",]
    
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
    tryCatch({RSQ = (cor(10^(p), pres_sp1$trunc_bm))^2}, error = function(e) RSQ <<- NA)
    
    Species = unique(melted_cov_bm$Species)[i]
    
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
    
  }
  
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
  
})}
