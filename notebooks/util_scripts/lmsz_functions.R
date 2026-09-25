# backtransform_centiles is adapted from Apply.Param function in the brainchart/Lifespan repo
# https://github.com/brainchart/Lifespan
backtransform_centiles <- function(NEWData, FITParam,
                                   Reference.Holder=NULL,
                                   MissingToZero=TRUE, NAToZero=TRUE, Prefix="",
                                   Pred.Name="PRED", Only.New=FALSE,
                                   verbose=FALSE ) {
  
  ##
  ## check arguments
  if( is.null(NEWData) || is.null(FITParam) ) { stop("Must supply valid NEWData and FITParam objects") }
  
  if( !is.null(Reference.Holder) ) {
    NEWData <- ValidateCleanInput(IN=NEWData,
                                  Reference.Subset=Reference.Holder$SUBSET,
                                  Reference.Model=Reference.Holder$MODEL,
                                  Reference.Param=Reference.Holder$FIT.EXTRACT$param )
    ## above line checks that input conforms to reference
    ##
  }
  
  ##
  ## vvv not sure this is necessary, but a precaution against mixing with model.frame(), model.matrix(), etc.
  Saved.Attributes <- attributes(NEWData)
  attributes(NEWData)[ !names(attributes(NEWData)) %in% c("names","row.names","class") ]  <- NULL
  
  if( Only.New ) {
    stop("This isn't quite right at the moment, need to fix this to perhaps Only.Modified?")
    IN.NAMES <- names(NEWData)
  }
  
  In.Model <- attr(FITParam,"model")
  
  Model.Columns <- unlist(In.Model$covariates[c("X","BY","OTHER","RANEF")])
  
  FAMILY <- get(FITParam$family)()
  
  Model.Parameters <- names(FAMILY$parameters)
  
  Any.Missing.Columns <- Model.Columns[!c(Model.Columns %in% names(NEWData))]
  if( length(Any.Missing.Columns) > 0 ) {
    NEWData[ , Any.Missing.Columns ] <- NA
    if( verbose ) { cat("Missing columns in NEWData from within Model. Set to NA. [", Any.Missing.Columns,"]\n") }
  }
  
  for( wPARAM in Model.Parameters ) {
    Model.Formula <- as.formula(sprintf(" ~ %s",FITParam[[ wPARAM ]]$equ$fixef))
    Model.Frame <- model.frame(formula=Model.Formula, data=NEWData, na.action=na.pass )
    Model.Matrix <- model.matrix( Model.Formula, Model.Frame, contrasts.arg=FITParam$contrasts, na.action=na.pass )
    
    if( MissingToZero==TRUE ) {
      lASSIGN <- attr(Model.Matrix,"assign")
      lAllNA <- apply(Model.Matrix,2,function(X){all(is.na(X))})
      if( length(lASSIGN)!=length(lAllNA) ) {stop("FATAL ERROR: lASSIGN different length from lAllNA")}
      if( sum(lAllNA) > 0 ) {
        lLABELS <- paste( attr(terms(Model.Formula),"term.labels")[unique(lASSIGN[lAllNA])], collapse=", " )
        if( verbose ) { 
          cat(wPARAM,"Some columns of the MODEL MATRIX are all NAs (",lLABELS,"), setting to zero in MODEL MATRIX.\n")
          cat(wPARAM,"NOTE: This is mainly acceptable for factors using the contr.sum() encoding, otherwise it may give unexpected results.\n")
        }
        Model.Matrix[ , which(lAllNA) ]  <- 0
      }
    }
    if( NAToZero==TRUE ) {
      local.sum <- sum(apply(Model.Matrix,1,function(X){any(is.na(X))}))
      if( local.sum > 0 ) {
        ARR.IND <- which( is.na(Model.Matrix), arr.ind=TRUE )
        
        lCOLUMNS <- paste( attr(terms(Model.Formula),"term.labels")[unique(attr(Model.Matrix,"assign")[unique(ARR.IND[,"col"])])], collapse=", " )
        STR <- sprintf("MODEL MATRIX has %i rows (across columns: %s ) with NAs, all set to zero",local.sum,lCOLUMNS)
        if( verbose ) { 
          cat(wPARAM,STR,"\n")
          cat(wPARAM,"NOTE: This is mainly acceptable for factors using the contr.sum() encoding, otherwise it may give unexpected results.\n")
        }
        Model.Matrix[ ARR.IND ]  <- 0
      }
    }
    if( !all( colnames(Model.Matrix) %in% names(FITParam[[ wPARAM ]]$fixef) ) ) {stop("FITParam inconsistent with INPUT")}
    Fit.Fixef <- matrix( FITParam[[ wPARAM ]]$fixef[colnames(Model.Matrix)], ncol=1, dimnames=list(colnames(Model.Matrix),"Beta") )
    
    NEWData[,sprintf("%s%s.pop",Prefix,wPARAM)] <- as.vector(Model.Matrix %*% Fit.Fixef)
    if( !is.na( FITParam[[ wPARAM ]]$equ$ranef ) ) {
      ##
      ## In this block we find the matching random-effect (from FITParam) and insert it into NEWData
      ## NOTE: This code block assumes only one random-effect (per gamlss-component) [ENHANCEMENT PLAN: expand to multiple random-effects]
      ## 
      POSITION.IN.MAP <- as.numeric(NEWData[, FITParam[[ wPARAM ]]$equ$ranef ])
      MAPPING <- match( levels(NEWData[, FITParam[[ wPARAM ]]$equ$ranef ]), names(FITParam[[ wPARAM ]]$ranef) )
      NEWData[,sprintf("%s%s.ranef",Prefix,wPARAM)] <- FITParam[[ wPARAM ]]$ranef[ MAPPING[POSITION.IN.MAP] ]
      
      NEWData[,sprintf("%s%s.wre",Prefix,wPARAM)] <- NEWData[,sprintf("%s%s.pop",Prefix,wPARAM)] + NEWData[,sprintf("%s%s.ranef",Prefix,wPARAM)]
    } else {
      NEWData[,sprintf("%s%s.ranef",Prefix,wPARAM)] <- Inf ## to differentiate from missing random-effect levels (NAs), use Inf if no ranef at all
      NEWData[,sprintf("%s%s.wre",Prefix,wPARAM)] <- NEWData[,sprintf("%s%s.pop",Prefix,wPARAM)]
    }
  }
  
  ARGUMENTS <- list()
  for( lTYPE in c("pop","wre") ) {
    ARGUMENTS[[ lTYPE ]] <- list()
    for( wPARAM in Model.Parameters ) {
      LinkInvFun <- FAMILY[[sprintf(sprintf("%s.linkinv",wPARAM))]]
      ComponentValues <- NEWData[,sprintf("%s%s.%s",Prefix,wPARAM,lTYPE)]
      
      ARGUMENTS[[ lTYPE ]][[ wPARAM ]] <- LinkInvFun( ComponentValues )
    }
  }
  
  ##
  ## Select rows (within arguments) that are valid, ie. not NAs
  KEEP <- lapply( ARGUMENTS, function(Y){ as.vector( Reduce(f=`&`,x=lapply(Y,function(X){!is.na(X)})), mode="logical") } )
  SHORT <- setNames(lapply(1:length(ARGUMENTS),function(IDX){ lapply(ARGUMENTS[[IDX]],function(X){ as.vector(X[ KEEP[[IDX]] ]) } ) } ), names(KEEP) )
  WHICH <- lapply( KEEP, which )
  
  ##
  ## in below we make versions that account for NAs in outcome (OUT.NAME) column
  OUT.NAME <- In.Model$covariates$Y
  if( (OUT.NAME %in% names(NEWData)) && any(!is.na(NEWData[,OUT.NAME])) ) {
    out.KEEP.ARGS <- lapply( ARGUMENTS, function(Y){ as.vector( Reduce(f=`&`,x=lapply(Y,function(X){!is.na(X)})), mode="logical") } )
    out.KEEP <- lapply( out.KEEP.ARGS, function(X) { (X) & (!is.na(NEWData[,OUT.NAME])) } )
    out.SHORT <- setNames(lapply(1:length(ARGUMENTS),function(IDX){ lapply(ARGUMENTS[[IDX]],function(X){ as.vector(X[ out.KEEP[[IDX]] ]) } ) } ), names(out.KEEP) )
    out.WHICH <- lapply( out.KEEP, which )
    
    KEEP.C <- Reduce(`&`,out.KEEP)
    SHORT.C <- lapply( ARGUMENTS, function(Y){ lapply(Y,function(X){ X[KEEP.C] } ) } )
    WHICH.C <- which( KEEP.C )
  }
  
  for( lTYPE in c("pop","wre") ) {
    if( length(WHICH[[lTYPE]])>0 ) {
      
      NEWData[WHICH[[lTYPE]],sprintf("%s%s.%s",Prefix,Pred.Name,lTYPE)] <- do.call(what=get(paste0("q",FAMILY$family[1])),
                                                                                   args=c(SHORT[[lTYPE]],list(p=NEWData$value)))
      if( (OUT.NAME %in% names(NEWData)) && any(!is.na(NEWData[,OUT.NAME])) ) {
        
        NEWData[out.WHICH[[lTYPE]],sprintf("%s%s.q.%s",Prefix,OUT.NAME,lTYPE)] <- do.call(what=get(paste0("p",FAMILY$family[1])),
                                                                                          args=c(out.SHORT[[lTYPE]],
                                                                                                 list(q=NEWData[out.WHICH[[lTYPE]],OUT.NAME,drop=TRUE])))
      }
    }
  }
  
  if( Only.New ) {
    NEWData <- NEWData[,names(NEWData)[!names(NEWData)%in%IN.NAMES]]
  }
  attributes(NEWData) <- c( attributes(NEWData), Saved.Attributes[!names(Saved.Attributes) %in% c("names","row.names","class")] )
  
  return( NEWData )
}


backtransform_centile_fans <- function(pheno, sim_list) {
  # Steps:
  # 1 - Load SLIP and LMSz growth charts
  # 2 - Generate fan lines for LMSz growth charts. These will give the centile 
  #     lines in brain z-score space (i.e., the starting point of the LMSz models)
  # 3 - Reformat the z-score space centile lines
  # 4 - Calculate the original SLIP centile fans just to have
  # 5 - Feed z-score space centile fans to backtransform_centiles. This function
  #     acts similar to Apply.Param in that it backtransforms a centile to the original
  #     measure, but while Apply.Param applies one or a selected number of centile fans
  #     across all ages, backtransform_centiles operates row-by-row, allowing for 
  #     each row (combination of age, sex, and centile fan) to have a different value
  # 6 - Combine all of the results
  # Step 1) Load SLIP and LMSz growth charts
  # Load data
  pheno_name <- glue("{pheno}Transformed.q.wre")
  norm_name <- glue("{pheno}Transformed.normalised")
  zscore_name <- glue("{pheno}Transformed.q.zscore")
  
  old_pheno <- slip_features$old_pheno[slip_features$new_pheno == pheno]
  if (length(old_pheno) == 0) {
    return()
  }
  growthChartModel <- readRDS(glue("../lmsz/lmsz_models_slip/SLIP_lmsz_{old_pheno}Transformed.rds"))
  # Load SLIP fit model
  if (grepl("CT",pheno)) {
    fit_file <- glue("../../../../braincharts/SLIP/recon-all-clinical/gamlss/RDS/mpr-{old_pheno}Transformed/FIT.EXTRACT.rds")
  } else {
    fit_file <- glue("../../../../braincharts/SLIP/recon-all-clinical/gamlss/RDS/rac-{old_pheno}Transformed/FIT.EXTRACT.rds")
  }
  
  orig_fit <- readRDS(fit_file)
  
  df_tmp <- df %>%
    select(participant, AgeTransformed, sex, !!pheno_name, !!norm_name) %>%
    rename(pheno = !!pheno_name,
           normalised_value = !!norm_name) %>%
    mutate(pheno = !!pheno_name,
           zscore = !!zscore_name) %>%
    na.omit()
  
  # Step 2) Generate fan lines for LMSz growth charts.
  # Simulate centiles
  lmsz_sim <- centile_predict(gamlssModel = growthChartModel, 
                              sim_df_list = sim_list, 
                              x_var = "AgeTransformed", 
                              desiredCentiles = c(0.025,0.250,0.50,0.750,0.975),
                              df = df_tmp,
                              average_over = FALSE,
                              resid_terms = FALSE)
  
  # Step 3) Reformatting
  # Reformat simulated values
  lmsz_sim$fanCentiles_Female <- lmsz_sim$fanCentiles_Female %>%
    pivot_longer(-c(AgeTransformed),
                 names_to = "centile_cat",
                 values_to = "value") %>%
    mutate(sex = "Female")
  
  lmsz_sim$fanCentiles_Male <- lmsz_sim$fanCentiles_Male %>%
    pivot_longer(-c(AgeTransformed),
                 names_to = "centile_cat",
                 values_to = "value") %>%
    mutate(sex = "Male")
  
  lmsz_df <- bind_rows(lmsz_sim$fanCentiles_Female,
                       lmsz_sim$fanCentiles_Male) %>%
    mutate(site = NA,
           sex = factor(sex, levels = c("Female","Male")),
           LMSz_input_fans = value,
           LMSz_input_centile = pnorm(value)) %>%
    mutate(line_size = case_match(centile_cat,
                                  c("cent_0.025","cent_0.975") ~ 0.3,
                                  c("cent_0.25","cent_0.75") ~ 0.3,
                                  c("cent_0.5") ~ 0.4),
           line_type = case_match(centile_cat,
                                  c("cent_0.025","cent_0.975") ~ "dotted",
                                  c("cent_0.25","cent_0.75") ~ "longdash",
                                  c("cent_0.5") ~ "solid"),
           .after = sex) %>%
    select(-c(value))
  # Step 4) Calculate the original SLIP centile fans just to have
  # Plot original SLIP Model
  NEWData <- lmsz_df %>%
    select(AgeTransformed, sex) %>% 
    distinct() %>%
    as.data.frame()
  
  original_fans <- Apply.Param(NEWData, orig_fit$param,
                               Pred.Set=c("l025"=0.025,"l250"=0.250,
                                          "m500"=0.5,"u750"=0.750,
                                          "u975"=0.975)) %>%
    pivot_longer(c(PRED.l025.pop,PRED.l250.pop,
                   PRED.m500.pop,PRED.u750.pop,
                   PRED.u975.pop),
                 names_to = "centile_cat", values_to = "SLIP_input_fans") %>%
    mutate(centile_cat = case_match(centile_cat,
                                    "PRED.l025.pop" ~ "cent_0.025",
                                    "PRED.l250.pop" ~ "cent_0.25",
                                    "PRED.m500.pop" ~ "cent_0.5",
                                    "PRED.u750.pop" ~ "cent_0.75",
                                    "PRED.u975.pop" ~ "cent_0.975")) %>%
    select(AgeTransformed, sex, centile_cat, SLIP_input_fans)
  
  # Step 5) Backtransform LMSz centile fans further into original measures
  # Plot backtransformed centiles
  NEWData <- lmsz_df %>%
    select(AgeTransformed, sex,
           centile_cat,
           value = LMSz_input_centile) %>% 
    distinct() %>%
    as.data.frame()
  
  backtransformed_fans <- backtransform_centiles(NEWData, orig_fit$param)
  
  backtransformed_fans <- backtransformed_fans %>%
    select(AgeTransformed, sex, centile_cat, 
           LMSz_fans = PRED.pop)
  
  # Step 6) Combine results
  # Join all results
  out_df <- lmsz_df %>%
    inner_join(original_fans, by= join_by(AgeTransformed, centile_cat, sex)) %>%
    inner_join(backtransformed_fans, by= join_by(AgeTransformed, centile_cat, sex)) %>%
    mutate(pheno = pheno, model = "SLIP", .before = AgeTransformed) %>%
    mutate(fit_file = fit_file)
  
  return(out_df)
}

# pred_og_centile_mod is adapted from the gamlssTools package
# https://github.com/BGDlab/gamlssTools

pred_og_centile_mod <- function(gamlssModel, og.data, get.std.scores = FALSE, new.data=NULL){
  pheno <- gamlssModel$mu.terms[[2]]
  
  #subset df cols just to predictors from model
  predictor_list <- list_predictors(gamlssModel)
  stopifnot("Dataframe columns and model covariates don't match" = 
              predictor_list %in% names(og.data))
  if (is.null(new.data)) {
    # newData <- subset(og.data, select = predictor_list)
    newData <- og.data
    predict_me <- og.data
  } else {
    stopifnot("Dataframe columns and model covariates don't match" = 
                predictor_list %in% names(new.data))
    #newData <- subset(new.data, select = predictor_list)
    newData <- new.data
    predict_me <- new.data
    #make sure all vals are within range of those originally modeled
    check_range(subset(og.data, select = predictor_list), newData)
  }
  
  #predict
  predModel <- predictAll(gamlssModel, newdata=newData, data=og.data, type= "response")
  
  #get dist type (e.g. GG, BCCG) and write out function
  fname <- gamlssModel$family[1]
  pfun <- paste0("p", fname)
  
  #look for moments
  has_sigma <- "sigma" %in% gamlssModel[[2]]
  has_nu <- "nu" %in% gamlssModel[[2]]
  has_tau <- "tau" %in% gamlssModel[[2]]
  
  centiles <- c()
  #iterate through participants
  for (i in 1:nrow(predict_me)){
    cent_args <- list(predict_me[[pheno]][[i]], predModel$mu[[i]])
    
    if (has_sigma){
      cent_args$sigma <- predModel$sigma[[i]]
    }
    if (has_nu){
      cent_args$nu <- predModel$nu[[i]]
    } 
    if (has_tau){
      cent_args$tau <- predModel$tau[[i]]
    } 
    
    centiles[i] <- do.call(pfun, cent_args)
    
    #don't let centile = 1 (for z-scores)!
    if (centiles[i] == 1) {
      centiles[i] <- 0.9999999999999999
    }
    #don't let centile = 0 (for z-scores)!
    if (centiles[i] == 0) {
      centiles[i] <- 0.0000000000000001
    }
    
  }
  if (get.std.scores == FALSE){
    return(centiles)
  } else {
    #get 'z scores' from normed centiles - how z.score() does it
    rqres <- qnorm(centiles)
    
    #return dataframe
    df <- data.frame("centile" = centiles,
                     "std_score" = rqres)
    return(df)
  } 
  
}

