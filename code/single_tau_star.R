library(data.table)
library(rmeta)

annot_cell = "/n/groups/price/kushal/CorShrink2/data/ANNOTATIONS"
results_cell = "/n/groups/price/kushal/CorShrink2/data/LDSC_RESULTS"
annot_names_1 = as.character(read.table("/n/groups/price/kushal/CorShrink2/data/annot.txt")[,1])
annot_names = annot_names_1
#annot_names = annot_names_1[c(32, 33)]
annot_idx = 1



all_traits = c('UKB_460K.body_BMIz','UKB_460K.cov_EDU_YEARS','UKB_460K.lung_FVCzSMOKE','UKB_460K.cov_SMOKING_STATUS',
                 'UKB_460K.mental_NEUROTICISM','UKB_460K.blood_WHITE_COUNT','PASS_Years_of_Education2','UKB_460K.bp_SYSTOLICadjMEDz',
                 'UKB_460K.body_HEIGHTz','UKB_460K.other_MORNINGPERSON','UKB_460K.body_WHRadjBMIz','UKB_460K.lung_FEV1FVCzSMOKE',
                 'UKB_460K.repro_MENARCHE_AGE','UKB_460K.blood_RED_COUNT','UKB_460K.blood_PLATELET_COUNT','UKB_460K.bmd_HEEL_TSCOREz',
                 'UKB_460K.blood_EOSINOPHIL_COUNT','PASS_Schizophrenia','UKB_460K.blood_RBC_DISTRIB_WIDTH','PASS_Height1','PASS_BMI1',
                 'UKB_460K.disease_T2D','PASS_AgeFirstBirth','UKB_460K.disease_RESPIRATORY_ENT','UKB_460K.body_BALDING1','UKB_460K.disease_HYPOTHYROIDISM_SELF_REP',
                 'UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED','UKB_460K.disease_HI_CHOL_SELF_REP','UKB_460K.repro_MENOPAUSE_AGE','PASS_HDL','UKB_460K.pigment_SUNBURN',
                 'PASS_NumberChildrenEverBorn','PASS_Anorexia','PASS_LDL','PASS_Crohns_Disease','PASS_DS','PASS_Ever_Smoked','UKB_460K.pigment_HAIR',
                 'PASS_Rheumatoid_Arthritis','PASS_Type_2_Diabetes','PASS_Autism','UKB_460K.pigment_TANNING','PASS_Ulcerative_Colitis',
                 'UKB_460K.disease_DERMATOLOGY','PASS_Coronary_Artery_Disease','UKB_460K.disease_AID_SURE','UKB_460K.pigment_SKIN')

blood_traits = c("UKB_460K.blood_RBC_DISTRIB_WIDTH", "UKB_460K.blood_RED_COUNT", "UKB_460K.blood_WHITE_COUNT",
                 "UKB_460K.blood_PLATELET_COUNT", "UKB_460K.blood_EOSINOPHIL_COUNT")

autoimmune_traits = c("UKB_460K.disease_AID_SURE", "PASS_Ulcerative_Colitis", "PASS_Crohns_Disease", "PASS_Rheumatoid_Arthritis",
                      "PASS_Type_2_Diabetes")

brain_traits = c("PASS_Ever_Smoked", "UKB_460K.cov_SMOKING_STATUS", "UKB_460K.mental_NEUROTICISM", "UKB_460K.repro_MENARCHE_AGE",
                  "PASS_Years_of_Education2", "PASS_DS", "PASS_Schizophrenia", "UKB_460K.body_WHRadjBMIz",
                 "PASS_BMI1", "UKB_460K.body_BMIz")



get_sd_annot = function(cell_path, annot_index = 1, flag=0){
  if(flag == 0){
    if(file.exists(paste0(cell_path, "/", "sd_annot_", annot_index, ".rda"))){
        sd_annot = get(load(paste0(cell_path, "/", "sd_annot_", annot_index, ".rda")))
        return(sd_annot)
    }else{
    	flag = 1
    }}
        
    if(flag == 1){
    	num = 0
        den = 0
        ll <- list.files(cell_path, pattern = ".annot.gz")
        for(m in 1:length(ll)){
            dat <- data.frame(fread(paste0("zcat ", cell_path, "/", ll[m])))
            num = num  + (nrow(dat)-1) * var(dat[,4+annot_index])
            den = den + (nrow(dat)-1)
            rm(dat)
       }   
    }
    
  estd_sd_annot = sqrt(num/den)
  save(estd_sd_annot, file = paste0(cell_path, "/", "sd_annot_", annot_index, ".rda"))
  return(estd_sd_annot)
}

run_single_tau_analysis = function(annot_cell,
                                   results_cell,
                                   annotations,
                                   traits,
                                   index_in_results=1,
                                   base_index = NULL,
                                   flag = 1){
    if(is.null(base_index)){base_index = index_in_results}
    tau_star_table = matrix(0, length(annotations), 3)
    for(annot_id in 1:length(annotations)){
        cell_path = paste0(annot_cell, "/", annotations[annot_id])
        sd_annot1=get_sd_annot(cell_path, annot_index=index_in_results, flag = flag)
        Mref = 5961159
        df = c()
        for(trait_id in 1:length(traits)){
            result.file=paste0(results_cell, "/", annotations[annot_id], "/", traits[trait_id], ".sumstats.part_delete")
            new_table=read.table(result.file,header=F)
            sc=c()
            logfile = paste(results_cell, "/", annotations[annot_id], "/", traits[trait_id],".sumstats.log", sep="")
            log = read.table(logfile,h=F,fill=T)
            h2g = as.numeric(as.character(log[which(log$V4=="h2:"),5]))
            coef1=sd_annot1*Mref/h2g
            for(i in 1:dim(new_table)[1]){
                  tau1=as.numeric(new_table[i,base_index])
                  taus1=tau1*coef1
                  sc=c(sc,taus1)
                  #cat("Block ", i, "\n")
            }
            mean_sc=mean(sc)
            se_sc=sqrt(199**2/200*var(sc))
            df = rbind(df, c(mean_sc,se_sc))
        }
        test_tauj=meta.summaries(df[,1],df[,2],method="random")
        tau=test_tauj$summary
        tau_se=test_tauj$se.summary
        z=tau/tau_se
        cat("Printing results for annotation:", annotations[annot_id], "\n")
        cat(tau, " ", tau_se, " ", 2*pnorm(-abs(z)), "\n")
        tau_star_table[annot_id, ] = c(tau, tau_se, 2*pnorm(-abs(z)))
    }
    rownames(tau_star_table) = annotations
    return(tau_star_table)
}





out1 = run_single_tau_analysis(annot_cell, results_cell, annotations = annot_names, traits = all_traits, 
                                index_in_results = annot_idx, flag = 0)
out2 = run_single_tau_analysis(annot_cell, results_cell, annotations = annot_names, traits = brain_traits, 
                                index_in_results = annot_idx, flag = 0)
out3 = run_single_tau_analysis(annot_cell, results_cell, annotations = annot_names, traits = c(blood_traits, autoimmune_traits), 
                                index_in_results = annot_idx, flag = 0)

ll <- list()
ll[["All"]] = out1
ll[["Brain"]] = out2
ll[["Blood"]] = out3


save(ll, file = "/n/groups/price/kushal/EXPECTOCPP-PROJECT/Figures/CNN_SELCTION_PREDICT.rda")







