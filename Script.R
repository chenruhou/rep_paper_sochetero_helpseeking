##############
# Replication script of 
# ”Is Social Heterogeneity in Classrooms Associated with Reduced Achievement 
# Inequality? The Role of Help-Seeking in Peer Networks"
# Authors: Hou, Chenru, Georg Lorenz, Camilla Rjosk
# Data set: Stanat, P., Schipolowski, S., Mahler, N., Weirich, S., Henschel, S., Holtmann, M., Becker, B. & Kölm, J. (2022). 
#           IQB-Bildungstrend Mathematik und Naturwissenschaften 2018 (IQB-BT 2018) (Version 2) [Datensatz]. 
#           Berlin: IQB – Institut zur Qualitätsentwicklung im Bildungswesen. http://doi.org/10.5159/IQB_BT_2018_v2
#############


library(dplyr)
library(psych)
library(lmerTest)
library(sjPlot)
library(tidyr)
library(datawizard)
library(ggplot2)


# 1. Data cleaning -------------
## sample restriction before imputation ----------
data_res<- subset(dat_com,dat_com$TR_SCHULFORM!=8) #special education school and school without type information 

data_res<- subset(data_res, data_res$TR_T_SFB !=0 &
                    data_res$TR_T_SFB !=8 &
                    data_res$TR_T_SFB !=9 &
                    data_res$SFB_bearbeitet !=0) #remove individuals do not have any information 
  ## TR_T_SFB: 0= abwesend/Ausfall, 8= Elterngenehmigung fehlt, 9=Schule/Klasse verlassen
  ## SFB_bearbeitet: 0 = nein

# 2. Multiple Imputation --------------
data_res_select <-  select(data_res, IDTESTGROUP_FDZ, IDSTUD_FDZ,HISEI, HISCED, hEGP_6kat, 
                           herkunft_FDZ,
                           TR_GESCHLECHT,
                           befki_wle_1, TR_NOTE_MAT_r,
                           TR_SCHULFORM)


data_res_select$HISCED <- as.factor(data_res_select$HISCED)
data_res_select$hEGP_6kat <- as.factor(data_res_select$hEGP_6kat)
data_res_select$TR_NOTE_MAT_r <- as.factor(data_res_select$TR_NOTE_MAT_r)
data_res_select$herkunft_FDZ <- as.factor(data_res_select$herkunft_FDZ)
data_res_select$TR_SCHULFORM <- as.factor(data_res_select$TR_SCHULFORM)
data_res_select$TR_GESCHLECHT <- as.factor(data_res_select$TR_GESCHLECHT)

  missing_counts <- colSums(is.na(data_res_select))
  max_na_column <- names(which.max(missing_counts))
  max_na_count <- missing_counts[max_na_column]
  cat("The column with the most missing values is:", max_na_column,
      "with", max_na_count, "missing values.\n")
  
  
  library(mice)
  analytical_vars <- c(data_res_select[, c("HISEI", "HISCED", "hEGP_6kat", 
                                           "herkunft_FDZ",
                                           "TR_GESCHLECHT",
                                           "befki_wle_1", "TR_NOTE_MAT_r",
                                           "TR_SCHULFORM")])
  pred_matrix <- quickpred(data_res_select, include = analytical_vars, mincor = 0.2, method = "pearson")
  impMethod <- make.method(data_res_select)
  impMethod[]<- ""
  
  impMethod [c("HISCED","hEGP_6kat", "TR_NOTE_MAT_r") ] <- "polr"
  impMethod [c("TR_SCHULFORM") ] <- "polyreg"
  impMethod [c("TR_GESCHLECHT") ] <- "logreg"
  impMethod [c("herkunft_FDZ","HISEI","befki_wle_1")] <- "pmm"
  
  exclude_vars <- c("IDTESTGROUP_FDZ", "IDSTUD_FDZ")
  pred_matrix[, exclude_vars] <- 0
  
  diag(pred_matrix) <- 0
  
  imp<- mice(
    data_res_select,
    m=15,
    maxit = 1,
    method=impMethod,
    predictorMatrix=pred_matrix,
    remove_lindep=F,
    seed=1234
  )
  
  imp_res <- mice::complete(imp,action=1,include=FALSE)
  sum(is.na(imp_res))
  

# 3. Sample restriction ---------------
  ## Without network data ---------
help_cols <- paste0("help_", 1:34)
  
only_na_network <- network_all[apply(network_all[, help_cols], 1, function(x) all(is.na(x))), ]
  only_na_network_id <- select(only_na_network, IDSTUD_FDZ)
  only_na_network_id$IDSTUD_FDZ <- as.character(only_na_network_id$IDSTUD_FDZ)

    imp_res2 <- anti_join(imp_res, only_na_network_id)  
  
isolate <- subset(network_all, network_all$total_tie ==0)
  isolate_id$IDSTUD_FDZ <- as.character(isolate_id$IDSTUD_FDZ)
  isolate_id <- select(isolate, IDSTUD_FDZ)
  
    imp_res2 <- anti_join(imp_res2, isolate_id)
    
  ## Classrooms with less than 10 students after restriction -------------
imp_res3 <- imp_res2%>%
      group_by(IDTESTGROUP_FDZ)%>%
      mutate(size_ana=n())%>%
      ungroup()
    
  imp_res3 <- subset(imp_res3, size_ana>=10) 
    
data_dc <- imp_res3
data_dc$IDTESTGROUP_FDZ <- as.character(data_dc$IDTESTGROUP_FDZ )
data_dc$IDSTUD_FDZ <- as.character(data_dc$IDSTUD_FDZ)
    
  
# 4. Variables generation after imputation---------
  ## Classroom heterogeneity with 3-category HISEI ------------
data_dc$ses_3c <-gtools::quantcut(data_dc$HISEI, 3)
  
  #define function: calculating Blau's index  
simps <- function(df, column_name) { 
  library(vegan)
  simpdat <- df[which(!(is.na(df[[column_name]]))), ]
  simpdat1 <- simpdat %>% group_by(IDTESTGROUP_FDZ) %>% count(!!sym(column_name))
  simpdat2 <- simpdat1 %>% tidyr::pivot_wider(names_from = IDTESTGROUP_FDZ, values_from = n)
  simpdat2 <- arrange(simpdat2, !!sym(column_name))
  simpdat3 <- data.frame(t(simpdat2[-1]))
  
  name <- rownames(simpdat2) 
  colnames(simpdat3) <- name
  simpdat3[is.na(simpdat3)] <- 0
  
  simps <- as.data.frame(diversity(simpdat3, index = "simpson"))
  
  names <- rownames(simps)
  rownames(simps) <- NULL
  simps<- cbind(names, simps)
  return(simps)
}

simps.3c <- simps(data_dc, "ses_3c")
  simps.3c<- rename(simps.3c, simps.3c=`diversity(simpdat3, index = "simpson")`)
  simps.3c <- rename(simps.3c, IDTESTGROUP_FDZ=names)
  simps.3c $IDTESTGROUP_FDZ <- as.character(simps.3c$IDTESTGROUP_FDZ)


  ## Proportion of cross-social-origin help-seeking -----------
# load network data with all individuals 
network_all <- select(dat_com, IDTESTGROUP_FDZ, IDinClass,IDSTUD_FDZ,help,HISEI)
  network_all$total_tie <- rowSums(network_all[, 4:37], na.rm = TRUE)

network_all$ses_3c_dc <- cut(network_all$HISEI,
                             breaks = c(11, 42.3, 61.8, 89),
                             right = TRUE,         # intervals are right-closed (i.e., include the upper limit)
                             include.lowest = TRUE,  # ensures the first interval includes the lowest value (11)
                             labels = c("[11,42.3]", "(42.3,61.8]", "(61.8,89]"))

# separate dataframe into classrooms
network_all_cl<- split(network_all, network_all$IDTESTGROUP_FDZ)
network_all_cl<- lapply(network_all_cl, function(df) {
  df[order(df$IDinClass), ]
}) #change the order of individuals within classrooms

# define function: calculate the proportion 
ratio_ses_3c_dc <- function(group) {
  individual_results <- list()
  help_columns <- paste0("help_", 1:34)
  
  for (i in 1:nrow(group)) {
    same_ties <- 0
    inter_ties <- 0
    total_ties <- group$total_tie[i]
    
    for (help_col in help_columns) {
      adjacency_value <- group[i, help_col]
      
      row_index <- which(colnames(group) == help_col) - 3  
      
      if (!is.na(adjacency_value) && adjacency_value == 1 && 
          !is.na(group$ses_3c_dc[i]) && !is.na(group$ses_3c_dc[row_index])) {
        

      if (group$ses_3c_dc[i] == group$ses_3c_dc[row_index]) {
          same_ties <- same_ties + 1
        }
      } 
    }
    
    # Calculate inter_ties
    inter_ties <- total_ties - same_ties
    
    # Calculate proportion, 
    ratio <- if (total_ties > 0) {
      if (is.na(inter_ties / total_ties)) {
        NA
      } else {
        inter_ties / total_ties
      }
    } else {
      NA 
    }
    
    individual_results[[i]] <- list(IDTESTGROUP_FDZ = group$IDTESTGROUP_FDZ[i], 
                                    IDSTUD_FDZ = group$IDSTUD_FDZ[i], 
                                    ratio_3c_dc = ratio, 
                                    inter_tie_3c_dc = inter_ties, 
                                    same_tie_3c_dc = same_ties,
                                    total_ties = total_ties)
  }
  
  return(individual_results)
}


ratio_all_dc<- lapply(network_all_cl, ratio_ses_3c_dc)
  ratio_all_dc2<- lapply(ratio_all_dc, bind_rows)%>%
  bind_rows(.)
  ratio_all_dc2$IDSTUD_FDZ <- as.character(ratio_all_dc2$IDSTUD_FDZ)
  ratio_all_dc2$IDTESTGROUP_FDZ <- as.character(ratio_all_dc2$IDTESTGROUP_FDZ)

  
  
data_dc<- left_join(data_dc, ratio_all_dc2, by=c("IDSTUD_FDZ", "IDTESTGROUP_FDZ"))

  ## Achievement --------------
replacements <- c("1"=6, "2"=5, "3"=4, "4"=3, "5"=2, "6"=1)
data_dc$TR_NOTE_MAT_recode <- replacements[as.factor(data_dc$TR_NOTE_MAT_r)]


  ## Consolidation of social and ethnic origins ----------

data_conso_eth_dc <- dplyr::select(data_dc, IDTESTGROUP_FDZ,herkunft_FDZ, ses_3c)%>%as.data.frame()
  data_conso_eth_dc$herkunft_FDZ <- as.factor(data_conso_eth_dc$herkunft_FDZ)

  library(asw.cluster)
shaw_eth_dc <- faultlines(data_conso_eth_dc, 
                          group.par = c("IDTESTGROUP_FDZ"), 
                          attr.type = c("nominal","nominal"), 
                          rescale = "sd",
                          method = "shaw",
                          usesghomo = FALSE)

  shaw_eth_imp_long_dc <- summary(shaw_eth_dc)$long
  shaw_eth_imp_indi_dc <- select(shaw_eth_imp_long_dc, c(team,fl.value))
  shaw_eth_imp_indi_dc <- rename(shaw_eth_imp_indi_dc, IDTESTGROUP_FDZ=team, shaw_eth_dc=fl.value)

shaw_eth_imp_cl_dc<- shaw_eth_imp_indi_dc%>%
  group_by(IDTESTGROUP_FDZ) %>%
  summarise(shaw_eth_dc = unique(shaw_eth_dc))

  data_dc <- left_join(data_dc, shaw_eth_imp_cl_dc)

  ## Consolidation of social and sex ----------

data_conso_gen_dc <- dplyr::select(data_dc, IDTESTGROUP_FDZ,TR_GESCHLECHT, ses_3c)%>%as.data.frame()
data_conso_gen_dc$TR_GESCHLECHT <- as.factor(data_conso_gen_dc$TR_GESCHLECHT)

shaw_gen_dc <- asw.cluster::faultlines(data_conso_gen_dc, 
                                       group.par = c("IDTESTGROUP_FDZ"), 
                                       attr.type = c("nominal","nominal"), 
                                       rescale = "sd",
                                       method = "shaw",
                                       usesghomo = FALSE)


  shaw_gen_imp_long_dc <- summary(shaw_gen_dc)$long
  shaw_gen_imp_indi_dc <- select(shaw_gen_imp_long_dc, c(team,fl.value))
  shaw_gen_imp_indi_dc <- rename(shaw_gen_imp_indi_dc, IDTESTGROUP_FDZ=team, shaw_gen_dc=fl.value)


shaw_gen_imp_cl_dc <- shaw_gen_imp_indi_dc%>%
  group_by(IDTESTGROUP_FDZ) %>%
  summarise(shaw_gen_dc= unique(shaw_gen_dc))

data_dc <- left_join(data_dc, shaw_gen_imp_cl_dc)


  ## Classroom level control variables -----------

# Average social origin within classrooms
data_dc <- data_dc%>%
  group_by(IDTESTGROUP_FDZ) %>%
  mutate(mean_ses = mean(HISEI))%>%
  ungroup()

# Classroom size (without restriction)
network_all<- network_all%>%
  group_by(IDTESTGROUP_FDZ)%>%
  mutate(size_orig=n())%>%
  ungroup()
size <- select(network_all, IDSTUD_FDZ,size_orig)
  size$IDSTUD_FDZ <- as.character(size$IDSTUD_FDZ )

  data_dc  <- left_join(data_dc,size)  
# Type of school
data_dc<- data_dc%>%
  mutate(gym=ifelse(TR_SCHULFORM==5,1,0))
  
    data_dc$gym <- as.factor(data_dc$gym)
    
  ## Individual level control variables ----------

# Immigrant background 
data_dc<- data_dc%>%
  mutate(immi=ifelse(herkunft_FDZ==1, 0, 1))
    
    data_dc$immi <- as.factor(data_dc$immi)  
    
# Number of same-social-origin peers within classroom
data_dc<- data_dc%>%
  group_by(IDTESTGROUP_FDZ, ses_3c) %>%
  mutate(same_hisei_3c = n()-1) %>%
  ungroup()


# 5. Descriptive statistics ------------

  ## figure ---------------
    ### Panel a: social origin inequality  ------------

graph1<- data_dc %>%
  group_by(ses_3c,TR_NOTE_MAT_recode) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(ses_3c) %>%
  mutate(total = sum(count)) %>%
  mutate(perc = count / total * 100)

graph1$ses_3c <- factor(graph1$ses_3c,
                        levels = c("[11,42.3]", "(42.3,61.8]", "(61.8,89]"),
                        labels = c("Low", "Intermediate", "High"))

library(ggplot2)
ineq_graph_3cat<- ggplot(graph1 , aes(x = ses_3c, y = perc, fill = as.factor(TR_NOTE_MAT_recode))) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Blues") + # Use a blue palette
  theme_minimal() +
  labs(tag = "a", x= "Social origin group", y = "Percentage", fill = "Grade Category") +
  guides(fill = guide_legend(title = "Reversed \n Grade", size=5))


      ### Panel b: proportion of cross-social-origin help-seeking ties ---------
ratio_graph_3cat <- select(data_dc, ratio_3c_dc)
ratio_graph_3cat <- rename(ratio_graph_3cat, value=ratio_3c_dc)
ratio_graph_3cat$label <- 'ratio'

graph_seg_3cat <- ggplot(ratio_graph_3cat, aes(x = value, fill = label)) +
  geom_density(position = "identity", alpha = 0.5) +
  scale_fill_manual(values = c("ratio" = "#99152b"), guide = "none") +
  theme_minimal() +
  labs(tag = "b", y = "Density",x="Proportion of cross-social-origin help-seeking ties")

      ### Panel c: classroom social heterogenity ----------
graph_div_3cat  <- dplyr::select(data_dc, simps.3c)
graph_div_3cat   <- rename(graph_div_3cat , value=simps.3c)
graph_div_3cat $label <- "simps.3c"

graph_div_3cat  <- ggplot(graph_div_3cat , aes(x = value, fill = label)) +
  geom_density(position = "identity") +
  scale_fill_manual(values = c("simps.3c" = "#405449"), guide = FALSE) +
  theme_minimal() +
  labs(tag = "c", y = "Density",x="Classroom social heterogeneity")+
  scale_x_continuous(limits = c(0, 1))

      ### Panel d: consolidation 
graph_conso_3cat <- select(data_dc, shaw_eth_dc)
graph_conso_3cat   <- rename(graph_conso_3cat ,value=shaw_eth_dc)
graph_conso_3cat$label <- "shaw_3c_eth"


graph_gen_3cat<- select(data_dc, shaw_gen_dc)
graph_gen_3cat  <- rename(graph_gen_3cat ,value=shaw_gen_dc)
graph_gen_3cat$label <- "shaw_3c_gen"

graph_conso_3cat<- rbind(graph_conso_3cat,graph_gen_3cat)

graph_eth_gen_3cat <- ggplot(graph_conso_3cat, aes(x = value, fill=label)) +
  geom_density(position = "identity", alpha = 0.5) +
  scale_fill_manual(values = c("shaw_3c_eth" = "#4e3227", "shaw_3c_gen" = "#dcbfa2"), guide = FALSE) +
  theme_minimal() +
  labs(tag = "d", x = "Consolidation",y="Density",
       caption = "Dark brown: Consolidation with social and ethnic origin
       Light brown: Consolidation with social origin and sex")


    ## Output 
  library(ggpubr)
  library(patchwork)

figure_3cat <- ggarrange(ineq_graph_3cat,graph_seg_3cat,graph_div_3cat,graph_eth_gen_3cat)

note_text <- paste(
  "Notes:",
    "          Panel a. Distribution of students' reversed grade by social origins (HISEI tertiles)",
    "          Panel b. Distribution of the proportion of cross-social-origin help-seeking ties",
    "          Panel c. Distribution of classroom social heterogeneity",
    "          Panel d. Distribution of social origin consolidation within classrooms",
  sep = "\n"
)

  library(cowplot)

title_row <- ggdraw() +
  draw_label("Figure 1. Distribution of selected variables",
             x = 0, hjust = 0, vjust = 1,size = 12) 

notes_row <- ggdraw() +
  draw_label(note_text,
             x = 0, hjust = 0, vjust = 0, size = 10, lineheight = 1.05)

figure_1 <- plot_grid(
  title_row,
  ggdraw(figure_3cat),
  notes_row,
  ncol = 1,
  rel_heights = c(0.07, 1, 0.25)
)

ggsave("figure 1.svg", figure_1, width = 12, height = 10, dpi = 1200)



  ## ICC -----------

  library(performance)
m0 <- lmer(ratio_3c_dc~ 1 + (1 | IDTESTGROUP_FDZ), data = data_dc, REML = TRUE)
  icc_tbl <- icc(m0)
  print(icc_tbl)

vc <- as.data.frame(VarCorr(m0))
  between_var <- vc$vcov[ vc$grp == "IDTESTGROUP_FDZ" ]
  within_var  <- attr(VarCorr(m0), "sc")^2


m1 <- lmer(TR_NOTE_MAT_recode~ 1 + (1 | IDTESTGROUP_FDZ), data = data_dc, REML = TRUE)
  icc_tbl <- icc(m1)
  print(icc_tbl)

vc <- as.data.frame(VarCorr(m1))
  between_var <- vc$vcov[ vc$grp == "IDTESTGROUP_FDZ" ]
  within_var  <- attr(VarCorr(m1), "sc")^2



# 6. Regression ---------
data_dc_st <- standardise(data_dc)
  
  ## Classroom social heterogeneity and consolidation predicting proportion of cross-social-origin ties ---------

div_ratio_3c_dc<-lmer(ratio_3c_dc~simps.3c+
                        TR_GESCHLECHT+immi+ses_3c+mean_ses+size_orig+same_hisei_3c+
                        (1|IDTESTGROUP_FDZ),
                      data=data_dc_st, REML=T)


div_ratio_3c_ethmod_dc <- lmer(ratio_3c_dc~simps.3c*shaw_eth_dc+
                                 size_orig+TR_GESCHLECHT+immi+ses_3c+mean_ses+same_hisei_3c+
                                 (1|IDTESTGROUP_FDZ),
                               data=data_dc_st, REML=T)


div_ratio_3c_genmod_dc <- lmer(ratio_3c_dc~simps.3c*shaw_gen_dc+
                                 size_orig+TR_GESCHLECHT+immi+ses_3c+mean_ses+same_hisei_3c+
                                 (1|IDTESTGROUP_FDZ),
                               data=data_dc_st, REML=T)



  ## Proportion of cross-social-origin ties predicting achievement inequality ---------

ratio_ineq_3c_dc_fix<- lmer(TR_NOTE_MAT_recode~ses_3c*ratio_3c_dc+
                              immi+TR_GESCHLECHT+gym+mean_ses+size_orig+
                              (ses_3c| IDTESTGROUP_FDZ), 
                            data = data_dc_st,REML=T)  


  ## Direct association between classroom social heterogeneity, consolidation, and achievement inequality ------------

div_ineq_3c_dc <-lmer(TR_NOTE_MAT_recode~ses_3c*simps.3c+
                        immi+TR_GESCHLECHT+gym+mean_ses+size_orig+
                        (ses_3c|IDTESTGROUP_FDZ),
                      data=data_dc_st, REML=T)


div_ineq_3c_ethmod_dc<-lmer(TR_NOTE_MAT_recode ~ ses_3c*simps.3c * shaw_eth_dc + 
                              immi +size_orig + TR_GESCHLECHT +  gym + mean_ses + 
                              (ses_3c | IDTESTGROUP_FDZ),
                            data = data_dc_st, REML = TRUE)


div_ineq_3c_genmod_dc <-lmer(TR_NOTE_MAT_recode ~ ses_3c*simps.3c*shaw_gen_dc+ 
                               immi+size_orig+TR_GESCHLECHT+gym+mean_ses + 
                               (ses_3c | IDTESTGROUP_FDZ),
                             data = data_dc_st, REML = TRUE)

  ### include the proportion as control 

div_ineq_3c_dc_ratio <-lmer(TR_NOTE_MAT_recode~ses_3c*simps.3c+ratio_3c_dc*ses_3c+ 
                              immi+TR_GESCHLECHT+gym+mean_ses+size_orig+
                              (ses_3c|IDTESTGROUP_FDZ),
                            data=data_dc_st, REML=T)


div_ineq_3c_ethmod_dc_ratio<-lmer(TR_NOTE_MAT_recode ~ ses_3c*simps.3c * shaw_eth_dc+ratio_3c_dc*ses_3c+ 
                                    immi +size_orig + TR_GESCHLECHT +  gym + mean_ses + 
                                    (ses_3c | IDTESTGROUP_FDZ),
                                  data = data_dc_st, REML = TRUE)


div_ineq_3c_genmod_dc_ratio <-lmer(TR_NOTE_MAT_recode ~ ses_3c*simps.3c*shaw_gen_dc+ratio_3c_dc*ses_3c+  
                                     immi+size_orig+TR_GESCHLECHT+gym+mean_ses + 
                                     (ses_3c | IDTESTGROUP_FDZ),
                                   data = data_dc_st, REML = TRUE)

  ### Results figure ----------
    #### Figure 2 --------

library(visreg)

ratio_achi_clean <- visreg(
  ratio_ineq_3c_dc_fix,          
  xvar    = "ratio_3c_dc",   
  by      = "ses_3c",        
  overlay = TRUE,          
  band    = TRUE,            
  partial = FALSE, 
  rug     = FALSE,          
  legend  = FALSE,           
  line.par = list(           
    col = c("#55A868","#4C72B0","#C44E52"),
    lwd = 2                  
  ),
  fill.par = list(          
    col = c("#55A86880","#4C72B080","#C44E5280")
  ),
  xlab    = "Proportion of inter-social-origin help-seeking",
  ylab    = "Achievement",
  main  ="Figure 2. Regression plot showing the association between proportion of 
  cross-social-origin help-seeking ties and achievement"
)

par(xpd = TRUE)

legend("topright",
       legend = c("low", "intermediate", "high"),
       title  = "Social Origin Group",
       col    = c("#55A868","#4C72B0",  "#C44E52"),
       lty    = 1,
       bty    = "o") 


library(ggeffects)

# Create predicted data with CI by group
pred <- ggpredict(ratio_ineq_3c_dc_fix, terms = c("ratio_3c_dc", "ses_3c"))

pred$group <- recode_factor(pred$group,
                            `[11,42.3]` = "low",
                            `(42.3,61.8]` = "intermediate",
                            `(61.8,89]` = "high"
)

div_ratio_figure <- ggplot(pred, aes(x = x, y = predicted, color = group, fill = group)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  scale_color_manual(
    values = c("low" ="#55A868" , "intermediate" ="#4C72B0" , "high" = "#C44E52"),
    labels = c("low", "intermediate", "high")
  ) +
  scale_fill_manual(
    values = c("low" = "#55A86880", "intermediate" ="#4C72B080" , "high" = "#C44E5280"),
    labels = c("low", "intermediate", "high")
  ) +
  labs(
    x = "Proportion of inter-social-origin help-seeking",
    y = "Achievement",
    title="Figure 2. Regression plot showing the association between proportion of cross-social-origin help-seeking ties and achievement",
    caption ="Notes: The horizontal axis displays the proportion of inter-social-origin help-seeking with standardized values. 
            The shaded ribbon the 95% confidence interval around the fitted regression lines for each SES group. 
            Full results are provided in the Appendix Table A3.",
    color = "Social origin group",
    fill = "Social origin group"
  ) +
  theme_minimal()+
  theme(
    plot.title.position   = "plot",
    plot.caption.position = "plot",
    plot.title   = element_text(hjust = 0, face = "plain", size = 12, margin = margin(b = 6)),
    plot.caption = element_text(hjust = 0, size = 9, lineheight = 1.05, margin = margin(t = 10)),
    plot.margin  = margin(t = 6, r = 10, b = 16, l = 10)
  )

ggsave("figure 2.svg", div_ratio_figure, width = 12, height = 8, dpi = 1200)
        #### Figure 3 -------------

      library(forcats)
      library(ggplot2)
      library(ggpubr)
      library(emmeans)

data_figure <- data_dc_st %>%
  transmute(
    ses_3c = fct_recode(
      ses_3c,
      low          = "[11,42.3]",
      intermediate = "(42.3,61.8]",
      high         = "(61.8,89]"
    ),
    shaw_gen= shaw_gen_dc,
    shaw_eth= shaw_eth_dc,
    simps_c = simps.3c,
    TR_NOTE_MAT_recode,
    immi, size_orig, TR_GESCHLECHT, gym, mean_ses, IDTESTGROUP_FDZ
  )

simps_seq <- seq(
  from = min(data_figure$simps_c, na.rm=TRUE),
  to   = max(data_figure$simps_c, na.rm=TRUE),
  length.out = 100
)


pred_grid <- expand.grid(
  simps_c = simps_seq,
  ses_3c  = levels(data_figure$ses_3c),
  shaw_eth  = c(-1, 0, 1),
  shaw_gen = c(-1, 0, 1),
  immi          = names(sort(table(data_figure$immi),          decreasing=TRUE))[1],
  TR_GESCHLECHT = names(sort(table(data_figure$TR_GESCHLECHT), decreasing=TRUE))[1],
  gym           = names(sort(table(data_figure$gym),           decreasing=TRUE))[1],
  size_orig     = mean(data_figure$size_orig, na.rm=TRUE),
  mean_ses      = mean(data_figure$mean_ses,  na.rm=TRUE)
)

  
    ##### Panel a: plot simple slopes of simps_c by ses_3c -----------
div_ineq_3c_dc_cat <-lmer(TR_NOTE_MAT_recode~ses_3c*simps_c+
                            immi+TR_GESCHLECHT+gym+mean_ses+size_orig+
                            (ses_3c|IDTESTGROUP_FDZ),
                          data=data_figure, REML=T)

emm_options(lmerTest.limit = 29597)
emm_options(pbkrtest.limit = 29597)

emm_direct <- emmeans(
  div_ineq_3c_dc_cat,
  specs = c("simps_c", "ses_3c"),
  at    = list(simps_c = simps_seq, ses_3c = levels(data_figure$ses_3c)),
  type  = "response"
)

emm_direct <- as.data.frame(emm_direct)

direct_plot <- ggplot(emm_direct, 
                      aes(x = simps_c, y = emmean, color = ses_3c, fill = ses_3c)) +
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL),
              alpha = 0.2, colour = NA) +
  geom_line(size = 1.5) +
  scale_color_manual(
    values = c("low" ="#55A868" , "intermediate" ="#4C72B0" , "high" = "#C44E52"),
    labels = c("low", "intermediate", "high")
  ) +
  scale_fill_manual(
    values = c("low" = "#55A86880", "intermediate" ="#4C72B080" , "high" = "#C44E5280"),
    labels = c("low", "intermediate", "high")
  ) +
  labs(
    x     = "Classroom Social Diversity",
    y     = "Achievement",
    title = "a. Classroom social diversity predicting achievement by social origin"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "right",
    strip.text      = element_text(face = "bold", size = 10),
    plot.title = element_text(size = 10, face = "plain", hjust = 0, margin = margin(b = 4))
  )


    ##### Panel b: ethnic consolidation interaction -----------

data_figure <- data_figure %>%
  mutate(
    shaw_bin_eth = case_when(
      shaw_eth<= -1 ~ "-1 SD",
      shaw_eth>=  1 ~ "+1 SD",
      TRUE         ~ "Mean"
    ) %>% factor(levels = c("-1 SD","Mean","+1 SD"))
  )

div_ineq_3c_ethmod_cat<-lmer(TR_NOTE_MAT_recode ~ ses_3c*simps_c * shaw_eth + 
                               immi +size_orig + TR_GESCHLECHT +  gym + mean_ses + 
                               (ses_3c | IDTESTGROUP_FDZ),
                             data = data_figure, REML = TRUE)

emm_eth <- emmeans(
  div_ineq_3c_ethmod_cat,
  specs = c("simps_c", "ses_3c", "shaw_eth"),
  at    = list(simps_c = simps_seq, 
               ses_3c = levels(data_figure$ses_3c), 
               shaw_eth = c(-1,0,1)),
  type  = "response"
)

emm_eth <- as.data.frame(emm_eth) %>%
  mutate(
    shaw_bin_eth = factor(shaw_eth, labels = c("-1 SD","Mean","+1 SD"))
  )


interaction_plot_eth <- ggplot(emm_eth, aes(x = simps_c, y = emmean, color = ses_3c, fill = ses_3c)) +
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL),
              alpha = 0.2, colour = NA) +
  geom_line(size = 1.5) +
  facet_wrap(~ shaw_bin_eth, nrow = 1) +
  scale_color_manual(
    values = c("low" ="#55A868" , "intermediate" ="#4C72B0" , "high" = "#C44E52"),
    labels = c("low", "intermediate", "high")
  ) +
  scale_fill_manual(
    values = c("low" = "#55A86880", "intermediate" ="#4C72B080" , "high" = "#C44E5280"),
    labels = c("low", "intermediate", "high")
  ) +
  labs(
    x = "Classroom Social Diversity",
    y = "Achievement",
    title = "b. Interaction of social origin, classroom social heterogeneity, and consolidation (with ethnic origin) in predicting achievement"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    strip.text      = element_text(face = "bold", size = 10),
    legend.position = "right",
    plot.title = element_text(size = 10, face = "plain", hjust = 0, margin = margin(b = 4))
  )


    ##### Panel c: gender consolidation interaction -----------------
data_figure <- data_figure %>%
  mutate(
    shaw_bin_gen = case_when(
      shaw_gen<= -1 ~ "-1 SD",
      shaw_gen>=  1 ~ "+1 SD",
      TRUE         ~ "Mean"
    ) %>% factor(levels = c("-1 SD","Mean","+1 SD"))
  ) 

div_ineq_3c_genmod_dc_cat<-lmer(TR_NOTE_MAT_recode ~ 
                                  ses_3c*simps_c*shaw_gen+ 
                                  immi+size_orig+TR_GESCHLECHT+gym+mean_ses+ 
                                  (ses_3c | IDTESTGROUP_FDZ),
                                data = data_figure, REML = TRUE)


emm_gen <- emmeans(
  div_ineq_3c_genmod_dc_cat,
  specs = c("simps_c", "ses_3c", "shaw_gen"),
  at    = list(simps_c = simps_seq, ses_3c = levels(data_figure$ses_3c), shaw_gen= c(-1,0,1)),
  type  = "response"
)

emm_gen <- as.data.frame(emm_gen) %>%
  mutate(
    shaw_bin_gen = factor(shaw_gen, labels = c("-1 SD","Mean","+1 SD"))
  )


interaction_plot_gen <- ggplot(emm_gen, aes(x = simps_c, y = emmean, color = ses_3c, fill = ses_3c)) +
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL),
              alpha = 0.2, colour = NA) +
  geom_line(size = 1.5) +
  facet_wrap(~ shaw_bin_gen, nrow = 1) +
  scale_color_manual(
    values = c("low" ="#55A868" , "intermediate" ="#4C72B0" , "high" = "#C44E52"),
    labels = c("low", "intermediate", "high")
  ) +
  scale_fill_manual(
    values = c("low" = "#55A86880", "intermediate" ="#4C72B080" , "high" = "#C44E5280"),
    labels = c("low", "intermediate", "high")
  ) +
  labs(
    x = "Classroom Social Diversity",
    y = "Achievement",
    title = "c. Interaction of social origin, classroom social heterogeneity, and consolidation (with sex) in predicting achievement"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    strip.text      = element_text(face = "bold", size = 10),
    legend.position = "right",
    plot.title = element_text(size = 10, face = "plain", hjust = 0, margin = margin(b = 4))
  )


  ### output
library(ggpubr)
interaction_plot <- ggarrange(direct_plot,interaction_plot_eth,interaction_plot_gen,
                              nrow=3,ncol=1,
                              common.legend= TRUE,
                              legend = "right")

interaction_plot<- annotate_figure(
  interaction_plot,
  top = text_grob(
    "Figure 3. Regression plots showing the direct associations between classroom social heterogeneity and achievement, 
including consolidation with ethnic origin and sex as interaction terms",,
    size = 12,
    hjust = 0,
    x=0,
    y=0.5
  ),
  bottom = text_grob(
    "Notes: All continues variables are z-standardized. Shaded bands are 95% CIs.
         -1 SD, Mean, and +1 SD refer to the value of consolidation.
          Full results table are provided in Appendix A4.",
    hjust = 0,
    x     = 0,
    size  = 10
  )
)
ggsave(interaction_plot, file = "Figure 3.svg", width = 12, height = 8, dpi = 1200)




# 7. Robustness check ------------------

## 7.1 Using different categorizations of HISEI --------------

      ### variables  -----------
library(gtools)
data_dc$ses_4c <- quantcut(data_dc$HISEI, 4)
data_dc$ses_5c <- quantcut(data_dc$HISEI, 5)
data_dc$ses_6c <- quantcut(data_dc$HISEI, 6)
data_dc$ses_7c <- quantcut(data_dc$HISEI, 7)
data_dc$ses_8c <- quantcut(data_dc$HISEI, 8)
data_dc$ses_9c <- quantcut(data_dc$HISEI, 9)
data_dc$ses_10c <- quantcut(data_dc$HISEI, 10)

breaks <- quantile(data_dc$HISEI, probs = seq(0, 1, length.out = 5), na.rm = TRUE)
  network_all$ses_4c<- cut(network_all$HISEI, breaks = breaks, include.lowest = T)

breaks <- quantile(data_dc$HISEI, probs = seq(0, 1, length.out = 6), na.rm = TRUE)
  network_all$ses_5c<- cut(network_all$HISEI, breaks = breaks, include.lowest = T)

breaks <- quantile(data_dc$HISEI, probs = seq(0, 1, length.out = 7), na.rm = TRUE)
  network_all$ses_6c<- cut(network_all$HISEI, breaks = breaks, include.lowest = T)

breaks <- quantile(data_dc$HISEI, probs = seq(0, 1, length.out = 8), na.rm = TRUE)
  network_all$ses_7c<- cut(network_all$HISEI, breaks = breaks, include.lowest = T)

breaks <- quantile(data_dc$HISEI, probs = seq(0, 1, length.out = 9), na.rm = TRUE)
  network_all$ses_8c<- cut(network_all$HISEI, breaks = breaks, include.lowest = T)

breaks <- quantile(data_dc$HISEI, probs = seq(0, 1, length.out = 10), na.rm = TRUE)
  network_all$ses_9c<- cut(network_all$HISEI, breaks = breaks, include.lowest = T)

breaks <- quantile(data_dc$HISEI, probs = seq(0, 1, length.out = 11), na.rm = TRUE)
  network_all$ses_10c<- cut(network_all$HISEI, breaks = breaks, include.lowest = T)

network_all_cl<- split(network_all, network_all$IDTESTGROUP_FDZ)
  network_all_cl<- lapply(network_all_cl, function(df) {
  df[order(df$IDinClass), ]
})
  
  
  #### 4 categories  ---------------------
simps.4c <- simps(data_dc, "ses_4c")
  simps.4c<- rename(simps.4c, simps.4c=`diversity(simpdat3, index = "simpson")`)
  simps.4c <- rename(simps.4c, IDTESTGROUP_FDZ=names)
    simps.4c $IDTESTGROUP_FDZ <- as.character(simps.4c $IDTESTGROUP_FDZ)

  data_dc <- left_join(data_dc, simps.4c)
  
  
ratio_ses_4c_dc <- function(group) {
    individual_results <- list()
    
    help_columns <- paste0("help_", 1:34)
    
    for (i in 1:nrow(group)) {
      same_ties <- 0
      inter_ties <- 0
      total_ties <- group$total_tie[i]
      
      for (help_col in help_columns) {
        adjacency_value <- group[i, help_col]
        
        row_index <- which(colnames(group) == help_col) - 3  
        
        if (!is.na(adjacency_value) && adjacency_value == 1 && 
            !is.na(group$ses_4c[i]) && !is.na(group$ses_4c[row_index])) {
          
        if (group$ses_4c[i] == group$ses_4c[row_index]) {
            same_ties <- same_ties + 1
          }
        } 
      }
      

      inter_ties <- total_ties - same_ties
      ratio <- if (total_ties > 0) {
        if (is.na(inter_ties / total_ties)) {
          NA
        } else {
          inter_ties / total_ties
        }
      } else {
        NA
      }
      

      individual_results[[i]] <- list(IDTESTGROUP_FDZ = group$IDTESTGROUP_FDZ[i], 
                                      IDSTUD_FDZ = group$IDSTUD_FDZ[i], 
                                      ratio_4c_dc = ratio, 
                                      inter_tie_4c_dc = inter_ties, 
                                      same_tie_4c_dc = same_ties)
    }
    
    return(individual_results)
  }
  
ratio_ses4cat<- lapply(network_all_cl, ratio_ses_4c_dc)
  ratio_ses4cat<<- lapply(ratio_ses4cat, bind_rows)%>%
    bind_rows(.)
    ratio_ses4cat$IDSTUD_FDZ <- as.character(ratio_ses4cat$IDSTUD_FDZ)
    ratio_ses4cat$IDTESTGROUP_FDZ <- as.character(ratio_ses4cat$IDTESTGROUP_FDZ)

  data_dc<- left_join(data_dc, ratio_ses4cat, by=c("IDSTUD_FDZ", "IDTESTGROUP_FDZ"))
  
data_dc$ses_4c <- as.factor(data_dc$ses_4c)
data_dc<- data_dc %>%
  group_by(IDTESTGROUP_FDZ, ses_4c) %>%
  mutate(same_hisei_4c = n() - 1) %>%
  ungroup()

data_conso_eth_4c <- dplyr::select(data_dc, IDTESTGROUP_FDZ,herkunft_FDZ, ses_4c)%>%as.data.frame()
  data_conso_eth_4c$herkunft_FDZ <- as.factor(data_conso_eth_4c$herkunft_FDZ)

shaw_eth_4c_dc <- asw.cluster::faultlines(data_conso_eth_4c, 
                                          group.par = c("IDTESTGROUP_FDZ"), 
                                          attr.type = c("nominal","nominal"), 
                                          rescale = "sd",
                                          method = "shaw",
                                          usesghomo = FALSE)


shaw_eth_4c_dc_long <- summary(shaw_eth_4c_dc)$long
  shaw_eth_4c_dc_indi <- select(shaw_eth_4c_dc_long, c(team,fl.value))
  shaw_eth_4c_dc_indi <- rename(shaw_eth_4c_dc_indi, IDTESTGROUP_FDZ=team, shaw_eth_4c=fl.value)

shaw_eth_4c_dc_cl <- shaw_eth_4c_dc_indi%>%
  group_by(IDTESTGROUP_FDZ) %>%
  summarise(shaw_eth_4c = unique(shaw_eth_4c))

data_dc<- left_join(data_dc, shaw_eth_4c_dc_cl)


data_conso_gen_4c <- dplyr::select(data_dc, IDTESTGROUP_FDZ,TR_GESCHLECHT, ses_4c)%>%as.data.frame()
  data_conso_gen_4c$TR_GESCHLECHT <- as.factor(data_conso_gen_4c$TR_GESCHLECHT)


shaw_gen_4c_dc<- asw.cluster::faultlines(data_conso_gen_4c, 
                                         group.par = c("IDTESTGROUP_FDZ"), 
                                         attr.type = c("nominal","nominal"), 
                                         rescale = "sd",
                                         method = "shaw",
                                         usesghomo = FALSE)

shaw_gen_4c_dc_long <- summary(shaw_gen_4c_dc)$long
shaw_gen_4c_dc_indi <- select(shaw_gen_4c_dc_long, c(team,fl.value))
shaw_gen_4c_dc_indi <- rename(shaw_gen_4c_dc_indi, IDTESTGROUP_FDZ=team, shaw_gen_4c=fl.value)

shaw_gen_4c_dc_cl <- shaw_gen_4c_dc_indi%>%
  group_by(IDTESTGROUP_FDZ) %>%
  summarise(shaw_gen_4c = unique(shaw_gen_4c))

data_dc<- left_join(data_dc, shaw_gen_4c_dc_cl)
  
    #### 5 categories  ------------

simps.5c <- simps(data_dc, "ses_5c")
  simps.5c<- rename(simps.5c, simps.5c=`diversity(simpdat3, index = "simpson")`)
  simps.5c <- rename(simps.5c, IDTESTGROUP_FDZ=names)
  simps.5c $IDTESTGROUP_FDZ <- as.character(simps.5c $IDTESTGROUP_FDZ)

  
  data_dc <- left_join(data_dc, simps.5c)
  
ratio_ses_5c_dc <- function(group) {
    individual_results <- list()
    
    help_columns <- paste0("help_", 1:34)
    
    for (i in 1:nrow(group)) {
      same_ties <- 0
      inter_ties <- 0
      total_ties <- group$total_tie[i]
      

      for (help_col in help_columns) {
        adjacency_value <- group[i, help_col]
        
        row_index <- which(colnames(group) == help_col) - 3  
        
        if (!is.na(adjacency_value) && adjacency_value == 1 && 
            !is.na(group$ses_5c[i]) && !is.na(group$ses_5c[row_index])) {
          
        if (group$ses_5c[i] == group$ses_5c[row_index]) {
            same_ties <- same_ties + 1
          }
        } 
      }
      
      inter_ties <- total_ties - same_ties
 
          ratio <- if (total_ties > 0) {
        if (is.na(inter_ties / total_ties)) {
          NA
        } else {
          inter_ties / total_ties
        }
      } else {
        NA
      }

      individual_results[[i]] <- list(IDTESTGROUP_FDZ = group$IDTESTGROUP_FDZ[i], 
                                      IDSTUD_FDZ = group$IDSTUD_FDZ[i], 
                                      ratio_5c_dc = ratio, 
                                      inter_tie_5c_dc = inter_ties, 
                                      same_tie_5c_dc = same_ties)
    }
    
    return(individual_results)
  }
  
ratio_ses5cat<- lapply(network_all_cl, ratio_ses_5c_dc)
  ratio_ses5cat<<- lapply(ratio_ses5cat, bind_rows)%>%
  bind_rows(.)
  ratio_ses5cat$IDSTUD_FDZ <- as.character(ratio_ses5cat$IDSTUD_FDZ)
  ratio_ses5cat$IDTESTGROUP_FDZ <- as.character(ratio_ses5cat$IDTESTGROUP_FDZ)


  data_dc<- left_join(data_dc, ratio_ses5cat, by=c("IDSTUD_FDZ", "IDTESTGROUP_FDZ"))
  
data_dc$ses_5c <- as.factor(data_dc$ses_5c)
data_dc<- data_dc %>%
  group_by(IDTESTGROUP_FDZ, ses_5c) %>%
  mutate(same_hisei_5cat = n() - 1) %>%
  ungroup()


data_conso_eth_5c <- dplyr::select(data_dc, IDTESTGROUP_FDZ,herkunft_FDZ, ses_5c)%>%as.data.frame()
  data_conso_eth_5c$herkunft_FDZ <- as.factor(data_conso_eth_5c$herkunft_FDZ)

shaw_eth_5c_dc <- asw.cluster::faultlines(data_conso_eth_5c, 
                                          group.par = c("IDTESTGROUP_FDZ"), 
                                          attr.type = c("nominal","nominal"), 
                                          rescale = "sd",
                                          method = "shaw",
                                          usesghomo = FALSE)


data_conso_gen_5c <- dplyr::select(data_dc, IDTESTGROUP_FDZ,TR_GESCHLECHT, ses_5c)%>%as.data.frame()
  data_conso_gen_5c$TR_GESCHLECHT <- as.factor(data_conso_gen_5c$TR_GESCHLECHT)


shaw_gen_5c_dc<- asw.cluster::faultlines(data_conso_gen_5c, 
                                         group.par = c("IDTESTGROUP_FDZ"), 
                                         attr.type = c("nominal","nominal"), 
                                         rescale = "sd",
                                         method = "shaw",
                                         usesghomo = FALSE)

shaw_eth_5c_dc_long <- summary(shaw_eth_5c_dc)$long
  shaw_eth_5c_dc_indi <- select(shaw_eth_5c_dc_long, c(team,fl.value))
  shaw_eth_5c_dc_indi <- rename(shaw_eth_5c_dc_indi, IDTESTGROUP_FDZ=team, shaw_eth_5c=fl.value)

shaw_eth_5c_dc_cl <- shaw_eth_5c_dc_indi%>%
  group_by(IDTESTGROUP_FDZ) %>%
  summarise(shaw_eth_5c = unique(shaw_eth_5c))

data_dc<- left_join(data_dc, shaw_eth_5c_dc_cl)

shaw_gen_5c_dc_long <- summary(shaw_gen_5c_dc)$long
  shaw_gen_5c_dc_indi <- select(shaw_gen_5c_dc_long, c(team,fl.value))
  shaw_gen_5c_dc_indi <- rename(shaw_gen_5c_dc_indi, IDTESTGROUP_FDZ=team, shaw_gen_5c=fl.value)

shaw_gen_5c_dc_cl <- shaw_gen_5c_dc_indi%>%
  group_by(IDTESTGROUP_FDZ) %>%
  summarise(shaw_gen_5c = unique(shaw_gen_5c))

  data_dc<- left_join(data_dc, shaw_gen_5c_dc_cl)

      #### 6 categories  ------------

simps.6c <- simps(data_dc, "ses_6c")
  simps.6c<- rename(simps.6c, simps.6c=`diversity(simpdat3, index = "simpson")`)
  simps.6c <- rename(simps.6c, IDTESTGROUP_FDZ=names)
  simps.6c $IDTESTGROUP_FDZ <- as.character(simps.6c $IDTESTGROUP_FDZ)
  
  data_dc <- left_join(data_dc, simps.6c)

ratio_ses_6c_dc <- function(group) {
    individual_results <- list()

    help_columns <- paste0("help_", 1:34)
    
    for (i in 1:nrow(group)) {
      same_ties <- 0
      inter_ties <- 0
      total_ties <- group$total_tie[i]
      
      for (help_col in help_columns) {
        adjacency_value <- group[i, help_col]
        
        row_index <- which(colnames(group) == help_col) - 3 
        
        if (!is.na(adjacency_value) && adjacency_value == 1 && 
            !is.na(group$ses_6c[i]) && !is.na(group$ses_6c[row_index])) {
          
        if (group$ses_6c[i] == group$ses_6c[row_index]) {
            same_ties <- same_ties + 1
          }
        } 
      }
      
      inter_ties <- total_ties - same_ties
      
      ratio <- if (total_ties > 0) {
        if (is.na(inter_ties / total_ties)) {
          NA
        } else {
          inter_ties / total_ties
        }
      } else {
        NA
      }
      
      individual_results[[i]] <- list(IDTESTGROUP_FDZ = group$IDTESTGROUP_FDZ[i], 
                                      IDSTUD_FDZ = group$IDSTUD_FDZ[i], 
                                      ratio_6c_dc = ratio, 
                                      inter_tie_6c_dc = inter_ties, 
                                      same_tie_6c_dc = same_ties)
    }
    
    return(individual_results)
}

ratio_ses6cat<- lapply(network_all_cl, ratio_ses_6c_dc)
  ratio_ses6cat<<- lapply(ratio_ses6cat, bind_rows)%>%
  bind_rows(.)

  ratio_ses6cat$IDSTUD_FDZ <- as.character(ratio_ses6cat$IDSTUD_FDZ)
  ratio_ses6cat$IDTESTGROUP_FDZ <- as.character(ratio_ses6cat$IDTESTGROUP_FDZ)

  data_dc<- left_join(data_dc, ratio_ses6cat, by=c("IDSTUD_FDZ", "IDTESTGROUP_FDZ"))

data_dc$ses_6c <- as.factor(data_dc$ses_6c)
data_dc<- data_dc %>%
  group_by(IDTESTGROUP_FDZ, ses_6c) %>%
  mutate(same_hisei_6cat = n() - 1) %>%
  ungroup()
  
data_conso_eth_6c <- dplyr::select(data_dc, IDTESTGROUP_FDZ,herkunft_FDZ, ses_6c)%>%as.data.frame()
  data_conso_eth_6c$herkunft_FDZ <- as.factor(data_conso_eth_6c$herkunft_FDZ)

shaw_eth_6c_dc <- asw.cluster::faultlines(data_conso_eth_6c, 
                                          group.par = c("IDTESTGROUP_FDZ"), 
                                          attr.type = c("nominal","nominal"), 
                                          rescale = "sd",
                                          method = "shaw",
                                          usesghomo = FALSE)


data_conso_gen_6c <- dplyr::select(data_dc, IDTESTGROUP_FDZ,TR_GESCHLECHT, ses_6c)%>%as.data.frame()
  data_conso_gen_6c$TR_GESCHLECHT <- as.factor(data_conso_gen_6c$TR_GESCHLECHT)


shaw_gen_6c_dc<- asw.cluster::faultlines(data_conso_gen_6c, 
                                         group.par = c("IDTESTGROUP_FDZ"), 
                                         attr.type = c("nominal","nominal"), 
                                         rescale = "sd",
                                         method = "shaw",
                                         usesghomo = FALSE)


shaw_eth_6c_dc_long <- summary(shaw_eth_6c_dc)$long
  shaw_eth_6c_dc_indi <- select(shaw_eth_6c_dc_long, c(team,fl.value))
  shaw_eth_6c_dc_indi <- rename(shaw_eth_6c_dc_indi, IDTESTGROUP_FDZ=team, shaw_eth_6c=fl.value)

shaw_eth_6c_dc_cl <- shaw_eth_6c_dc_indi%>%
  group_by(IDTESTGROUP_FDZ) %>%
  summarise(shaw_eth_6c = unique(shaw_eth_6c))

  data_dc<- left_join(data_dc, shaw_eth_6c_dc_cl)

shaw_gen_6c_dc_long <- summary(shaw_gen_6c_dc)$long
  shaw_gen_6c_dc_indi <- select(shaw_gen_6c_dc_long, c(team,fl.value))
  shaw_gen_6c_dc_indi <- rename(shaw_gen_6c_dc_indi, IDTESTGROUP_FDZ=team, shaw_gen_6c=fl.value)

shaw_gen_6c_dc_cl <- shaw_gen_6c_dc_indi%>%
  group_by(IDTESTGROUP_FDZ) %>%
  summarise(shaw_gen_6c = unique(shaw_gen_6c))

  data_dc<- left_join(data_dc, shaw_gen_6c_dc_cl)

  
      #### 7 categories ----------
simps.7c <- simps(data_dc, "ses_7c")
  simps.7c<- rename(simps.7c, simps.7c=`diversity(simpdat3, index = "simpson")`)
  simps.7c <- rename(simps.7c, IDTESTGROUP_FDZ=names)
  simps.7c $IDTESTGROUP_FDZ <- as.character(simps.7c $IDTESTGROUP_FDZ)
  
  data_dc <- left_join(data_dc, simps.7c)
  
  ratio_ses_7c_dc <- function(group) {
    individual_results <- list()
    
    help_columns <- paste0("help_", 1:34)
    
    for (i in 1:nrow(group)) {
      same_ties <- 0
      inter_ties <- 0
      total_ties <- group$total_tie[i]
      
      for (help_col in help_columns) {
        adjacency_value <- group[i, help_col]
        
        row_index <- which(colnames(group) == help_col) - 3  
        
        if (!is.na(adjacency_value) && adjacency_value == 1 && 
            !is.na(group$ses_7c[i]) && !is.na(group$ses_7c[row_index])) {
          

        if (group$ses_7c[i] == group$ses_7c[row_index]) {
            same_ties <- same_ties + 1
          }
        } 
      }
      
      inter_ties <- total_ties - same_ties
      
      ratio <- if (total_ties > 0) {
        if (is.na(inter_ties / total_ties)) {
          NA
        } else {
          inter_ties / total_ties
        }
      } else {
        NA
      }
      
      individual_results[[i]] <- list(IDTESTGROUP_FDZ = group$IDTESTGROUP_FDZ[i], 
                                      IDSTUD_FDZ = group$IDSTUD_FDZ[i], 
                                      ratio_7c_dc = ratio, 
                                      inter_tie_7c_dc = inter_ties, 
                                      same_tie_7c_dc = same_ties)
    }
    
    return(individual_results)
  }
  
ratio_ses7cat<- lapply(network_all_cl, ratio_ses_7c_dc)
  ratio_ses7cat<<- lapply(ratio_ses7cat, bind_rows)%>%
    bind_rows(.)
  
  ratio_ses7cat$IDSTUD_FDZ <- as.character(ratio_ses7cat$IDSTUD_FDZ)
  ratio_ses7cat$IDTESTGROUP_FDZ <- as.character(ratio_ses7cat$IDTESTGROUP_FDZ)
 
  data_dc<- left_join(data_dc, ratio_ses7cat, by=c("IDSTUD_FDZ", "IDTESTGROUP_FDZ"))
  
data_dc$ses_7c <- as.factor(data_dc$ses_7c)
data_dc<- data_dc %>%
    group_by(IDTESTGROUP_FDZ, ses_7c) %>%
    mutate(same_hisei_7cat = n() - 1) %>%
    ungroup()


data_conso_eth_7c <- dplyr::select(data_dc, IDTESTGROUP_FDZ,herkunft_FDZ, ses_7c)%>%as.data.frame()
  data_conso_eth_7c$herkunft_FDZ <- as.factor(data_conso_eth_7c$herkunft_FDZ)

shaw_eth_7c_dc <- asw.cluster::faultlines(data_conso_eth_7c, 
                                          group.par = c("IDTESTGROUP_FDZ"), 
                                          attr.type = c("nominal","nominal"), 
                                          rescale = "sd",
                                          method = "shaw",
                                          usesghomo = FALSE)


data_conso_gen_7c <- dplyr::select(data_dc, IDTESTGROUP_FDZ,TR_GESCHLECHT, ses_7c)%>%as.data.frame()
  data_conso_gen_7c$TR_GESCHLECHT <- as.factor(data_conso_gen_7c$TR_GESCHLECHT)


shaw_gen_7c_dc<- asw.cluster::faultlines(data_conso_gen_7c, 
                                         group.par = c("IDTESTGROUP_FDZ"), 
                                         attr.type = c("nominal","nominal"), 
                                         rescale = "sd",
                                         method = "shaw",
                                         usesghomo = FALSE)


shaw_eth_7c_dc_long <- summary(shaw_eth_7c_dc)$long
  shaw_eth_7c_dc_indi <- select(shaw_eth_7c_dc_long, c(team,fl.value))
  shaw_eth_7c_dc_indi <- rename(shaw_eth_7c_dc_indi, IDTESTGROUP_FDZ=team, shaw_eth_7c=fl.value)

shaw_eth_7c_dc_cl <- shaw_eth_7c_dc_indi%>%
  group_by(IDTESTGROUP_FDZ) %>%
  summarise(shaw_eth_7c = unique(shaw_eth_7c))

data_dc<- left_join(data_dc, shaw_eth_7c_dc_cl)

shaw_gen_7c_dc_long <- summary(shaw_gen_7c_dc)$long
  shaw_gen_7c_dc_indi <- select(shaw_gen_7c_dc_long, c(team,fl.value))
  shaw_gen_7c_dc_indi <- rename(shaw_gen_7c_dc_indi, IDTESTGROUP_FDZ=team, shaw_gen_7c=fl.value)

shaw_gen_7c_dc_cl <- shaw_gen_7c_dc_indi%>%
  group_by(IDTESTGROUP_FDZ) %>%
  summarise(shaw_gen_7c = unique(shaw_gen_7c))

  data_dc<- left_join(data_dc, shaw_gen_7c_dc_cl)

  
      #### 8 categories ------------
simps.8c <- simps(data_dc, "ses_8c")
  simps.8c<- rename(simps.8c, simps.8c=`diversity(simpdat3, index = "simpson")`)
  simps.8c <- rename(simps.8c, IDTESTGROUP_FDZ=names)
  simps.8c $IDTESTGROUP_FDZ <- as.character(simps.8c $IDTESTGROUP_FDZ)
  describe(simps.8c) 
  
  data_dc <- left_join(data_dc, simps.8c)
  
  
ratio_ses_8c_dc <- function(group) {

    individual_results <- list()
    
    help_columns <- paste0("help_", 1:34)
    
    for (i in 1:nrow(group)) {
      same_ties <- 0
      inter_ties <- 0
      total_ties <- group$total_tie[i]
      
      for (help_col in help_columns) {
        adjacency_value <- group[i, help_col]
        
      row_index <- which(colnames(group) == help_col) - 3  
      
        if (!is.na(adjacency_value) && adjacency_value == 1 && 
            !is.na(group$ses_8c[i]) && !is.na(group$ses_8c[row_index])) {
          
        if (group$ses_8c[i] == group$ses_8c[row_index]) {
            same_ties <- same_ties + 1
          }
        } 
      }
      
      inter_ties <- total_ties - same_ties
      
      ratio <- if (total_ties > 0) {
        if (is.na(inter_ties / total_ties)) {
          NA
        } else {
          inter_ties / total_ties
        }
      } else {
        NA
      }
      
      individual_results[[i]] <- list(IDTESTGROUP_FDZ = group$IDTESTGROUP_FDZ[i], 
                                      IDSTUD_FDZ = group$IDSTUD_FDZ[i], 
                                      ratio_8c_dc = ratio, 
                                      inter_tie_8c_dc = inter_ties, 
                                      same_tie_8c_dc = same_ties)
    }
    
    return(individual_results)
}


ratio_ses8cat<- lapply(network_all_cl, ratio_ses_8c_dc)
  ratio_ses8cat<<- lapply(ratio_ses8cat, bind_rows)%>%
  bind_rows(.)

  ratio_ses8cat$IDSTUD_FDZ <- as.character(ratio_ses8cat$IDSTUD_FDZ)
  ratio_ses8cat$IDTESTGROUP_FDZ <- as.character(ratio_ses8cat$IDTESTGROUP_FDZ)

  data_dc<- left_join(data_dc, ratio_ses8cat, by=c("IDSTUD_FDZ", "IDTESTGROUP_FDZ"))

data_dc$ses_8c <- as.factor(data_dc$ses_8c)
data_dc<- data_dc %>%
  group_by(IDTESTGROUP_FDZ, ses_8c) %>%
  mutate(same_hisei_8cat = n() - 1) %>%
  ungroup() 


data_conso_eth_8c <- dplyr::select(data_dc, IDTESTGROUP_FDZ,herkunft_FDZ, ses_8c)%>%as.data.frame()
  data_conso_eth_8c$herkunft_FDZ <- as.factor(data_conso_eth_8c$herkunft_FDZ)

shaw_eth_8c_dc <- asw.cluster::faultlines(data_conso_eth_8c, 
                                          group.par = c("IDTESTGROUP_FDZ"), 
                                          attr.type = c("nominal","nominal"), 
                                          rescale = "sd",
                                          method = "shaw",
                                          usesghomo = FALSE)


data_conso_gen_8c <- dplyr::select(data_dc, IDTESTGROUP_FDZ,TR_GESCHLECHT, ses_8c)%>%as.data.frame()
  data_conso_gen_8c$TR_GESCHLECHT <- as.factor(data_conso_gen_8c$TR_GESCHLECHT)


shaw_gen_8c_dc<- asw.cluster::faultlines(data_conso_gen_8c, 
                                         group.par = c("IDTESTGROUP_FDZ"), 
                                         attr.type = c("nominal","nominal"), 
                                         rescale = "sd",
                                         method = "shaw",
                                         usesghomo = FALSE)


shaw_eth_8c_dc_long <- summary(shaw_eth_8c_dc)$long
  shaw_eth_8c_dc_indi <- select(shaw_eth_8c_dc_long, c(team,fl.value))
  shaw_eth_8c_dc_indi <- rename(shaw_eth_8c_dc_indi, IDTESTGROUP_FDZ=team, shaw_eth_8c=fl.value)

shaw_eth_8c_dc_cl <- shaw_eth_8c_dc_indi%>%
  group_by(IDTESTGROUP_FDZ) %>%
  summarise(shaw_eth_8c = unique(shaw_eth_8c))

  data_dc<- left_join(data_dc, shaw_eth_8c_dc_cl)


shaw_gen_8c_dc_long <- summary(shaw_gen_8c_dc)$long
  shaw_gen_8c_dc_indi <- select(shaw_gen_8c_dc_long, c(team,fl.value))
shaw_gen_8c_dc_indi <- rename(shaw_gen_8c_dc_indi, IDTESTGROUP_FDZ=team, shaw_gen_8c=fl.value)

shaw_gen_8c_dc_cl <- shaw_gen_8c_dc_indi%>%
  group_by(IDTESTGROUP_FDZ) %>%
  summarise(shaw_gen_8c = unique(shaw_gen_8c))

  data_dc<- left_join(data_dc, shaw_gen_8c_dc_cl)
  
  
      #### 9 categories -------------
simps.9c <- simps(data_dc, "ses_9c")
  simps.9c<- rename(simps.9c, simps.9c=`diversity(simpdat3, index = "simpson")`)
  simps.9c <- rename(simps.9c, IDTESTGROUP_FDZ=names)
  simps.9c $IDTESTGROUP_FDZ <- as.character(simps.9c $IDTESTGROUP_FDZ)
  
  data_dc <- left_join(data_dc, simps.9c)
  

ratio_ses_9c_dc <- function(group) {

    individual_results <- list()
    
    help_columns <- paste0("help_", 1:34)
    
    for (i in 1:nrow(group)) {
      same_ties <- 0
      inter_ties <- 0
      total_ties <- group$total_tie[i]
      
      for (help_col in help_columns) {
        adjacency_value <- group[i, help_col]
        
        row_index <- which(colnames(group) == help_col) - 3  
        
        if (!is.na(adjacency_value) && adjacency_value == 1 && 
            !is.na(group$ses_9c[i]) && !is.na(group$ses_9c[row_index])) {
          
        if (group$ses_9c[i] == group$ses_9c[row_index]) {
            same_ties <- same_ties + 1
          }
        } 
      }
      
      inter_ties <- total_ties - same_ties
      
      ratio <- if (total_ties > 0) {
        if (is.na(inter_ties / total_ties)) {
          NA
        } else {
          inter_ties / total_ties
        }
      } else {
        NA
      }
      
      individual_results[[i]] <- list(IDTESTGROUP_FDZ = group$IDTESTGROUP_FDZ[i], 
                                      IDSTUD_FDZ = group$IDSTUD_FDZ[i], 
                                      ratio_9c_dc = ratio, 
                                      inter_tie_9c_dc = inter_ties, 
                                      same_tie_9c_dc = same_ties)
    }
    
    return(individual_results)
}

ratio_ses9cat<- lapply(network_all_cl, ratio_ses_9c_dc)
ratio_ses9cat<<- lapply(ratio_ses9cat, bind_rows)%>%
  bind_rows(.)

ratio_ses9cat$IDSTUD_FDZ <- as.character(ratio_ses9cat$IDSTUD_FDZ)
ratio_ses9cat$IDTESTGROUP_FDZ <- as.character(ratio_ses9cat$IDTESTGROUP_FDZ)
data_dc <- select(data_dc,-c(ratio_9c_dc))
data_dc<- left_join(data_dc, ratio_ses9cat, by=c("IDSTUD_FDZ", "IDTESTGROUP_FDZ"))
data_dc$ses_9c <- as.factor(data_dc$ses_9c)

data_dc<- data_dc %>%
  group_by(IDTESTGROUP_FDZ, ses_9c) %>%
  mutate(same_hisei_9cat = n() - 1) %>%
  ungroup()


data_conso_eth_9c <- dplyr::select(data_dc, IDTESTGROUP_FDZ,herkunft_FDZ, ses_9c)%>%as.data.frame()
data_conso_eth_9c$herkunft_FDZ <- as.factor(data_conso_eth_9c$herkunft_FDZ)

shaw_eth_9c_dc <- asw.cluster::faultlines(data_conso_eth_9c, 
                                          group.par = c("IDTESTGROUP_FDZ"), 
                                          attr.type = c("nominal","nominal"), 
                                          rescale = "sd",
                                          method = "shaw",
                                          usesghomo = FALSE)


data_conso_gen_9c <- dplyr::select(data_dc, IDTESTGROUP_FDZ,TR_GESCHLECHT, ses_9c)%>%as.data.frame()
data_conso_gen_9c$TR_GESCHLECHT <- as.factor(data_conso_gen_9c$TR_GESCHLECHT)


shaw_gen_9c_dc<- asw.cluster::faultlines(data_conso_gen_9c, 
                                         group.par = c("IDTESTGROUP_FDZ"), 
                                         attr.type = c("nominal","nominal"), 
                                         rescale = "sd",
                                         method = "shaw",
                                         usesghomo = FALSE)

shaw_eth_9c_dc_long <- summary(shaw_eth_9c_dc)$long
  shaw_eth_9c_dc_indi <- select(shaw_eth_9c_dc_long, c(team,fl.value))
  shaw_eth_9c_dc_indi <- rename(shaw_eth_9c_dc_indi, IDTESTGROUP_FDZ=team, shaw_eth_9c=fl.value)

shaw_eth_9c_dc_cl <- shaw_eth_9c_dc_indi%>%
  group_by(IDTESTGROUP_FDZ) %>%
  summarise(shaw_eth_9c = unique(shaw_eth_9c))

data_dc<- left_join(data_dc, shaw_eth_9c_dc_cl)



shaw_gen_9c_dc_long <- summary(shaw_gen_9c_dc)$long
  shaw_gen_9c_dc_indi <- select(shaw_gen_9c_dc_long, c(team,fl.value))
  shaw_gen_9c_dc_indi <- rename(shaw_gen_9c_dc_indi, IDTESTGROUP_FDZ=team, shaw_gen_9c=fl.value)

shaw_gen_9c_dc_cl <- shaw_gen_9c_dc_indi%>%
  group_by(IDTESTGROUP_FDZ) %>%
  summarise(shaw_gen_9c = unique(shaw_gen_9c))

data_dc<- left_join(data_dc, shaw_gen_9c_dc_cl)


      #### 10 categories  -----
simps.10c <- simps(data_dc, "ses_10c")
  simps.10c<- rename(simps.10c, simps.10c=`diversity(simpdat3, index = "simpson")`)
  simps.10c <- rename(simps.10c, IDTESTGROUP_FDZ=names)
  simps.10c $IDTESTGROUP_FDZ <- as.character(simps.10c $IDTESTGROUP_FDZ)

  data_dc <- left_join(data_dc, simps.10c)


ratio_ses_10c_dc <- function(group) {
  
  individual_results <- list()
  
  help_columns <- paste0("help_", 1:34)
  
  for (i in 1:nrow(group)) {
    same_ties <- 0
    inter_ties <- 0
    total_ties <- group$total_tie[i]
    
    for (help_col in help_columns) {
      adjacency_value <- group[i, help_col]
      
      row_index <- which(colnames(group) == help_col) - 3  
      
      if (!is.na(adjacency_value) && adjacency_value == 1 && 
          !is.na(group$ses_10c[i]) && !is.na(group$ses_10c[row_index])) {
        
        if (group$ses_10c[i] == group$ses_10c[row_index]) {
          same_ties <- same_ties + 1
        }
      } 
    }
    
    inter_ties <- total_ties - same_ties

    ratio <- if (total_ties > 0) {
      if (is.na(inter_ties / total_ties)) {
        NA
      } else {
        inter_ties / total_ties
      }
    } else {
      NA
    }
    
    individual_results[[i]] <- list(IDTESTGROUP_FDZ = group$IDTESTGROUP_FDZ[i], 
                                    IDSTUD_FDZ = group$IDSTUD_FDZ[i], 
                                    ratio_10c_dc = ratio, 
                                    inter_tie_10c_dc = inter_ties, 
                                    same_tie_10c_dc = same_ties)
  }
  
  return(individual_results)
}


ratio_ses10cat<- lapply(network_all_cl, ratio_ses_10c_dc)
  ratio_ses10cat<<- lapply(ratio_ses10cat, bind_rows)%>%
  bind_rows(.)

  ratio_ses10cat$IDSTUD_FDZ <- as.character(ratio_ses10cat$IDSTUD_FDZ)
  ratio_ses10cat$IDTESTGROUP_FDZ <- as.character(ratio_ses10cat$IDTESTGROUP_FDZ)

  data_dc<- left_join(data_dc, ratio_ses10cat, by=c("IDSTUD_FDZ", "IDTESTGROUP_FDZ"))

data_dc$ses_10c <- as.factor(data_dc$ses_10c)
data_dc<- data_dc %>%
  group_by(IDTESTGROUP_FDZ, ses_10c) %>%
  mutate(same_hisei_10cat = n() - 1) %>%
  ungroup()


data_conso_eth_10c <- dplyr::select(data_dc, IDTESTGROUP_FDZ,herkunft_FDZ, ses_10c)%>%as.data.frame()
  data_conso_eth_10c$herkunft_FDZ <- as.factor(data_conso_eth_10c$herkunft_FDZ)

shaw_eth_10c_dc <- asw.cluster::faultlines(data_conso_eth_10c, 
                                           group.par = c("IDTESTGROUP_FDZ"), 
                                           attr.type = c("nominal","nominal"), 
                                           rescale = "sd",
                                           method = "shaw",
                                           usesghomo = FALSE)


data_conso_gen_10c <- dplyr::select(data_dc, IDTESTGROUP_FDZ,TR_GESCHLECHT, ses_10c)%>%as.data.frame()
  data_conso_gen_10c$TR_GESCHLECHT <- as.factor(data_conso_gen_10c$TR_GESCHLECHT)


shaw_gen_10c_dc<- asw.cluster::faultlines(data_conso_gen_10c, 
                                          group.par = c("IDTESTGROUP_FDZ"), 
                                          attr.type = c("nominal","nominal"), 
                                          rescale = "sd",
                                          method = "shaw",
                                          usesghomo = FALSE)


shaw_eth_10c_dc_long <- summary(shaw_eth_10c_dc)$long
  shaw_eth_10c_dc_indi <- select(shaw_eth_10c_dc_long, c(team,fl.value))
  shaw_eth_10c_dc_indi <- rename(shaw_eth_10c_dc_indi, IDTESTGROUP_FDZ=team, shaw_eth_10c=fl.value)

shaw_eth_10c_dc_cl <- shaw_eth_10c_dc_indi%>%
  group_by(IDTESTGROUP_FDZ) %>%
  summarise(shaw_eth_10c = unique(shaw_eth_10c))

data_dc<- left_join(data_dc, shaw_eth_10c_dc_cl)

shaw_gen_10c_dc_long <- summary(shaw_gen_10c_dc)$long
  shaw_gen_10c_dc_indi <- select(shaw_gen_10c_dc_long, c(team,fl.value))
  shaw_gen_10c_dc_indi <- rename(shaw_gen_10c_dc_indi, IDTESTGROUP_FDZ=team, shaw_gen_10c=fl.value)

shaw_gen_10c_dc_cl <- shaw_gen_10c_dc_indi%>%
  group_by(IDTESTGROUP_FDZ) %>%
  summarise(shaw_gen_10c = unique(shaw_gen_10c))

  data_dc<- left_join(data_dc, shaw_gen_10c_dc_cl)


    ### Regression ------------

      #### Classroom social heterogeneity and consolidation predicting proportion of cross-social-origin ties  --------------

  #4c
div_ratio_4c_dc<-lmer(ratio_4c_dc~simps.4c+
                        TR_GESCHLECHT+immi+ses_4c+mean_ses+size_orig+same_hisei_4c+
                        (1|IDTESTGROUP_FDZ),
                      data=data_dc_st, REML=T)

div_ratio_4c_ethmod_dc <- lmer(ratio_4c_dc~simps.4c*shaw_eth_dc+
                                 size_orig+TR_GESCHLECHT+immi+ses_4c+mean_ses+same_hisei_4c+
                                 (1|IDTESTGROUP_FDZ),
                               data=data_dc_st, REML=T)


div_ratio_4c_genmod_dc <- lmer(ratio_4c_dc~simps.4c*shaw_gen_dc+
                                 size_orig+TR_GESCHLECHT+immi+ses_4c+mean_ses+same_hisei_4c+
                                 (1|IDTESTGROUP_FDZ),
                               data=data_dc_st, REML=T)

  #5c 

div_ratio_5c_dc<-lmer(ratio_5c_dc~simps.5c+
                        TR_GESCHLECHT+immi+ses_5c+mean_ses+size_orig+same_hisei_5cat+
                        (1|IDTESTGROUP_FDZ),
                      data=data_dc_st, REML=T)

div_ratio_5c_ethmod_dc <- lmer(ratio_5c_dc~simps.5c*shaw_eth_dc+
                                 size_orig+TR_GESCHLECHT+immi+ses_5c+mean_ses+same_hisei_5cat+
                                 (1|IDTESTGROUP_FDZ),
                               data=data_dc_st, REML=T)


div_ratio_5c_genmod_dc <- lmer(ratio_5c_dc~simps.5c*shaw_gen_dc+
                                 size_orig+TR_GESCHLECHT+immi+ses_5c+mean_ses+same_hisei_5cat+
                                 (1|IDTESTGROUP_FDZ),
                               data=data_dc_st, REML=T)

  # 6c 
div_ratio_6c_dc<-lmer(ratio_6c_dc~simps.6c+
                        TR_GESCHLECHT+immi+ses_6c+mean_ses+size_orig+same_hisei_6cat+
                        (1|IDTESTGROUP_FDZ),
                      data=data_dc_st, REML=T)

div_ratio_6c_ethmod_dc <- lmer(ratio_6c_dc~simps.6c*shaw_eth_dc+
                                 size_orig+TR_GESCHLECHT+immi+ses_6c+mean_ses+same_hisei_6cat+
                                 (1|IDTESTGROUP_FDZ),
                               data=data_dc_st, REML=T)


div_ratio_6c_genmod_dc <- lmer(ratio_6c_dc~simps.6c*shaw_gen_dc+
                                 size_orig+TR_GESCHLECHT+immi+ses_6c+mean_ses+same_hisei_6cat+
                                 (1|IDTESTGROUP_FDZ),
                               data=data_dc_st, REML=T)


  # 7c 
div_ratio_7c_dc<-lmer(ratio_7c_dc~simps.7c+
                        TR_GESCHLECHT+immi+ses_7c+mean_ses+size_orig+same_hisei_7cat+
                        (1|IDTESTGROUP_FDZ),
                      data=data_dc_st, REML=T)

div_ratio_7c_ethmod_dc <- lmer(ratio_7c_dc~simps.7c*shaw_eth_dc+
                                 size_orig+TR_GESCHLECHT+immi+ses_7c+mean_ses+same_hisei_7cat+
                                 (1|IDTESTGROUP_FDZ),
                               data=data_dc_st, REML=T)


div_ratio_7c_genmod_dc <- lmer(ratio_7c_dc~simps.7c*shaw_gen_dc+
                                 size_orig+TR_GESCHLECHT+immi+ses_7c+mean_ses+same_hisei_7cat+
                                 (1|IDTESTGROUP_FDZ),
                               data=data_dc_st, REML=T)

  # 8c 
div_ratio_8c_dc<-lmer(ratio_8c_dc~simps.8c+
                        TR_GESCHLECHT+immi+ses_8c+mean_ses+size_orig+same_hisei_8cat+
                        (1|IDTESTGROUP_FDZ),
                      data=data_dc_st, REML=T)

div_ratio_8c_ethmod_dc <- lmer(ratio_8c_dc~simps.8c*shaw_eth_dc+
                                 size_orig+TR_GESCHLECHT+immi+ses_8c+mean_ses+same_hisei_8cat+
                                 (1|IDTESTGROUP_FDZ),
                               data=data_dc_st, REML=T)


div_ratio_8c_genmod_dc <- lmer(ratio_8c_dc~simps.8c*shaw_gen_dc+
                                 size_orig+TR_GESCHLECHT+immi+ses_8c+mean_ses+same_hisei_8cat+
                                 (1|IDTESTGROUP_FDZ),
                               data=data_dc_st, REML=T)
  
  # 9c 
div_ratio_9c_dc<-lmer(ratio_9c_dc~simps.9c+
                        TR_GESCHLECHT+immi+ses_9c+mean_ses+size_orig+same_hisei_9cat+
                        (1|IDTESTGROUP_FDZ),
                      data=data_dc_st, REML=T)

div_ratio_9c_ethmod_dc <- lmer(ratio_9c_dc~simps.9c*shaw_eth_dc+
                                 size_orig+TR_GESCHLECHT+immi+ses_9c+mean_ses+same_hisei_9cat+
                                 (1|IDTESTGROUP_FDZ),
                               data=data_dc_st, REML=T)


div_ratio_9c_genmod_dc <- lmer(ratio_9c_dc~simps.9c*shaw_gen_dc+
                                 size_orig+TR_GESCHLECHT+immi+ses_9c+mean_ses+same_hisei_9cat+
                                 (1|IDTESTGROUP_FDZ),
                               data=data_dc_st, REML=T)
  #10c 
div_ratio_10c_dc<-lmer(ratio_10c_dc~simps.10c+
                        TR_GESCHLECHT+immi+ses_10c+mean_ses+size_orig+same_hisei_10cat+
                        (1|IDTESTGROUP_FDZ),
                      data=data_dc_st, REML=T)

div_ratio_10c_ethmod_dc <- lmer(ratio_10c_dc~simps.10c*shaw_eth_dc+
                                 size_orig+TR_GESCHLECHT+immi+ses_10c+mean_ses+same_hisei_10cat+
                                 (1|IDTESTGROUP_FDZ),
                               data=data_dc_st, REML=T)


div_ratio_10c_genmod_dc <- lmer(ratio_10c_dc~simps.10c*shaw_gen_dc+
                                 size_orig+TR_GESCHLECHT+immi+ses_10c+mean_ses+same_hisei_10cat+
                                 (1|IDTESTGROUP_FDZ),
                               data=data_dc_st, REML=T)


    #### Proportion of cross-social-origin ties predicting achievement inequality --------------

ratio_ineq_4c_dc<- lmer(TR_NOTE_MAT_recode~ses_4c*ratio_4c_dc+
                          immi+TR_GESCHLECHT+gym+mean_ses+size_orig+
                          (ses_4c| IDTESTGROUP_FDZ), 
                        data = data_dc_st,REML=T)

ratio_ineq_5c_dc<- lmer(TR_NOTE_MAT_recode~ses_5c*ratio_5c_dc+
                          immi+TR_GESCHLECHT+gym+mean_ses+size_orig+
                          (ses_5c| IDTESTGROUP_FDZ), 
                        data = data_dc_st,REML=T)

ratio_ineq_6c_dc<- lmer(TR_NOTE_MAT_recode~ses_6c*ratio_6c_dc+
                          immi+TR_GESCHLECHT+gym+mean_ses+size_orig+
                          (ses_6c+TR_GESCHLECHT+immi| IDTESTGROUP_FDZ), 
                        data = data_dc_st,REML=T)

ratio_ineq_7c_dc<- lmer(TR_NOTE_MAT_recode~ses_7c*ratio_7c_dc+
                          immi+TR_GESCHLECHT+gym+mean_ses+size_orig+
                          (ses_7c| IDTESTGROUP_FDZ), 
                        data = data_dc_st,REML=T)

ratio_ineq_8c_dc<- lmer(TR_NOTE_MAT_recode~ses_8c*ratio_8c_dc+
                          immi+TR_GESCHLECHT+gym+mean_ses+size_orig+
                          (ses_8c+TR_GESCHLECHT+immi| IDTESTGROUP_FDZ), 
                        data = data_dc_st,REML=T)

ratio_ineq_9c_dc<- lmer(TR_NOTE_MAT_recode~ses_9c*ratio_9c_dc+
                          immi+TR_GESCHLECHT+gym+mean_ses+size_orig+
                          (ses_9c| IDTESTGROUP_FDZ), 
                        data = data_dc_st,REML=T)

ratio_ineq_10c_dc<- lmer(TR_NOTE_MAT_recode~ses_10c*ratio_10c_dc+
                           immi+TR_GESCHLECHT+gym+mean_ses+size_orig+
                           (ses_10c| IDTESTGROUP_FDZ), 
                         data = data_dc_st,REML=T)
    

      #### Direct associations ------------
div_ineq_4c_dc <-lmer(TR_NOTE_MAT_recode~ses_4c*simps.4c+
                        immi+TR_GESCHLECHT+gym+mean_ses+size_orig+
                        (ses_4c|IDTESTGROUP_FDZ),
                      data=data_dc_st, REML=T)

div_ineq_5c_dc <-lmer(TR_NOTE_MAT_recode~ses_5c*simps.5c+
                        immi+TR_GESCHLECHT+gym+mean_ses+size_orig+
                        (ses_5c|IDTESTGROUP_FDZ),
                      data=data_dc_st, REML=T)

div_ineq_6c_dc <-lmer(TR_NOTE_MAT_recode~ses_6c*simps.6c+
                        immi+TR_GESCHLECHT+gym+mean_ses+size_orig+
                        (ses_6c|IDTESTGROUP_FDZ),
                      data=data_dc_st, REML=T)

div_ineq_7c_dc <-lmer(TR_NOTE_MAT_recode~ses_7c*simps.7c+
                        immi+TR_GESCHLECHT+gym+mean_ses+size_orig+
                        (ses_7c|IDTESTGROUP_FDZ),
                      data=data_dc_st, REML=T)


div_ineq_8c_dc <-lmer(TR_NOTE_MAT_recode~ses_8c*simps.8c+
                        immi+TR_GESCHLECHT+gym+mean_ses+size_orig+
                        (ses_8c|IDTESTGROUP_FDZ),
                      data=data_dc_st, REML=T)

div_ineq_9c_dc <-lmer(TR_NOTE_MAT_recode~ses_9c*simps.9c+
                        immi+TR_GESCHLECHT+gym+mean_ses+size_orig+
                        (ses_9c|IDTESTGROUP_FDZ),
                      data=data_dc_st, REML=T)

div_ineq_10c_dc <-lmer(TR_NOTE_MAT_recode~ses_10c*simps.10c+
                         immi+TR_GESCHLECHT+gym+mean_ses+size_orig+
                         (ses_10c|IDTESTGROUP_FDZ),
                       data=data_dc_st, REML=T)

    ### AICBIC figure -------------
      #### Figure A1 ----------------
rm(models)
models <- base::list(M_3cat = div_ratio_3c_dc,
                     M_4cat = div_ratio_4c_dc, 
                     M_5cat = div_ratio_5c_dc,
                     M_6cat = div_ratio_6c_dc,
                     M_7cat = div_ratio_7c_dc,
                     M_8cat = div_ratio_8c_dc,
                     M_9cat = div_ratio_9c_dc,
                     M_10cat = div_ratio_10c_dc)


extract_model_info <- function(model, model_name) {
  tidy_data <- broom.mixed::tidy(model) %>%
    filter(effect == "fixed") %>%
    mutate(model = model_name)
  
  aic_value <- AIC(model)
  bic_value <- BIC(model)
  
  tidy_data <- tidy_data %>%
    mutate(AIC = aic_value, BIC = bic_value)
  
  return(tidy_data)
}


model_info <- purrr::map2_df(models, names(models), extract_model_info)

iv_names <- c("simps.3c", "simps.4c", "simps.5c", 
              "simps.6c", "simps.7c","simps.8c","simps.9c", "simps.10c")

coefficients_data <- model_info %>%
  filter(term %in% iv_names) %>%
  rename(coefficient = term)

model_names <- c("3", "4", "5", "6", "7", "8", "9", "10")
coefficients_data <- coefficients_data %>%
  mutate(model_numeric = as.numeric(factor(model, levels = names(models))),
         model_label = factor(model_numeric, labels = model_names))

min_aic <- min(coefficients_data$AIC, na.rm = TRUE)
min_bic <- min(coefficients_data$BIC, na.rm = TRUE)

library(ggplot2)
plot <- ggplot() +
  geom_point(data = coefficients_data, aes(x = model_label, y = estimate),
             position = position_dodge(width = 0.4), size = 3, shape = 16) +
  geom_errorbar(data = coefficients_data, aes(x = model_label, ymin = estimate - std.error, ymax = estimate + std.error),
                width = 0.2, position = position_dodge(width = 0.4)) +
  geom_text(data = coefficients_data, aes(x = model_label, y = max(estimate + std.error) * 1.05, label = paste0("AIC: ", round(AIC, 2)),
                                          color = ifelse(AIC == min_aic, "red", "black")),
            position = position_dodge(width = 0.4), vjust = -1, size = 2.7 ) +
  geom_text(data = coefficients_data, aes(x = model_label, y = max(estimate + std.error) * 1.2, label = paste0("BIC: ", round(BIC, 2)),
                                          color = ifelse(BIC == min_bic, "red", "black")),
            position = position_dodge(width = 0.4), vjust = -1, size = 2.7) +
  scale_color_manual(values = c("red" = "red", "black" = "black"), guide = "none") +
  # Axis labels
  labs(
    x = "Number of HISEI categories",
    y = "Coefficient"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels
    legend.position = "none" 
  )

ggsave("Figure A1", plot, width = 10, height = 7, dpi = 600)


  ### Figure A2  ----------------
rm(models)
models <- list(
  M_3cat = ratio_ineq_3c_dc,
  M_4cat = ratio_ineq_4c_dc, 
  M_5cat = ratio_ineq_5c_dc,
  M_6cat = ratio_ineq_6c_dc,
  M_7cat = ratio_ineq_7c_dc,
  M_8cat = ratio_ineq_8c_dc,
  M_9cat = ratio_ineq_9c_dc,
  M_10cat = ratio_ineq_10c_dc)

extract_model_info <- function(model, model_name) {
  tidy_data <- broom.mixed::tidy(model) %>%
    filter(effect == "fixed") %>%
    mutate(model = model_name)
  
  aic_value <- AIC(model)
  bic_value <- BIC(model)
  
  tidy_data <- tidy_data %>%
    mutate(AIC = aic_value, BIC = bic_value)
  
  return(tidy_data)
}

model_info <- purrr::map2_df(models, names(models), extract_model_info)

iv_names <- model_info %>%
  filter(grepl("[:*]", term)) %>%
  pull(term) %>%
  unique()

coefficients_data <- model_info %>%
  filter(term %in% iv_names) %>%
  rename(coefficient = term)

model_names <- c("3", "4", "5", "6", "7", "8", "9", "10")
coefficients_data <- coefficients_data %>%
  mutate(model_numeric = as.numeric(factor(model, levels = names(models))),
         model_label = factor(model_numeric, labels = model_names))

min_aic <- min(coefficients_data$AIC, na.rm = TRUE)
min_bic <- min(coefficients_data$BIC, na.rm = TRUE)

plot <- ggplot() +
   geom_point(data = coefficients_data, aes(x = model_label, y = estimate),
             position = position_dodge(width = 0.4), size = 3, shape = 16) +
   geom_errorbar(data = coefficients_data, aes(x = model_label, ymin = estimate - std.error, ymax = estimate + std.error),
                width = 0.2, position = position_dodge(width = 0.4)) +
   geom_text(data = coefficients_data, aes(x = model_label, y = max(estimate + std.error) * 1.05, label = paste0("AIC: ", round(AIC, 2)),
                                          color = ifelse(AIC == min_aic, "red", "black")),
            position = position_dodge(width = 0.4), vjust = -1, size = 2.7 ) +
   geom_text(data = coefficients_data, aes(x = model_label, y = max(estimate + std.error) * 1.2, label = paste0("BIC: ", round(BIC, 2)),
                                          color = ifelse(BIC == min_bic, "red", "black")),
            position = position_dodge(width = 0.4), vjust = -1, size = 2.7) +
   scale_color_manual(values = c("red" = "red", "black" = "black"), guide = "none") +

  labs(
    x = "Number of HISEI categories",
    y = "Coefficient"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels
    legend.position = "none" 
  )

print(plot)

ggsave("Figure A2.jpg", plot, width = 10, height = 7, dpi = 600)

      #### Figure A3 -------------
models <- list(
  M_3cat = div_ineq_3c_dc,
  M_4cat = div_ineq_4c_dc, 
  M_5cat = div_ineq_5c_dc,
  M_6cat = div_ineq_6c_dc,
  M_7cat = div_ineq_7c_dc,
  M_8cat = div_ineq_8c_dc,
  M_9cat = div_ineq_9c_dc,
  M_10cat = div_ineq_10c_dc)

extract_model_info <- function(model, model_name) {
  tidy_data <- broom.mixed::tidy(model) %>%
    filter(effect == "fixed") %>%
    mutate(model = model_name)
  
  aic_value <- AIC(model)
  bic_value <- BIC(model)
  
  tidy_data <- tidy_data %>%
    mutate(AIC = aic_value, BIC = bic_value)
  
  return(tidy_data)
}

model_info <- purrr::map2_df(models, names(models), extract_model_info)

iv_names <- model_info %>%
  filter(grepl("[:*]", term)) %>%
  pull(term) %>%
  unique()

coefficients_data <- model_info %>%
  filter(term %in% iv_names) %>%
  rename(coefficient = term)

model_names <- c("3", "4", "5", "6", "7", "8", "9", "10")
coefficients_data <- coefficients_data %>%
  mutate(model_numeric = as.numeric(factor(model, levels = names(models))),
         model_label = factor(model_numeric, labels = model_names))

min_aic <- min(coefficients_data$AIC, na.rm = TRUE)
min_bic <- min(coefficients_data$BIC, na.rm = TRUE)

plot <- ggplot() +
  geom_point(data = coefficients_data, aes(x = model_label, y = estimate),
             position = position_dodge(width = 0.4), size = 3, shape = 16) +
  geom_errorbar(data = coefficients_data, aes(x = model_label, ymin = estimate - std.error, ymax = estimate + std.error),
                width = 0.2, position = position_dodge(width = 0.4)) +
  geom_text(data = coefficients_data, aes(x = model_label, y = max(estimate + std.error) * 1.05, label = paste0("AIC: ", round(AIC, 2)),
                                          color = ifelse(AIC == min_aic, "red", "black")),
            position = position_dodge(width = 0.4), vjust = -1, size = 2.7) +
  geom_text(data = coefficients_data, aes(x = model_label, y = max(estimate + std.error) * 1.2, label = paste0("BIC: ", round(BIC, 2)),
                                          color = ifelse(BIC == min_bic, "red", "black")),
            position = position_dodge(width = 0.4), vjust = -1, size = 2.7) +
  scale_color_manual(values = c("red" = "red", "black" = "black"), guide = "none") +
  labs(
    x = "Number of HISEI categories",
    y = "Coefficient"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels
    legend.position = "none" 
  )

print(plot)
ggsave("Figure A3.jpg", plot, width = 10, height = 7, dpi = 600)




  ## 7.2 Using HISCED replace 3-category HISEI --------------

# classroom heterogeneity and consolidation
simps.hisced <- simps(data_dc, "HISCED")
  simps.hisced<- rename(simps.hisced, simps.hisced=`diversity(simpdat3, index = "simpson")`)
  simps.hisced <- rename(simps.hisced, IDTESTGROUP_FDZ=names)
  simps.hisced $IDTESTGROUP_FDZ <- as.character(simps.hisced $IDTESTGROUP_FDZ)

  data_dc<- left_join(data_dc, simps.hisced)




data_conso_eth_hisced <- dplyr::select(data_dc, IDTESTGROUP_FDZ,herkunft_FDZ, HISCED)%>%as.data.frame()
  data_conso_eth_hisced$herkunft_FDZ <- as.factor(data_conso_eth_hisced$herkunft_FDZ)


library(asw.cluster)
shaw_eth_isced <- asw.cluster::faultlines(data_conso_eth_hisced, 
                                          group.par = c("IDTESTGROUP_FDZ"), 
                                          attr.type = c("nominal","nominal"), 
                                          rescale = "sd",
                                          method = "shaw",
                                          usesghomo = FALSE)

shaw_eth_isced_long <- summary(shaw_eth_isced )$long
  shaw_eth_isced_long <- select(shaw_eth_isced_long, c(team,fl.value))
  shaw_eth_isced_long <- rename(shaw_eth_isced_long, IDTESTGROUP_FDZ=team, shaw_eth_isced_dc=fl.value)

shaw_eth_isced_cl <- shaw_eth_isced_long %>%
  group_by(IDTESTGROUP_FDZ) %>%
  summarise(shaw_eth_isced_dc = unique(shaw_eth_isced_dc))

  data_dc<- left_join(data_dc, shaw_eth_isced_cl)


  
data_conso_gen_hisced <- dplyr::select(data_dc, IDTESTGROUP_FDZ,TR_GESCHLECHT, HISCED)%>%as.data.frame()
  data_conso_gen_hisced$TR_GESCHLECHT <- as.factor(data_conso_gen_hisced$TR_GESCHLECHT)



shaw_gen_isced <- asw.cluster::faultlines(data_conso_gen_hisced, 
                                          group.par = c("IDTESTGROUP_FDZ"), 
                                          attr.type = c("nominal","nominal"), 
                                          rescale = "sd",
                                          method = "shaw",
                                          usesghomo = FALSE)

shaw_gen_isced_long <- summary(shaw_gen_isced )$long
  shaw_gen_isced_long <- select(shaw_gen_isced_long, c(team,fl.value))
  shaw_gen_isced_long <- rename(shaw_gen_isced_long, IDTESTGROUP_FDZ=team, shaw_gen_isced_dc=fl.value)

shaw_gen_isced_cl <- shaw_gen_isced_long %>%
  group_by(IDTESTGROUP_FDZ) %>%
  summarise(shaw_gen_isced_dc = unique(shaw_gen_isced_dc))

  data_dc<- left_join(data_dc, shaw_gen_isced_cl)
  
  
# proportion of cross-social-origin help-seeking ties 
  
ratio_hisced_dc <- function(group) {
  individual_results <- list()
  
  help_columns <- paste0("help_", 1:34)
  
  for (i in 1:nrow(group)) {
    same_ties <- 0
    inter_ties <- 0
    total_ties <- group$total_tie[i]
    
    for (help_col in help_columns) {
      adjacency_value <- group[i, help_col]
      
      row_index <- which(colnames(group) == help_col) - 3  
      
      if (!is.na(adjacency_value) && adjacency_value == 1 && 
          !is.na(group$HISCED[i]) && !is.na(group$HISCED[row_index])) {
        
      if (group$HISCED[i] == group$HISCED[row_index]) {
          same_ties <- same_ties + 1
        }
      } 
    }
    
    inter_ties <- total_ties - same_ties
    
   ratio <- if (total_ties > 0) {
      if (is.na(inter_ties / total_ties)) {
        NA
      } else {
        inter_ties / total_ties
      }
    } else {
      NA
    }
    
    individual_results[[i]] <- list(IDTESTGROUP_FDZ = group$IDTESTGROUP_FDZ[i], 
                                    IDSTUD_FDZ = group$IDSTUD_FDZ[i], 
                                    ratio_hisced_dc = ratio, 
                                    inter_tie_hisced_dc = inter_ties, 
                                    same_tie_hisced_dc = same_ties)
  }
  
  return(individual_results)
}

  hisced <- select(dat_com, IDSTUD_FDZ, HISCED)
    hisced$IDSTUD_FDZ <- as.numeric(hisced$IDSTUD_FDZ )
    network_all <- left_join(network_all, hisced)
    network_all_cl <- split(network_all, network_all$IDTESTGROUP_FDZ)

ratio_hisced_dc<- lapply(network_all_cl,  ratio_hisced_dc)
ratio_hisced_dc<<- lapply(ratio_hisced_dc, bind_rows)%>%
  bind_rows(.)
    ratio_hisced_dc$IDSTUD_FDZ <- as.character(ratio_hisced_dc$IDSTUD_FDZ)
    ratio_hisced_dc$IDTESTGROUP_FDZ <- as.character(ratio_hisced_dc$IDTESTGROUP_FDZ)

  
  data_dc<- left_join(data_dc, ratio_hisced_dc, by=c("IDSTUD_FDZ", "IDTESTGROUP_FDZ"))

  
  data_dc<- data_dc%>%
    group_by(IDTESTGROUP_FDZ, HISCED) %>%
    mutate(same_hisced = n()-1) %>%
    ungroup()
  
  
    ### regression -----------
  data_dc_st <- standardise(data_dc)
     #### Classroom social heterogeneity and consolidation predicting proportion of cross-social-origin ties  --------------

div_ratio_hisced_dc<-lmer(ratio_hisced_dc~simps.hisced+
                              TR_GESCHLECHT+immi+HISCED+mean_ses+size_orig+same_hisced+
                              (1|IDTESTGROUP_FDZ),
                            data=data_dc_st, REML=T)
  
div_ratio_hisced_ethmod_dc <- lmer(ratio_hisced_dc~simps.hisced*shaw_eth_isced_dc+
                                       size_orig+TR_GESCHLECHT+immi+HISCED+mean_ses+same_hisced+
                                       (1|IDTESTGROUP_FDZ),
                                     data=data_dc_st, REML=T)
  
  
div_ratio_hisced_genmod_dc <- lmer(ratio_hisced_dc~simps.hisced*shaw_gen_isced_dc+
                                       size_orig+TR_GESCHLECHT+immi+HISCED+mean_ses+same_hisced+
                                       (1|IDTESTGROUP_FDZ),
                                     data=data_dc_st, REML=T)
  
      #### Proportion of cross-social-origin ties predicting achievement inequality --------------
ratio_ineq_hisced_dc<- lmer(TR_NOTE_MAT_recode~HISCED*ratio_hisced_dc+
                              immi+TR_GESCHLECHT+gym+mean_ses+size_orig+
                              (HISCED| IDTESTGROUP_FDZ), 
                            data = data_dc_st,REML=T)


      #### Direct associations --------------

div_ineq_hisced_dc<-lmer(TR_NOTE_MAT_recode~HISCED*simps.hisced+
                           immi+TR_GESCHLECHT+gym+mean_ses+size_orig+
                           (HISCED|IDTESTGROUP_FDZ),
                         data=data_dc_st, REML=T)

div_ineq_hisced_ethmod_dc<-lmer(TR_NOTE_MAT_recode ~ HISCED*simps.hisced * shaw_eth_isced_dc + 
                                  immi +size_orig + TR_GESCHLECHT +  gym + mean_ses + 
                                  (HISCED| IDTESTGROUP_FDZ),
                                data = data_dc_st, REML = TRUE)


div_ineq_hisced_genmod_dc <-lmer(TR_NOTE_MAT_recode ~HISCED*simps.hisced*shaw_gen_isced_dc+ 
                                   immi+size_orig+TR_GESCHLECHT+gym+mean_ses + 
                                   (HISCED| IDTESTGROUP_FDZ),
                                 data = data_dc_st, REML = TRUE)





  ## 7.3 Cross-social-origin help-seeking ties with and without friendship -------------

# combine friendship network data
var <- dat_com[,grep("friend_", names(dat_com))]
var <- rbind(names(var), var)
friend <- as.character(var[1,])

network_friend <- select(dat_com,IDTESTGROUP_FDZ, IDSTUD_FDZ,IDinClass,HISEI,friend)
network_friend$HISEI <- ifelse(is.na(network_friend $HISEI), data_imp_all$HISEI, network_friend$HISEI)
sum(is.na(network_friend $HISEI))

network_friend$ses_3c_dc <- cut(network_friend$HISEI,
                                breaks = c(11, 42.3, 61.8, 89),
                                right = TRUE,         
                                include.lowest = TRUE,  
                                labels = c("[11,42.3]", "(42.3,61.8]", "(61.8,89]"))

network_friend$IDSTUD_FDZ <- as.character(network_friend$IDSTUD_FDZ)
network_friend$IDTESTGROUP_FDZ <- as.character(network_friend$IDTESTGROUP_FDZ)

network_friend <- network_friend%>%
  mutate(total_tie=rowSums(network_friend[, 5:38],na.rm=T))


# define calculation function
ratio_help_friend_dc<- function(group) {
   individual_results <- vector("list", length = nrow(group))
  
  help_columns <- paste0("help_", 1:34)
  friend_columns<- paste0("friend_", 1:34)
  
  for (i in seq_len(nrow(group))) {
    total_ties <- group$total_tie[i]
    friend_inter_ses_ties <- 0
    non_friend_inter_ses_ties <- 0
    
    
  for (j in seq_along(help_columns)) {
      help_col <- help_columns[j]
      friend_col <- friend_columns[j]
      
       if (!is.na(group[[help_col]][i]) && group[[help_col]][i] == 1) {
       if (!is.na(group$ses_3c_dc[i]) && !is.na(group$ses_3c_dc[j]) && group$ses_3c_dc[i] != group$ses_3c_dc[j]) {
       if (!is.na(group[[friend_col]][i]) && group[[friend_col]][i] == 1) {
            friend_inter_ses_ties <- friend_inter_ses_ties + 1
          } else {
            non_friend_inter_ses_ties <- non_friend_inter_ses_ties + 1
          }
        }
      }
    }
    ratio_friend_inter_ses <- ifelse(total_ties > 0, friend_inter_ses_ties / total_ties, NA)
    ratio_non_friend_inter_ses <- ifelse(total_ties > 0, non_friend_inter_ses_ties / total_ties, NA)
    
    
    individual_results[[i]] <- list(
      IDTESTGROUP_FDZ = group$IDTESTGROUP_FDZ[i],
      IDSTUD_FDZ = group$IDSTUD_FDZ[i],
      ratio_help_friend_dc = ratio_friend_inter_ses,
      ratio_help_nonfriend_dc= ratio_non_friend_inter_ses,
      n_help_non_friend_dc=non_friend_inter_ses_ties,
      n_help_friend_dc=friend_inter_ses_ties
    )
  }
  
  return(individual_results)
}

network_help_friend <- select(dat_com, IDTESTGROUP_FDZ, IDSTUD_FDZ, IDinClass, HISEI, help, friend)

breaks <- quantile(data_rc_imp$HISEI, probs = seq(0, 1, length.out = 4), na.rm = TRUE)
  network_help_friend$ses_3c_dc<- cut(network_all$HISEI, breaks = breaks, include.lowest = T)
  network_help_friend$total_tie <- rowSums(network_help_friend[, 5:38], na.rm = TRUE)
  network_help_friend$total_friend_tie <- rowSums(network_help_friend[, 39:72], na.rm = TRUE)

network_help_friend_cl<- split(network_help_friend, network_help_friend$IDTESTGROUP_FDZ)
network_help_friend_cl<- lapply(network_help_friend_cl, function(df) {
  df[order(df$IDinClass), ]
})

ratio_help_friend<- lapply(network_help_friend_cl, ratio_help_friend_dc)
ratio_help_friend2<- lapply(ratio_help_friend, bind_rows)%>%
  bind_rows(.)

  ratio_help_friend2$IDTESTGROUP_FDZ <- as.character(ratio_help_friend2$IDTESTGROUP_FDZ)
  ratio_help_friend2$IDSTUD_FDZ <- as.character(ratio_help_friend2$IDSTUD_FDZ )

  data_dc <- left_join(data_dc, ratio_help_friend2, by=c("IDTESTGROUP_FDZ", "IDSTUD_FDZ"))


  ### regression -----------
  data_dc_st <- standardise(data_dc)

    #### Classroom social heterogeneity and consolidation predicting proportion of cross-social-origin ties ---------
div_ratio_nonfriend_dc<- lmer(ratio_help_nonfriend_dc~simps.3c+
                             size_orig+TR_GESCHLECHT+immi+same_hisei_3c+ses_3c+mean_ses+
                             (1|IDTESTGROUP_FDZ),data=data_dc_st, REML=T)

div_ratio_nonfriend_eth_dc<- lmer(ratio_help_nonfriend_dc~simps.3c*shaw_eth_dc+
                                 size_orig+TR_GESCHLECHT+immi+same_hisei_3c+ses_3c+mean_ses+
                                 (1|IDTESTGROUP_FDZ),data=data_dc_st, REML=T)

div_ratio_nonfriend_gen_dc<- lmer(ratio_help_nonfriend_dc~simps.3c*shaw_gen_dc+
                                 size_orig+TR_GESCHLECHT+immi+same_hisei_3c+ses_3c+mean_ses+
                                 (1|IDTESTGROUP_FDZ),data=data_dc_st, REML=T)

div_ratio_friend_dc<- lmer(ratio_help_friend_dc~simps.3c+
                          size_orig+TR_GESCHLECHT+immi+same_hisei_3c+ses_3c+mean_ses+
                          (1|IDTESTGROUP_FDZ),data=data_dc_st , REML=T)

div_ratio_friend_eth_dc<- lmer(ratio_help_friend_dc~simps.3c*shaw_eth_dc+
                              size_orig+TR_GESCHLECHT+immi+same_hisei_3c+ses_3c+mean_ses+
                              (1|IDTESTGROUP_FDZ),data=data_dc_st , REML=T)

div_ratio_friend_gen_dc<- lmer(ratio_help_friend_dc~simps.3c*shaw_gen_dc+
                              size_orig+TR_GESCHLECHT+immi+same_hisei_3c+ses_3c+mean_ses+
                              (1|IDTESTGROUP_FDZ),data=data_dc_st , REML=T)
      


    #### Proportion of cross-social-origin ties predicting achievement inequality ---------
ratio_ineq_nonfriend_dc <- lmer(TR_NOTE_MAT_recode~ses_3c*ratio_help_nonfriend_dc+simps.3c+
                                  immi+ size_orig + TR_GESCHLECHT+ gym + same_hisei_3c+mean_ses+
                                  (ses_3c| IDTESTGROUP_FDZ), 
                                data = data_dc_st , REML=T)

ratio_ineq_friend_dc <- lmer(TR_NOTE_MAT_recode~ses_3c*ratio_help_friend_dc+simps.3c+
                               immi+ size_orig+TR_GESCHLECHT+ gym + same_hisei_3c+mean_ses+
                               (ses_3c| IDTESTGROUP_FDZ), 
                             data = data_dc_st , REML=T)



      #### Direct associations--------------

div_ineq_nonfriend_dc<-lmer(TR_NOTE_MAT_recode~ses_3c*simps.3c+ses_3c*ratio_help_nonfriend_dc+
                              immi+size_orig+TR_GESCHLECHT+gym+mean_ses+same_hisei_3c+
                              (immi+TR_GESCHLECHT+ses_3c+same_hisei_3c|IDTESTGROUP_FDZ),
                            data=data_dc_st , REML=T)

div_ineq_friend_dc<-lmer(TR_NOTE_MAT_recode~ses_3c*simps.3c+ses_3c*ratio_help_friend_dc+
                           immi+size_orig+TR_GESCHLECHT+gym+mean_ses+same_hisei_3c+
                           (ses_3c|IDTESTGROUP_FDZ),
                         data=data_dc_st , REML=T)



  ## 7.4 Curvilinear --------------
div_ratio_3c_cov_dc<- lmer(ratio_3c_dc~simps.3c+I(simps.3c^2)+
                             immi+size_orig+TR_GESCHLECHT+same_hisei_3c+ses_3c+mean_ses+
                             (1|IDTESTGROUP_FDZ),data=data_dc_st, REML=T)

  X_star <- - b["simps.3c"] / (2 * b["I(simps.3c^2)"])
  X_star

  mu   <- mean(data_dc$simps.3c,   na.rm=TRUE)
  sigma<- sd( data_dc$simps.3c,   na.rm=TRUE)
  X_star_orig <- X_star * sigma + mu

  test <- unique(data_dc[, c("IDTESTGROUP_FDZ", "simps.3c")])
  n_under_0.2 <- sum(test$simps.3c < 0.2, na.rm = TRUE)
