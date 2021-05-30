
# Packages ----------------------------------------------------------------
library(webshot)
library(metafor)
library(tidyverse)
library(rmarkdown)
library(magrittr)
library(grid)
library(forestplot)
library(glue)
library(Gmisc)
library(Matrix)


# Cleaning data -----------------------------------------------------------
music <- read.csv("6 Extracted Data/outcomes.csv")

music$study <- gsub(",", " et al.,", music$study)
music$comparison <- gsub("vs", "vs.", music$comparison)

names(music)[2] <- "author"
names(music)[33] <- "mean_post1"
names(music)[34] <- "sd_post1"
names(music)[40] <- "mean_post2"
names(music)[41] <- "sd_post2"


music[music$outcome=='Quality of life - Overall', c(31,33,35,38,40,42)] <- 
  (-1 * music[music$outcome=='Quality of life - Overall', c(31,33,35,38,40,42)])

names(music)[28] <- 'int_time_weeks'
music$int_time_weeks <- ifelse(is.na(music$int_time_weeks), "NR", music$int_time_weeks)


# Functions ---------------------------------------------------------------
md <- function(database){
  
  db <- database
  db <- escalc(measure="SMD", 
               m1i=mean_post1, m2i = mean_post2, 
               sd1i  = sd_post1, sd2i = sd_post2, 
               n1i = n1, n2i= n2, data=database)
  db$vi <- ifelse(is.na(db$vi), 
                  ((db$mdul_post-db$mdll_post)/((2*abs(qt(0.05/2, db$total-1)))^2)), db$vi)
  db <- db[order(db$yi),]
  db <- summary(db)
  
  db$md <- paste(formatC(db$yi, format='f', digits =2)," ",
                 "(",formatC(db$ci.lb, format='f', digits =2),
                 ",",formatC(db$ci.ub, format='f', digits=2),")") 
  db$postmean1 <- paste(formatC(db$mean_post1, format='f', digits=1),'(', formatC(db$sd_post1, format='f', digits=1),')')
  db$postmean2 <- paste(formatC(db$mean_post2, format='f', digits=1),'(', formatC(db$sd_post2, format='f', digits=1),')')
  
  ma <- rma(yi, vi, measure='SMD', data=db, method='REML', test='knha')
  
  db$weights <- weights(ma)
  db$w <- paste(formatC(db$weights, format='f', digits = 0),'%')
  
  list(pre = db, ma = ma)
}

subgroup_interaction <- function(group1, group2, name1, name2){
  
  x <- data.frame(estimate = c(coef(group1$ma), coef(group2$ma)), 
                  stderror = c(group1$ma$se, group2$ma$se),
                  meta = c(name1, name2), 
                  tau2 = round(c(group1$ma$tau2, group2$ma$tau2),3))
  
  subgroup_test <- rma(estimate, sei=stderror, mods = ~ meta, 
                       method="FE", data=x, digits=3)
}


table_md <- function(analysis, nstudies, int, comp, outcome, col,
                     follow=FALSE){
  
  authors <- c("Author", analysis$pre$author, 
               paste("Overall for ", analysis$ma$k, "studies"), 
               paste("(Tau^2 = ", (formatC(analysis$ma$tau2, digits=2, format="f")), ", df = ", 
                     (analysis$ma$k - analysis$ma$p),
                     ", p ", (ifelse(analysis$ma$QEp < 0.001, 
                                     paste("< 0.001"),
                                     paste("= ", formatC(analysis$ma$QEp, digits=3, format="f")))),
                     "; ", "I^2", " = ", (formatC(analysis$ma$I2, digits=1, format="f")), "%)"))
  comparison <- c("Comparison", analysis$pre$comparison, NA, NA)
  fu         <- c("Intervention\nTime", analysis$pre$int_time_weeks, NA, NA)
  int_pop    <- c(paste(int), analysis$pre$n1, sum(analysis$pre$n1), NA)
  int_cont   <- c(paste(col), analysis$pre$postmean1, NA, NA)
  comp_pop   <- c(paste(comp), analysis$pre$n2,sum(analysis$pre$n2), NA)
  comp_cont  <- c(paste(col), analysis$pre$postmean2, NA, NA)
  smd        <- c(paste0("Standardized Mean","\n","Difference (95% CI)"),
                  analysis$pre$md, 
                  paste(formatC(analysis$ma$b, format='f', digits =2), 
                        " (",formatC(analysis$ma$ci.lb, format='f', digits=2),
                        ",", formatC(analysis$ma$ci.ub, format='f', digits=2), ")"), NA)
  weight     <- c("Weight (%)", analysis$pre$w, "100 %", NA)
  
  
  ifelse(follow==T,
         b <- cbind(authors, comparison, fu, int_pop, int_cont, comp_pop, comp_cont, smd, weight),
         b <- cbind(authors, comparison, int_pop, int_cont, comp_pop, comp_cont, smd, weight))
  
  ifelse(nstudies>1,
         b <- rbind(b[1,], NA, NA, b[2:(nrow(b)-1),], b[nrow(b),]),
         b <- rbind(b[1,], NA, NA, b[2:(nrow(b)-2),], NA))
  
  ifelse(nstudies > 1,
         (c <- structure(list(
           mean = c(rep(NA, 3), analysis$pre$yi, analysis$ma$b,NA),
           lower = c(rep(NA, 3), analysis$pre$ci.lb, analysis$ma$ci.lb, NA),
           upper = c(rep(NA, 3), analysis$pre$ci.ub, analysis$ma$ci.ub, NA)),
           .Names = c("mean", "lower", "upper"),
           row.names = c(NA, -1L*nrow(b)),
           class = "data.frame")),
         (c <- structure(list(
           mean = c(rep(NA, 3), analysis$pre$yi, NA),
           lower = c(rep(NA, 3), analysis$pre$ci.lb, NA),
           upper = c(rep(NA, 3), analysis$pre$ci.ub, NA)),
           .Names = c("mean", "lower", "upper"),
           row.names = c(NA, -1L*nrow(b)),
           class = "data.frame")))
  
  c <- as_tibble(c)
  
  list(b = b, c = c)
} 



plotmd_single_studies <- function(words, numbers,  
                                  xtick, sizes, bolding, aligning, fpPosition,
                                  lines) {
  (forestplot(words,
              graph.pos = fpPosition,
              zero = 0,
              numbers,
              new_page = TRUE,
              colgap = unit(5, "mm"),
              hrzl_lines = lines,
              lineheight=unit(0.7,'cm'),
              boxsize = sizes,
              line.margin = 2,
              is.summary = bolding,
              align = aligning,
              ci.vertices = TRUE,
              txt_gp = fpTxtGp(label =gpar (cex=0.9), 
                               ticks = gpar(cex = 0.9, fontface="bold"),
                               summary = gpar(cex = 0.9, fontsize=13),
                               xlab = gpar(cex=0.9)),
              xticks = xtick,
              xlog=FALSE,
              clip = c(0,1),
              grid = xtick,
              lwd.xaxis = 1,
              lwd.ci = 2.2,
              lwd.zero = 1.5,
              graphwidth = unit(10,"cm"),
              col=fpColors(box=c("black"),line="grey", zero = 'dodgerblue4', axes="grey20", summary='black')))
  
}


# Meta-analyses primary outcomes -----------------------------------------------------------

#Speech
speech <- subset(music, outcome=='Speech improvement')
speech1 <- subset(speech, sensitivity1==1)
speech1_md <- md(speech1)
tspeech1 <- table_md(analysis = speech1_md, nstudies = 5, 
                     int = "Arm 1 (n)", 
                     comp = "Arm 2 (n)",
                     col = "Mean (SD)",
                     follow = FALSE)

# comprehension_analysis
comprehension <- subset(music, outcome=='Comprehension')
comp1 <- subset(comprehension, sensitivity1==1)
comp1_md <- md(comp1)
tcomp1 <- table_md(analysis = comp1_md, nstudies = 5, 
                   int = "Arm 1 (n)", 
                   comp = "Arm 2 (n)",
                   col = "Mean (SD)",
                   follow = FALSE)
# mood_analysis
mood <- subset(music, outcome=='Mood')
mood1 <- subset(mood, sensitivity1==1)
mood1_md <- md(mood1)
tmood1 <- table_md(analysis = mood1_md, nstudies = 5, 
                   int = "Arm 1 (n)", 
                   comp = "Arm 2 (n)",
                   col = "Mean (SD)",
                   follow = FALSE)
# social_analysis
social <- subset(music, outcome=='Social abilities')
social1 <- subset(social, sensitivity1==1)
social1_md <- md(social1)
tsocial1 <- table_md(analysis = social1_md, nstudies = 1, 
                     int = "Arm 1 (n)", 
                     comp = "Arm 2 (n)",
                     col = "Mean (SD)",
                     follow = FALSE)
# qol_analysis
qol <- subset(music, outcome=='Quality of life - Overall')
qol1 <- subset(qol, sensitivity1==1)
qol1_md <- md(qol1)
tqol1 <- table_md(analysis = qol1_md, nstudies = 1, 
                  int = "Arm 1 (n)", 
                  comp = "Arm 2 (n)",
                  col = "Mean (SD)",
                  follow = FALSE)


# Figure - overall main outcomes ------------------------------------------
overall1 <- rbind(
  tspeech1$b[1:3,],
  c('Outcome: Speech', rep(NA, 7)),
  tspeech1$b[4:12,], NA, NA,
  c('Outcome: Comprehension', rep(NA, 7)),
  tcomp1$b[4:10,], NA, NA,
  c('Outcome: Mood', rep(NA, 7)),
  tmood1$b[4:7,], NA, NA,
  c('Outcome: Social Abilities', rep(NA, 7)),
  tsocial1$b[4,], NA, NA,
  c('Outcome: Overall Quality of Life', rep(NA, 7)),
  tqol1$b[4,], NA)

sub1 <- c(5:13, 17:23, 27:30, 34,38)
overall1[,1][sub1] <- paste("       ",overall1[,1][sub1])

sub2 <- c(4, 16, 26, 33, 37)
overall1[,1][sub2] <- paste("   ",overall1[,1][sub2])


c_overall1 <- rbind(tspeech1$c[1:3,],
                    NA,
                    tspeech1$c[4:12,], NA, NA,
                    NA,
                    tcomp1$c[4:10,], NA, NA,
                    NA,
                    tmood1$c[4:7,], NA, NA,
                    NA,
                    tsocial1$c[4,], NA, NA,
                    NA,
                    tqol1$c[4,], NA)

plot <- plotmd_single_studies(words = overall1,
                              numbers = c_overall1,
                              fpPosition = ncol(overall1)-1,
                              xtick = c(-7,-4,-2,0,2,4,7), 
                              sizes = c(rep(NA, 4),
                                        0.009*(speech1_md$pre$weights+50),
                                        1,rep(NA, 4),
                                        0.009*(comp1_md$pre$weights+50),
                                        1, rep(NA, 4),
                                        0.008*(mood1_md$pre$weights+50),
                                        0.6, rep(NA, 4),
                                        1, NA, NA, NA,
                                        1, NA),
                              bolding = c(T, rep(F, 2), 
                                          T, rep(F, length(speech1_md$pre$comparison)), T, T,
                                          F, F,
                                          T, rep(F, length(comp1_md$pre$comparison)), T, T,
                                          F, F,
                                          T, rep(F, length(mood1_md$pre$comparison)), T, T,
                                          F, F,
                                          T, F, F, F,
                                          T, F, F),
                              lines = 
                                list('3' = gpar (lwd=1, columns=c(1:(ncol(tspeech1$b)+1)), col="black"),
                                     '12' = gpar (lwd=1, columns=c(3,5), col="#404040"),
                                     '15' = gpar (lwd=1, columns=c(1:(ncol(tspeech1$b)+1)), col="grey"),
                                     '22' = gpar (lwd=1, columns=c(3,5), col="#404040"),
                                     '25' = gpar (lwd=1, columns=c(1:(ncol(tspeech1$b)+1)), col="grey"),
                                     '29' = gpar (lwd=1, columns=c(3,5), col="#404040"),
                                     '32' = gpar (lwd=1, columns=c(1:(ncol(tspeech1$b)+1)), col="grey"),
                                     '36' = gpar (lwd=1, columns=c(1:(ncol(tspeech1$b)+1)), col="grey")),
                              aligning = c('l','l', 'c','l','c','l','r','c' ))


# Figure for subgroup analyses --------------------------------------------

# Supplementary Figure 1: </b> Subgroup analysis for speech by stroke time 
speech_subacute <- subset(music, outcome=='Speech improvement' & 
                            stroketime=='Subacute' & 
                            sensitivity1==1)
speech_chronic <- subset(music, outcome=='Speech improvement' & 
                           stroketime=='Chronic' &
                           sensitivity1==1)
speech_subacute_md <- md(speech_subacute)
speech_chronic_md <- md(speech_chronic)
tspeech_subacute <- table_md(analysis = speech_subacute_md, nstudies = speech_subacute_md$ma$k.all, 
                             int = "Arm 1 (n)", 
                             comp = "Arm 2 (n)",
                             col = "Mean (SD)",
                             follow = FALSE)
tspeech_chronic <- table_md(analysis = speech_chronic_md, nstudies = speech_chronic_md$ma$k.all, 
                            int = "Arm 1 (n)", 
                            comp = "Arm 2 (n)",
                            col = "Mean (SD)",
                            follow = FALSE)
subgroup_speech <- subgroup_interaction(group1 = speech_subacute_md,
                                        group2 = speech_chronic_md,
                                        name1 = "subacute",
                                        name2 = "chronic")
sub_speech_words <- rbind(
  tspeech_subacute$b[1:3,],
  c('Subgroup: <6 months since stroke (subacute)', rep(NA, 7)),
  tspeech_subacute$b[4:7,], NA, NA,
  c('Subgroup: >6 months since stroke (chronic)', rep(NA, 7)),
  tspeech_chronic$b[4:8,], NA)
sub1 <- c(5:8, 12:16)
sub_speech_words[,1][sub1] <- paste("       ",sub_speech_words[,1][sub1])
sub2 <- c(4, 11)
sub_speech_words[,1][sub2] <- paste("   ",sub_speech_words[,1][sub2])

sub_speech_num <- rbind(
  tspeech_subacute$c[1:3,],
  NA,
  tspeech_subacute$c[4:7,], NA, NA,
  NA,
  tspeech_chronic$c[4:8,], NA)
plot <- plotmd_single_studies(words = sub_speech_words,
                              numbers = sub_speech_num,
                              fpPosition = ncol(overall1)-1,
                              xtick = c(-7,-4,-2,0,2,4,7), 
                              sizes = c(rep(NA, 4),
                                        0.006*(speech_subacute_md$pre$weights+50),
                                        1,rep(NA, 4),
                                        0.006*(speech_chronic_md$pre$weights+50), 
                                        1, NA, NA),
                              bolding = c(T, rep(F, 2), 
                                          T, rep(F, 2), T, T,
                                          F, F,
                                          T, rep(F,3), T, T,
                                          F),
                              lines = list('3' = gpar (lwd=1, columns=c(1:(ncol(tspeech1$b)+1)), col="black"),
                                           '7' = gpar (lwd=1, columns=c(3,5), col="#404040"),
                                           '10' = gpar (lwd=1, columns=c(1:(ncol(tspeech1$b)+1)), col="grey"),
                                           '15' = gpar (lwd=1, columns=c(3,5), col="#404040")),
                              aligning = c('l','l', 'c','l','c','l','r','c' ))


x <- unit(.148, 'npc')
y <- unit(.12, 'npc')
grid.text('Test for subgroup differences: df=1, p = 0.617, I2 = 0%', x, y, gp = gpar(fontsize=13, fontface = 1, fontfamily='sans'))

# Supplementary Figure 2: Subgroup analysis for speech by musical experience 
speech_mus_exp <- subset(music, outcome=='Speech improvement' & 
                           musicalexperience=='Patients with musical experience' &
                           sensitivity1==1)
speech_mus_noexp <- subset(music, outcome=='Speech improvement' & 
                             musicalexperience=='Patients without musical experience' &
                             sensitivity1==1)
speech_mus_exp_md <- md(speech_mus_exp)
speech_mus_noexp_md <- md(speech_mus_noexp)
tspeech_mus_exp <- table_md(analysis = speech_mus_exp_md, nstudies = speech_mus_exp_md$ma$k.all, 
                            int = "Arm 1 (n)", 
                            comp = "Arm 2 (n)",
                            col = "Mean (SD)",
                            follow = FALSE)
tspeech_mus_noexp <- table_md(analysis = speech_mus_noexp_md, nstudies = speech_mus_noexp_md$ma$k.all, 
                              int = "Arm 1 (n)", 
                              comp = "Arm 2 (n)",
                              col = "Mean (SD)",
                              follow = FALSE)
speech_subg_musexp <- subgroup_interaction(group1 = speech_mus_exp_md,
                                           group2 = speech_mus_noexp_md,
                                           name1 = "musical_exp",
                                           name2 = "no_musical_exp")

subs_peech_musexp_words <- rbind(
  tspeech_mus_exp$b[1:3,],
  c('Subgroup: Patients with musical experience', rep(NA, 7)),
  tspeech_mus_exp$b[4,], NA, NA,
  c('Subgroup: Patients without musical experience', rep(NA, 7)),
  tspeech_mus_noexp$b[4:11,], NA)

sub1 <- c(5, 9:16)
subs_peech_musexp_words[,1][sub1] <- paste("       ",subs_peech_musexp_words[,1][sub1])
sub2 <- c(4, 8)
subs_peech_musexp_words[,1][sub2] <- paste("   ",subs_peech_musexp_words[,1][sub2])

subs_peech_musexp_num <- rbind(
  tspeech_mus_exp$c[1:3,],
  NA,
  tspeech_mus_exp$c[4,], NA, NA,
  NA,
  tspeech_mus_noexp$c[4:11,], NA)

plot <- plotmd_single_studies(words = subs_peech_musexp_words,
                              numbers = subs_peech_musexp_num,
                              fpPosition = ncol(overall1)-1,
                              xtick = c(-7,-4,-2,0,2,4,7), 
                              sizes = c(rep(NA, 4),
                                        1,
                                        rep(NA, 3),
                                        0.006*(speech_mus_noexp_md$pre$weights+50), 
                                        1, NA, NA),
                              bolding = c(T, rep(F, 2), 
                                          T, rep(F, 3), 
                                          T, rep(F, 6), T, T, F),
                              lines = list('3' = gpar (lwd=1, columns=c(1:9), col="black"),
                                           '7' = gpar (lwd=1, columns=c(1:9), col="grey")),
                              aligning = c('l','l', 'c','l','c','l','r','c' ))
x <- unit(.148, 'npc')
y <- unit(.12, 'npc')
grid.text('Test for subgroup differences: df=1, p = 0.192, I2 = 0%', x, y, gp = gpar(fontsize=13, fontface = 1, fontfamily='sans'))

# Supplementary Figure 3:</b> Subgroup analysis for speech by interventions based on rythm 
speech_rythm <- subset(music, outcome=='Speech improvement' & 
                         rythmbasedtherapy=='Rhythm-based therapy' & 
                         sensitivity1==1)
speech_norhythm <- subset(music, outcome=='Speech improvement' & 
                            rythmbasedtherapy=='No rhythm-based therapy' &
                            sensitivity1==1)

speech_rythm_md <- md(speech_rythm)
speech_norhythm_md <- md(speech_norhythm)

tspeech_rythm <- table_md(analysis = speech_rythm_md, nstudies = speech_rythm_md$ma$k.all, 
                          int = "Arm 1 (n)", 
                          comp = "Arm 2 (n)",
                          col = "Mean (SD)",
                          follow = FALSE)
tspeech_norhythm <- table_md(analysis = speech_norhythm_md, nstudies = speech_norhythm_md$ma$k.all, 
                             int = "Arm 1 (n)", 
                             comp = "Arm 2 (n)",
                             col = "Mean (SD)",
                             follow = FALSE)

speech_subg_rhythm <- subgroup_interaction(group1 = speech_rythm_md,
                                           group2 = speech_norhythm_md,
                                           name1 = "rythm",
                                           name2 = "norhythm")

sub_speech_rythm_words <- rbind(
  tspeech_rythm$b[1:3,],
  c('Subgroup: Interventions based on rythm', rep(NA, 7)),
  tspeech_rythm$b[4:10,], NA, NA,
  c('Subgroup: Interventions not based on rythm', rep(NA, 7)),
  tspeech_norhythm$b[4:7,], NA)

sub1 <- c(5:11, 15:18)
sub_speech_rythm_words[,1][sub1] <- paste("       ",sub_speech_rythm_words[,1][sub1])

sub2 <- c(4, 14)
sub_speech_rythm_words[,1][sub2] <- paste("   ",sub_speech_rythm_words[,1][sub2])

sub_speech_rythm_num <- rbind(
  tspeech_rythm$c[1:3,],
  NA,
  tspeech_rythm$c[4:10,], NA, NA,
  NA,
  tspeech_norhythm$c[4:7,], NA)

plot <- plotmd_single_studies(words = sub_speech_rythm_words,
                              numbers = sub_speech_rythm_num,
                              fpPosition = ncol(overall1)-1,
                              xtick = c(-7,-4,-2,0,2,4,7), 
                              sizes = c(rep(NA, 4),
                                        0.006*(speech_rythm_md$pre$weights+50), 1,
                                        rep(NA, 4),
                                        0.006*(speech_norhythm_md$pre$weights+50), 0.6,
                                        NA, NA),
                              bolding = c(T, rep(F, 2), 
                                          T, rep(F, length(speech_rythm_md$pre$comparison)), T, T,  
                                          F, F,
                                          T, rep(F, length(speech_norhythm_md$pre$comparison)), T, T,
                                          F),
                              lines = list('3' = gpar (lwd=1, columns=c(1:9), col="black"),
                                           '13' = gpar (lwd=3, columns=c(1:9), col="grey")),
                              aligning = c('l','l', 'c','l','c','l','r','c' ))
# Supplementary Figure 4: Subgroup analysis for speech by interventions based in melodic intonation 
speech_mit <- subset(music, outcome=='Speech improvement' & 
                       mit=='Melodic intonation therapy' & 
                       sensitivity1==1)
speech_nomit <- subset(music, outcome=='Speech improvement' & 
                         mit=='No melodic intonation therapy' &
                         sensitivity1==1)

speech_mit_md <- md(speech_mit)
speech_nomit_md <- md(speech_nomit)

tspeech_mit <- table_md(analysis = speech_mit_md, nstudies = speech_mit_md$ma$k.all, 
                        int = "Arm 1 (n)", 
                        comp = "Arm 2 (n)",
                        col = "Mean (SD)",
                        follow = FALSE)
tspeech_nomit <- table_md(analysis = speech_nomit_md, nstudies = speech_nomit_md$ma$k.all, 
                          int = "Arm 1 (n)", 
                          comp = "Arm 2 (n)",
                          col = "Mean (SD)",
                          follow = FALSE)

speech_subg_mit <- subgroup_interaction(group1 = speech_mit_md,
                                        group2 = speech_nomit_md,
                                        name1 = "mit",
                                        name2 = "nomit")

sub_speech_mit_words <- rbind(
  tspeech_mit$b[1:3,],
  c('Subgroup: Interventions based on Melodic Intonation', rep(NA, 7)),
  tspeech_mit$b[4:9,], NA, NA,
  c('Subgroup: Interventions not based on Melodic Intonation', rep(NA, 7)),
  tspeech_nomit$b[4:8,], NA)

sub1 <- c(5:10, 14:18)
sub_speech_mit_words[,1][sub1] <- paste("       ",sub_speech_mit_words[,1][sub1])

sub2 <- c(4, 13)
sub_speech_mit_words[,1][sub2] <- paste("   ",sub_speech_mit_words[,1][sub2])

sub_speech_mit_num <- rbind(
  tspeech_mit$c[1:3,],
  NA,
  tspeech_mit$c[4:9,], NA, NA,
  NA,
  tspeech_nomit$c[4:8,], NA)

plot <- plotmd_single_studies(words = sub_speech_mit_words,
                              numbers = sub_speech_mit_num,
                              fpPosition = ncol(overall1)-1,
                              xtick = c(-7,-4,-2,0,2,4,7), 
                              sizes = c(rep(NA, 4),
                                        0.006*(speech_mit_md$pre$weights+50), 1,
                                        rep(NA, 4),
                                        0.006*(speech_nomit_md$pre$weights+50), 0.8,
                                        NA, NA),
                              bolding = c(T, rep(F, 2), 
                                          T, rep(F, length(speech_mit_md$pre$comparison)), T, T,  
                                          F, F,
                                          T, rep(F, length(speech_nomit_md$pre$comparison)), T, T,
                                          F),
                              lines = list('3' = gpar (lwd=1, columns=c(1:9), col="black"),
                                           '12' = gpar (lwd=1, columns=c(1:9), col="grey")),
                              aligning = c('l','l', 'c','l','c','l','r','c' ))

# Supplementary Figure 5: Subgroup analysis for comprehension by stroke time 
comp_subacute <- subset(music, outcome=='Comprehension' & 
                          stroketime=='Subacute' & 
                          sensitivity1==1)
comp_chronic <- subset(music, outcome=='Comprehension' & 
                         stroketime=='Chronic' &
                         sensitivity1==1)

comp_subacute_md <- md(comp_subacute)
comp_chronic_md <- md(comp_chronic)

tcomp_subacute <- table_md(analysis = comp_subacute_md, nstudies = comp_subacute_md$ma$k.all, 
                           int = "Arm 1 (n)", 
                           comp = "Arm 2 (n)",
                           col = "Mean (SD)",
                           follow = FALSE)
tcomp_chronic <- table_md(analysis = comp_chronic_md, nstudies = comp_chronic_md$ma$k.all, 
                          int = "Arm 1 (n)", 
                          comp = "Arm 2 (n)",
                          col = "Mean (SD)",
                          follow = FALSE)

sub_comp_stroke <- subgroup_interaction(group1 = comp_subacute_md,
                                        group2 = comp_chronic_md,
                                        name1 = "subacute",
                                        name2 = "chronic")

sub_comp_words <- rbind(
  tcomp_subacute$b[1:3,],
  c('Subgroup: <6 months since stroke (subacute)', rep(NA, 7)),
  tcomp_subacute$b[4,], NA, NA,
  c('Subgroup: >6 months since stroke (chronic)', rep(NA, 7)),
  tcomp_chronic$b[4:8,], NA)

sub1 <- c(5, 9:13)
sub_comp_words[,1][sub1] <- paste("       ",sub_comp_words[,1][sub1])

sub2 <- c(4, 8)
sub_comp_words[,1][sub2] <- paste("   ",sub_comp_words[,1][sub2])


sub_comp_num <- rbind(
  tcomp_subacute$c[1:3,],
  NA,
  tcomp_subacute$c[4,], NA, NA,
  NA,
  tcomp_chronic$c[4:8,], NA)

plot <- plotmd_single_studies(words = sub_comp_words,
                              numbers = sub_comp_num,
                              fpPosition = ncol(overall1)-1,
                              xtick = c(-7,-4,-2,0,2,4,7), 
                              sizes = c(rep(NA, 4),
                                        1,rep(NA, 3),
                                        0.006*(comp_chronic_md$pre$weights+50), 
                                        1, NA, NA),
                              bolding = c(T, rep(F, 2), 
                                          T, rep(F, 3), 
                                          T, rep(F, 3), T, T, F),
                              lines = list('3' = gpar (lwd=1, columns=c(1:(ncol(tcomp1$b)+1)), col="black"),
                                           '7' = gpar (lwd=1, columns=c(1:(ncol(tcomp1$b)+1)), col="grey")),
                              aligning = c('l','l', 'c','l','c','l','r','c' ))
# Supplementary Figure 6: Subgroup analysis for comprehension by musical experience 
comp_mus_exp <- subset(music, outcome=='Comprehension' & 
                         musicalexperience=='Patients with musical experience' &
                         sensitivity1==1)
comp_mus_noexp <- subset(music, outcome=='Comprehension' & 
                           musicalexperience=='Patients without musical experience' &
                           sensitivity1==1)
comp_mus_exp_md <- md(comp_mus_exp)
comp_mus_noexp_md <- md(comp_mus_noexp)
tcomp_mus_exp <- table_md(analysis = comp_mus_exp_md, nstudies = comp_mus_exp_md$ma$k.all, 
                          int = "Arm 1 (n)", 
                          comp = "Arm 2 (n)",
                          col = "Mean (SD)",
                          follow = FALSE)
tcomp_mus_noexp <- table_md(analysis = comp_mus_noexp_md, nstudies = comp_mus_noexp_md$ma$k.all, 
                            int = "Arm 1 (n)", 
                            comp = "Arm 2 (n)",
                            col = "Mean (SD)",
                            follow = FALSE)
comp_subg_musexp <- subgroup_interaction(group1 = comp_mus_exp_md,
                                         group2 = comp_mus_noexp_md,
                                         name1 = "musical_exp",
                                         name2 = "no_musical_exp")


subs_peech_musexp_words <- rbind(
  tcomp_mus_exp$b[1:3,],
  c('Subgroup: Patients with musical experience', rep(NA, 7)),
  tcomp_mus_exp$b[4,], NA, NA,
  c('Subgroup: Patients without musical experience', rep(NA, 7)),
  tcomp_mus_noexp$b[4:9,], NA)

sub1 <- c(5, 9:14)
subs_peech_musexp_words[,1][sub1] <- paste("       ",subs_peech_musexp_words[,1][sub1])

sub2 <- c(4, 8)
subs_peech_musexp_words[,1][sub2] <- paste("   ",subs_peech_musexp_words[,1][sub2])


subs_peech_musexp_num <- rbind(
  tcomp_mus_exp$c[1:3,],
  NA,
  tcomp_mus_exp$c[4,], NA, NA,
  NA,
  tcomp_mus_noexp$c[4:9,], NA)

plot <- plotmd_single_studies(words = subs_peech_musexp_words,
                              numbers = subs_peech_musexp_num,
                              fpPosition = ncol(overall1)-1,
                              xtick = c(-7,-4,-2,0,2,4,7), 
                              sizes = c(rep(NA, 4),
                                        1,
                                        rep(NA, 3),
                                        0.006*(comp_mus_noexp_md$pre$weights+50), 
                                        1, NA, NA),
                              bolding = c(T, rep(F, 2), 
                                          T, rep(F, 3), 
                                          T, rep(F, 4), T, T, F),
                              lines = list('3' = gpar (lwd=1, columns=c(1:9), col="black"),
                                           '7' = gpar (lwd=1, columns=c(1:9), col="grey")),
                              aligning = c('l','l', 'c','l','c','l','r','c' ))
# Supplementary Figure 7: Subgroup analysis for comprehension by interventions based on rythm 
comp_rythm <- subset(music, outcome=='Comprehension' & 
                       rythmbasedtherapy=='Rhythm-based therapy' & 
                       sensitivity1==1)
comp_norhythm <- subset(music, outcome=='Comprehension' & 
                          rythmbasedtherapy=='No rhythm-based therapy' &
                          sensitivity1==1)

comp_rythm_md <- md(comp_rythm)
comp_norhythm_md <- md(comp_norhythm)

tcomp_rythm <- table_md(analysis = comp_rythm_md, nstudies = comp_rythm_md$ma$k.all, 
                        int = "Arm 1 (n)", 
                        comp = "Arm 2 (n)",
                        col = "Mean (SD)",
                        follow = FALSE)
tcomp_norhythm <- table_md(analysis = comp_norhythm_md, nstudies = comp_norhythm_md$ma$k.all, 
                           int = "Arm 1 (n)", 
                           comp = "Arm 2 (n)",
                           col = "Mean (SD)",
                           follow = FALSE)

comp_subg_rhythm <- subgroup_interaction(group1 = comp_rythm_md,
                                         group2 = comp_norhythm_md,
                                         name1 = "rythm",
                                         name2 = "norhythm")

sub_comp_rythm_words <- rbind(
  tcomp_rythm$b[1:3,],
  c('Subgroup: Interventions based on rythm', rep(NA, 7)),
  tcomp_rythm$b[4:8,], NA, NA,
  c('Subgroup: Interventions not based on rythm', rep(NA, 7)),
  tcomp_norhythm$b[4:7,], NA)

sub1 <- c(5:9, 13:16)
sub_comp_rythm_words[,1][sub1] <- paste("       ",sub_comp_rythm_words[,1][sub1])

sub2 <- c(4, 12)
sub_comp_rythm_words[,1][sub2] <- paste("   ",sub_comp_rythm_words[,1][sub2])

sub_comp_rythm_num <- rbind(
  tcomp_rythm$c[1:3,],
  NA,
  tcomp_rythm$c[4:8,], NA, NA,
  NA,
  tcomp_norhythm$c[4:7,], NA)

plot <- plotmd_single_studies(words = sub_comp_rythm_words,
                              numbers = sub_comp_rythm_num,
                              fpPosition = ncol(overall1)-1,
                              xtick = c(-7.5,-4,-2,0,2,4,7.5), 
                              sizes = c(rep(NA, 4),
                                        0.006*(comp_rythm_md$pre$weights+50), 1,
                                        rep(NA, 4),
                                        0.006*(comp_norhythm_md$pre$weights+50), 0.5,
                                        NA, NA),
                              bolding = c(T, rep(F, 2), 
                                          T, rep(F, length(comp_rythm_md$pre$comparison)), T, T,  
                                          F, F,
                                          T, rep(F, length(comp_norhythm_md$pre$comparison)), T, T,
                                          F),
                              lines = list('3' = gpar (lwd=1, columns=c(1:9), col="black"),
                                           '11' = gpar (lwd=1, columns=c(1:9), col="grey")),
                              aligning = c('l','l', 'c','l','c','l','r','c' ))

# Supplementary Figure 8: Subgroup analysis for comprehension by interventions based on melodic intonation 
comp_mit <- subset(music, outcome=='Comprehension' & 
                     mit=='Melodic intonation therapy' & 
                     sensitivity1==1)
comp_nomit <- subset(music, outcome=='Comprehension' & 
                       mit=='No melodic intonation therapy' &
                       sensitivity1==1)

comp_mit_md <- md(comp_mit)
comp_nomit_md <- md(comp_nomit)

tcomp_mit <- table_md(analysis = comp_mit_md, nstudies = comp_mit_md$ma$k.all, 
                      int = "Arm 1 (n)", 
                      comp = "Arm 2 (n)",
                      col = "Mean (SD)",
                      follow = FALSE)
tcomp_nomit <- table_md(analysis = comp_nomit_md, nstudies = comp_nomit_md$ma$k.all, 
                        int = "Arm 1 (n)", 
                        comp = "Arm 2 (n)",
                        col = "Mean (SD)",
                        follow = FALSE)

comp_subg_mit <- subgroup_interaction(group1 = comp_mit_md,
                                      group2 = comp_nomit_md,
                                      name1 = "mit",
                                      name2 = "nomit")

sub_comp_mit_words <- rbind(
  tcomp_mit$b[1:3,],
  c('Subgroup: Interventions based on Melodic Intonation', rep(NA, 7)),
  tcomp_mit$b[4:7,], NA, NA,
  c('Subgroup: Interventions not based on Melodic Intonation', rep(NA, 7)),
  tcomp_nomit$b[4:8,], NA)

sub1 <- c(5:8, 12:16)
sub_comp_mit_words[,1][sub1] <- paste("       ",sub_comp_mit_words[,1][sub1])

sub2 <- c(4, 11)
sub_comp_mit_words[,1][sub2] <- paste("   ",sub_comp_mit_words[,1][sub2])

sub_comp_mit_num <- rbind(
  tcomp_mit$c[1:3,],
  NA,
  tcomp_mit$c[4:7,], NA, NA,
  NA,
  tcomp_nomit$c[4:8,], NA)

plot <- plotmd_single_studies(words = sub_comp_mit_words,
                              numbers = sub_comp_mit_num,
                              fpPosition = ncol(overall1)-1,
                              xtick = c(-7,-4,-2,0,2,4,7), 
                              sizes = c(rep(NA, 4),
                                        0.006*(comp_mit_md$pre$weights+50), 0.7,
                                        rep(NA, 4),
                                        0.006*(comp_nomit_md$pre$weights+50), 1,
                                        NA, NA),
                              bolding = c(T, rep(F, 2), 
                                          T, rep(F, length(comp_mit_md$pre$comparison)), T, T,  
                                          F, F,
                                          T, rep(F, length(comp_nomit_md$pre$comparison)), T, T,
                                          F),
                              lines = list('3' = gpar (lwd=1, columns=c(1:9), col="black"),
                                           '10' = gpar (lwd=1, columns=c(1:9), col="grey")),
                              aligning = c('l','l', 'c','l','c','l','r','c' ))

# Supplementary Figure 9: Subgroup analysis for mood by stroke time 
mood_mus_exp <- subset(music, outcome=='Mood' & 
                         musicalexperience=='Patients with musical experience' &
                         sensitivity1==1)
mood_mus_noexp <- subset(music, outcome=='Mood' & 
                           musicalexperience=='Patients without musical experience' &
                           sensitivity1==1)
mood_mus_exp_md <- md(mood_mus_exp)
mood_mus_noexp_md <- md(mood_mus_noexp)
tmood_mus_exp <- table_md(analysis = mood_mus_exp_md, nstudies = mood_mus_exp_md$ma$k.all, 
                          int = "Arm 1 (n)", 
                          comp = "Arm 2 (n)",
                          col = "Mean (SD)",
                          follow = FALSE)
tmood_mus_noexp <- table_md(analysis = mood_mus_noexp_md, nstudies = mood_mus_noexp_md$ma$k.all, 
                            int = "Arm 1 (n)", 
                            comp = "Arm 2 (n)",
                            col = "Mean (SD)",
                            follow = FALSE)
mood_subg_musexp <- subgroup_interaction(group1 = mood_mus_exp_md,
                                         group2 = mood_mus_noexp_md,
                                         name1 = "musical_exp",
                                         name2 = "no_musical_exp")


subs_peech_musexp_words <- rbind(
  tmood_mus_exp$b[1:3,],
  c('Subgroup: Patients with musical experience', rep(NA, 7)),
  tmood_mus_exp$b[4,], NA, NA,
  c('Subgroup: Patients without musical experience', rep(NA, 7)),
  tmood_mus_noexp$b[4,], NA)

sub1 <- c(5, 9)
subs_peech_musexp_words[,1][sub1] <- paste("       ",subs_peech_musexp_words[,1][sub1])

sub2 <- c(4, 8)
subs_peech_musexp_words[,1][sub2] <- paste("   ",subs_peech_musexp_words[,1][sub2])


subs_peech_musexp_num <- rbind(
  tmood_mus_exp$c[1:3,],
  NA,
  tmood_mus_exp$c[4,], NA, NA,
  NA,
  tmood_mus_noexp$c[4,], NA)

plot <- plotmd_single_studies(words = subs_peech_musexp_words,
                              numbers = subs_peech_musexp_num,
                              fpPosition = ncol(overall1)-1,
                              xtick = c(-7,-4,-2,0,2,4,7), 
                              sizes = c(rep(NA, 4),
                                        1,
                                        rep(NA, 3),
                                        1, NA),
                              bolding = c(T, rep(F, 2), 
                                          T, rep(F, 3), 
                                          T, rep(F, 2)),
                              lines = list('3' = gpar (lwd=1, columns=c(1:9), col="black"),
                                           '7' = gpar (lwd=1, columns=c(1:9), col="grey")),
                              aligning = c('l','l', 'c','l','c','l','r','c' ))
# Supplementary Figure 10: Subgroup analysis for mood by interventions based on rythm 
mood_rythm <- subset(music, outcome=='Mood' & 
                       rythmbasedtherapy=='Rhythm-based therapy' & 
                       sensitivity1==1)
mood_norhythm <- subset(music, outcome=='Mood' & 
                          rythmbasedtherapy=='No rhythm-based therapy' &
                          sensitivity1==1)

mood_rythm_md <- md(mood_rythm)
mood_norhythm_md <- md(mood_norhythm)

tmood_rythm <- table_md(analysis = mood_rythm_md, nstudies = mood_rythm_md$ma$k.all, 
                        int = "Arm 1 (n)", 
                        comp = "Arm 2 (n)",
                        col = "Mean (SD)",
                        follow = FALSE)
tmood_norhythm <- table_md(analysis = mood_norhythm_md, nstudies = mood_norhythm_md$ma$k.all, 
                           int = "Arm 1 (n)", 
                           comp = "Arm 2 (n)",
                           col = "Mean (SD)",
                           follow = FALSE)

mood_subg_rhythm <- subgroup_interaction(group1 = mood_rythm_md,
                                         group2 = mood_norhythm_md,
                                         name1 = "rythm",
                                         name2 = "norhythm")

sub_mood_rythm_words <- rbind(
  tmood_rythm$b[1:3,],
  c('Subgroup: Interventions based on rythm', rep(NA, 7)),
  tmood_rythm$b[4,], NA, NA,
  c('Subgroup: Interventions not based on rythm', rep(NA, 7)),
  tmood_norhythm$b[4,], NA)

sub1 <- c(5, 9)
sub_mood_rythm_words[,1][sub1] <- paste("       ",sub_mood_rythm_words[,1][sub1])

sub2 <- c(4, 8)
sub_mood_rythm_words[,1][sub2] <- paste("   ",sub_mood_rythm_words[,1][sub2])

sub_mood_rythm_num <- rbind(
  tmood_rythm$c[1:3,],
  NA,
  tmood_rythm$c[4,], NA, NA,
  NA,
  tmood_norhythm$c[4,], NA)

plot <- plotmd_single_studies(words = sub_mood_rythm_words,
                              numbers = sub_mood_rythm_num,
                              fpPosition = ncol(overall1)-1,
                              xtick = c(-7.5,-4,-2,0,2,4,7.5), 
                              sizes = c(rep(NA, 4), 1,rep(NA, 3),1,NA),
                              bolding = c(T, rep(F, 2), T, rep(F,3), T, rep(F,2)),
                              lines = list('3' = gpar (lwd=1, columns=c(1:9), col="black"),
                                           '7' = gpar (lwd=1, columns=c(1:9), col="grey")),
                              aligning = c('l','l', 'c','l','c','l','r','c' ))



# Meta analyses surrogate outcomes ----------------------------------------
# motorspeech_analysis
motorsp <- subset(music, surrogate=='Motor-speech skills')
motorsp1 <- subset(motorsp, sensitivity1==1)
motorsp1_md <- md(motorsp1)
tmotorsp1 <- table_md(analysis = motorsp1_md, nstudies = 1, 
                      int = "Arm 1 (n)", 
                      comp = "Arm 2 (n)",
                      col = "Mean (SD)",
                      follow = FALSE)
# naming_analysis
naming <- subset(music, surrogate=='Naming')
naming1 <- subset(naming, sensitivity1==1)
naming1_md <- md(naming1)
tnaming1 <- table_md(analysis = naming1_md, nstudies = naming1_md$ma$k, 
                     int = "Arm 1 (n)", 
                     comp = "Arm 2 (n)",
                     col = "Mean (SD)",
                     follow = FALSE)
# repetition_analysis
repetition <- subset(music, surrogate=='Repetition')
repetition1 <- subset(repetition, sensitivity1==1)
repetition1_md <- md(repetition1)
trepetition1 <- table_md(
  analysis = repetition1_md, 
  nstudies = repetition1_md$ma$k, 
  int = "Arm 1 (n)", 
  comp = "Arm 2 (n)",
  col = "Mean (SD)",
  follow = FALSE)
# spontaneous_speech_analysis
spontSpeech1 <- subset(music, surrogate=='Spontanous speech')
spontSpeech1_md <- md(spontSpeech1)
tspontSpeech1 <- table_md(
  analysis = spontSpeech1_md, 
  nstudies = spontSpeech1_md$ma$k, 
  int = "Arm 1 (n)", 
  comp = "Arm 2 (n)",
  col = "Mean (SD)",
  follow = FALSE)
# story_telling_analysis
storyTell1 <- subset(music, surrogate=='Story telling')
storyTell1_md <- md(storyTell1)
tstoryTell1 <- table_md(
  analysis = storyTell1_md, 
  nstudies = storyTell1_md$ma$k, 
  int = "Arm 1 (n)", 
  comp = "Arm 2 (n)",
  col = "Mean (SD)",
  follow = FALSE)
# reading_comprehension_analysis
readComp <- subset(music, surrogate=='Reading comprehension')
readComp1 <- subset(readComp, sensitivity1==1)
readComp1_md <- md(readComp1)
treadComp1 <- table_md(analysis = readComp1_md, nstudies = 1, 
                       int = "Arm 1 (n)", 
                       comp = "Arm 2 (n)",
                       col = "Mean (SD)",
                       follow = FALSE)
# auditory_comprehension_analysis
audComp1 <- subset(music, surrogate=='Auditory comprehension')
audComp1_md <- md(audComp1)
taudComp1 <- table_md(
  analysis = audComp1_md, 
  nstudies = audComp1_md$ma$k, 
  int = "Arm 1 (n)", 
  comp = "Arm 2 (n)",
  col = "Mean (SD)",
  follow = FALSE)
# emotional_stability_analysis
emotionStab1 <- subset(music, surrogate=='Emotional stability')
emotionStab1_md <- md(emotionStab1)
temotionStab1 <- table_md(analysis = emotionStab1_md, nstudies = 1, 
                          int = "Arm 1 (n)", 
                          comp = "Arm 2 (n)",
                          col = "Mean (SD)",
                          follow = FALSE)
# energy_analysis
energy1 <- subset(music, surrogate=='Energy')
energy1_md <- md(energy1)
tenergy1 <- table_md(analysis = energy1_md, nstudies = 1, 
                     int = "Arm 1 (n)", 
                     comp = "Arm 2 (n)",
                     col = "Mean (SD)",
                     follow = FALSE)
# qol_mental_analysis
mh1 <- subset(music, outcome=='Quality of life - Mental health')
mh1_md <- md(mh1)
tmh1 <- table_md(analysis = mh1_md, nstudies = 1, 
                 int = "Arm 1 (n)", 
                 comp = "Arm 2 (n)",
                 col = "Mean (SD)",
                 follow = FALSE)
# qol_physical_analysis
ph1 <- subset(music, outcome=='Quality of life - Physical health')
ph1_md <- md(ph1)
tph1 <- table_md(analysis = ph1_md, nstudies = 1, 
                 int = "Arm 1 (n)", 
                 comp = "Arm 2 (n)",
                 col = "Mean (SD)",
                 follow = FALSE)
# qol_perception_analysis
perception1 <- subset(music, outcome=='Quality of life - General health')
perception1_md <- md(perception1)
tperception1 <- table_md(analysis = perception1_md, nstudies = 1, 
                         int = "Arm 1 (n)", 
                         comp = "Arm 2 (n)",
                         col = "Mean (SD)",
                         follow = FALSE)


# Figure surrogate outcomes 1 ---------------------------------------------
overall2 <- rbind(
  tnaming1$b[1:3,],
  c('Outcome: Naming', rep(NA, 7)),
  tnaming1$b[4:10,], NA, NA,
  c('Outcome: Spontanous speech', rep(NA, 7)),
  tspontSpeech1$b[4:11,], NA, NA,
  c('Outcome: Story telling', rep(NA, 7)),
  tstoryTell1$b[4:8,], NA, NA,
  c('Outcome: Motor-speech skills', rep(NA, 7)),
  tmotorsp1$b[4,], NA)

sub3 <- c(5:11, 15:22,26:30,34)
overall2[,1][sub3] <- paste("       ",overall2[,1][sub3])

sub4 <- c(4, 14, 25, 33)
overall2[,1][sub4] <- paste("   ",overall2[,1][sub4])

c_overall2 <- rbind(
  tnaming1$c[1:3,],
  NA,
  tnaming1$c[4:10,], NA, NA,
  NA,
  tspontSpeech1$c[4:11,], NA, NA,
  NA,
  tstoryTell1$c[4:8,], NA, NA,
  NA,
  tmotorsp1$c[4,], NA)

plot <- plotmd_single_studies(words = overall2,
                              numbers = c_overall2,
                              fpPosition = ncol(overall2)-1,
                              xtick = c(-3,-2,0,2,3), 
                              sizes = c(rep(NA, 4),
                                        0.009*(naming1_md$pre$weights+50),
                                        1,rep(NA, 4),
                                        0.009*(spontSpeech1_md$pre$weights+50),
                                        1, rep(NA, 4),
                                        0.007*(storyTell1_md$pre$weights+50),
                                        1, rep(NA, 4), 
                                        1, NA),
                              bolding = c(T, rep(F, 2), 
                                          T, rep(F, length(naming1_md$pre$comparison)), T, T,
                                          F, F, 
                                          T, rep(F,length(spontSpeech1_md$pre$comparison)), T, T,
                                          F, F,
                                          T, rep(F, length(storyTell1_md$pre$comparison)), T, T,
                                          F, F, 
                                          T, rep(F, 2)),
                              lines = list('3' = gpar (lwd=1, columns=c(1:9), col="black"),
                                           '10' = gpar (lwd=1, columns=c(3,5), col="#404040"),
                                           '13' = gpar (lwd=1, columns=c(1:9), col="grey"),
                                           '21' = gpar (lwd=1, columns=c(3,5), col="#404040"),
                                           '24' = gpar (lwd=1, columns=c(1:9), col="grey"),
                                           
                                           
                                           '29' = gpar (lwd=1, columns=c(3,5), col="#404040"),
                                           '32' = gpar (lwd=1, columns=c(1:9), col="grey")),
                              aligning = c('l','l', 'c','l','c','l','r','c' ))



# Figure surrogate outcomes 2 ---------------------------------------------
overall3 <- rbind(
  taudComp1$b[1:3,],
  c('Outcome: Auditory comprehension', rep(NA, 7)),
  taudComp1$b[4:7,], NA, NA,
  c('Outcome: Reading comprehension', rep(NA, 7)),
  treadComp1$b[4,], NA, NA,
  c('Outcome: Emotional stability', rep(NA, 7)),
  temotionStab1$b[4,], NA, NA,
  c('Outcome: Energy', rep(NA, 7)),
  tenergy1$b[4,], NA, NA,
  c('Outcome: Mental health', rep(NA, 7)),
  tmh1$b[4,], NA, NA,
  c('Outcome: Physical health', rep(NA, 7)),
  tph1$b[4,], NA, NA,
  c('Outcome: Perception of health', rep(NA, 7)),
  tperception1$b[4,], NA)

sub5 <- c(5:8, 12, 16, 20, 24, 28, 32)
overall3[,1][sub5] <- paste("       ",overall3[,1][sub5])
sub6 <- c(4, 11, 15, 19, 23, 27, 31)
overall3[,1][sub6] <- paste("   ",overall3[,1][sub6])

c_overall3 <- rbind(
  taudComp1$c[1:3,],
  NA,
  taudComp1$c[4:7,], NA, NA,
  NA,
  treadComp1$c[4,], NA, NA,
  NA,
  temotionStab1$c[4,], NA, NA,
  NA,
  tenergy1$c[4,], NA, NA,
  NA,
  tmh1$c[4,], NA, NA,
  NA,
  tph1$c[4,], NA, NA,
  NA,
  tperception1$c[4,], NA)


plot <- plotmd_single_studies(words = overall3,
                              numbers = c_overall3,
                              fpPosition = ncol(overall3)-1,
                              xtick = c(-3,-2,0,2,3), 
                              sizes = c(rep(NA, 4),
                                        0.007*(audComp1_md$pre$weights+50),
                                        0.6, rep(NA, 4),
                                        1, rep(NA, 3),
                                        1, rep(NA, 3),
                                        1, rep(NA, 3),
                                        1, rep(NA, 3),
                                        1, rep(NA, 3),
                                        1, NA),
                              bolding = c(T, F, F, 
                                          T, rep(F, 2), T, T, F, F,
                                          T, rep(F, 3),
                                          T, rep(F, 3), 
                                          T, rep(F, 3),
                                          T, rep(F, 3),
                                          T, rep(F, 3),
                                          T, F, F),
                              lines = list('3' = gpar (lwd=1, columns=c(1:9), col="black"),
                                           '7' = gpar (lwd=1, columns=c(3,5), col="#404040"),
                                           '10' = gpar (lwd=1, columns=c(1:9), col="grey"),
                                           '14' = gpar (lwd=1, columns=c(1:9), col="grey"),
                                           '18' = gpar (lwd=1, columns=c(1:9), col="grey"),
                                           '22' = gpar (lwd=1, columns=c(1:9), col="grey"),
                                           '26' = gpar (lwd=1, columns=c(1:9), col="grey"),
                                           '30' = gpar (lwd=1, columns=c(1:9), col="grey")),
                              aligning = c('l','l', 'c','l','c','l','r','c' ))


# Supplementary figures sensitivity analysis----------------------------

# Outcome: Speech 
speech2 <- subset(speech, sensitivity2==2)
speech2_md <- md(speech2)
tspeech2 <- table_md(analysis = speech2_md, nstudies = 5, 
                     int = "Arm 1 (n)", 
                     comp = "Arm 2 (n)",
                     col = "Mean (SD)",
                     follow = FALSE)
plot <- plotmd_single_studies(words = tspeech2$b,
                              numbers = tspeech2$c,
                              fpPosition = ncol(tspeech2$b)-1,
                              xtick = c(-2,-1,0,1,2), 
                              sizes = c(rep(NA,3), 0.011*(speech2_md$pre$weights+50), 1, 1),
                              bolding = c(T, F, F,rep(F, length(speech2_md$pre$comparison)), T, T),
                              lines = list('3' = gpar (lwd=1, columns=c(1:9), col="black"),
                                           '11' = gpar (lwd=1, columns=c(3,5), col="#404040")),
                              aligning = c('l','l', rep('c', 4), 'l', 'l'))
# Outcome: Comprehension 
comp2 <- subset(comprehension, sensitivity2==2)
comp2_md <- md(comp2)
tcomp2 <- table_md(analysis = comp2_md, nstudies = 5, 
                   int = "Arm 1 (n)", 
                   comp = "Arm 2 (n)",
                   col = "Mean (SD)",
                   follow = FALSE)
plot <- plotmd_single_studies(words = tcomp2$b,
                              numbers = tcomp2$c,
                              fpPosition = ncol(tcomp2$b)-1,
                              xtick = c(-2,-1,0,1,2), 
                              sizes = 
                                c(rep(NA,3), 0.011*(comp2_md$pre$weights+50), 0.8,NA),
                              bolding = c(T, F, F, rep(F, length(comp2_md$pre$comparison)), T, T),
                              lines = list('3' = gpar (lwd=1, columns=c(1:9), col="black"),
                                           '9' = gpar (lwd=1, columns=c(3,5), col="#404040")),
                              aligning = c('l','l', rep('c', 4), 'l', 'l'))
# Outcome: Mood  
mood2 <- subset(mood, sensitivity2==2)
mood2_md <- md(mood2)
tmood2 <- table_md(analysis = mood2_md, nstudies = 5, 
                   int = "Arm 1 (n)", 
                   comp = "Arm 2 (n)",
                   col = "Mean (SD)",
                   follow = FALSE)
plot <- plotmd_single_studies(words = tmood2$b,
                              numbers = tmood2$c,
                              fpPosition = ncol(tmood2$b)-1,
                              xtick = c(-12, -6, 0, 6, 12), 
                              sizes = 
                                c(rep(NA,3), 0.008*(mood2_md$pre$weights+50), 0.5,NA),
                              bolding = c(T, F, F, rep(F, length(mood2_md$pre$comparison)), T, T),
                              lines = list('3' = gpar (lwd=1, columns=c(1:9), col="black"),
                                           '6' = gpar (lwd=1, columns=c(3,5), col="#404040")),
                              aligning = c('l','l', rep('c', 4), 'l', 'l'))
# Outcome: Overall quality of life 
qol2 <- subset(qol, sensitivity2==2)

qol2_md <- md(qol2)
tqol2 <- table_md(analysis = qol2_md, nstudies = 1, 
                  int = "Arm 1 (n)", 
                  comp = "Arm 2 (n)",
                  col = "Mean (SD)",
                  follow = FALSE)

plot <- plotmd_single_studies(words = tqol2$b,
                              numbers = tqol2$c,
                              fpPosition = ncol(tqol2$b)-1,
                              xtick = c(-3,-2,-1,0,1,2,3), 
                              sizes = 
                                c(rep(NA,3), 1, NA),
                              bolding = c(T, F, F, rep(F, 2), T, F),
                              lines = list('3' = gpar (lwd=1, columns=c(1:9), col="black")),
                              aligning = c('l', rep('c', 4), 'c', 'l', 'l'))
# Surrogate outcome: Motor speech 
motorsp2 <- subset(motorsp, sensitivity2==2)
motorsp2_md <- md(motorsp2)
tmotorsp2 <- table_md(analysis = motorsp2_md, nstudies = 1, 
                      int = "Arm 1 (n)", 
                      comp = "Arm 2 (n)",
                      col = "Mean (SD)",
                      follow = FALSE)
plot <- plotmd_single_studies(words = tmotorsp2$b,
                              numbers = tmotorsp2$c,
                              fpPosition = ncol(tmotorsp2$b),
                              xtick = c(-3,-2,-1,0,1,2,3), 
                              sizes = 
                                c(rep(NA,3), 1, NA),
                              bolding = c(T, F, F, rep(F, 2), T, F),
                              lines = list('3' = gpar (lwd=1, columns=c(1:9), col="black")),
                              aligning = c('l', rep('c', 4), 'c', 'l', 'l'))
# Surrogate outcome: Naming 
naming2 <- subset(naming, sensitivity2==2)
naming2_md <- md(naming2)
tnaming2 <- table_md(analysis = naming2_md, nstudies = naming2_md$ma$k, 
                     int = "Arm 1 (n)", 
                     comp = "Arm 2 (n)",
                     col = "Mean (SD)",
                     follow = FALSE)
plot <- plotmd_single_studies(words = tnaming2$b,
                              numbers = tnaming2$c,
                              fpPosition = ncol(tnaming2$b)-1,
                              xtick = c(-3,-2,-1,0,1,2,3), 
                              sizes = 
                                c(rep(NA,3), 0.007*(naming2_md$pre$weights+50), 0.9,NA),
                              bolding = c(T, rep(F, 2), rep(F, length(naming2_md$pre$weights)), T, T),
                              lines = list('3' = gpar (lwd=1, columns=c(1:9), col="black"),
                                           '9' = gpar (lwd=1, columns=c(3,5), col="#404040")),
                              aligning = c('l','l', rep('c', 4), 'c', 'l', 'l'))
# Surrogate outcome: Repetition 
repetition2 <- subset(repetition, sensitivity2==2)
repetition2_md <- md(repetition2)
trepetition2 <- table_md(analysis = repetition2_md, nstudies = repetition2_md$ma$k, 
                         int = "Arm 1 (n)", 
                         comp = "Arm 2 (n)",
                         col = "Mean (SD)",
                         follow = FALSE)
plot <- plotmd_single_studies(words = trepetition2$b,
                              numbers = trepetition2$c,
                              fpPosition = ncol(trepetition2$b)-1,
                              xtick = c(-3,-2,-1,0,1,2,3), 
                              sizes = 
                                c(rep(NA,3), 0.007*(repetition2_md$pre$weights+50),  0.9,NA),
                              bolding = c(T, rep(F, 2), rep(F, length(repetition2_md$pre$weights)), T, T),
                              lines = list('3' = gpar (lwd=1, columns=c(1:9), col="black"),
                                           '10' = gpar (lwd=1, columns=c(3,5), col="#404040")),
                              aligning = c('l','l', rep('c', 4), 'c', 'l', 'l'))
# Surrogate outcome: Reading comprehension 
readComp2 <- subset(readComp, sensitivity2==2)
readComp2_md <- md(readComp2)
treadComp2 <- table_md(analysis = readComp2_md, nstudies = 1, 
                       int = "Arm 1 (n)", 
                       comp = "Arm 2 (n)",
                       col = "Mean (SD)",
                       follow = FALSE)
plot <- plotmd_single_studies(words = treadComp2$b,
                              numbers = treadComp2$c,
                              fpPosition = ncol(treadComp2$b)-1,
                              xtick = c(-3,-2,-1,0,1,2,3), 
                              sizes = 
                                c(rep(NA,3), 1, NA),
                              bolding = c(T, rep(F, 2), F, F),
                              lines = list('3' = gpar (lwd=1, columns=c(1:9), col="black")),
                              aligning = c('l', rep('c', 4), 'c', 'l', 'l'))