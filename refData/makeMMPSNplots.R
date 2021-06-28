library('tidyverse')
library('RColorBrewer')
library("survival")
library("survminer")
library('grid')
surv_psn = readRDS("overall_survival_PSN_groups.RDS")
colnames(surv_psn) = c("censos" ,"ttcos","SubGroup")
fit <- survfit(Surv(ttcos, censos) ~ SubGroup, data = surv_psn)
subgroups <- unique(surv_psn$SubGroup)

# theme <- theme()
psnCurveList <- list()
for (group in subgroups){

  tmp <- surv_psn
  # move factor level for group of interest to very end so it's on top
  tmp$SubGroup <- factor(tmp$SubGroup, c(subgroups[subgroups != group], group))
  fit <- survfit(Surv(ttcos, censos) ~ SubGroup, data = tmp)
  plt <- ggsurvplot(fit, data = tmp, risk.table = FALSE, pval = FALSE, censor=FALSE,
                   palette = c(rep("gray80",length(subgroups)-1), "#00AEEF"),
                   ylab="Overall Survival", legend.title="SubGroup", pval.size = 18, legend="none",
                   ggtheme = theme_classic())

  psnCurveList[[group]] <- plt
}

saveRDS(psnCurveList, file='mmPSN_survivalCurvePlotList.rds')


# #1a
# pdf(file="survival_SubGroups_OS_sample_1a.pdf", width = 1.8, height =1.3)
# ## use these colors if need colored plot, the colors are in sequnce 1a,1b,1c,1d,2a,2b,2c,2d,2e,3a,3b,3c
# #ggsurvplot(fit, data = s, risk.table = TRUE, pval = TRUE,palette =)
# aa=ggsurvplot(
#   fit, data = surv_psn,
#   risk.table = FALSE, pval = FALSE,palette = c("black", "grey",  "grey","grey","grey",  "grey",  "grey","grey","grey", "grey","grey", "grey"),ylab="Overall Survival"
#   ,legend.title="SubGroup",pval.size = 18,size = 0.25,legend="none",censor.size=0.2,

#   ggtheme = theme_survminer(

#     font.x = c(5, "bold", "black"),
#     font.y = c(5, "bold", "black"),
#     font.legend = c(5, "bold"), 
#     font.tickslab = c(5, "bold", "black"),
#   )
# )
# #aa$plot <- aa$plot + theme(legend.key.height =unit(3, "line"),legend.key.size = unit(3,"line"))
# aa
# dev.off()


# #1b

# pdf(file="survival_SubGroups_OS_sample_1b.pdf", width = 1.8, height =1.3)
# ## use these colors if need colored plot, the colors are in sequnce 1a,1b,1c,1d,2a,2b,2c,2d,2e,3a,3b,3c
# #ggsurvplot(fit, data = s, risk.table = TRUE, pval = TRUE,palette = c( "#1B9E77", "aquamarine2",  "limegreen","#666666","#F46D43",  "darkred",  "#FA9FB5","#E6AB02","magenta1", "purple3",  "deepskyblue",  "gray4"))
# aa=ggsurvplot(
#   fit, data = surv_psn,
#   risk.table = FALSE, pval = FALSE,palette = c("grey", "black",  "grey","grey","grey",  "grey",  "grey","grey","grey", "grey","grey", "grey"),ylab="Overall Survival"
#   ,legend.title="SubGroup",pval.size = 18,size = 0.25,legend="none",censor.size=0.2,

#   ggtheme = theme_survminer(

#     font.x = c(5, "bold", "black"),
#     font.y = c(5, "bold", "black"),
#     font.legend = c(5, "bold"), 
#     font.tickslab = c(5, "bold", "black"),
#   )
# )
# #aa$plot <- aa$plot + theme(legend.key.height =unit(3, "line"),legend.key.size = unit(3,"line"))
# aa
# dev.off()


# #1c

# pdf(file="survival_SubGroups_OS_sample_1c.pdf", width = 1.8, height =1.3)
# ## use these colors if need colored plot, the colors are in sequnce 1a,1b,1c,1d,2a,2b,2c,2d,2e,3a,3b,3c
# #ggsurvplot(fit, data = s, risk.table = TRUE, pval = TRUE,palette = c( "#1B9E77", "aquamarine2",  "limegreen","#666666","#F46D43",  "darkred",  "#FA9FB5","#E6AB02","magenta1", "purple3",  "deepskyblue",  "gray4"))
# aa=ggsurvplot(
#   fit, data = surv_psn,
#   risk.table = FALSE, pval = FALSE,palette = c("grey", "grey",  "black","grey","grey",  "grey",  "grey","grey","grey", "grey","grey", "grey"),ylab="Overall Survival"
#   ,legend.title="SubGroup",pval.size = 18,size = 0.25,legend="none",censor.size=0.2,

#   ggtheme = theme_survminer(

#     font.x = c(5, "bold", "black"),
#     font.y = c(5, "bold", "black"),
#     font.legend = c(5, "bold"), 
#     font.tickslab = c(5, "bold", "black"),
#   )
# )
# #aa$plot <- aa$plot + theme(legend.key.height =unit(3, "line"),legend.key.size = unit(3,"line"))
# aa
# dev.off()


# #1d

# pdf(file="survival_SubGroups_OS_sample_1d.pdf", width = 1.8, height =1.3)
# ## use these colors if need colored plot, the colors are in sequnce 1a,1b,1c,1d,2a,2b,2c,2d,2e,3a,3b,3c
# #ggsurvplot(fit, data = s, risk.table = TRUE, pval = TRUE,palette = c( "#1B9E77", "aquamarine2",  "limegreen","#666666","#F46D43",  "darkred",  "#FA9FB5","#E6AB02","magenta1", "purple3",  "deepskyblue",  "gray4"))
# aa=ggsurvplot(
#   fit, data = surv_psn,
#   risk.table = FALSE, pval = FALSE,palette = c("grey", "grey",  "grey","black","grey",  "grey",  "grey","grey","grey", "grey","grey", "grey"),ylab="Overall Survival"
#   ,legend.title="SubGroup",pval.size = 18,size = 0.25,legend="none",censor.size=0.2,

#   ggtheme = theme_survminer(

#     font.x = c(5, "bold", "black"),
#     font.y = c(5, "bold", "black"),
#     font.legend = c(5, "bold"), 
#     font.tickslab = c(5, "bold", "black"),
#   )
# )
# #aa$plot <- aa$plot + theme(legend.key.height =unit(3, "line"),legend.key.size = unit(3,"line"))
# aa
# dev.off()


# #2a

# pdf(file="survival_SubGroups_OS_sample_2a.pdf", width = 1.8, height =1.3)
# ## use these colors if need colored plot, the colors are in sequnce 1a,1b,1c,1d,2a,2b,2c,2d,2e,3a,3b,3c
# #ggsurvplot(fit, data = s, risk.table = TRUE, pval = TRUE,palette = c( "#1B9E77", "aquamarine2",  "limegreen","#666666","#F46D43",  "darkred",  "#FA9FB5","#E6AB02","magenta1", "purple3",  "deepskyblue",  "gray4"))
# aa=ggsurvplot(
#   fit, data = surv_psn,
#   risk.table = FALSE, pval = FALSE,palette = c("grey", "grey",  "grey","grey","black",  "grey",  "grey","grey","grey", "grey","grey", "grey"),ylab="Overall Survival"
#   ,legend.title="SubGroup",pval.size = 18,size = 0.25,legend="none",censor.size=0.2,

#   ggtheme = theme_survminer(

#     font.x = c(5, "bold", "black"),
#     font.y = c(5, "bold", "black"),
#     font.legend = c(5, "bold"), 
#     font.tickslab = c(5, "bold", "black"),
#   )
# )
# #aa$plot <- aa$plot + theme(legend.key.height =unit(3, "line"),legend.key.size = unit(3,"line"))
# aa
# dev.off()


# #2b

# pdf(file="survival_SubGroups_OS_sample_2b.pdf", width = 1.8, height =1.3)
# ## use these colors if need colored plot, the colors are in sequnce 1a,1b,1c,1d,2a,2b,2c,2d,2e,3a,3b,3c
# #ggsurvplot(fit, data = s, risk.table = TRUE, pval = TRUE,palette = c( "#1B9E77", "aquamarine2",  "limegreen","#666666","#F46D43",  "darkred",  "#FA9FB5","#E6AB02","magenta1", "purple3",  "deepskyblue",  "gray4"))
# aa=ggsurvplot(
#   fit, data = surv_psn,
#   risk.table = FALSE, pval = FALSE,palette = c("grey", "grey",  "grey","grey","grey",  "black",  "grey","grey","grey", "grey","grey", "grey"),ylab="Overall Survival"
#   ,legend.title="SubGroup",pval.size = 18,size = 0.25,legend="none",censor.size=0.2,

#   ggtheme = theme_survminer(

#     font.x = c(5, "bold", "black"),
#     font.y = c(5, "bold", "black"),
#     font.legend = c(5, "bold"), 
#     font.tickslab = c(5, "bold", "black"),
#   )
# )
# #aa$plot <- aa$plot + theme(legend.key.height =unit(3, "line"),legend.key.size = unit(3,"line"))
# aa
# dev.off()


# #2c

# pdf(file="survival_SubGroups_OS_sample_2c.pdf", width = 1.8, height =1.3)
# ## use these colors if need colored plot, the colors are in sequnce 1a,1b,1c,1d,2a,2b,2c,2d,2e,3a,3b,3c
# #ggsurvplot(fit, data = s, risk.table = TRUE, pval = TRUE,palette = c( "#1B9E77", "aquamarine2",  "limegreen","#666666","#F46D43",  "darkred",  "#FA9FB5","#E6AB02","magenta1", "purple3",  "deepskyblue",  "gray4"))
# aa=ggsurvplot(
#   fit, data = surv_psn,
#   risk.table = FALSE, pval = FALSE,palette = c("grey", "grey",  "grey","grey","grey",  "grey",  "black","grey","grey", "grey","grey", "grey"),ylab="Overall Survival"
#   ,legend.title="SubGroup",pval.size = 18,size = 0.25,legend="none",censor.size=0.2,

#   ggtheme = theme_survminer(

#     font.x = c(5, "bold", "black"),
#     font.y = c(5, "bold", "black"),
#     font.legend = c(5, "bold"), 
#     font.tickslab = c(5, "bold", "black"),
#   )
# )
# #aa$plot <- aa$plot + theme(legend.key.height =unit(3, "line"),legend.key.size = unit(3,"line"))
# aa
# dev.off()


# #2d

# pdf(file="survival_SubGroups_OS_sample_2d.pdf", width = 1.8, height =1.3)
# ## use these colors if need colored plot, the colors are in sequnce 1a,1b,1c,1d,2a,2b,2c,2d,2e,3a,3b,3c
# #ggsurvplot(fit, data = s, risk.table = TRUE, pval = TRUE,palette = c( "#1B9E77", "aquamarine2",  "limegreen","#666666","#F46D43",  "darkred",  "#FA9FB5","#E6AB02","magenta1", "purple3",  "deepskyblue",  "gray4"))
# aa=ggsurvplot(
#   fit, data = surv_psn,
#   risk.table = FALSE, pval = FALSE,palette = c("grey", "grey",  "grey","grey","grey",  "grey",  "grey","black","grey", "grey","grey", "grey"),ylab="Overall Survival"
#   ,legend.title="SubGroup",pval.size = 18,size = 0.25,legend="none",censor.size=0.2,

#   ggtheme = theme_survminer(

#     font.x = c(5, "bold", "black"),
#     font.y = c(5, "bold", "black"),
#     font.legend = c(5, "bold"), 
#     font.tickslab = c(5, "bold", "black"),
#   )
# )
# #aa$plot <- aa$plot + theme(legend.key.height =unit(3, "line"),legend.key.size = unit(3,"line"))
# aa
# dev.off()



# #2e

# pdf(file="survival_SubGroups_OS_sample_2e.pdf", width = 1.8, height =1.3)
# ## use these colors if need colored plot, the colors are in sequnce 1a,1b,1c,1d,2a,2b,2c,2d,2e,3a,3b,3c
# #ggsurvplot(fit, data = s, risk.table = TRUE, pval = TRUE,palette = c( "#1B9E77", "aquamarine2",  "limegreen","#666666","#F46D43",  "darkred",  "#FA9FB5","#E6AB02","magenta1", "purple3",  "deepskyblue",  "gray4"))
# aa=ggsurvplot(
#   fit, data = surv_psn,
#   risk.table = FALSE, pval = FALSE,palette = c("grey", "grey",  "grey","grey","grey",  "grey",  "grey","grey","black", "grey","grey", "grey"),ylab="Overall Survival"
#   ,legend.title="SubGroup",pval.size = 18,size = 0.25,legend="none",censor.size=0.2,

#   ggtheme = theme_survminer(

#     font.x = c(5, "bold", "black"),
#     font.y = c(5, "bold", "black"),
#     font.legend = c(5, "bold"), 
#     font.tickslab = c(5, "bold", "black"),
#   )
# )
# #aa$plot <- aa$plot + theme(legend.key.height =unit(3, "line"),legend.key.size = unit(3,"line"))
# aa
# dev.off()

# #3a

# pdf(file="survival_SubGroups_OS_sample_3a.pdf", width = 1.8, height =1.3)
# ## use these colors if need colored plot, the colors are in sequnce 1a,1b,1c,1d,2a,2b,2c,2d,2e,3a,3b,3c
# #ggsurvplot(fit, data = s, risk.table = TRUE, pval = TRUE,palette = c( "#1B9E77", "aquamarine2",  "limegreen","#666666","#F46D43",  "darkred",  "#FA9FB5","#E6AB02","magenta1", "purple3",  "deepskyblue",  "gray4"))
# aa=ggsurvplot(
#   fit, data = surv_psn,
#   risk.table = FALSE, pval = FALSE,palette = c("grey", "grey",  "grey","grey","grey",  "grey",  "grey","grey","grey", "black","grey", "grey"),ylab="Overall Survival"
#   ,legend.title="SubGroup",pval.size = 18,size = 0.25,legend="none",censor.size=0.2,

#   ggtheme = theme_survminer(

#     font.x = c(5, "bold", "black"),
#     font.y = c(5, "bold", "black"),
#     font.legend = c(5, "bold"), 
#     font.tickslab = c(5, "bold", "black"),
#   )
# )
# #aa$plot <- aa$plot + theme(legend.key.height =unit(3, "line"),legend.key.size = unit(3,"line"))
# aa
# dev.off()


# #3b

# pdf(file="survival_SubGroups_OS_sample_3b.pdf", width = 1.8, height =1.3)
# ## use these colors if need colored plot, the colors are in sequnce 1a,1b,1c,1d,2a,2b,2c,2d,2e,3a,3b,3c
# #ggsurvplot(fit, data = s, risk.table = TRUE, pval = TRUE,palette = c( "#1B9E77", "aquamarine2",  "limegreen","#666666","#F46D43",  "darkred",  "#FA9FB5","#E6AB02","magenta1", "purple3",  "deepskyblue",  "gray4"))
# aa=ggsurvplot(
#   fit, data = surv_psn,
#   risk.table = FALSE, pval = FALSE,palette = c("grey", "grey",  "grey","grey","grey",  "grey",  "grey","grey","grey", "grey","black", "grey"),ylab="Overall Survival"
#   ,legend.title="SubGroup",pval.size = 18,size = 0.25,legend="none",censor.size=0.2,

#   ggtheme = theme_survminer(

#     font.x = c(5, "bold", "black"),
#     font.y = c(5, "bold", "black"),
#     font.legend = c(5, "bold"), 
#     font.tickslab = c(5, "bold", "black"),
#   )
# )
# #aa$plot <- aa$plot + theme(legend.key.height =unit(3, "line"),legend.key.size = unit(3,"line"))
# aa
# dev.off()


# #3c

# pdf(file="survival_SubGroups_OS_sample_3c.pdf", width = 1.8, height =1.3)
# ## use these colors if need colored plot, the colors are in sequnce 1a,1b,1c,1d,2a,2b,2c,2d,2e,3a,3b,3c
# #ggsurvplot(fit, data = s, risk.table = TRUE, pval = TRUE,palette = c( "#1B9E77", "aquamarine2",  "limegreen","#666666","#F46D43",  "darkred",  "#FA9FB5","#E6AB02","magenta1", "purple3",  "deepskyblue",  "gray4"))
# aa=ggsurvplot(
#   fit, data = surv_psn,
#   risk.table = FALSE, pval = FALSE,palette = c("grey", "grey",  "grey","grey","grey",  "grey",  "grey","grey","grey", "grey","grey", "black"),ylab="Overall Survival"
#   ,legend.title="SubGroup",pval.size = 18,size = 0.25,legend="none",censor.size=0.2,

#   ggtheme = theme_survminer(

#     font.x = c(5, "bold", "black"),
#     font.y = c(5, "bold", "black"),
#     font.legend = c(5, "bold"), 
#     font.tickslab = c(5, "bold", "black"),
#   )
# )
# #aa$plot <- aa$plot + theme(legend.key.height =unit(3, "line"),legend.key.size = unit(3,"line"))
# aa
# dev.off()