setwd("D:/ChIP/")
library(ggplot2)
library(ggrepel)
library(GenomicRanges)
library(tidyverse)
library(xlsx)
library(patchwork)
library(scales)
library(tidyr)
library(gridExtra)
library(cowplot)
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(ggpubr)
library(ComplexHeatmap)
############################
######### figure 1 #########

##### c distance ######

# USE BASH: Find the nearest TSS for each interval
#

# awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2, $3, $1"_"$2"_"$3}' /home/castillaa/chip/revisions/input/sp6_idrOptPeaks_pool.bed > /home/castillaa/chip/revisions/input/sp6_pool_4c.bed
# awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2, $3, $1"_"$2"_"$3}' /home/castillaa/chip/revisions/input/sp8_idrOptPeaks_pool.bed > /home/castillaa/chip/revisions/input/sp8_pool_4c.bed
# sort -k1,1 -k2,2n /home/castillaa/chip/revisions/input/sp6_pool_4c.bed > /home/castillaa/chip/revisions/input/sp6_pool_4cs.bed
# sort -k1,1 -k2,2n /home/castillaa/chip/revisions/input/sp8_pool_4c.bed > /home/castillaa/chip/revisions/input/sp8_pool_4cs.bed
# ##1 gene
# bedtools closest -D ref -t first -a /home/castillaa/chip/revisions/input/sp6_pool_4cs.bed -b /home/castillaa/chip/revisions/input/input.tss.sorted.bed > /home/castillaa/chip/revisions/output/distance_to_tss_sp6_pool.bed #
# bedtools closest -D ref -t first -a /home/castillaa/chip/revisions/input/sp8_pool_4cs.bed -b /home/castillaa/chip/revisions/input/input.tss.sorted.bed > /home/castillaa/chip/revisions/output/distance_to_tss_sp8_pool.bed 


# --- Read input BED-like file ---
sp6_tss <- read.table("revisions/input/distance_to_tss_sp6_pool.bed", header=FALSE)
sp8_tss <- read.table("revisions/input/distance_to_tss_sp8_pool.bed", header=FALSE)

# scale distances from V9 (bp) into kb
sp6_tss$V10 <- sp6_tss$V9 / 1000
sp8_tss$V10 <- sp8_tss$V9 / 1000

# --- Proximal vs Distal (±2 kb only) ---
sp6_tss$distance_label <- cut(sp6_tss$V10, 
                              breaks=c(-Inf, -2, 2, Inf), 
                              labels=c("Distal","Proximal","Distal"))
sp8_tss$distance_label <- cut(sp8_tss$V10, 
                              breaks=c(-Inf, -2, 2, Inf), 
                              labels=c("Distal","Proximal","Distal"))

# --- Discrete bins with ±2 kb proximal ---
breaks_discrete <- c(-Inf,-500, -50, -2, 0, 2, 50, 500, Inf)
labels_discrete <- c("< -500","-500 to -50", "-50 to -2",
                     "-2 to 0", "0 to 2", "2 to 50", "50 to 500",
                     "> 500")

sp6_tss$Distance_Group <- cut(sp6_tss$V10, breaks=breaks_discrete, labels=labels_discrete)
sp8_tss$Distance_Group <- cut(sp8_tss$V10, breaks=breaks_discrete, labels=labels_discrete)

# --- Name columns
sp6_tss_sd1 <- sp6_tss
colnames(sp6_tss_sd1) <- c("chr","start","end","peak_id","gene_chr",
                       "gene_start","gene_end","gene","distance",
                       "distance_kb","distance_label", "Distance_Group")
sp8_tss_sd1 <- sp8_tss
colnames(sp8_tss_sd1) <- c("chr","start","end","peak_id","gene_chr",
                           "gene_start","gene_end","gene","distance",
                           "distance_kb","distance_label", "Distance_Group")
# --- Make the Source data 1 with IDR values
sp6_sd1 <- read.table("revisions/idr_redo/results/sp6_idrOptPeaks_pool.bed", header=FALSE)
sp8_sd1 <- read.table("revisions/idr_redo/results/sp8_idrOptPeaks_pool.bed", header=FALSE)
colnames(sp6_sd1) <- c("chr","start","end","def_name","IDR_score_sc",
                           "strand","signal","p-val","q-val",
                           "summit","local_IDR", "global_IDR")
colnames(sp8_sd1) <- c("chr","start","end","def_name","IDR_score_sc",
                       "strand","signal","p-val","q-val",
                       "summit","local_IDR", "global_IDR")

# --- Merge info
sp6_sd1 <- sp6_sd1 %>%
  mutate(peak_id= paste0(chr,"_",start,"_",end) )
sp8_sd1 <- sp8_sd1 %>%
  mutate(peak_id= paste0(chr,"_",start,"_",end) )

sp6_tss_sd1 <- merge(sp6_tss_sd1, sp6_sd1, by="peak_id")
sp8_tss_sd1 <- merge(sp8_tss_sd1, sp8_sd1, by="peak_id")

# --- Save annotated tables ---
xlsx::write.xlsx(sp8_tss_sd1[,c(1:12,17,19:24)], "revisions/output/Source_data1_sp8_sp6_peaks_2kb.xlsx", sheetName="Sp8_peaks")
xlsx::write.xlsx(sp8_tss_sd1[,c(1:12,17,19:24)], "revisions/output/Source_data1_sp8_sp6_peaks_2kb.xlsx", sheetName="Sp6_peaks", append=TRUE)

# --- Colors ---
fill_sp8 <- c("mediumorchid1","mediumorchid1","mediumorchid1",
              "mediumorchid4","mediumorchid4","mediumorchid1","mediumorchid1",
              "mediumorchid1")
fill_sp6 <- c("olivedrab1","olivedrab1","olivedrab1",
              "olivedrab4","olivedrab4","olivedrab1","olivedrab1",
              "olivedrab1")


###SP6
signed_log_distance <- sign(sp6_tss$V10) * log10(abs(sp6_tss$V10) + 1)
ggplot(data.frame(distance = signed_log_distance), aes(x = distance)) +
  geom_density(fill = "olivedrab1", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +     # TSS
  geom_vline(xintercept = sign(c(-2, 2)) * log10(abs(c(-2, 2)) + 1), linetype = "dashed") +
  labs(
    x = "Log10 Distance to TSS (kb)",
    y = "Density"
  ) +
  annotate("rect", xmin = -log10(2+1), xmax = log10(2+1), ymin = 0, ymax = Inf, 
           fill = "olivedrab4", alpha = 0.3) +
  scale_x_continuous(
    breaks = sign(c(-1000, -100, -10, -2, 0, 2, 10, 100, 1000)) *
      log10(abs(c(-1000, -100, -10, -2, 0, 2, 10, 100, 1000)) + 1),
    labels = c("-1000 kb", "-100 kb", "-10 kb", "-2 kb", 
               "TSS", "2 kb", "10 kb", "100 kb", "1000 kb")
  ) +geom_rug(sides = "b", color = "olivedrab4", alpha = 0.3, linewidth = 0.2) +
  theme_minimal(base_size = 16)+
  theme( text = element_text(color = "black"),
         axis.text = element_text(color = "black"),
         axis.title = element_text(color = "black"),
         plot.title = element_text(color = "black"),
         plot.subtitle = element_text(color = "black"),
         legend.text = element_text(color = "black"),
         legend.title = element_text(color = "black"),
         strip.text = element_text(color = "black")
         ,axis.text.x=element_text(angle = -45, vjust = -0.5),
         legend.position = "top")
ggsave("D:/ChIP/revisions/Plots/density_plots_Sp6.png",dpi = 300, height = 4, width=5)

###SP8
signed_log_distance8 <- sign(sp8_tss$V10) * log10(abs(sp8_tss$V10) + 1)
ggplot(data.frame(distance = signed_log_distance8), aes(x = distance)) +
    geom_density(fill = "mediumorchid1", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "black") +     # TSS
    geom_vline(xintercept = sign(c(-2, 2)) * log10(abs(c(-2, 2)) + 1), linetype = "dashed") +
  labs(
    x = "Log10 Distance to TSS (kb)",
    y = "Density"
  ) +
    annotate("rect", xmin = -log10(2+1), xmax = log10(2+1), ymin = 0, ymax = Inf, 
             fill = "mediumorchid4", alpha = 0.3) +
    scale_x_continuous(
      breaks = sign(c(-1000, -100, -10, -2, 0, 2, 10, 100, 1000)) *
        log10(abs(c(-1000, -100, -10, -2, 0, 2, 10, 100, 1000)) + 1),
      labels = c("-1000 kb", "-100 kb", "-10 kb", "-2 kb", 
                 "TSS", "2 kb", "10 kb", "100 kb", "1000 kb")
    ) +geom_rug(sides = "b", color = "mediumorchid4", alpha = 0.3, linewidth = 0.2) +
  theme_minimal(base_size = 16)+
  theme( text = element_text(color = "black"),
         axis.text = element_text(color = "black"),
         axis.title = element_text(color = "black"),
         plot.title = element_text(color = "black"),
         plot.subtitle = element_text(color = "black"),
         legend.text = element_text(color = "black"),
         legend.title = element_text(color = "black"),
         strip.text = element_text(color = "black")
         ,axis.text.x=element_text(angle = -45, vjust = -0.5),
         legend.position = "top")
ggsave("D:/ChIP/revisions/Plots/density_plots_Sp8.png",dpi = 300, height = 4, width=5)

# ##BOTH
# # --- Transform both datasets ---
# signed_log_distance6 <- data.frame(distance = sign(sp6_tss$V10) * log10(abs(sp6_tss$V10) + 1),
#                    sample = "Sp6")
# signed_log_distance8 <- data.frame(distance = sign(sp8_tss$V10) * log10(abs(sp8_tss$V10) + 1),
#                     sample = "Sp8")
# # Combine them
# df_combined_signed_log <- rbind(signed_log_distance6, signed_log_distance8)
# # --- Plot both in one ggplot ---
# distribution_sp6_sp8 <- ggplot(df_combined_signed_log, aes(x = distance, fill = sample, color = sample)) +
#     
#     # Density curves
#     geom_density(alpha = 0.4, linewidth = 1) +
#     geom_vline(xintercept = 0, linetype = "dotted", color = "black") +      # Add vertical reference lines 
#     geom_vline(xintercept = sign(c(-2, 2)) * log10(abs(c(-2, 2)) + 1),# TSS
#                linetype = "dashed", color = "black") +
#     annotate("rect", xmin = -log10(2 + 1), xmax = log10(2 + 1), # Promoter-proximal shading
#              ymin = 0, ymax = Inf,
#              fill = "gray70", alpha = 0.2) +
#     geom_rug(   # Add rug plots (tiny dots for individual data points)
#       data = subset(df_combined_signed_log, sample == "Sp6"),
#       aes(color = sample), sides = "b", alpha = 0.3, linewidth = 0.3
#     ) +
#     geom_rug(
#       data = subset(df_combined_signed_log, sample == "Sp8"),
#       aes(color = sample), sides = "t", alpha = 0.3, linewidth = 0.3
#     ) +
#     # Axes and titles
#     scale_x_continuous(
#       breaks = sign(c(-1000, -100, -10, -2, 0, 2, 10, 100, 1000)) *
#         log10(abs(c(-1000, -100, -10, -2, 0, 2, 10, 100, 1000)) + 1),
#       labels = c("-1000 kb", "-100 kb", "-10 kb", "-2 kb", 
#                  "TSS", "2 kb", "10 kb", "100 kb", "1000 kb")
#     ) +
#     labs(
#       x = "Distance to TSS (kb)",
#       y = "Density",
#       fill = "Sample",
#       color = "Sample"
#     ) +
#     scale_fill_manual(values = c("Sp6" = "olivedrab1", "Sp8" = "mediumorchid1")) +
#     scale_color_manual(values = c("Sp6" = "olivedrab4", "Sp8" = "mediumorchid4")) +
#   theme_minimal(base_size = 16)+
#   theme( text = element_text(color = "black"),
#          axis.text = element_text(color = "black"),
#          axis.title = element_text(color = "black"),
#          plot.title = element_text(color = "black"),
#          plot.subtitle = element_text(color = "black"),
#          legend.text = element_text(color = "black"),
#          legend.title = element_text(color = "black"),
#          strip.text = element_text(color = "black")
#          ,axis.text.x=element_text(angle = -45, vjust = -0.5),
#          legend.position = "top")
# distribution_sp6_sp8_noleg <- distribution_sp6_sp8+theme(legend.position = "none")
# ggsave("D:/ChIP/revisions/Plots/density_plots_distance.png",distribution_sp6_sp8,dpi = 300)
# ggsave("D:/ChIP/revisions/Plots/density_plots_distance_noleg.png",distribution_sp6_sp8_noleg,dpi = 300)

###################################### Split into proximal/distal using ±2 kb
# USE BASH: Split Sp6 into proximal/distal using ±2 kb
# awk '{if ($9 >= -2000 && $9 <= 2000) 
#          print > "/home/castillaa/chip/revisions/output/sp6_tss_pool_proximal_2kb.bed"; 
#      else 
#          print > "/home/castillaa/chip/revisions/output/sp6_tss_pool_distal_2kb.bed"}' \
# /home/castillaa/chip/revisions/output/distance_to_tss_sp6_pool.bed
# 
# # Split Sp8 into proximal/distal using ±2 kb
# awk '{if ($9 >= -2000 && $9 <= 2000) 
#          print > "/home/castillaa/chip/revisions/output/sp8_tss_pool_proximal_2kb.bed"; 
#      else 
#          print > "/home/castillaa/chip/revisions/output/sp8_tss_pool_distal_2kb.bed"}' \
# /home/castillaa/chip/revisions/output/distance_to_tss_sp8_pool.bed


##### e motives and go terms obtained from GREAT  ##### 
#rownumber from excel file -4 (as skip=4) # 
sp6_distal <-read.table("D:/ChIP/revisions/great_bpe/greatExportAll_Sp6all.tsv", skip = 4, sep = "\t") 
sp8_distal <-   read.table("D:/ChIP/revisions/great_bpe/2kb/greatExportAll_sp8_dist_2kb.tsv", skip = 4, sep = "\t") 
sp8_proximal <- read.table("D:/ChIP/revisions/great_bpe/2kb/greatExportAll_sp8_prox_2kb.tsv", skip = 4, sep = "\t") 

sp6_distal <- sp6_distal %>% filter(V1=="GO Biological Process" & V16 < 0.05) %>% mutate(FDR=V16) %>% mutate(Count=V19) 
sp8_distal <-   sp8_distal %>% filter(V1=="GO Biological Process" & V16 < 0.05) %>% mutate(FDR=V16) %>% mutate(Count=V19)
sp8_proximal <- sp8_proximal %>% filter(V1=="GO Biological Process" & V16 < 0.05) %>% mutate(FDR=V16) %>% mutate(Count=V19)

write.xlsx(sp8_proximal, "D:/ChIP/revisions/output/Source_data_2_GREAT_BP_terms.xlsx", sheetName="Sp8_proximal")
write.xlsx(sp8_distal, "D:/ChIP/revisions/output/Source_data_2_GREAT_BP_terms.xlsx", sheetName="Sp8_distal", append=TRUE)
write.xlsx(sp6_distal, "D:/ChIP/revisions/output/Source_data_2_GREAT_BP_terms.xlsx", sheetName="Sp6_distal", append=TRUE) 


####################### GO terms fig 1 ##############
#dotplot #########
# TOP 30 SELECTED

p5 <- ggplot(data = sp6_distal[c(1,2,3,7,9), ], aes(y=reorder(V3, -log10(FDR)), x=V19/V20))+
  geom_point(stat = "identity", color="black", pch= 21, aes(size=Count, fill=-log10(FDR)))+ 
  scale_y_discrete(labels = wrap_format(34))+ scale_x_continuous(position="top")+ 
  xlim(c(0,0.75))+
  theme_bw(base_size = 10)+
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size =14, colour="black"),
        axis.text.x = element_text(colour="black"), panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())+ labs(x="Gene Ratio")+ 
  scale_fill_gradient(low="olivedrab1",high = "olivedrab4")+ 
  scale_size_continuous(range=c(3,6)) 


p2 <- ggplot(data = sp8_proximal[c(2,4,6,9,11), ], aes(y=reorder(V3, -log10(FDR)), x=V19/V20))+ 
  geom_point(stat = "identity", color="black", pch= 21, aes(size=Count, fill=-log10(FDR)))+ 
  scale_y_discrete(labels = wrap_format(34))+ scale_x_continuous(position="top")+ xlim(c(0,0.75))+ 
  theme_bw(base_size = 10)+ 
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 14, colour="black"), 
                    axis.text.x = element_text(colour="black"), panel.grid.major.x = element_blank(),
                    panel.grid.minor.x = element_blank())+ labs(x="Gene Ratio")+ 
  scale_fill_gradient(low="mediumorchid1",high = "mediumorchid4")+ scale_size_continuous(range=c(3,6)) 

p1 <-ggplot(data = sp8_distal[c(1,2, 4, 6,8), ], aes(y=reorder(V3, -log10(FDR)), x=V19/V20))+ 
  geom_point(stat = "identity", color="black", pch= 21, aes(size=Count, fill=-log10(FDR)))+ 
  scale_y_discrete(labels = wrap_format(34))+ scale_x_continuous(position="top")+ xlim(c(0,0.75))+
  scale_size_continuous(range=c(3,6))+ theme_bw(base_size = 10)+ 
  theme(axis.title.y = element_blank(),
                                                       axis.text.y = element_text(size = 14, colour="black"),
                                                         axis.text.x = element_text(colour="black"),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())+ 
  labs(x="Gene Ratio")+ scale_fill_gradient(low="mediumorchid1",high = "mediumorchid4")


combined_plot <- plot_grid( p2, p1, p5,p5, nrow = 2, 
                            # 2x2 grid 
                            align = "hv", # Align horizontally and vertically 
                            axis = "tblr" # Align axes across all plots 
                            ) 
ggsave("D:/ChIP/revisions/Plots/Great_bpe_Fig1_aligned_plots.pdf", combined_plot, width = 11, height = 6.5)


##### d heatmaps proximal distal ######
#data from cistrome galaxy server- input sp6 and sp8 optimal pool peaks and run default in mm10 genome

x<-c(-1500,-1470,-1440,-1410,-1380,-1350,-1320,-1290,-1260,-1230,-1200,-1170,-1140,-1110,-1080,-1050,-1020,-990,-960,-930,-900,-870,-840,-810,-780,-750,-720,-690,-660,-630,-600,-570,-540,-510,-480,-450,-420,-390,-360,-330,-300,-270,-240,-210,-180,-150,-120,-90,-60,-30,0,30,60,90,120,150,180,210,240,270,300,330,360,390,420,450,480,510,540,570,600,630,660,690,720,750,780,810,840,870,900,930,960,990,1020,1050,1080,1110,1140,1170,1200,1230,1260,1290,1320,1350,1380,1410,1440)
y0<-c(0.160766497501,0.164041549061,0.163343856336,0.163371322351,0.162497199194,0.162579346467,0.164006151255,0.162499669865,0.170795815569,0.170636457466,0.163493506242,0.160343827381,0.161736034425,0.163089603685,0.160396683951,0.166250470426,0.169639128553,0.17400551237,0.176809805877,0.176584080441,0.178723803869,0.170582335936,0.167186799062,0.169017576301,0.179256424287,0.169469438098,0.175875934047,0.175711908291,0.174466693999,0.171224526737,0.172701216607,0.170026834885,0.16934493445,0.171951535911,0.182365404585,0.181192496072,0.191752454942,0.194282680604,0.202586501938,0.211308129231,0.217322043271,0.215701444605,0.224938424972,0.24203029653,0.250250666168,0.270824722222,0.28970941915,0.306375244765,0.326612108343,0.338562020565,0.351338649861,0.329414476881,0.327479004079,0.300546369653,0.280130705416,0.257751889601,0.235941953473,0.225132151644,0.215114857585,0.202779377896,0.201261642577,0.192130636427,0.196277021703,0.186158636671,0.191562514722,0.194650118695,0.188553457327,0.18332324135,0.186272762637,0.179191371207,0.182166350695,0.179928002594,0.179385692512,0.168412936209,0.171411157757,0.171170513763,0.172547686762,0.16790188897,0.168191808407,0.165233518182,0.166351378158,0.161390213716,0.167808028408,0.170823746466,0.173270141448,0.171503295105,0.17034358605,0.161526965416,0.172102032769,0.173033490808,0.164456426933,0.175750221136,0.174654369849,0.180051124138,0.170625255874,0.170927755862,0.17217627514,0.166777588999,0.157105222928)
y1<-c(0.154334638313,0.145202095821,0.154213891277,0.15627609759,0.150656161319,0.164839790646,0.15043392844,0.150459310363,0.160017212552,0.163978082917,0.166555920112,0.171819385477,0.168006111896,0.166234839542,0.162806857905,0.159636749092,0.168901515854,0.15566109658,0.1605124659,0.160440086967,0.162646485536,0.161849261543,0.15596161948,0.149484310241,0.145623999632,0.143798069619,0.135735537291,0.140149067159,0.153056767681,0.151691429833,0.159435103239,0.166059861901,0.158223527814,0.157028375508,0.158844542916,0.167537663889,0.175669793858,0.183325653156,0.185013866867,0.19105107637,0.191645819566,0.201135324395,0.21258539143,0.21959896378,0.230303535846,0.236324820216,0.240955536785,0.272948897553,0.284044637727,0.294039525977,0.308232254007,0.302085212873,0.292331468282,0.285458834746,0.253626849096,0.237887691177,0.228880239135,0.222186847235,0.199484092767,0.202623020837,0.187141259105,0.184411941949,0.179967004512,0.179823594217,0.170636400643,0.162103739411,0.154152796307,0.154115293936,0.167342677474,0.172410153661,0.166969617288,0.163235122766,0.160123612084,0.163236658568,0.159657447757,0.159994392348,0.168655253598,0.165981737215,0.164044235503,0.160589512583,0.148077876431,0.150818812262,0.150797922824,0.152402994937,0.156035775776,0.152918090148,0.159296507527,0.148131485367,0.163019544756,0.169211740828,0.171343946202,0.172670572738,0.165093206547,0.169205985668,0.169189210454,0.175720506824,0.170043320943,0.165402663833,0.163983080966)
ymax <- max(y0,y1)
ymin <- min(y0,y1)
yquart <- (ymax-ymin)/4
pdf("revisions/Plots/cistrome_cons_opt_sp8.pdf",height=5,width=5)
plot(x,y0,type="l",col="darkorchid1", lwd=3,
     main="Phastcons around the Center of Sites",
     xlab="Distance from the Center (bp)",ylab="Average Phastcons",ylim=c(ymin-yquart,ymax+yquart),
     cex.main = 1.5, cex.lab = 1.4, cex.axis = 1.3)
abline(v=0, lwd=1)
legend("topright",c('Sp8'),col=c("darkorchid1"),lty=c(1,1), lwd=c(3,3),cex = 1.2)
dev.off()

pdf("revisions/Plots/cistrome_cons_opt_sp6.pdf",height=5,width=5)
plot(x,y1,type="l",col="chartreuse2", lwd=3,
     main="Phastcons around the Center of Sites",
     xlab="Distance from the Center (bp)",ylab="Average Phastcons",ylim=c(ymin-yquart,ymax+yquart),
     cex.main = 1.5, cex.lab = 1.4, cex.axis = 1.3)
abline(v=0, lwd=1)
legend("topright",c('Sp6'),col=c("chartreuse2"),lty=c(1), lwd=c(3),cex = 1.2)
dev.off()


##################### HEATMAPS FIG S2 ######################
# #######USE BASH: heatmaps-------------------------------------
# 
# ##atac--------------------------------
# computeMatrix reference-point -S /home/castillaa/chip/desanlis_2020_atac_ect/GSM4323465_WT.8_ectoderm_treat_pileup.bw \
# -R /home/castillaa/chip/revisions/input/sp6_pool_4cs.bed \
# -a 3000 -b 3000 -p 16 \
# -o /home/castillaa/chip/bigwigs_bamcoverage/matrix/sp6_pool_peaks_in_atac.gz --referencePoint center 
# 
# plotHeatmap -m /home/castillaa/chip/bigwigs_bamcoverage/matrix/sp6_pool_peaks_in_atac.gz \
# -o /home/castillaa/chip/heatmaps/sp6_pool_peaks_in_atac.pdf --startLabel center \
# --colorList  blue,white,red --legendLocation none --zMax  1 
# 
# computeMatrix reference-point -S /home/castillaa/chip/desanlis_2020_atac_ect/GSM4323465_WT.8_ectoderm_treat_pileup.bw \
# -R /home/castillaa/chip/revisions/input/sp8_pool_4cs.bed \
# -a 3000 -b 3000 -p 16 \
# -o /home/castillaa/chip/bigwigs_bamcoverage/matrix/sp8_pool_peaks_in_atac.gz --referencePoint center 
# 
# plotHeatmap -m /home/castillaa/chip/bigwigs_bamcoverage/matrix/sp8_pool_peaks_in_atac.gz \
# -o /home/castillaa/chip/heatmaps/sp8_pool_peaks_in_atac.pdf --startLabel center \
# --colorList  blue,white,red --legendLocation none --zMax  1 
# 
# 
# ##me2--------------------------------
# computeMatrix reference-point -S /home/castillaa/chip/bigwigs_bamcoverage/rochp_me2_rep1_sorted_filt_chr.bw \
# -R /home/castillaa/chip/revisions/input/sp8_pool_4cs.bed \
# -a 3000 -b 3000 -p 16 \
# -o /home/castillaa/chip/bigwigs_bamcoverage/matrix/sp8_pool_peaks_in_me2.gz --referencePoint center 
# 
# plotHeatmap -m /home/castillaa/chip/bigwigs_bamcoverage/matrix/sp8_pool_peaks_in_me2.gz \
# -o /home/castillaa/chip/heatmaps/sp8_pool_peaks_in_me2.pdf --startLabel center \
# --colorList  blue,white,red --legendLocation none --zMax  4 
# 
# computeMatrix reference-point -S /home/castillaa/chip/bigwigs_bamcoverage/rochp_me2_rep1_sorted_filt_chr.bw \
# -R /home/castillaa/chip/revisions/input/sp6_pool_4cs.bed \
# -a 3000 -b 3000 -p 16 \
# -o /home/castillaa/chip/bigwigs_bamcoverage/matrix/sp6_pool_peaks_in_me2.gz --referencePoint center 
# 
# plotHeatmap -m /home/castillaa/chip/bigwigs_bamcoverage/matrix/sp6_pool_peaks_in_me2.gz \
# -o /home/castillaa/chip/heatmaps/sp6_pool_peaks_in_me2.pdf --startLabel center \
# --colorList  blue,white,red --legendLocation none --zMax  4 
# 
# 
# ##signal--------------------------------
# computeMatrix reference-point -S /home/castillaa/chip/bigwigs_bamcoverage/rochp_sp8_rep1_sorted_filt_chr.bw \
# /home/castillaa/chip/bigwigs_bamcoverage/rochp_sp8_rep2_sorted_filt_chr.bw \
# -R //home/castillaa/chip/revisions/input/sp8_pool_4cs.bed \
# -a 3000 -b 3000 -p 16 \
# -o /home/castillaa/chip/bigwigs_bamcoverage/matrix/sp8_pool_peaks_in_sp8.gz --referencePoint center 
# 
# plotHeatmap -m /home/castillaa/chip/bigwigs_bamcoverage/matrix/sp8_pool_peaks_in_sp8.gz \
# -o /home/castillaa/chip/heatmaps/sp8_pool_peaks_in_sp8.pdf --startLabel center \
# --colorMap  bwr --legendLocation none 
# 
# # --zMax 7.5 7.5 4 4 
# 
# 
# computeMatrix reference-point -S /home/castillaa/chip/bigwigs_bamcoverage/chipd_sp6_rep1_sorted_filt_chr.bw \
# /home/castillaa/chip/bigwigs_bamcoverage/chipd_sp6_rep2_sorted_filt_chr.bw \
# -R /home/castillaa/chip/revisions/input/sp6_pool_4cs.bed \
# -a 3000 -b 3000 -p 16 \
# -o /home/castillaa/chip/bigwigs_bamcoverage/matrix/sp6_pool_peaks_in_sp6.gz --referencePoint center 
# 
# plotHeatmap -m /home/castillaa/chip/bigwigs_bamcoverage/matrix/sp6_pool_peaks_in_sp6.gz \
# -o /home/castillaa/chip/heatmaps/sp6_pool_peaks_in_sp6.pdf --startLabel center \
# --colorMap  bwr --legendLocation none 
# 
# # --zMax 7.5 7.5 4 4 
# 
# #rep3
# computeMatrix reference-point -S /home/castillaa/chip/output_sp6_rep3/bwa/merged_library/bigwig/SP6_REP1.mLb.clN.bigWig \
# /home/castillaa/chip/output_sp6_rep3/bwa/merged_library/bigwig/SP6_input_REP1.mLb.clN.bigWig \
# -R /home/castillaa/chip/revisions/input/sp6_pool_4cs.bed \
# -a 3000 -b 3000 -p 16 \
# -o /home/castillaa/chip/bigwigs_bamcoverage/matrix/sp6_pool_peaks_in_sp6_r3.gz --referencePoint center 
# 
# plotHeatmap -m /home/castillaa/chip/bigwigs_bamcoverage/matrix/sp6_pool_peaks_in_sp6_r3.gz \
# -o /home/castillaa/chip/heatmaps/sp6_pool_peaks_in_sp6_r3.pdf --startLabel center \
# --colorMap  bwr --legendLocation none --zMax 0.8
# 
# ##signal in dlx5 --------------------------------
# computeMatrix reference-point -S /home/castillaa/chip/revisions/input/DLX5_HL_11.bw \
# /home/castillaa/chip/revisions/input/DLX5_FL_11.bw \
# /home/castillaa/chip/revisions/input/DLX5_BA_11.bw \
# -R //home/castillaa/chip/revisions/input/sp8_pool_4cs.bed \
# -a 3000 -b 3000 -p 16 \
# -o /home/castillaa/chip/bigwigs_bamcoverage/matrix/sp8_pool_peaks_in_dlx5.gz --referencePoint center 
# 
# plotHeatmap -m /home/castillaa/chip/bigwigs_bamcoverage/matrix/sp8_pool_peaks_in_dlx5.gz \
# -o /home/castillaa/chip/heatmaps/sp8_pool_peaks_in_dlx5.pdf --startLabel center \
# --colorMap  bwr --legendLocation none 
# 
# # --zMax 7.5 7.5 4 4 
# 
# 
# computeMatrix reference-point -S /home/castillaa/chip/revisions/input/DLX5_HL_11.bw \
# /home/castillaa/chip/revisions/input/DLX5_FL_11.bw \
# /home/castillaa/chip/revisions/input/DLX5_BA_11.bw \
# -R //home/castillaa/chip/revisions/input/sp6_pool_4cs.bed \
# -a 3000 -b 3000 -p 16 \
# -o /home/castillaa/chip/bigwigs_bamcoverage/matrix/sp6_pool_peaks_in_dlx5.gz --referencePoint center 
# 
# plotHeatmap -m /home/castillaa/chip/bigwigs_bamcoverage/matrix/sp6_pool_peaks_in_dlx5.gz \
# -o /home/castillaa/chip/heatmaps/sp6_pool_peaks_in_dlx5.pdf --startLabel center \
# --colorMap  bwr --legendLocation none 
# 
# # --zMax 7.5 7.5 4 4 
# 
# ##signal in dlx5 common peaks--------------------------------
# computeMatrix reference-point -S /home/castillaa/chip/revisions/input/DLX5_HL_11.bw \
# /home/castillaa/chip/revisions/input/DLX5_FL_11.bw \
# /home/castillaa/chip/revisions/input/DLX5_BA_11.bw \
# -R //home/castillaa/chip/revisions/input/common_with_gb_for_great_s.bed \
# -a 3000 -b 3000 -p 16 \
# -o /home/castillaa/chip/bigwigs_bamcoverage/matrix/common_pool_peaks_in_dlx5.gz --referencePoint center 
# 
# plotHeatmap -m /home/castillaa/chip/bigwigs_bamcoverage/matrix/common_pool_peaks_in_dlx5.gz \
# -o /home/castillaa/chip/heatmaps/common_pool_peaks_in_dlx5.pdf --startLabel center \
# --colorMap  bwr --legendLocation none 
# 
# # --zMax 7.5 7.5 4 4 
# 
# 
# computeMatrix reference-point -S /home/castillaa/chip/revisions/input/DLX5_HL_11.bw \
# /home/castillaa/chip/revisions/input/DLX5_FL_11.bw \
# /home/castillaa/chip/revisions/input/DLX5_BA_11.bw \
# -R //home/castillaa/chip/revisions/input/sp8_bpe_for_great_s.bed \
# -a 3000 -b 3000 -p 16 \
# -o /home/castillaa/chip/bigwigs_bamcoverage/matrix/sp8_spe_pool_peaks_in_dlx5.gz --referencePoint center 
# 
# plotHeatmap -m /home/castillaa/chip/bigwigs_bamcoverage/matrix/sp8_spe_pool_peaks_in_dlx5.gz \
# -o /home/castillaa/chip/heatmaps/sp8_spe_pool_peaks_in_dlx5.pdf --startLabel center \
# --colorMap  bwr --legendLocation none 
# 
# 
# bedtools intersect -a DLX5_HL_11_merge_p001_peaks.narrowPeak \
# -b /home/castillaa/chip/revisions/input/sp8_bpe_for_great_s.bed \
# > intersect_dlx5_hl_sp8_spe.bed
# awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2, $3, $4}' /home/castillaa/chip/revisions/input/intersect_dlx5_hl_sp8_spe.bed > /home/castillaa/chip/revisions/input/intersect_dlx5_hl_sp8_spe.bed_4c.bed
# 
# # --zMax 7.5 7.5 4 4 
# #####correlation--------------------------------------
# 
# multiBamSummary bins --bamfiles /home/castillaa/chip/IDR/idr_redo/sp8/rochp_sp8_rep1_sorted_filt_chr.bam \
# /home/castillaa/chip/IDR/idr_redo/sp8/rochp_sp8_rep2_sorted_filt_chr.bam \
# /home/castillaa/chip/IDR/idr_redo/sp8/roinp_sp8_rep1_sorted_filt_chr.bam \
# /home/castillaa/chip/IDR/idr_redo/sp8/roinp_sp8_rep2_sorted_filt_chr.bam \
# -o /home/castillaa/chip/IDR/idr_redo/sp8/sp8_inp_multibamsummary.npz -p 10
# 
# plotCorrelation -in /home/castillaa/chip/IDR/idr_redo/sp8/sp8_inp_multibamsummary.npz \
# -c pearson -p scatterplot --removeOutliers --log1p \
# --labels Sp8_rep1 Sp8_rep2 Sp8_inp_rep1 Sp8_inp_rep2 \
# -o /home/castillaa/chip/heatmaps/plotCorrelation_sp8_inp_pearson.png
# 
# plotCorrelation -in /home/castillaa/chip/IDR/idr_redo/sp8/sp8_inp_multibamsummary.npz \
# -c spearman -p scatterplot --removeOutliers --log1p \
# --labels Sp8_rep1 Sp8_rep2 Sp8_inp_rep1 Sp8_inp_rep2 \
# -o /home/castillaa/chip/heatmaps/plotCorrelation_sp8_inp_spearman.png
# 
# 
# multiBamSummary bins --bamfiles /home/castillaa/chip/IDR/idr_redo/sp6/chipd_sp6_rep1_sorted_filt_chr.bam \
# /home/castillaa/chip/IDR/idr_redo/sp6/chipd_sp6_rep2_sorted_filt_chr.bam \
# /home/castillaa/chip/IDR/idr_redo/sp6/chipd_sp6_rep1_sorted_filt_chr.bam \
# /home/castillaa/chip/IDR/idr_redo/sp6/chipd_sp6_rep2_sorted_filt_chr.bam \
# -o /home/castillaa/chip/IDR/idr_redo/sp6/sp6_inp_multibamsummary.npz -p 10
# 
# plotCorrelation -in /home/castillaa/chip/IDR/idr_redo/sp6/sp6_inp_multibamsummary.npz \
# -c pearson -p scatterplot --removeOutliers --log1p \
# --labels sp6_rep1 sp6_rep2 sp6_inp_rep1 sp6_inp_rep2 \
# -o /home/castillaa/chip/heatmaps/plotCorrelation_sp6_inp_pearson.png
# 
# plotCorrelation -in /home/castillaa/chip/IDR/idr_redo/sp6/sp6_inp_multibamsummary.npz \
# -c spearman -p scatterplot --removeOutliers --log1p \
# --labels sp6_rep1 sp6_rep2 sp6_inp_rep1 sp6_inp_rep2 \
# -o /home/castillaa/chip/heatmaps/plotCorrelation_sp6_inp_spearman.png
# 
# 
# ##### d heatmaps of common and specific ######
# 
# computeMatrix reference-point -S /home/castillaa/chip/bigwigs_bamcoverage/rochp_sp8_rep1_sorted_filt_chr.bw \
# /home/castillaa/chip/bigwigs_bamcoverage/rochp_sp8_rep2_sorted_filt_chr.bw  \
# /home/castillaa/chip/bigwigs_bamcoverage/chipd_sp6_rep1_sorted_filt_chr.bw \
# /home/castillaa/chip/bigwigs_bamcoverage/chipd_sp6_rep2_sorted_filt_chr.bw  \
# -R /home/castillaa/chip/revisions/input/common_with_gb_for_great_s.bed \
# /home/castillaa/chip/revisions/input/sp8_bpe_for_great_s.bed \
# /home/castillaa/chip/revisions/input/sp6_bpe_for_great_s.bed \
# -a 3000 -b 3000 -p 16 \
# -o /home/castillaa/chip/bigwigs_bamcoverage/matrix/common_specif_pool_peaks_in_sp6-sp8.gz --referencePoint center 
# 
# plotHeatmap -m /home/castillaa/chip/bigwigs_bamcoverage/matrix/common_specif_pool_peaks_in_sp6-sp8.gz \
# -o /home/castillaa/chip/revisions/output/common_sp6-8_pool_peaks_in_sp6-8.pdf --startLabel center \
# --colorMap  bwr --legendLocation none  --zMax 3.5 3.5 5 5
# 
# #atac----------------
# computeMatrix reference-point -S /home/castillaa/chip/desanlis_2020_atac_ect/GSM4323465_WT.8_ectoderm_treat_pileup.bw \
# -R /home/castillaa/chip/revisions/input/common_with_gb_for_great_s.bed \
# /home/castillaa/chip/revisions/input/sp8_bpe_for_great_s.bed \
# /home/castillaa/chip/revisions/input/sp6_bpe_for_great_s.bed \
# -a 3000 -b 3000 -p 16 \
# -o /home/castillaa/chip/bigwigs_bamcoverage/matrix/kmeans_pool_peaks_in_atac.gz --referencePoint center 
# 
# plotHeatmap -m /home/castillaa/chip/bigwigs_bamcoverage/matrix/kmeans_pool_peaks_in_atac.gz \
# -o /home/castillaa/chip/revisions/output/kmeans_pool_peaks_in_atac.pdf --startLabel center \
# --colorList  blue,white,red --legendLocation none --zMax  1 
# 
# #me2----------------
# computeMatrix reference-point -S /home/castillaa/chip/bigwigs_bamcoverage/rochp_me2_rep1_sorted_filt_chr.bw \
# -R /home/castillaa/chip/revisions/input/common_with_gb_for_great_s.bed \
# /home/castillaa/chip/revisions/input/sp8_bpe_for_great_s.bed \
# /home/castillaa/chip/revisions/input/sp6_bpe_for_great_s.bed \
# -a 3000 -b 3000 -p 16 \
# -o /home/castillaa/chip/bigwigs_bamcoverage/matrix/kmeans_pool_peaks_in_me2.gz --referencePoint center 
# 
# plotHeatmap -m /home/castillaa/chip/bigwigs_bamcoverage/matrix/kmeans_pool_peaks_in_me2.gz \
# -o /home/castillaa/chip/revisions/output/kmeans_pool_peaks_in_me2.pdf --startLabel center \
# --colorList  blue,white,red --legendLocation none --zMax  4 
# 
# 
# #k27ac----------------
# computeMatrix reference-point -S /home/castillaa/chip/bigwigs_bamcoverage/chipd_ach_rep1_sorted_filt_chr_norm.bw \
# /home/castillaa/chip/bigwigs_bamcoverage/chipd_ach_rep2_sorted_filt_chr_norm.bw \
# /home/castillaa/chip/bigwigs_bamcoverage/chipd_acf_rep2_sorted_filt_chr_norm.bw \
# -R /home/castillaa/chip/revisions/input/common_with_gb_for_great_s.bed \
# /home/castillaa/chip/revisions/input/sp8_bpe_for_great_s.bed \
# /home/castillaa/chip/revisions/input/sp6_bpe_for_great_s.bed \
# -a 3000 -b 3000 -p 16 \
# -o /home/castillaa/chip/bigwigs_bamcoverage/matrix/kmeans_pool_peaks_in_k27ac.gz --referencePoint center 
# 
# plotHeatmap -m /home/castillaa/chip/bigwigs_bamcoverage/matrix/kmeans_pool_peaks_in_k27ac.gz \
# -o /home/castillaa/chip/revisions/output/kmeans_pool_peaks_in_k27ac.pdf --startLabel center \
# --colorList  blue,white,red --legendLocation none --zMax 4.5 4.5 9
# 


####################################### FIGURE 2:RNASEQ 3 MUTANTS AND dtgs each TF separately #################################

############## rna-seq SP6 SP8 DEGS ###########################

countDa_sp8 <- as.matrix(read.csv("input/gene_count_matrix_no_c_sp8.csv", row.names="gene_id"))
countDa_sp8 <- countDa_sp8[,c(4,5,6,1,2,3)] #### I NEEDED TO FIX THIS WAY BECAUSE THE WT AND KO WHERE SWAPPED
colData_sp8 <- read.csv("input/sampleplan_sp8.csv", sep=",", row.names=1)
countDa_sp6 <- as.matrix(read.csv("input/gene_count_matrix_no_c_sp6.csv", row.names="gene_id"))
colData_sp6 <- read.csv("input/sampleplan_sp6.csv", sep=",", row.names=1)

colnames(countDa_sp8) <- as.matrix(c("sp8_KO_rep1", "sp8_KO_rep2", "sp8_KO_rep3", "sp8_WT_rep1", "sp8_WT_rep2", "sp8_WT_rep3"))
countDa_sp8 <- as.data.frame(countDa_sp8)

colnames(countDa_sp6) <- gsub("RNASEQ_ect_", "", colnames(countDa_sp6))
countDa_sp6 <- as.data.frame(countDa_sp6)

#Check all sample IDs in colData are also in CountData and match their orders

all(rownames(colData_sp8) %in% colnames(countDa_sp8))

countDa_sp8 <- countDa_sp8[, rownames(colData_sp8)]
all(rownames(colData_sp8) == colnames(countDa_sp8))

all(rownames(colData_sp6) %in% colnames(countDa_sp6))

countDa_sp6 <- countDa_sp6[, rownames(colData_sp6)]
all(rownames(colData_sp6) == colnames(countDa_sp6))

### prepare deseq2 object and run deseq2

dds_sp8 <- DESeqDataSetFromMatrix(countData = na.omit(countDa_sp8),
                                  colData = colData_sp8, design = ~ genotype)

dds_sp8 <- DESeq(dds_sp8)
res_sp8 <- results(dds_sp8)
(resOrdered_sp8 <- as.data.frame(res_sp8[order(res_sp8$padj), ]))

#transfer rownames into columns
resOrdered_sp8$names <- rownames(resOrdered_sp8)
resOrdered_sp8 <- resOrdered_sp8 %>% 
  separate(names, into = c("Gene.ID", "Gene.name"), sep = "\\|")
rownames(resOrdered_sp8) <- NULL

### prepare deseq2 object and run deseq2
dds_sp6 <- DESeqDataSetFromMatrix(countData = countDa_sp6,
                                  colData = colData_sp6, design = ~ genotype)

dds_sp6 <- DESeq(dds_sp6)
res_sp6 <- results(dds_sp6)
(resOrdered_sp6 <- as.data.frame(res_sp6[order(res_sp6$padj), ]))

#transfer rownames into columns
resOrdered_sp6$names <- rownames(resOrdered_sp6)
resOrdered_sp6 <- resOrdered_sp6 %>% separate(names, into = c("Gene.ID", "Gene.name"), sep = "\\|")
rownames(resOrdered_sp6) <- NULL


### QC SP8

# Apply rlog transformation
vst_sp8 <- vst(dds_sp8, blind = T)

#pca
vst_sp8_assay <- as.data.frame(assay(vst_sp8))

# Perform PCA on the gene expression matrix
gene_pca_sp8 <- prcomp(t(vst_sp8_assay), scale. = F)

# Create a dataframe with the principal components
pca_df_sp8 <- as.data.frame(gene_pca_sp8$x)

# Add sample names or conditions for coloring if available
# Here, we're just using sample names
pca_df_sp8$Sample <- rownames(pca_df_sp8)
pca_df_sp8$group <- c("sp8_KO","sp8_KO","sp8_KO",
                      "sp8_WT","sp8_WT","sp8_WT")

# Plot PC1 vs PC2
ggplot(pca_df_sp8, aes(x = PC1, y = PC2, label = Sample, color=group)) +
  geom_point(size=5) +
  geom_text_repel(aes(label=Sample), size=5) +
  labs(title = "PCA: PC1 vs PC2",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  ylim(c(-20,20))

# ggsave("revisions/Plots/PCA_sp8_1.png")
ggsave("revisions/Plots/PCA_sp8_1_prev.png")

# Extract the loadings (rotation) matrix
loadings_sp8 <- gene_pca_sp8$rotation

# Determine the top 100 genes that impact PC1 the most
top_100_genes_pc1_sp8 <- head(order(abs(loadings_sp8[, 1]), decreasing = TRUE), 250)

# Display the top 100 genes
top_genes_pc1_sp8 <- rownames(loadings_sp8)[top_100_genes_pc1_sp8]
top_genes_pc1_sp8

# If you are interested in other principal components, e.g., PC2
top_100_genes_pc2_sp8 <- head(order(abs(loadings_sp8[, 2]), decreasing = TRUE), 250)
top_genes_pc2_sp8 <- rownames(loadings_sp8)[top_100_genes_pc2_sp8]
top_genes_pc2_sp8

# Get the names of the top 100 genes
top_genes_pc1_sp8 <- rownames(loadings_sp8)[top_100_genes_pc1_sp8]
top_genes_pc2_sp8 <- rownames(loadings_sp8)[top_100_genes_pc2_sp8]
# Subset the original gene expression matrix to include only the top 100 genes
top_genes_matrix_sp8 <- countDa_sp8[top_genes_pc1_sp8,c(1:6) ]
top_genes_matrix_pc2_sp8 <- countDa_sp8[top_genes_pc2_sp8,c(1:6) ]

# Create a heatmap of the top 100 genes
pdf("heatmaps/top250_counts_sp8_f.pdf",height = 30)
pheatmap(top_genes_matrix_sp8, scale = "row", clustering_distance_rows = "correlation"
         , show_rownames = TRUE, show_colnames = TRUE,
         main = "Heatmap of Top 250 Genes Impacting PC1 _sp8",cluster_cols = F)
dev.off()

############## QC SP6 ###########

vst_sp6 <- vst(dds_sp6, blind=T)

# Plot PCA
plotPCA(vst_sp6, intgroup = "genotype")

#other pcs
vst_sp6_assay <- as.data.frame(assay(vst_sp6))

# Perform PCA on the gene expression matrix
gene_pca_sp6 <- prcomp(t(vst_sp6_assay), scale. = F)

# Create a dataframe with the principal components
pca_df_sp6 <- as.data.frame(gene_pca_sp6$x)

# Add sample names or conditions for coloring if available
# Here, we're just using sample names
pca_df_sp6$Sample <- rownames(pca_df_sp6)
pca_df_sp6$group <- c("sp6_KO",
                      "sp6_KO",
                      "sp6_KO",
                      "sp6_WT",
                      "sp6_WT",
                      "sp6_WT")

# Plot PC1 vs PC2
ggplot(pca_df_sp6, aes(x = PC1, y = PC2, label = Sample, color=group)) +
  geom_point(size=5) +
  geom_text_repel(aes(label=Sample), size=5) +
  labs(title = "PCA: PC1 vs PC2",
       x = "Principal Component 1",
       y = "Principal Component 2")  
# ggsave("revisions/Plots/PCA_sp6_1.png")
ggsave("revisions/Plots/PCA_sp6_1_after.png")

# Extract the loadings (rotation) matrix
loadings_sp6 <- gene_pca_sp6$rotation

# Determine the top 100 genes that impact PC1 the most
top_100_genes_pc1_sp6 <- head(order(abs(loadings_sp6[, 1]), decreasing = TRUE), 250)

# Display the top 100 genes
top_genes_pc1_sp6 <- rownames(loadings_sp6)[top_100_genes_pc1_sp6]
top_genes_pc1_sp6

# If you are interested in other principal components, e.g., PC2
top_100_genes_pc2_sp6 <- head(order(abs(loadings_sp6[, 2]), decreasing = TRUE), 250)
top_genes_pc2_sp6 <- rownames(loadings_sp6)[top_100_genes_pc2_sp6]
top_genes_pc2_sp6

# Get the names of the top 100 genes
top_genes_pc1_sp6 <- rownames(loadings_sp6)[top_100_genes_pc1_sp6]
top_genes_pc2_sp6 <- rownames(loadings_sp6)[top_100_genes_pc2_sp6]
# Subset the original gene expression matrix to include only the top 100 genes
top_genes_matrix_sp6 <- countDa_sp6[top_genes_pc1_sp6,c(1:6) ]
top_genes_matrix_pc2_sp6 <- countDa_sp6[top_genes_pc2_sp6,c(1:6) ]

# Create a heatmap of the top 100 genes
pdf("heatmaps/top250_counts_sp6_f_f.pdf",height = 30)
pheatmap(top_genes_matrix_sp6, scale = "row", clustering_distance_rows = "correlation"
         , show_rownames = TRUE, show_colnames = TRUE,
         main = "Heatmap of Top 250 Genes Impacting PC1 _sp6",cluster_cols = F)
dev.off()


################ POST-FILTERING ANALYSIS #############

## start fresh
countDa_sp8 <- as.matrix(read.csv("input/gene_count_matrix_no_c_sp8.csv", row.names="gene_id"))
countDa_sp8 <- countDa_sp8[,c(4,5,6,1,2,3)] #### I NEEDED TO FIX THIS WAY BECAUSE THE WT AND KO WHERE SWAPPED
colData_sp8 <- read.csv("input/sampleplan_sp8.csv", sep=",", row.names=1)

colnames(countDa_sp8) <- as.matrix(c("sp8_KO_rep1", "sp8_KO_rep2", "sp8_KO_rep3", "sp8_WT_rep1", "sp8_WT_rep2", "sp8_WT_rep3"))
countDa_sp8 <- as.data.frame(countDa_sp8)


## remove rep2 KO as it contains contamination
countDa_sp8 <- countDa_sp8[,c(1,3,4,5,6)]
colData_sp8 <- colData_sp8[c(1,3,4,5,6),]
#coldata 
# sample,genotype,replicate
# sp8_KO_rep1,KO,1
#
# sp8_KO_rep3,KO,3
# sp8_WT_rep1,WT,1
# sp8_WT_rep2,WT,2
# sp8_WT_rep3,WT,1

RNASEQ_ect_sp8_KO_rep1 <- read.table("input/RNASEQ_ect_sp8_KO_rep1.tab",sep = "\t", header=T)
RNASEQ_ect_sp8_KO_rep3 <- read.table("input/RNASEQ_ect_sp8_KO_rep3.tab",sep = "\t", header=T)
RNASEQ_ect_sp8_WT_rep1 <- read.table("input/RNASEQ_ect_sp8_WT_rep1.tab",sep = "\t", header=T)
RNASEQ_ect_sp8_WT_rep2 <- read.table("input/RNASEQ_ect_sp8_WT_rep2.tab",sep = "\t", header=T)
RNASEQ_ect_sp8_WT_rep3 <- read.table("input/RNASEQ_ect_sp8_WT_rep3.tab",sep = "\t", header=T)

RNASEQ_ect_sp8_KO_rep1 <- RNASEQ_ect_sp8_KO_rep1 %>% dplyr::rename(TPM_sp8_KO_rep1 = TPM)

RNASEQ_ect_sp8_KO_rep3 <- RNASEQ_ect_sp8_KO_rep3 %>% dplyr::rename(TPM_sp8_KO_rep3 = TPM)
RNASEQ_ect_sp8_WT_rep1 <- RNASEQ_ect_sp8_WT_rep1 %>% dplyr::rename(TPM_sp8_WT_rep1 = TPM)
RNASEQ_ect_sp8_WT_rep2 <- RNASEQ_ect_sp8_WT_rep2 %>% dplyr::rename(TPM_sp8_WT_rep2 = TPM)
RNASEQ_ect_sp8_WT_rep3 <- RNASEQ_ect_sp8_WT_rep3 %>% dplyr::rename(TPM_sp8_WT_rep3 = TPM)

TPM_sp8 <- RNASEQ_ect_sp8_KO_rep1 %>%
  
  inner_join(RNASEQ_ect_sp8_KO_rep3, by = c("Gene.ID", "Gene.Name")) %>%
  inner_join(RNASEQ_ect_sp8_WT_rep1, by = c("Gene.ID", "Gene.Name")) %>%
  inner_join(RNASEQ_ect_sp8_WT_rep2, by = c("Gene.ID", "Gene.Name")) %>%
  inner_join(RNASEQ_ect_sp8_WT_rep3, by = c("Gene.ID", "Gene.Name"))

TPM_sp8 <- TPM_sp8[,c(1,2,9,16,23,30,37)]


###simplyfy column names rep3 is now rep2 KO
colnames(countDa_sp8) <- as.matrix(c("sp8_KO_rep1", "sp8_KO_rep2", "sp8_WT_rep1", "sp8_WT_rep2", "sp8_WT_rep3"))
rownames(colData_sp8) <- as.matrix(c("sp8_KO_rep1", "sp8_KO_rep2", "sp8_WT_rep1", "sp8_WT_rep2", "sp8_WT_rep3"))
countDa_sp8 <- as.data.frame(countDa_sp8)

all(rownames(colData_sp8) %in% colnames(countDa_sp8))

countDa_sp8 <- countDa_sp8[, rownames(colData_sp8)]
all(rownames(colData_sp8) == colnames(countDa_sp8))

### prepare deseq2 object and run deseq2

dds_sp8 <- DESeqDataSetFromMatrix(countData = na.omit(countDa_sp8),
                                  colData = colData_sp8, design = ~ genotype)

dds_sp8 <- DESeq(dds_sp8)
res_sp8 <- results(dds_sp8)
resOrdered_sp8 <- as.data.frame(res_sp8[order(res_sp8$padj), ])
saveRDS(dds_sp8, "output/dds_sp8_filtered.RDS")

vst_sp8 <- vst(dds_sp8, blind = T)

##################
####pca sp8#######
##################
##put the names back again in the rownames
row.names(countDa_sp8) <- countDa_sp8$gene_id
countDa_sp8 <- countDa_sp8[,c(1:5)]

# Plot PCA
plotPCA(vst_sp8, intgroup = "genotype")

#other pcs
vst_sp8_assay <- as.data.frame(assay(vst_sp8))

# Perform PCA on the gene expression matrix
gene_pca_sp8 <- prcomp(t(vst_sp8_assay), scale. = F)

# Create a dataframe with the principal components
pca_df_sp8 <- as.data.frame(gene_pca_sp8$x)

# Add sample names or conditions for coloring if available
# Here, we're just using sample names
pca_df_sp8$Sample <- rownames(pca_df_sp8)
pca_df_sp8$group <- c("sp8_KO","sp8_KO",
                      "sp8_WT","sp8_WT","sp8_WT")

# Plot PC1 vs PC2
ggplot(pca_df_sp8, aes(x = PC1, y = PC2, label = Sample, color=group)) +
  geom_point(size=5) +
  geom_text_repel(aes(label=Sample), size=5) +
  labs(title = "PCA: PC1 vs PC2",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  ylim(c(-20,20))

# ggsave("revisions/Plots/PCA_sp8_1.png")
ggsave("revisions/Plots/PCA_sp8_1_after.png")

# Extract the loadings (rotation) matrix
loadings_sp8 <- gene_pca_sp8$rotation

# Determine the top 100 genes that impact PC1 the most
top_100_genes_pc1_sp8 <- head(order(abs(loadings_sp8[, 1]), decreasing = TRUE), 250)

# Display the top 100 genes
top_genes_pc1_sp8 <- rownames(loadings_sp8)[top_100_genes_pc1_sp8]

# Subset the original gene expression matrix to include only the top 100 genes
top_genes_matrix_sp8 <- countDa_sp8[top_genes_pc1_sp8,c(1:5) ]

# Create a heatmap of the top 100 genes
pdf("heatmaps/top250_counts_sp8_f.pdf",height = 30)
pheatmap(top_genes_matrix_sp8, scale = "row", clustering_distance_rows = "correlation"
         , show_rownames = TRUE, show_colnames = TRUE,
         main = "Heatmap of Top 250 Genes Impacting PC1 _sp8",cluster_cols = F)
dev.off()

#transfer rownames into columns
resOrdered_sp8$names <- rownames(resOrdered_sp8)
resOrdered_sp8 <- resOrdered_sp8 %>% separate(names, into = c("Gene.ID", "Gene.name"), sep = "\\|")
rownames(resOrdered_sp8) <- NULL

#merge with tpms
resord_sp8 <- merge(resOrdered_sp8, TPM_sp8, by = "Gene.ID", all.x=T)
resord_sp8 <- resord_sp8[,c(1, 8,2,3,4,5,6,7,10,11,12,13,14)]

###simplyfy column names rep3 is now rep2 KO
colnames(resord_sp8)[10] <- as.matrix(c("TPM_sp8_KO_rep2"))

# writexl::write_xlsx(resord_sp8, "output/results_sp8_filtered.xlsx")

# resord_sp8 <- readxl::read_xlsx("output/results_sp8_filtered.xlsx")

#### assign sp6 name
resord_sp8[28924,]$Gene.name <- "Sp6"

## calculate degs

sp8_deg <- subset(resord_sp8, subset = padj <= 0.05 & (log2FoldChange > 0.585 | log2FoldChange < -0.585))
sp8_deg <- sp8_deg[sp8_deg$Gene.name!="Sp8",] #remove sp8
sp8_deg <- sp8_deg %>%
  filter(Gene.ID != "ENSMUSG00000117694.1")
table(duplicated(sp8_deg$Gene.name))


##########################
########sp6 ##############
##########################

countDa_sp6 <- as.matrix(read.csv("input/gene_count_matrix_no_c_sp6.csv", row.names="gene_id"))
colData_sp6 <- read.csv("input/sampleplan_sp6.csv", sep=",", row.names=1)

colnames(countDa_sp6) <- gsub("RNASEQ_ect_", "", colnames(countDa_sp6))
countDa_sp6 <- as.data.frame(countDa_sp6)


#coldata 
# sample,genotype,replicate

# sp6_KO_rep2,KO,2
# sp6_KO_rep3,KO,3
# sp6_WT_rep1,WT,1
# sp6_WT_rep2,WT,2
# sp6_WT_rep3,WT,1

## remove rep1 KO and rep2 WT as it contains contamination
countDa_sp6 <- countDa_sp6[,c(2,
                              3,
                              4,
                              # 5,
                              6)]
colData_sp6 <- colData_sp6[c(2,
                             3,
                             4,
                             # 5,
                             6),]
colData_sp6[4,2] <- "3" #fix error


# countDa_sp6 <- countDa_sp6[-grep("rpl|mt-|Hbb|Hba", rownames(countDa_sp6), ignore.case = TRUE),]


RNASEQ_ect_sp6_KO_rep2 <- read.table("input/RNASEQ_ect_sp6_KO_rep2.tab",sep = "\t", header=T)
RNASEQ_ect_sp6_KO_rep3 <- read.table("input/RNASEQ_ect_sp6_KO_rep3.tab",sep = "\t", header=T)
RNASEQ_ect_sp6_WT_rep1 <- read.table("input/RNASEQ_ect_sp6_WT_rep1.tab",sep = "\t", header=T)
# RNASEQ_ect_sp6_WT_rep2 <- read.table("input/RNASEQ_ect_sp6_WT_rep2.tab",sep = "\t", header=T)
RNASEQ_ect_sp6_WT_rep3 <- read.table("input/RNASEQ_ect_sp6_WT_rep3.tab",sep = "\t", header=T)



RNASEQ_ect_sp6_KO_rep2 <- RNASEQ_ect_sp6_KO_rep2 %>% dplyr::rename(TPM_sp6_KO_rep2 = TPM)
RNASEQ_ect_sp6_KO_rep3 <- RNASEQ_ect_sp6_KO_rep3 %>% dplyr::rename(TPM_sp6_KO_rep3 = TPM)
RNASEQ_ect_sp6_WT_rep1 <- RNASEQ_ect_sp6_WT_rep1 %>% dplyr::rename(TPM_sp6_WT_rep1 = TPM)
# RNASEQ_ect_sp6_WT_rep2 <- RNASEQ_ect_sp6_WT_rep2 %>% dplyr::rename(TPM_sp6_WT_rep2 = TPM)
RNASEQ_ect_sp6_WT_rep3 <- RNASEQ_ect_sp6_WT_rep3 %>% dplyr::rename(TPM_sp6_WT_rep3 = TPM)


TPM_sp6 <- RNASEQ_ect_sp6_KO_rep2 %>%
  inner_join(RNASEQ_ect_sp6_KO_rep3, by = c("Gene.ID", "Gene.Name")) %>%
  inner_join(RNASEQ_ect_sp6_WT_rep1, by = c("Gene.ID", "Gene.Name")) %>%
  # inner_join(RNASEQ_ect_sp6_WT_rep2, by = c("Gene.ID", "Gene.Name")) %>%
  inner_join(RNASEQ_ect_sp6_WT_rep3, by = c("Gene.ID", "Gene.Name"))


TPM_sp6 <- TPM_sp6[,c(1,2,9,16,23,30
                      # ,37
)]


colnames(countDa_sp6) <- gsub("RNASEQ_ect_", "", colnames(countDa_sp6))
countDa_sp6 <- as.data.frame(countDa_sp6)

dds_sp6 <- DESeqDataSetFromMatrix(countData = countDa_sp6,
                                  colData = colData_sp6, design = ~ genotype)

dds_sp6 <- DESeq(dds_sp6)

res_sp6 <- results(dds_sp6)

resOrdered_sp6 <- as.data.frame(res_sp6[order(res_sp6$padj), ])

saveRDS(dds_sp6,"output/dds_sp6_filtered.RDS")

#transfer rownames into columns
resOrdered_sp6$names <- rownames(resOrdered_sp6)
resOrdered_sp6 <- resOrdered_sp6 %>% separate(names, into = c("Gene.ID", "Gene.name"), sep = "\\|")
rownames(resOrdered_sp6) <- NULL

#merge with tpms
resord_sp6 <- merge(resOrdered_sp6, TPM_sp6, by = "Gene.ID", all.x=T)
resord_sp6 <- resord_sp6[,c(1, 8,2,3,4,5,6,7,10,11,12,13
                            # ,14
)]
###simplyfy column names rep3 is now rep2 KO
colnames(resord_sp6)[9]  <- as.matrix(c("TPM_sp6_KO_rep1"))
colnames(resord_sp6)[10] <- as.matrix(c("TPM_sp6_KO_rep2"))
colnames(resord_sp6)[11] <- as.matrix(c("TPM_sp6_WT_rep1"))
colnames(resord_sp6)[12] <- as.matrix(c("TPM_sp6_WT_rep2"))
## calculate degs 

sp6_deg <- subset(resord_sp6, subset = padj <= 0.05 & (log2FoldChange > 0.585 | log2FoldChange < -0.585))
sp6_deg <- sp6_deg[sp6_deg$Gene.name!="Sp6",] #remove sp6

table(duplicated(sp6_deg$Gene.name))
# writexl::write_xlsx(resord_sp6, "output/results_sp6_filtered.xlsx")
# 
# write.csv(sp6_deg_pval, "output/sp6_degs_pval_for_david.txt")





########################################
############ DESeq2 - DKO  #############
########################################
countDa <- as.matrix(read.csv("revisions/input/gene_count_matrix_DKO.csv", row.names = 1))

# Sum duplicate rows (same gene_id)
countDa <- as.data.frame(countDa)
countDa$gene_id <- rownames(countDa)

countDa <- countDa %>%
  group_by(gene_id) %>%
  summarise(across(where(is.numeric), sum)) %>%
  as.data.frame()

rownames(countDa) <- countDa$gene_id
countDa$gene_id <- NULL
countDa <- as.matrix(countDa)

cat("Genes after collapsing duplicates:", nrow(countDa), "\n")

# Sample order: WT(3),  DKO(3)
# colData: one row per sample
colData <- data.frame(
  row.names = c("RNA_WT_32", "RNA_WT_53", "RNA_WT_56",
                "RNA_DKO_17", "RNA_DKO_40", "RNA_DKO_41"),
  genotype = factor(c("WT","WT","WT",
                      "DKO","DKO","DKO"),
                    levels = c("WT","Sp6","Sp8","DKO"))
)

# Verify samples match
all(rownames(colData) %in% colnames(countDa))
countDa <- countDa[, rownames(colData)]
all(rownames(colData) == colnames(countDa))

########################################
##### Pre-filtering low-count genes ####
########################################

# Keep genes with at least 5 counts in at least N samples
# N = size of smallest group (Sp8 has 2, so use 2 as minimum)
min_counts <- 5
min_samples <- 2

keep <- rowSums(countDa >= min_counts) >= min_samples
countDa <- countDa[keep, ]

cat("Genes before filtering:", nrow(countDa) + sum(!keep), "\n")
cat("Genes after filtering:", nrow(countDa), "\n")

########################################
##### Load and merge TPM tables ########
########################################

# Genes to exclude (high-variance contaminants / non-informative)
exclude_genes <- c("Hbb-bs", "Hbb-bt", "Hbb-bh1", "Hbb-bh2", "Hbb-y",
                   "Hba-a1", "Hba-a2", "Hba-x")  # full globin family

samples <- c("RNA_WT_32", "RNA_WT_53", "RNA_WT_56",
             "RNA_DKO_17", "RNA_DKO_40", "RNA_DKO_41")

tpm_list <- lapply(samples, function(s) {
  df <- read.table(paste0("revisions/input/", s, ".tab"), sep = "\t", header = TRUE)
  
  # Collapse duplicate Gene.ID entries by summing TPM
  df <- df %>%
    group_by(Gene.ID, Gene.Name) %>%
    summarise(TPM = sum(TPM), .groups = "drop")
  
  # Remove globin genes
  # df <- df %>% filter(!Gene.Name %in% exclude_genes)
  
  df <- df %>% dplyr::rename(!!paste0("TPM_", s) := TPM)
  df
})

TPM_all <- tpm_list[[1]]
for (i in 2:length(tpm_list)) {
  TPM_all <- inner_join(TPM_all, tpm_list[[i]], by = c("Gene.ID", "Gene.Name"))
}
TPM_all <- TPM_all[, c("Gene.ID", "Gene.Name", paste0("TPM_", samples))]

########################################
###### Helper: format DESeq2 results ###
########################################

format_results <- function(res, TPM) {
  df <- as.data.frame(res[order(res$padj), ])
  df$names <- rownames(df)
  df <- df %>% separate(names, into = c("Gene.ID", "Gene.name"), sep = "\\|")
  rownames(df) <- NULL
  df <- merge(df, TPM, by = "Gene.ID", all.x = TRUE)
  df
}


########################################
###### pairwise analysis ######
########################################

run_pairwise <- function(count_matrix, coldata_full, groups, label) {
  keep <- rownames(coldata_full)[coldata_full$genotype %in% groups]
  cd <- droplevels(coldata_full[keep, , drop = FALSE])
  cm <- count_matrix[, keep]
  keep_genes <- rowSums(cm >= 10) >= 2
  cm <- cm[keep_genes, ]
  dds <- DESeqDataSetFromMatrix(countData = na.omit(cm),
                                colData = cd,
                                design = ~ genotype)
  ko_group <- groups[groups != "WT"]
  dds$genotype <- relevel(dds$genotype, ref = ko_group)
  dds <- DESeq(dds)
  saveRDS(dds, paste0("output/dds_", label, ".RDS"))
  res <- results(dds)
  list(dds = dds, res = res)
}


# WT vs DKO
pw_DKO <- run_pairwise(countDa, colData, c("WT","DKO"), "WTvsDKO")
resord_DKO <- format_results(pw_DKO$res, TPM_all)
writexl::write_xlsx(resord_DKO, "output/results_WTvsDKO.xlsx")
deg_DKO      <- subset(resord_DKO, padj <= 0.05 & (log2FoldChange > 0.585 | log2FoldChange < -0.585))
deg_DKO <- deg_DKO[deg_DKO$Gene.name!="Sp8",] #remove sp8

# writexl::write_xlsx(deg_DKO, "output/deg_DKO.xlsx")
########################################
########## PCA - all samples ###########
########################################

vst_all <- vst(dds_all, blind = TRUE)
vst_assay <- as.data.frame(assay(vst_all))

gene_pca <- prcomp(t(vst_assay), scale. = FALSE)
pca_df <- as.data.frame(gene_pca$x)
pca_df$Sample <- rownames(pca_df)
pca_df$group <- colData[pca_df$Sample, "genotype"]
pct_var <- setNames(round(100 * gene_pca$sdev^2 / sum(gene_pca$sdev^2), 1), paste0("PC", seq_along(gene_pca$sdev)))

plot_pca <- function(df, pcx, pcy, pct_var = NULL) {
  xlab <- if (!is.null(pct_var)) paste0(pcx, " (", pct_var[pcx], "%)") else pcx
  ylab <- if (!is.null(pct_var)) paste0(pcy, " (", pct_var[pcy], "%)") else pcy
  ggplot(df, aes_string(x = pcx, y = pcy, color = "group", label = "Sample")) +
    geom_point(size = 5) +
    geom_text_repel(size = 4) +
    labs(title = paste("PCA:", pcx, "vs", pcy),
         x = xlab, y = ylab) +
    theme_minimal()
}

plot_pca(pca_df, "PC1", "PC2", pct_var); ggsave("revisions/Plots/PCA_DKO_PC1vsPC2.png")


########################################
###### Heatmaps - top PC genes #########
########################################

loadings <- gene_pca$rotation

top_genes_pc <- function(loadings, pc, n = 250) {
  rownames(loadings)[head(order(abs(loadings[, pc]), decreasing = TRUE), n)]
}

top_pc1 <- top_genes_pc(loadings, 1)

# Use VST matrix for heatmaps (more appropriate than raw counts)
vst_mat <- as.data.frame(assay(vst_all))

pdf("heatmaps/top250_vst_PC1_DKO.pdf", height = 30)
pheatmap(vst_mat[top_pc1, ], scale = "row",
         clustering_distance_rows = "correlation",
         show_rownames = TRUE, show_colnames = TRUE,
         main = "Top 250 genes impacting PC1 - DKO project",
         cluster_cols = FALSE)
dev.off()



##########VOLCANO PLOTS ###############

#rna-seq degs only significant

Genes_sp8 <- unique(c("Fgf8", "Rspo2","En1","Wnt7a", "Msx1","Msx2","Fzd1",
                      "Sost", "Sp6", "Lgr6","Crabp1", "Cyp26a1", "Dlx1"
                      # , "Bambi","Bmp4","Hoxc9","Crabp1",
))



resord_sp8$Category <- with(resord_sp8, ifelse(padj < 0.05 & log2FoldChange > 0.585, "Positive",
                                               ifelse(padj < 0.05 & log2FoldChange < -0.585, "Negative", "Neutral")))

# Volcano plot with border for labeled genes
ggplot(resord_sp8, aes(x = log2FoldChange, y = -log10(padj))) +
  
  # Base layer: all points
  geom_point(aes(color = Category), size = 1.5) +
  
  # Outline points for labeled genes (overlay with larger black-bordered shape)
  geom_point(data = subset(resord_sp8, Gene.name %in% Genes_sp8),
             shape = 21, stroke = 0.8, color = "black", 
             aes(fill = Category), size = 2.5) +
  
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey80") +
  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed", color = "grey80") +
  # Labels
  geom_text_repel(data  = subset(resord_sp8, Gene.name %in% Genes_sp8),
                  aes(label = Gene.name),
                  box.padding = 0.3,
                  point.padding = 0.3,
                  segment.color = 'black',
                  segment.size = 0.75,
                  max.overlaps = Inf,
                  size = 7.5,
                  min.segment.length = 0,
                  seed = 12) +
  # Axes limits
  xlim(c(-10, 10)) +
  # ylim(c(0, 125)) +
  # Colors
  scale_color_manual(values = c("Positive" = "orange1", 
                                "Negative" = "dodgerblue", 
                                "Neutral" = "darkgrey")) +
  scale_fill_manual(values = c("Positive" = "orange1", 
                               "Negative" = "dodgerblue", 
                               "Neutral" = "darkgrey")) +
  # Theme
  theme_classic(base_size = 16) +
  labs(x = "Log2 Fold Change",
       y = "-Log10(P-value)",
       title = NULL) +
  theme(legend.position = "none")
ggsave("Revisions/Plots/volcano_plot_sp8.png",units = "cm",width = 15,height = 15, dpi = 300)



#broken volcano plot 
Genes_sp6 <- unique(c(
  # "Fgf8", "Rspo2","En1","Wnt7a", "Msx1","Msx2","Fzd1",
  "Sost", "Sp6", "Lgr6","Crabp1", "Cyp26a1", "Dlx1"))

resord_sp6$Category <- with(resord_sp6, ifelse(padj < 0.05 & log2FoldChange > 0.585, "Positive",
                                               ifelse(padj < 0.05 & log2FoldChange < -0.585, "Negative", "Neutral")))

# volcano2a <- 
ggplot(resord_sp6, aes(x = log2FoldChange, y = -log10(padj))) +
  
  # Base layer: all points
  geom_point(aes(color = Category), size = 1.5) +
  
  # Outline points for labeled genes (overlay with larger black-bordered shape)
  geom_point(data = subset(resord_sp6, Gene.name %in% Genes_sp6),
             shape = 21, stroke = 0.8, color = "black", 
             aes(fill = Category), size = 2.5) +
  
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey80") +
  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed", color = "grey80") +
  # Labels
  geom_text_repel(data  = subset(resord_sp6, Gene.name %in% Genes_sp6),
                  aes(label = Gene.name),
                  box.padding = 0.3,
                  point.padding = 0.3,
                  segment.color = 'black',
                  segment.size = 0.75,
                  max.overlaps = Inf,
                  size = 7.5,
                  min.segment.length = 0,
                  seed = 12) +
  # Axes limits
  xlim(c(-10, 10)) +
  ylim(c(0, 25)) +
  # Colors
  scale_color_manual(values = c("Positive" = "orange1", 
                                "Negative" = "dodgerblue", 
                                "Neutral" = "darkgrey")) +
  scale_fill_manual(values = c("Positive" = "orange1", 
                               "Negative" = "dodgerblue", 
                               "Neutral" = "darkgrey")) +
  # Theme
  theme_classic(base_size = 16) +
  labs(x = "Log2 Fold Change",
       y = "-Log10(P-value)",
       title = NULL) +
  theme(legend.position = "none")
ggsave("Revisions/Plots/volcano_plot_sp6.png",dpi=300, unit="cm",width = 15,height = 15)


Genes_sp8 <- unique(c("Fgf8", "Rspo2","En1","Msx1","Fzd1",
                      "Sost",  "Lgr6", "Cyp26a1", "Dlx1"
                      # , "Bambi","Bmp4","Hoxc9","Crabp1",
))
resord_DKO$Category <- with(resord_DKO, ifelse(padj < 0.05 & log2FoldChange > 0.585, "Positive",
                                               ifelse(padj < 0.05 & log2FoldChange < -0.585, "Negative", "Neutral")))
ggplot(resord_DKO, aes(x = log2FoldChange, y = -log10(padj))) +
  
  # Base layer: all points
  geom_point(aes(color = Category), size = 1.5) +
  
  # Outline points for labeled genes (overlay with larger black-bordered shape)
  geom_point(data = subset(resord_DKO, Gene.name %in% Genes_sp8),
             shape = 21, stroke = 0.8, color = "black", 
             aes(fill = Category), size = 2.5) +
  
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey80") +
  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed", color = "grey80") +
  # Labels
  geom_text_repel(data  = subset(resord_DKO, Gene.name %in% Genes_sp8),
                  aes(label = Gene.name),
                  box.padding = 0.3,
                  point.padding = 0.3,
                  segment.color = 'black',
                  segment.size = 0.75,
                  max.overlaps = Inf,
                  size = 7.5,
                  min.segment.length = 0,
                  seed = 12) +
  # Axes limits
  xlim(c(-10, 10)) +
  ylim(c(0, 75)) +
  # Colors
  scale_color_manual(values = c("Positive" = "orange1", 
                                "Negative" = "dodgerblue", 
                                "Neutral" = "darkgrey")) +
  scale_fill_manual(values = c("Positive" = "orange1", 
                               "Negative" = "dodgerblue", 
                               "Neutral" = "darkgrey")) +
  # Theme
  theme_classic(base_size = 16) +
  labs(x = "Log2 Fold Change",
       y = "-Log10(P-value)",
       title = NULL) +
  theme(legend.position = "none")
ggsave("Revisions/Plots/volcano_plot_DKO.png",dpi=300, unit="cm",width = 15,height = 15)


#--------------------------
############# SPIA FOR SPS 
#-------------------------

# we need to make a database for SPIA with all pathways:
# At the time of analysis, KEGG pathways were obtained by:
# curl -s http://rest.kegg.jp/list/pathway/mmu | awk '{split($1,a,":"); print "curl -s http://rest.kegg.jp/get/"a[2]"/kgml -o extdata/keggxml/mmu/"a[2]".xml"}' | bash
# They are now in input/KEGGPATHWAYS

library(SPIA)
library(org.Mm.eg.db)
library(msigdbr)
library(KEGGREST)

makeSPIAdata(kgml.path="input/KEGGPATHWAYS/", organism="mmu", out.path=system.file("extdata/",package="SPIA"))

# Load data
# Map SYMBOL → ENTREZ
resord_sp8$ENTREZ <- mapIds(org.Mm.eg.db,
                            keys = resord_sp8$Gene.name,
                            column = "ENTREZID",
                            keytype = "SYMBOL",
                            multiVals = "first")

# Remove duplicates and NAs
resord_sp8 <- resord_sp8[!duplicated(resord_sp8$ENTREZ), ]
resord_sp8 <- resord_sp8[!is.na(resord_sp8$ENTREZ), ]

# Define DE gene list
res_filt_sp8_pval <- subset(resord_sp8,
                            padj <= 0.05 & (log2FoldChange > 0.585 | log2FoldChange < -0.585))

input_sp8 <- sp8_deg$log2FoldChange
names(input_sp8) <- res_filt_sp8_pval$ENTREZ

all_genes <- res_filt_sp8$ENTREZ  # now proper Entrez IDs

# Run SPIA
# res <- spia(de = input_sp8, all = all_genes, organism = "mmu",
#             nB = 2000, plots = TRUE, combine = "fisher")
# 
# 
writexl::write_xlsx(res, "revisions/output/SPIA_sp8_raw_results.xlsx")
res_spia_sp8 <- readxl::read_excel("output/SPIA_sp8_raw_results.xlsx")

# Pathways of interest

pathways_to_check <- c("04310", "04350", "04010", 
                       "04512", 
                       # "04510","04530", 
                       # "04540",
                       "04330"
                       # "04015",
                       # "04810"
)
names(pathways_to_check) <- c("Wnt", "TGFB","Mapk" ,
                              "ECM-receptor interaction", 
                              # "Focal adhesion", "Tight junction", 
                              # "Gap junction",
                              "Notch"
                              # "Rap1 signaling pathway",
                              # "04180"
)


res_spia_sp8$status <- "NS"
res_spia_sp8$status[res_spia_sp8$pG < 0.05 & res_spia_sp8$tA > 0] <- "Inhibited in KO"
res_spia_sp8$status[res_spia_sp8$pG < 0.05 & res_spia_sp8$tA < 0] <- "Activated in KO"
fixedColors <- list('status'=c('Activated in KO'="dodgerblue",'Inhibited in KO'="orange1",'NS'="grey"))
res_spia_sp8$custom_label <- names(pathways_to_check)[match(res_spia_sp8$ID, pathways_to_check)]
ggplot(res_spia_sp8, aes(x = tA, y = -log(pG))) +
  geom_point(aes(color = status)) +
  geom_label_repel(data = subset(res_spia_sp8, ID %in% pathways_to_check),
                   aes(label = custom_label,
                       fill = status),        # colored box background
                   color = "black",           # black text + border
                   min.segment.length = 0,
                   show.legend = F,
                   box.padding = 0.6,
                   point.padding = 0.4,
                   force = 10,
                   alpha = 0.8,              # slight transparency so fill isn't too heavy
                   size = 6) +
  xlim(c(-100, 100)) +
  ylim(c(0, 30)) +
  scale_color_manual("", values = fixedColors[["status"]]) +
  scale_fill_manual(values = fixedColors[["status"]], guide = "none") +
  theme_classic(base_size = 16) +
  geom_hline(yintercept = -log(0.05), linetype = "dashed", color = "grey80") +
  theme(legend.position = "bottom")

ggsave("revisions/Plots/spia_sp8.png",units = "cm",width = 15,height = 15, dpi = 300)

#----------------------------------------------------------------------
#----------------------- sp6---------------------
#-------------------------------------------------------------------

# Map SYMBOL → ENTREZ
resord_sp6$ENTREZ <- mapIds(org.Mm.eg.db,
                            keys = resord_sp6$Gene.name,
                            column = "ENTREZID",
                            keytype = "SYMBOL",
                            multiVals = "first")

# Remove duplicates and NAs
resord_sp6 <- resord_sp6[!duplicated(resord_sp6$ENTREZ), ]
resord_sp6 <- resord_sp6[!is.na(resord_sp6$ENTREZ), ]

# Define DE gene list
res_filt_sp6_pval <- subset(resord_sp6,
                            padj <= 0.05 & (log2FoldChange > 0.585 | log2FoldChange < -0.585))

input_sp6 <- res_filt_sp6_pval$log2FoldChange
names(input_sp6) <- res_filt_sp6_pval$ENTREZ

all_genes <- resord_sp6$ENTREZ  # now proper Entrez IDs

# Run SPIA
# res_sp6 <- spia(de = input_sp6, all = all_genes, organism = "mmu",
#             nB = 2000, plots = TRUE, combine = "fisher")
# 
# 
writexl::write_xlsx(res_sp6, "output/SPIA_Sp6_raw_results.xlsx")
res_spia_sp6 <- readxl::read_excel("output/SPIA_sp6_raw_results.xlsx")



res_spia_sp6$status <- "NS"
res_spia_sp6$status[res_spia_sp6$pG < 0.05 & res_spia_sp6$tA > 0] <- "Inhibited in KO"
res_spia_sp6$status[res_spia_sp6$pG < 0.05 & res_spia_sp6$tA < 0] <- "Activated in KO"
fixedColors <- list('status'=c('Activated in KO'="dodgerblue",'Inhibited in KO'="orange2",'NS'="grey"))
res_spia_sp6$custom_label <- names(pathways_to_check)[match(res_spia_sp6$ID, pathways_to_check)]
ggplot(res_spia_sp8, aes(x = tA, y = -log(pG))) +
  geom_point(aes(color = status)) +
  geom_label_repel(data = subset(res_spia_sp6, ID %in% pathways_to_check),
                   aes(label = custom_label,
                       fill = status),        # colored box background
                   color = "black",           # black text + border
                   min.segment.length = 0,
                   show.legend = F,
                   box.padding = 0.6,
                   point.padding = 0.4,
                   force = 10,
                   alpha = 0.8,              # slight transparency so fill isn't too heavy
                   size = 6) +
  xlim(c(-100, 100)) +
  ylim(c(0, 30)) +
  scale_color_manual("", values = fixedColors[["status"]]) +
  scale_fill_manual(values = fixedColors[["status"]], guide = "none") +
  theme_classic(base_size = 16) +
  geom_hline(yintercept = -log(0.05), linetype = "dashed", color = "grey80") +
  theme(legend.position = "bottom")

ggsave("revisions/Plots/spia_sp6.png",units = "cm",width = 15,height = 15, dpi = 300)

#----------------------------------------------------------------------
#----------------------- DKO---------------------
#-------------------------------------------------------------------

# Map SYMBOL → ENTREZ
resord_DKO$ENTREZ <- mapIds(org.Mm.eg.db,
                            keys = resord_DKO$Gene.name,
                            column = "ENTREZID",
                            keytype = "SYMBOL",
                            multiVals = "first")

# Remove duplicates and NAs
resord_DKO <- resord_DKO[!duplicated(resord_DKO$ENTREZ), ]
resord_DKO <- resord_DKO[!is.na(resord_DKO$ENTREZ), ]

# Define DE gene list
res_filt_DKO <- subset(resord_DKO,
                       padj <= 0.05 & (log2FoldChange > 0.585 | log2FoldChange < -0.585))

input_DKO <- res_filt_DKO$log2FoldChange
names(input_DKO) <- res_filt_DKO$ENTREZ

all_genes <- resord_DKO$ENTREZ  # now proper Entrez IDs

# Run SPIA
res_DKO <- spia(de = input_DKO, all = all_genes, organism = "mmu",
                nB = 2000, plots = TRUE, combine = "fisher")


writexl::write_xlsx(res_DKO, "output/SPIA_DKO_raw_results.xlsx")
res_spia_DKO <- readxl::read_excel("output/SPIA_DKO_raw_results.xlsx")



res_spia_DKO$status <- "NS"
res_spia_DKO$status[res_spia_DKO$pG < 0.05 & res_spia_DKO$tA > 0] <- "Inhibited in KO"
res_spia_DKO$status[res_spia_DKO$pG < 0.05 & res_spia_DKO$tA < 0] <- "Activated in KO"
fixedColors <- list('status'=c('Activated in KO'="dodgerblue",'Inhibited in KO'="orange2",'NS'="grey"))
res_spia_DKO$custom_label <- names(pathways_to_check)[match(res_spia_DKO$ID, pathways_to_check)]
ggplot(res_spia_DKO, aes(x = tA, y = -log(pG))) +
  geom_point(aes(color = status)) +
  geom_label_repel(data = subset(res_spia_DKO, ID %in% pathways_to_check),
                   aes(label = custom_label,
                       fill = status),        # colored box background
                   color = "black",           # black text + border
                   min.segment.length = 0,
                   show.legend = F,
                   box.padding = 0.6,
                   point.padding = 0.4,
                   force = 10,
                   alpha = 0.8,              # slight transparency so fill isn't too heavy
                   size = 6) +
  xlim(c(-100, 100)) +
  ylim(c(0, 30)) +
  scale_color_manual("", values = fixedColors[["status"]]) +
  scale_fill_manual(values = fixedColors[["status"]], guide = "none") +
  theme_classic(base_size = 16) +
  geom_hline(yintercept = -log(0.05), linetype = "dashed", color = "grey80") +
  theme(legend.position = "bottom")

ggsave("revisions/Plots/spia_DKO.png",units = "cm",width = 15,height = 15, dpi = 300)

############################### GO TERMS SINGLEKO DEGS ######################


sp8_david <- read.table("david/sp8_degs_david_gos.txt", sep = "\t", header = T)
sp6_david <- read.table("david/sp6_degs_david_v2.txt", sep = "\t", header = T)

sp6_david <- sp6_david %>%
  filter(FDR < 0.05) 
sp6_david$Term <- sub(".*~", "", sp6_david$Term)

sp8_david <- sp8_david %>%
  filter(FDR < 0.05) 
sp8_david$Term <- sub(".*~", "", sp8_david$Term)

sp8_david_BP <- sp8_david %>%
  filter(Category=="GOTERM_BP_DIRECT" & PValue < 0.05) %>%
  mutate(Term = gsub(".*~", "", Term))


p3 <- ggplot(data = sp6_david[c(12,17,19,20,23), ], aes(y=reorder(Term, -log10(FDR)), 
                                                        x=Count/Pop.Hits))+
  geom_point(stat = "identity", color="black", pch= 21, aes(size=Count, fill=-log10(FDR)))+ 
  scale_y_discrete(labels = wrap_format(30))+
  
  scale_x_continuous(position="top")+
  xlim(c(0,0.5))+
  labs(x="Gene Ratio")+
  scale_fill_gradient(low="olivedrab1",high =  "olivedrab4")+
  scale_size_continuous(range=c(3,6))+
  theme_bw(base_size = 10)+
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size =14, colour="black"),
        axis.text.x = element_text(colour="black"), panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())




p4 <- ggplot(data = sp8_david_BP[c(3, 7,10, 13, 18), ], aes(y=reorder(Term, -log10(FDR)), 
                                                            x=Count/Pop.Hits))+
  geom_point(stat = "identity", color="black", pch= 21, aes(size=Count, fill=-log10(FDR)))+
  scale_y_discrete(labels = wrap_format(30))+
  
  scale_x_continuous(position="top")+
  xlim(c(0,0.50))+
  
  labs(x="Gene Ratio")+
  scale_fill_gradient(low="mediumorchid1",high =  "mediumorchid4")+
  scale_size_continuous(range=c(3,6))+
  theme_bw(base_size = 10)+
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size =14, colour="black"),
        axis.text.x = element_text(colour="black"), panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())



combined_plot <- plot_grid(
  p4,p3,
  nrow = 1, # 2x2 grid
  align = "hv", # Align horizontally and vertically
  axis = "tblr" # Align axes across all plots
)

ggsave("revisions/Plots/aligned_plots_DAVID_degs.pdf", combined_plot, width = 11, height = 3.5)

write.xlsx(sp8_deg,      "revisions/output/Source_data_3_DEGs_and_DAVID.xlsx", sheetName = "Sp8_deg")
write.xlsx(sp6_deg,      "revisions/output/Source_data_3_DEGs_and_DAVID.xlsx", sheetName = "Sp6_DEGs", append = T)
write.xlsx(deg_DKO,      "revisions/output/Source_data_3_DEGs_and_DAVID.xlsx", sheetName = "DKO_DEGs", append = T)
write.xlsx(sp8_david_BP, "revisions/output/Source_data_3_DEGs_and_DAVID.xlsx", sheetName = "Sp8_DAVID", append = T)
write.xlsx(sp6_david,    "revisions/output/Source_data_3_DEGs_and_DAVID.xlsx", sheetName = "Sp6_DAVID", append = T)
write.xlsx(res_spia_sp8, "revisions/output/Source_data_3_DEGs_and_DAVID.xlsx", sheetName = "res_spia_sp8", append = T)
write.xlsx(res_spia_sp6, "revisions/output/Source_data_3_DEGs_and_DAVID.xlsx", sheetName = "res_spia_sp6", append = T)
write.xlsx(res_spia_DKO, "revisions/output/Source_data_3_DEGs_and_DAVID.xlsx", sheetName = "res_spia_DKO", append = T)

############################### DIRECT TARGETS SP8 SP6 BY THEMSELVES ######################

library(GenomicRanges)
library(xlsx)


#genes annotated to peaks, basal plus extension
sp6_pool_bpe <- read.table("revisions/great_bpe/gene_assignment/sp6_pool_20260216-public-4.0.4-k56pj0-mm10-all-region.txt",skip = 1, sep = "\t")
df_split_sp6 <- separate(sp6_pool_bpe, col = V2, into = c("Gene_1", "Gene_2", "Gene_3"), sep = ", ")
df_numbers_sp6 <- extract(df_split_sp6, col = "Gene_1", into = c("Gene_1", "Dist_1"), regex = "(.*)\\s\\((.*)\\)")
df_numbers_sp6 <- extract(df_numbers_sp6, col = "Gene_2", into = c("Gene_2", "Dist_2"), regex = "(.*)\\s\\((.*)\\)")
df_numbers_sp6 <- extract(df_numbers_sp6, col = "Gene_3", into = c("Gene_3", "Dist_3"), regex = "(.*)\\s\\((.*)\\)")
df_numbers_1_sp6 <- data.frame("peak"=df_numbers_sp6[,1], "gene"=df_numbers_sp6[,2], "dist"=df_numbers_sp6[,3])
df_numbers_2_sp6 <- data.frame("peak"=df_numbers_sp6[,1], "gene"=df_numbers_sp6[,4], "dist"=df_numbers_sp6[,5])
df_numbers_3_sp6 <- data.frame("peak"=df_numbers_sp6[,1], "gene"=df_numbers_sp6[,6], "dist"=df_numbers_sp6[,7])
df_combined_sp6 <- rbind(df_numbers_1_sp6, df_numbers_2_sp6,df_numbers_3_sp6)
df_combined_sp6 <- na.omit(df_combined_sp6)

##direct targets
#with bpe assignation, fold change filter
sp6_dtg_bpe <- as.data.frame(intersect(sp6_deg$Gene.name, df_combined_sp6$gene))
sp6_dtg_bpe <- merge(sp6_deg, sp6_dtg_bpe, by.x = "Gene.name", by.y = colnames(sp6_dtg_bpe)[1]) ##DEG info

sum(sp6_dtg_bpe$log2FoldChange < 0) # 12 repressed
sum(sp6_dtg_bpe$log2FoldChange > 0) # 26 activated

sp6_dtg_bpe_peak <- as.data.frame(intersect(sp6_dtg_bpe$Gene.name, df_combined_sp6$gene))
sp6_dtg_bpe_peak <- merge(sp6_dtg_bpe_peak, df_combined_sp6, by.x=colnames(sp6_dtg_bpe_peak[1]), by.y="gene")
sp6_dtg_bpe_peak <- merge(sp6_tss, sp6_dtg_bpe_peak, by.x = "V4", by.y = "peak") ##PEAK INFO
colnames(sp6_dtg_bpe_peak) <- c("peak","chr","start","end","genechr","genestart","geneend","gene_tss","dist_tss","dist_tss_kb","category","dist_group","gene_bpe","dist_bpe")

#### peaks in direct targets for study
sp6_dtg_bpe_peak$dist_bpe_kb <- as.numeric(sp6_dtg_bpe_peak$dist_bpe)/1000 

#peak distance assignment
sp6_dtg_bpe_peak$Distance_Group_bpe <- cut(sp6_dtg_bpe_peak$dist_bpe_kb,
                                       breaks = c(-Inf,-500, -50, -5, 0, 5, 50, 500, Inf),
                                       labels = c("< -500","-500 to -50", "-50 to -5",
                                                  "-5 to 0", "0 to 5", "5 to 50", "50 to 500",
                                                  "> 500"))#overrides basal plus extension info with tss assignation info

sp6_dtg_bpe_peak_deg <- merge(sp6_dtg_bpe_peak, sp6_deg, by.x="gene_bpe", by.y="Gene.name")

###### peaks activated/repressed distal/proximal sp6
sp6_dtg_bpe_peak_deg$status <- ifelse(sp6_dtg_bpe_peak_deg$log2FoldChange > 0, "Down in KO", "Up in KO")


#genes annotated to peaks, basal plus extension
sp8_pool_bpe <- read.table("revisions/great_bpe/gene_assignment/sp8_pool_20260216-public-4.0.4-iqDEeB-mm10-all-region.txt",skip = 1, sep = "\t")
df_split_sp8 <- separate(sp8_pool_bpe, col = V2, into = c("Gene_1", "Gene_2", "Gene_3"), sep = ", ")
df_numbers_sp8 <- extract(df_split_sp8, col = "Gene_1", into = c("Gene_1", "Dist_1"), regex = "(.*)\\s\\((.*)\\)")
df_numbers_sp8 <- extract(df_numbers_sp8, col = "Gene_2", into = c("Gene_2", "Dist_2"), regex = "(.*)\\s\\((.*)\\)")
df_numbers_sp8 <- extract(df_numbers_sp8, col = "Gene_3", into = c("Gene_3", "Dist_3"), regex = "(.*)\\s\\((.*)\\)")
df_numbers_1_sp8 <- data.frame("peak"=df_numbers_sp8[,1], "gene"=df_numbers_sp8[,2], "dist"=df_numbers_sp8[,3])
df_numbers_2_sp8 <- data.frame("peak"=df_numbers_sp8[,1], "gene"=df_numbers_sp8[,4], "dist"=df_numbers_sp8[,5])
df_numbers_3_sp8 <- data.frame("peak"=df_numbers_sp8[,1], "gene"=df_numbers_sp8[,6], "dist"=df_numbers_sp8[,7])
df_combined_sp8 <- rbind(df_numbers_1_sp8, df_numbers_2_sp8,df_numbers_3_sp8)
df_combined_sp8 <- na.omit(df_combined_sp8)

##direct targets
#with bpe assignation, fold change filter
sp8_dtg_bpe <- as.data.frame(intersect(sp8_deg$Gene.name, df_combined_sp8$gene))
sp8_dtg_bpe <- merge(sp8_deg, sp8_dtg_bpe, by.x = "Gene.name", by.y = colnames(sp8_dtg_bpe)[1]) ##DEG info

sum(sp8_dtg_bpe$log2FoldChange < 0) # 72 repressed
sum(sp8_dtg_bpe$log2FoldChange > 0) # 160 activated

sp8_dtg_bpe_peak <- as.data.frame(intersect(sp8_dtg_bpe$Gene.name, df_combined_sp8$gene))
sp8_dtg_bpe_peak <- merge(sp8_dtg_bpe_peak, df_combined_sp8, by.x=colnames(sp8_dtg_bpe_peak[1]), by.y="gene")
sp8_dtg_bpe_peak <- merge(sp8_tss, sp8_dtg_bpe_peak, by.x = "V4", by.y = "peak") ##PEAK INFO
colnames(sp8_dtg_bpe_peak) <- c("peak","chr","start","end","genechr","genestart","geneend","gene_tss","dist_tss","dist_tss_kb","category","dist_group","gene_bpe","dist_bpe")

#### peaks in direct targets for study
sp8_dtg_bpe_peak$dist_bpe_kb <- as.numeric(sp8_dtg_bpe_peak$dist_bpe)/1000 

#peak distance assignment
sp8_dtg_bpe_peak$Distance_Group_bpe <- cut(sp8_dtg_bpe_peak$dist_bpe_kb,
                                           breaks = c(-Inf,-500, -50, -5, 0, 5, 50, 500, Inf),
                                           labels = c("< -500","-500 to -50", "-50 to -5",
                                                      "-5 to 0", "0 to 5", "5 to 50", "50 to 500",
                                                      "> 500"))#overrides basal plus extension info with tss assignation info

sp8_dtg_bpe_peak_deg <- merge(sp8_dtg_bpe_peak, sp8_deg, by.x="gene_bpe", by.y="Gene.name")

###### peaks activated/repressed distal/proximal sp8
sp8_dtg_bpe_peak_deg$status <- ifelse(sp8_dtg_bpe_peak_deg$log2FoldChange > 0, "Down in KO", "Up in KO")

#####input for meme
write.table(sp8_dtg_bpe_peak_deg[,c(3,4,5,2)],"revisions/input/sp8_dtg_for_meme.bed",  quote = F, row.names = F, col.names = F, sep = "\t")
write.table(sp6_dtg_bpe_peak_deg[,c(3,4,5,2)],"revisions/input/sp6_dtg_for_meme.bed",  quote = F, row.names = F, col.names = F, sep = "\t")

write.xlsx(sp8_dtg_bpe_peak_deg,"revisions/output/Source_data_4_direct_targets_sKO.xlsx", sheetName="Sp8_dir_tar")
write.xlsx(sp6_dtg_bpe_peak_deg,"revisions/output/Source_data_4_direct_targets_sKO.xlsx", sheetName="Sp6_dir_tar", append=TRUE)

#####peak distribution
signed_log_distancedtg8 <- sign(sp8_dtg_bpe_peak_deg$dist_tss_kb) * log10(abs(sp8_dtg_bpe_peak_deg$dist_tss_kb) + 1)
sum(sp8_dtg_bpe_peak_deg$category=="Proximal") #39
sum(sp8_dtg_bpe_peak_deg$category=="Distal")#288
ggplot(data.frame(distance = signed_log_distancedtg8), aes(x = distance)) +
  geom_density(fill = "mediumorchid1", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +     # TSS
  geom_vline(xintercept = sign(c(-2, 2)) * log10(abs(c(-2, 2)) + 1), linetype = "dashed") +
  labs(
    x = "Log10 Distance to TSS (kb)",
    y = "Density"
  ) +
  annotate("rect", xmin = -log10(2+1), xmax = log10(2+1), ymin = 0, ymax = Inf, 
           fill = "mediumorchid4", alpha = 0.3) +
  scale_x_continuous(
    breaks = sign(c(-1000, -100, -10, -2, 0, 2, 10, 100, 1000)) *
      log10(abs(c(-1000, -100, -10, -2, 0, 2, 10, 100, 1000)) + 1),
    labels = c("-1000 kb", "-100 kb", "-10 kb", "-2 kb", 
               "TSS", "2 kb", "10 kb", "100 kb", "1000 kb")
  ) +geom_rug(sides = "b", color = "mediumorchid4", alpha = 0.3, linewidth = 0.2) +
  theme_minimal(base_size = 18)+
  theme( text = element_text(color = "black"),
         axis.text = element_text(color = "black"),
         axis.title = element_text(color = "black"),
         plot.title = element_text(color = "black"),
         plot.subtitle = element_text(color = "black"),
         legend.text = element_text(color = "black"),
         legend.title = element_text(color = "black"),
         strip.text = element_text(color = "black")
         ,axis.text.x=element_text(angle = -45, vjust = -0.5),
         legend.position = "top")
ggsave("revisions/Plots/density_plots_Sp8_dtg.png",dpi = 300, height = 4, width=5)

signed_log_distancedtg6 <- sign(sp6_dtg_bpe_peak_deg$dist_tss_kb) * log10(abs(sp6_dtg_bpe_peak_deg$dist_tss_kb) + 1)
sum(sp6_dtg_bpe_peak_deg$category=="Proximal") #5
sum(sp6_dtg_bpe_peak_deg$category=="Distal") #43
ggplot(data.frame(distance = signed_log_distancedtg6), aes(x = distance)) +
  geom_density(fill = "olivedrab1", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +     # TSS
  geom_vline(xintercept = sign(c(-2, 2)) * log10(abs(c(-2, 2)) + 1), linetype = "dashed") +
  labs(
    x = "Log10 Distance to TSS (kb)",
    y = "Density"
  ) +
  annotate("rect", xmin = -log10(2+1), xmax = log10(2+1), ymin = 0, ymax = Inf, 
           fill = "olivedrab4", alpha = 0.3) +
  scale_x_continuous(
    breaks = sign(c(-1000, -100, -10, -2, 0, 2, 10, 100, 1000)) *
      log10(abs(c(-1000, -100, -10, -2, 0, 2, 10, 100, 1000)) + 1),
    labels = c("-1000 kb", "-100 kb", "-10 kb", "-2 kb", 
               "TSS", "2 kb", "10 kb", "100 kb", "1000 kb")
  ) +geom_rug(sides = "b", color = "black", alpha = 0.3, linewidth = 0.2) +
  theme_minimal(base_size = 18)+
  theme( text = element_text(color = "black"),
         axis.text = element_text(color = "black"),
         axis.title = element_text(color = "black"),
         plot.title = element_text(color = "black"),
         plot.subtitle = element_text(color = "black"),
         legend.text = element_text(color = "black"),
         legend.title = element_text(color = "black"),
         strip.text = element_text(color = "black")
         ,axis.text.x=element_text(angle = -45, vjust = -0.5),
         legend.position = "top")
ggsave("D:/ChIP/revisions/Plots/density_plots_Sp6_dtg.png",dpi = 300, height = 4, width=5)


############overlap of direct targets ##########
length(unique(df_combined_sp8$gene))
length(unique(sp8_deg[ sp8_deg$log2FoldChange>0,]$Gene.name))
length(unique(sp8_deg[ sp8_deg$log2FoldChange<0,]$Gene.name))
length(unique(df_combined_sp6$gene))
length(unique(sp6_deg[ sp6_deg$log2FoldChange>0,]$Gene.name))
length(unique(sp6_deg[ sp6_deg$log2FoldChange<0,]$Gene.name))
length(unique(deg_DKO[ deg_DKO$log2FoldChange>0,]$Gene.name))
length(unique(deg_DKO[ deg_DKO$log2FoldChange<0,]$Gene.name))
length(intersect(unique(df_combined_sp8$gene), 
                 unique(sp8_deg$Gene.name)))
length(intersect(unique(df_combined_sp6$gene), 
                 unique(sp6_deg$Gene.name)))


length(unique(sp8_dtg_bpe_peak_deg[ sp8_dtg_bpe_peak_deg$log2FoldChange>0,]$gene_bpe))+
length(unique(sp8_dtg_bpe_peak_deg[ sp8_dtg_bpe_peak_deg$log2FoldChange<0,]$gene_bpe))
length(unique(sp6_dtg_bpe_peak_deg[ sp6_dtg_bpe_peak_deg$log2FoldChange>0,]$gene_bpe))
length(unique(sp6_dtg_bpe_peak_deg[ sp6_dtg_bpe_peak_deg$log2FoldChange<0,]$gene_bpe))

list_euler_sp8 <- list("Peaks Sp8" = unique(df_combined_sp8$gene),
                       "Sp8"= unique(sp8_deg$Gene.name))
list_euler_sp6 <- list("Peaks Sp6" = unique(df_combined_sp6$gene),
                       "Sp6"= unique(sp6_deg$Gene.name))

sp8_eulerr <- eulerr::euler(combinations = list_euler_sp8)
sp6_eulerr <- eulerr::euler(combinations = list_euler_sp6)

pdf("revisions/Plots/euler_direct_targets.pdf", width = 10, height = 10)
plot(sp8_eulerr, fills = c("#e499cd", "#af99e4","mediumorchid1"), labels=list(cex = 2),quantities=list(cex = 2))
plot(sp6_eulerr, fills = c("#94df84", "#dfd084","olivedrab1"), labels=list(cex = 2),quantities=list(cex = 2))
plot(sp8_eulerr, fills = c("#e499cd", "#af99e4","mediumorchid1"), labels=list(cex = 0),quantities=list(cex = 0))
plot(sp6_eulerr, fills = c("#94df84", "#dfd084","olivedrab1"), labels=list(cex = 0),quantities=list(cex = 0))
dev.off()


intersect(intersect(unique(df_combined_sp8$gene),
          unique(sp8_deg$Gene.name)),
intersect(unique(df_combined_sp6$gene),
                 unique(sp6_deg$Gene.name)))


################################################################################################
################################################################################################
################################################################################################
#### FIG3 OVERLAPS
################################################################################################
################################################################################################

### PREPARE BPE assignment
##INPUT:SEQMINER CLUSTERS
commongb_bpe <- read.table("revisions/seqminer_redo/output/common_pool.bed")
sp8_bpe <- read.table("revisions/seqminer_redo/output/sp8_spe_pool.bed")
sp6_bpe <- read.table("revisions/seqminer_redo/output/sp6_spe_pool.bed")

commongb_bpe <- commongb_bpe %>%
  mutate(peak=paste0(V1, "_",V2,"_",V3)) %>% 
  select(V1, V2, V3, peak)
sp8_bpe <- sp8_bpe %>%
  mutate(peak=paste0(V1, "_",V2,"_",V3)) %>% 
  select(V1, V2, V3, peak)
sp6_bpe <- sp6_bpe %>%
  mutate(peak=paste0(V1, "_",V2,"_",V3)) %>% 
  select(V1, V2, V3, peak)

#GENES ASSIGNED TO PEAKS 
#INPUT:BPE ASSIGNMENT FROM GREAT
commongb_bpe <- read.table("revisions/great_bpe/gene_assignment/common_20250901-public-4.0.4-g9mdsr-mm10-all-region.txt",skip = 1, sep = "\t")
sp8_bpe <- read.table("revisions/great_bpe/gene_assignment/sp8_20250901-public-4.0.4-gRlw9f-mm10-all-region.txt",skip = 1, sep = "\t")
sp6_bpe <- read.table("revisions/great_bpe/gene_assignment/sp6_20250901-public-4.0.4-06bnHY-mm10-all-region.txt",skip = 1, sep = "\t")

df_split_com <- separate(commongb_bpe, col = V2, into = c("Gene_1", "Gene_2", "Gene_3"), sep = ", ")
df_numbers_com <- extract(df_split_com, col = "Gene_1", into = c("Gene_1", "Dist_1"), regex = "(.*)\\s\\((.*)\\)")
df_numbers_com <- extract(df_numbers_com, col = "Gene_2", into = c("Gene_2", "Dist_2"), regex = "(.*)\\s\\((.*)\\)")
df_numbers_com <- extract(df_numbers_com, col = "Gene_3", into = c("Gene_3", "Dist_3"), regex = "(.*)\\s\\((.*)\\)")
df_numbers_1_com <- data.frame("peak"=df_numbers_com[,1], "gene"=df_numbers_com[,2], "dist"=df_numbers_com[,3])
df_numbers_2_com <- data.frame("peak"=df_numbers_com[,1], "gene"=df_numbers_com[,4], "dist"=df_numbers_com[,5])
df_numbers_3_com <- data.frame("peak"=df_numbers_com[,1], "gene"=df_numbers_com[,6], "dist"=df_numbers_com[,7])
df_combined_com <- rbind(df_numbers_1_com, df_numbers_2_com,df_numbers_3_com)
df_combined_com <- na.omit(df_combined_com)

df_split_sp8 <- tidyr::separate(sp8_bpe, col = V2, into = c("Gene_1", "Gene_2", "Gene_3"), sep = ", ")
df_numbers_sp8 <- extract(df_split_sp8, col = "Gene_1", into = c("Gene_1", "Dist_1"), regex = "(.*)\\s\\((.*)\\)")
df_numbers_sp8 <- extract(df_numbers_sp8, col = "Gene_2", into = c("Gene_2", "Dist_2"), regex = "(.*)\\s\\((.*)\\)")
df_numbers_sp8 <- extract(df_numbers_sp8, col = "Gene_3", into = c("Gene_3", "Dist_3"), regex = "(.*)\\s\\((.*)\\)")
df_numbers_1_sp8 <- data.frame("peak"=df_numbers_sp8[,1], "gene"=df_numbers_sp8[,2], "dist"=df_numbers_sp8[,3])
df_numbers_2_sp8 <- data.frame("peak"=df_numbers_sp8[,1], "gene"=df_numbers_sp8[,4], "dist"=df_numbers_sp8[,5])
df_numbers_3_sp8 <- data.frame("peak"=df_numbers_sp8[,1], "gene"=df_numbers_sp8[,6], "dist"=df_numbers_sp8[,7])
df_combined_sp8 <- rbind(df_numbers_1_sp8, df_numbers_2_sp8,df_numbers_3_sp8)
df_combined_sp8 <- na.omit(df_combined_sp8)

df_split_sp6 <- tidyr::separate(sp6_bpe, col = V2, into = c("Gene_1", "Gene_2", "Gene_3"), sep = ", ")
df_numbers_sp6 <- extract(df_split_sp6, col = "Gene_1", into = c("Gene_1", "Dist_1"), regex = "(.*)\\s\\((.*)\\)")
df_numbers_sp6 <- extract(df_numbers_sp6, col = "Gene_2", into = c("Gene_2", "Dist_2"), regex = "(.*)\\s\\((.*)\\)")
df_numbers_sp6 <- extract(df_numbers_sp6, col = "Gene_3", into = c("Gene_3", "Dist_3"), regex = "(.*)\\s\\((.*)\\)")
df_numbers_1_sp6 <- data.frame("peak"=df_numbers_sp6[,1], "gene"=df_numbers_sp6[,2], "dist"=df_numbers_sp6[,3])
df_numbers_2_sp6 <- data.frame("peak"=df_numbers_sp6[,1], "gene"=df_numbers_sp6[,4], "dist"=df_numbers_sp6[,5])
df_numbers_3_sp6 <- data.frame("peak"=df_numbers_sp6[,1], "gene"=df_numbers_sp6[,6], "dist"=df_numbers_sp6[,7])
df_combined_sp6 <- rbind(df_numbers_1_sp6, df_numbers_2_sp6,df_numbers_3_sp6)
df_combined_sp6 <- na.omit(df_combined_sp6)


# RUN IN BASH
#### distribution of peaks 
# sort -k1,1 -k2,2n /home/castillaa/chip/revisions/input/sp6_bpe_for_great.bed        > /home/castillaa/chip/revisions/input/sp6_bpe_for_great_s.bed       
# sort -k1,1 -k2,2n /home/castillaa/chip/revisions/input/sp8_bpe_for_great.bed        > /home/castillaa/chip/revisions/input/sp8_bpe_for_great_s.bed       
# sort -k1,1 -k2,2n /home/castillaa/chip/revisions/input/common_with_gb_for_great.bed > /home/castillaa/chip/revisions/input/common_with_gb_for_great_s.bed
# 
# bedtools closest -D ref -t first -a /home/castillaa/chip/revisions/input/sp6_bpe_for_great_s.bed -b input.tss.sorted.bed >        /home/castillaa/chip/revisions/input/distance_to_tss_sp6_enriched_pool.bed
# bedtools closest -D ref -t first -a /home/castillaa/chip/revisions/input/sp8_bpe_for_great_s.bed -b input.tss.sorted.bed >        /home/castillaa/chip/revisions/input/distance_to_tss_sp8_enriched_pool.bed #
# bedtools closest -D ref -t first -a /home/castillaa/chip/revisions/input/common_with_gb_for_great_s.bed -b input.tss.sorted.bed > /home/castillaa/chip/revisions/input/distance_to_tss_common_with__pool.bed

# INPUT: DISTRIBUTION OF PEAKS. ANNOTATION TO NEAREST TSS
commongb_tss <-   read.table("revisions/input/distance_to_tss_common_with__pool.bed", sep = "\t")
sp8_enrich_tss <- read.table("revisions/input/distance_to_tss_sp8_enriched_pool.bed", sep = "\t")
sp6_enrich_tss <- read.table("revisions/input/distance_to_tss_sp6_enriched_pool.bed", sep = "\t")


commongb_tss$V12 <- commongb_tss$V9/1000
sp8_enrich_tss$V12 <- sp8_enrich_tss$V9/1000
sp6_enrich_tss$V12 <- sp6_enrich_tss$V9/1000


commongb_tss$Distance_Group <- cut(commongb_tss$V12,
                                   breaks = c(-Inf,-500, -50, -2, 0, 2, 50, 500, Inf),
                                   labels = c("< -500","-500 to -50", "-50 to -5",
                                              "-2 to 0", "0 to 2", "2 to 50", "50 to 500",
                                              "> 500"))
sp8_enrich_tss$Distance_Group <- cut(sp8_enrich_tss$V12,
                                     breaks = c(-Inf,-500, -50, -2, 0, 2, 50, 500, Inf),
                                     labels = c("< -500","-500 to -50", "-50 to -5",
                                                "-2 to 0", "0 to 2", "2 to 50", "50 to 500",
                                                "> 500"))
sp6_enrich_tss$Distance_Group <- cut(sp6_enrich_tss$V12,
                                     breaks = c(-Inf,-500, -50, -2, 0, 2, 50, 500, Inf),
                                     labels = c("< -500","-500 to -50", "-50 to -5",
                                                "-2 to 0", "0 to 2", "2 to 50", "50 to 500",
                                                "> 500"))

commongb_tss$distance_label <- cut(commongb_tss$V12, 
                                   breaks = c(-Inf, -2, 2, Inf),
                                   labels= c( "Distal","Proximal","Distal"))
sp8_enrich_tss$distance_label <- cut(sp8_enrich_tss$V12, 
                                     breaks = c(-Inf, -2, 2, Inf),
                                     labels= c( "Distal","Proximal","Distal"))
sp6_enrich_tss$distance_label <- cut(sp6_enrich_tss$V12, 
                                     breaks = c(-Inf, -2, 2, Inf),
                                     labels= c( "Distal","Proximal","Distal"))

df_summary_commongb_tss <- commongb_tss %>%  
  group_by(Distance_Group) %>% 
  summarise(n = n()) %>% 
  mutate(pct =n/sum(n)*100) %>% 
  mutate(label = paste0(format(round(pct,0), nsmall=0), "%"))  %>% 
  ungroup()

fill_common<-c("rosybrown1","rosybrown1","rosybrown1",
       "rosybrown","rosybrown","rosybrown1","rosybrown1",
       "rosybrown1")

signed_log_distance_com <- sign(commongb_tss$V12) * log10(abs(commongb_tss$V12) + 1)
common_continuous <- ggplot(data.frame(distance = signed_log_distance_com), aes(x = distance)) +
  geom_density(fill = "rosybrown1", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +     # TSS
  geom_vline(xintercept = sign(c(-2, 2)) * log10(abs(c(-2, 2)) + 1), linetype = "dashed") +
  labs(       x = "Distance to TSS (kb)",
       y = "Density") +
  annotate("rect", xmin = -log10(2+1), xmax = log10(2+1), ymin = 0, ymax = Inf, 
           fill = "rosybrown", alpha = 0.3) +
  scale_x_continuous(
    breaks = sign(c(-1000, -100, -10, -2, 0, 2, 10, 100, 1000)) *
      log10(abs(c(-1000, -100, -10, -2, 0, 2, 10, 100, 1000)) + 1),
    labels = c("-1000 kb", "-100 kb", "-10 kb", "-2 kb", 
               "TSS", "2 kb", "10 kb", "100 kb", "1000 kb")
  ) +geom_rug(sides = "b", color = "black", alpha = 0.3, linewidth = 0.2) +
  theme_minimal(base_size = 16)+
  theme( text = element_text(color = "black"),
         axis.text = element_text(color = "black"),
         axis.title = element_text(color = "black"),
         plot.title = element_text(color = "black"),
         plot.subtitle = element_text(color = "black"),
         legend.text = element_text(color = "black"),
         legend.title = element_text(color = "black"),
         strip.text = element_text(color = "black")
         ,axis.text.x=element_text(angle = -45, vjust = -0.5))
ggsave("revisions/Plots/continuous_common_2kb.png", common_continuous, dpi = 300, height = 4, width = 5)


signed_log_distance_s8s <- sign(sp8_enrich_tss$V12) * log10(abs(sp8_enrich_tss$V12) + 1)
sp8_spe_continuous <- ggplot(data.frame(distance = signed_log_distance_s8s), aes(x = distance)) +
  geom_density(fill = "mediumorchid1", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +     # TSS
  geom_vline(xintercept = sign(c(-2, 2)) * log10(abs(c(-2, 2)) + 1), linetype = "dashed") +
  labs(       x = "Distance to TSS (kb)",
       y = "Density") +
  annotate("rect", xmin = -log10(2+1), xmax = log10(2+1), ymin = 0, ymax = Inf, 
           fill = "mediumorchid", alpha = 0.3) +
  scale_x_continuous(
    breaks = sign(c(-1000, -100, -10, -2, 0, 2, 10, 100, 1000)) *
      log10(abs(c(-1000, -100, -10, -2, 0, 2, 10, 100, 1000)) + 1),
    labels = c("-1000 kb", "-100 kb", "-10 kb", "-2 kb", 
               "TSS", "2 kb", "10 kb", "100 kb", "1000 kb")
  ) +geom_rug(sides = "b", color = "black", alpha = 0.3, linewidth = 0.2) +
  theme_minimal(base_size = 18)+
  theme( text = element_text(color = "black"),
         axis.text = element_text(color = "black"),
         axis.title = element_text(color = "black"),
         plot.title = element_text(color = "black"),
         plot.subtitle = element_text(color = "black"),
         legend.text = element_text(color = "black"),
         legend.title = element_text(color = "black"),
         strip.text = element_text(color = "black")
         ,axis.text.x=element_text(angle = -45, vjust = -0.5))
ggsave("revisions/output/continuous_sp8_spe_2kb.png", sp8_spe_continuous, dpi = 300)

signed_log_distance_s6s <- sign(sp6_enrich_tss$V12) * log10(abs(sp6_enrich_tss$V12) + 1)
sp6_spe_continuous <- ggplot(data.frame(distance = signed_log_distance_s6s), aes(x = distance)) +
  geom_density(fill = "olivedrab1", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +     # TSS
  geom_vline(xintercept = sign(c(-2, 2)) * log10(abs(c(-2, 2)) + 1), linetype = "dashed") +
  labs(       x = "Distance to TSS (kb)",
       y = "Density") +
  annotate("rect", xmin = -log10(2+1), xmax = log10(2+1), ymin = 0, ymax = Inf, 
           fill = "olivedrab4", alpha = 0.3) +
  scale_x_continuous(
    breaks = sign(c(-1000, -100, -10, -2, 0, 2, 10, 100, 1000)) *
      log10(abs(c(-1000, -100, -10, -2, 0, 2, 10, 100, 1000)) + 1),
    labels = c("-1000 kb", "-100 kb", "-10 kb", "-2 kb", 
               "TSS", "2 kb", "10 kb", "100 kb", "1000 kb")
  ) +geom_rug(sides = "b", color = "black", alpha = 0.3, linewidth = 0.2)  +
  theme_minimal(base_size = 18)+
  theme( text = element_text(color = "black"),
         axis.text = element_text(color = "black"),
         axis.title = element_text(color = "black"),
         plot.title = element_text(color = "black"),
         plot.subtitle = element_text(color = "black"),
         legend.text = element_text(color = "black"),
         legend.title = element_text(color = "black"),
         strip.text = element_text(color = "black")
         ,axis.text.x=element_text(angle = -45, vjust = -0.5))
ggsave("revisions/output/continuous_sp6_spe_2kb.png", sp6_spe_continuous, dpi = 300)

commongb_tss_ids <- commongb_tss[,c(1,2,3,11,12)]
sp8_enrich_tss_ids <- sp8_enrich_tss[,c(1,2,3,11,12)]
sp6_enrich_tss_ids <- sp6_enrich_tss[,c(1,2,3,11,12)]

commongb_tss_ids <- commongb_tss_ids %>%
  mutate(peak=paste0(V1, "_", V2, "_", V3))
sp8_enrich_tss_ids <- sp8_enrich_tss_ids %>%
  mutate(peak=paste0(V1, "_", V2, "_", V3))
sp6_enrich_tss_ids <- sp6_enrich_tss_ids %>%
  mutate(peak=paste0(V1, "_", V2, "_", V3))

################################### motif content ###############################################

### input: FIMO RUN WITH MOTIFS OBTAINED FROM MEME IN .MEME FORMAT, DEFAULT PARAM

common_TAATTAYAGGBTWB      <- read.table("revisions/fimo_final/com_at.tsv",header = T)
common_GGSSSGGGRGGGGGCGGGGC<- read.table("revisions/fimo_final/com_gc.tsv",header = T)

sp8_TAATTAYAGGBTWB      <- read.table("revisions/fimo_final/sp8_at.tsv",header = T)
sp8_GGSSSGGGRGGGGGCGGGGC<- read.table("revisions/fimo_final/sp8_gc.tsv",header = T)

sp6_TAATTAYAGGBTWB      <- read.table("revisions/fimo_final/sp6_at.tsv",header = T)
sp6_GGSSSGGGRGGGGGCGGGGC<- read.table("revisions/fimo_final/sp6_gc.tsv",header = T)


####################################################################################################################
############################################################### with meme motives ##################################
####################################################################################################################

common_TAATTAYAGGBTWB_uniq <- common_TAATTAYAGGBTWB %>%  distinct(sequence_name, .keep_all = TRUE)
common_GGSSSGGGRGGGGGCGGGGC_uniq <- common_GGSSSGGGRGGGGGCGGGGC %>%  distinct(sequence_name, .keep_all = TRUE)

#using AT from meme and GC from meme

common_GC_AT_peaks <- intersect(common_TAATTAYAGGBTWB_uniq$sequence_name, common_GGSSSGGGRGGGGGCGGGGC_uniq$sequence_name)
common_AT_peaks <- setdiff(common_TAATTAYAGGBTWB_uniq$sequence_name, common_GGSSSGGGRGGGGGCGGGGC_uniq$sequence_name)
common_GC_peaks <- setdiff(common_GGSSSGGGRGGGGGCGGGGC_uniq$sequence_name, common_TAATTAYAGGBTWB_uniq$sequence_name)

### get peak IDs properly
common_TAATTAYAGGBTWB_uniq$sequence_name <- gsub("chr(\\w+):(\\d+)-(\\d+)\\(\\+\\)", "chr\\1_\\2_\\3", 
                                                 common_TAATTAYAGGBTWB_uniq$sequence_name)

common_GGSSSGGGRGGGGGCGGGGC_uniq$sequence_name <- gsub("chr(\\w+):(\\d+)-(\\d+)\\(\\+\\)", "chr\\1_\\2_\\3", 
                                                       common_GGSSSGGGRGGGGGCGGGGC_uniq$sequence_name)


common_GC_AT_peaks <-  gsub("chr(\\w+):(\\d+)-(\\d+)\\(\\+\\)", "chr\\1_\\2_\\3", common_GC_AT_peaks)
common_GC_AT_peaks <- as.data.frame(common_GC_AT_peaks)

common_GC_peaks<-  gsub("chr(\\w+):(\\d+)-(\\d+)\\(\\+\\)", "chr\\1_\\2_\\3", common_GC_peaks)
common_GC_peaks <- as.data.frame(common_GC_peaks)

common_AT_peaks<-  gsub("chr(\\w+):(\\d+)-(\\d+)\\(\\+\\)", "chr\\1_\\2_\\3", common_AT_peaks)
common_AT_peaks <- as.data.frame(common_AT_peaks)

# sp8_dtg_tss_peak <- sp8_dtg_tss_peak %>%
#   mutate(peak_id=paste0(chr.x, "_", start.x, "_", end.x))
# common_AT_peaks <- merge(common_AT_peaks, sp8_dtg_tss_peak, by.x=colnames(common_AT_peaks)[1], by.y="peak_id" )

#### abundance analysis
common_with_gb.bed <- read.table("revisions/seqminer_redo/output/common_pool.bed")
common_with_gb.bed <- common_with_gb.bed  %>%  
  mutate(V4=paste0(V1,"_",V2,"_",V3))
common_with_gb.bed <- common_with_gb.bed %>%  distinct(V4, .keep_all = TRUE)

common_orphans <- setdiff(unique(common_with_gb.bed$V4), 
                          unique(append(common_TAATTAYAGGBTWB_uniq$sequence_name,
                                        common_GGSSSGGGRGGGGGCGGGGC_uniq$sequence_name)))

# Create Data

datacom <- data.frame(
  group = c(
    paste0("GC_AT_peaks: ", length(common_GC_AT_peaks[, 1]), " (", round(length(common_GC_AT_peaks[, 1]) / length(common_with_gb.bed$V4) * 100), "%)"),
    paste0("AT_peaks: ", length(common_AT_peaks[, 1]), " (", round(length(common_AT_peaks[, 1]) / length(common_with_gb.bed$V4) * 100), "%)"),
    paste0("GC_peaks: ", length(common_GC_peaks[, 1]), " (", round(length(common_GC_peaks[, 1]) / length(common_with_gb.bed$V4) * 100), "%)"),
    paste0("other: ", length(common_orphans), " (", round(length(common_orphans) / length(common_with_gb.bed$V4) * 100), "%)")
  ),
  value = c(
    nrow(common_GC_AT_peaks),
    nrow(common_AT_peaks),
    nrow(common_GC_peaks),
    length(common_orphans)
  )
)

# Compute the position of labels
datacom <- datacom %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(datacom$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )



####################################################################################################################
############################################################### with meme motives ##################################
####################################################################################################################

sp8_TAATTAYAGGBTWB_uniq <- sp8_TAATTAYAGGBTWB %>%  distinct(sequence_name, .keep_all = TRUE)
sp8_GGSSSGGGRGGGGGCGGGGC_uniq <- sp8_GGSSSGGGRGGGGGCGGGGC %>%  distinct(sequence_name, .keep_all = TRUE)

#using AT from meme and GC from meme

sp8_GC_AT_peaks <- intersect(sp8_TAATTAYAGGBTWB_uniq$sequence_name, sp8_GGSSSGGGRGGGGGCGGGGC_uniq$sequence_name)
sp8_AT_peaks <- setdiff(sp8_TAATTAYAGGBTWB_uniq$sequence_name, sp8_GGSSSGGGRGGGGGCGGGGC_uniq$sequence_name)
sp8_GC_peaks <- setdiff(sp8_GGSSSGGGRGGGGGCGGGGC_uniq$sequence_name, sp8_TAATTAYAGGBTWB_uniq$sequence_name)

### get peak IDs properly
sp8_TAATTAYAGGBTWB_uniq$sequence_name <- gsub("chr(\\w+):(\\d+)-(\\d+)\\(\\+\\)", "chr\\1_\\2_\\3", 
                                              sp8_TAATTAYAGGBTWB_uniq$sequence_name)

sp8_GGSSSGGGRGGGGGCGGGGC_uniq$sequence_name <- gsub("chr(\\w+):(\\d+)-(\\d+)\\(\\+\\)", "chr\\1_\\2_\\3", 
                                                    sp8_GGSSSGGGRGGGGGCGGGGC_uniq$sequence_name)


sp8_GC_AT_peaks <-  gsub("chr(\\w+):(\\d+)-(\\d+)\\(\\+\\)", "chr\\1_\\2_\\3", sp8_GC_AT_peaks)
sp8_GC_AT_peaks <- as.data.frame(sp8_GC_AT_peaks)

sp8_GC_peaks<-  gsub("chr(\\w+):(\\d+)-(\\d+)\\(\\+\\)", "chr\\1_\\2_\\3", sp8_GC_peaks)
sp8_GC_peaks <- as.data.frame(sp8_GC_peaks)

sp8_AT_peaks<-  gsub("chr(\\w+):(\\d+)-(\\d+)\\(\\+\\)", "chr\\1_\\2_\\3", sp8_AT_peaks)
sp8_AT_peaks <- as.data.frame(sp8_AT_peaks)


#### abundance analysis
sp8_with_gb.bed <- read.table("revisions/seqminer_redo/output/sp8_spe_pool.bed")
sp8_with_gb.bed <- sp8_with_gb.bed  %>%  
  mutate(V4=paste0(V1,"_",V2,"_",V3))
sp8_with_gb.bed <- sp8_with_gb.bed %>%  distinct(V4, .keep_all = TRUE)

sp8_orphans <- setdiff(unique(sp8_with_gb.bed$V4), 
                       unique(append(sp8_TAATTAYAGGBTWB_uniq$sequence_name,
                                     sp8_GGSSSGGGRGGGGGCGGGGC_uniq$sequence_name)))

# Create Data

data8 <- data.frame(
  group = c(
    paste0("GC_AT_peaks: ", length(sp8_GC_AT_peaks[, 1]), " (", round(length(sp8_GC_AT_peaks[, 1]) / length(sp8_with_gb.bed$V4) * 100), "%)"),
    paste0("AT_peaks: ", length(sp8_AT_peaks[, 1]), " (", round(length(sp8_AT_peaks[, 1]) / length(sp8_with_gb.bed$V4) * 100), "%)"),
    paste0("GC_peaks: ", length(sp8_GC_peaks[, 1]), " (", round(length(sp8_GC_peaks[, 1]) / length(sp8_with_gb.bed$V4) * 100), "%)"),
    paste0("other: ", length(sp8_orphans), " (", round(length(sp8_orphans) / length(sp8_with_gb.bed$V4) * 100), "%)")
  ),
  value = c(
    nrow(sp8_GC_AT_peaks),
    nrow(sp8_AT_peaks),
    nrow(sp8_GC_peaks),
    length(sp8_orphans)
  )
)

# Compute the position of labels
data8 <- data8 %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(data8$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )


####################################################################################################################
############################################################### with meme motives ##################################
####################################################################################################################

sp6_TAATTAYAGGBTWB_uniq <- sp6_TAATTAYAGGBTWB %>%  distinct(sequence_name, .keep_all = TRUE)
sp6_GGSSSGGGRGGGGGCGGGGC_uniq <- sp6_GGSSSGGGRGGGGGCGGGGC %>%  distinct(sequence_name, .keep_all = TRUE)

#using AT from meme and GC from meme

sp6_GC_AT_peaks <- intersect(sp6_TAATTAYAGGBTWB_uniq$sequence_name, sp6_GGSSSGGGRGGGGGCGGGGC_uniq$sequence_name)
sp6_AT_peaks <- setdiff(sp6_TAATTAYAGGBTWB_uniq$sequence_name, sp6_GGSSSGGGRGGGGGCGGGGC_uniq$sequence_name)
sp6_GC_peaks <- setdiff(sp6_GGSSSGGGRGGGGGCGGGGC_uniq$sequence_name, sp6_TAATTAYAGGBTWB_uniq$sequence_name)

### get peak IDs properly
sp6_TAATTAYAGGBTWB_uniq$sequence_name <- gsub("chr(\\w+):(\\d+)-(\\d+)\\(\\+\\)", "chr\\1_\\2_\\3", 
                                              sp6_TAATTAYAGGBTWB_uniq$sequence_name)

sp6_GGSSSGGGRGGGGGCGGGGC_uniq$sequence_name <- gsub("chr(\\w+):(\\d+)-(\\d+)\\(\\+\\)", "chr\\1_\\2_\\3", 
                                                    sp6_GGSSSGGGRGGGGGCGGGGC_uniq$sequence_name)


sp6_GC_AT_peaks <-  gsub("chr(\\w+):(\\d+)-(\\d+)\\(\\+\\)", "chr\\1_\\2_\\3", sp6_GC_AT_peaks)
sp6_GC_AT_peaks <- as.data.frame(sp6_GC_AT_peaks)

sp6_GC_peaks<-  gsub("chr(\\w+):(\\d+)-(\\d+)\\(\\+\\)", "chr\\1_\\2_\\3", sp6_GC_peaks)
sp6_GC_peaks <- as.data.frame(sp6_GC_peaks)

sp6_AT_peaks<-  gsub("chr(\\w+):(\\d+)-(\\d+)\\(\\+\\)", "chr\\1_\\2_\\3", sp6_AT_peaks)
sp6_AT_peaks <- as.data.frame(sp6_AT_peaks)


#### abundance analysis
sp6_with_gb.bed <- read.table("revisions/seqminer_redo/output/sp6_spe_pool.bed")
sp6_with_gb.bed <- sp6_with_gb.bed  %>%  
  mutate(V4=paste0(V1,"_",V2,"_",V3))
sp6_with_gb.bed <- sp6_with_gb.bed %>%  distinct(V4, .keep_all = TRUE)

sp6_orphans <- setdiff(unique(sp6_with_gb.bed$V4), 
                       unique(append(sp6_TAATTAYAGGBTWB_uniq$sequence_name,
                                     sp6_GGSSSGGGRGGGGGCGGGGC_uniq$sequence_name)))

# Create Data

data6 <- data.frame(
  group = c(
    paste0("GC_AT_peaks: ", length(sp6_GC_AT_peaks[, 1]), " (", round(length(sp6_GC_AT_peaks[, 1]) / length(sp6_with_gb.bed$V4) * 100), "%)"),
    paste0("AT_peaks: ", length(sp6_AT_peaks[, 1]), " (", round(length(sp6_AT_peaks[, 1]) / length(sp6_with_gb.bed$V4) * 100), "%)"),
    paste0("GC_peaks: ", length(sp6_GC_peaks[, 1]), " (", round(length(sp6_GC_peaks[, 1]) / length(sp6_with_gb.bed$V4) * 100), "%)"),
    paste0("other: ", length(sp6_orphans), " (", round(length(sp6_orphans) / length(sp6_with_gb.bed$V4) * 100), "%)")
  ),
  value = c(
    nrow(sp6_GC_AT_peaks),
    nrow(sp6_AT_peaks),
    nrow(sp6_GC_peaks),
    length(sp6_orphans)
  )
)

# Compute the position of labels
data6 <- data6 %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(data6$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )


# Combine them into a single data frame
dataall <- rbind(datacom, data8, data6)
dataall$dataset <- c("Shared","Shared","Shared","Shared",
                     "Sp8","Sp8","Sp8","Sp8",
                     "Sp6","Sp6","Sp6","Sp6")
dataall$group_clean <- factor(c("Other","only GC peaks","GC+AT peaks","only AT peaks",
                                "Other","only GC peaks","GC+AT peaks","only AT peaks",
                                "Other","only GC peaks","GC+AT peaks","only AT peaks"), levels = c("only GC peaks","only AT peaks","GC+AT peaks","Other"))
cat = factor(c('No', 'Yes', 'NA'), levels = c('NA', 'Yes', 'No'))
# Change order of categories (bottom to top)
dataall$Category <- factor(dataall$group_clean, levels = c("only GC peaks", "only AT peaks", "GC+AT peaks", "Other"))  # Custom order
dataall[10,3] <- dataall[10,3] - 0.01  #fix 100% max

# n per dataset for x-axis labels
n_labels <- c(
  Shared = nrow(common_with_gb.bed),
  Sp8    = nrow(sp8_with_gb.bed),
  Sp6    = nrow(sp6_with_gb.bed)
)
dataset_levels <- c("Shared", "Sp8", "Sp6")
dataall$x_label <- factor(
  paste0(dataall$dataset, "\n(n=", n_labels[dataall$dataset], ")"),
  levels = paste0(dataset_levels, "\n(n=", n_labels[dataset_levels], ")")
)

motif_colors <- c(
  "GC+AT peaks" = "#59A14F",
  "only AT peaks"    = "#F28E2B",
  "only GC peaks"    = "#4E79A7",
  "Other"       = "grey75"
)

# Create the stacked bar plot
ggplot(dataall, aes(x = x_label, y = prop, fill = group_clean)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = motif_colors,
                    breaks = names(motif_colors),
                    guide  = guide_legend(reverse = TRUE)) +
  labs(y = "Percentage (%)", x = element_blank()) +
  geom_text(aes(label = paste0(sprintf("%1.1f", prop), "%")),
            position = position_stack(vjust = 0.5),
            size = 8.25, color = "black") +
  ylim(c(0, 100)) +
  theme_classic(base_size = 25) +
  theme(legend.title = element_blank(),
        axis.text.x  = element_text(size = 27, colour = "black"),
        axis.text.y  = element_text(size = 27, colour = "black"),
        axis.title.y = element_text(size = 27, colour = "black"))

ggsave("revisions/Plots/motif_distribution_stacked_plot.png", width = 10, height = 6, dpi = 600)


#### DIRECT TARGET INFERENCE
##INPUT: GENES ASSIGNED TO PEAKS
##INPUT: DEGS

#### gene intersect ### genes that are not common direct targets 
sp8_deg$status <- ifelse(sp8_deg$log2FoldChange > 0, "Down in KO", "Up in KO")

commongb_dtg8_tss4 <- as.data.frame(intersect(sp8_deg$Gene.name, df_combined_com$gene))
commongb_dtg8_tss4 <- as.data.frame(merge(df_combined_com, commongb_dtg8_tss4, by.x="gene", by.y=colnames(commongb_dtg8_tss4)[1]))
commongb_dtg8_tss4 <- as.data.frame(merge(commongb_dtg8_tss4, commongb_tss_ids, by="peak"))
commongb_dtg8_tss4 <- as.data.frame(merge(commongb_dtg8_tss4, sp8_deg, by.x="gene", by.y="Gene.name"))

sp6_deg$status <- ifelse(sp6_deg$log2FoldChange > 0, "Down in KO", "Up in KO")

commongb_dtg6_tss4 <- as.data.frame(intersect(sp6_deg$Gene.name, df_combined_com$gene))
commongb_dtg6_tss4 <- as.data.frame(merge(df_combined_com, commongb_dtg6_tss4, by.x="gene", by.y=colnames(commongb_dtg6_tss4)[1]))
commongb_dtg6_tss4 <- as.data.frame(merge(commongb_dtg6_tss4, commongb_tss_ids, by="peak"))
commongb_dtg6_tss4 <- as.data.frame(merge(commongb_dtg6_tss4, sp6_deg, by.x="gene", by.y="Gene.name"))

commongb_dtg8_tss4 <- as.data.frame(merge(commongb_dtg8_tss4, resord_sp6, by.x="gene", by.y="Gene.name"))
commongb_dtg6_tss4 <- as.data.frame(merge(commongb_dtg6_tss4, resord_sp8, by.x="gene", by.y="Gene.name"))


commongb_tss4_common_dtg <- as.data.frame(intersect(commongb_dtg8_tss4$gene, commongb_dtg6_tss4$gene))
# [1] "Cd36"     "Crabp1"   "Cyp26a1"  "Diras2"   "Dlx1"     "Etv3"     "Fam189a2" "Grm4"     "Irx2"     "Lclat1"   "Lgr6"     "Lrrtm1"  
# [13] "Paqr8"    "Plcxd3"   "Rasgef1b" "Snap91"   "Sost"  
commongb_tss4_common_dtg <- merge(commongb_tss4_common_dtg, df_combined_com, by.x=colnames(commongb_tss4_common_dtg)[1], by.y="gene")
commongb_tss4_common_dtg <- as.data.frame(merge(commongb_tss4_common_dtg, commongb_tss_ids, by="peak"))
commongb_tss4_common_dtg <- as.data.frame(merge(commongb_tss4_common_dtg, sp8_deg, by.x=colnames(commongb_tss4_common_dtg[2]), by.y="Gene.name"))
commongb_tss4_common_dtg <- as.data.frame(merge(commongb_tss4_common_dtg, sp6_deg, by.x=colnames(commongb_tss4_common_dtg[1]), by.y="Gene.name"))
commongb_tss4_common_dtg <- commongb_tss4_common_dtg[,c(2:6,1,7,8,11,21,15:20,24,34,28:32)]
colnames(commongb_tss4_common_dtg) <- c("peak","distance_to_genetss", "chr","start","end","gene","Distance_group","Distance_label",
                                        "L2FC_Sp8","status_Sp8","padj_Sp8","TPM_sp8_KO_rep1","TPM_sp8_KO_rep2","TPM_sp8_WT_rep1","TPM_sp8_WT_rep2",
                                        "TPM_sp8_WT_rep3", "L2FC_Sp6","status_Sp6","padj_Sp6","TPM_sp6_KO_rep1","TPM_sp6_KO_rep2","TPM_sp6_WT_rep1","TPM_sp6_WT_rep2")
write.table(unique(commongb_tss4_common_dtg$gene), "revisions/david/input/common_for_david.txt", quote = F, row.names = F, col.names = F)

ls_eul_common_binding <- list("common_binding_associated_genes"= unique(df_combined_com$gene),
                              "DEG_sp8"= unique(sp8_deg$Gene.name),
                              "DEG_sp6"= unique(sp6_deg$Gene.name))

eul_common_binding <- eulerr::euler(combinations = ls_eul_common_binding)

pdf("revisions/Plots/eul_common_binding.pdf", width = 10, height = 10)
plot(eul_common_binding, fills = c("lightcoral","mediumorchid1","olivedrab1"), labels=list(cex = 2),quantities=list(cex = 2))
plot(eul_common_binding, fills = c("lightcoral","mediumorchid1","olivedrab1"), labels=list(cex = 0),quantities=list(cex = 0))
dev.off()

commongb_redundant <- df_combined_com %>% 
  filter(!gene %in% commongb_dtg8_tss4$gene ) %>%
  filter(!gene %in% commongb_dtg6_tss4$gene )

commongb_redundant <- merge(commongb_redundant, resord_sp8, by.x = "gene", by.y = "Gene.name")
commongb_redundant <- merge(commongb_redundant, resord_sp6, by.x = "gene", by.y = "Gene.name")

ls_eul_redund <- list("common_binding_associated_genes"= unique(commongb_redundant$gene),
                              "DEG_DKO"=unique(deg_DKO$Gene.name))

eul_redun_binding <- eulerr::euler(combinations = ls_eul_redund)

pdf("revisions/Plots/eul_redun_binding.pdf", width = 10, height = 10)
plot(eul_redun_binding, fills = c("lightcoral","darkorange"), labels=list(cex = 2),quantities=list(cex = 2))
plot(eul_redun_binding, fills = c("lightcoral","darkorange"), labels=list(cex = 0),quantities=list(cex = 0))
dev.off()

commongb_redundant_degdko <- merge(commongb_redundant,deg_DKO, by.x = "gene", by.y="Gene.name")
commongb_redundant_degdko$status <- ifelse(commongb_redundant_degdko$log2FoldChange > 0, "Down in KO", "Up in KO")

commongb_redundant_degdko <- commongb_redundant_degdko[,c(2:3,1,43,31,35,37:42,6,16,10:15,19,28,23:27)]

colnames(commongb_redundant_degdko) <- c("peak" ,            "dist" ,            "gene"      ,       "status_dko"        ,   "log2FoldChange_dko"  , "padj_dko"     ,        "TPM_RNA_WT_rep1" ,   "TPM_RNA_WT_rep2"   ,
                                         "TPM_RNA_WT_rep3",    "TPM_RNA_DKO_rep1"  , "TPM_RNA_DKO_rep2"  , "TPM_RNA_DKO_rep3"  , "log2FoldChange.sp8" ,"Category.sp8"  ,     "padj.sp8"       ,    "TPM_sp8_KO_rep1" ,
                                          "TPM_sp8_KO_rep2",  "TPM_sp8_WT_rep1",  "TPM_sp8_WT_rep2" , "TPM_sp8_WT_rep3"  ,"log2FoldChange.sp6", "Category.sp6"  ,     "padj.sp6"       ,    "TPM_sp6_KO_rep1" ,
                                          "TPM_sp6_KO_rep2",  "TPM_sp6_WT_rep1",  "TPM_sp6_WT_rep2" )

write.table(unique(commongb_redundant_degdko$gene), "revisions/david/input/redundant_degdko_genes_for_david.txt", quote = F, row.names = F, col.names = F)

commongb_dtg8_tss4 <- commongb_dtg8_tss4 %>%
  filter(!gene %in% commongb_tss4_common_dtg$gene)
commongb_dtg8_tss4 <- commongb_dtg8_tss4[,c(2:6,1,7,8,11,21,15:20,24,33,28:32)]
colnames(commongb_dtg8_tss4) <- c("peak","distance_to_genetss", "chr","start","end","gene","Distance_group","Distance_label",
                                  "L2FC_Sp8","status_Sp8","padj_Sp8","TPM_sp8_KO_rep1","TPM_sp8_KO_rep2","TPM_sp8_WT_rep1","TPM_sp8_WT_rep2",
                                  "TPM_sp8_WT_rep3", "L2FC_Sp6","status_Sp6","padj_Sp6","TPM_sp6_KO_rep1","TPM_sp6_KO_rep2","TPM_sp6_WT_rep1","TPM_sp6_WT_rep2")
write.table(unique(commongb_dtg8_tss4$gene), "revisions/david/input/common_dtg8.txt", quote = F, row.names = F, col.names = F)

commongb_dtg6_tss4 <- commongb_dtg6_tss4 %>%
  filter(!gene %in% commongb_tss4_common_dtg$gene)
commongb_dtg6_tss4 <- commongb_dtg6_tss4[,c(2:6,1,7,8,11,21,15:19,24,34,28:33)]
colnames(commongb_dtg6_tss4) <- c("peak","distance_to_genetss", "chr","start","end","gene","Distance_group","Distance_label",
                                  "L2FC_Sp6","status_Sp6","padj_Sp6","TPM_sp6_KO_rep1","TPM_sp6_KO_rep2","TPM_sp6_WT_rep1","TPM_sp6_WT_rep2",
                                  "L2FC_Sp8","status_Sp8","padj_Sp8","TPM_sp8_KO_rep1","TPM_sp8_KO_rep2","TPM_sp8_WT_rep1","TPM_sp8_WT_rep2",
                                  "TPM_sp8_WT_rep3")
write.table(unique(commongb_dtg6_tss4$gene), "revisions/david/input/common_dtg6.txt", quote = F, row.names = F, col.names = F)

###############sp8 specific

sp8_rich_dtg_tss4 <- as.data.frame(intersect(sp8_deg$Gene.name, df_combined_sp8$gene))
sp8_rich_dtg_tss4 <- as.data.frame(merge(df_combined_sp8, sp8_rich_dtg_tss4, by.x="gene", by.y=colnames(sp8_rich_dtg_tss4)[1]))
sp8_rich_dtg_tss4 <- as.data.frame(merge(sp8_rich_dtg_tss4, sp8_enrich_tss_ids, by="peak"))
sp8_rich_dtg_tss4 <- as.data.frame(merge(sp8_rich_dtg_tss4, sp8_deg, by.x="gene", by.y="Gene.name"))
sp8_rich_dtg_tss4 <- sp8_rich_dtg_tss4[,c(2:6,1,7,8,11,21,15:20)]
colnames(sp8_rich_dtg_tss4) <- c("peak","distance_to_genetss", "chr","start","end","gene","Distance_group","Distance_label",
                                 "L2FC_Sp8","status_Sp8","padj_Sp8","TPM_sp8_KO_rep1","TPM_sp8_KO_rep3","TPM_sp8_WT_rep1","TPM_sp8_WT_rep2",
                                 "TPM_sp8_WT_rep3")
write.table(unique(sp8_rich_dtg_tss4$gene), "revisions/david/input/spe8_dtg8.txt", quote = F, row.names = F, col.names = F)

eul_sp8_rich_dtg_tss4 <- list("sp8_rich_binding_associated_genes"= unique(df_combined_sp8$gene),
                              "DEG_sp8"= unique(sp8_deg$Gene.name))
eul_sp8_rich_dtg_tss4 <- eulerr::euler(combinations = eul_sp8_rich_dtg_tss4)
pdf("revisions/Plots/eul_sp8_rich_binding.pdf", width = 10, height = 10)
plot(eul_sp8_rich_dtg_tss4, fills = c("mediumpurple","mediumorchid1"), labels=list(cex = 2),quantities=list(cex = 2))
plot(eul_sp8_rich_dtg_tss4, fills = c("mediumpurple","mediumorchid1"), labels=list(cex = 0),quantities=list(cex = 0))
dev.off()

sp6_rich_dtg_tss4 <- as.data.frame(intersect(sp6_deg$Gene.name, df_combined_sp6$gene))
sp6_rich_dtg_tss4 <- as.data.frame(merge(df_combined_sp6, sp6_rich_dtg_tss4, by.x="gene", by.y=colnames(sp6_rich_dtg_tss4)[1]))
sp6_rich_dtg_tss4 <- as.data.frame(merge(sp6_rich_dtg_tss4, sp6_enrich_tss_ids, by="peak"))
sp6_rich_dtg_tss4 <- as.data.frame(merge(sp6_rich_dtg_tss4, sp6_deg, by.x="gene", by.y="Gene.name"))
sp6_rich_dtg_tss4 <- sp6_rich_dtg_tss4[,c(2:6,1,7,8,11,21,15:19)]
colnames(sp6_rich_dtg_tss4) <- c("peak","distance_to_genetss", "chr","start","end","gene","Distance_group","Distance_label",
                                 "L2FC_Sp6","status_Sp6","padj_Sp6","TPM_sp6_KO_rep1","TPM_sp6_KO_rep2","TPM_sp6_WT_rep1","TPM_sp6_WT_rep2")
write.table(unique(sp6_rich_dtg_tss4$gene), "revisions/david/input/spe6_dtg6.txt", quote = F, row.names = F, col.names = F)

eul_sp6_rich_dtg_tss4 <- list("sp6_rich_binding_associated_genes"= unique(df_combined_sp6$gene),
                              "DEG_sp6"= unique(sp6_deg$Gene.name))
eul_sp6_rich_dtg_tss4 <- eulerr::euler(combinations = eul_sp6_rich_dtg_tss4)
pdf("revisions/Plots/eul_sp6_rich_binding.pdf", width = 10, height = 10)
plot(eul_sp6_rich_dtg_tss4, fills = c("#d5fec2","olivedrab1"), labels=list(cex = 2),quantities=list(cex = 2))
plot(eul_sp6_rich_dtg_tss4, fills = c("#d5fec2","olivedrab1"), labels=list(cex = 0),quantities=list(cex = 0))
dev.off()


intersect(sp6_rich_dtg_tss4$gene, commongb_dtg6_tss4$gene)
# "Nfe2l3"
intersect(sp8_rich_dtg_tss4$gene, commongb_dtg8_tss4$gene)
# [1] "Ankrd44"  "Aox3"     "Arrdc4"   "Barx2"    "Creb3l1"  "Dact1"    "Ddhd1"    "Dlx2"     "Dlx4"     "Dnajc1"   "Dok5"     "Epha4"    "Haao"     "Irx4"    
# [15] "Irx5"     "Klf14"    "Lnx1"     "Mycn"     "Nfil3"    "Nr3c1"    "Nrxn3"    "Palld"    "Pcsk6"    "Pdzd2"    "Pla2r1"   "Pou4f1"   "Sall3"    "Sdk2"    
# [29] "Skap2"    "Slc39a11" "Sorbs2"   "Sox9"     "Spry1"    "Tfap2a"   "Tiam2"    "Ttc28"    "Txnrd1"   "Ust"      "Wnt5a"    "Zfp322a"  "Zfp36l2" 
intersect(sp8_rich_dtg_tss4$gene, sp6_rich_dtg_tss4$gene) 


# ── Build motif annotation per peak ──────────────────────────────────────────
build_motif_lookup <- function(at_uniq, gc_uniq, all_bed) {
  
  # Standardize peak IDs in bed file
  all_bed <- all_bed %>%
    mutate(peak_id = paste0(V1, "_", V2, "_", V3)) %>%
    distinct(peak_id)
  
  # Classify each peak
  has_AT <- at_uniq$sequence_name
  has_GC <- gc_uniq$sequence_name
  
  all_bed %>%
    mutate(motif_content = case_when(
      peak_id %in% has_AT & peak_id %in% has_GC ~ "GC+AT",
      peak_id %in% has_AT                        ~ "AT_only",
      peak_id %in% has_GC                        ~ "GC_only",
      TRUE                                       ~ "No_motif"
    ))
}

# Build one lookup per binding category
motif_common <- build_motif_lookup(common_TAATTAYAGGBTWB_uniq, 
                                   common_GGSSSGGGRGGGGGCGGGGC_uniq,
                                   common_with_gb.bed)

motif_sp8 <- build_motif_lookup(sp8_TAATTAYAGGBTWB_uniq,
                                sp8_GGSSSGGGRGGGGGCGGGGC_uniq,
                                sp8_with_gb.bed)

motif_sp6 <- build_motif_lookup(sp6_TAATTAYAGGBTWB_uniq,
                                sp6_GGSSSGGGRGGGGGCGGGGC_uniq,
                                sp6_with_gb.bed)

# Combine into one master motif lookup
motif_lookup <- bind_rows(motif_common, motif_sp8, motif_sp6) %>%
  distinct(peak_id, .keep_all = TRUE) %>%
  select(peak_id, motif_content)

# ── Add peak_id column and merge motif info into each source dataframe ────────
add_motif_to_df <- function(df) {
  # Build peak_id from coordinates — adjust column names if different
  df %>%
    mutate(peak_id = gsub(":", "_", gsub("-", "_", peak))) %>%  # adapt if your peak col uses chr:start-end format
    left_join(motif_lookup, by = "peak_id") %>%
    mutate(motif_content = replace_na(motif_content, "No_motif"))
}

commongb_tss4_common_dtg  <- add_motif_to_df(commongb_tss4_common_dtg)
commongb_redundant_degdko <- add_motif_to_df(commongb_redundant_degdko)
commongb_dtg8_tss4        <- add_motif_to_df(commongb_dtg8_tss4)
commongb_dtg6_tss4        <- add_motif_to_df(commongb_dtg6_tss4)
sp8_rich_dtg_tss4         <- add_motif_to_df(sp8_rich_dtg_tss4)
sp6_rich_dtg_tss4         <- add_motif_to_df(sp6_rich_dtg_tss4)

######################## GO TERMs PLOTS BY DAVID TOOL #########################

#INPUT: DAVID WEB ANALYSIS WITH EACH CATEGORIES' GENES
com_david <- read.table("revisions/david/output/common_both.txt", sep = "\t", header = T)
cmr_david <- read.table("revisions/david/output/DAVIDChartReport_User List_2026-06-19_redunDKO.csv", sep = ",", header = T)
cm8_david <- read.table("revisions/david/output/common_sp8.txt", sep = "\t", header = T)
cm6_david <- read.table("revisions/david/output/common_sp6.txt", sep = "\t", header = T)
sp8_david <- read.table("revisions/david/output/sp8_spe.txt", sep = "\t", header = T)
# sp6_david <- read.table("revisions/david/output/spe_sp6.txt", sep = "\t", header = T)

com_david_BP <- com_david %>%
  filter( PValue < 0.05) %>%
  mutate(Term = gsub(".*~", "", Term))

cmr_david_BP <- cmr_david %>%
  filter(P.Value < 0.05) %>%
  mutate(Term = gsub(".*~", "", Term))

cm8_david_BP <- cm8_david %>%
  filter(PValue < 0.05) %>%
  mutate(Term = gsub(".*~", "", Term))

cm6_david_BP <- cm6_david %>%
  filter(PValue < 0.05) %>%
  mutate(Term = gsub(".*~", "", Term))

# sp6_david_BP <- sp6_david %>%
  # filter( PValue < 0.05) %>%
  # mutate(Term = gsub(".*~", "", Term))

sp8_david_BP <- sp8_david %>%
  filter( PValue < 0.05) %>%
  mutate(Term = gsub(".*~", "", Term))

###### TOP XX selected

# Set a reusable theme
my_theme <- theme_bw(base_size = 10) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_text(colour = "black"),
    axis.text.x  = element_text(colour = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

# Plot 1
p1 <- ggplot(com_david_BP[c(1,3,7,6), ],
             aes(y = reorder(Term, -log10(PValue)),
                 x = Count/Pop.Hits,
                 size = Count,
                 fill = -log10(PValue))) +
  geom_point(stat = "identity", color = "black", pch = 21) +
  scale_y_discrete(labels = wrap_format(30)) +
  scale_x_continuous(position = "top",
                     limits = c(0, 0.2),
                     breaks = c(0.0, 0.1, 0.2),
                     expand = expansion(mult = c(0, 0.06))) +
  labs(x = "Gene Ratio") +
  scale_fill_continuous(low = "lightcoral", high = "lightsalmon4") +
  scale_size_continuous(range = c(3, 6)) +
  theme_bw(base_size = 10) +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(colour = "black"), panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

# Plot 2
p2 <- ggplot(cm8_david_BP[c(4, 15, 25, 26, 27), ],
             aes(y = reorder(Term, -log10(PValue)),
                 x = Count/Pop.Hits,
                 size = Count,
                 fill = -log10(PValue))) +
  geom_point(stat = "identity", color = "black", pch = 21) +
  scale_y_discrete(labels = wrap_format(30)) +
  scale_x_continuous(position = "top",
                     limits = c(0, 0.5),
                     breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5),
                     expand = expansion(mult = c(0, 0.06))) +
  labs(x = "Gene Ratio") +
  scale_fill_continuous(low = "lightcoral", high = "lightsalmon4") +
  scale_size_continuous(range = c(3, 6)) +
  theme_bw(base_size = 10) +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(colour = "black"), panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

# Plot 3
p3 <- ggplot(cm6_david_BP[c(3,6,8,12), ],
             aes(y = reorder(Term, -log10(PValue)),
                 x = Count/Pop.Hits,
                 size = Count,
                 fill = -log10(PValue))) +
  geom_point(stat = "identity", color = "black", pch = 21) +
  scale_y_discrete(labels = wrap_format(30)) +
  scale_x_continuous(position = "top",
                     limits = c(0, 0.2),
                     breaks = c(0.0, 0.1, 0.2),
                     expand = expansion(mult = c(0, 0.06))) +
  labs(x = "Gene Ratio") +
  scale_fill_continuous(low = "lightcoral", high = "lightsalmon4") +
  scale_size_continuous(range = c(3, 6)) +
  theme_bw(base_size = 10) +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(colour = "black"), panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

# Plot 4
p4 <- ggplot(sp8_david_BP[c(4, 19, 24, 28, 87), ],
             aes(y = reorder(Term, -log10(PValue)),
                 x = Count/Pop.Hits,
                 size = Count,
                 fill = -log10(PValue))) +
  geom_point(stat = "identity", color = "black", pch = 21) +
  scale_y_discrete(labels = wrap_format(30)) +
  scale_x_continuous(position = "top",
                     limits = c(0, 0.3),
                     breaks = c(0.0, 0.1, 0.2, 0.3),
                     expand = expansion(mult = c(0, 0.06))) +
  labs(x = "Gene Ratio") +
  scale_fill_continuous(low = "mediumorchid1", high = "mediumorchid4") +
  scale_size_continuous(range = c(3, 6)) +
  theme_bw(base_size = 10) +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(colour = "black"), panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

# Plot 5
p5 <- ggplot(cmr_david[c(12,14,15,17,28), ],
             aes(y = reorder(Term, -log10(P.Value)),
                 x = Count/Pop.Hits,
                 size = Count,
                 fill = -log10(P.Value))) +
  geom_point(stat = "identity", color = "black", pch = 21) +
  scale_y_discrete(labels = wrap_format(30)) +
  scale_x_continuous(position = "top",
                     limits = c(0, 0.3),
                     breaks = c(0.0, 0.1, 0.2, 0.3),
                     expand = expansion(mult = c(0, 0.06))) +
  labs(x = "Gene Ratio") +
  scale_fill_continuous(low = "darkorange", high = "darkorange4") +
  scale_size_continuous(range = c(3, 6)) +
  theme_bw(base_size = 10) +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(colour = "black"), panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

combined_plot <- (p2 | p4) / (p1 | p3) / (p5 | p5)

ggsave("revisions/Plots/DAVID_Fig3_aligned_plots.jpg", combined_plot, dpi = 300, width = 10, height = 9.75)

xlsx::write.xlsx(com_david_BP, "output/Source_data_6_Targets_DAVID.xlsx", sheetName = "Shared_Sp8_Sp6_Dependent_DAVID")
xlsx::write.xlsx(cmr_david_BP, "output/Source_data_6_Targets_DAVID.xlsx", sheetName = "Redundant_DAVID", append = T)
xlsx::write.xlsx(cm8_david_BP, "output/Source_data_6_Targets_DAVID.xlsx", sheetName = "Shared_Sp8_Dependent_DAVID", append = T)
xlsx::write.xlsx(cm6_david_BP, "output/Source_data_6_Targets_DAVID.xlsx", sheetName = "Shared_Sp6_Dependent_DAVID", append = T)
xlsx::write.xlsx(sp8_david_BP, "output/Source_data_6_Targets_DAVID.xlsx", sheetName = "Sp6-specific_DAVID", append = T)

#################################################
###### intersects with prev lists fig s4e #########
#################################################


intersect(unique(deg_DKO$Gene.name), unique(sp6_deg$Gene.name))
intersect(unique(deg_DKO$Gene.name), unique(sp8_deg$Gene.name))
intersect(unique(deg_DKO$Gene.name), unique(commongb_tss4_common_dtg$gene)); length(unique(commongb_tss4_common_dtg$gene))
length(intersect(unique(deg_DKO$Gene.name), unique(commongb_tss4_common_dtg$gene)))/length(unique(commongb_tss4_common_dtg$gene)) *100
intersect(unique(deg_DKO$Gene.name), unique(commongb_dtg8_tss4$gene)) ; length(unique(commongb_dtg8_tss4$gene))
length(intersect(unique(deg_DKO$Gene.name), unique(commongb_dtg8_tss4$gene)))/length(unique(commongb_dtg8_tss4$gene)) *100
intersect(unique(deg_DKO$Gene.name), unique(commongb_dtg6_tss4$gene)); length(unique(commongb_dtg6_tss4$gene))
length(intersect(unique(deg_DKO$Gene.name), unique(commongb_dtg6_tss4$gene)))/length(unique(commongb_dtg6_tss4$gene)) *100
intersect(unique(deg_DKO$Gene.name), unique(sp8_rich_dtg_tss4$gene)); length(unique(sp8_rich_dtg_tss4$gene))
length(intersect(unique(deg_DKO$Gene.name), unique(sp8_rich_dtg_tss4$gene)))/length(unique(sp8_rich_dtg_tss4$gene)) *100
intersect(unique(deg_DKO$Gene.name), unique(sp6_rich_dtg_tss4$gene)); length(unique(sp6_rich_dtg_tss4$gene))
length(intersect(unique(deg_DKO$Gene.name), unique(sp6_rich_dtg_tss4$gene)))/length(unique(sp6_rich_dtg_tss4$gene)) *100

gene_sets <- list(
  "Shared Sp8/Sp6-dependent" = commongb_tss4_common_dtg,
  "Shared Sp8-dependent"     = commongb_dtg8_tss4,
  "Shared Sp6-dependent"     = commongb_dtg6_tss4,
  "Sp8-specific"             = sp8_rich_dtg_tss4,
  "Sp6-specific"             = sp6_rich_dtg_tss4
)

overlap_table <- do.call(rbind, lapply(names(gene_sets), function(name) {
  gs    <- unique(gene_sets[[name]]$gene)
  dko   <- unique(deg_DKO$Gene.name)
  n_int <- length(intersect(dko, gs))
  n_tot <- length(gs)
  data.frame(
    "Gene category"  = name,
    "Overlap / Total" = paste0(n_int, "/", n_tot),
    "Overlap (%)"    = paste0(round(n_int / n_tot * 100, 1), "%"),
    check.names = FALSE
  )
}))

ggtexttable(overlap_table, rows = NULL,
            theme = ttheme(
              colnames.style = colnames_style(fill = "#2c2c2c", color = "white",
                                              face = "bold", size = 10),
              tbody.style    = tbody_style(fill = c("#f0f0f0", "white"),
                                           color = "black", size = 10)
            ))

ggsave("revisions/Plots/table_overlaps.jpeg", dpi = 300, width = 5)


##################################################
##### sc RNAseq gene annotations -------------------------
##################################################

library(Seurat)
data10<-readRDS(file = "input/E105_Allou_Kelly.RDS")
# load("Rdata_sc105.RData")
FeaturePlot(data10, "Fgf8", label = T, pt.size = 1.5, order = T, cols = c("grey","red"))
FeaturePlot(data10, "Epcam", label = T, pt.size = 1.5, order = T, cols = c("grey","red"))
FeaturePlot(data10, "Pecam", label = T, pt.size = 1.5, order = T, cols = c("grey","red"))
FeaturePlot(data10, "Prrx1", label = T, pt.size = 1.5, order = T, cols = c("grey","red"))
FeaturePlot(data10, "Sox10", label = T, pt.size = 1.5, order = T, cols = c("grey","red"))
FeaturePlot(data10, "Pax3", label = T, pt.size = 1.5, order = T, cols = c("grey","red"))
FeaturePlot(data10, "Hba-a2", label = T, pt.size = 1.5, order = T, cols = c("grey","red"))

# Idents(data10) <- "RNA_snn_res.0.3"
data10 <- RenameIdents(data10,
                       "0" = "Mesoderm_1",
                       "1" = "Mesoderm_2",
                       "2" = "Mesoderm_3",
                       "3"="Mesoderm_4",
                       "4"="Mesoderm_5",
                       "5"="Muscular",
                       "6"="nonAER_ECT",
                       "7"="Endothelial",
                       "8"="Blood",
                       "9"="AER",
                       "10"="NC_lineage")


FeaturePlot(data10, "Cyp26b1", label = T, pt.size = 1, order = T, cols = c("yellow","purple"), repel = T)

############################ fully annotated pathways ###################################


celltype_map <- c("0"  = "Mesoderm_1",
                  "1"  = "Mesoderm_2",
                  "2"  = "Mesoderm_3",
                  "3"  = "Mesoderm_4",
                  "4"  = "Mesoderm_5",
                  "5"  = "Muscular",
                  "6"  = "nonAER_ECT",
                  "7"  = "Endothelial",
                  "8"  = "Blood",
                  "9"  = "AER",
                  "10" = "Schwann.Melanoblast")

#  New x-axis labels: celltype + cluster number
new_xlabels <- paste0(celltype_map, "\n(", names(celltype_map), ")")
names(new_xlabels) <- names(celltype_map)
#  Category rename map
category_map <- c(
  "Enriched_AER"          = "AER_enriched",
  "Enriched_ECT"          = "nonAER_ECT_enriched",
  "Enriched_both_ECT_AER" = "panECT_enriched",
  "AER_ECT_broad"         = "panECT_broad",
  "AER_broad"             = "AER_broad",
  "ECT_broad"             = "nonAER_ECT_broad",
  "Absent_AER_ECT"        = "ECT_absent"
)

#  Color palette
# Organized by biological hierarchy:
# AER (purples) ??? pan-ectoderm (greens) ??? non-AER ECT (teals) ??? broad/absent (grays)
ect_colors <- c(
  "AER_enriched"       = "#7B2D8B",   # deep purple  - AER specific, high confidence
  "AER_broad"          = "#C084D4",   # light purple - AER present, not exclusive
  "panECT_enriched"    = "#1B7A4A",   # deep green   - enriched across all ectoderm
  "panECT_broad"       = "#74C496",   # light green  - broad pan-ectodermal
  "nonAER_ECT_enriched"= "#0E6B7A",   # deep teal    - non-AER ectoderm specific
  "nonAER_ECT_broad"   = "#6DC0CC",   # light teal   - non-AER ectoderm broad
  "ECT_absent"         = "#D3D3D3",   # light gray   - not ectodermal
  "Not_in_scRNA"       = "#2C2C2C"    # near black   - not detected in scRNA
)

# Directory setup
dir.create("revisions/FeaturePlots/Common",        recursive = TRUE, showWarnings = FALSE)
dir.create("revisions/FeaturePlots/Common_sp8",    recursive = TRUE, showWarnings = FALSE)
dir.create("revisions/FeaturePlots/Common_sp6",    recursive = TRUE, showWarnings = FALSE)
dir.create("revisions/FeaturePlots/spe6_dtg6",     recursive = TRUE, showWarnings = FALSE)
dir.create("revisions/FeaturePlots/spe8_dtg8",     recursive = TRUE, showWarnings = FALSE)
dir.create("revisions/FeaturePlots/redun",         recursive = TRUE, showWarnings = FALSE)
dir.create("revisions/FeaturePlots/notbound_DEGs", recursive = TRUE, showWarnings = FALSE)

geneset_config <- list(
  "Common"       = list(data = data10, output_dir = "revisions/FeaturePlots/Common"),
  "Common_sp8"   = list(data = data10, output_dir = "revisions/FeaturePlots/Common_sp8"),
  "Common_sp6"   = list(data = data10, output_dir = "revisions/FeaturePlots/Common_sp6"),
  "spe6_dtg6"    = list(data = data10, output_dir = "revisions/FeaturePlots/spe6_dtg6"),
  "spe8_dtg8"    = list(data = data10, output_dir = "revisions/FeaturePlots/spe8_dtg8"),
  "redun"        = list(data = data10, output_dir = "revisions/FeaturePlots/redun"),
  "notbound_DEGs"= list(data = data10, output_dir = "revisions/FeaturePlots/notbound_DEGs")
)

# Gene sets
# deg_DKO <- read.xlsx("output/deg_DKO.xlsx", sheetIndex = 1)
GS_dtr_sKO <- list(
  "Sp8_dtr" = unique(sp8_deg$Gene.name),
  "Sp6_dtr" = unique(sp6_deg$Gene.name),
  "DKO_dtr" = unique(deg_DKO$Gene.name)
)
GS_dtr_sKO <- unlist(GS_dtr_sKO)

common     <- as.vector(unique(commongb_tss4_common_dtg$gene))
common_sp8 <- as.vector(unique(setdiff(commongb_dtg8_tss4$gene, commongb_tss4_common_dtg$gene)))
common_sp6 <- as.vector(unique(setdiff(commongb_dtg6_tss4$gene, commongb_tss4_common_dtg$gene)))
spe_sp8    <- as.vector(unique(sp8_rich_dtg_tss4$gene))
spe_sp6    <- as.vector(unique(sp6_rich_dtg_tss4$gene))
redun      <- as.vector(unique(intersect(commongb_redundant$gene, deg_DKO$Gene.name)))
genes_to_remove <- c(common, common_sp8, common_sp6, spe_sp6, spe_sp8, redun)
notbound_DEGs   <- as.vector(GS_dtr_sKO[!(GS_dtr_sKO %in% genes_to_remove)])

GS <- list(
  "Common"        = common,
  "Common_sp8"    = common_sp8,
  "Common_sp6"    = common_sp6,
  "spe6_dtg6"     = spe_sp6,
  "spe8_dtg8"     = spe_sp8,
  "redun"         = redun,
  "notbound_DEGs" = notbound_DEGs
)

# Main loop
exp_threshold <- 0.2
all_results   <- list()


####### very long run, caution #######
for (gs_name in names(GS)) {
  
  gene_list  <- GS[[gs_name]]
  seurat_obj <- geneset_config[[gs_name]]$data
  output_dir <- geneset_config[[gs_name]]$output_dir
  
  message("\n========== Processing gene set: ", gs_name, " ==========")
  
  results_list <- list()
  
  for (gene in gene_list) {
    
    if (gene %in% rownames(seurat_obj)) {
      
      # 1. Fetch data
      plot_data <- FetchData(seurat_obj, vars = c(gene, "ident"))
      colnames(plot_data) <- c("Expression", "Cluster")
      
      # 2. Per-cluster stats
      label_data <- plot_data %>%
        group_by(Cluster) %>%
        summarize(
          pct       = mean(Expression > exp_threshold) * 100,
          n_cells   = n(),
          mean_expr = mean(Expression),
          y_pos     = max(Expression) + (max(plot_data$Expression) * 0.05),
          .groups   = "drop"
        ) %>%
        mutate(label = paste0(round(pct, 1), "%"))
      
      # 3. Extract C6 and C9
      pct_C6 <- label_data %>% filter(Cluster == "nonAER_ECT") %>% pull(pct)
      pct_C9 <- label_data %>% filter(Cluster == "AER") %>% pull(pct)
      pct_C6 <- ifelse(length(pct_C6) == 0, 0, pct_C6)
      pct_C9 <- ifelse(length(pct_C9) == 0, 0, pct_C9)
      
      # And update the others filter accordingly
      pct_others <- label_data %>% filter(!Cluster %in% c("nonAER_ECT", "AER")) %>% pull(pct)
      
      # 4. Others
      max_pct_others  <- ifelse(length(pct_others) == 0, 0, max(pct_others))
      mean_pct_others <- ifelse(length(pct_others) == 0, 0, mean(pct_others))
      
      # 5. Ratios
      pseudocount <- 0.1
      ratio_C6 <- (pct_C6 + pseudocount) / (max_pct_others + pseudocount)
      ratio_C9 <- (pct_C9 + pseudocount) / (max_pct_others + pseudocount)
      
      # 6. Classify with new category names
      absent_threshold <- 2
      ratio_threshold  <- 1.5
      
      classification <- case_when(
        pct_C6 < absent_threshold & pct_C9 < absent_threshold                       ~ "ECT_absent",
        ratio_C6 > ratio_threshold & ratio_C9 > ratio_threshold &
          (pct_C6 >= absent_threshold | pct_C9 >= absent_threshold)                 ~ "panECT_enriched",
        ratio_C6 > ratio_threshold & pct_C6 >= absent_threshold                     ~ "nonAER_ECT_enriched",
        ratio_C9 > ratio_threshold & pct_C9 >= absent_threshold                     ~ "AER_enriched",
        pct_C9 >= absent_threshold & pct_C6 >= absent_threshold                     ~ "panECT_broad",
        pct_C9 >= absent_threshold & pct_C6 <  absent_threshold                     ~ "AER_broad",
        pct_C6 >= absent_threshold & pct_C9 <  absent_threshold                     ~ "nonAER_ECT_broad",
        TRUE                                                                          ~ "ECT_absent"
      )
      
      # 7. Store metrics
      results_list[[gene]] <- data.frame(
        gene_set        = gs_name,
        gene            = gene,
        pct_C6          = round(pct_C6, 2),
        pct_C9          = round(pct_C9, 2),
        max_pct_others  = round(max_pct_others, 2),
        mean_pct_others = round(mean_pct_others, 2),
        ratio_C6        = round(ratio_C6, 3),
        ratio_C9        = round(ratio_C9, 3),
        classification  = classification,
        stringsAsFactors = FALSE
      )
      
      # 8. Violin plot with cell type labels and new colors
      q <- VlnPlot(seurat_obj, features = gene, pt.size = 0.75) +
        geom_text(data = label_data,
                  aes(x = Cluster, y = y_pos, label = label),
                  inherit.aes = FALSE, size = 3.5, fontface = "bold") +
        scale_x_discrete(labels = new_xlabels) +
        scale_fill_manual(values = ect_colors) +
        ggtitle(paste0("[", gs_name, "] ", gene,
                       " (% > ", exp_threshold, ") | ", classification)) +
        ylim(NA, max(label_data$y_pos) * 1.1) +
        theme(axis.text.x = element_text(size = 8))
      
      file_name <- file.path(output_dir, paste0(gene, "_vln_labeled.png"))
      ggsave(filename = file_name, plot = q, width = 7, height = 6)
      
    } else {
      message("Gene ", gene, " not found in ", gs_name)
    }
  }
  
  # 9. Combine results
  if (length(results_list) > 0) {
    gs_df <- bind_rows(results_list) %>%
      arrange(desc(pmax(ratio_C6, ratio_C9)))
    all_results[[gs_name]] <- gs_df
    
    # 10. Scatter plot with new category names and colors
    p <- ggplot(gs_df, aes(x = ratio_C6, y = ratio_C9, color = classification)) +
      geom_point(size = 3, alpha = 0.8) +
      geom_vline(xintercept = ratio_threshold, linetype = "dashed", color = "gray50") +
      geom_hline(yintercept = ratio_threshold, linetype = "dashed", color = "gray50") +
      geom_text(aes(label = gene), size = 2.5, vjust = -0.8, check_overlap = TRUE) +
      scale_color_manual(values = ect_colors) +
      labs(title = paste0("Enrichment ratios: nonAER_ECT vs AER - ", gs_name),
           x     = "Ratio nonAER_ECT / max(non-ECT clusters)",
           y     = "Ratio AER / max(non-ECT clusters)",
           color = "scRNA class") +
      theme_bw() +
      theme(legend.position = "right",axis.text.x = element_text(size = 11, angle = 45, hjust = 1))
    
    ggsave(filename = file.path(output_dir, paste0("0_", gs_name, "_ratio_scatter.png")),
           plot = p, width = 8, height = 6)
  }
}

# Single master dataframe with all gene sets
master_df <- bind_rows(all_results)

# Optional: save to CSV
write.csv(master_df, file = "gene_classification_scRNAseq.csv", row.names = FALSE)

master_df <- read.csv("gene_classification_scRNAseq.csv", header = T)

scrna_lookup <- master_df %>%
  select(gene, classification) %>%
  distinct(gene, .keep_all = TRUE)

#add this info to source data
commongb_tss4_common_dtg  <- merge(commongb_tss4_common_dtg  , scrna_lookup, by = "gene")
commongb_redundant_degdko <- merge(commongb_redundant_degdko , scrna_lookup, by = "gene")
commongb_dtg8_tss4        <- merge(commongb_dtg8_tss4        , scrna_lookup, by = "gene")
commongb_dtg6_tss4        <- merge(commongb_dtg6_tss4        , scrna_lookup, by = "gene")
sp8_rich_dtg_tss4         <- merge(sp8_rich_dtg_tss4         , scrna_lookup, by = "gene")
sp6_rich_dtg_tss4         <- merge(sp6_rich_dtg_tss4         , scrna_lookup, by = "gene")


table(master_df$gene_set,master_df$classification)

#  Raw counts table
count_table <- master_df %>%
  group_by(gene_set, classification) %>%
  summarise(n = n(), .groups = "drop")

#  Percentage within each gene_set (rows sum to 100
pct_table <- count_table %>%
  group_by(gene_set) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  ungroup()

#  Wide format for easy reading
pct_wide <- pct_table %>%
  select(gene_set, classification, pct) %>%
  tidyr::pivot_wider(names_from  = classification,
                     values_from = pct,
                     values_fill = 0)

library(ggpubr)


# Visualization: stacked bar 
# Gene set display name map
geneset_labels <- c(
  "Common"        = "Shared_Sp8&Sp6Dependent",
  "Common_sp8"    = "Shared_Sp8Dependent",
  "Common_sp6"    = "Shared_Sp6Dependent",
  "redun"         = "Redundant",
  "spe8_dtg8"     = "Sp8_specific",
  "notbound_DEGs" = "DEGs_Not_Bound"
)

scrna_colors <- c(
  "AER_enriched"        = "#7B2D8B",
  "AER_broad"           = "#C084D4",
  "panECT_enriched"     = "#1B7A4A",
  "panECT_broad"        = "#74C496",
  "nonAER_ECT_enriched" = "#0E6B7A",
  "nonAER_ECT_broad"    = "#6DC0CC",
  "ECT_absent"          = "#D3D3D3",
  "Not_in_scRNA"        = "#2C2C2C"
)

ggplot(pct_table %>%
         filter(gene_set %in% names(geneset_labels)) %>%
         mutate(
           gene_set       = factor(geneset_labels[gene_set], levels = geneset_labels),
           classification = factor(classification, levels = rev(names(scrna_colors)))
         ),
       aes(x = gene_set, y = pct, fill = classification)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  geom_text(aes(label = ifelse(pct > 3, paste0(pct, "%"), "")),
            position = position_stack(vjust = 0.5),
            size = 3, fontface = "bold") +
  scale_fill_manual(values = scrna_colors, breaks = names(scrna_colors),
                    guide  = guide_legend(reverse = TRUE)) +
  labs(title = "scRNA classification per gene set",
       x     = "Gene set",
       y     = "Percentage (%)",
       fill  = "Classification") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("revisions/Plots/classification_pct_barplot.png",
       width = 8, height = 5, dpi = 200)

rm(geneset_config)
##################################  limb phenotype annotation


## --- Step 1: Map human -> mouse genes ---


# Start engine UPHENO:0002945 LIMB PHENOTYPE
library(monarchr)
gf <- monarch_engine()
# Step 1: fetch the UPHENO node
gf <- gf |>
  fetch_nodes(query_ids = c("UPHENO:0002945"))
# Step 2: expand to all descendant phenotypes (subclass_of, transitive closure)
gf <- gf |>
  expand(predicates = "biolink:subclass_of", direction = "in", transitive = TRUE)
# Step 3: expand to genes associated with those phenotypes
gf <- gf |>
  expand(categories = "biolink:Gene")
# Extract nodes and edges from igraph
edgesf <- igraph::as_data_frame(gf, what = "edges")
edgesf_f <- edgesf |>
  filter(predicate !="biolink:subclass_of")
# Convert list columns to strings
edgesf_f2 <- data.frame(lapply(edgesf_f, function(x) {
  if (is.list(x)) {
    sapply(x, function(y) paste(y, collapse = ";"))
  } else {
    x
  }
}), stringsAsFactors = FALSE)

write.table(edgesf_f2,
            file = "revisions/monarch_database_limb-fin_phenotype.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

### V2 getting the limb annotaiton from monarch, mgi and hpo

# Start engine ###UPHENO:0047354 disgenesis of the ectoderm
gf <- monarch_engine() 
# Step 1: fetch the UPHENO node
gf_ect <- gf |>
  fetch_nodes(query_ids = c("UPHENO:0047354")) #### disgenesis of the ectoderm
# Step 2: expand to all descendant phenotypes (subclass_of, transitive closure)
gf_ect <- gf_ect |>
  expand(predicates = "biolink:subclass_of", direction = "in", transitive = TRUE)
# Step 3: expand to genes associated with those phenotypes
gf_ect <- gf_ect |>
  expand(categories = "biolink:Gene")
# Extract nodes and edges from igraph
edges_ect <- igraph::as_data_frame(gf_ect, what = "edges")
edges_ect <- edges_ect |>
  filter(predicate !="biolink:subclass_of")
edges_ect_2 <- data.frame(lapply(edges_ect, function(x) {
  if (is.list(x)) {
    sapply(x, function(y) paste(y, collapse = ";"))
  } else {
    x
  }
}), stringsAsFactors = FALSE)

write.table(edges_ect_2,
            file = "monarch_database_ectodermal_phenotype.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

library(biomaRt)


# Assuming limb_df has human gene symbols in "from"
human_genes <- unique(edgesf_f$from)
human_genes_ect <- unique(edges_ect$from)
library(homologene)
mapping <- homologene(edgesf_f$from, inTax = 9606, outTax = 10090)
mapping_ect <- homologene(edges_ect$from, inTax = 9606, outTax = 10090)
# 9606 = human, 10090 = mouse
limb_monarch <- unique(mapping$`10090`)
ect_monarch <- unique(mapping_ect$`10090`)


# mapping has: hgnc_symbol | mgi_symbol
# limb_mouse_genes <- unique(mapping$`10090_ID`)
# download the latest file (large: ~200 MB)
url <- "http://www.informatics.jax.org/downloads/reports/MGI_DiseaseGeneModel.rpt"
dest <- "MGI_PhenoGenoMP.rpt"
download.file(url, dest, mode = "wb")

library(data.table)

# Load file into R
mgi <- fread(dest, sep = "\t", header = FALSE, quote = "")[,c(1,4,7)] 

# The file has many columns; important ones are:
# V1 = MGI Allele ID
# V5 = MGI Gene Symbol
# V6 = MGI Gene ID
# V9 = MP term ID
# V10 = MP term name
# V11 = MP term synonyms

colnames(mgi)<- c("Genesymbol","Disease","Mousesymbol")

# Filter for limb-related phenotypes
limb_mgi <- mgi[grep("limb|digit|arm|leg|paw|hand|foot|dact", Disease, ignore.case = TRUE)]
ect_mgi <- mgi[grep("ectodermal", Disease, ignore.case = TRUE)]
# Unique list of mouse genes with limb phenotypes
limb_mouse_genes <- unique(c(limb_mgi$Mousesymbol, limb_monarch))
limb_ect_genes <- unique(c(ect_monarch, ect_mgi$Mousesymbol))

# Construct the URL to the latest release of HPO
base <- "https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2025-09-01/genes_to_phenotype.txt"


url <- paste0(base)
dest_hpo <- "genes_to_phenotype.txt"

if (!file.exists(dest_hpo)) {
  download.file(url, dest_hpo, mode = "wb")
}

hpo <- fread(dest_hpo, sep = "\t", header = TRUE, quote = "")

# Columns:
# 1: entrez-gene-id
# 2: gene-symbol
# 3: HPO-ID
# 4: HPO-Term-Name
# 5+: disease info

# Filter for limb-related phenotypes
limb_hpo <- hpo[grep("limb|digit|arm|leg|hand|foot|dact", hpo_name, ignore.case = TRUE)]
ect_hpo <- hpo[grep("ectodermal", hpo_name, ignore.case = TRUE)]
# Unique human genes associated with limb phenotypes
limb_human_hpo <- unique(limb_hpo$gene_symbol)
ect_human_hpo <- unique(ect_hpo$gene_symbol)
# Option A: homologene
library(homologene)
mapping_hpo <- homologene(limb_human_hpo, inTax = 9606, outTax = 10090)
mapping_hpo_ect <- homologene(ect_human_hpo, inTax = 9606, outTax = 10090)
# 9606 = human, 10090 = mouse
limb_human_hpo <- unique(mapping_hpo$`10090`)
#merge
limb_mouse_genes <- unique(c(limb_mouse_genes, limb_human_hpo))
limb_ect_genes <- unique(c(limb_ect_genes, mapping_hpo_ect$`10090`))

# Phenotype lookup from your vectors
phenotype_lookup <- data.frame(
  gene = unique(c(unlist(lapply(list(commongb_tss4_common_dtg, commongb_redundant_degdko,
                                     commongb_dtg8_tss4, commongb_dtg6_tss4,
                                     sp8_rich_dtg_tss4, sp6_rich_dtg_tss4),
                                function(df) df$gene)))),
  stringsAsFactors = FALSE
) %>%
  mutate(
    limb_phenotype = gene %in% limb_mouse_genes,
    ect_phenotype  = gene %in% limb_ect_genes,
    phenotype_annot = case_when(
      limb_phenotype & ect_phenotype ~ "Limb & Ectodermal",
      limb_phenotype                 ~ "Limb",
      ect_phenotype                  ~ "Ectodermal",
      TRUE                           ~ "None"
    )
  ) %>%
  dplyr::select(gene, phenotype_annot)

#ADD ANNOTATION TO SOURCE DATA 5
commongb_tss4_common_dtg  <- merge(commongb_tss4_common_dtg  , phenotype_lookup, by = "gene")
commongb_redundant_degdko <- merge(commongb_redundant_degdko , phenotype_lookup, by = "gene")
commongb_dtg8_tss4        <- merge(commongb_dtg8_tss4        , phenotype_lookup, by = "gene")
commongb_dtg6_tss4        <- merge(commongb_dtg6_tss4        , phenotype_lookup, by = "gene")
sp8_rich_dtg_tss4         <- merge(sp8_rich_dtg_tss4         , phenotype_lookup, by = "gene")
sp6_rich_dtg_tss4         <- merge(sp6_rich_dtg_tss4         , phenotype_lookup, by = "gene")

#rewrite
xlsx::write.xlsx(commongb_tss, "revisions/output/Source_data_5_Overlap_analysis.xlsx", sheetName = "Shared_peaks", row.names = F)
xlsx::write.xlsx(sp8_enrich_tss, "revisions/output/Source_data_5_Overlap_analysis.xlsx", sheetName = "Sp8_specific_peaks", append = T, row.names = F)
xlsx::write.xlsx(sp6_enrich_tss, "revisions/output/Source_data_5_Overlap_analysis.xlsx", sheetName = "Sp6_specific_peaks", append = T, row.names = F)
xlsx::write.xlsx(commongb_tss4_common_dtg, "revisions/output/Source_data_5_Overlap_analysis.xlsx", sheetName = "Shared_Sp8_Sp6_Dependent.xlsx",append = T, row.names = F)
xlsx::write.xlsx(commongb_redundant_degdko, "revisions/output/Source_data_5_Overlap_analysis.xlsx", sheetName = "Shared_redundant.xlsx",append = T, row.names = F)
xlsx::write.xlsx(commongb_dtg8_tss4, "revisions/output/Source_data_5_Overlap_analysis.xlsx", sheetName = "Shared_Sp8_Dependent.xlsx",append = T, row.names = F)
xlsx::write.xlsx(commongb_dtg6_tss4, "revisions/output/Source_data_5_Overlap_analysis.xlsx", sheetName = "Shared_Sp6_Dependent.xlsx",append = T, row.names = F)
xlsx::write.xlsx(sp8_rich_dtg_tss4, "revisions/output/Source_data_5_Overlap_analysis.xlsx", sheetName = "Sp8_specific.xlsx",append = T,row.names = F)
xlsx::write.xlsx(sp6_rich_dtg_tss4, "revisions/output/Source_data_5_Overlap_analysis.xlsx", sheetName = "Sp6_specific.xlsx",append = T,row.names = F)


############################ ############################ ############################ ############################ 
############################ ############################ TABLE FOR PHENOTYPES ############################ ############################ 
############################ ############################ ############################ ############################ 

# Optional: Convert GeneSet to a factor to preserve your specific list order
plot_data$GeneSet <- factor(plot_data$GeneSet, levels = names(GS))
# ── Phenotype overlap by gene_set and classification ─────────────────────────
phenotype_table <- master_df %>%
  mutate(
    has_limb = gene %in% limb_mouse_genes,
    has_ecto = gene %in% limb_ect_genes,
    has_both = has_limb & has_ecto,
    phenotype = case_when(
      has_both  ~ "Limb & Ectodermal",
      has_limb  ~ "Limb",
      has_ecto  ~ "Ectodermal",
      TRUE      ~ "None"
    )
  ) %>%
  group_by(gene_set) %>%
  summarise(
    n_total          = n(),
    n_limb           = sum(has_limb),
    n_ecto           = sum(has_ecto),
    n_both           = sum(has_both),
    n_none           = sum(phenotype == "None"),
    pct_limb         = round(n_limb / n_total * 100, 1),
    pct_ecto         = round(n_ecto / n_total * 100, 1),
    pct_both         = round(n_both / n_total * 100, 1),
    pct_none         = round(n_none / n_total * 100, 1),
    .groups = "drop"
  )

print(phenotype_table)

display_df <- phenotype_table %>%
  dplyr::select(gene_set, n_total,
         pct_limb, pct_ecto, pct_both, pct_none) %>%
  arrange(gene_set) %>%
  rename(
    "Gene set"      = gene_set,
    "N"             = n_total,
    "Limb (%)"      = pct_limb,
    "Ecto (%)"      = pct_ecto,
    "Both (%)"      = pct_both,
    "None (%)"      = pct_none
  )

phen_table <- ggtexttable(display_df,
                          rows  = NULL,
                          theme = ttheme(
                            colnames.style = colnames_style(fill  = "#2c2c2c",
                                                            color = "white",
                                                            face  = "bold",
                                                            size  = 16),
                            tbody.style    = tbody_style(fill  = c("#f9f9f9", "white"),
                                                         color = "black",
                                                         size  = 16)
                          ))
phen_table

# ── Stacked barplot: phenotype overlap per gene set ──────────────────────────
geneset_labels <- c(
  "Common"        = "Shared_Sp8&Sp6Dependent",
  "Common_sp8"    = "Shared_Sp8Dependent",
  "Common_sp6"    = "Shared_Sp6Dependent",
  "redun"         = "Redundant",
  "spe8_dtg8"     = "Sp8_specific",
  "notbound_DEGs" = "DEGs_Not_Bound"
)

phenotype_colors <- c(
  "Limb & Ectodermal" = "orange",
  "Limb"              = "gold",
  "Ectodermal"        = "mediumseagreen",
  "None"              = "grey85"
)

# Long-format table with per-gene-set n embedded in x-axis labels
phen_long <- master_df %>%
  mutate(
    has_limb  = gene %in% limb_mouse_genes,
    has_ecto  = gene %in% limb_ect_genes,
    phenotype = case_when(
      has_limb & has_ecto ~ "Limb & Ectodermal",
      has_limb            ~ "Limb",
      has_ecto            ~ "Ectodermal",
      TRUE                ~ "None"
    )
  ) %>%
  filter(gene_set %in% names(geneset_labels)) %>%
  group_by(gene_set) %>%
  mutate(n_total = n()) %>%
  group_by(gene_set, phenotype, n_total) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(
    pct      = round(n / n_total * 100, 1),
    gs_label = factor(
      paste0(geneset_labels[gene_set], "\n(n=", n_total, ")"),
      levels = paste0(geneset_labels[names(geneset_labels)], "\n(n=",
                      n_total[match(names(geneset_labels), gene_set)], ")")
    ),
    phenotype = factor(phenotype, levels = rev(names(phenotype_colors)))
  )

ggplot(phen_long, aes(x = gs_label, y = pct, fill = phenotype)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  geom_text(aes(label = ifelse(pct > 0, paste0(pct, "%"), "")),
            position = position_stack(vjust = 0.5),
            size = 3, fontface = "bold") +
  scale_fill_manual(values = phenotype_colors,
                    breaks = names(phenotype_colors),
                    guide  = guide_legend(reverse = TRUE)) +
  labs(title = "Phenotype overlap per gene set",
       x     = "Gene set",
       y     = "Percentage (%)",
       fill  = "Phenotype") +
  coord_cartesian(ylim = c(0, 100)) +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("revisions/Plots/phenotype_pct_barplot.png",
       width = 8, height = 5, dpi = 600)

############################ ############################ ############################ ############################ 
############################  annotated pathways heatmaps ###################################
############################ ############################ ############################ ############################ 

library(KEGGREST)
library(ComplexHeatmap)
library(circlize)

padj_threshold <- 0.05
lfc_threshold  <- 0.585

# ── 1. GET KEGG GENES ────────────────────────────────────────────────────────

pathways_to_check <- c("04310", "04350", "04010", 
                       "04512","04330")
names(pathways_to_check) <- c("Wnt", "TGFB","Mapk" ,
                              "ECM-receptor interaction","Notch signaling patwhay")
kegg_gene_list <- list()
for (i in seq_along(pathways_to_check)) {
  path_id   <- trimws(as.character(pathways_to_check[i]))
  path_name <- names(pathways_to_check)[i]
  path_id_clean <- sub("^mmu", "", path_id)
  kegg_id       <- paste0("mmu", path_id_clean)
  message("Querying: ", kegg_id)
  pathway      <- keggGet(kegg_id)[[1]]
  gene_info    <- pathway$GENE[seq(2, length(pathway$GENE), 2)]
  gene_symbols <- sub(";.*", "", gene_info)
  if (path_name == "Wnt") {
    # Check to make sure we aren't duplicating it if KEGG updates later
    if (!"Wls" %in% gene_symbols) {
      gene_symbols <- c(gene_symbols, "Wls")
    }
  }
  kegg_gene_list[[path_name]] <- gene_symbols
  message(path_name, ": ", length(gene_symbols), " genes retrieved")
}

# ── 2. PREP DEG TABLES ───────────────────────────────────────────────────────
resord_sp8 <- resord_sp8 %>%
  arrange(padj) %>%
  filter(!duplicated(Gene.name)) %>%
  as.data.frame()
rownames(resord_sp8) <- resord_sp8$Gene.name

resord_sp6 <- resord_sp6 %>%
  arrange(padj) %>%
  filter(!duplicated(Gene.name)) %>%
  as.data.frame()
rownames(resord_sp6) <- resord_sp6$Gene.name

resord_DKO <- resord_DKO %>%
  arrange(padj) %>%
  filter(!duplicated(Gene.name)) %>%
  as.data.frame()
rownames(resord_DKO) <- resord_DKO$Gene.name


deg_sp8_genes <- resord_sp8$Gene.name[resord_sp8$padj < padj_threshold &
                                        abs(resord_sp8$log2FoldChange) > lfc_threshold &
                                        !is.na(resord_sp8$padj)]
deg_sp6_genes <- resord_sp6$Gene.name[resord_sp6$padj < padj_threshold &
                                        abs(resord_sp6$log2FoldChange) > lfc_threshold &
                                        !is.na(resord_sp6$padj)]
deg_DKO_genes <- resord_DKO$Gene.name[resord_DKO$padj < padj_threshold &
                                        abs(resord_DKO$log2FoldChange) > lfc_threshold &
                                        !is.na(resord_DKO$padj)]
any_deg_genes <- union(deg_sp8_genes, union(deg_sp6_genes, deg_DKO_genes))

# ── 3. UPDATE GS WITH notbound_DEGs ──────────────────────────────────────────
# ── Define the six source dataframes and their names ─────────────────────────
source_data_list <- list(
  "Common"        = commongb_tss4_common_dtg,
  "Redundant_DEG" = commongb_redundant_degdko,
  "Common_Sp8"    = commongb_dtg8_tss4,
  "Common_Sp6"    = commongb_dtg6_tss4,
  "Specific_Sp8"  = sp8_rich_dtg_tss4,
  "Specific_Sp6"  = sp6_rich_dtg_tss4
)
# Genes that are DEG but not in any existing GS group

all_GS_genes <- unique(unlist(GS))

# Priority order (highest priority LAST — used to pick when gene can't be repeated)
gs_priority <- c("notbound_DEGs",  "spe6_dtg6", "spe8_dtg8", 
                 "redun","Common_sp6", "Common_sp8", "Common")

# ── 4. COLORS ────────────────────────────────────────────────────────────────
col_lfc <- colorRamp2(c(-3, 0, 3), c("tomato", "white", "steelblue"))

gs_colors <- list(
  Common        = c("0" = "white", "1" = "black"),
  Common_sp8    = c("0" = "white", "1" = "dodgerblue"),
  Common_sp6    = c("0" = "white", "1" = "tomato"),
  spe6_dtg6     = c("0" = "white", "1" = "salmon"),
  spe8_dtg8     = c("0" = "white", "1" = "skyblue"),
  redun         = c("0" = "white", "1" = "grey60"),
  notbound_DEGs = c("0" = "white", "1" = "gold")
)

ann_colors <- list(
  DEG_Sp8 = c("Down in KO" = "steelblue",
              "Up in KO"   = "tomato",
              NS           = "grey85",
              Not_tested   = "white"),
  DEG_Sp6 = c("Down in KO"  = "steelblue",
              "Up in KO"    = "tomato",
              NS            = "grey85",
              Not_tested    = "white",
              Excluded_self = "white"),
  DEG_DKO = c("Down in KO"  = "steelblue",
              "Up in KO"    = "tomato",
              NS            = "grey85",
              Not_tested    = "white",
              Excluded_self = "white"),
  scRNA_class = c(
    AER_enriched        = "#7B2D8B",
    AER_broad          = "#C084D4",
    panECT_enriched     = "#1B7A4A",
    panECT_broad        = "#74C496",
    nonAER_ECT_enriched = "#0E6B7A",
    nonAER_ECT_broad    = "#6DC0CC",
    ECT_absent          = "#D3D3D3",
    Not_in_scRNA        = "#2C2C2C"
  ),

  Phenotype = c("Limb & Ectodermal" = "darkorange",
                Limb                = "gold",
                Ectodermal          = "mediumseagreen",
                None                = "white"),
  Motif = c(
    "GC+AT"    = "#59A14F",
    "AT_only"  = "#F28E2B",
    "GC_only"  = "#4E79A7",
    "No_motif" = "grey75"
  )
)
gs_cat_colors <- c(
  "Common"        = "salmon",
  "Common_sp8"    = "salmon2",
  "Common_sp6"    = "salmon3",
  "spe8_dtg8"     = "mediumorchid4",
  "spe6_dtg6"     = "olivedrab",
  "redun"         = "orange",
  "notbound_DEGs" = "gray80"
)

# Display labels for gene set categories
geneset_labels <- c(
  "Common"        = "Shared_Sp8&Sp6Dependent",
  "Common_sp8"    = "Shared_Sp8Dependent",
  "Common_sp6"    = "Shared_Sp6Dependent",
  "spe6_dtg6"     = "Sp6_specific",
  "spe8_dtg8"     = "Sp8_specific",
  "redun"         = "Redundant",
  "notbound_DEGs" = "DEGs_Not_Bound"
)

# gs_cat_colors mapped to display labels
gs_cat_colors_labeled <- setNames(gs_cat_colors, geneset_labels[names(gs_cat_colors)])

# ── PUBLICATION FONTS & GEOMETRY ─────────────────────────────────────────────
FONT_SIZE_TITLE  <- 14
FONT_SIZE_AXIS   <- 13
FONT_SIZE_ROW    <- 10
FONT_SIZE_LEGEND <- 11
MATRIX_WIDTH_CM  <- 3.5
ROW_HEIGHT_MM    <- 4

# ── 5. LOOP OVER PATHWAYS ────────────────────────────────────────────────────
for (path_name in names(kegg_gene_list)) {
  
  message("\n===== Heatmap for: ", path_name, " =====")
  
  pathway_genes        <- kegg_gene_list[[path_name]]
  genes_in_path_and_GS <- intersect(pathway_genes, all_GS_genes)
  
  if (length(genes_in_path_and_GS) < 2) {
    message("Skipping ", path_name, " — fewer than 2 genes overlap with GS")
    next
  }
  
  # ── 5a. EXPAND GENES FOR ROW SPLIT ────────────────────────────────────────
  gene_rows <- do.call(rbind, lapply(genes_in_path_and_GS, function(g) {
    hits <- names(GS)[sapply(GS, function(x) g %in% x)]
    hits <- hits[hits %in% gs_priority]
    if (length(hits) == 0)
      return(data.frame(gene = g, gs_group = "notbound_DEGs", stringsAsFactors = FALSE))
    data.frame(gene = g, gs_group = hits, stringsAsFactors = FALSE)
  }))
  
  gene_rows$gs_group <- factor(gene_rows$gs_group, levels = rev(gs_priority))
  gene_rows <- gene_rows %>% arrange(gs_group)
  
  # Relabel factor levels for display in row titles and left annotation
  levels(gene_rows$gs_group) <- geneset_labels[levels(gene_rows$gs_group)]
  
  # ── Motif annotation ──────────────────────────────────────────────────────
  # motif_per_gene <- sapply(gene_rows$gene, function(g) {
  #   df <- rbind(source_data_list[[1]][,c(1,6)],
  #               source_data_list[[2]][,c(1,6)],
  #               source_data_list[[3]][,c(1,6)],
  #               source_data_list[[4]][,c(1,6)],
  #               source_data_list[[5]][,c(1,6)],
  #               source_data_list[[6]][,c(1,6)])
  #   peaks_for_gene <- df$peak[df$gene == g]
  #   motifs <- motif_lookup$motif_content[motif_lookup$peak_id %in% peaks_for_gene]
  #   if (length(motifs) == 0) return("No_motif")
  #   names(sort(table(motifs), decreasing = TRUE))[1]
  # })
  
  gene_counts <- table(gene_rows$gene)
  gene_rows$row_label <- ifelse(
    gene_counts[gene_rows$gene] > 1,
    paste0(gene_rows$gene, " (", gene_rows$gs_group, ")"),
    gene_rows$gene
  )
  
  # ── 5b. LFC MATRIX ────────────────────────────────────────────────────────
  lfc_sp8 <- resord_sp8[gene_rows$gene, "log2FoldChange"]
  lfc_sp6 <- resord_sp6[gene_rows$gene, "log2FoldChange"]
  lfc_DKO <- resord_DKO[gene_rows$gene, "log2FoldChange"]
  
  lfc_matrix <- data.frame(
    WTvsSp8KO = ifelse(is.na(lfc_sp8), 0, lfc_sp8),
    WTvsSp6KO = ifelse(is.na(lfc_sp6), 0, lfc_sp6),
    WTvsDKO   = ifelse(is.na(lfc_DKO), 0, lfc_DKO),
    row.names = gene_rows$row_label
  )
  
  # ── 5c. ANNOTATIONS ───────────────────────────────────────────────────────
  deg_sp8_status <- sapply(gene_rows$gene, function(g) {
    if (!g %in% rownames(resord_sp8)) return("Not_tested")
    r <- resord_sp8[g, ]
    if (is.na(r$padj)) return("NS")
    if (r$padj < padj_threshold & r$log2FoldChange >  lfc_threshold) return("Down in KO")
    if (r$padj < padj_threshold & r$log2FoldChange < -lfc_threshold) return("Up in KO")
    return("NS")
  })
  
  deg_sp6_status <- sapply(gene_rows$gene, function(g) {
    if (!g %in% rownames(resord_sp6)) return("Not_tested")
    r <- resord_sp6[g, ]
    if (is.na(r$padj)) return("NS")
    if (r$padj < padj_threshold & r$log2FoldChange >  lfc_threshold) return("Down in KO")
    if (r$padj < padj_threshold & r$log2FoldChange < -lfc_threshold) return("Up in KO")
    return("NS")
  })
  
  deg_dko_status <- sapply(gene_rows$gene, function(g) {
    if (!g %in% rownames(resord_DKO)) return("Not_tested")
    r <- resord_DKO[g, ]
    if (is.na(r$padj)) return("NS")
    if (r$padj < padj_threshold & r$log2FoldChange >  lfc_threshold) return("Down in KO")
    if (r$padj < padj_threshold & r$log2FoldChange < -lfc_threshold) return("Up in KO")
    return("NS")
  })
  
  limb_annot <- ifelse(gene_rows$gene %in% limb_mouse_genes, "Limb", "")
  ecto_annot <- ifelse(gene_rows$gene %in% limb_ect_genes,   "Ectodermal", "")
  phenotype_annot <- case_when(
    limb_annot != "" & ecto_annot != "" ~ "Limb & Ectodermal",
    limb_annot != ""                    ~ "Limb",
    ecto_annot != ""                    ~ "Ectodermal",
    TRUE                                ~ "None"
  )
  
  # ── 5d. ANNOTATIONS ───────────────────────────────────────────────────────
  row_ann <- rowAnnotation(
    DEG_Sp8   = deg_sp8_status,
    DEG_Sp6   = deg_sp6_status,
    DEG_DKO   = deg_dko_status,
    Phenotype = phenotype_annot,
    col = c(gs_colors, ann_colors[c("DEG_Sp8", "DEG_Sp6", "DEG_DKO", "Phenotype")]),
    border = TRUE,
    annotation_legend_param = list(
      DEG_Sp8   = list(title = "DEG Sp8KO"),
      DEG_Sp6   = list(title = "DEG Sp6KO"),
      DEG_DKO   = list(title = "DEG DKO"),
      Phenotype = list(title = "Phenotype")
    )
  )
  
  left_ann <- rowAnnotation(
    Category = as.character(gene_rows$gs_group),
    col      = list(Category = gs_cat_colors_labeled),
    show_legend = TRUE,
    show_annotation_name = FALSE,
    annotation_legend_param = list(Category = list(title = "Gene category")),
    width = unit(3, "mm")
  )
  
  # ── 5e. HEATMAP ───────────────────────────────────────────────────────────
  ht <- Heatmap(
    as.matrix(lfc_matrix),
    name             = "L2FC",
    col              = col_lfc,
    cluster_rows     = TRUE,
    cluster_columns  = FALSE,
    show_row_names   = TRUE,
    left_annotation  = left_ann,
    right_annotation = row_ann,
    row_split        = gene_rows$gs_group,
    cluster_row_slices = FALSE,
    column_title     = paste0(path_name, " — L2FC: Sp8KO vs WT  |  Sp6KO vs WT  |  DKO vs WT"),
    border = TRUE,
    heatmap_legend_param = list(
      title      = "log2 FC",
      title_gp   = gpar(fontsize = FONT_SIZE_LEGEND, fontface = "bold"),
      labels_gp  = gpar(fontsize = FONT_SIZE_LEGEND),
      at         = c(-2.5, 0, 2.5)
    ),
    width           = unit(MATRIX_WIDTH_CM, "cm"),
    height          = unit(ROW_HEIGHT_MM * nrow(lfc_matrix), "mm"),
    row_names_gp    = gpar(fontsize = FONT_SIZE_ROW),
    column_names_gp = gpar(fontsize = FONT_SIZE_AXIS, fontface = "bold"),
    column_names_rot = 45,
    row_title_gp    = gpar(fontsize = FONT_SIZE_AXIS, fontface = "bold"),
    column_title_gp = gpar(fontsize = FONT_SIZE_TITLE, fontface = "bold")
  )
  
  # ── 5f. SAVE ──────────────────────────────────────────────────────────────
  safe_name   <- gsub("[^A-Za-z0-9_]", "_", path_name)
  out_file    <- paste0("revisions/Plots/KEGG_heatmap-2.5_sc", safe_name, ".pdf")
  num_rows    <- nrow(lfc_matrix)
  calc_height <- (num_rows * ROW_HEIGHT_MM) / 25.4 + 3.5
  
  pdf(out_file, height = max(7, calc_height))
  draw(ht, merge_legend = TRUE)
  dev.off()
  message("Saved: ", out_file)
}

# ############################################### ANNOTATE genes in all CATEGORIES ################
# 
# 
# # ── Loop ─────────────────────────────────────────────────────────────────────
# for (df_name in names(source_data_list)) {
#   
#   message("\n===== Heatmap for: ", df_name, " =====")
#   
#   # Get genes from the dataframe
#   genes_in_set <- unique(source_data_list[[df_name]]$gene)
#   genes_in_set <- genes_in_set[!is.na(genes_in_set)]
#   
#   if (length(genes_in_set) < 2) {
#     message("Skipping ", df_name, " — fewer than 2 genes")
#     next
#   }
#   
#   # ── EXPAND GENES FOR ROW SPLIT ────────────────────────────────────────────
#   gene_rows <- do.call(rbind, lapply(genes_in_set, function(g) {
#     hits <- names(GS)[sapply(GS, function(x) g %in% x)]
#     hits <- hits[hits %in% gs_priority]
#     if (length(hits) == 0) {
#       return(data.frame(gene = g, gs_group = "notbound_DEGs", stringsAsFactors = FALSE))
#     }
#     data.frame(gene = g, gs_group = hits, stringsAsFactors = FALSE)
#   }))
#   
#   gene_rows$gs_group <- factor(gene_rows$gs_group, levels = rev(gs_priority))
#   gene_rows <- gene_rows %>% arrange(gs_group)
#   
#   gene_counts <- table(gene_rows$gene)
#   gene_rows$row_label <- ifelse(
#     gene_counts[gene_rows$gene] > 1,
#     paste0(gene_rows$gene, " (", gene_rows$gs_group, ")"),
#     gene_rows$gene
#   )
#   motif_per_gene <- sapply(gene_rows$gene, function(g) {
#     # Get all peaks associated with this gene in the current source df
#     
#     df <- rbind(source_data_list[[1]][,c(1,6)],
#                 source_data_list[[2]][,c(1,3)],
#                 source_data_list[[3]][,c(1,6)],
#                 source_data_list[[4]][,c(1,6)],
#                 source_data_list[[5]][,c(1,6)],
#                 source_data_list[[6]][,c(1,6)])
#     peaks_for_gene <- df$peak[df$gene == g]
#     motifs <- motif_lookup$motif_content[motif_lookup$peak_id %in% peaks_for_gene]
#     if (length(motifs) == 0) return("No_motif")
#     # Return most frequent motif type
#     names(sort(table(motifs), decreasing = TRUE))[1]
#   })
#   # ── LFC MATRIX ────────────────────────────────────────────────────────────
#   lfc_sp8 <- resord_sp8[gene_rows$gene, "log2FoldChange"]
#   lfc_sp6 <- resord_sp6[gene_rows$gene, "log2FoldChange"]
#   lfc_DKO <- resord_DKO[gene_rows$gene, "log2FoldChange"]
#   
#   lfc_matrix <- data.frame(
#     WTvsSp8KO = ifelse(is.na(lfc_sp8), 0, lfc_sp8),
#     WTvsSp6KO = ifelse(is.na(lfc_sp6), 0, lfc_sp6),
#     WTvsDKO   = ifelse(is.na(lfc_DKO), 0, lfc_DKO),
#     row.names = gene_rows$row_label
#   )
#   
#   # ── ANNOTATIONS ───────────────────────────────────────────────────────────
#   deg_sp8_status <- sapply(gene_rows$gene, function(g) {
#     if (!g %in% rownames(resord_sp8)) return("Not_tested")
#     r <- resord_sp8[g, ]
#     if (is.na(r$padj)) return("NS")
#     if (r$padj < padj_threshold & r$log2FoldChange >  lfc_threshold) return("Activated")
#     if (r$padj < padj_threshold & r$log2FoldChange < -lfc_threshold) return("Repressed")
#     return("NS")
#   })
#   
#   deg_sp6_status <- sapply(gene_rows$gene, function(g) {
#     if (!g %in% rownames(resord_sp6)) return("Not_tested")
#     r <- resord_sp6[g, ]
#     if (is.na(r$padj)) return("NS")
#     if (r$padj < padj_threshold & r$log2FoldChange >  lfc_threshold) return("Activated")
#     if (r$padj < padj_threshold & r$log2FoldChange < -lfc_threshold) return("Repressed")
#     return("NS")
#   })
#   
#   deg_dko_status <- sapply(gene_rows$gene, function(g) {
#     if (!g %in% rownames(resord_DKO)) return("Not_tested")
#     r <- resord_DKO[g, ]
#     if (is.na(r$padj)) return("NS")
#     if (r$padj < padj_threshold & r$log2FoldChange >  lfc_threshold) return("Activated")
#     if (r$padj < padj_threshold & r$log2FoldChange < -lfc_threshold) return("Repressed")
#     return("NS")
#   })
#   
#   scrna_class <- master_df$classification[match(gene_rows$gene, master_df$gene)]
#   scrna_class[is.na(scrna_class)] <- "Not_in_scRNA"
#   
#   limb_annot <- ifelse(gene_rows$gene %in% limb_mouse_genes, "Limb", "")
#   ecto_annot <- ifelse(gene_rows$gene %in% limb_ect_genes,   "Ectodermal", "")
#   phenotype_annot <- case_when(
#     limb_annot != "" & ecto_annot != "" ~ "Limb & Ectodermal",
#     limb_annot != ""                    ~ "Limb",
#     ecto_annot != ""                    ~ "Ectodermal",
#     TRUE                                ~ "None"
#   )
#   
#   # Binary GS bars
#   gs_binary_chr <- as.data.frame(lapply(gs_priority, function(gs) {
#     as.character(as.integer(gene_rows$gene %in% GS[[gs]]))
#   }))
#   colnames(gs_binary_chr) <- gs_priority
#   rownames(gs_binary_chr) <- gene_rows$row_label
#   
#   # ── ROW ANNOTATION ────────────────────────────────────────────────────────
#   row_ann <- rowAnnotation(
#     DEG_Sp8     = deg_sp8_status,
#     DEG_Sp6     = deg_sp6_status,
#     DEG_DKO     = deg_dko_status,
#     scRNA_class = scrna_class,
#     Phenotype   = phenotype_annot,
#     Motif       = motif_per_gene, 
#     col = c(gs_colors, ann_colors[c("DEG_Sp8", "DEG_Sp6", "DEG_DKO", "scRNA_class", "Phenotype","Motif")]),
#     border = TRUE,
#     annotation_legend_param = list(
#       DEG_Sp8     = list(title = "DEG Sp8KO"),
#       DEG_Sp6     = list(title = "DEG Sp6KO"),
#       DEG_DKO     = list(title = "DEG DKO"),
#       scRNA_class = list(title = "scRNA class"),
#       Phenotype   = list(title = "Phenotype"),
#       Motif       = list(title = "Motif content")
#     )
#   )
#   
#   # ── HEATMAP ───────────────────────────────────────────────────────────────
#   ht <- Heatmap(
#     as.matrix(lfc_matrix),
#     name             = "L2FC",
#     col              = col_lfc,
#     cluster_rows     = TRUE,
#     cluster_columns  = FALSE,
#     show_row_names   = TRUE,
#     right_annotation = row_ann,
#     row_split        = gene_rows$gs_group,
#     cluster_row_slices = FALSE,
#     column_title     = paste0(df_name, " — L2FC: Sp8KO vs WT  |  Sp6KO vs WT  |  DKO vs WT"),
#     
#     border = TRUE,
#     heatmap_legend_param = list(title = "log2 FC", title_gp = gpar(fontsize = FONT_SIZE_LEGEND, fontface = "bold"),
#                                 labels_gp = gpar(fontsize = FONT_SIZE_LEGEND),at =  c(-2.5, 0, 2.5)),
#     # Font Size Updates
#     row_height      = unit(ROW_HEIGHT_MM, "mm"), 
#     row_names_gp    = gpar(fontsize = FONT_SIZE_ROW),
#     column_names_gp = gpar(fontsize = FONT_SIZE_AXIS, fontface = "bold"),
#     row_title_gp    = gpar(fontsize = FONT_SIZE_AXIS, fontface = "bold"),
#     column_title_gp = gpar(fontsize = FONT_SIZE_TITLE, fontface = "bold")
#   )
#   
#   # ── SAVE ──────────────────────────────────────────────────────────────────
#   safe_name <- gsub("[^A-Za-z0-9_]", "_", df_name)
#   out_file  = paste0("revisions/Plots/Binding_heatmap-2.5_", safe_name, ".png")
#   
#   num_rows <- nrow(lfc_matrix)
#   calc_height <- (num_rows * ROW_HEIGHT_MM) / 25.4 + 3.5 
#   
#   png(out_file,
#       width  = 12,
#       height = max(7, calc_height),
#       units  = "in", res = 300) # High DPI output
#   draw(ht, merge_legend = TRUE)
#   dev.off()
#   message("Saved: ", out_file)
# }


################################################################################
############################# Custom gene heatmap ##############################
################################################################################

custom_genes <- c("Sost", "Dlx1", "Cyp26b1", "Rspo2", "Fzd1", "Bmp4",
                  "En1", "Wnt7a", "Bhlha9", "Fgf8", "Sall1",
                  "Msx1","Msx2","Bambi","Wnt5a","Sp6","Dlx2","Dlx3","Dlx4","Dlx5","Dlx6","Jag2")

# Build gene_rows (same logic as other loops)
gene_rows <- do.call(rbind, lapply(custom_genes, function(g) {
  hits <- names(GS)[sapply(GS, function(x) g %in% x)]
  hits <- hits[hits %in% gs_priority]
  if (length(hits) == 0) {
    return(data.frame(gene = g, gs_group = "notbound_DEGs", stringsAsFactors = FALSE))
  }
  data.frame(gene = g, gs_group = hits, stringsAsFactors = FALSE)
}))

gene_rows$gs_group <- factor(gene_rows$gs_group, levels = rev(gs_priority))
gene_rows <- gene_rows %>% arrange(gs_group)

gene_counts <- table(gene_rows$gene)
gene_rows$row_label <- ifelse(
  gene_counts[gene_rows$gene] > 1,
  paste0(gene_rows$gene, " (", gene_rows$gs_group, ")"),
  gene_rows$gene
)

motif_per_gene <- sapply(gene_rows$gene, function(g) {
  for (sdf in source_data_list) {
    if (!"motif_content" %in% colnames(sdf) || !"gene" %in% colnames(sdf)) next
    motifs <- sdf$motif_content[sdf$gene == g]
    motifs <- motifs[!is.na(motifs) & motifs != "No_motif"]
    if (length(motifs) > 0)
      return(names(sort(table(motifs), decreasing = TRUE))[1])
  }
  "No_motif"
})

# LFC matrix
lfc_sp8 <- resord_sp8[gene_rows$gene, "log2FoldChange"]
lfc_sp6 <- resord_sp6[gene_rows$gene, "log2FoldChange"]
lfc_DKO <- resord_DKO[gene_rows$gene, "log2FoldChange"]

lfc_matrix <- data.frame(
  WTvsSp8KO = ifelse(is.na(lfc_sp8), 0, lfc_sp8),
  WTvsSp6KO = ifelse(is.na(lfc_sp6), 0, lfc_sp6),
  WTvsDKO   = ifelse(is.na(lfc_DKO), 0, lfc_DKO),
  row.names = gene_rows$row_label
)

# Annotations
deg_sp8_status <- sapply(gene_rows$gene, function(g) {
  if (!g %in% rownames(resord_sp8)) return("Not_tested")
  r <- resord_sp8[g, ]
  if (is.na(r$padj)) return("NS")
  if (r$padj < padj_threshold & r$log2FoldChange >  lfc_threshold) return("Down in KO")
  if (r$padj < padj_threshold & r$log2FoldChange < -lfc_threshold) return("Up in KO")
  return("NS")
})

deg_sp6_status <- sapply(gene_rows$gene, function(g) {
  if (!g %in% rownames(resord_sp6)) return("Not_tested")
  r <- resord_sp6[g, ]
  if (is.na(r$padj)) return("NS")
  if (r$padj < padj_threshold & r$log2FoldChange >  lfc_threshold) return("Down in KO")
  if (r$padj < padj_threshold & r$log2FoldChange < -lfc_threshold) return("Up in KO")
  return("NS")
})

deg_dko_status <- sapply(gene_rows$gene, function(g) {
  if (!g %in% rownames(resord_DKO)) return("Not_tested")
  r <- resord_DKO[g, ]
  if (is.na(r$padj)) return("NS")
  if (r$padj < padj_threshold & r$log2FoldChange >  lfc_threshold) return("Down in KO")
  if (r$padj < padj_threshold & r$log2FoldChange < -lfc_threshold) return("Up in KO")
  return("NS")
})
deg_sp6_status[gene_rows$gene == "Sp6"] <- "Excluded_self"
deg_dko_status[gene_rows$gene == "Sp6"] <- "Excluded_self"

scrna_class <- master_df$classification[match(gene_rows$gene, master_df$gene)]
scrna_class[is.na(scrna_class)] <- "Not_in_scRNA"

limb_annot <- ifelse(gene_rows$gene %in% limb_mouse_genes, "Limb", "")
ecto_annot <- ifelse(gene_rows$gene %in% limb_ect_genes,   "Ectodermal", "")
phenotype_annot <- case_when(
  limb_annot != "" & ecto_annot != "" ~ "Limb & Ectodermal",
  limb_annot != ""                    ~ "Limb",
  ecto_annot != ""                    ~ "Ectodermal",
  TRUE                                ~ "None"
)

deg_col_self <- c("Down in KO"  = "steelblue",
                  "Up in KO"    = "tomato",
                  NS            = "grey85",
                  Not_tested    = "white",
                  Excluded_self = "white")

row_ann <- rowAnnotation(
  DEG_Sp8     = deg_sp8_status,
  DEG_Sp6     = deg_sp6_status,
  DEG_DKO     = deg_dko_status,
  scRNA_class = scrna_class,
  Phenotype   = phenotype_annot,
  Motif       = motif_per_gene,
  col = c(gs_colors,
          list(DEG_Sp8     = ann_colors$DEG_Sp8,
               DEG_Sp6     = deg_col_self,
               DEG_DKO     = deg_col_self,
               scRNA_class = ann_colors$scRNA_class,
               Phenotype   = ann_colors$Phenotype,
               Motif       = ann_colors$Motif)),
  border = TRUE,
  annotation_legend_param = list(
    DEG_Sp8     = list(title = "DEG Sp8KO"),
    DEG_Sp6     = list(title = "DEG Sp6KO"),
    DEG_DKO     = list(title = "DEG DKO"),
    scRNA_class = list(title = "scRNA class"),
    Phenotype   = list(title = "Phenotype"),
    Motif       = list(title = "Motif content")
  )
)

left_ann <- rowAnnotation(
  Category = geneset_labels[as.character(gene_rows$gs_group)],
  col      = list(Category = gs_cat_colors_labeled),
  show_legend = TRUE,
  show_annotation_name = FALSE,
  annotation_legend_param = list(Category = list(title = "Gene category")),
  width = unit(3, "mm")
)

ht_custom <- Heatmap(
  as.matrix(lfc_matrix),
  name             = "L2FC",
  col              = col_lfc,
  cluster_rows     = TRUE,
  cluster_columns  = FALSE,
  show_row_names   = TRUE,
  left_annotation  = left_ann,
  right_annotation = row_ann,
  row_split        = gene_rows$gs_group,
  cluster_row_slices = FALSE,
  column_title     = "Custom genes — L2FC: WT vs Sp8KO  |  WT vs Sp6KO  |  WT vs DKO",
  border = TRUE,
  heatmap_legend_param = list(title = "log2 FC", title_gp = gpar(fontsize = FONT_SIZE_LEGEND, fontface = "bold"),
                              labels_gp = gpar(fontsize = FONT_SIZE_LEGEND), at = c(-2.5, 0, 2.5)),
  width           = unit(MATRIX_WIDTH_CM, "cm"),
  height          = unit(ROW_HEIGHT_MM * nrow(lfc_matrix), "mm"),
  row_names_gp    = gpar(fontsize = FONT_SIZE_ROW),
  column_names_gp = gpar(fontsize = FONT_SIZE_AXIS, fontface = "bold"),
  column_names_rot = 45,
  row_title_gp    = gpar(fontsize = FONT_SIZE_AXIS, fontface = "bold"),
  column_title_gp = gpar(fontsize = FONT_SIZE_TITLE, fontface = "bold")
)

num_rows    <- nrow(lfc_matrix)
calc_height <- (num_rows * ROW_HEIGHT_MM) / 25.4 + 3.5

pdf("revisions/Plots/custom_genes_heatmap.pdf",
    width  = 12,
    height = max(8, calc_height + 1))
draw(ht_custom, merge_legend = TRUE)
dev.off()
message("Saved: revisions/Plots/custom_genes_heatmap.pdf")

################################################################################
#################### Source Data 7 — KEGG pathway heatmap data ################
################################################################################

# ── Helper: classify DEG status ──────────────────────────────────────────────
deg_status <- function(lfc, pval,
                       lfc_thr = lfc_threshold, p_thr = padj_threshold) {
  case_when(
    is.na(pval)                              ~ "Not_tested",
    pval < p_thr & lfc >  lfc_thr           ~ "Down in KO",
    pval < p_thr & lfc < -lfc_thr          ~ "Up in KO",
    TRUE                                     ~ "NS"
  )
}

# ── Normalise each source df to a common schema ───────────────────────────────
# Common / Common_Sp8 / Specific_Sp8/Sp6 already have standardised names;
# commongb_redundant_degdko uses its own naming and already has DKO info.

norm_source <- function(df, nm) {
  d <- df

  # commongb_redundant_degdko: rename to standard names
  if ("dist" %in% colnames(d))
    d <- dplyr::rename(d, distance_to_genetss = dist)
  if ("log2FoldChange.sp8" %in% colnames(d))
    d <- dplyr::rename(d, L2FC_Sp8 = log2FoldChange.sp8, padj_Sp8 = padj.sp8) %>%
      dplyr::select(-any_of("Category.sp8")) %>%
      mutate(status_Sp8 = deg_status(L2FC_Sp8, padj_Sp8))
  if ("log2FoldChange.sp6" %in% colnames(d))
    d <- dplyr::rename(d, L2FC_Sp6 = log2FoldChange.sp6, padj_Sp6 = padj.sp6) %>%
      dplyr::select(-any_of("Category.sp6")) %>%
      mutate(status_Sp6 = deg_status(L2FC_Sp6, padj_Sp6))
  # DKO info already present in redundant df
  if ("log2FoldChange_dko" %in% colnames(d))
    d <- dplyr::rename(d,
      L2FC_DKO   = log2FoldChange_dko,
      padj_DKO   = padj_dko) %>%
      dplyr::select(-any_of("status_dko")) %>%
      mutate(status_DKO = deg_status(L2FC_DKO, padj_DKO))

  if ("L2FC_Sp8" %in% colnames(d) && "padj_Sp8" %in% colnames(d))
    d <- d %>% mutate(status_Sp8 = deg_status(L2FC_Sp8, padj_Sp8))
  if ("L2FC_Sp6" %in% colnames(d) && "padj_Sp6" %in% colnames(d))
    d <- d %>% mutate(status_Sp6 = deg_status(L2FC_Sp6, padj_Sp6))
  if ("L2FC_DKO" %in% colnames(d) && "padj_DKO" %in% colnames(d))
    d <- d %>% mutate(status_DKO = deg_status(L2FC_DKO, padj_DKO))

  d %>% mutate(GS_category = nm)
}

sp8_deg_info <- resord_sp8 %>%
  dplyr::select(Gene.name, log2FoldChange, padj) %>%
  dplyr::rename(gene = Gene.name, L2FC_Sp8 = log2FoldChange, padj_Sp8 = padj) %>%
  distinct(gene, .keep_all = TRUE) %>%
  mutate(status_Sp8 = deg_status(L2FC_Sp8, padj_Sp8))

sp6_deg_info <- resord_sp6 %>%
  dplyr::select(Gene.name, log2FoldChange, padj) %>%
  dplyr::rename(gene = Gene.name, L2FC_Sp6 = log2FoldChange, padj_Sp6 = padj) %>%
  distinct(gene, .keep_all = TRUE) %>%
  mutate(status_Sp6 = deg_status(L2FC_Sp6, padj_Sp6))

dko_deg_info <- resord_DKO %>%
  dplyr::select(Gene.name, log2FoldChange, padj) %>%
  dplyr::rename(gene = Gene.name, L2FC_DKO = log2FoldChange, padj_DKO = padj) %>%
  distinct(gene, .keep_all = TRUE) %>%
  mutate(status_DKO = deg_status(L2FC_DKO, padj_DKO))

all_source_combined <- bind_rows(
  lapply(names(source_data_list), function(nm)
    norm_source(source_data_list[[nm]], nm))
) %>%
  # Fill missing Sp8 DEG columns (sp6_rich has none)
  left_join(sp8_deg_info %>% dplyr::select(gene, L2FC_Sp8, status_Sp8, padj_Sp8),
            by = "gene", suffix = c("", ".fill")) %>%
  mutate(
    L2FC_Sp8   = dplyr::coalesce(L2FC_Sp8,   L2FC_Sp8.fill),
    status_Sp8 = dplyr::coalesce(status_Sp8, status_Sp8.fill),
    padj_Sp8   = dplyr::coalesce(padj_Sp8,   padj_Sp8.fill)
  ) %>%
  dplyr::select(-ends_with(".fill")) %>%
  # Fill missing Sp6 DEG columns (sp8_rich has none)
  left_join(sp6_deg_info %>% dplyr::select(gene, L2FC_Sp6, status_Sp6, padj_Sp6),
            by = "gene", suffix = c("", ".fill")) %>%
  mutate(
    L2FC_Sp6   = dplyr::coalesce(L2FC_Sp6,   L2FC_Sp6.fill),
    status_Sp6 = dplyr::coalesce(status_Sp6, status_Sp6.fill),
    padj_Sp6   = dplyr::coalesce(padj_Sp6,   padj_Sp6.fill)
  ) %>%
  dplyr::select(-ends_with(".fill")) %>%
  # Fill missing DKO columns (all except redundant_degdko)
  left_join(dko_deg_info %>% dplyr::select(gene, L2FC_DKO, status_DKO, padj_DKO),
            by = "gene", suffix = c("", ".fill")) %>%
  mutate(
    L2FC_DKO   = dplyr::coalesce(L2FC_DKO,   L2FC_DKO.fill),
    status_DKO = dplyr::coalesce(status_DKO, status_DKO.fill),
    padj_DKO   = dplyr::coalesce(padj_DKO,   padj_DKO.fill)
  ) %>%
  dplyr::select(-ends_with(".fill"))

# ── Preferred column order ────────────────────────────────────────────────────
preferred_cols <- c(
  "KEGG_pathway", "gene", "GS_category", "GS_category_label",
  "peak", "chr", "start", "end", "distance_to_genetss",
  "Distance_group", "Distance_label",
  "L2FC_Sp8", "status_Sp8", "padj_Sp8",
  "L2FC_Sp6", "status_Sp6", "padj_Sp6",
  "L2FC_DKO", "status_DKO", "padj_DKO",
  "motif_content", "classification", "phenotype_annot"
)

# ── One sheet per KEGG pathway ────────────────────────────────────────────────
out_file_sd7 <- "revisions/output/Source_Data_7_KEGG_heatmaps.xlsx"
first_sheet  <- TRUE

for (path_name in names(kegg_gene_list)) {
  pathway_genes        <- kegg_gene_list[[path_name]]
  genes_in_path_and_GS <- intersect(pathway_genes, all_GS_genes)
  if (length(genes_in_path_and_GS) < 2) next

  gene_gs <- do.call(rbind, lapply(genes_in_path_and_GS, function(g) {
    hits <- names(GS)[sapply(GS, function(x) g %in% x)]
    hits <- hits[hits %in% gs_priority]
    if (length(hits) == 0)
      return(data.frame(gene = g, gs_group = "notbound_DEGs", stringsAsFactors = FALSE))
    data.frame(gene = g, gs_group = hits, stringsAsFactors = FALSE)
  })) %>%
    mutate(GS_category_label = geneset_labels[gs_group])

  bound_genes    <- unique(gene_gs$gene[gene_gs$gs_group != "notbound_DEGs"])
  notbound_genes <- unique(gene_gs$gene[gene_gs$gs_group == "notbound_DEGs"])

  # Bound genes: full source data rows
  sheet_bound <- all_source_combined %>%
    filter(gene %in% bound_genes) %>%
    left_join(
      gene_gs %>% dplyr::select(gene, GS_category_label) %>% distinct(),
      by = "gene"
    ) %>%
    mutate(KEGG_pathway = path_name)

  # notbound_DEGs: DEG info only, no peak data
  sheet_notbound <- if (length(notbound_genes) > 0) {
    sp8_deg_info %>%
      filter(gene %in% notbound_genes) %>%
      left_join(sp6_deg_info, by = "gene") %>%
      left_join(dko_deg_info, by = "gene") %>%
      left_join(
        master_df %>% dplyr::select(gene, classification) %>% distinct(),
        by = "gene"
      ) %>%
      left_join(
        phenotype_lookup %>% dplyr::rename(phenotype_annot = phenotype_annot),
        by = "gene"
      ) %>%
      mutate(
        GS_category       = "notbound_DEGs",
        GS_category_label = "DEGs_Not_Bound",
        KEGG_pathway      = path_name
      )
  } else NULL

  sheet_data <- bind_rows(sheet_bound, sheet_notbound) %>%
    dplyr::select(any_of(preferred_cols))

  safe_name <- substr(gsub("[^A-Za-z0-9_]", "_", path_name), 1, 31)
  xlsx::write.xlsx(as.data.frame(sheet_data),
                   out_file_sd7,
                   sheetName = safe_name,
                   append    = !first_sheet,
                   row.names = FALSE)
  first_sheet <- FALSE
  message("Written sheet: ", safe_name)
}

message("Saved: ", out_file_sd7)


###################### vista ###############################
mm10_db_info <- read.delim("revisions/experiments.tsv",header = T, na.strings = c("","NA") )

mm10_db_info <- mm10_db_info %>%
  mutate(
    chr = str_extract(coordinate_mm10, "^[^:]+"),
    start = as.integer(str_extract(coordinate_mm10, "(?<=:)[0-9]+")),
    end = as.integer(str_extract(coordinate_mm10, "(?<=-)[0-9]+"))
  )

mm10_db_info<- mm10_db_info[!is.na(mm10_db_info$start) | !is.na(mm10_db_info$end),]

mm10_db_gr <- GenomicRanges::GRanges(seqnames = mm10_db_info$chr, ranges = IRanges(start = mm10_db_info$start, end = mm10_db_info$end), 
                                     id=mm10_db_info$vista_id, curation = mm10_db_info$curation_status )

common <- read.table("revisions/seqminer_redo/output/common_pool.bed")
common_gr <- GenomicRanges::GRanges(seqnames = common$V1, ranges = IRanges(start =common$V2, end = common$V3) )

# Find overlapping regions
hits <- findOverlaps(mm10_db_gr, common_gr)
# Extract intersected ranges
overlapping_ranges <- pintersect(mm10_db_gr[queryHits(hits)], common_gr[subjectHits(hits)])
# Preserve metadata from gr1
mcols(overlapping_ranges) <- mcols(mm10_db_gr)[queryHits(hits), , drop = FALSE]
# Print result
print(overlapping_ranges)
overlapping_ranges <- as.data.frame(overlapping_ranges)
overlapping_ranges$source <- "Common"

sp8_spe <- read.table("revisions/seqminer_redo/output/sp8_spe_pool.bed")
sp8_spe_gr <- GenomicRanges::GRanges(seqnames = sp8_spe$V1, ranges = IRanges(start =sp8_spe$V2, end = sp8_spe$V3) )


# Find overlapping regions
hits8 <- findOverlaps(mm10_db_gr, sp8_spe_gr)
# Extract intersected ranges
overlapping_ranges8 <- pintersect(mm10_db_gr[queryHits(hits8)], sp8_spe_gr[subjectHits(hits8)])
# Preserve metadata from gr1
mcols(overlapping_ranges8) <- mcols(mm10_db_gr)[queryHits(hits8), , drop = FALSE]
# Print result
print(overlapping_ranges8)
overlapping_ranges8 <- as.data.frame(overlapping_ranges8)
overlapping_ranges8$source <- "Sp8"


sp6_spe <- read.table("revisions/seqminer_redo/output/sp6_spe_pool.bed")
sp6_spe_gr <- GenomicRanges::GRanges(seqnames = sp6_spe$V1, ranges = IRanges(start =sp6_spe$V2, end = sp6_spe$V3) )


# Find overlapping regions
hits6 <- findOverlaps(mm10_db_gr, sp6_spe_gr)
# Extract intersected ranges
overlapping_ranges6 <- pintersect(mm10_db_gr[queryHits(hits6)], sp6_spe_gr[subjectHits(hits6)])
# Preserve metadata from gr1
mcols(overlapping_ranges6) <- mcols(mm10_db_gr)[queryHits(hits6), , drop = FALSE]
# Print result
print(overlapping_ranges6)
overlapping_ranges6 <- as.data.frame(overlapping_ranges6)
overlapping_ranges6$source <- "Sp6"


all_hits <- rbind(overlapping_ranges, overlapping_ranges8,overlapping_ranges6)
all_hits <- all_hits %>%
  mutate(peakid=paste0(seqnames,":",start,"-",end))
all_hits_db <- merge(all_hits, mm10_db_info, by.x="id", by.y="vista_id" )

writexl::write_xlsx(all_hits_db, "vista_hits.xlsx")

# ==============================================================================
# 1. Standardize and Flatten the Gene List (GS) for Hierarchy
# ==============================================================================
target_categories <- c("Common", "Common_sp8", "Common_sp6","redun", "spe8_dtg8", "spe6_dtg6")

gs_lookup <- imap_dfr(GS, ~ tibble(gene_symbol = .x, raw_category = .y)) %>%
  filter(raw_category %in% target_categories) %>%
  mutate(raw_category = factor(raw_category, levels = target_categories)) %>%
  arrange(gene_symbol, raw_category) %>%
  distinct(gene_symbol, .keep_all = TRUE) %>%
  mutate(raw_category = as.character(raw_category))

# Manual curation reference data
manual_lookup <- tibble::tribble(
  ~vista_id, ~manual_genes,     ~Lacz_AERorECT,
  "mm1663",  "Slc29a3,Unc5b",   "Yes",
  "mm264",   "Tex14",           "No",
  "mm876",   "Copz2,Cdk5rap3",  "No",
  "mm265",   "Sox4,Prl5a1",     "Yes?",
  "mm84",    "Ttll12,Scube1",   "Yes",
  "mm1599",  "Spry1,Ankrd50",   "Yes",
  "mm1616",  "Cd1d1,Kirrel",    "Yes",
  "mm1537",  "Slc16a1,Lrig2",   "Yes",
  "mm426",   "Msx1",            "Yes",
  "mm723",   "Pxn,Sirt4",       "negative",
  "mm1539",  "Nr2f2",           "negative",
  "mm907",   "Erich1,Tdrp",     "negative",
  "mm1598",  "Zfp703,Thap1",    "Yes",
  "mm2141",  "Smad3,Smad6",     "Yes",
  "mm698",   "Bcl11a",          "negative",
  "mm1996",  "Rec114,Hcn4",     "negative",
  "hs422",   "Dlx1",            "Yes",
  "hs484",   "Zfp503",          "Yes"
)

# ==============================================================================
# 2. Integrated Data Processing Pipeline
# ==============================================================================
processed_hits <- all_hits_db %>%
  left_join(manual_lookup, by = c("id" = "vista_id")) %>%
  
  mutate(cleaned_genes = bracket_mm10) %>%
  mutate(cleaned_genes = str_remove_all(cleaned_genes, "\\s*\\([^)]+\\)")) %>%
  mutate(final_genes = coalesce(manual_genes, cleaned_genes)) %>%
  
  mutate(final_genes = str_replace_all(final_genes, "-", ",")) %>%
  separate_rows(final_genes, sep = ",") %>%
  mutate(final_genes = str_trim(final_genes)) %>%
  
  mutate(Lacz_AERorECT = coalesce(Lacz_AERorECT, "Untested/Other"))

# Extract Base Peak Categories from great bpe dataframes
peaks_com_genes <- unique(df_combined_com$gene)
peaks_sp8_genes <- unique(df_combined_sp8$gene)
peaks_sp6_genes <- unique(df_combined_sp6$gene)

# ==============================================================================
# 3. Integrated Peak Processing and Classification Pipeline (Including Non-DEGs)
# ==============================================================================
annotated_hits_updated <- processed_hits %>%
  distinct(id, .keep_all = TRUE) %>%

  left_join(gs_lookup, by = c("final_genes" = "gene_symbol")) %>%

  mutate(
    source = case_when(
      source == "Common" ~ "Shared",
      source == "Sp8"    ~ "Sp8-specific",
      source == "Sp6"    ~ "Sp6-specific",
      TRUE               ~ source
    ),
    category = case_when(
      raw_category == "Common"     ~ "shared_sp8&sp6dependent",
      raw_category == "Common_sp8" ~ "shared_sp8dependent",
      raw_category == "redun"      ~ "redundant",
      raw_category == "Common_sp6" ~ "shared_sp6dependent",
      is.na(raw_category)          ~ "non_DEG",
      TRUE                         ~ raw_category
    )
  ) %>%
  select(-raw_category) %>%
  relocate(final_genes, .after = source) %>%
  relocate(category,    .after = source)

writexl::write_xlsx(annotated_hits_updated, "revisions/output/Source_Data_8_VISTA_hits.xlsx")
