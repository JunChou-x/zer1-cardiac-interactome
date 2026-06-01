# ============================================================
# 03_visualization.R
#
# Visualization of prioritized candidates
# ============================================================

library(readxl)
library(dplyr)
library(ggplot2)
library(ggrepel)

# ------------------------------------------------------------
# Lollipop plot
# ------------------------------------------------------------

tier1_data <- read_excel(
  "results/Tier1_candidates.xlsx"
)

plot_data <- tier1_data %>%
  arrange(desc(log2FC)) %>%
  mutate(
    
    Rank=row_number(),
    
    Percentile=
      if(n()>1)
        100*(1-(Rank-1)/(n()-1))
    else
      100,
    
    Color_Group=
      case_when(
        Gene.Symbol=="Dvl2" ~ "Dvl2",
        Rank<=5 ~ "Control",
        TRUE ~ "Other"
      ),
    
    Is_Highlight=
      Gene.Symbol=="Dvl2" |
      Rank<=5
  )

n_tier1 <- nrow(plot_data)

p_lollipop <- ggplot(
  plot_data,
  aes(Rank,Percentile)
)+
  
  geom_segment(
    aes(
      xend=Rank,
      y=0,
      yend=Percentile,
      color=Color_Group
    )
  )+
  
  geom_point(
    aes(
      fill=Color_Group,
      size=Color_Group
    ),
    shape=21
  )+
  
  geom_label_repel(
    data=
      subset(
        plot_data,
        Is_Highlight
      ),
    aes(
      label=Gene.Symbol
    )
  )+
  
  theme_bw()

ggsave(
  "results/Lollipop_plot.pdf",
  p_lollipop,
  width=8,
  height=5
)

# ------------------------------------------------------------
# Volcano plot
# ------------------------------------------------------------

volcano_df <- read_excel(
  "data/Volcano_input.xlsx"
)

volcano_df <- volcano_df %>%
  mutate(
    
    negLogP=
      -log10(`p-value`+1e-10),
    
    Significance=
      case_when(
        
        log2FC>1 &
          `p-value`<0.05
        ~ "Up-regulated",
        
        log2FC< -1 &
          `p-value`<0.05
        ~ "Down-regulated",
        
        TRUE
        ~ "Not Significant"
      )
  )

p_volcano <- ggplot(
  
  volcano_df,
  
  aes(
    log2FC,
    negLogP,
    color=Significance
  )
  
)+
  
  geom_point() +
  
  geom_text_repel(
    
    data=
      subset(
        volcano_df,
        Significance!="Not Significant"
      ),
    
    aes(
      label=`Gene Symbol`
    )
  )

ggsave(
  "results/Volcano_plot.pdf",
  p_volcano,
  width=8,
  height=6
)

cat("Visualization completed.\n")