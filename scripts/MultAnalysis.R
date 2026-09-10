############# Analysis and figures for Evolution of skotomorphogenic development in seedlings of Zea and Tripsacum


##Clean session 

rm(list = ls())
graphics.off()
cat("\014")


## Load libraries

library(ggplot2)
library(dplyr)
library(ggfortify)
library(MASS)
library(patchwork)
library(tidyr)
library(ape)
library(agricolae)

# All results generated in this script will be saved in the "results" directory

dir.create("results", showWarnings = FALSE)


##Load Data

DF <- read.csv("data/EvolEscoDB.csv")



DF$Group <- factor(DF$Group, levels = c("Maize", "Teosinte", "Tripsacum"))

corOrgan <- DF %>%
  dplyr::select(Group, Stem, Mesocotyl, Coleoptile, FPL, FPLCol, PR, ARCN, SR, RDSD) %>%
  tidyr::drop_na()

#Name Correction

colnames(corOrgan)[colnames(corOrgan) == "FPLCol"] <- "FPL-Col"

# PCA analysis

corOrgansp <- corOrgan %>%
  dplyr::select(-Group)

respca <- prcomp(corOrgansp, scale. = TRUE, center = TRUE)

# Explain variance

var_exp_pca <- summary(respca)$importance[2, 1:2] * 100

# LOADINGS (variables)

loadings <- as.data.frame(respca$rotation)

# Save loading's

write.csv(
  loadings,
  "results/Table_S13_PCA_loadings.csv",
  row.names = TRUE
)

print(loadings)

# Save variance explain 

var_exp <- as.data.frame(summary(respca)$importance)
var_exp <- t(var_exp)
var_exp <- as.data.frame(var_exp)

write.csv(
  var_exp,
  "results/Table_S13_PCA_variance_explained.csv",
  row.names = TRUE
)

print(var_exp)

# PCA plot

g1 <- autoplot(
  respca,
  data = corOrgan,
  colour = "Group",
  loadings = TRUE,
  loadings.label = TRUE,
  loadings.label.repel = TRUE,
  
  loadings.colour = "#6A3D9A",
  loadings.label.colour = "#6A3D9A",
  
  loadings.size = 2,
  loadings.label.size = 5,
  x = 1,
  y = 2
) +
  
  scale_color_manual(values = c(
    "Maize" = "orange",
    "Teosinte" = "darkgreen",
    "Tripsacum" = "blue"
  )) +
  
  labs(
    x = paste0("PC1 (", round(var_exp_pca[1], 1), "%)"),
    y = paste0("PC2 (", round(var_exp_pca[2], 1), "%)")
  ) +
  
  theme(
    panel.background = element_rect(fill = "white"),
    
    panel.grid.major = element_line(
      color = "gray",
      linewidth = 0.5
    ),
    
    panel.grid = element_blank(),
    
    legend.title = element_blank(),
    legend.position = "top",
    
    axis.text = element_text(
      size = 16,
      face = "bold"
    ),
    
    axis.title = element_text(
      size = 18,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 14,
      face = "bold"
    ),
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1.5
    )
  )

# Save PCA plot

ggsave(
  "results/Figure_6.tiff",
  plot = g1,
  width = 7,
  height = 6,
  dpi = 300,
  compression = "lzw"
)

# LDA Analysis 


Data <- DF %>%
  dplyr::select(Group, Stem, Mesocotyl, Coleoptile, FPL, FPLCol, PR, ARCN, SR, RDSD) %>%
  tidyr::drop_na()

modelo_lda <- lda(Group ~ ., data = Data)

# Predictions

pred <- predict(modelo_lda)

lda_df <- data.frame(
  LD1 = pred$x[,1],
  LD2 = pred$x[,2],
  Group = Data$Group
)

# Eplain variace LDA

var_exp_lda <- modelo_lda$svd^2 / sum(modelo_lda$svd^2)
var_exp_lda <- round(var_exp_lda * 100, 1)

# Coefficients LDA

coef_df <- as.data.frame(modelo_lda$scaling)
coef_df$Variable <- rownames(coef_df)

# Sort

coef_df <- coef_df %>%
  dplyr::select(Variable, everything())

# Save coefficients
write.csv(
  coef_df,
  "results/Table_S14_LDA_coefficients.csv",
  row.names = FALSE
)

print(coef_df)

# Scatter plot LDA

g7a <- ggplot(
  lda_df,
  aes(x = LD1, y = LD2, color = Group)
) +
  
  geom_point(size = 3, alpha = 0.8) +
  
  scale_color_manual(
    values = c(
      "Maize" = "orange",
      "Teosinte" = "darkgreen",
      "Tripsacum" = "blue"
    ),
    
    labels = c(
      "Maize",
      expression(bolditalic("Teosinte")),
      expression(bolditalic("Tripsacum"))
    )
  )+
  
  labs(
    x = paste0("LD1 (", round(var_exp_lda[1], 1), "%)"),
    y = paste0("LD2 (", round(var_exp_lda[2], 1), "%)")
  ) +
  
  theme_minimal() +
  
  theme(
    
    axis.text = element_text(
      size = 16,
      face = "bold"
    ),
    
    axis.title = element_text(
      size = 18,
      face = "bold"
    ),
    
    legend.title = element_blank(),
    
    legend.text = element_text(
      size = 14,
      face = "bold"
    ),
    
    legend.position = "top",
    
    panel.grid.major = element_line(
      color = "gray80",
      linewidth = 0.5
    ),
    
    panel.grid.minor = element_line(
      color = "gray90",
      linewidth = 0.25
    ),
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1.5
    )
  )


# Guardar
ggsave(
  "results/Figure_7.tiff",
  plot = g7a,
  width = 7,
  height = 6,
  dpi = 300,
  compression = "lzw"
)


#========================================
#Figure_S14


# LD1

g6b <- coef_df %>%
  
  mutate(
    Variable = ifelse(
      Variable == "FPLCol",
      "FPL-Col",
      Variable
    )
  ) %>%
  
  arrange(desc(abs(LD1))) %>%
  
  mutate(
    Variable = factor(
      Variable,
      levels = Variable
    )
  ) %>%
  
  ggplot(aes(x = Variable, y = LD1)) +
  
  geom_col(fill = "black") +
  
  labs(
    y = "Weighting on LD1",
    x = ""
  ) +
  
  theme_minimal() +
  
  theme(
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 12,
      face = "bold"
    ),
    
    axis.text.y = element_text(
      size = 12,
      face = "bold"
    ),
    
    axis.title.y = element_text(
      size = 14,
      face = "bold"
    ),
    
    panel.grid.major = element_line(
      color = "gray85",
      linewidth = 0.5
    ),
    
    panel.grid.minor = element_blank(),
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1.5
    )
  )


# LD2

g6c <- coef_df %>%
  
  mutate(
    Variable = ifelse(
      Variable == "FPLCol",
      "FPL-Col",
      Variable
    )
  ) %>%
  
  arrange(desc(abs(LD2))) %>%
  
  mutate(
    Variable = factor(
      Variable,
      levels = Variable
    )
  ) %>%
  
  ggplot(aes(x = Variable, y = LD2)) +
  
  geom_col(fill = "black") +
  
  labs(
    y = "Weighting on LD2",
    x = ""
  ) +
  
  theme_minimal() +
  
  theme(
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 12,
      face = "bold"
    ),
    
    axis.text.y = element_text(
      size = 12,
      face = "bold"
    ),
    
    axis.title.y = element_text(
      size = 14,
      face = "bold"
    ),
    
    panel.grid.major = element_line(
      color = "gray85",
      linewidth = 0.5
    ),
    
    panel.grid.minor = element_blank(),
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1.5
    )
  )


# Joint plot

final_plot <- (g6b + g6c) +
  
  plot_layout(ncol = 2) +
  
  plot_annotation(tag_levels = "A") &
  
  theme(
    
    plot.tag = element_text(
      size = 18,
      face = "bold"
    ),
    
    plot.tag.position = c(0, 1)
  )

# Show

print(final_plot)

# Save

ggsave(
  "results/Figure_S14.tiff",
  plot = final_plot,
  width = 12,
  height = 5,
  dpi = 300,
  compression = "lzw"
)

#=================================================================================  
#Hierachical Cluster analysis 


# All results generated in this script will be saved in the "results" directory

dir.create("results", showWarnings = FALSE)

#Load data 

DF <- read.csv("data/EvolEscoMean.csv")

# Replace NA

DF[is.na(DF)] <- 0

# The "Code" column is used as row names

rownames(DF) <- DF$Code

# Variables

data <- DF %>%
  dplyr::select(
    Stem, Mesocotyl, Coleoptile,
    FPL, FPLCol, PR,
    ARCN, SR, BreakCol, RDSD
  )

# Scale

data_scaled <- scale(data)

# Distance matrix

distancias <- dist(data_scaled)

# Clustering

hc <- hclust(distancias, method = "ward.D2")
phy <- as.phylo(hc)

# Calculation of the cophenetic coefficient

dist_cof <- cophenetic(hc)
cofenetico <- cor(distancias, dist_cof)

print(paste("Coeficiente cofenético:", round(cofenetico, 3)))

# CLUSTERS

k <- 4
clusters <- cutree(hc, k = k)

# Cluster table

tabla_clusters <- data.frame(
  Code = rownames(DF),
  Cluster = clusters
)

# Cluster assignation

tabla_clusters$Cluster <- factor(
  tabla_clusters$Cluster,
  levels = 1:4,
  labels = c("I", "II", "IV", "III")
)

#Save cluster table

dir.create("results", showWarnings = FALSE)

write.csv(
  tabla_clusters,
  "results/Table_Clusters_Assignation.csv",
  row.names = FALSE
)

#Assign cluster Colors


colores <- c("#CDCD00", "#D55E00",  "#6A3D9A", "#009E73")
tip_colors <- colores[as.numeric(tabla_clusters$Cluster)]

# Save tree diagram

tiff("results/Figure_8.tiff",
     width = 7, 
     height = 7, 
     units = "in",
     res = 300, 
     compression = "lzw")

plot(phy,
     type = "fan",
     tip.color = tip_colors,
     cex = 0.6,
     no.margin = TRUE,
     label.offset = 0.5)

# tree cordinates
pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)

# Cluster tips assigment
cluster_tip <- tabla_clusters$Cluster[
  match(phy$tip.label, tabla_clusters$Code)
]
x <- pp$xx[1:length(phy$tip.label)]
y <- pp$yy[1:length(phy$tip.label)]

# Tips angle
angulos <- atan2(y, x)

orden <- order(angulos)

angulos_ord <- angulos[orden]
cluster_ord <- cluster_tip[orden]
cambios <- which(
  cluster_ord[-length(cluster_ord)] !=
    cluster_ord[-1]
)

#Radius max

radio_max <- max(sqrt(x^2 + y^2))

# Dotted line

for(i in cambios){
  
  theta1 <- angulos_ord[i]
  theta2 <- angulos_ord[i + 1]
  
  if(abs(theta2 - theta1) > pi){
    
    if(theta1 < theta2){
      theta1 <- theta1 + 2*pi
    } else {
      theta2 <- theta2 + 2*pi
    }
  }
  
  theta <- (theta1 + theta2) / 2
  
  # Dotted line  center to tips
  radio_inicio <- radio_max * 0.4
  radio_final  <- radio_max * 1.15
  
  segments(
    radio_inicio * cos(theta),
    radio_inicio * sin(theta),
    radio_final * cos(theta),
    radio_final * sin(theta),
    lty = 2,
    lwd = 1,
    col = "black"
  )
}

# Tags

text(-4.8,  3.5,  "III",   cex = 1.5, font = 2)
text( 4.8,  3.5,  "IV",  cex = 1.5, font = 2)
text( 3.8, -4.5,  "II", cex = 1.5, font = 2)
text(-2.8, -4.5,  "I",  cex = 1.5, font = 2)

dev.off()


# TABLA S15 Accesions assigned to  Cluster

# Variables en orden fijo

##Clean session 

rm(list = ls())
graphics.off()
cat("\014")


# This script analyzes phenotypic data from an 8-day skotomorphogenic experiment in 
# Zea and Tripsacum seedlings, where each row corresponds to an individual assigned to a 
# cluster como se obtuvo en el analsisi pasado.
# All tables and figures are derived from EvolEscoDBHCJ.csv dataset


DF <- read.csv("data/EvolEscoDBHCJ.csv")


vars <- c(
  "Stem","Mesocotyl","Coleoptile",
  "FPL","FPLCol","PR",
  "ARCN","SR","BreakCol","RDSD"
)

#Cluster labels


DF$Cluster <- factor(DF$Cluster, levels = c("I","II","III","IV"))

# Function to compute summary statistics for Table S15 (by cluster)


calc_summary_cluster <- function(data, cluster_name){
  
  data %>%
    dplyr::select(all_of(vars)) %>%
    pivot_longer(cols = everything(),
                 names_to = "Variable",
                 values_to = "Value") %>%
    
    mutate(Variable = factor(Variable, levels = vars)) %>%
    
    group_by(Variable) %>%
    summarise(
      Min = min(Value, na.rm = TRUE),
      Max = max(Value, na.rm = TRUE),
      Mean = mean(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE),
      CV = ifelse(Mean == 0, NA, (SD / Mean) * 100),
      .groups = "drop"
    ) %>%
    
    arrange(Variable) %>%
    mutate(Cluster = cluster_name)
}


# Subset selection


clusters <- levels(DF$Cluster)

tabla_S15 <- lapply(clusters, function(cl){
  
  df_sub <- DF %>% filter(Cluster == cl)
  
  calc_summary_cluster(df_sub, cl)
})

tabla_S15 <- bind_rows(tabla_S15) %>%
  dplyr::select(Cluster, Variable, Min, Max, Mean, SD, CV)

# Save Table S15

write.csv(tabla_S15,
          "results/Table_S15_Summary_Cluster.csv",
          row.names = FALSE)

# ANOVA + TUKEY by (Cluster)


run_anova_tukey_cluster <- function(data, var){
  
  form <- as.formula(paste(var, "~ Cluster"))
  model <- aov(form, data = data)
  
  tuk <- HSD.test(model, "Cluster", group = TRUE)
  
  groups <- tuk$groups
  groups$Cluster <- rownames(groups)
  
  colnames(groups)[1] <- "Mean"
  colnames(groups)[2] <- "Letter"
  
  sm <- summary(model)[[1]]
  
  data.frame(
    Variable = var,
    F_value = sm["Cluster", "F value"],
    p_value = sm["Cluster", "Pr(>F)"],
    Cluster = groups$Cluster,
    Mean = groups$Mean,
    Letter = groups$Letter
  )
}

anova_results_cluster <- lapply(vars, function(v){
  
  if(!v %in% colnames(DF)){
    return(data.frame(
      Variable = v,
      F_value = NA,
      p_value = NA,
      Cluster = NA,
      Mean = NA,
      Letter = NA
    ))
  }
  
  run_anova_tukey_cluster(DF, v)
})

anova_table_cluster <- do.call(rbind, anova_results_cluster)

#### order

anova_table_cluster$Variable <- factor(anova_table_cluster$Variable, levels = vars)

anova_table_cluster <- anova_table_cluster %>%
  arrange(Variable, Cluster)

###Save Table S15 anova

write.csv(anova_table_cluster,
          "results/Table_S15_ANOVA_Tukey.csv",
          row.names = FALSE)



#===============================================================================

#####  Figure 11 Violin plot and dot plot by clusters


DF_cluster <- DF %>%
  filter(!is.na(Cluster))

DF_cluster$Cluster <- factor(DF_cluster$Cluster,
                             levels = c("I", "II", "III", "IV"))

DF_cluster$Group <- DF_cluster$Cluster

#### Function violin by cluster


# Custom titles

custom_titles <- c(
  Stem = "Stem",
  Mesocotyl = "Mesocotyl",
  Coleoptile = "Coleoptile",
  FPL = "FPL",
  FPLCol = "FPL-Col",
  PR = "PR",
  SR = "SR",
  RDSD = "RDSD",
  ARCN = "ARCN",
  SR = "SR"
)



# Function violin by cluster

plot_violin_cluster <- function(data,
                                variable,
                                y_label = "Length (cm)",
                                custom_title = NULL) {
  
  # ANOVA
  model <- aov(as.formula(paste(variable, "~ Group")), data = data)
  
  # Tukey groups
  res <- agricolae::HSD.test(model, "Group", group = TRUE)
  
  letras <- res$groups
  letras$Group <- rownames(letras)
  colnames(letras)[2] <- "letra"
  
  # Position letters
  pos <- data %>%
    group_by(Group) %>%
    summarise(
      ymax = max(.data[[variable]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(letras, by = "Group")
  
  # Plot
  ggplot(data,
         aes(
           x = Group,
           y = .data[[variable]],
           color = Group
         )) +
    
    geom_violin(
      trim = TRUE,
      linewidth = 1.5,
      scale = "count"
    ) +
    
    geom_boxplot(
      width = 0.05,
      linewidth = 1,
      outlier.size = 0.5
    ) +
    
    stat_summary(fun = mean,
                 geom = "point",
                 shape = 16,
                 size = 1,
                 color = "black") +
    
    geom_text(
      data = pos,
      aes(x = Group,
          y = ymax + 10,
          label = letra),
      inherit.aes = FALSE,
      size = 7,
      fontface = "bold"
    ) +
    
    scale_color_manual(values = c(
      "I" = "#CDCD00",
      "II" = "#D55E00",
      "III" = "#009E73",
      "IV" = "#6A3D9A"
    )) +
    
    labs(
      title = custom_title,
      x = "",
      y = y_label
    ) +
    
    theme_minimal() +  
    theme(
      
      plot.title = element_text(
        hjust = 0.5,
        size = 18,
        face = "bold"
      ),
      
      axis.text.x = element_text(
        size = 18,
        face = "bold"
      ),
      
      axis.title.y = element_text(
        size = 18,
        face = "bold"
      ),
      
      axis.text.y = element_text(
        size = 18,
        face = "bold"
      ),
      
      legend.position = "none",
      
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 1.5
      )
    )
}

# Function bubble by cluster

plot_bubble_cluster <- function(data,
                                variable,
                                y_label = "Number per seedling",
                                custom_title = NULL) {
  
  # ANOVA
  model <- aov(as.formula(paste(variable, "~ Group")), data = data)
  
  # Tukey groups
  res <- agricolae::HSD.test(model, "Group", group = TRUE)
  
  letras <- res$groups
  letras$Group <- rownames(letras)
  colnames(letras)[2] <- "letra"
  
  # Position letters
  pos <- data %>%
    group_by(Group) %>%
    summarise(
      ymax = max(.data[[variable]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(letras, by = "Group")
  
  # Plot
  ggplot(data,
         aes(
           x = Group,
           y = .data[[variable]],
           color = Group
         )) +
    
    geom_count(alpha = 0.7) +
    
    scale_size_area(max_size = 8) +
    scale_y_continuous(
      breaks = function(x)
        seq(floor(min(x)), ceiling(max(x)), by = 1)
    ) +
    
    geom_text(
      data = pos,
      aes(x = Group,
          y = ymax + 1.5,
          label = letra),
      inherit.aes = FALSE,
      size = 8,
      color = "black",
      fontface = "bold"
    ) +
        scale_color_manual(values = c(
      "I" = "#0072B2",
      "II" = "#D55E00",
      "III" = "#009E73",
      "IV" = "#6A3D9A"
    )) +
    
    labs(
      title = custom_title,
      x = "",
      y = y_label
    ) +
    
    theme_minimal() +
    
    theme(
      
      plot.title = element_text(
        hjust = 0.5,
        size = 18,
        face = "bold"
      ),
      
      axis.text.x = element_text(
        size = 18,
        face = "bold"
      ),
      
      axis.title.y = element_text(
        size = 18,
        face = "bold"
      ),
      
      axis.text.y = element_text(
        size = 18,
        face = "bold"
      ),
      
      legend.position = "none",
      
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 1.5
      )
    )
}

# Variables

vars <- c(
  "Stem",
  "Mesocotyl",
  "Coleoptile",
  "FPL",
  "FPLCol",
  "PR",
  "RDSD"
)

vars <- vars[vars %in% colnames(DF_cluster)]

vars_bubble <- c("ARCN", "SR")

vars_bubble <- vars_bubble[vars_bubble %in% colnames(DF_cluster)]

# Generate plots

plots <- list()

for (v in vars) {
  
  plots[[length(plots) + 1]] <- plot_violin_cluster(
    DF_cluster,
    variable = v,
    custom_title = custom_titles[v]
  )
}

for (v in vars_bubble) {
  
  plots[[length(plots) + 1]] <- plot_bubble_cluster(
    DF_cluster,
    variable = v,
    custom_title = custom_titles[v]
  )
}

# Add labels A-I

letras <- LETTERS[1:length(plots)]

plots <- lapply(seq_along(plots), function(i) {
  
  plots[[i]] +
    
    labs(tag = letras[i]) +
    
    theme(
      plot.tag = element_text(
        size = 16,
        face = "bold"
      ),
      plot.tag.position = c(0.02, 0.98)
    )
})

# Final figure

fig_final <- wrap_plots(plots, ncol = 3)
# Save figure


ggsave(
  "results/Figure_9.tiff",
  fig_final,
  width = 14,
  height = 12,
  dpi = 300,
  compression = "lzw"
)










