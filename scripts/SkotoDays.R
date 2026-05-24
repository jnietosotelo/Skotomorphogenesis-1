############# Analysis and figures for Evolution of skotomorphogenic development in seedlings of Zea and Tripsacum


##Clean session 

rm(list = ls())
graphics.off()
cat("\014")


## Load libraries

library(ggplot2)
library(dplyr)
library(agricolae)
library(patchwork)
library(cowplot)
library(tidyr)


# All results generated in this script will be saved in the "results" directory

dir.create("results", showWarnings = FALSE)


# This script analyzes phenotypic data from an 8-day skotomorphogenic experiment in 
# Zea and Tripsacum seedlings, where each row corresponds to an individual.
# All tables and figures are derived from EvolEscoDB.csv dataset


DF <- read.csv("data/SDDBD.csv")


# Traits analyses

vars <- c(
  "Stem","Mesocotyl","Coleoptile",
  "FPL","FPLCol","PR",
  "ARCN","SR","BreakCol","RDSD", "SRV", "SRD"
)

#Data subsets


DF_maize <- DF %>% filter(Group == "Maize")
DF_teo   <- DF %>% filter(Group == "Teosinte")
DF_trip  <- DF %>% filter(Group == "Tripsacum")

#==============================================================================

#Table S10

# Function to compute summary statistics (Min, Max, Mean, SD, and CV) for Maize, Teosinte, and Tripsacum

calc_summary <- function(data, group_name){
  
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
    mutate(Group = group_name)
}

# Combine results into a single table

table_Sko<- bind_rows(
  calc_summary(DF_maize, "Maize"),
  calc_summary(DF_teo, "Teosinte"),
  calc_summary(DF_trip, "Tripsacum")
) %>%
  dplyr::select(Group, Variable, Min, Max, Mean, SD, CV)


# Save table

dir.create("results", showWarnings = FALSE)

write.csv(table_Sko,
          "results/Table_SkotoDays_Summary_.csv",
          row.names = FALSE)



# Function to perform two-way ANOVA and Tukey's test for all traits

run_anova_tukey <- function(data, var){
  
  form <- as.formula(paste(var, "~ Group"))
  model <- aov(form, data = data)
  
  tuk <- HSD.test(model, "Group", group = TRUE)
  
  groups <- tuk$groups
  groups$Group <- rownames(groups)
  
  colnames(groups)[1] <- "Mean"
  colnames(groups)[2] <- "Letter"
  
  sm <- summary(model)[[1]]
  
  data.frame(
    Variable = var,
    F_value = sm["Group", "F value"],
    p_value = sm["Group", "Pr(>F)"],
    Group = groups$Group,
    Mean = groups$Mean,
    Letter = groups$Letter
  )
}

# Perform one-way ANOVA for all traits

anova_results <- lapply(vars, function(v){
  
  if(!v %in% colnames(DF)){
    return(data.frame(
      Variable = v,
      F_value = NA,
      p_value = NA,
      Group = NA,
      Mean = NA,
      Letter = NA
    ))
  }
  
  run_anova_tukey(DF, v)
})

anova_table <- do.call(rbind, anova_results)

anova_table$Variable <- factor(anova_table$Variable, levels = vars)

anova_table <- anova_table %>%
  arrange(Variable, Group)

#Save table

write.csv(anova_table,
          "results/Table_Sko_ANOVA_Tukey.csv",
          row.names = FALSE)

#===============================================================================

# Violin plots showing trait variation across groups

# Custom titles

custom_titles <- c(
  Stem = "Stem",
  Mesocotyl = "Mesocotyl",
  Coleoptile = "Coleoptile",
  FPL = "FPL",
  FPLCol = "FPL-Col",
  PR = "PR",
  SR = "SR",
  SRV = "SRV",
  SRD = "SRD",
  RDSD = "RDSD",
  ARCN = "ARCN"
  
)

# Function to generate violin plots

plot_violin <- function(data, variable,
                        y_label = "Length (cm)",
                        custom_title = NULL) {
  
  # ANOVA
  formula <- as.formula(paste(variable, "~ Group"))
  model <- aov(formula, data = data)
  
  # Tukey letters
  res <- HSD.test(model, "Group", group = TRUE)
  letras <- res$groups
  letras$Group <- rownames(letras)
  colnames(letras)[1] <- "media"
  colnames(letras)[2] <- "letra"
  
  # Position for letters
  pos <- data %>%
    group_by(Group) %>%
    summarise(
      ymax = max(.data[[variable]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(letras, by = "Group")
  
  # Plot
  p <- ggplot(data, aes(
    x = factor(Group, levels = c("Maize","Teosinte","Tripsacum")),
    y = .data[[variable]],
    color = Group
  )) +
    
    geom_violin(trim = TRUE,size = 1.5, scale = "count" ) +
    
    geom_boxplot(  width = 0.05,size = 1,   outlier.size = 0.5) +
    
    stat_summary(  fun = mean,geom = "point", shape = 16, size = 3,color = "black") +
    
    geom_text(data = pos,  aes(x = Group, y = ymax + 7, label = letra),
              inherit.aes = FALSE,size = 8,  color = "black",   fontface = "bold"
    ) +
    
    scale_color_manual(values = c(
      "Maize" = "orange",
      "Teosinte" = "darkgreen",
      "Tripsacum" = "blue"
    )) +
    
    labs(  title = custom_title,x = "",  y = y_label) +
    
    theme_minimal() +
    
    theme(plot.title = element_text(hjust = 0.5,size = 18,face = "bold"  ),
          
          axis.text.x = element_text( size = 18,face = "bold"  ),
          
          axis.title.y = element_text(size = 18,face = "bold"  ),
          
          axis.text.y = element_text( size = 18,  face = "bold"  ),
          
          legend.position = "none",  panel.border = element_rect(  color = "black",
                                                                   fill = NA,  size = 1.5)
    )+
    
    scale_x_discrete(labels = c(
      "Maize" = expression(bold("Maize, 7d")),
      "Teosinte" = expression(bolditalic("Teosinte, 7d")),
      "Tripsacum" = expression(bolditalic("Tripsacum, 17d"))
    ))
  
  return(p)
}

# Variables for plotting

vars <- c(
  "Stem",
  "Mesocotyl",
  "Coleoptile",
  "FPL",
  "FPLCol",
  "PR",
  "RDSD"
)

# Check variables

vars <- vars[vars %in% colnames(DF)]

print(vars)

# Generate plots

for (v in vars) {
  
  message("Processing: ", v)
  
  tryCatch({
    
    p <- plot_violin(
      DF,
      variable = v,
      custom_title = custom_titles[v]
    )

    
  }, error = function(e) {
    
    message("Error in: ", v)
    message(e)
  })
}

# Bubble plot function

plot_bubble <- function(data, variable,
                        y_label = "Number per seedling",
                        custom_title = NULL) {
  
  # ANOVA
  formula <- as.formula(paste(variable, "~ Group"))
  model <- aov(formula, data = data)
  
  # Tukey letters
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
  p <- ggplot(
    data,
    aes(
      x = Group,
      y = .data[[variable]],
      color = Group
    )
  ) +
    
    geom_count(alpha = 0.7) +
    
    scale_size_area(max_size = 8) +
    scale_y_continuous(
      breaks = function(x) seq(
        floor(min(x, na.rm = TRUE)),
        ceiling(max(x, na.rm = TRUE)),
        by = 1
      )
    )+
    
    geom_text(
      data = pos,
      aes(x = Group, y = ymax + 1, label = letra),
      inherit.aes = FALSE,
      size = 8,
      color = "black",
      fontface = "bold"
    ) +
    
    scale_color_manual(values = c(
      "Maize" = "orange",
      "Teosinte" = "darkgreen",
      "Tripsacum" = "blue"
    )) +
    
    labs(
      title = custom_title,
      x = "",
      y = y_label
    ) +
    
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18),
      text = element_text(face = "bold"),
      axis.text.x = element_text(size = 18, face = "bold"),
      axis.title.x  = element_text(size = 18, face = "bold"),
      axis.title.y = element_text(size = 18, face = "bold"),
      axis.text.y = element_text(size = 18, face = "bold"),
      legend.position = "none",
      panel.border = element_rect(
        color = "black",
        fill = NA,
        size = 1.5
      )
    ) +
    
    scale_x_discrete(labels = c(
      "Maize" = expression(bold("Maize, 7d")),
      "Teosinte" = expression(bolditalic("Teosinte, 7d")),
      "Tripsacum" = expression(bolditalic("Tripsacum, 17d"))
    ))
  
  return(p)
}

# Bubble variables

vars_bubble <- c("ARCN","SRV", "SRD","SR")

vars_bubble <- vars_bubble[vars_bubble %in% colnames(DF)]

print(vars_bubble)

# Generate bubble plots

for (v in vars_bubble) {
  
  message("Processing: ", v)
  
  tryCatch({
    
    p <- plot_bubble(
      DF,
      variable = v,
      custom_title = custom_titles[v]
    )

    
  }, error = function(e) {
    
    message("Error in bubble: ", v)
    message(e)
  })
}

# Combine plots

plots <- list()

# Preserve order
for (v in vars) {
  
  plots[[length(plots) + 1]] <- plot_violin(
    DF,
    variable = v,
    custom_title = custom_titles[v]
  )
}

for (v in vars_bubble) {
  
  plots[[length(plots) + 1]] <- plot_bubble(
    DF,
    variable = v,
    custom_title = custom_titles[v]
  )
}

# Add panel labels

letras <- LETTERS[1:length(plots)]

plots <- lapply(seq_along(plots), function(i) {
  
  plots[[i]] +
    
    labs(tag = letras[i]) +
    
    theme(
      plot.tag = element_text(size = 20, face = "bold"),
      plot.tag.position = c(0, 1)
    )
})

# Final figure

fig_final <- wrap_plots(plots, ncol = 3)

# Save figure

ggsave(
  "results/Figure_S15.tiff",
  fig_final,
  width = 20,
  height = 14,
  dpi = 300,
  compression = "lzw"
)


#===============================================================================

