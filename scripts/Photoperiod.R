
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
library(tidyr)
library(tibble)

# All results generated in this script will be saved in the "results" directory

dir.create("results", showWarnings = FALSE)


# Data

DF <- read.csv("data/Photoperiod.csv", stringsAsFactors = TRUE)

# Traits selected

vars <- c(
  "Stem","Mesocotyl","Coleoptile",
  "FPL","FPLCol","PR",
  "ARCN","SR","BreakCol","RDSD"
)

#Sort of experiment

DF$Experiment <- factor(DF$Experiment,
                        levels = c("8D","8L","15D","15L"))

#Subsets by accession


DF_maize    <- DF %>% filter(Code == "M52")
DF_teocinte <- DF %>% filter(Code == "Teo2")

#Function Summary data


calc_summary_exp <- function(data, code_name){
  
  data %>%
    dplyr::select(Experiment, all_of(vars)) %>%
    
    pivot_longer(
      cols = all_of(vars),
      names_to = "Variable",
      values_to = "Value"
    ) %>%
    
    mutate(Variable = factor(Variable, levels = vars)) %>%
    
    group_by(Experiment, Variable) %>%
    summarise(
      Min = min(Value, na.rm = TRUE),
      Max = max(Value, na.rm = TRUE),
      Mean = mean(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE),
      CV = ifelse(Mean == 0, NA, (SD / Mean) * 100),
      .groups = "drop"
    ) %>%
    
    arrange(Variable) %>%
    
    mutate(Code = code_name)
}

#Function to maize and teocinte 

tab_maize <- calc_summary_exp(DF_maize, "M52")
tab_teo   <- calc_summary_exp(DF_teocinte, "Teo2")

##Joint table


tabla_S9 <- bind_rows(tab_maize, tab_teo) %>%
  dplyr::select(Code, Experiment, Variable, Min, Max, Mean, SD, CV)

#Save table 

dir.create("results", showWarnings = FALSE)

write.csv(
  tabla_S9,
  "results/Table_S8_Summary_Phototoperiod.csv",
  row.names = FALSE
)

# Function ANOVA + TUKEY (for subset)


run_anova_tukey_exp <- function(data, var, code_name){
  
  form <- as.formula(paste(var, "~ Experiment"))
  model <- aov(form, data = data)
  
  tuk <- HSD.test(model, "Experiment", group = TRUE)
  
  groups <- tuk$groups
  groups$Experiment <- rownames(groups)
  
  colnames(groups)[1] <- "Mean"
  colnames(groups)[2] <- "Letter"
  
  sm <- summary(model)[[1]]
  
  data.frame(
    Code = code_name,
    Variable = var,
    F_value = sm["Experiment", "F value"],
    p_value = sm["Experiment", "Pr(>F)"],
    Experiment = groups$Experiment,
    Mean = groups$Mean,
    Letter = groups$Letter
  )
}

#Anova analysis

anova_maize <- lapply(vars, function(v){
  
  if(!v %in% colnames(DF_maize)){
    return(data.frame(
      Code = "M52",
      Variable = v,
      F_value = NA,
      p_value = NA,
      Experiment = NA,
      Mean = NA,
      Letter = NA
    ))
  }
  
  run_anova_tukey_exp(DF_maize, v, "M52")
})

anova_teo <- lapply(vars, function(v){
  
  if(!v %in% colnames(DF_teocinte)){
    return(data.frame(
      Code = "Teo2",
      Variable = v,
      F_value = NA,
      p_value = NA,
      Experiment = NA,
      Mean = NA,
      Letter = NA
    ))
  }
  
  run_anova_tukey_exp(DF_teocinte, v, "Teo2")
})


anova_S9 <- bind_rows(
  do.call(rbind, anova_maize),
  do.call(rbind, anova_teo)
)

anova_S9$Variable <- factor(anova_S9$Variable, levels = vars)

anova_S9 <- anova_S9 %>%
  arrange(Code, Variable, Experiment)

#save anova result

write.csv(
  anova_S9,
  "results/Table_S8ANOVA_Tukey_Photoperiod.csv",
  row.names = FALSE
)
#============================================================================
####Figure 

# Subsets

DF_maize <- DF %>% filter(Group == "Maize")
DF_teocinte <- DF %>% filter(Group == "Teocinte")

# Function violín plot

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

plot_violin <- function(data, variable, y_label = "Length (cm)",
                        custom_title = NULL) {
  # ANOVA
  model <- aov(as.formula(paste(variable, "~ Experiment")), data = data)
  # Tukey letters 
  res <- HSD.test(model, "Experiment", group = TRUE)
  letras <- res$groups
  letras$Experiment <- rownames(letras)
  colnames(letras)[2] <- "letra"
  
  # Position for letters  
  
  pos <- data %>%
    group_by(Experiment) %>%
    summarise(ymax = max(.data[[variable]], na.rm = TRUE), .groups = "drop") %>%
    left_join(letras, by = "Experiment")
  
  # Plot
  
  ggplot(data, aes(
    x = factor(Experiment, levels = c("8D","8L","15D","15L")),
    y = .data[[variable]],
    color = Experiment
  )) +
    geom_violin(trim = TRUE,size = 1.5, scale = "count" ) +
    geom_boxplot(  width = 0.05,size = 1,   outlier.size = 0.5) +
    stat_summary(  fun = mean,geom = "point", shape = 16, size = 3,color = "black") +
    
    geom_text(data = pos,  aes(x = Experiment, y = ymax + 7, label = letra),
              inherit.aes = FALSE,size = 8,  color = "black",   fontface = "bold"
    ) +
    
    scale_color_manual(values = c(
      "8D" = "#0072B2", "15D" = "#0072B2",
      "8L" = "#D55E00",   "15L" = "#D55E00"
    )) +
    
    labs(  title = custom_title,x = "",  y = y_label) +
    
    theme_minimal() +
    
    theme(plot.title = element_text(hjust = 0.5,size = 18,face = "bold"  ),
          
          axis.text.x = element_text( size = 18,face = "bold"  ),
          
          axis.title.y = element_text(size = 18,face = "bold"  ),
          
          axis.text.y = element_text( size = 18,  face = "bold"  ),
          
          legend.position = "none",  panel.border = element_rect(  color = "black",
                                                                   fill = NA,  size = 1.5)
    )
}

# Función bubble (dot plot)


plot_bubble <- function(data, variable, y_label = "Number per seedling",
                        custom_title = NULL) {
  data$Experiment <- factor(data$Experiment,
                            levels = c("8D", "8L", "15D", "15L"))
  
  model <- aov(as.formula(paste(variable, "~ Experiment")), data = data)
  
  # Tukey letters  
  res <- HSD.test(model, "Experiment", group = TRUE)
  letras <- res$groups
  letras$Experiment <- rownames(letras)
  colnames(letras)[2] <- "letra"
  
  # Position letters
  
  pos <- data %>%
    group_by(Experiment) %>%
    summarise(ymax = max(.data[[variable]], na.rm = TRUE), .groups = "drop") %>%
    left_join(letras, by = "Experiment")
  
  ggplot(data, aes(x = Experiment, y = .data[[variable]], color = Experiment)) +
    geom_count(alpha = 0.7) +
    
    scale_size_area(max_size = 8) +
    scale_y_continuous(
      breaks = scales::pretty_breaks(n = 6)
    )+
    
    geom_text(
      data = pos,
      aes(x = Experiment, y = ymax + 1, label = letra),
      inherit.aes = FALSE,
      size = 8,
      color = "black",
      fontface = "bold"
    ) +
    
    scale_color_manual(values = c(
      "8D" = "#0072B2", "15D" = "#0072B2",
      "8L" = "#D55E00",   "15L" = "#D55E00"
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
      ))
}


# Variables

vars <- c("Stem","Mesocotyl","Coleoptile","FPL","FPLCol","PR","RDSD")
vars <- vars[vars %in% colnames(DF)]

vars_bubble <- c("ARCN","SR")
vars_bubble <- vars_bubble[vars_bubble %in% colnames(DF)]

# Function plots


make_figure <- function(data) {
  
  plots <- list()
  
  # Violin plots
  
  for (v in vars) {
    
    plots[[length(plots) + 1]] <- plot_violin(
      data,
      variable = v,
      custom_title = custom_titles[v]
    )
  }
  
  # Bubble plots
  
  for (v in vars_bubble) {
    
    plots[[length(plots) + 1]] <- plot_bubble(
      data,
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
          size = 20,
          face = "bold"
        ),
        
        plot.tag.position = c(0.02, 0.98)
      )
  })
  
  # Final figure
  
  wrap_plots(plots, ncol = 3)
}
# Figure 4 (Teocinte)

fig6 <- make_figure(DF_teocinte)

ggsave(
  "results/Figure_3.tiff",
  fig6,
  width = 14,
  height = 12,
  dpi = 300,
  compression = "lzw"
)

# Figure 5 (Maize)
fig5 <- make_figure(DF_maize)

ggsave(
  "results/Figure_4.tiff",
  fig5,
  width = 14,
  height = 12,
  dpi = 300,
  compression = "lzw"
)

