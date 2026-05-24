

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


DF <- read.csv("data/EvolEscoDB.csv")


# Traits analyses

vars <- c(
  "Stem","Mesocotyl","Coleoptile",
  "FPL","FPLCol","PR",
  "ARCN","SR","BreakCol","RDSD"
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

table_S10 <- bind_rows(
  calc_summary(DF_maize, "Maize"),
  calc_summary(DF_teo, "Teosinte"),
  calc_summary(DF_trip, "Tripsacum")
) %>%
  dplyr::select(Group, Variable, Min, Max, Mean, SD, CV)


# Save table

dir.create("results", showWarnings = FALSE)

write.csv(table_S10,
          "results/Table_S10_Summary_by_Group.csv",
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
          "results/Table_S10_ANOVA_Tukey.csv",
          row.names = FALSE)

#==============================================================================

# Tables S5, S6, and S9. Mean trait values per accession including all groups
# Coleoptile breakage was recorded as a binary trait 
# (0 = no breakage, 1 = breakage). Therefore, the percentage of coleoptile 
# breakage was calculated as the mean value multiplied by 100.

DF <- DF %>%
  mutate(BreakCol = BreakCol * 100)

# Ensure consistent variable order


vars <- vars[vars %in% colnames(DF)]

tabla_S6 <- DF %>%
  
  group_by(Code) %>%
  summarise(n = sum(complete.cases(across(all_of(vars)))), .groups = "drop") %>%
  
  left_join(
    DF %>%
      dplyr::select(Code, all_of(vars)) %>%
      pivot_longer(
        cols = all_of(vars),
        names_to = "Variable",
        values_to = "Value"
      ) %>%
      mutate(Variable = factor(Variable, levels = vars)) %>%
      group_by(Code, Variable) %>%
      summarise(Mean = mean(Value, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(
        names_from = Variable,
        values_from = Mean
      ),
    by = "Code"
  ) %>%
  
  dplyr::select (Code, n, everything())


# Save table containing the data for Tables S6,S8 and S9

dir.create("results", showWarnings = FALSE)

write.csv(
  tabla_S6,
  "results/Table_S679.csv",
  row.names = FALSE
)

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
      "Maize" = expression(bold("Maize")),
      "Teosinte" = expression(bolditalic("Teosinte")),
      "Tripsacum" = expression(bolditalic("Tripsacum"))
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
      breaks = scales::pretty_breaks(n = 6)
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
      "Maize" = expression(bold("Maize")),
      "Teosinte" = expression(bolditalic("Teosinte")),
      "Tripsacum" = expression(bolditalic("Tripsacum"))
    ))
  
  return(p)
}

# Bubble variables

vars_bubble <- c("ARCN", "SR")

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
  "results/Figure_1.tiff",
  fig_final,
  width = 14,
  height = 12,
  dpi = 300,
  compression = "lzw"
)


#===============================================================================


# Figure S5. Percentage of coleoptile tip breakage for maize and teosinte

## Figure S5A

mean_maize_code <- DF %>%
  filter(Group == "Maize") %>%
  group_by(Code) %>%
  summarise(
    mean = mean(BreakCol, na.rm = TRUE) * 1,
    sd = sd(BreakCol, na.rm = TRUE) * 1,
    n = n(),
    se = sd / sqrt(n),
    .groups = "drop"
  )


mean_maize_code <- mean_maize_code %>%
  arrange(desc(mean))

mean_maize_code$Code <- factor(mean_maize_code$Code,
                               levels = mean_maize_code$Code)

p5 <- ggplot(mean_maize_code,
             aes(x = Code, y = mean)) +
  
  geom_bar(
    stat = "identity",
    fill = "black",
    color = "white",
    linewidth = 0.3
  ) +
  
  coord_flip() +
  
  labs(
    title = "",
    x = "",
    y = "% Coleoptile tip breakage"
  ) +
  
  theme_minimal() +
  
  theme(
    
    plot.title = element_text(
      hjust = 0.5,
      size = 18,
      face = "bold"
    ),
    
    axis.text.x = element_text(
      size = 14,
      face = "bold"
    ),
    
    axis.text.y = element_text(
      size = 10,
      face = "bold"
    ),
    
    axis.title = element_text(
      size = 16,
      face = "bold"
    ),
    
    panel.grid.major.y = element_blank(),
    
    panel.grid.major.x = element_line(
      color = "gray85",
      linewidth = 0.5
    ),
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1.5
    )
  )


##S5B

mean_Teosinte_code <- DF %>%
  filter(Group == "Teosinte") %>%
  group_by(Code) %>%
  summarise(
    mean = mean(BreakCol, na.rm = TRUE) * 1,
    sd = sd(BreakCol, na.rm = TRUE) * 1,
    n = n(),
    se = sd / sqrt(n),
    .groups = "drop"
  )


mean_Teosinte_code <- mean_Teosinte_code %>%
  arrange(desc(mean))

mean_Teosinte_code$Code <- factor(mean_Teosinte_code$Code,
                               levels = mean_Teosinte_code$Code)

p6 <- ggplot(
  mean_Teosinte_code,
  aes(x = Code, y = mean)
) +
  
  geom_bar(
    stat = "identity",
    fill = "black",
    color = "white",
    linewidth = 0.3
  ) +
  
  coord_flip() +
  
  labs(
    title = "",
    x = "",
    y = "% Coleoptile tip breakage"
  ) +
  
  theme_minimal() +
  
  theme(
    
    plot.title = element_text(
      hjust = 0.5,
      size = 18,
      face = "bold"
    ),
    
    axis.text.x = element_text(
      size = 14,
      face = "bold"
    ),
    
    axis.text.y = element_text(
      size = 10,
      face = "bold"
    ),
    
    axis.title = element_text(
      size = 16,
      face = "bold"
    ),
    
    panel.grid.major.y = element_blank(),
    
    panel.grid.major.x = element_line(
      color = "gray85",
      linewidth = 0.5
    ),
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1.5
    )
  )

p5 <- p5 + theme(axis.title.y = element_blank())
p6 <- p6 + theme(axis.title.y = element_blank())


fig_S5 <- (p5 + p6) +
  
  plot_layout(ncol = 2) +
  
  plot_annotation(tag_levels = "A") &
  
  theme(
    plot.tag = element_text(
      size = 18,
      face = "bold"
    )
  )

ggsave(
  "results/Figure_S5.tiff",
  fig_S5,
  width = 9,
  height = 12,
  dpi = 300,
  compression = "lzw"
)

#===============================================================================

# Figure S6 and S7, Number of roots for seedling

#ARCN

DF_maize <- DF %>%
  filter(Group == "Maize")

freq_ARCN <- DF_maize %>%
  group_by(ARCN) %>%
  summarise(n = n(), .groups = "drop")


freq_ARCN <- freq_ARCN %>%
  mutate(
    percent = (n / sum(n)) * 100
  )

freq_ARCN$ARCN <- factor(freq_ARCN$ARCN, levels = sort(unique(freq_ARCN$ARCN)))


p1 <- ggplot(
  freq_ARCN,
  aes(x = ARCN, y = percent)
) +
  
  geom_bar(
    stat = "identity",
    fill = "black",
    color = "white",
    linewidth = 0.3
  ) +
  
  labs(
    title = "",
    x = "Number of ARCN per seedling",
    y = "Seedlings (%)"
  ) +
  
  theme_minimal() +
  
  theme(
    
    plot.title = element_text(
      hjust = 0.5,
      size = 18,
      face = "bold"
    ),
    
    axis.text = element_text(
      size = 16,
      face = "bold"
    ),
    
    axis.title = element_text(
      size = 18,
      face = "bold"
    ),
    
    panel.grid.major.x = element_blank(),
    
    panel.grid.major.y = element_line(
      color = "gray85",
      linewidth = 0.5
    ),
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1.5
    )
  )

ggsave(
  "results/Figure_S6.tiff",
  plot = p1,
  width = 6,
  height = 5,
  dpi = 300,
  compression = "lzw"
)


#Seminal Roots

freq_SR <- DF_maize %>%
  group_by(SR) %>%
  summarise(n = n(), .groups = "drop")


freq_SR <- freq_SR %>%
  mutate(
    percent = (n / sum(n)) * 100
  )

freq_SR$SR <- factor(freq_SR$SR, levels = sort(unique(freq_SR$SR)))


p2 <- ggplot(
  freq_SR,
  aes(x = SR, y = percent)
) +
  
  geom_bar(
    stat = "identity",
    fill = "black",
    color = "white",
    linewidth = 0.3
  ) +
  
  labs(
    title = "",
    x = "Number of seminal roots per seedling",
    y = "Seedlings (%)"
  ) +
  
  theme_minimal() +
  
  theme(
    
    plot.title = element_text(
      hjust = 0.5,
      size = 18,
      face = "bold"
    ),
    
    axis.text = element_text(
      size = 16,
      face = "bold"
    ),
    
    axis.title = element_text(
      size = 18,
      face = "bold"
    ),
    
    panel.grid.major.x = element_blank(),
    
    panel.grid.major.y = element_line(
      color = "gray85",
      linewidth = 0.5
    ),
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1.5
    )
  )


ggsave(
  "results/Figure_S7.tiff",
  plot = p2,
  width = 6,
  height = 5,
  dpi = 300,
  compression = "lzw"
)

#===============================================================================

# Figure S8 Breake coleoptile tip percent for groups 

mean_BreakCol <- DF %>%
  group_by(Group) %>%
  summarise(
    mean = mean(BreakCol, na.rm = TRUE) * 1,
    sd = sd(BreakCol, na.rm = TRUE) * 1,
    n = n(),
    se = sd / sqrt(n),
    .groups = "drop"
  )


mean_BreakCol$Group <- factor(mean_BreakCol$Group,
                              levels = c("Maize", "Teosinte", "Tripsacum"))


model <- aov(BreakCol ~ Group, data = DF)

res <- agricolae::HSD.test(model, "Group", group = TRUE)

letras <- res$groups
letras$Group <- rownames(letras)
colnames(letras)[2] <- "letra"

pos <- mean_BreakCol %>%
  left_join(letras, by = "Group")

p8 <- ggplot(
  mean_BreakCol,
  aes(x = Group, y = mean, fill = Group)
) +
  
  geom_bar(
    stat = "identity",
    color = "black",
    linewidth = 0.5
  ) +
  
  geom_text(
    data = pos,
    aes(
      x = Group,
      y = mean + se + 10,
      label = letra
    ),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 7
  ) +
  
  scale_fill_manual(values = c(
    "Maize" = "orange",
    "Teosinte" = "darkgreen",
    "Tripsacum" = "blue"
  )) +
  
  scale_x_discrete(labels = c(
    "Maize" = expression(bold("Maize")),
    "Teosinte" = expression(bolditalic("Teosinte")),
    "Tripsacum" = expression(bolditalic("Tripsacum"))
  )) +
  
  labs(
    title = "",
    x = "",
    y = "Coleoptile tip breakage (%)"
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
      size = 16,
      face = "bold"
    ),
    
    panel.grid.major.x = element_blank(),
    
    panel.grid.major.y = element_line(
      color = "gray85",
      linewidth = 0.5
    ),
    
    legend.position = "none",
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1.5
    )
  )

ggsave(
  filename = "results/Figure_S8.tiff",
  plot = p8,
  width = 5,
  height = 4,
  dpi = 300,
  compression = "lzw"
)

