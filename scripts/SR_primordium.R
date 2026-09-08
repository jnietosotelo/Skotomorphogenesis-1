

############# Analysis and figures for Evolution of skotomorphogenic development in seedlings of Zea and Tripsacum

##Clean session 

rm(list = ls())
graphics.off()
cat("\014")


##Load libraries 

library(ggplot2)
library(dplyr)
library(agricolae)
library(patchwork)
library(cowplot)
library(tidyr)
library(tibble)

# All results generated in this script will be saved in the "results" directory

dir.create("results", showWarnings = FALSE)

#Load data 

DF <- read.csv("data/PrimordiumSR.csv", stringsAsFactors = TRUE)


DF <- DF %>%
  filter(!is.na(AccEvolID), !is.na(Origen), !is.na(Rasem))


#Sort

DF$Origen <- factor(
  DF$Origen,
  levels = c("Primordia", "Darkness", "Photoperiod")
)

# FUNCTION: ANOVA + TUKEY + LETERS

get_anova_full <- function(data, accesion_id){
  
  df_sub <- data %>%
    filter(AccEvolID == accesion_id)
  
  if(nrow(df_sub) < 3){
    return(NULL)
  }
  

  # ANOVA
  
  model <- aov(Rasem ~ Origen, data = df_sub)
  
  anova_sum <- summary(model)[[1]]
  
  F_val <- anova_sum["Origen", "F value"]
  p_val <- anova_sum["Origen", "Pr(>F)"]
  

  # TUKEY 

  
  tuk <- agricolae::HSD.test(model, "Origen", group = FALSE)
  comp <- data.frame(tuk$comparison)
  comp$Comparison <- rownames(comp)
  rownames(comp) <- NULL
  
 
  # LETERS
  
  tuk_groups <- agricolae::HSD.test(model, "Origen", group = TRUE)
  letras <- tuk_groups$groups
  letras$Origen <- rownames(letras)
  
  colnames(letras)[1] <- "Mean"
  colnames(letras)[2] <- "Letter"
  
  rownames(letras) <- NULL

  
  comp <- comp %>%
    mutate(
      Origen1 = sub(" - .*", "", Comparison),
      Origen2 = sub(".*- ", "", Comparison)
    ) %>%
    
    left_join(letras %>% dplyr::select(Origen, Letter), 
              by = c("Origen1" = "Origen")) %>%
    rename(Letter1 = Letter) %>%
    
    left_join(letras %>% dplyr::select(Origen, Letter), 
              by = c("Origen2" = "Origen")) %>%
    rename(Letter2 = Letter)
  
#General info
  
  comp$AccEvolID <- accesion_id
  comp$F_value <- F_val
  comp$p_value_anova <- p_val
  
  # sort
  
  comp <- comp %>%
    dplyr::select(
      AccEvolID,
      F_value,
      p_value_anova,
      Comparison,
      Origen1, Origen2,
      Letter1, Letter2,
      everything()
    )
  
  return(comp)
}

##ID

accesiones <- levels(as.factor(DF$AccEvolID))


#RUN

anova_list <- lapply(accesiones, function(acc){
  
  message("Procesando ANOVA: ", acc)
  
  tryCatch({
    get_anova_full(DF, acc)
  }, error = function(e){
    message("Error en: ", acc)
    print(e)
    traceback()
    return(NULL)
  })
})

anova_list <- anova_list[!sapply(anova_list, is.null)]

##Table

tabla_anova_final <- do.call(rbind, anova_list)

#Save table 

dir.create("results", showWarnings = FALSE)

write.csv(
  tabla_anova_final,
  "results/Table_ANOVA_Tukey_Primordia.csv",
  row.names = FALSE
)


plot_bubble_accesion <- function(data, accesion_id) {
  
  # Filter by ID
  
  df_sub <- data %>%
    filter(AccEvolID == accesion_id)
  
  # Anova
  
  model <- aov(Rasem ~ Origen, data = df_sub)
  
  # Tukey
  
  tukey <- agricolae::HSD.test(model, "Origen", group = TRUE)
  letras <- tukey$groups
  letras$Origen <- rownames(letras)
  colnames(letras)[1] <- "mean"
  colnames(letras)[2] <- "letra"
  
  pos <- df_sub %>%
    group_by(Origen) %>%
    summarise(ymax = max(Rasem, na.rm = TRUE), .groups = "drop") %>%
    left_join(letras, by = "Origen")
  
  # Plot
  
  p <- ggplot(df_sub, aes(x = Origen, y = Rasem, color = Origen)) +
    
    geom_count(alpha = 0.7) +
    scale_size_area(max_size = 8) +
    
    geom_text(
      data = pos,
      aes(x = Origen, y = ymax + 1.5, label = letra),
      inherit.aes = FALSE,
      size = 8,     
      color = "black",
      fontface = "bold"
    ) +
    
    scale_color_manual(values = c(
      "Primordia" = "darkgreen",
      "Darkness" = "black",
      "Photoperiod" = "orange"
    )) +
    
    labs(
      title = accesion_id,
      x = "",
      y = "No. Seminal roots"
    ) +
    
    theme_minimal() +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 18,
        face = "bold"),
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
  
  return(p)
}

##Accesion ID 

accesiones <- levels(as.factor(DF$AccEvolID))

## Plots whit tags

tags <- LETTERS[1:length(accesiones)]

plots <- lapply(seq_along(accesiones), function(i) {
  
  acc <- accesiones[i]
  tag <- tags[i]
  
  message("Procesando: ", acc)
  
  tryCatch({
    
    plot_bubble_accesion(DF, acc) +
      labs(tag = tag) +
      theme(
        plot.tag = element_text(size = 18, face = "bold"),
        plot.tag.position = c(0.02, 0.98)
      )
    
  }, error = function(e) {
    message("Error en: ", acc)
    print(e)
    NULL
  })
})

# Remove NULL

plots <- plots[!sapply(plots, is.null)]



final_plot <- plot_grid(plotlist = plots, ncol = 3)


##Save plot

ggsave(
  "results/Figure_2.tiff",
  plot = final_plot,
  dpi = 300,
  width = 17,
  height = 5,
  units = "in",
  bg = "white",
  compression = "lzw"
)

