############# Analysis and figures for Evolution of skotomorphogenic development in seedlings of Zea and Tripsacum


##Clean session 

rm(list = ls())
graphics.off()
cat("\014")


## Load libraries

library(dplyr)
library(ggplot2)
library(cowplot)


select <- dplyr::select
filter <- dplyr::filter


# All results generated in this script will be saved in the "results" directory

dir.create("results", showWarnings = FALSE)


# Data: mean trait values used for correlation with collection altitude



DF <- read.csv("data/EvolEscoMeanAlt.csv", stringsAsFactors = TRUE)




# Variables
vars_cor <- c("Stem","Mesocotyl","Coleoptile","FPL","FPLCol",
              "PR","ARCN","SR","BreakCol","RDSD")

# # Function: clean dataset by removing missing values (NA) for each variable pair used in correlation analysis


clean_data <- function(data, x, y){
  data %>%
    dplyr::select(all_of(c(x, y))) %>%
    na.omit()
}

#Functoion correlation test

get_cor_stats <- function(data, x, y) {
  
  df_clean <- clean_data(data, x, y)
  
  if(nrow(df_clean) < 3){
    return(list(r = NA, p = NA, data = df_clean))
  }
  
  test <- cor.test(df_clean[[x]], df_clean[[y]], method = "pearson")
  
  r <- round(test$estimate, 2)
  
  p <- ifelse(test$p.value < 0.0001,
              "<0.0001",
              round(test$p.value, 4))
  
  return(list(r = r, p = p, data = df_clean))
}


# Function_ X Y plot

plot_corr <- function(data, x, y, ylab,
                      y_limits = NULL, y_breaks = NULL,
                      label_pos = "top_left") {
  
  stats <- get_cor_stats(data, x, y)
  df_plot <- stats$data
  
  if(nrow(df_plot) < 3){
    return(ggplot() + theme_void())
  }
  
  # Posición etiqueta
  if(label_pos == "top_left"){
    xpos <- -Inf; ypos <- Inf; h <- -0.1; v <- 1.1
  } else {
    xpos <- -Inf; ypos <- -Inf; h <- -0.1; v <- -0.5
  }
  
  p <- ggplot(df_plot, aes(x = .data[[x]], y = .data[[y]])) +
    
    geom_point(color = "black", size = 2, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    
    annotate(
      "text",
      x = xpos,
      y = ypos,
      label = paste0("r = ", stats$r, "\np = ", stats$p),
      hjust = h,
      vjust = v,
      size = 6,
      fontface = "bold"
    ) +
    
    labs(
      x = "Altitude (m a.s.l)",
      y = ylab
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
        size = 1.5)
    )
  
  if(!is.null(y_limits)){
    p <- p + scale_y_continuous(breaks = y_breaks, limits = y_limits)
  }
  
  p <- p + scale_x_continuous(breaks = scales::pretty_breaks(n = 5))
  
  return(p)
}

# SUBSETS 

DF_maize <- DF %>%
  dplyr::filter(Group == "maize") %>%
  dplyr::select(Altitude, dplyr::all_of(vars_cor))

DF_trip <- DF %>%
  dplyr::filter(Group == "Tripsacum") %>%
  dplyr::select(Altitude, dplyr::all_of(vars_cor))


# FIGURE (Maize Correlation)

g1 <- plot_corr(DF_maize, "Altitude", "Mesocotyl",
                "Mesocotyl length (cm)", c(0,25), seq(0,25,5))

g2 <- plot_corr(DF_maize, "Altitude", "Coleoptile",
                "Coleoptile length (cm)", c(0,10), seq(0,10,2))

g3 <- plot_corr(DF_maize, "Altitude", "FPL",
                "FPL length (cm)", c(0,25), seq(0,25,5))

g4 <- plot_corr(DF_maize, "Altitude", "FPLCol",
                "FPL-Col length (cm)", c(0,17), seq(0,17,2))

g5 <- plot_corr(DF_maize, "Altitude", "ARCN",
                "Number of ARCN per seedling", c(0,4), seq(0,4,1))

g6 <- plot_corr(DF_maize, "Altitude", "RDSD",
                "RDSD length (cm)", c(0,25), seq(0,25,5))

g7 <- plot_corr(DF_maize, "Altitude", "BreakCol",
                "Coleoptile tip breakage (%)",
                c(0,120), seq(0,100,20), "bottom_left")

plots <- list(
  g1 + labs(tag="A")+
    theme(plot.tag = element_text(size = 18, face = "bold")),
  g2 + labs(tag="B")+
    theme(plot.tag = element_text(size = 18, face = "bold")),
  g3 + labs(tag="C")+
    theme(plot.tag = element_text(size = 18, face = "bold")),
  g4 + labs(tag="D")+
    theme(plot.tag = element_text(size = 18, face = "bold")),
  g5 + labs(tag="E")+
    theme(plot.tag = element_text(size = 18, face = "bold")),
  g6 + labs(tag="F")+
    theme(plot.tag = element_text(size = 18, face = "bold")),
  g7 + labs(tag="G")
)

# Completar paneles vacíos
while(length(plots) < 9){
  plots[[length(plots)+1]] <- ggplot() + theme_void()
}

final_plot <- plot_grid(plotlist = plots, ncol = 3)

# Save figure
ggsave("results/Figure_5.tiff",
       plot = final_plot,
       dpi = 300,
       width = 16,
       height = 13,
       units = "in",
       bg = "white",
       compression = "lzw")

# TABLE MAIZE altitudeXTraits Correlation

cor_maize <- lapply(vars_cor, function(var) {
  stats <- get_cor_stats(DF_maize, "Altitude", var)
  data.frame(Variable = var, r = stats$r, p = stats$p)
})

cor_maize <- do.call(rbind, cor_maize)

write.csv(cor_maize,
          "results/Table_S12_Correlation_Altitude_Maize.csv",
          row.names = FALSE)

# Table tripsacum correlation

cor_trip <- lapply(vars_cor, function(var) {
  stats <- get_cor_stats(DF_trip, "Altitude", var)
  data.frame(Variable = var, r = stats$r, p = stats$p)
})

cor_trip <- do.call(rbind, cor_trip)

write.csv(cor_trip,
          "results/Table_S12_Correlation_Altitude_Tripsacum.csv",
          row.names = FALSE)


