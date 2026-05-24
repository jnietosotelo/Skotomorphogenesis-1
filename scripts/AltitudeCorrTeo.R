############# Analysis and figures for Evolution of skotomorphogenic development in seedlings of Zea and Tripsacum


##Clean session 

rm(list = ls())
graphics.off()
cat("\014")


## Load libraries

library(dplyr)
library(ggplot2)
library(cowplot)



# All results generated in this script will be saved in the "results" directory

dir.create("results", showWarnings = FALSE)


# Data: mean trait values used for correlation with collection altitude

DF <- read.csv("data/EvolEscoMeanAlt.csv", stringsAsFactors = TRUE)

#Filter of teocinte 

DF <- DF %>%
  filter(Group == "Teosinte") %>%
  filter(!is.na(Altitude))


#Color of teocinte specie

col_subSp <- c(
  "Z. mexicana" = "blue",
  "Z. parviglumis" = "brown",
  "Z. sp" = "purple",
  "Z. diploperennis" = "black",
  "Z. huehuetenanguensis" = "gray"
)



vars_cor <- c(
  "Stem","Mesocotyl","Coleoptile",
  "FPL", "FPLCol", "PR",
  "ARCN","SR","BreakCol","RDSD"
)

# Función correlación
# =========================
get_cor_stats <- function(data, x, y) {
  
  df_clean <- data %>%
    dplyr::select(all_of(c(x, y, "SubSp"))) %>% 
    na.omit()
  
  test <- cor.test(df_clean[[x]], df_clean[[y]], method = "pearson")
  
  r <- round(test$estimate, 2)
  
  p <- ifelse(test$p.value < 0.0001,
              "<0.0001",
              round(test$p.value, 4))
  
  return(list(r = r, p = p, data = df_clean))
}

#Correlation analysis altitude x variables

cor_results <- lapply(vars_cor, function(var) {
  

  if(!var %in% colnames(DF)){
    return(data.frame(
      Variable = var,
      r = NA,
      p = NA,
      Note = "Variable no encontrada"
    ))
  }
  
  stats <- get_cor_stats(DF, "Altitude", var)
  
  data.frame(
    Variable = var,
    r = stats$r,
    p = stats$p,
    Note = NA
  )
})

# Joint results

cor_results_df <- do.call(rbind, cor_results)

#Save table 

dir.create("results", showWarnings = FALSE)

write.csv(
  cor_results_df,
  "results/Table_S12_Correlation_Altitude_Teosinte.csv",
  row.names = FALSE
)

# =========================
# Ver resultados
# =========================
print(cor_results_df)
# =========================




# =========================
# Función gráfica
# =========================
plot_corr_teo <- function(data, x, y, ylab,
                          y_limits = NULL,
                          y_breaks = NULL,
                          label_pos = "top_left") {
  
  stats <- get_cor_stats(data, x, y)
  df_plot <- stats$data
  
  if(label_pos == "top_left"){
    
    xpos <- -Inf
    ypos <- Inf
    h <- -0.1
    v <- 1.1
    
  } else {
    
    xpos <- -Inf
    ypos <- -Inf
    h <- -0.1
    v <- -0.5
  }
  
  p <- ggplot(
    df_plot,
    aes(
      x = .data[[x]],
      y = .data[[y]],
      color = SubSp
    )
  ) +
    
    geom_point(
      size = 3,
      alpha = 0.7
    ) +
    
    geom_smooth(
      method = "lm",
      se = FALSE,
      color = "black",
      linewidth = 1.2
    ) +
    
    annotate(
      "text",
      x = xpos,
      y = ypos,
      label = paste0(
        "r = ", stats$r,
        "\np = ", stats$p
      ),
      hjust = h,
      vjust = v,
      size = 6,
      fontface = "bold"
    ) +
    
    scale_color_manual(values = col_subSp) +
    
    labs(
      x = "Altitude (m a.s.l)",
      y = ylab
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
      
      panel.grid.major = element_line(
        color = "gray85",
        linewidth = 0.5
      ),
      
      panel.grid.minor = element_blank(),
      
      legend.position = "none",
      
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 1.5
      )
    )
  
  if(!is.null(y_limits)){
    
    p <- p + scale_y_continuous(
      breaks = y_breaks,
      limits = y_limits
    )
  }
  
  p <- p + scale_x_continuous(
    breaks = seq(0, 2800, 500),
    limits = c(0, 2800)
  )
  
  return(p)
}

# Figure A

g_alt <- DF %>%
  
  dplyr::arrange(Altitude) %>%
  
  dplyr::mutate(
    Name = factor(Name, levels = Name)
  ) %>%
  
  ggplot(aes(
    x = Name,
    y = Altitude,
    color = SubSp
  )) +
  
  geom_point(
    size = 3
  ) +
  
  scale_color_manual(
    values = col_subSp,
    
    labels = c(
      "Z. mexicana" = expression(bold(italic("Z. mexicana"))),
      "Z. parviglumis" = expression(bold(italic("Z. parviglumis"))),
      "Z. sp" = expression(bold(italic("Z. sp"))),
      "Z. diploperennis" = expression(bold(italic("Z. diploperennis"))),
      "Z. huehuetenanguensis" = expression(bold(italic("Z. huehuetenanguensis")))
    )
  ) +
  
  labs(
    x = "",
    y = "Altitude (m a.s.l)",
    color = NULL
  ) +
  
  scale_y_continuous(
    limits = c(0, 3100),
    breaks = seq(0, 3100, 500)
  ) +
  
  theme_minimal() +
  
  theme(
    
    axis.text.x = element_text(
      angle = 70,
      hjust = 1,
      size = 12,
      face = "bold"
    ),
    
    axis.text.y = element_text(
      size = 16,
      face = "bold"
    ),
    
    axis.title.y = element_text(
      size = 18,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 12,
      face = "bold"
    ),
    
    panel.grid.major = element_line(
      color = "gray85",
      linewidth = 0.5
    ),
    
    panel.grid.minor = element_blank(),
    
    legend.position = c(0.95, 0.05),
    
    legend.justification = c(
      "right",
      "bottom"
    ),
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1.5
    )
  )

# Correlation plots

g1 <- plot_corr_teo(DF, "Altitude", "FPL", "FPL length (cm)", c(0,15), seq(0,15,2))
g2 <- plot_corr_teo(DF, "Altitude", "FPLCol", "FPL-Col length (cm)", c(0,8), seq(0,8,2))
g3 <- plot_corr_teo(DF, "Altitude", "PR", "PR length (cm)", c(0,5), seq(0,5,1))
g4 <- plot_corr_teo(DF, "Altitude", "BreakCol", "Coleoptile tip breakage (%)",
                    c(0,110), seq(0,100,20), label_pos = "bottom_left")



tag_theme <- theme(
  plot.tag = element_text(size = 18, face = "bold"),
  plot.tag.position = c(0.01, .99),
  plot.margin = margin(15, 15, 15, 15)
)

g_alt <- g_alt + labs(tag = "A") + tag_theme
g1 <- g1 + labs(tag = "B") + tag_theme
g2 <- g2 + labs(tag = "C") + tag_theme
g3 <- g3 + labs(tag = "D") + tag_theme
g4 <- g4 + labs(tag = "E") + tag_theme

g1 <- g1 + theme(
  plot.tag.position = c(-0.01, 1.1)   
)

g2 <- g2 + theme(
  plot.tag.position = c(-0.03, 1.1)   
)

g3 <- g3 + theme(
  plot.tag.position = c(-0.01, 1.1)   
)

g4 <- g4 + theme(
  plot.tag.position = c(-0.03, 1.1)   
)

bottom_panel <- plot_grid(
  g1, g2,
  g3, g4,
  ncol = 2
)


final_plot <- plot_grid(
  g_alt,
  bottom_panel,
  ncol = 1,
  rel_heights = c(2, 1.5)
)

final_plot <- final_plot +
  theme(plot.background = element_rect(fill = "white", color = NA))


##Save plot

ggsave(
  "results/Figure_S13.tiff",
  plot = final_plot,
  dpi = 300,
  compression = "lzw",
  width = 10,
  height = 18,
  units = "in",
  bg = "white"
)