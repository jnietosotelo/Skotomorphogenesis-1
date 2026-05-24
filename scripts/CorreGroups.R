

############# Analysis and figures for Evolution of skotomorphogenic development in seedlings of Zea and Tripsacum

##Clean session 

rm(list = ls())
graphics.off()
cat("\014")



## Load libraries

library(dplyr)
library(tidyr)


# All results generated in this script will be saved in the "results" directory

dir.create("results", showWarnings = FALSE)


# Data
DF <- read.csv("data/EvolEscoDB.csv")



# Variables 

vars <- c(
  "Stem","Mesocotyl","Coleoptile",
  "FPL","FPLCol","PR",
  "ARCN","SR","RDSD"
)

# SUBSETS for group

DF_maize <- DF %>% filter(Group == "Maize")
DF_teo   <- DF %>% filter(Group == "Teosinte")
DF_trip  <- DF %>% filter(Group == "Tripsacum")

# FUNCTION FUL MATRIZ (r y p value) 

cor_matrix <- function(data, group_name){
  
  results <- list()
  
  for(i in 1:length(vars)){
    for(j in 1:length(vars)){
      
      var1 <- vars[i]
      var2 <- vars[j]
      
      if(!(var1 %in% colnames(data)) | !(var2 %in% colnames(data))){
        r <- NA
        p <- NA
        
      } else {
        
        df_pair <- data %>%
          dplyr::select(all_of(c(var1, var2))) %>%
          na.omit()
        
        if(ncol(df_pair) < 2 || nrow(df_pair) < 3){
          r <- NA
          p <- NA
          
        } else if(sd(df_pair[[1]]) == 0 | sd(df_pair[[2]]) == 0){
          r <- NA
          p <- NA
          
        } else {
          test <- cor.test(df_pair[[1]], df_pair[[2]], method = "pearson")
          r <- round(test$estimate, 3)
          p <- round(test$p.value, 5)
        }
      }
      
      results[[length(results) + 1]] <- data.frame(
        Group = group_name,
        Var1 = var1,
        Var2 = var2,
        r = r,
        p_value = p
      )
    }
  }
  
  do.call(rbind, results)
}

# Full matrix analysis

cor_maize <- cor_matrix(DF_maize, "Maize")
cor_teo   <- cor_matrix(DF_teo, "Teosinte")
cor_trip  <- cor_matrix(DF_trip, "Tripsacum")

tabla_cor_full <- bind_rows(cor_maize, cor_teo, cor_trip)


tabla_cor_full$Var1 <- factor(tabla_cor_full$Var1, levels = vars)
tabla_cor_full$Var2 <- factor(tabla_cor_full$Var2, levels = vars)

tabla_cor_full <- tabla_cor_full %>%
  arrange(Group, Var1, Var2)

# Save table

write.csv(
  tabla_cor_full,
  "results/Table_Correlation_Group_Full_Matrix.csv",
  row.names = FALSE
)


# Function matrix r 

cor_matrix_r <- function(data){
  
  n <- length(vars)
  
  mat <- matrix("", nrow = n, ncol = n)
  colnames(mat) <- vars
  rownames(mat) <- vars
  
  for(i in 1:n){
    for(j in 1:n){
      
      if(i == j){
        mat[i,j] <- "1"
        
      } else if(i > j){
        
        var1 <- vars[i]
        var2 <- vars[j]
        
        if(!(var1 %in% colnames(data)) | !(var2 %in% colnames(data))){
          r <- NA
          
        } else {
          
          df_pair <- data %>%
            dplyr::select(all_of(c(var1, var2))) %>%
            na.omit()
          
          if(ncol(df_pair) < 2 || nrow(df_pair) < 3){
            r <- NA
            
          } else if(sd(df_pair[[1]]) == 0 | sd(df_pair[[2]]) == 0){
            r <- NA
            
          } else {
            r <- cor(df_pair[[1]], df_pair[[2]])
          }
        }
        
        mat[i,j] <- ifelse(is.na(r), "", round(r, 2))
        
      } else {
        mat[i,j] <- ""
      }
    }
  }
  
  mat_df <- as.data.frame(mat)
  mat_df <- cbind(Variable = rownames(mat_df), mat_df)
  rownames(mat_df) <- NULL
  
  return(mat_df)
}

# Correlation matrix r value

mat_maize <- cor_matrix_r(DF_maize)
mat_teo   <- cor_matrix_r(DF_teo)
mat_trip  <- cor_matrix_r(DF_trip)

# Save matrix correlation by groups

write.csv(mat_maize, "results/Table_Correlation_r_Maize.csv", row.names = FALSE)
write.csv(mat_teo,   "results/Table_Correlation_r_Teosinte.csv", row.names = FALSE)
write.csv(mat_trip,  "results/Table_Correlation_r_Tripsacum.csv", row.names = FALSE)
