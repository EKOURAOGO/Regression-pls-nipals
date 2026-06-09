# =========================================================
# 9. MODÉLISATION UNIVARIÉE 
# =========================================================

library(pls)
library(ggplot2)
library(knitr)
library(reshape2)
library(gridExtra)
library(grid)
library(MASS)

run_univariate_models <- function(X_matrix, Y_matrix, label){
  
  cat("\n==============================\n")
  cat("ANALYSE :", label, "\n")
  cat("==============================\n")
  
  # Retrait contrainte compositionnelle
  if("p_blanc1" %in% colnames(X_matrix)) {
    X_matrix <- X_matrix[, colnames(X_matrix) != "p_blanc1"]
    cat("\n⚠️  Variable 'p_blanc1' retirée (contrainte compositionnelle)\n")
  }
  
  X_scaled <- scale(X_matrix)
  n <- nrow(X_scaled)
  p <- ncol(X_scaled)
  max_comp <- min(10, p-1, n-2)
  
  Y_vars <- colnames(Y_matrix)
  results_list <- list()
  plot_list <- list()
  
  # =====================================================
  # DIAGNOSTICS MULTICOLINÉARITÉ
  # =====================================================
  
  cat("\n--- DIAGNOSTICS MULTICOLINÉARITÉ ---\n")
  
  cor_matrix <- cor(X_scaled)
  print(round(cor_matrix,2))
  
  condition_number <- kappa(X_scaled)
  cat(sprintf("\nIndice de conditionnement : %.2f\n", condition_number))
  
  vif_values <- tryCatch({
    diag(solve(cor_matrix))
  }, error=function(e){
    diag(MASS::ginv(cor_matrix))
  })
  
  vif_df <- data.frame(
    Variable = colnames(X_matrix),
    VIF = round(vif_values,2)
  )
  print(vif_df)
  
  # =====================================================
  # MODÉLISATION
  # =====================================================
  
  for(Y in Y_vars){
    
    cat("\n--- Variable :", Y, "---\n")
    
    y <- Y_matrix[,Y]
    
    # ----------------------
    # MCO
    # ----------------------
    model_mco <- lm(y ~ X_scaled)
    r2_mco <- summary(model_mco)$r.squared
    r2_adj_mco <- summary(model_mco)$adj.r.squared
    coef_mco <- coef(model_mco)[-1]
    
    cat(sprintf("MCO  : R² = %.4f | R² adj = %.4f\n",
                r2_mco, r2_adj_mco))
    
    # ----------------------
    # PCR
    # ----------------------
    model_pcr <- pcr(y ~ X_scaled,
                     validation="LOO",
                     scale=FALSE,
                     ncomp=max_comp)
    
    q2_pcr_vec <- as.numeric(drop(R2(model_pcr, estimate="CV")$val))[-1]
    ncomp_pcr <- which.max(q2_pcr_vec)
    q2_pcr <- q2_pcr_vec[ncomp_pcr]
    
    cat(sprintf("PCR  : %d composantes | Q² = %.4f\n",
                ncomp_pcr, q2_pcr))
    
    # ----------------------
    # PLS
    # ----------------------
    model_pls <- plsr(y ~ X_scaled,
                      validation="LOO",
                      scale=FALSE,
                      ncomp=max_comp)
    
    q2_pls_vec <- as.numeric(drop(R2(model_pls, estimate="CV")$val))[-1]
    ncomp_pls <- which.max(q2_pls_vec)
    q2_pls <- q2_pls_vec[ncomp_pls]
    
    cat(sprintf("PLS  : %d composantes | Q² = %.4f\n",
                ncomp_pls, q2_pls))
    
    # =====================================================
    # INTERPRÉTATION PLS
    # =====================================================
    
    coef_pls <- as.vector(coef(model_pls,
                               ncomp=ncomp_pls,
                               intercept=FALSE))
    
    # -------- VIP robuste
    W <- as.matrix(model_pls$loading.weights[,1:ncomp_pls,drop=FALSE])
    T_scores <- as.matrix(model_pls$scores[,1:ncomp_pls,drop=FALSE])
    SS <- colSums(T_scores^2)
    
    vip <- numeric(p)
    for(j in 1:p){
      vip[j] <- sqrt(p * sum(SS * W[j,]^2) / sum(SS))
    }
    
    vip_df <- data.frame(
      Variable = colnames(X_matrix),
      VIP = round(vip,3)
    )
    vip_df <- vip_df[order(-vip_df$VIP),]
    
    cat("\nVIP (Top 5) :\n")
    print(head(vip_df,5))
    
    # -------- Loadings
    loadings <- model_pls$loadings[,1:min(2,ncomp_pls),drop=FALSE]
    
    loadings_df <- data.frame(
      Variable = colnames(X_matrix),
      Comp1 = round(loadings[,1],3)
    )
    
    if(ncol(loadings)>=2){
      loadings_df$Comp2 <- round(loadings[,2],3)
    }
    
    cat("\nLoadings (Comp 1 & 2) :\n")
    print(loadings_df)
    
    # -------- Comparaison coefficients
    cat("\n--- Comparaison MCO vs PLS ---\n")
    
    comp_coef <- data.frame(
      Variable = colnames(X_matrix),
      MCO = round(coef_mco,4),
      PLS = round(coef_pls,4),
      Diff = round(abs(coef_mco - coef_pls),4)
    )
    
    print(head(comp_coef[order(-comp_coef$Diff),],5))
    
    # =====================================================
    # GRAPHIQUE COMPARATIF
    # =====================================================
    
    df_plot <- data.frame(
      ncomp = 1:max_comp,
      PCR = q2_pcr_vec,
      PLS = q2_pls_vec
    )
    
    df_long <- melt(df_plot, id.vars="ncomp")
    
    y_min <- min(c(q2_pcr_vec,q2_pls_vec,r2_adj_mco), na.rm=TRUE)
    y_max <- max(c(q2_pcr_vec,q2_pls_vec,r2_adj_mco), na.rm=TRUE)
    
    g <- ggplot(df_long, aes(x=ncomp, y=value, color=variable)) +
      geom_line(linewidth=1.2) +
      geom_point(size=2) +
      
      geom_hline(yintercept=r2_adj_mco,
                 linetype="dashed",
                 color="darkgreen",
                 linewidth=1) +
      
      annotate("text",
               x=max_comp*0.65,
               y=r2_adj_mco+0.02,
               label="R² adj MCO",
               color="darkgreen",
               size=3) +
      
      annotate("point",
               x=ncomp_pls,
               y=q2_pls,
               color="blue",
               size=4,
               shape=8) +
      
      annotate("point",
               x=ncomp_pcr,
               y=q2_pcr,
               color="red",
               size=4,
               shape=8) +
      
      scale_color_manual(values=c("PCR"="red","PLS"="blue")) +
      
      coord_cartesian(ylim=c(max(0,y_min-0.05),
                             min(1,y_max+0.05))) +
      
      labs(title=gsub("p_","",Y),
           x="Nombre de composantes",
           y=expression(Q^2~"(LOO)"),
           color="Méthode") +
      
      theme_minimal(base_size=12) +
      theme(legend.position="bottom",
            plot.title=element_text(face="bold"))
    
    plot_list[[Y]] <- g
    
    # =====================================================
    # STOCKAGE
    # =====================================================
    
    results_list[[Y]] <- data.frame(
      Variable=Y,
      ncomp_PCR=ncomp_pcr,
      ncomp_PLS=ncomp_pls,
      R2_MCO=round(r2_mco,4),
      R2adj_MCO=round(r2_adj_mco,4),
      Q2_PCR=round(q2_pcr,4),
      Q2_PLS=round(q2_pls,4),
      VIP_max=round(max(vip),3),
      VIP_var=vip_df$Variable[1]
    )
  }
  
  # =====================================================
  # AFFICHAGE COMBINÉ
  # =====================================================
  
  grid.arrange(
    grobs=plot_list,
    ncol=2,
    top=textGrob(paste("Évolution du Q² –",label),
                 gp=gpar(fontsize=16,fontface="bold"))
  )
  
  results <- do.call(rbind,results_list)
  rownames(results) <- NULL
  
  cat("\n==============================\n")
  cat("TABLEAU RÉCAPITULATIF -",label,"\n")
  cat("==============================\n")
  print(knitr::kable(results,digits=4))
  
  return(invisible(results))
}

# =========================================================
# EXÉCUTION
# =========================================================

cat("\nDÉMARRAGE DES MODÈLES...\n")

results_all <- run_univariate_models(
  data_all_xy$X,
  data_all_xy$Y,
  "107 DÉPARTEMENTS"
)

results_metro <- run_univariate_models(
  data_metro_xy$X,
  data_metro_xy$Y,
  "MÉTROPOLE"
)

# =========================================================
# 10. MODÉLISATION MULTIVARIÉE PLS2 - VERSION CORRIGÉE ET FINALISÉE
# =========================================================

library(pls)
library(ggplot2)

cat("\n=== MODÉLISATION PLS2 (MULTIVARIÉE) ===\n")

# =========================================================
# PRÉPARATION DES DONNÉES
# =========================================================

X_pls2 <- data_all_xy$X
Y_pls2 <- data_all_xy$Y

if("p_blanc1" %in% colnames(X_pls2)) {
  X_pls2 <- X_pls2[, -which(colnames(X_pls2) == "p_blanc1")]
  cat("Variable 'p_blanc1' retirée\n")
}

n <- nrow(X_pls2)
p <- ncol(X_pls2)
q <- ncol(Y_pls2)
max_comp <- min(12, p-1, n-2)

cat("Dimensions :", n, "obs,", p, "X,", q, "Y\n")

# =========================================================
# MODÈLE PLS2
# =========================================================

model_pls2_all <- plsr(Y_pls2 ~ X_pls2, 
                       scale = TRUE, 
                       ncomp = max_comp,
                       validation = "CV")

# =========================================================
# CHOIX DU NOMBRE DE COMPOSANTES
# =========================================================

validationplot(model_pls2_all, val.type = "R2",
               main = "Choix du nombre de composantes – PLS2")

Q2_array <- R2(model_pls2_all, estimate = "CV")$val
Q2_mean <- numeric(max_comp)

for(k in 1:max_comp) {
  Q2_mean[k] <- mean(Q2_array[1, , k+1])
}

best_n_pls2 <- which.max(Q2_mean)
best_q2_pls2 <- Q2_mean[best_n_pls2]

cat("\nQ² moyen par composante:\n")
print(round(Q2_mean, 4))
cat("\nOptimal:", best_n_pls2, "composantes, Q² =", round(best_q2_pls2, 4), "\n")

# Graphique Q²
df_q2 <- data.frame(ncomp = 1:max_comp, Q2 = Q2_mean)

ggplot(df_q2, aes(x = ncomp, y = Q2)) +
  geom_line(color = "blue", linewidth = 1.2) +
  geom_point(color = "blue", size = 3) +
  geom_point(data = df_q2[best_n_pls2, ], aes(x = ncomp, y = Q2),
             color = "red", size = 5, shape = 8) +
  labs(title = "Choix du nombre de composantes PLS2",
       subtitle = paste("Optimal =", best_n_pls2, "comp., Q² =", round(best_q2_pls2, 4)),
       x = "Nombre de composantes", y = "Q² moyen (CV)") +
  theme_minimal()

# =========================================================
# MODÈLE FINAL
# =========================================================

pls2_final <- plsr(Y_pls2 ~ X_pls2, 
                   ncomp = best_n_pls2, 
                   scale = TRUE)

var_X <- cumsum(pls2_final$Xvar / pls2_final$Xtotvar * 100)
var_Y <- drop(R2(pls2_final, estimate = "train")$val)[-1]

cat("\nVariance X:", round(var_X[best_n_pls2], 1), "%\n")
cat("R² Y moyen:", round(mean(var_Y[best_n_pls2]), 4), "\n")

# =========================================================
# SCORES ET DÉPARTEMENTS EXTRÊMES
# =========================================================

cat("\n=== DÉPARTEMENTS EXTRÊMES ===\n")

scores_pls <- scores(pls2_final)

# PLS1
pls1 <- scores_pls[, 1]
idx_min <- order(pls1)[1:5]
idx_max <- order(pls1, decreasing = TRUE)[1:5]

cat("\nPLS1 minimum (5 départements):\n")
for(i in idx_min) {
  cat(sprintf("  %3d: %.3f\n", i, pls1[i]))
}

cat("\nPLS1 maximum (5 départements):\n")
for(i in idx_max) {
  cat(sprintf("  %3d: %.3f\n", i, pls1[i]))
}

# PLS2
if(ncol(scores_pls) >= 2) {
  pls2 <- scores_pls[, 2]
  idx_min2 <- order(pls2)[1:5]
  idx_max2 <- order(pls2, decreasing = TRUE)[1:5]
  
  cat("\nPLS2 minimum (5 départements):\n")
  for(i in idx_min2) {
    cat(sprintf("  %3d: %.3f\n", i, pls2[i]))
  }
  
  cat("\nPLS2 maximum (5 départements):\n")
  for(i in idx_max2) {
    cat(sprintf("  %3d: %.3f\n", i, pls2[i]))
  }
}

# =========================================================
# CERCLE DE CORRÉLATION
# =========================================================

cat("\n=== VISUALISATION ===\n")

loadings <- cor(X_pls2, scores_pls[, 1:2])

plot(loadings, xlim = c(-1, 1), ylim = c(-1, 1), asp = 1,
     main = "Cercle de corrélation X / PLS", pch = 19, col = "steelblue")
text(loadings, labels = colnames(X_pls2), pos = 3, cex = 0.7)
symbols(0, 0, circles = 1, add = TRUE, inches = FALSE, fg = "gray")
abline(h = 0, v = 0, lty = 2)
arrows(0, 0, loadings[, 1]*0.9, loadings[, 2]*0.9, 
       length = 0.1, col = "steelblue", lwd = 1.5)

# =========================================================
# PERFORMANCE PAR VARIABLE - CORRECTION
# =========================================================

cat("\n=== PERFORMANCE PAR VARIABLE ===\n")

# R² ajustement
R2_fit <- drop(R2(pls2_final, estimate = "train")$val)[-1]

# Q² prédiction
Q2_cv <- Q2_array[1, , best_n_pls2 + 1]

# CORRECTION : Prédictions en array 3D -> matrice 2D
Y_pred_array <- predict(pls2_final, ncomp = best_n_pls2)

# Si c'est un array 3D [n, q, 1], on convertit en matrice [n, q]
if(length(dim(Y_pred_array)) == 3) {
  Y_pred <- Y_pred_array[, , 1]
} else {
  Y_pred <- Y_pred_array
}

# Vérification des dimensions
cat("Dimensions Y_pls2:", dim(Y_pls2), "\n")
cat("Dimensions Y_pred:", dim(Y_pred), "\n")

# RMSE et corrélation
rmse <- sqrt(colMeans((Y_pls2 - Y_pred)^2))
cor_obs_pred <- diag(cor(Y_pls2, Y_pred))

results <- data.frame(
  Variable = colnames(Y_pls2),
  R2 = round(R2_fit[best_n_pls2], 4),
  Q2 = round(Q2_cv, 4),
  RMSE = round(rmse, 4),
  Cor = round(cor_obs_pred, 4)
)

print(results, row.names = FALSE)

# =========================================================
# VIP
# =========================================================

cat("\n=== VIP ===\n")

vip <- matrix(0, nrow = p, ncol = q)
for(j in 1:q) {
  pls_temp <- plsr(Y_pls2[, j] ~ X_pls2, ncomp = best_n_pls2, scale = TRUE)
  W <- pls_temp$loading.weights[, 1:best_n_pls2, drop = FALSE]
  ss <- colSums(pls_temp$scores[, 1:best_n_pls2]^2)
  vip[, j] <- sqrt(p * rowSums((W^2) * matrix(rep(ss, each = p), nrow = p)) / sum(ss))
}

vip_mean <- rowMeans(vip)
vip_df <- data.frame(Variable = colnames(X_pls2), VIP = round(vip_mean, 3))
vip_df <- vip_df[order(-vip_df$VIP), ]

cat("Top 5 VIP:\n")
print(head(vip_df, 5))

# =========================================================
# DÉPARTEMENTS ATYPIQUES
# =========================================================

cat("\n=== ATYPIQUES ===\n")

residus <- Y_pls2 - Y_pred
res_norm <- sqrt(rowSums(residus^2))
res_std <- scale(res_norm)

atyp <- which(abs(res_std) > 2)
if(length(atyp) > 0) {
  cat("Départements atypiques (|résidu std| > 2):", length(atyp), "\n")
  atyp_df <- data.frame(
    Index = atyp,
    Residu_Std = round(res_std[atyp], 2),
    Residu_Norm = round(res_norm[atyp], 4)
  )
  print(atyp_df)
} else {
  cat("Pas de départements atypiques détectés\n")
}

# Top 5 meilleures et pires
cat("\n5 meilleures prédictions:\n")
best_idx <- order(res_norm)[1:5]
best_df <- data.frame(
  Index = best_idx,
  Residu_Norm = round(res_norm[best_idx], 4)
)
print(best_df, row.names = FALSE)

cat("\n5 pires prédictions:\n")
worst_idx <- order(res_norm, decreasing = TRUE)[1:5]
worst_df <- data.frame(
  Index = worst_idx,
  Residu_Norm = round(res_norm[worst_idx], 4)
)
print(worst_df, row.names = FALSE)

cat("\n=== FIN PLS2 ===\n")

# =========================================================
# 11. DÉTECTION DES ATYPIQUES (Hotelling T² et résidus)
# =========================================================

library(pls)
library(ggplot2)
library(ggrepel)

cat("\n=== DÉTECTION DES ATYPIQUES PLS2 ===\n")

# =========================================================
# PRÉPARATION
# =========================================================

X_pls2 <- data_all_xy$X
Y_pls2 <- data_all_xy$Y

if("p_blanc1" %in% colnames(X_pls2)) {
  X_pls2 <- X_pls2[, colnames(X_pls2) != "p_blanc1", drop = FALSE]
}

# Noms des départements / territoires
dept_names <- data$Nom_dept

# Vérification cohérence
if(length(dept_names) != nrow(X_pls2)) {
  stop("Le nombre de noms de départements ne correspond pas au nombre de lignes de X_pls2.")
}

# =========================================================
# MODÈLE PLS2 POUR L’ANALYSE DES ATYPIQUES
# =========================================================

n_comp_T2 <- 2
pls2_T2 <- plsr(Y_pls2 ~ X_pls2, ncomp = n_comp_T2, scale = TRUE)

scores_T2 <- scores(pls2_T2)[, 1:n_comp_T2, drop = FALSE]

# =========================================================
# CALCUL T² DE HOTELLING
# =========================================================

scores_scaled <- scale(scores_T2)
T2 <- rowSums(scores_scaled^2)

seuil_T2 <- qchisq(0.95, df = n_comp_T2)
cat("Seuil T² (95%, df =", n_comp_T2, ") :", round(seuil_T2, 2), "\n")

# =========================================================
# RÉSIDUS PLS2
# =========================================================

Y_pred_T2 <- predict(pls2_T2, ncomp = n_comp_T2)
if(length(dim(Y_pred_T2)) == 3) {
  Y_pred_T2 <- Y_pred_T2[, , 1]
}

residus_T2 <- Y_pls2 - Y_pred_T2
resid_norm_T2 <- sqrt(rowSums(residus_T2^2))

seuil_resid <- as.numeric(quantile(resid_norm_T2, 0.95, na.rm = TRUE))

cat("Norme moyenne des résidus :", round(mean(resid_norm_T2), 4), "\n")
cat("Seuil résidus (95%) :", round(seuil_resid, 4), "\n")

# =========================================================
# TABLEAU DE SYNTHÈSE
# =========================================================

results_atyp <- data.frame(
  Index = 1:nrow(X_pls2),
  Departement = dept_names,
  T2 = as.numeric(T2),
  Residu_Norm = as.numeric(resid_norm_T2),
  Atypique_T2 = T2 > seuil_T2,
  Atypique_Resid = resid_norm_T2 > seuil_resid,
  stringsAsFactors = FALSE
)

results_atyp$Type <- "Normal"
results_atyp$Type[results_atyp$Atypique_T2 & !results_atyp$Atypique_Resid] <- "Atypique T²"
results_atyp$Type[!results_atyp$Atypique_T2 & results_atyp$Atypique_Resid] <- "Atypique résidu"
results_atyp$Type[results_atyp$Atypique_T2 & results_atyp$Atypique_Resid] <- "Atypique T² + résidu"

# =========================================================
# AFFICHAGE TEXTE
# =========================================================

cat("\n--- Départements atypiques selon T² ---\n")
atyp_T2_df <- results_atyp[results_atyp$Atypique_T2, c("Departement", "Index", "T2", "Residu_Norm")]
atyp_T2_df <- atyp_T2_df[order(-atyp_T2_df$T2), ]
print(atyp_T2_df, row.names = FALSE)

cat("\n--- Départements atypiques selon les résidus ---\n")
atyp_resid_df <- results_atyp[results_atyp$Atypique_Resid, c("Departement", "Index", "T2", "Residu_Norm")]
atyp_resid_df <- atyp_resid_df[order(-atyp_resid_df$Residu_Norm), ]
print(atyp_resid_df, row.names = FALSE)

cat("\n--- Départements atypiques sur les deux critères ---\n")
double_atyp_df <- results_atyp[
  results_atyp$Atypique_T2 & results_atyp$Atypique_Resid,
  c("Departement", "Index", "T2", "Residu_Norm")
]
double_atyp_df <- double_atyp_df[order(-double_atyp_df$T2), ]
print(double_atyp_df, row.names = FALSE)

# =========================================================
# PARAMÈTRES GRAPHIQUES COMMUNS
# =========================================================

# Espace supplémentaire à droite pour éviter que les labels soient coupés
x_max_extra <- max(results_atyp$Index) + 12

theme_atyp <- theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold"),
    plot.margin = margin(t = 15, r = 90, b = 15, l = 15)
  )

# =========================================================
# GRAPHE 1 : DISTANCE DE HOTELLING
# =========================================================

plot_T2 <- ggplot(results_atyp, aes(x = Index, y = T2, color = Atypique_T2)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_hline(
    yintercept = seuil_T2,
    linetype = "dashed",
    color = "red",
    linewidth = 1
  ) +
  ggrepel::geom_text_repel(
    data = subset(results_atyp, Atypique_T2),
    aes(label = Departement),
    color = "red",
    size = 5,
    direction = "y",
    hjust = 0,
    nudge_x = 2,
    box.padding = 0.35,
    point.padding = 0.25,
    segment.color = "grey50",
    segment.size = 0.4,
    min.segment.length = 0,
    max.overlaps = Inf,
    seed = 123
  ) +
  scale_color_manual(values = c("FALSE" = "steelblue", "TRUE" = "red")) +
  scale_x_continuous(
    limits = c(1, x_max_extra),
    breaks = seq(0, max(results_atyp$Index), by = 10),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    title = "Distance de Hotelling par département",
    x = "Indice des départements",
    y = "Distance de Hotelling (T²)",
    color = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_atyp

# =========================================================
# GRAPHE 2 : RÉSIDUS DE LA PLS2
# =========================================================

plot_resid <- ggplot(results_atyp, aes(x = Index, y = Residu_Norm, color = Atypique_Resid)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_hline(
    yintercept = seuil_resid,
    linetype = "dashed",
    color = "red",
    linewidth = 1
  ) +
  ggrepel::geom_text_repel(
    data = subset(results_atyp, Atypique_Resid),
    aes(label = Departement),
    color = "red",
    size = 5,
    direction = "y",
    hjust = 0,
    nudge_x = 2,
    box.padding = 0.35,
    point.padding = 0.25,
    segment.color = "grey50",
    segment.size = 0.4,
    min.segment.length = 0,
    max.overlaps = Inf,
    seed = 123
  ) +
  scale_color_manual(values = c("FALSE" = "steelblue", "TRUE" = "red")) +
  scale_x_continuous(
    limits = c(1, x_max_extra),
    breaks = seq(0, max(results_atyp$Index), by = 10),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    title = "Résidus de la PLS2 par département",
    x = "Indice des départements",
    y = "Distance euclidienne des résidus",
    color = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_atyp

print(plot_T2)
print(plot_resid)

# =========================================================
# PHRASES D’INTERPRÉTATION AUTOMATIQUES
# =========================================================

top_T2 <- head(results_atyp[order(-results_atyp$T2), "Departement"], 3)
top_resid <- head(results_atyp[order(-results_atyp$Residu_Norm), "Departement"], 3)

cat("\n=== INTERPRÉTATION SYNTHÉTIQUE ===\n")

cat(
  paste0(
    "La distance de Hotelling met en évidence des territoires au profil électoral particulièrement atypique dans l’espace latent de la PLS2, notamment ",
    paste(top_T2, collapse = ", "),
    ". "
  )
)

cat(
  paste0(
    "Les résidus les plus élevés concernent surtout ",
    paste(top_resid, collapse = ", "),
    ", ce qui indique que la reconstitution du second tour par le modèle est moins bonne pour ces territoires. "
  )
)

cat(
  "La confrontation des deux indicateurs permet ainsi de distinguer les départements simplement atypiques de ceux qui sont à la fois atypiques et mal expliqués par la relation entre premier et second tour.\n"
)

cat("\n=== FIN DÉTECTION ATYPIQUES ===\n")