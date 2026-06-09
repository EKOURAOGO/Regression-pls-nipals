# =========================================================
# PROJET PLS – PRÉSIDENTIELLE 2022
# Master 2 IMSD – 2026
# =========================================================

packages <- c("tidyverse","corrplot","FactoMineR","factoextra",
              "pls","plsdepot","car","knitr","reshape2","gridExtra")

# Vérification version plsdepot
ver <- packageVersion("plsdepot")
cat("Version plsdepot:", as.character(ver), "\n")
if(ver < "0.3.0") {
  cat("⚠️  Version obsolète - risque de bugs Q²/R²\n")
}

for(p in packages){
  if(!require(p, character.only=TRUE)){
    install.packages(p)
    library(p, character.only=TRUE)
  }
}

theme_set(theme_minimal())
set.seed(123)

# =========================================================
# 1. IMPORTATION
# =========================================================

chemin_fichier <- "C:/Users/Pc/OneDrive/Bureau/Dossier/IMSD/PLS/Projet PLS/PLS/Résultats présidentielles 2022.txt"

# ⚠️ Encodage corrigé (fichier probablement en latin1)
data <- read.delim(chemin_fichier,
                   sep = "\t",
                   header = TRUE,
                   encoding = "latin1",
                   stringsAsFactors = FALSE)

# Correction encodage colonnes texte
data$Nom_dept <- iconv(data$Nom_dept, from = "latin1", to = "UTF-8")

# =========================================================
# 2. RENOMMAGE ET NETTOYAGE
# =========================================================

cat("\nNoms originaux des colonnes:\n")
print(names(data))

# Renommage adaptatif sécurisé
names(data) <- gsub("^MACRON$", "MACRON_1", names(data))
names(data) <- gsub("^LE_PEN$|^LE\\.PEN$|^LE PEN$", "LE_PEN_1", names(data))

cat("\nNoms après renommage:\n")
print(names(data))

cat("\nDimensions :", nrow(data), "départements x", 
    ncol(data), "variables\n")

# =========================================================
# 2. CONTRÔLE SIMPLE DES DONNÉES (Conforme TP)
# =========================================================

cat("\n=== CONTRÔLE DES DONNÉES ===\n")

# Nombre d'individus
cat("Nombre de départements :", nrow(data), "\n")

# Doublons éventuels
if(anyDuplicated(data$Num_dept) > 0){
  cat("Doublons détectés\n")
} else {
  cat("Aucun doublon détecté\n")
}

# Valeurs manquantes
na_count <- sum(is.na(data))
cat("Nombre total de valeurs manquantes :", na_count, "\n")

# Cohérence Tour 1
check_t1 <- data$Abstentions_1 + data$Votants_1 - data$Inscrits_1
cat("Écart max cohérence T1 :", max(abs(check_t1)), "\n")

# Cohérence Tour 2
check_t2 <- data$Abstentions_2 + data$Votants_2 - data$Inscrits_2
cat("Écart max cohérence T2 :", max(abs(check_t2)), "\n")

# =========================================================
# 3. TRANSFORMATION EN PROPORTIONS
# =========================================================

candidats_t1 <- c("ARTHAUD","ROUSSEL","MACRON_1","LASSALLE",
                  "LE_PEN_1","ZEMMOUR","MELENCHON","HIDALGO",
                  "JADOT","PECRESSE","POUTOU","DUPONT_AIGNAN")

# Vérification colonnes
candidats_existants <- candidats_t1[candidats_t1 %in% names(data)]

if(length(candidats_existants) < length(candidats_t1)){
  warning("Candidats manquants : ",
          paste(setdiff(candidats_t1, candidats_existants), collapse=", "))
}

# Proportions T1 (rapport aux inscrits)
for(c in candidats_existants){
  data[[paste0("p_", c)]] <- data[[c]] / data$Inscrits_1
}

data <- data %>%
  mutate(
    p_abst1  = Abstentions_1 / Inscrits_1,
    p_blanc1 = (Blancs_1 + Nuls_1) / Inscrits_1
  )

# Vérification somme T1 = 1
somme_t1 <- data$p_abst1 + data$p_blanc1 +
  rowSums(data[paste0("p_", candidats_existants)])

cat("\nVérification T1 (somme ≈ 1):\n")
print(summary(somme_t1))


# =========================================================
# Proportions T2
# =========================================================

data <- data %>%
  mutate(
    p_abst2   = Abstentions_2 / Inscrits_2,
    p_blanc2  = (Blancs_2 + Nuls_2) / Inscrits_2,
    p_macron2 = MACRON_2 / Inscrits_2,
    p_lepen2  = LE_PEN_2 / Inscrits_2
  )

# Vérification somme T2 = 1
somme_t2 <- data$p_abst2 + data$p_blanc2 +
  data$p_macron2 + data$p_lepen2

cat("\nVérification T2 (somme ≈ 1):\n")
print(summary(somme_t2))


# =========================================================
# Filtrage Métropole
# =========================================================

data_metropole <- data %>%
  filter(Num_dept %in% c(1:95, 201, 202))

cat("\nDépartements métropole :", nrow(data_metropole), "\n")



# =========================================================
# 4. MATRICES X ET Y
# =========================================================

build_XY <- function(data_in){
  
  X_loc <- data_in %>%
    select(starts_with("p_")) %>%
    select(-p_abst2, -p_blanc2, -p_macron2, -p_lepen2)
  
  Y_loc <- data_in %>%
    select(p_abst2, p_blanc2, p_macron2, p_lepen2)
  
  list(X = as.matrix(X_loc),
       Y = as.matrix(Y_loc))
}

data_all_xy   <- build_XY(data)
data_metro_xy <- build_XY(data_metropole)

cat("\nDimensions X (all):", dim(data_all_xy$X))
cat("\nDimensions X (metro):", dim(data_metro_xy$X), "\n")

# =========================================================
# 5. STATISTIQUES DESCRIPTIVES
# =========================================================

cat("\n=== STATISTIQUES DESCRIPTIVES ===\n")

# Définition explicite (sécurité)
X <- data_all_xy$X
Y <- data_all_xy$Y

# ---------------------------------------------------------
# 1) Premier tour – PROPORTIONS
# ---------------------------------------------------------

cat("\n--- Premier Tour (proportions des inscrits) ---\n")

desc_t1 <- data.frame(
  Variable = colnames(X),
  Moyenne = round(colMeans(X), 4),
  Ecart_type = round(apply(X, 2, sd), 4),
  Min = round(apply(X, 2, min), 4),
  Max = round(apply(X, 2, max), 4)
)

desc_t1$CV_pct <- round(desc_t1$Ecart_type / desc_t1$Moyenne * 100, 1)

print(knitr::kable(desc_t1, digits = 4))


# ---------------------------------------------------------
# 2) Second tour – PROPORTIONS
# ---------------------------------------------------------

cat("\n--- Second Tour (proportions des inscrits) ---\n")

desc_t2 <- data.frame(
  Variable = colnames(Y),
  Moyenne = round(colMeans(Y), 4),
  Ecart_type = round(apply(Y, 2, sd), 4),
  Min = round(apply(Y, 2, min), 4),
  Max = round(apply(Y, 2, max), 4)
)

desc_t2$CV_pct <- round(desc_t2$Ecart_type / desc_t2$Moyenne * 100, 1)

print(knitr::kable(desc_t2, digits = 4))


# ---------------------------------------------------------
# 3) Comparaison avec les NOMBRES BRUTS (exigé par le TP)
# ---------------------------------------------------------

cat("\n--- Comparaison avec les nombres bruts (Tour 1) ---\n")

vars_bruts_t1 <- c("Abstentions_1","Blancs_1","Nuls_1",
                   candidats_existants)

desc_brut_t1 <- data.frame(
  Variable = vars_bruts_t1,
  Moyenne = round(colMeans(data[vars_bruts_t1]), 0),
  Ecart_type = round(apply(data[vars_bruts_t1], 2, sd), 0),
  Min = round(apply(data[vars_bruts_t1], 2, min), 0),
  Max = round(apply(data[vars_bruts_t1], 2, max), 0)
)

print(knitr::kable(desc_brut_t1))

# =========================================================
# 5. COMPARAISON BRUTS VS PROPORTIONS 
# =========================================================

par(mfrow=c(2,2))

# ---------------------------
# MACRON
# ---------------------------

hist(data$MACRON_1,
     main="Macron T1 - Brut (voix)",
     xlab="Nombre de voix",
     col="lightblue",
     border="white")

abline(v=mean(data$MACRON_1), col="blue", lwd=2)

hist(data$p_MACRON_1,
     main="Macron T1 - Proportion",
     xlab="Proportion des inscrits",
     col="lightgreen",
     border="white")

abline(v=mean(data$p_MACRON_1), col="darkgreen", lwd=2)


# ---------------------------
# LE PEN
# ---------------------------

hist(data$LE_PEN_1,
     main="Le Pen T1 - Brut (voix)",
     xlab="Nombre de voix",
     col="lightcoral",
     border="white")

abline(v=mean(data$LE_PEN_1), col="red", lwd=2)

hist(data$p_LE_PEN_1,
     main="Le Pen T1 - Proportion",
     xlab="Proportion des inscrits",
     col="mistyrose",
     border="white")

abline(v=mean(data$p_LE_PEN_1), col="darkred", lwd=2)

par(mfrow=c(1,1))


cat("\nMoyennes comparées :\n")
cat("Macron brut :", round(mean(data$MACRON_1),0), "\n")
cat("Macron proportion :", round(mean(data$p_MACRON_1),4), "\n")
cat("Le Pen brut :", round(mean(data$LE_PEN_1),0), "\n")
cat("Le Pen proportion :", round(mean(data$p_LE_PEN_1),4), "\n")

cat("\n=== STATISTIQUES DESCRIPTIVES – MÉTROPOLE ===\n")

# Définition explicite (sécurité)
X_metro <- data_metro_xy$X
Y_metro <- data_metro_xy$Y

# ---------------------------------------------------------
# 1) Premier tour – PROPORTIONS
# ---------------------------------------------------------

cat("\n--- Premier Tour (proportions des inscrits) – MÉTROPOLE ---\n")

desc_t1_metro <- data.frame(
  Variable = colnames(X_metro),
  Moyenne = round(colMeans(X_metro), 4),
  Ecart_type = round(apply(X_metro, 2, sd), 4),
  Min = round(apply(X_metro, 2, min), 4),
  Max = round(apply(X_metro, 2, max), 4)
)

desc_t1_metro$CV_pct <- round(desc_t1_metro$Ecart_type / 
                                desc_t1_metro$Moyenne * 100, 1)

print(knitr::kable(desc_t1_metro, digits = 4))


# ---------------------------------------------------------
# 2) Second tour – PROPORTIONS
# ---------------------------------------------------------

cat("\n--- Second Tour (proportions des inscrits) – MÉTROPOLE ---\n")

desc_t2_metro <- data.frame(
  Variable = colnames(Y_metro),
  Moyenne = round(colMeans(Y_metro), 4),
  Ecart_type = round(apply(Y_metro, 2, sd), 4),
  Min = round(apply(Y_metro, 2, min), 4),
  Max = round(apply(Y_metro, 2, max), 4)
)

desc_t2_metro$CV_pct <- round(desc_t2_metro$Ecart_type / 
                                desc_t2_metro$Moyenne * 100, 1)

print(knitr::kable(desc_t2_metro, digits = 4))


# ---------------------------------------------------------
# 3) Comparaison avec les NOMBRES BRUTS (Tour 1) – MÉTROPOLE
# ---------------------------------------------------------

cat("\n--- Comparaison avec les nombres bruts (Tour 1) – MÉTROPOLE ---\n")

desc_brut_t1_metro <- data.frame(
  Variable = vars_bruts_t1,
  Moyenne = round(colMeans(data_metropole[vars_bruts_t1]), 0),
  Ecart_type = round(apply(data_metropole[vars_bruts_t1], 2, sd), 0),
  Min = round(apply(data_metropole[vars_bruts_t1], 2, min), 0),
  Max = round(apply(data_metropole[vars_bruts_t1], 2, max), 0)
)

print(knitr::kable(desc_brut_t1_metro))


# =========================================================
# 5 BIS. COMPARAISON BRUTS VS PROPORTIONS – MÉTROPOLE
# =========================================================

par(mfrow=c(2,2))

# ---------------------------
# MACRON – MÉTROPOLE
# ---------------------------

hist(data_metropole$MACRON_1,
     main="Macron T1 - Brut (voix) – Métropole",
     xlab="Nombre de voix",
     col="lightblue",
     border="white")

abline(v=mean(data_metropole$MACRON_1), col="blue", lwd=2)

hist(data_metropole$p_MACRON_1,
     main="Macron T1 - Proportion – Métropole",
     xlab="Proportion des inscrits",
     col="lightgreen",
     border="white")

abline(v=mean(data_metropole$p_MACRON_1), col="darkgreen", lwd=2)


# ---------------------------
# LE PEN – MÉTROPOLE
# ---------------------------

hist(data_metropole$LE_PEN_1,
     main="Le Pen T1 - Brut (voix) – Métropole",
     xlab="Nombre de voix",
     col="lightcoral",
     border="white")

abline(v=mean(data_metropole$LE_PEN_1), col="red", lwd=2)

hist(data_metropole$p_LE_PEN_1,
     main="Le Pen T1 - Proportion – Métropole",
     xlab="Proportion des inscrits",
     col="mistyrose",
     border="white")

abline(v=mean(data_metropole$p_LE_PEN_1), col="darkred", lwd=2)

par(mfrow=c(1,1))

cat("\nMoyennes comparées – MÉTROPOLE :\n")
cat("Macron brut :", round(mean(data_metropole$MACRON_1),0), "\n")
cat("Macron proportion :", round(mean(data_metropole$p_MACRON_1),4), "\n")
cat("Le Pen brut :", round(mean(data_metropole$LE_PEN_1),0), "\n")
cat("Le Pen proportion :", round(mean(data_metropole$p_LE_PEN_1),4), "\n")

# =========================================================
# 6. SCATTERPLOTS POLITIQUES
# =========================================================

# Corrélations
cor_t1 <- cor(data$p_MELENCHON, data$p_LE_PEN_1)
cor_t2 <- cor(data$p_macron2, data$p_lepen2)
cor_brut_t2 <- cor(data$MACRON_2, data$LE_PEN_2)

# ---------------------------------------------------------
# 1) Opposition idéologique T1 (proportions)
# ---------------------------------------------------------

p1 <- ggplot(data, aes(x = p_MELENCHON, y = p_LE_PEN_1)) +
  geom_point(color = "darkblue", size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
  coord_cartesian(xlim = c(0, max(data$p_MELENCHON)),
                  ylim = c(0, max(data$p_LE_PEN_1))) +
  labs(title = "Opposition gauche / droite – T1 (proportions)",
       subtitle = paste0("Corrélation r = ", round(cor_t1, 2)),
       x = "Mélenchon (proportion des inscrits)",
       y = "Le Pen (proportion des inscrits)") +
  theme_minimal(base_size = 12)

# ---------------------------------------------------------
# 2) Bipolarisation T2 (proportions)
# ---------------------------------------------------------

p2 <- ggplot(data, aes(x = p_macron2, y = p_lepen2)) +
  geom_point(color = "firebrick2", size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
  coord_cartesian(xlim = c(0, max(data$p_macron2)),
                  ylim = c(0, max(data$p_lepen2))) +
  labs(title = "Bipolarisation – T2 (proportions)",
       subtitle = paste0("Corrélation r = ", round(cor_t2, 2)),
       x = "Macron (proportion des inscrits)",
       y = "Le Pen (proportion des inscrits)") +
  theme_minimal(base_size = 12)

# ---------------------------------------------------------
# 3) Brut vs Brut – Effet taille T2
# ---------------------------------------------------------

p3 <- ggplot(data, aes(x = MACRON_2, y = LE_PEN_2)) +
  geom_point(color = "darkgreen", size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(title = "Macron vs Le Pen – T2 (nombres bruts)",
       subtitle = paste0("Corrélation r = ", round(cor_brut_t2, 2)),
       x = "Macron (voix brutes)",
       y = "Le Pen (voix brutes)") +
  theme_minimal(base_size = 12)

# ---------------------------------------------------------
# 4) Brut vs Proportion – Macron T1
# ---------------------------------------------------------

p4 <- ggplot(data, aes(x = MACRON_1, y = p_MACRON_1)) +
  geom_point(color = "steelblue", size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(title = "Brut vs Proportion – Macron T1",
       x = "Voix brutes",
       y = "Proportion des inscrits") +
  theme_minimal(base_size = 12)

# Affichage
gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)

# =========================================================
# 7. MULTICOLINÉARITÉ (VERSION MÉTHODOLOGIQUE CORRIGÉE)
# =========================================================

cat("\n=== ANALYSE DE LA MULTICOLINÉARITÉ ===\n")

# ---------------------------------------------------------
# 1) MATRICES AVEC ET SANS CONTRAINTE
# ---------------------------------------------------------

X_full <- data_all_xy$X

# Suppression d'une variable pour lever la contrainte de somme = 1
X_reduced <- X_full[, colnames(X_full) != "p_blanc1"]

cat("Dimensions X complet :", dim(X_full), "\n")
cat("Dimensions X réduit  :", dim(X_reduced), "\n")

# ---------------------------------------------------------
# 2) MATRICE DE CORRÉLATION (SANS CONTRAINTE)
# ---------------------------------------------------------

cor_matrix <- cor(X_reduced)

corrplot::corrplot(
  cor_matrix,
  method = "color",
  type = "upper",
  addCoef.col = "black",
  number.cex = 0.6,
  tl.cex = 0.7,
  tl.col = "black",
  tl.srt = 45,
  diag = FALSE,
  mar = c(0,0,2,0),
  title = "Matrice de corrélation – T1 (sans contrainte de somme)"
)

# Heatmap ggplot
cor_melt <- reshape2::melt(cor_matrix)

ggplot(cor_melt, aes(Var1, Var2, fill=value)) +
  geom_tile(color="white") +
  geom_text(aes(label=round(value,2)), size=2.8) +
  scale_fill_gradient2(low="navy", high="darkred", mid="white", 
                       midpoint=0, limits=c(-1,1)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, size=7),
        axis.text.y = element_text(size=7)) +
  labs(title="Heatmap des corrélations – T1 (sans contrainte)",
       x="", y="", fill="Corr")

# ---------------------------------------------------------
# 3) VIF (Variance Inflation Factor)
# ---------------------------------------------------------

# --- VIF avec contrainte (complet) ---
vif_full <- sapply(1:ncol(X_full), function(i){
  r2 <- summary(lm(X_full[,i] ~ X_full[,-i]))$r.squared
  ifelse(r2 < 0.999, 1/(1-r2), Inf)
})

# --- VIF sans contrainte ---
vif_reduced <- sapply(1:ncol(X_reduced), function(i){
  r2 <- summary(lm(X_reduced[,i] ~ X_reduced[,-i]))$r.squared
  1/(1-r2)
})

vif_df <- data.frame(
  Variable = colnames(X_reduced),
  VIF_SansContrainte = round(vif_reduced, 2)
)

print(knitr::kable(vif_df, digits = 2))

cat("\nVIF max (complet) :", max(vif_full, na.rm=TRUE))
cat("\nVIF max (sans contrainte) :", round(max(vif_reduced), 2), "\n")

cat("\nInterprétation VIF :\n")
if(max(vif_reduced) > 10){
  cat("→ Multicolinéarité sévère détectée (VIF > 10).\n")
} else if(max(vif_reduced) > 5){
  cat("→ Multicolinéarité modérée.\n")
} else {
  cat("→ Multicolinéarité faible.\n")
}

# ---------------------------------------------------------
# 4) CONDITION NUMBER (κ)
# ---------------------------------------------------------

X_scaled <- scale(X_reduced)
eigen_vals <- eigen(cor(X_scaled))$values
kappa <- sqrt(max(eigen_vals) / min(eigen_vals))

cat("\nNombre de condition (kappa) :", round(kappa, 2), "\n")

cat("Interprétation kappa :\n")
if(kappa < 10){
  cat("→ Pas de problème majeur.\n")
} else if(kappa < 30){
  cat("→ Multicolinéarité modérée.\n")
} else if(kappa < 100){
  cat("→ Multicolinéarité forte.\n")
} else {
  cat("→ Multicolinéarité très sévère.\n")
}

# ---------------------------------------------------------
# 5) CONCLUSION AUTOMATIQUE
# ---------------------------------------------------------

cat("\nConclusion :\n")
cat("Les variables du premier tour présentent une forte interdépendance.\n")
cat("Les indicateurs VIF et le nombre de condition confirment\n")
cat("la présence d'une multicolinéarité importante.\n")
cat("L'utilisation de la PLS est donc méthodologiquement justifiée.\n")


# =========================================================
# 8. ANALYSE EN COMPOSANTES PRINCIPALES (EXPLORATOIRE)
# =========================================================

cat("\n=== ANALYSE EN COMPOSANTES PRINCIPALES ===\n")

# IMPORTANT : on utilise la matrice sans contrainte
X_acp <- X_reduced
Y_acp <- data_all_xy$Y

# =========================================================
# 1) ACP PREMIER TOUR (X)
# =========================================================

acp_t1 <- PCA(X_acp, scale.unit = TRUE, graph = FALSE)

cat("\n--- ACP Premier Tour ---\n")
cat("Variance expliquée (%):\n")
print(round(acp_t1$eig[1:6, 2], 2))

# Scree plot
fviz_eig(acp_t1, addlabels = TRUE,
         title = "Scree plot – ACP Premier Tour")

# Cercle des corrélations
fviz_pca_var(acp_t1,
             axes = c(1, 2),
             repel = TRUE,
             col.var = "steelblue",
             title = "ACP T1 – Cercle des corrélations")

# Contributions principales
contrib_t1 <- acp_t1$var$contrib[,1:2]
print(round(contrib_t1[order(-contrib_t1[,1]),],2))

cat("\nInterprétation suggérée :\n")
cat("Dim1 semble structurer un axe idéologique majeur.\n")
cat("Dim2 pourrait capter une dynamique centre/périphérie ou participation.\n")

# =========================================================
# 2) ACP SECOND TOUR (Y)
# =========================================================

acp_t2 <- PCA(Y_acp, scale.unit = TRUE, graph = FALSE)

cat("\n--- ACP Second Tour ---\n")
cat("Variance expliquée (%):\n")
print(round(acp_t2$eig[, 2], 2))

fviz_eig(acp_t2, addlabels = TRUE,
         title = "Scree plot – ACP Second Tour")

fviz_pca_var(acp_t2,
             axes = c(1,2),
             repel = TRUE,
             title = "ACP T2 – Structuration du vote")

cat("\nInterprétation attendue :\n")
cat("Dim1 devrait opposer Macron et Le Pen.\n")
cat("Dim2 pourrait isoler l'abstention et les votes blancs.\n")

# =========================================================
# 3) CORRÉLATION ENTRE AXES T1 ET T2
# =========================================================

scores_t1 <- acp_t1$ind$coord[, 1:3]
scores_t2 <- acp_t2$ind$coord[, 1:2]

cor_axes <- cor(scores_t1, scores_t2)

cat("\nCorrélations entre axes ACP T1 et T2:\n")
print(round(cor_axes, 3))

cat("\nInterprétation :\n")
cat("Une forte corrélation entre Dim1_T1 et Dim1_T2\n")
cat("indique que la structuration idéologique du premier tour\n")
cat("se retrouve au second tour.\n")

# =========================================================
# 3) ACP CONJOINTE (T1 + T2)
# =========================================================

cat("\n--- ACP Conjointe (T1 + T2) ---\n")

# On combine la matrice X réduite (sans contrainte)
# et la matrice Y
XY_joint <- cbind(X_acp, Y_acp)

acp_joint <- PCA(XY_joint,
                 scale.unit = TRUE,
                 graph = FALSE)

cat("Variance expliquée (%):\n")
print(round(acp_joint$eig[1:6, 2], 2))

# Scree plot
fviz_eig(acp_joint,
         addlabels = TRUE,
         title = "Scree plot – ACP conjointe (T1 + T2)")

# Cercle des corrélations
fviz_pca_var(acp_joint,
             axes = c(1,2),
             repel = TRUE,
             title = "ACP conjointe – Structure globale T1 + T2")

# # =========================================================
# graphe
# =========================================================
cat("\n=== CARTOGRAPHIE BRUTE vs STANDARDISÉE ===\n")

# -----------------------------------------
# 1) CARTOGRAPHIE AVANT STANDARDISATION
# -----------------------------------------

par(mfrow=c(1,1))
image(t(X_full),
      col=heat.colors(60),
      axes=FALSE,
      main="Cartographie brute – T1 (proportions)")
axis(1, at=seq(0,1,length.out=ncol(X_full)),
     labels=colnames(X_full), las=2, cex.axis=0.6)
axis(2, at=seq(0,1,length.out=nrow(X_full)),
     labels=1:nrow(X_full), las=1, cex.axis=0.6)

# -----------------------------------------
# 2) CARTOGRAPHIE APRÈS STANDARDISATION
# -----------------------------------------

image(t(scale(X_full)),
      col=heat.colors(60),
      axes=FALSE,
      main="Cartographie standardisée – T1")
axis(1, at=seq(0,1,length.out=ncol(X_full)),
     labels=colnames(X_full), las=2, cex.axis=0.6)
axis(2, at=seq(0,1,length.out=nrow(X_full)),
     labels=1:nrow(X_full), las=1, cex.axis=0.6)

par(mfrow=c(1,1))

#MATRICE COMPLÈTE DES SCATTERPLOTS
cat("\n=== MATRICE SCATTERPLOTS T1 ===\n")

pairs(X_reduced,
      main="Matrice des scatterplots – T1",
      pch=19,
      col=rgb(0,0,0,0.5))

#MATRICE AVEC CORRÉLATIONS ET P-VALUES
cat("\n=== MATRICE CORRÉLATION + P-VALUE ===\n")

panel.cor <- function(x, y, digits=2, prefix="", cex.cor){
  usr <- par("usr"); on.exit(par(usr))
  par(usr=c(0,1,0,1))
  r <- cor(x, y)
  test <- cor.test(x, y)
  txt <- paste0("r=", format(c(r, 0.123456789), digits=digits)[1])
  pval <- paste0("p=", format.pval(test$p.value, digits=2))
  text(0.5, 0.6, txt, cex=1.2)
  text(0.5, 0.4, pval, cex=0.9)
}

pairs(X_reduced,
      upper.panel=panel.cor,
      lower.panel=panel.smooth,
      main="Corrélations linéaires et p-values – T1")


cat("\n=== HEATMAP GLOBALE T1 + T2 ===\n")

XY_joint_scaled <- scale(cbind(X_acp, Y_acp))

heatmap(XY_joint_scaled,
        scale="none",
        col=colorRampPalette(c("navy","white","darkred"))(60),
        margins=c(6,6),
        main="Cartographie globale standardisée – T1 + T2")


library(GGally)

cat("\n=== GGPairs – MATRICE MODERNE ===\n")

GGally::ggpairs(
  as.data.frame(X_reduced),
  upper = list(continuous = wrap("cor", size=3)),
  lower = list(continuous = "points"),
  diag  = list(continuous = "densityDiag")
)
