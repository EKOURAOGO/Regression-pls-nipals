# Régression PLS-NIPALS en R

> Régression PLS (Partial Least Squares) via l'algorithme NIPALS  
> Application sur les résultats de l'élection présidentielle 2022  
> Master 2 IMSD · Paris-Saclay · 2025-2026

---

## Contexte

Ce projet applique la **régression PLS-NIPALS** à des données socio-économiques départementales
pour expliquer les résultats du premier tour de l'**élection présidentielle française 2022**.

Problématique centrale : comment les caractéristiques socio-économiques des départements
(chômage, revenus, densité, niveau d'éducation...) expliquent-elles les scores des candidats ?

---

## Structure du projet

```
Regression-pls-nipals/
├── PLS1_2.R                   # PLS1 et PLS2 — algorithme NIPALS, composantes, R²
├── PLS_3_4.R                  # PLS3 et PLS4 — validation croisée, scores optimaux
├── PROJET_PLS_EMMANUEL.pdf    # Rapport complet
└── README.md
```

---

## Méthodes implémentées

### PLS1 — Une seule variable réponse
- Algorithme NIPALS implémenté manuellement en R
- Déflation des matrices X et Y
- Sélection du nombre de composantes par **validation croisée**
- R² cumulé, coefficients de régression, importance des variables (VIP)

### PLS2 — Plusieurs variables réponses simultanées
- Extension PLS1 à plusieurs candidats simultanément
- Cercle des corrélations (variables explicatives vs scores des candidats)
- Biplot composantes / départements

### Analyse comparative
- Comparaison PLS vs régression linéaire classique
- Gestion de la **multicolinéarité** entre variables socio-économiques
- Stabilité des composantes par bootstrap

---

## Dataset — Élection présidentielle 2022

| Variable | Type | Description |
|----------|------|-------------|
| Score candidats | Réponse (Y) | % voix par département au 1er tour |
| Chômage | Explicative (X) | Taux de chômage départemental |
| Revenu médian | Explicative (X) | Revenu médian par UC |
| Densité | Explicative (X) | Densité de population |
| Niveau éducation | Explicative (X) | % diplômés supérieur |
| ... | ... | Variables socio-économiques départementales |

---

## Installation

```r
install.packages(c("pls", "ggplot2", "corrplot", "dplyr", "tidyr", "MASS"))
source("PLS1_2.R")
source("PLS_3_4.R")
```

---

## Stack technique

![R](https://img.shields.io/badge/R-276DC3?style=flat-square&logo=r&logoColor=white)
![pls](https://img.shields.io/badge/pls-NIPALS-blue?style=flat-square)
![ggplot2](https://img.shields.io/badge/ggplot2-visualization-red?style=flat-square)
![tidyverse](https://img.shields.io/badge/tidyverse-276DC3?style=flat-square)

---

## Auteur

**KOURAOGO Emmanuel** — M2 IMSD · Paris-Saclay  
Data Scientist & Data Analyst · DREES

[![GitHub](https://img.shields.io/badge/GitHub-EKOURAOGO-181717?style=flat-square&logo=github)](https://github.com/EKOURAOGO)
