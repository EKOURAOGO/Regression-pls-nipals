# Dynamiques Électorales & Régression PLS - Présidentielle 2022

> Analyse départementale par régression PLS (Partial Least Squares)  
> Anticipation du second tour à partir des résultats du premier tour  


---

## Contexte

L'élection présidentielle française de 2022 a mis en évidence de fortes **disparités territoriales** dans les comportements électoraux. Ce projet analyse dans quelle mesure les résultats du **premier tour** permettent d'expliquer et d'anticiper ceux du **second tour** à l'échelle départementale.

**Question centrale :** Le premier tour structure-t-il fortement le second, ou observe-t-on des recompositions plus complexes selon les territoires ?

En raison de la nature **compositionnelle** des données électorales (multicolinéarité structurelle), ce projet compare trois approches : **MCO**, **PCR** et **PLS**, avant de conduire une modélisation multivariée par **PLS2**.

> Encadrant : Christian DERQUENNE

---

## Données

| Caractéristique | Valeur |
|---|---|
| Unité d'observation | Département |
| Nombre d'unités | 107 (métropole + DOM-TOM + Français de l'étranger) |
| Candidats T1 | 12 (Macron, Le Pen, Mélenchon, Zemmour, Pécresse...) |
| Variables réponses T2 | p_abst2, p_blanc2, p_macron2, p_lepen2 |

---

## Structure du projet

```
Regression-pls-nipals/
├── PLS1_2.R                   # PLS1 + PLS2, ACP, MCO, PCR, VIP, validation croisée LOO
├── PLS_3_4.R                  # PLS2 multivariée, Hotelling, départements atypiques
├── PROJET_PLS_EMMANUEL.pdf    # Rapport complet (25 pages)
└── README.md
```

---

## Méthodes implémentées

| Méthode | Principe |
|---------|----------|
| MCO | Régression linéaire standard |
| PCR | Composantes principales de X puis régression |
| PLS1 | Composantes maximisant cov(X, y) - validation croisée LOO |
| PLS2 | Réponse simultanée sur 4 variables T2 |
| ACP | T1, T2 et conjointe T1+T2 |
| VIF | Détection multicolinéarité structurelle |
| Hotelling | Détection des départements atypiques |

---

## Principaux résultats

### Multicolinéarité structurelle

| Variable | VIF |
|---|---|
| p_abst1 | **1950.62** |
| p_LE_PEN_1 | 465.33 |
| p_MACRON_1 | 316.54 |

→ Justifie pleinement le recours à PLS.

### Comparaison MCO / PCR / PLS1

| Variable | R² MCO | Q² PCR | Q² PLS1 | Optimal |
|---|---|---|---|---|
| p_abst2 | 0.984 | 0.980 | 0.979 | PCR ≈ PLS |
| p_blanc2 | 0.946 | 0.865 | 0.828 | **PCR** |
| p_macron2 | 0.981 | 0.952 | **0.963** | **PLS1** |
| p_lepen2 | 0.941 | 0.907 | **0.909** | **PLS1** |

### PLS2 multivariée (Q²)

| Variable | Q² |
|---|---|
| p_abst2 | 0.9799 |
| p_macron2 | 0.9574 |
| p_lepen2 | 0.9116 |
| p_blanc2 | 0.7995 |

**Composante 1** → axe de mobilisation électorale (abstention vs vote)  
**Composante 2** → axe de différenciation politique (Le Pen vs Mélenchon/Macron/Jadot)

---

## Installation

```r
install.packages(c("pls", "FactoMineR", "factoextra", "corrplot",
                   "car", "ggplot2", "dplyr", "tidyr", "ggrepel"))

source("PLS1_2.R")
source("PLS_3_4.R")
```

---

## Stack technique

![R](https://img.shields.io/badge/R-276DC3?style=flat-square&logo=r&logoColor=white)
![pls](https://img.shields.io/badge/pls-PLS1%20·%20PLS2%20·%20PCR-blue?style=flat-square)
![FactoMineR](https://img.shields.io/badge/FactoMineR-ACP-orange?style=flat-square)
![ggplot2](https://img.shields.io/badge/ggplot2-visualization-red?style=flat-square)
![car](https://img.shields.io/badge/car-VIF-green?style=flat-square)

---

## Auteur

**KOURAOGO Emmanuel** 
Data Scientist & Data Analyst 

[![GitHub](https://img.shields.io/badge/GitHub-EKOURAOGO-181717?style=flat-square&logo=github)](https://github.com/EKOURAOGO)

*Encadrant : Christian DERQUENNE · 2025-2026*
