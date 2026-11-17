# Indexation automatique (recherche de maille) multi-méthodes pour GSAS-II

*[English version: README.md](README.md)*

Extension pour **GSAS-II** qui cherche la **maille d'un diagramme de poudre** en lançant plusieurs méthodes d'indexation côte à côte. Elle compare leurs solutions avec les mêmes figures de mérite, puis envoie les mailles que vous choisissez dans la *Unit Cells List* de votre projet.

C'est la première étape d'une chaîne simple :
**Indexation** (cette extension) → [**gsasii-cellsearch**](https://github.com/yrousselin/gsasii-cellsearch) (cette maille est-elle déjà connue ?) → [**gsasii-phasematch**](https://github.com/yrousselin/gsasii-phasematch) (quelles phases contient l'échantillon ?).

## Ce que fait l'extension

- Ajoute à GSAS-II le menu **Export → Entire project as → (user) Auto-index (find unit cell)**.
- Lit les pics de l'histogramme choisi : l'*Index Peak List* si elle est remplie, sinon la *Peak List*.
- Lance **en parallèle** les méthodes que vous cochez :

| Méthode | Ce que c'est | Systèmes cristallins |
|---|---|---|
| **Conograph** | programme d'indexation externe (conseillé) | tous, y compris triclinique |
| **Dichotomie** | réécriture en Python du principe de DICVOL (dichotomie par tranches de volume) | cubique à monoclinique |
| **Essai-erreur** | réécriture en Python du principe de TREOR (essais sur les raies de bas angle) | cubique à monoclinique |
| **GSAS-II natif** | le moteur d'indexation intégré à GSAS-II | tous (lent en triclinique) |

- Affine chaque maille candidate sur les pics, puis recalcule **M20 et F20 de la même façon pour toutes les méthodes**. Leurs solutions peuvent donc être comparées équitablement.
- Affiche un compte rendu avec la comparaison des méthodes et un classement général. Vous cochez les mailles à garder.
- Écrit les mailles cochées dans la **Unit Cells List** du projet ouvert. Rien n'est enregistré sur le disque tant que vous ne sauvegardez pas le projet.
- Le calcul tourne dans un processus séparé : GSAS-II reste utilisable, et un bouton **Annuler** arrête la recherche.

Les moteurs « dichotomie » et « essai-erreur » sont des réécritures du principe de DICVOL et de TREOR, **pas les programmes originaux**.

## Prérequis

| Quoi | Pourquoi | Où le trouver |
|---|---|---|
| **GSAS-II** | accueille l'extension ; fournit numpy, scipy et wxPython | [GSAS-II](https://advancedphotonsource.github.io/GSAS-II-tutorials/install.html) |
| **Conograph** *(conseillé)* | la méthode la plus robuste, la seule avec GSAS-II pour les mailles tricliniques | [téléchargements Z-Code](https://z-code-software.com/downloads/) |
| Une liste de pics | l'entrée de toutes les méthodes | *Peak List* ou *Index Peak List* dans GSAS-II |

Sans Conograph, les trois autres méthodes fonctionnent quand même.

**Conseil :** des pics ajustés dans GSAS-II (*Peak List → Peak fitting*) donnent des mailles plus précises que des pics simplement pointés.

## Installation

### Étape 1 : copier les fichiers dans GSAS-II

**macOS (et Linux)**
1. Téléchargez ce dépôt : bouton vert **Code** → **Download ZIP**, puis double-cliquez sur le ZIP pour le décompresser.
2. Ouvrez le **Terminal**. Tapez `bash ` (avec un espace après), faites glisser `install_autoindex_plugin.sh` dans la fenêtre, puis appuyez sur Entrée. Les deux fichiers sont copiés dans `~/.GSASII/exports/`.

**Windows**
1. Téléchargez le ZIP comme ci-dessus, puis clic droit dessus → **Extraire tout**.
2. Double-cliquez sur **`install_windows.bat`**. Les deux fichiers sont copiés dans `C:\Users\<votre nom>\.GSASII\exports\`.

### Étape 2 : installer Conograph (conseillé)

1. Téléchargez Conograph depuis la [page de téléchargement de Z-Code](https://z-code-software.com/downloads/) et décompressez-le dans le dossier de votre choix. On peut aussi le compiler à partir de son code source.
2. Notez le chemin complet de la version **en ligne de commande** de Conograph. Sous Windows, l'installeur de Z-Code fournit plusieurs `.exe` : utilisez **`Conograph CUI.exe`** (par défaut `C:\Program Files (x86)\Conograph\Conograph CUI.exe`, trouvé automatiquement), et non `Conograph GUI.exe` ni `PeakSearch CUI.exe`.
3. L'extension le trouve toute seule s'il est dans le chemin du système ou si la variable d'environnement `CONOGRAPH` pointe vers lui. Sinon, indiquez son chemin la première fois dans le champ **Exécutable Conograph** de la fenêtre. Il est mémorisé ensuite.

### Étape 3 : redémarrer GSAS-II

## Utilisation

1. Ouvrez un projet GSAS-II avec un diagramme de poudre, et remplissez la *Peak List* (idéalement avec des pics ajustés).
2. Ouvrez **Export → Entire project as → (user) Auto-index (find unit cell)**.
3. Dans la fenêtre, choisissez l'histogramme et les méthodes. Si besoin, ajustez les réglages :
   - nombre maximal de pics utilisés ;
   - décalage de zéro ;
   - tolérance sur 2θ ;
   - durée maximale de chaque recherche ;
   - recherche Conograph rapide ou poussée.
4. Validez. Les méthodes tournent en parallèle ; vous pouvez annuler à tout moment.
5. Lisez le compte rendu, cochez les mailles à garder et validez. Elles apparaissent dans la **Unit Cells List** de l'histogramme.
6. Enregistrez le projet.

**En ligne de commande** (avec le Python de GSAS-II). Les mailles sont alors écrites dans une **copie** du projet :

```
python g2_autoindex.py projet.gpx                        # demande les méthodes à lancer
python g2_autoindex.py projet.gpx --methods c,d,t        # Conograph, dichotomie, essai-erreur
python g2_autoindex.py projet.gpx --methods c --maxpeaks 25
python g2_autoindex.py --help                            # toutes les options
```

## Limites

- Les moteurs dichotomie et essai-erreur ne traitent pas le triclinique. Utilisez Conograph ou GSAS-II dans ce cas.
- L'indexation dépend de la qualité des positions de pics. Des raies d'impureté, ou trop peu de raies, peuvent donner des solutions fausses ou aucune solution. Vérifiez toujours la maille proposée, par exemple par un affinement de Le Bail.
- Développé et utilisé sur macOS. Windows devrait fonctionner, mais n'a pas encore été testé. Vos retours sont les bienvenus.

## Crédits et références

- **GSAS-II :** B. H. Toby, R. B. Von Dreele, *J. Appl. Cryst.* **2013**, 46, 544–549. La méthode native utilise le module d'indexation de GSAS-II.
- **Conograph** (indexation de poudre) : « New functions and graphical user interface attached to powder indexing software », *J. Appl. Cryst.* **2017**, 50, 651–659.
- **DICVOL** (méthode par dichotomie) : D. Louër, M. Louër, *J. Appl. Cryst.* **1972**, 5, 271–275 ; A. Boultif, D. Louër, *J. Appl. Cryst.* **2004**, 37, 724–731.
- **TREOR** (méthode par essai-erreur) : P.-E. Werner, L. Eriksson, M. Westdahl, *J. Appl. Cryst.* **1985**, 18, 367–370.
- **Figures de mérite :** M20 (de Wolff, *J. Appl. Cryst.* **1968**, 1, 108–113) et F20 (Smith & Snyder, *J. Appl. Cryst.* **1979**, 12, 60–65).
- Écrit par Yoann Rousselin avec l'aide de Claude (Anthropic) comme assistant de programmation.

## Licence

MIT. Voir [LICENSE](LICENSE). Conograph et GSAS-II sont distribués sous leurs propres licences.

## Auteur

**Yoann Rousselin**, ingénieur de recherche CNRS (cristallographie), ICMUB, Dijon.
Questions ou signalement d'un problème : ouvrez une *Issue* sur cette page.
