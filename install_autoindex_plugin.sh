#!/bin/bash
# Installe l'indexation automatique (recherche de maille) dans GSAS-II (menu Export) -- macOS/Linux
# Installs the auto-indexing (unit-cell search) plugin into GSAS-II (Export menu) -- macOS/Linux
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
DEST="$HOME/.GSASII/exports"
mkdir -p "$DEST"
cp "$HERE/G2export_AutoIndex.py" "$HERE/g2_autoindex.py" "$DEST/"
echo "Fichiers copies dans / Files copied to: $DEST"
echo
if command -v Conograph >/dev/null 2>&1 || [ -n "$CONOGRAPH" ] || [ -x "$HOME/ProjectPowderConograph/Conograph/src/Conograph" ]; then
  echo "Conograph trouve / found."
else
  echo "Conograph non trouve automatiquement : indiquez son chemin dans la fenetre du plugin."
  echo "Conograph not found automatically: give its path in the plugin window."
  echo "(Les methodes dichotomie, essai-erreur et GSAS-II fonctionnent sans lui / other methods work without it.)"
fi
echo
echo "Redemarrez GSAS-II / Restart GSAS-II: Export -> Entire project as -> (user) Auto-index (find unit cell)"
