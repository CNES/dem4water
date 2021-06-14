# dem4water

## Description

The aim of dem4water is to estimate water surface elevation (Z), surface area (S) and volume (V) relationships of reservoirs. It uses as input a digital elevation model (DEM) acquired after the reservoir construction and a water occurrence map (watermap). 

## Execution

The main script used to run the chain is [camp.sh](helper/camp.sh), gathering all steps together to process a list of dams. Using the script requires adapting the input definition section (first part of camp.sh file) as well as the dam list.
The script [camp.sh](helper/camp.sh) run a whole campaign and is in charge of setting the right environment and dependencies.
Each application embbed a documentation that can be accessed using --help.

## Examples

The following examples are directly extracted from [camp.sh](helper/camp.sh)

### Step 1 - area_mapping

This application extracts area specific data from the input watermap and DEM.   Then, it uses the water map and concentric cercles on the DEM to identify upstream/downstream areas (= orthogonal to flow direction) and the dam bottom location. 

``` sh
python3 area_mapping.py --debug \
  --id       "${DAMID}" \
  --infile   "${DB_PATH}" \
  --watermap "${WMAP_PATH}" \
  --dem      "${DEM_PATH}" \
  --radius   "$RADIUS" \
  --out      "$EXTR_DIR/${DAM}_${RADIUS}"
```

### Step 2 - szi_from_contourline

This application uses the dam bottom location to derive Z0 as well as the dam cut line (imaginary line cutting the valley on the dam position). It also extracts the contour lines from the DEM.

```sh
python3 szi_from_contourline.py --debug \
  --id           "${DAMID}" \
  --infile       "${DB_PATH}" \
  --watermap     "$EXTR_DIR/${DAM}_${RADIUS}/wmap_extract-$DAM.tif" \
  --dem          "$EXTR_DIR/${DAM}_${RADIUS}/dem_extract-$DAM.tif"  \
  --radius       "$RADIUS" \
  --pdbstep      5 \
  --pdbradius    500 \
  --elevsampling 1 \
  --elevoffset   60 \
  --tmp          "$ROOT_DIR/${DAM}_${RADIUS}/tmp" \
  --out          "$ROOT_DIR/${DAM}_${RADIUS}"
```

### Step 3 - cut_contourlines

Using the dam cutting line as well as the contour lines, this application estimates the virtual surfaces and the associated S(Zi).

```sh
python3 cut_contourlines.py --debug \
  --dem      "$EXTR_DIR/${DAM}_${RADIUS}/dem_extract-$DAM.tif"  \
  --info     "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_daminfo.json" \
  --cut      "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_cutline.json" \
  --level    "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_contourlines@1m.json" \
  --out      "$ROOT_DIR/${DAM}_${RADIUS}"
```

### Step 4 - szi_to_model

This application computes model parameters of the bathymetric relationship based on the best S(Zi) subset.

```sh
python3 szi_to_model.py --debug \
  --infile     "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_SZi.dat" \
  --daminfo    "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_daminfo.json" \
  --zminoffset 10 \
  --zmaxoffset 30 \
  --maemode    "first" \
  --outfile    "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_model.png"
```

### Step 5 - val_report

This application compares the estimated model to a ground truth model to evaluate its quality.

```sh
python3 val_report.py --debug \
  -i "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_model.json" \
  -r "${GT_PATH}" \
  -o "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_report.json"
```

# dem4water [FR]

## Approche

1. Caractérisation de la zone (step 1)

  + Hypotèse : localisation de l'ouvrage
  + Cercles concentriques pour identification [amont/aval] / ouvrage (= transverse à la direction d'écoulement)
  + [amont/aval] + carte d'occurence d'eau -> amont / aval -> localisation pied du barage -> Z_0

2. S(Z_i)  (step 2 et 3)

  + Génération des courbes de niveau intersectées par la transverse à la direction d'écoulement (~ ouvrage)
  + 2ème méthode _Analyse locale du MNT au sein du masque polygone fourni_ sur la base du masque d'occurrence (from SurfWater)


3. S(Z)    (step 4)

  + Suppressions des outliers
  + Fit du modèle: détermination de alpha, beta et Z_0

## Prise en main

Le point d'entré privilégié est le fichier d'exemple [camp.sh](helper/camp.sh) qui met en oeuvre la chaîne de traitement de bout en bout et propose un méchanisme de traitement en lot simple. Il est bien évidemment à adapter, notament la section de définition des données d'entré.
Le script [camp.sh](helper/camp.sh) de lancement d'une campagne se charge de modifier l'environement pour permettre l'accès aux dépendances et exécuter la chaîne dans les conditions optimales. Par défault, le script [camp.sh](helper/camp.sh) lance les applications en mode débogage qui affiche des informations supplémentaires sur les conditions d'opérations de la chaîne lors de sont exécution.
Chaque application contient la définition de ces paramètres (accessible par --help).

## Exemples d'utilisation

Les exemples sont donnés à titre d'illustration et sont directement extraitrs de [camp.sh](helper/camp.sh)

### area_mapping

Cette application extrait la zone d'intérêt autour de la localisaiton du barrage.

``` sh
python3 area_mapping.py --debug \
  --id       "${DAMID}" \
  --infile   "${DB_PATH}" \
  --watermap "${WMAP_PATH}" \
  --dem      "${DEM_PATH}" \
  --radius   "$RADIUS" \
  --out      "$EXTR_DIR/${DAM}_${RADIUS}"
```

### szi_from_contourline

Cette application détecte le pied de barrage pour estimer Z0 ainsi que la ligne de coupure. Elle prépare également les lignes de niveau de la vallée.

```sh
python3 szi_from_contourline.py --debug \
  --id           "${DAMID}" \
  --infile       "${DB_PATH}" \
  --watermap     "$EXTR_DIR/${DAM}_${RADIUS}/wmap_extract-$DAM.tif" \
  --dem          "$EXTR_DIR/${DAM}_${RADIUS}/dem_extract-$DAM.tif"  \
  --radius       "$RADIUS" \
  --pdbstep      5 \
  --pdbradius    500 \
  --elevsampling 1 \
  --elevoffset   60 \
  --tmp          "$ROOT_DIR/${DAM}_${RADIUS}/tmp" \
  --out          "$ROOT_DIR/${DAM}_${RADIUS}"
```

### cut_contourlines

Cette application exploite la ligne de coupure et les lignes de niveau afin d'estimer les surfaces virtuelles S(Zi).

```sh
python3 cut_contourlines.py --debug \
  --dem      "$EXTR_DIR/${DAM}_${RADIUS}/dem_extract-$DAM.tif"  \
  --info     "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_daminfo.json" \
  --cut      "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_cutline.json" \
  --level    "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_contourlines@1m.json" \
  --out      "$ROOT_DIR/${DAM}_${RADIUS}"
```

### szi_to_model

Cette application estime les paramètres du modèles en optimisant le choix des S(Zi) utilisés pour cette estimation.

```sh
python3 szi_to_model.py --debug \
  --infile     "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_SZi.dat" \
  --daminfo    "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_daminfo.json" \
  --zminoffset 10 \
  --zmaxoffset 30 \
  --maemode    "first" \
  --outfile    "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_model.png"
```

### val_report

Cette application permet d'évaluer la qualité du modèle estimé en le comparant à un modèle de référence.

```sh
python3 val_report.py --debug \
  -i "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_model.json" \
  -r "${GT_PATH}" \
  -o "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_report.json"
```
