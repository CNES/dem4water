# dem4water

## Approche

1. Caractérisation de la zone

  + Hypotèse : localisation de l'ouvrage
  + Cercles concentriques pour identification [amont/aval] / ouvrage (= transverse à la direction d'écoulement)
  + [amont/aval] + carte d'occurence d'eau -> amont / aval -> localisation pied du barage -> Z_0

2. S(Z_i)

  + Génération des courbes de niveau intersectées par la transverse à la direction d'écoulement (~ ouvrage)
  + 2ème méthode _Analyse locale du MNT au sein du masque polygone fourni_ sur la base du masque d'occurrence (from SurfWater)


3. S(Z)

  + Suppressions des outliers
  + Fit du modèle: détermination de alpha, beta et Z_0

## Prise en main

Le point d'entré privilégié est le fichier d'exemple [camp.sh](helper/camp.sh) qui met en oeuvre la chaîne de traitement de bout en bout et propose un méchanisme de traitement en lot simple. Il est bien évidemment à adapter, notament la section de définition des données d'entré.
Le script [camp.sh](helper/camp.sh) de lancement d'une campagne se charge de modifier l'environement pour permettre l'accès aux dépendances et exécuter la chaîne dans les conditions optimales. Par défault, le script [camp.sh](helper/camp.sh) lance les applications en mode débogage qui affiche des informations supplémentaires sur les conditions d'opérations de la chaîne lors de sont exécution.
Chaque application contient la définition de ces paramètres (accessible par --help).
