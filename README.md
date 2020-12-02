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
