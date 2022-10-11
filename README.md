+ je recupere le nom de l'échantillon dans l'en-tête du BAM.
+ j'extrais les reads avec une MAPQ>=30 et qui entourent l'intervalle `chr7:27_135_310-27_135_339` (pas de chevauchement partiel).
+ Pour chaque read
    + je récupère la sequence d'ADN du read
    + je récupère la CIGAR string du read qui décrit les matchs/dels/ins
    + je traverse la cigar strin et quand je suis uniquement dans l'interval chr7:27_135_310-27_135_339, je reconstitue la sequence d'ADN de cet intervalle
    + l'ADN est reverse-complenté et traduit en peptide
    + je scanne le peptide pour le transformer en 'symbole'
+ chaque sequence d'ADN différente est comptée
+ le séquences presentes <=2 fois sont ignorées
+ les deux sequences restantes les plus fréquentes sont affichées