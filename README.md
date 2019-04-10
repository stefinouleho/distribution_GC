# DISTRIBUTION DES TAILLES DES GRAPHES DE CYCLES & GRAPHES MOLECULAIRES DANS LA BASE DE DONNÉE chebi 
Similarité sur le graphe de cycles des molécules

- La base de données de molecules utilisée est CHEBI, il faut télécharger le fichier CHebi_Lite.pdf(https://www.dropbox.com/s/7k0ef9rmxvao70q/ChEBI_lite.sdf?dl=0) et le mettre dans le dossier.


Il faut lancer: make run  qui va generer :

+ le fichier /stats/distribution_cycles.data qui contient la distribution des tailles de cycles dans les GC.
+ le fichier /stats/distribution_atomes.data qui contient la distribution des tailles d'atomes dans les graphes molécuaires
+ le fichier /stats/moyenne_cycles.data qui contient la moyenne  des tailles de cycles dans les GC.
+ le fichier /stats/moyenne_atomes.data qui contient la moyenne  des tailles d'atomes dans les graphes molécuaires.


Ensuite excecuter le fichier gplt : gnuplot distribution.gplt pour affichier l'histogramme de distribution des tailles de cycles dans les GC.


Version à jour!
