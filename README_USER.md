# Bienvenue dans le readme utilisateur 

## Objectif :

Cet algorithme d'assemblage de type OLC permet de faire du "Gap Filling", avec en données, un jeu de reads de type Illumina et deux k-mers (start et stop). Le but est de reconstruire une séquence, construite à partir de reads chevauchants, débutant par un k-mer start et terminant par un k-mer stop.
 
 Voici ci dessous un graphe qui represente l'assemblage de type olc et que nous avons essayé de reprosuire dans cet outil
 ![Image description](https://github.com/m1-bioinfo-amu-01/Algo-OLC/blob/master/OLC%2COverlap-Layout-consensus.png)

## Avant de démarrer :

Vérifiez la présence des fichiers :
- OLC.py
- un fichier .fasta /.fa contenant les reads à assembler 
- un fichier .fasta /.fa contenant les k-mers start puis stop dans cet ordre

Assurez-vous de bien avoir les packages "time" et "ressource".

## Lancement :
Dans le terminal veuillez exécuter la commande comme suit :
```
python3 olc.py '/Users/ninamenet/PycharmProjects/olc/ecoli_2kb_perfect_forward_reads.fasta' '/Users/ninamenet/PycharmProjects/olc/start_stop_2kb.fa' 31 -p
```
Avec en argument 1 : le chemin complet du fichier contenant les reads

En argument 2 : le chemin complet du fichier contenant les k-mers start et stop

En argument 3 : la taille du k-mer que vous souhaitez utiliser pour faire l'extension 

En argument 4 : l'option qui peut être soit -p ou --perfect pour faire une extension sans mismatch
 
## Sortie :

Apres l'exécution, le programme affiche dans le terminal un message pour confirmer la bonne exécution : "code had run without encountering error", avec le temps d'exécution et la mémoire utilisée. 
Par ailleurs, cela a également créér localement un fichier texte "resultat.txt" contenant les résultats sous la forme :

path :['read0,read1']
extended sequence : ['ATCTGAATAACATCCT']

path :['read1', 'read3']
extended sequence : ['ATCCT', 'ATCCTACCGTC']

## Warning :

Les messages "WARNING unknown option" "WARNING kmer length is not a valid number" "FileNotFoundError" peuvent apparaître notamment si les chemins contiennent des espaces et sont donnés sans être entre guillemets, ou si les arguments passés sont invalides.

Attention, pour le fichier fasta contenant le start et stop, le premier k-mer doit être le start, et le second sera le stop et non l'inverse. 

La taille du k-mer doit être inférieure à la taille des reads utilisée.

# Résultats :

## Temps de calculs et mémoire :

Nous avons testé le programme sur 2 jeux de données disponibles dans le dossier data_alg.

Un jeu de donnée de 2kb et un de 100kb dont voici les graphiques représentant l'évolution des temps utilisés et de la mémoire en fonction de la taille de k-mer utilisée pour réaliser l'assemblage.

Gaphique pour le jeu de donnée de 2Kb:
![Image description](https://github.com/m1-bioinfo-amu-01/Algo-OLC/blob/master/2kb.png)

Gaphique pour le jeu de donnée de 100Kb:
![Image description](https://github.com/m1-bioinfo-amu-01/Algo-OLC/blob/master/100kb.png)


Pour le jeu de données de 100kb, nous n'avons pas pu effectué les tests de mémoire et de temps avec une taille de k-mer supérieure à 70 nucléotides, pour des raisons de temps, de mémoire et de performance de nos ordinateurs. Ce qui est une des limites les plus importantes de notre programme.Nous avons egalement fait tourner avec le jeu de donné de 5000kb nous avons obtenu un resultat pour un kmer de taille 31 :1524.0472841262817 secondes et 809066496 kilobytes ce qui nous a disuader de tenter ce graphique avec ce jeu de donné.

Nous remarquons dans le jeu de données de 2kb, avec une taille de k-mer de 150 nucléotides, que le temps d'exécution est beaucoup plus rapide et que cela prend moins d'espace. Ceci peut s'expliquer du fait de la taille importante du k-mer, qui ne permet que très peu voire aucune extension, donc la fin du programme est vite atteinte. 

Comme l'on peut le voir ici, le compromis idéal entre le temps et la mémoire est aux alentours de 20 nucléotides. 
Sachant que nos reads ont une longueur de 150 nucléotides, on estime qu'une bonne taille de k-mer correspond à 15% de la longueur des reads utilisés.

Cependant, il est important de prendre en compte la complexité des données d'entrées. En effet, il est recommandé, sur des données biologiques, d'utiliser des nombres impairs afin d'éviter les palindromes qui peuvent être fréquents. Ici nous ne l'avons pas fait car ce sont des données artificielles et que la différence n'était pas significative.

## Qualité :

Pour le fichier de 100 kb, nous trouvons un alignement de taille 10013 et 10 chemins possibles, ce qui est relativement proche de ce qui est attendu (10001). Les 12 nucléotides de différences s'expliquent par le fait que l'alignement ne s'arrête directement au dernier nucléotide du kmer stop, mais continue jusqu'à la fin du read contenant ce même k-mer.

Pour le fichier 2kb, nous trouvons une séquence étendue de 1065 nucléotides, ce qui correspond au résultat attendu (1065).

On peut conclure que notre outil réalise un alignement de qualité.

Cependant, les différents chemins produisent la même séquence, ce qui nous a fait douter de la véracité de la séquence. Nous avons alors testé notre programme sur un toy exemple afin de comprendre les raisons de ce résultat, et nous observons bien le comportement attendu. Les résultats obtenus sont présentés ci-dessous :

path :['read0,read1']
extended sequence : ['ATCTGAATAACATCCT']

path :['read1', 'read3']
extended sequence : ['ATCCT', 'ATCCTACCGTC']

Nous n'avons pas eu assez de temps pour vérifier manuellement tous les assemblages, mais étant donné que cela fonctionne sur notre toy example, nous estimons que cela devrait en théorie être correct.

## Discussion :

Nos temps de calculs ainsi que la mémoire utilisée restent relativement élevée, ce qui rend l'application de notre programme pour de grands jeux de données non envisageable.

Cela est dû entre autre à notre utilisation de multiple "deep-copy" lors de notre récursion et de l'état encore un peu "brut" de notre code, qui pourrait encore être bien plus optimisé.

De plus, nous avons décidé de conserver et d'afficher les chemins utilisés pour l'assembage, même si cela n'est pas indispensable. Il nous semblait intéressant, dans la mesure où nous avons realisé l'extension à partir de tous les starts possibles, de les garder afin de pouvoir les différencier, mais cela nous a également servi pour voir nos erreurs au moment de la programmation.

Notre code ne réalise pour l'instant que l'extension pour des matchs parfaits. Nous avons manqué de temps pour réaliser une option qui permettrait de faire l'extension avec des mismatch. Cependant, nous y avions réfléchis et nous avions prévus de faire cela avec un score de seuil d'erreur, en dessous duquel la comparaison de la séquence chevauchante serait considérée comme valable et continuerait la récursion.

