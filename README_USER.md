# Bienvenue dans le readme utilisateur 

## Objectif :

Cet algorithme d'assemblage de type OLC permet de faire du "Gap Filling".
Avec en données un jeu de reads de type Illumina et deux k-mers (start et stop), le but est de
reconstruire une séquence, construite à partir des reads chevauchants, débutant par le k-mer start
 et terminant par le k-mer stop .
 
 ![Image description](https://github.com/m1-bioinfo-amu-01/Algo-OLC/blob/master/OLC%2COverlap-Layout-consensus.png)

## Avant de démarrer :

Verifiez la presence des fichier :
- OLC.py
- un fichier .fasta /.fa contenant les reads a assembler 
- un fichier .fasta /.fa contenant les k-mers start puis stop dans cet ordre

Assurez vous de bien avoir les packages time et resource

## Lancement :
Dans le terminal veuillez executer la commande comme suit :
```
python3 olc.py '/Users/ninamenet/PycharmProjects/olc/ecoli_2kb_perfect_forward_reads.fasta' '/Users/ninamenet/PycharmProjects/olc/start_stop_2kb.fa' 31 -p
```
Avec en agrgument 1 : le chemin complet du fichier contenant les reads

En argument 2 : le chemin complet du fichier contenant les reads start et stop

En argument 3 : la taille du k-mer que vous souhaitez utiliser pour faire l'extension 

En argument 4 : l'option qui peut etre -p ou --perfect pour faire une extension sans mismatch
 
## Sortie :
Apres exécution le programme affiche dans le terminal un message pour confirmer la bonne exécution "code had run without encountering error" , le temps d'exécution et la memoire utilisée. 
Par ailleurs créé localement un fichier texte "resultat.txt" contenant les résultats sous la forme :

path :['read0,read1']
extended sequence : ['ATCTGAATAACATCCT']

path :['read1', 'read3']
extended sequence : ['ATCCT', 'ATCCTACCGTC']

## Warning :
"WARNING unknown option" 'WARNING kmer length is not a valid number' FileNotFoundError peuvent apparaitre notamment si les chemins contiennent des espaces et sont donnés sans être entre '' ou si les arguments passés sont invalides.

Attention le fichier fasta contenant le start et stop doit les présenter de facon à ce que le premier k-mer soit toujours le start 

La taille de k-mer doit être inferieure à la taille des reads utilisés.

# Résultats :

## Temps de calculs et mémoire :
Nous avons testé les programmes sur 2 jeux de données disponible dans le dossier data_alg.
Un jeu de donnée de 2kb et un de 100kb dont voici un graphique présentant l'évolution des temps utilisés et de la mémoire en fonction de la taille de k-mer utilisé pour réaliser l'assemblage.

Gaphique pour le jeu de donnée de 2Kb:
![Image description](https://github.com/m1-bioinfo-amu-01/Algo-OLC/blob/master/2kb.png)

Gaphique pour le jeu de donnée de 100Kb:
![Image description](https://github.com/m1-bioinfo-amu-01/Algo-OLC/blob/master/100kb.png)
Nous n'avons pas pu réaliser exactement les mêmes tests que pour le fichier précédent pour des raisons de temps, de mémoire et de performance de nos ordinateurs. Ce qui soulève une des limites les plus importantes de notre programme.

Nous avons été jusque 150 dans nos test si l'exécution est plus rapide cela s'explique parce qu'il s'agit de la taille des reads donc aucune extension ne peut être faite.

Comme l'on peut le voir ici que le compromis idéal entre temps et mémoire est aux alentours de 20-30 nt 
Sachant que nos reads ont une longeur de 150 nt on estime qu'un bon point de départ semble être aux alentours de 20% de la taille de k-mer cependant il est important de prendre en compte la complexitées des données entrées. Il est recommandé sur des données biologique d'utiliser des nombres impairs afin d'éviter les palindromes qui peuvent être fréquents. Ici nous ne l'avons pas fait car ce sont des données artificielles et que la différence n'était pas significative.

## Qualité :
Pour le fichier 100 kb nous trouvons un alignement de taille 10013 et 10 chemins possibles ce qui est relativement proche de ce qui est attendu (10001) les 12 nucléotides de différences s'expliquent par le fait que la séquence ne s'arrête pas pile après le codon stop.

Pour le fichier 2kb nous trouvons une séquence étendue de 1065 nucléotides ce qui est très proche de ce qui est attendu (1065).

On peut conclure que notre outil réalise un alignement de qualité.

Cependant les différents chemins produisent la même séquence ce qui nous a fait douter de la véracité de la séquence cepandant sur un toy exemple nous avons bien observé le comportement attendu. Ci-dessous :

path :['read0,read1']
extended sequence : ['ATCTGAATAACATCCT']

path :['read1', 'read3']
extended sequence : ['ATCCT', 'ATCCTACCGTC']

Nous n'avons pas eu assez de temps pour vérififer manuellement tout les assemblages mais vu que cela fonctionne sur les toy example nous estimons que cela devrait en théorie être correct.

## Discussion :

Nos temps de calculs et la mémoire utilisée restent relativement elevée ce qui rend sont application a de grands jeux de données non envisageable.
Cela est dû entre autre à notre utilisation de multiple deep-copy lors de notre récursion et de l'état encore un peu "brut" de notre code qui pourrait encore être bien plus optimisé.

De plus nous avons décider de conserver et d'afficher les chemins utilisés pour l'assembage même si cela n'est pas forcement indispensable il nous semblait intéressant dans la mesure où nous avons realisé l'extention à partir de tout les starts possible pour pouvoir les diferenncier mais également car cela nous a beaucoup aidé à voir nos erreurs lorsque nous codions. 

Notre code ne réalise pour l'instant que l'extension pour des matchs parfaits nous avons manqués de temps pour realiser une option qui permettrait de faire l'extension avec des mismatch. Cependant nous y avions reflechis et nous avions prévu de faire cela avec un score seuil d'erreur en dessous duquel la comparaison de la séquence chevauchante serait considérée comme valable et continuerait la récursion.

