# Bienvenue dans le readme utilisateur 

## Objectif 

Cet algorithme d'assemblage de type OLC permet de faire du "Gap Filling".
Avec en données un jeu de reads de type Illumina et deux k-mers (start et stop), le but est de
reconstruire une séquence, construite à partir des reads chevauchants, débutant par le k-mer start
 et terminant par le k-mer stop .
 
 ![Image description](https://github.com/m1-bioinfo-amu-01/Algo-OLC/blob/master/OLC%2COverlap-Layout-consensus.png)

## Avant de demarrer

Verifiez la presence des fichier :
- OLC.py
- un fichier .fasta /.fa contenant les reads a assembler 
- un fichier .fasta /.fa contenant les k-mers start puis stop dans cet ordre

Assurez vous de bien avoir les packages time et resource

## Lancement 
Dans le terminal veuillez executer la commande comme suit :
```
python3 olc.py '/Users/ninamenet/PycharmProjects/olc/ecoli_2kb_perfect_forward_reads.fasta' '/Users/ninamenet/PycharmProjects/olc/start_stop_2kb.fa' 31 -p
```
Avec en agrgument 1 : le chemin complet du fichier contenant les reads

En argument 2 : le chemin complet du fichier contenant les reads start et stop

En argument 3 : la taille du k-mer que vous souhaitez utiliser pour faire l'extension 

En argument 4 : l'option qui peut etre -p ou --perfect pour faire une extension sans mismatch
 
## Sortie
Apres exécution le programme affiche dans le terminal un message pour confirmer la bonne exécution "code had run without encountering error" , le temps d'exécution et la memoire utilisée. 
Par ailleurs créé localement un fichier texte "resultat.txt" contenant les résultats sous la forme :

path :['read0,read1']
extended sequence : ['ATCTGAATAACATCCT']

path :['read1', 'read3']
extended sequence : ['ATCCT', 'ATCCTACCGTC']

## Warning
"WARNING unknown option" 'WARNING kmer length is not a valid number' FileNotFoundError peuvent apparaitre notamment si les chemins contiennent des espaces et sont donnés sans être entre '' ou si les arguments passés sont invalides.

Attention le fichier fasta contenant le start et stop doit les présenter de facon à ce que le premier kmer soit toujours le start 

La taille de kmer doit etre inferieure à la taille des read utilisé

# resultats

## temps de calculs et memoire
nous avons tester les programme sur 2 jeu de données disponible dans le dossier data_alg 
un jeu de donnée de 2kb et un de 100kb dont voici un graph presentant l'evolution des temps utilisé et de la memoire en fonction de la taille de kmer utilisé pour realiser l'assemblage 
![Image description](link-to-image)


![Image description](link-to-image)
Nous n'avons pas pu réaliser exactement les mêmes tests que pour le fichier précédent pour des raisons de temps, de mémoire et de performance de nos ordinateurs. Ce qui soulèvent une des limites les plus importantes de notre programme

Nous avons été jusque 150 dans nos test si l'execution est plus rapide cela s'explique car il s'agit de la taille des reads donc aucune extention ne peut etre faite.

Comme l'on peut le voir ici que le compromis ideal entre temps et memoire est au alentour de 20-30 nt 
Sachant que nos read on une longeur de 150 nt on estime qu'un bon point de depart semble etre au alentours de 20% de la taille de kmer cependant il est important de prendre en compte la complexitées des données entrée. Il est recommandé sur des données biologique d'utiliser des nombres impair afin d'eviter les palindromes qui peuvent etre frequents. Ici nous ne l'avons pas fait car ce sont des données artificielles et que la différence n'était pas significative.

## qualité
pour le fichier 100 kb nous trouvons un alignement de taille 10013 et 10 chemin possibles ce qui est relativement proche de ce qui est attendu (10001) les 12 nucléotides de différences s'explique par le fait que la séquence ne s'arrete pas pile apres le codon stop 

pour le fichier 2kb nous trouvons une séquence etendue de 1065 nucléotides ce qui est tres proche de ce qui est attendu (1065).

on peut conclure que notre outil realise un alignement de qualité

cependant les différent chemin produisent la meme sequence se qui nous a fait douter de la veracité de la séquence cepandant sur un toy exemple nous avons bien observé le comportement attendu. Ci-dessous :

path :['read0,read1']
extended sequence : ['ATCTGAATAACATCCT']

path :['read1', 'read3']
extended sequence : ['ATCCT', 'ATCCTACCGTC']

nous n'avons pas eu le temps suffisant pour verififer manuellement tout les assemblages mais vu que cela fonctionne sur les toy example nous estimons que cela devrait en théorie etre correct.

## disscusion 

Nos temps de calculs et la memoire utilisée restent relativement elevé ce qui rend sont application a de grand jeu de données non envisable.
Cela est du entre autre de notre utilisation de multiple deep-copy lors de notre recursion et de l'etat encore un peu "brut" de notre code qui pourrait encore etre bien plus optimisé.

De plus nous avons decider de conserver et d'afficher les chemin utilisé pour l'assembage meme si cela n'est pas forcement indispensable il nous semblais interessant la mesure ou nous avons realiser l'extention a partir de tout les starts possible pour pouvoir les diferenncier mais egalement car cela nous a beaucoup aider a voir nos erreurs lorsque nous codions. 

notre code ne realise pour l'instant malheureusement seulement l'extension pour des match parfaits nous avons manqué de temps pour realiser une option qui permettrait de faire l'extention avec des mismatch. Cependant nous y avions reflechis et nous avions prevu de faire cela avec une score seuil d'erreur en dessous du quel la comparaison de la sequence chevauchante serait considérée comme valable et continuerait la recursion.

