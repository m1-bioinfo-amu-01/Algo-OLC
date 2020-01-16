# Bienvenue dans le readme utilisateur 

## objectif 

Cet algorithme d'assemblage de type OLC permet de faire du "Gap Filling".
Avec en donnés un jeu de reads de type Illumina et deux -mers (start et stop) , le but est de
reconstruire une séquence, construite à partir des reads chevauchants, débutant par le kmer start
 et terminant par le kmer stop .
 
 ![Image description](link-to-image)

## avant de demarrer

verifiez la presence des fichier :
- OLC.py
- un fichier .fasta /.fa contenant les read a assembler 
- un fichier .fasta /.fa contenant les kmers start puis stop dans cet ordre

assurez vous de bien avoir les packages time et resource

## lancement 
dans le terminal veuillez executer la commande comme suit :
```
python3 olc.py '/Users/ninamenet/PycharmProjects/olc/ecoli_2kb_perfect_forward_reads.fasta' '/Users/ninamenet/PycharmProjects/olc/start_stop_2kb.fa' 31 -p
```
avec en agrgument 1 : le chemin complet du fichier contenant les reads 
en argument 2 : le chemin complet du fichier contenant les read start et stop
en argument 3 : la taille du kmer que vous souhaitez utiliser pour faire l'extention 
en argument 4 : l'option qui peut etre -p ou --perfect pour faire une extention sans mismatch

## Sortie
apres execution le programme affiche dans le terminal un message pour confirmer la bonne execution "code had run without encountering error" , le temps d'execution et la memoire utilisée 
par ailleur cree localement un fichier texte "resultat.txt" contenant les resultats sous la forme :

path :['read0,read1']
extended sequence : ['ATCTGAATAACATCCT']

path :['read1', 'read3']
extended sequence : ['ATCCT', 'ATCCTACCGTC']

## Warning
"WARNING unknown option" 'WARNING kmer length is not a valid number' FileNotFoundError peuvent apparaitre notamment si les chemin contiennent des espaces et sont donnés sans etre entre '' ou si les arguments passé sont invalides.

Attention le fichier fasta contenant le start et stop doit les presenter de facon à ce que le premier kmer soit toujours le start 

La taille de kmer doit etre inferieure à la taille des read utilisé

# resultats

## temps de calculs et memoire
nous avons tester les programme sur 2 jeu de données disponible dans le dossier data_alg 
un jeu de donnée de 2kb et un de 100kb dont voici un graph presentant l'evolution des temps utilisé et de la memoire en fonction de la taille de kmer utilisé pour realiser l'assemblage 
![Image description](link-to-image)

comme l'on peut le voir ici que le compromis ideal entre temps et memoire est au alentour de 20-30 nt 
--> psq c'est de la taille du start/stop ?? 
nous avons été jusque 150 dans nos test si c'est plus rapide cela s'explique car il s'agit de la taille des reads donc aucune extention ne peut etre faite 
![Image description](link-to-image)

## qualité

## disscusion 
nos temps de calculs et la memoire utilisée restent relativement elevé ce qui rend sont application a de grand jeu de données inenvisable.
Cela est du entre autre de notre utilisation de multiple deep-copy lors de notre recursion et de l'etat encore un peu "brut" de notre code qui pourrait encore etre bien plus optimisé. deplus nous avons decider de conserver et d'afficher les chemin utilisé pour l'assembage meme si cela n'est pas forcement indispensable il nous semblais interessant la mesure ou nous avons realiser l'extention a partir de tout les starts possible pour pouvoir les diferenncier mais egalement car cela nous a beaucoup aider a voir nos erreurs lorsque nous codions. 



