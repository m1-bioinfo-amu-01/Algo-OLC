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

# resultats

## temps de calculs

## utilisation memoire 

## qualité

## conclusion

