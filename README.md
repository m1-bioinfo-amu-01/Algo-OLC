# Algo-OLC

## Bienvenue dans le read-me developpeur 

L'objectif de cet algorithme d'assemblage OLC est de faire du " Gap Filling " : reconstruire des séquences à partir des fichiers d'entrées. Pour ce faire, l'agorithme doit trouver, dans les reads données en entrée, un kmer start à partir duquel il essaye d'étendre une séquence nucléotidique. Il concatène donc les reads chevauchants de façon exacte et ce, jusqu'à ce qu'il trouve un kmer Stop. 

Ce programme nécessite 2 fichiers fasta en entrée : 1 contenant les reads type Illumina que l'on souhaite étendre et un fichier contenant les codons Start et Stop.  

## Structure principale 

### - import des modules nécessaires: 
  - sys 
  - copy
  - time 
  - resource

### - fonctions:

#### parserMultiFASTA : 

  Cette fonction permet d'ouvrir le fichier FASTA et de récupérer un dictionnaire nommé "tab_id_seq_pos" de la forme suivante 
   { identifiant des read : [sequence du read : position kmer start s'il y a ] }
   
   ex :
   
   |read3 | [ ATTGG... ] |  --> ici read sans kmer start dedans
   
   

   |read4 | [ CGTCA... , 90]| --> ici read avec un codon start commençant au 90 ème nucléotide


#### parser_start_stop

   lis le fichier et stocke la séquence du kmer start dans la variable start et la séquence de stop dans la variable stop. 
   
#### rev_comp(sequence):

   renvoie le reverse complément de la séquence nucléotidique passée en argument.


#### try_stop_start(matrice,id_read,kmer_to_find):

   Cette fonction permet de chercher un kmer donné, que ce soit dans un read forward ou dans son reverse complement.
   Il retourne la position du kmer ou -1 si ce dernier est absent dans un read forward. 
   Il retourne la position negative du kmer (-pos) ou None si ce dernier est absent dans un read reverse complement.

#### extend(id_read, tab_premiers_kmers, tab_id_seq_pos, taille_kmer, stop, result):

  Fonction récursive qui va permettre d'étendre les séquences et retourner une matrice resultat contenant : 
  - le path actuel qui est en train de s'étendre
  - une liste all_path qui va stocker l'ensemble des chemins trouvés pour les différents start ce qui va permettre de vider à chaque fois la liste "path"
  - une chaîne de caractère "seq" qui va contenir la séquence étendue du path actuel
  - une liste all_seq stockant toutes les séquencences des paths trouvés
  
  Structure : 
  
```

Une copie profonde du dictionnaire est utilisée afin de ne pas rendre indisponible certaines données lors du parcours de l'ensemble des chemins.
  
Si le read passé en argument contient un kmer stop la récursion s'arrête. Les données de path et de sequence contenues dans la copie sont alors "sauvegardées" dans le dictionnaire initial sous les clefs : all_path et all_seq

Sinon :
  Boucle for qui parcours l'ensemble des positions de départ possible pour notre kmer à tester :
  
     Vérifie si l'id du read est en reverse :
     
          Si c'est le cas, le kmer_test sera en reverse, en utilisant cette commande : rev-comp(seq du read)
          Sinon, le kmer_test utilisé sera bien un forward
      
     On regarde si le kmer test est dans le tableau tab_premiers_kmers:
        
          Si oui, on vérifie que le read que l'on souhaite étendre n'est pas égal à lui-même:
          
         	    S'il est différent de lui-même, on définit les séquences pour les chevauchements et extensions en fonction de   s'ils sont des reverse ou pas.  
              
              On test si le reste de la partie chevauchante est aussi égale: 
            
            	      Si oui, on vérifie que la séquence s'est bien étendue en vérifiant sa longueur :
            
            			        Si oui, on ajoute la partie non chevauchante dans la séquence etendue 
                          Appel de la fonction extend grâce à la récursion   

```

#### fill_tab_premier_kmers(matrice_all_reads,taille_kmer,start):

Fonction qui permet de créer un dictionnaire tab_premiers_kmer contenant toutes les séquences des premiers kmer associés au "id read",si la séquence est en forward, ou "-id read" si c'est le début d'un reverse complément.
tab_premiers_kmer est de la forme : {sequenceA: id_read1, -id_read5, sequenceB: id_read2 ... }

#### all_extentions (path_for_all_read, path_file_start_stop,length_kmer):

Fonction qui permet de realiser l'extension parfaite et de créer un fichier texte contenant les resultats obtenus.

```
Chargement des données des fichier dans les variables

Boucle qui parcours tous les reads: 
      Si le read contient un kmer start:
            Initialisation de la récursion pour recommencer à partir de ce start
            Résulat de la recursion est ajouté à une liste all_result

Création et écriture du fichier de sortie contenant : 
      - chemin 
      - séquence étendue
Fermeture du fichier

```
        
### -  variables : 

tab_id_seq_pos : dictionnaire contenant {identifiant des read : [sequence du read : position kmer start s'il y a ] }

start : séquence du kmer start

pos_start = position du kmer start

stop : séquence du kmer stop

pos_stop = position du kmer stop

taille_kmer : longueur du kmer à tester 

tab_premiers_kmers : dictionnaire contenant tous les read qui partagent le meme premier kmer {sequence : id des read}

result : tableau contenant le chemin (liste des read utilisés pour étendre le read de départ) et la séquence concatenée obtenue après extention

kmer_test = premier kmer du read dans tab_id_seq_pos, c'est le kmer à tester

id_read_for_extend = id du read utilisé pour l'extension

path_file_read : deuxième argument passé en ligne de commande (premier après le nom du fichier) lors de l'execution. Il correspond au chemin complet, sans espaces, du fichier contenant l'ensemble des reads

path_file_start_stop: troisième argument passé en ligne de commande  lors de l'execution. Il correspond au chemin complet,sans espaces, du fichier contenant les reads start et stop

kmer_len : quatrième argument passé en ligne de commande lors de l'execution. Il correspond à la taille des kmer à utiliser pour l'extension

option : cinquième argument passé en ligne de commande lors de l'execution. Il correspond à l'option choisie pour l'extention à realiser: peut être -p ou --perfect (TODO -m ou --mismatch)







 
