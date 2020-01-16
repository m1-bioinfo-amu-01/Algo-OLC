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

###- fonctions:

#### parserMultiFASTA : 

  Cette fonction permet d'ouvrir le fichier FASTA et de récupérer un dictionnaire nommé "tab_id_seq_pos" de la forme suivante 
   { identifiant des read : [sequence du read : position kmer start s'il y a ] }
   
   ex :
   
   |read3 | [ ATTGG... ] |  --> ici read sans kmer start dedans
   
   

   |read4 | [ CGTCA... , 90]| --> ici read avec un codon start commençant au 90 ème nucléotide


#### parser_start_stop

   lis le fichier et stock la séquence du kmer start dans la variable start et la séquence de stop dans la variable stop. 
   
#### rev_comp(sequence):
renvoie le reverse complement de la séquence nucléotidique passée en argument.


#### try_stop_start(matrice,id_read,kmer_to_find):

cette fonction permet de chercher un kmer donné que ce soit dans un read forward ou dans son reverse complement. Retourne la position du kmer ou -1 si ce dernier est absent dans un read forward. 
Retourne la position negative du kmer (-pos) ou None si ce dernier est absent dans un read reverse complement.



#### extend(read,seed,matrice,kmer,stop,result):

  Fonction récursive qui va permettre d'étendre les séquences et retourner une matrice resultat contenant : 
  - le path actuel qui est en train de s'étendre
  - une liste all_path qui va stocké l'ensemble des chemins trouvés pour les différents start ce qui va permettre de vider à chaque fois la liste "path"
  - une chaîne de caractère "seq" qui va contenir la séquence étendu du path actuel
  - une liste all_seq stockant toutes les séquencences des paths trouvés
  
  Structure : 
  
```

  une copie profonde du dictionnaire est utilisée afin de ne pas rendre indisponible certaine donées lors du parcours de l'ensemble des chemins
  
 si le read passé en argument contient un codons stop la recursion s'arrete les donées de path et de sequence contenue dans la copie sont alors "sauvegardées" dans le dictionnaire initial sous les clefs: all_path et all_seq

sinon :
  boucle for qui parcours l'ensemble des positions de départ possible pour notre kmer a tester :
    verifie si l'id du read est en forward ou en reverse :
      si en reverse il faut rev-comp(seq du read)
      
      	regarde si le test dans le tableau seed (liste des 1 kmer):
        
        pour chaque element test que ce n'est pas déjà lui même :
          
         		 test si le reste de la partie chevauchante et aussi égale: 
            
            			si oui on enlève cette référence de seed pour éviter une boucle infinie 
            
            			test le read permet bien d'etendre la séquence: 
              
             				ajout de la partie non chevauchante dans la séquence etendue 
              
              				recursion 
```

#### fill_tab_premier_kmers(matrice_all_reads,taille_kmer,start):
fonction qui permet de creer un dictionnaire contenant toutes les séquences les premier kmer associé au id read ou -id read si c'est le debut du reverse complement 
de la forme : {sequenceA: id_read1, -id_read5, sequenceB: id_read2 ... }

#### all_extentions (path_for_all_read, path_file_start_stop,length_kmer):
fonction qui permet de realiser l'extension parfaite et de creer un fichier texte contenant les resultats obtenus
```
chargement des données des fichier dans les variables

boucle qui parcours tous les reads: 
  si il contiennent un codon start:
  initialisation de la recursion 
  resulat de la recursion ajouté a une liste all_result

creation et ecriture du fichier de sortie 
  chemin 
  sequence etendue
fermeture du fichier

```
        
### -  variables : 

tab_id_seq_pos : dictionnaire contenant {identifiant des read : [sequence du read : position kmer start s'il y a ] }

start : séquence du kmer start

pos_start = position du kmer start

stop : sequence kmer stop

pos_stop = position du kmer stop

taille_kmer : longueur du kmer à tester 

tab_premiers_kmers : dictionnaire contenant tous les read qui partagent le meme premier kmer {sequence : id des read}

result : tableau contenant le chemin (liste des read utilisés pour étendre le read de départ) et la séquence concatenée obtenue après extention

kmer_test = premier kmer du read dans tab_id_seq_pos, c'est le kmer à tester

id_read_for_extend = read utilisé pour l'extension

path_file_read : deuxieme argument passé en ligne de commande (premier apres le nom du fichier) lors de l'execution c'est le chemin complet sans espaces du fichier contenant l'ensemble des reads

path_file_start_stop: troisieme argument passé en ligne de commande  lors de l'execution c'est le chemin complet sans espaces du fichier contenant les read start et stop

kmer_len : quatrieme argument passé en ligne de commande  lors de l'execution c'est la taille des kmer a utiliser pour l'extension

option : cinquieme argument passé en ligne de commande  lors de l'execution c'est l'option choisie pour l'extention a realiser peut etre -p ou --perfect (TODO -m ou --mismatch)







 
