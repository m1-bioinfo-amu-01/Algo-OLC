# Algo-OLC
## Bienvenue dans le read-me developpeur 

L'objectif de cet algorithme d'assemblage OLC est de faire du " Gap Filling " : reconstruire toutes les séquences à partir des fichiers d'entrées. Pour ce faire, l'agorithme doi trouver, dans les reads données en entrée,un kmer start à partir duquel il essaye d'étendre une séquence nucléotidique. Il concatène donc les reads chevauchants de façon exacte et ce, jusqu'à ce qu'il retourne un kmer Stop. 

Ce programme nécessite 2 fichiers fasta en entrée : 1 contenant les reads type Illumina que l'on souhaite étendre et un fichier contenant les codons Start et Stop.  

### Structure principale 

- import des modules nécessaires 

- fonctions:

  #### parserMultiFASTA : 

  Cette fonction permet d'ouvrir le fichier FASTA et de récupérer un dictionnaire nommé "tab_id_seq_pos" de la forme suivante 
   { identifiant des read : [sequence du read : position kmer start s'il y a ] }
   
   ex :
   
   |read3 | [ ATTGG... ] |  --> ici read sans kmer start dedans
   
   

   |read4 | [ CGTCA... , 90]| --> ici read avec un codon start commençant au 90 ème nucléotide
   

   #### parser_start_stop

   lis le fichier et stock la séquence du kmer start dans la variable start et la séquence de stop dans la variable stop. 


#### try_stop(file path):
not used yet

#### extend(read,seed,matrice,kmer,stop,result):

  Fonction récursive qui va permettre d'étendre les séquences et retourner une matrice resultat contenant : 
  - le path actuel qui est en train de s'étendre
  - une liste all_path qui va stocké l'ensemble des chemins trouvés pour les différents start ce qui va permettre de vider à chaque fois la liste "path"
  - une chaîne de caractère "seq" qui va contenir la séquence étendu du path actuel
  - une liste all_seq stockant tutes les séquencences des paths trouvés
  
  Structure : 
  
  boucle for qui parcours l'ensemble des positions de départ possible pour kmer_test. 
   
   verifie si les kmer test n'est pas un kmer stop : 
      
      si oui : 
      
      	return resulat 
      
      sinon : 
      
      	regarde si le test dans le tableau seed(liste des 1 kmer):
        
        	test que ce n'est pas déjà lui même :
          
         		 test si le reste de la partie chevauchante et aussi égale: 
            
            			si oui on enlève cette référence de seed pour éviter une boucle infinie 
            
            			test le read permet bien d'etendre la séquence: 
              
             				ajout de la partie non chevauchante dans la séquence etendue 
              
              				recursion 



TODO fonction pour le remplissage de tab_premiers_kmers ( pour l'instant boucle qui parcours la matrice et répertorie les n premier nt + id read associe dans seed & ajout de l'élément position dans matrice). Seed permet de récupérer tous les reads qui partagent le premier kmer. 
        
-  variables : 

tab_id_seq_pos : tableau contenant {identifiant des read : [sequence du read : position kmer start s'il y a ] }

start : séquence kmer start

pos_start = position du kmer start

stop : sequence kmer stop

pos_stop = position du kmer stop

taille_kmer : longueur du kmer à tester 

tab_premiers_kmers : tableau contenant tous les read qui partagent le meme premier kmer

result : tableau contenant le chemin (liste des read utilisés pour étendre le read de départ) et la séquence concatenée obtenue après extention

kmer_test = premier kmer du read dans tab_id_seq_pos, c'est le kmer à tester

id_read_for_extend = read utilisé pour l'extension







 
