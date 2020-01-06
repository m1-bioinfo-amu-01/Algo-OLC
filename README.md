# Algo-OLC
## Bienvenue dans le read-me developpeur 


Ce programme nécessite 2 fichiers fasta en entrée : 1 contenant les reads que l'on souhaite étendre et un fichier contenant les codons START et STOP.  


Il a pour objectif de trouver dans les reads un codon start à partir duquel il essaye d'étendre une séquence nucléotidique. Pour cela il concatène les reads chevauchants de façon exacte et ce, jusqu'à ce qu'il retourne un codon STOP. 

### Structure principale 

- import des modules nécessaires 

- fonctions:

  #### parserMultiFASTA : 

  cette fonction permet d'ouvrir le fichier FASTA et de récupérer un dictionnaire nommé "matrice" de la forme suivante 
   | ID read (read+n°)| liste de 1 ou 2 éléments : séquence et éventuellement position du start |
   
   ex :
   
   |read3 | [ ATTGG... ] |  --> ici read sans kmer start dedans
   
   

   |read4 | [ CGTCA... , 90]| --> ici read avec un codon start commençant au 90 ème nucléotide
   

   #### parser_start_stop

   lis le fichier et stock la séquence du codon start dans la variable start et la séquence de stop dans la variable stop. 


TODO fonction pour le remplissage de seed ( pour l'instant boucle qui parcours la matrice et répertorie les n premier nt + id read associe dans seed & ajout de l'élément position dans matrice). Seed permet de récupérer tous les reads qui partagent le premier kmer. 

#### try_stop(file path):
not used yet

#### extend(read,seed,matrice,kmer,stop,result):

  boucle for qui parcours chaque nucléotide:(pour décaler de 1 le kmer a tester):
   
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


        
-  variables : 

matrice : tableau contenant  identifiant des read | séquence du read 
eventuellement une troisieme "colone" avec si il y a un kmer start sa position de depart 

start : séquence kmer start

stop : sequence kmer stop

kmer :Longueur du kmer seed

seed : tableau contenant tout les read qui partagent le meme premier kmer

result : tableau contenant le chemin (liste des read utilises pour etendre le read de depart )et la séquence concatenee obtenue apres extention

test = read dans matrice que l'on veut tester 

essai = read dans seed[test] quand le read test est présent dans seed

 
