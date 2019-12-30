# Algo-OLC
## Bienvenue dans le read-me developpeur 

## Fichier OLC.py 
Ce programme necessite 2 fichier fasta en entrée 1 contenant les read que l'on souhaite etendre et un fichier contenant les condon START et STOP. 

Il a pour objectif de trouver dans les read un codon start a partir du quel il essaye d'entendre une sequence nucleotidique. Pour cela il concatene les read chevauchants de facon exacte et ce tant qu'il ne rencontre pas de codon STOP. 

### Structure principale 
- import des modules necessaires 

- fonctions:
  #### parserMultiFASTA : 
  cette fonction permet d'ouvrir le fichier fasta et de recuperer un dictionnaire nommé "matrice" de la forme suivante 
   IMAGE 

   #### parser_start_stop
   lis le fichier et stock la séquence du codon start dans la variable start et la séquence de stop dans la variable stop. 


TODO remplissage de seed dans une fonction ( parcours de la matrice et repertorie les n premier nt + id read associe)

#### try_stop(file path):
not used yet

#### extend(read,seed,matrice,kmer,stop,result):
  boucle for qui parcours chaque nucléotide:(pour decaler de 1 le kmer a tester)
    verifie si les kmer test n'est pas un codon stop : 
      si oui : 
      return resulat 
      sinon : 
      regarde si le test dans le tableau seed(liste des 1 kmer):
        test que ce n'est pas deja lui meme :
          test si le reste de la partie chevauchante et aussi egale: 
            si oui on enleve cette reference de seed pour eviter une boucle infinie 
            test le read permet bien d'etendre la séquence: 
              ajout de la partie non chevauchante dans la séquence etendue 
              recursion 
        
-  variables : 

matrice : tableau contenant  identifiant des read | séquence du read 
eventuellement une troisieme "colone" avec si il y a un codon start sa position de depart 

start : séquence codon start

stop : sequence codon stop

kmer :Longueur du kmer seed

seed : tableau contenant tout les read qui partagent le meme premier kmer

result : tableau contenant le chemin (liste des read utilises pour etendre le read de depart )et la séquence concatenee obtenue apres extention


 
