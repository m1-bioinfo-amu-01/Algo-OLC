# Algo-OLC
## Bienvenue dans le read-me developpeur 

## Fichier OLC.py 
Ce programme necessite 2 fichier fasta en entrée 1 contenant les read que l'on souhaite etendre et un fichier contenant les condon START et STOP. 
Il a pour objectif de trouver dans les read un codon start a partir du quel il essaye d'entendre une sequence nucleotidique. Pour cela il concatene les read chevauchants de facon exacte et ce tant qu'il ne rencontre pas de codon STOP. 

### Structure principale 
- import des modules necessaires 
- ouvertures des fichiers fasta 
