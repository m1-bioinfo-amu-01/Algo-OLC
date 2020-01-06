# Bienvenue dans le readme utilisateur 

## objectif 

Cet algorithme d'assemblage de type OLC permet de faire du "Gap Filling".
Avec en donnés un jeu de reads de type Illumina et deux -mers (start et stop) , le but est de
reconstruire une séquence, construite à partir des reads chevauchants, débutant par le kmer start
 et terminant par le kmer stop .
 
 ![Image description](link-to-image)

## getting started 

verifiez la presence des fichier :
- OLC.py
- un fichier .fasta /.fa contenant les read a assembler 
- un fichier .fasta /.fa contenant les kmers start puis stop dans cet ordre

## lancement 
veuillez executer l'ensemble du code present dans OLC.py 
--> TODO fonction / option et detail 
'' ## extention exacte : 
 veuillez executer :  nom_fonction(arguments)
 
## extention avec mismatch : 
 veuillez executer :  nom_fonction(arguments)

## resultat 
apres execution le programme affiche un message pour confirmer la bonne execution / echec  et en parallèle 
est cree un nouveau fichier fasta contenant la séquence "solution"
