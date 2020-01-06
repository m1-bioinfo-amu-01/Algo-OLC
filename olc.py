#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 30 14:58:00 2019

@author: anissachibani
"""

import sys

sys.setrecursionlimit(50000) # augmente le maximum de récursion


''' cette fonction permet d'ouvrir le fichier fasta et de recuperer un dictionnaire de la forme :
{ ID(read+numero) : ["sequence"]}'''

def parserMultiFASTA(nom_fichier):
	seq = []
	resultat = {}
	ID = ''
	f = open(nom_fichier, "r")
	while True:
		uneLigne = f.readline()
		if len(uneLigne) == 0:
			resultat[ID] = [seq_str.join(seq)]
			break
		if uneLigne.startswith(">"):
			seq_str = ''
			resultat[ID] = [seq_str.join(seq)]
			seq = []
			ID = uneLigne[1:-1]
			ID = ID.split('_')[0]
		else:
			seq.append(uneLigne.strip(' \n'))
	# on est dans la sequence ou sur une ligne vide
	f.close()

	return resultat

''' cette fonction permet d'ouvrir le fichier contenant la séquence stop et la séquence start :
    la premiere séquence est sauvegardée comme etant le stop et la deuxieme le start.
    retourne un tuple de 2 variables
'''
def parser_start_stop(nom_fichier):
	seq = []
	start = ''
	stop = ''
	f = open(nom_fichier, "r")
	while True:
		uneLigne = f.readline()
		if len(uneLigne) == 0:
			stop = stop.join(seq)
			break
		if uneLigne.startswith(">"):
			start = start.join(seq)
			seq = []
		else:
			seq.append(uneLigne.strip(' \n'))
	# on est dans la sequence ou sur une ligne vide
	f.close()

	return start, stop

''' fonction  qui permet de determinier si il y a un codon stop dans le read d'interet'''
# actuellement on ne l'utilise pas encore

def try_stop(matrice,read,stop):
	try:
		pos = matrice[read][0].find(stop)
	except IndexError :
		pos = -1 
	return pos

''' fonction recursive qui permet d'etendre la sequence du read initial
    - si le read contient un codon stop return le dictionnaire resultat
    - sinon compare le read avec seed pour trouver d'autres read qui partagent le meme kmer
    - si on essaye pas de comparer le read avec lui meme ob regarde si bien pareil sur tout la longeur
    - si oui on relance la fonction avec le nouveau read et on etend la sequence'''

def extend(read,tab_premiers_kmers,matrice,kmer,stop,result):
	result['path'].append(read) # ajout du read dans le chemin
	
	for pos_read_matrice in range(0,len(matrice[read][0])-kmer+1): # pos_read_matrice parcours chaque nucléotide des reads dans matrice (cf parserMultiFASTA(nom_fichier))
		
		test= matrice[read][0][pos_read_matrice:pos_read_matrice+kmer] # test récupère le premier kmer du read dans matrice = kmer à tester
		
		try:
			pos_stop = matrice[read][0].find(stop) # find(stop) cherche s'il y a un stop dans les read, et si oui, pos prend l'indice du stop 
			
		except IndexError :
			pos_stop = -1  # pas de kmer stop trouvé
			
	# on commence tout d'abord par vérifier s'il y a un stop pour pouvoir terminer le programme 
	
		if pos_stop != -1:    # si il y a un kmer stop
			
			print("stop")
			print("result",result)
			return result  # fin programme + retourne dictionnaire resultat
		
	# si le kmer stop est différent du read à tester (test), alors on regarde si test est dans tab_premiers_kmers (=tableau contenant tous les reads partageant le meme premier kmer)  
	
		if test in tab_premiers_kmers:  # compare avec le dicco des premiers kmers des reads
			
		for essai in tab_premiers_kmers[test]: # essai parcours tab_premiers_kmers[test] pour tester tous les reads qui ont le même premier kmer

			if essai != read: # on vérifie que essai ne retombe pas sur lui-même
				print("read ", read, "essai ", essai)
				
				# on vérifie si le reste de la partie chevauchante est la même --> 
				# /!\ MIEUX EXPLIQUER LE IF CI-DESSOUS ! 
				seq_read_to_extend = matrice[read][0][pos_read_matrice:]
				seq_new_read_to_test =matrice[essai][0][0:len(matrice[read][0][pos_read_matrice:])]
				if seq_read_to_extend  == seq_new_read_to_test: # on cherche dans matrice de pos_read_matrice(debut de la seq a tester) jusqu'à la fin pour read a étendre VS le nouveau read du debut a la fin de la partie chevauchante
					# /!\ MIEUX EXPLIQUER LE indice_test 					
					indice_test=tab_premiers_kmers[test].index(essai) #indice_test récupère l'indice de l'élément utilisé pour étendre dans la liste des 1ers kmer commun
					
					tab_premiers_kmers[test].pop(indice_test) # on enlève indice_test de tab_premiers_kmer[test] pour ne pas boucler infiniment
					
					if len((matrice[essai][0][len(matrice[read][0][pos_read_matrice:]):])) !=0: #verifie que l'extention apportée par le read apporte bien de nouveaux nucléotides ( extention supérieure à 0)
						
						result['seq']+=(matrice[essai][0][len(matrice[read][0][pos_read_matrice:]):]) # ajout de la partie non chevauchante dans result['seq']
						print("seq :",len(result['seq']) ) # affiche la longueur de la séquence pour pouvoir suivre son extension
						
						extend(essai,tab_premiers_kmers,matrice,kmer,stop,result) # recursion


'''chargement des données '''

matrice = parserMultiFASTA('/Users/anissachibani/Desktop/M2/ALG/Projet/data/ecoli_2kb_perfect_forward_reads.fasta')
start, stop = parser_start_stop('/Users/anissachibani/Desktop/M2/ALG/Projet/data/start_stop_2kb.fa')

'''def des variables'''

kmer=len(start) # taille de notre kmer
tab_premiers_kmers={} # matrice contant le 1kmer de tout les read de la forme {seq: id reads qui commencent par ça}  /!\ CA ?? 

 
for read in matrice: # remplissage de tab_premiers_kmers boucle for qui parcours l'ensemble des reads
	premier_kmer= matrice[read][0][0:kmer] # premier_kmer prend le premier kmer dans matrice
	if premier_kmer not in tab_premiers_kmers and premier_kmer!= '': # vérifie que premier_kmer n'est pas deja dans tab_premiers_kmers et qu'il n'est pas vide 
		tab_premiers_kmers[premier_kmer]=[read] # on ajoute l'id de ce read a une liste contenant l'enemble des reads partageant ce premier kmer
	else: # si il est deja dedans 
		if premier_kmer != '': # le kmer n'est pas vide 
			tab_premiers_kmers[premier_kmer].append(read) # on ajoute le read  a la liste  deja existante pour ce kmer
	try:
		pos_start = matrice[read][0].find(start) # find(start) cherche s'il y a un start dans les read, et si oui, pos_start prend l'indice du start 
	except IndexError : # quand il ne trouve pas de start find renvoi un IndexError
		pos_start=-1 # on defini alors la valeur de "absence de kmer start dans ce read" comme etant -1
	if pos_start != -1: 
		matrice[read].append(pos_start) #si il y a un start, la matrice prend la forme {id read : [sequence , position start]}


result = {'path': [], 'seq': ''} # creation du dictionnaire qui va nous permettre se stocker le "chemin" d'assemblage, et la séquence étendue

'''code principal : appel des fonction '''
for i in matrice:
	if len(matrice[i])>1 :
		pos=matrice[i][1]
		result['seq']=matrice[i][0][pos:] # initialisation de la séquence avec le premier read de la matrice contenant un start
		extend(i,tab_premiers_kmers,matrice,kmer,stop,result)
