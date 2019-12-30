#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 30 14:58:00 2019

@author: anissachibani
"""

import sys

sys.setrecursionlimit(50000)


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
def extend(read,seed,matrice,kmer,stop,result):
	result['path'].append(read) # ajout du read dans le chemin
	for NT in range(0,len(matrice[read][0])-kmer+1):
		test= matrice[read][0][NT:NT+kmer] # le premier kmer du read
		try:
			pos = matrice[read][0].find(stop) # cherche un codon stop dans le read d'interet si il y en a l'indice de sa position sera contenu dans pos
		except IndexError :
			pos = -1  # pas de codon stop trouvé
		if pos != -1:    # si il y a un codon stop
			print("stop")
			print("result",result)
			return result  # devrais terminer le programme
		if test in seed:  # compare avec le dico des premier kmers des reads
		for essai in seed[test]: # on test tout les read qui ont le meme premier kmer
			if essai != read: # si on ne retombe pas sur lui meme
				print("read ", read, "essai ", essai)
				if matrice[read][0][NT:] == matrice[essai][0][0:len(matrice[read][0][NT:])] : # si tout le reste de la seq est egalement communue
					indextest=seed[test].index(essai) #on recupere l'index de l'element utiliser pour etendre
					seed[test].pop(indextest) # on l'enleve de la liste des read partagant le premier kmer pour ne pas boucler infiniment
					if len((matrice[essai][0][len(matrice[read][0][NT:]):])) !=0: #verifie que l'extention apportée par le read apporte bien des nouveau nucléotides ( extention suppérieure a 0)
						result['seq']+=(matrice[essai][0][len(matrice[read][0][NT:]):]) # extention dxe la séquence
						print("seq :",len(result['seq']) ) # print la longueur de la séquence pour pouvoir suivre son extension
						extend(essai,seed,matrice,kmer,stop,result) # recursion


'''chargement des données '''
matrice = parserMultiFASTA('/Users/anissachibani/Desktop/M2/ALG/Projet/data/ecoli_2kb_perfect_forward_reads.fasta')

start, stop = parser_start_stop('/Users/anissachibani/Desktop/M2/ALG/Projet/data/start_stop_2kb.fa')

'''def des variables'''
kmer=len(start) # taille de notre kmer

seed={} # matrice contant le 1kmer de tout les read de la forme {seq: id reads qui commencent par ca}

for read in matrice: # remplissage de seed
	first_kmer= matrice[read][0][0:kmer]
	if first_kmer not in seed and first_kmer!= '':
		seed[first_kmer]=[read]
	else:
		if first_kmer != '':
			seed[first_kmer].append(read)
	try:
		pos = matrice[read][0].find(start) # en meme temps test dans chaque read si il y a un start et si ou l'indice de sa position est contenu dans pos
	except IndexError :
		pos=-1
	if pos != -1:
		matrice[read].append(pos) #si il y en a un matrice prend la forme {id read : [sequence , position start]}


result = {'path': [], 'seq': ''} # creation du dictionnaire qui va nous permettre se stocker le "chemin" d'assemblage, et la séquence étendue



'''code principal : appel des fonction '''
for i in matrice:
	if len(matrice[i])>1 :
		pos=matrice[i][1]
		result['seq']=matrice[i][0][pos:] # initialisation de la séquence avec le premier read de la matrice contenant un start
		extend(i,seed,matrice,kmer,stop,result)
