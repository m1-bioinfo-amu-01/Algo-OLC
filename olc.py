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

def try_stop(matrice,id_read,stop):
	try:
		pos = matrice[id_read][0].find(stop)
	except IndexError :
		pos = -1 
	return pos

''' fonction recursive qui permet d'etendre la sequence du read initial
    - si le read contient un codon stop return le dictionnaire resultat
    - sinon compare le read avec seed pour trouver d'autres read qui partagent le meme kmer
    - si on essaye pas de comparer le read avec lui meme ob regarde si bien pareil sur tout la longeur
    - si oui on relance la fonction avec le nouveau read et on etend la sequence'''

def extend(id_read, tab_premiers_kmers, tab_id_seq_pos, taille_kmer, stop, result):
	result['path'].append(id_read) # ajout du read dans le chemin
	
	for pos_read_tab in range(0,len(tab_id_seq_pos[id_read][0])-taille_kmer+1): # pos_read_matrice parcours chaque nucléotide des reads dans matrice (cf parserMultiFASTA(nom_fichier))
		
		kmer_test= tab_id_seq_pos[id_read][0][pos_read_tab:pos_read_tab+taille_kmer] # test récupère le premier kmer du read dans matrice = kmer à tester
		
		try:
			pos_stop = tab_id_seq_pos[id_read][0].find(stop) # find(stop) cherche s'il y a un stop dans les read, et si oui, pos prend l'indice du stop 
			
		except IndexError :
			pos_stop = -1  # pas de kmer stop trouvé
			
	# on commence tout d'abord par vérifier s'il y a un stop pour pouvoir terminer le programme 
	
		if pos_stop != -1:    # si il y a un kmer stop
			
			print("stop")
			path_actuel = ','.join(result['path']
			result['path']=[]
			result['all_path'].append(path_actuel)
			seq_actuelle = result['seq']
			result['seq'] = []result['all_seq'].append(seq_actuelle)
			print("result", result)
			return result  # fin programme + retourne dictionnaire resultat
		
	# si le kmer stop est différent du read à tester (test), alors on regarde si test est dans tab_premiers_kmers (=tableau contenant tous les reads partageant le meme premier kmer)  
	
		if kmer_test in tab_premiers_kmers:  # compare avec le dicco des premiers kmers des reads
			
		for id_read_ in tab_premiers_kmers[kmer_test]: # id_read_for_extend parcours tab_premiers_kmers[test] pour tester tous les reads qui ont le même premier kmer

			if id_read_for_extend != id_read: # on vérifie que id_read_for_extend ne retombe pas sur lui-même
				print("read ", id_read, "id_read_for_extend ", id_read_for_extend)
				

				# on vérifie si le reste de la partie chevauchante est la même --> on cherche dans matrice de pos_read_matrice jusqu'à la fin pour read VS  ?????
				# /!\ MIEUX EXPLIQUER LE IF CI-DESSOUS ! 
				if tab_id_seq_pos[id_read][0][pos_read_tab:] == tab_id_seq_pos[id_read_for_extend][0][0:len(tab_id_seq_pos[id_read][0][pos_read_tab:])] :
					# /!\ MIEUX EXPLIQUER LE indice_test 					
					indice_test=tab_premiers_kmers[kmer_test].index(id_read_for_extend) #indice_test récupère l'indice de l'élément utilisé pour étendre 
					
					tab_premiers_kmers[kmer_test].pop(indice_test) # on enlève indice_test de tab_premiers_kmer[test] pour ne pas boucler infiniment

					if len((tab_id_seq_pos[id_read_for_extend][0][len(tab_id_seq_pos[id_read][0][pos_read_tab:]):])) !=0: #verifie que l'extention apportée par le read apporte bien de nouveaux nucléotides ( extention supérieure à 0)
						
						result['seq']+=(tab_id_seq_pos[id_read_for_extend][0][len(tab_id_seq_pos[id_read][0][pos_read_tab:]):]) # ajout de la partie non chevauchante dans result['seq']
						print("seq :",len(result['seq']) ) # affiche la longueur de la séquence pour pouvoir suivre son extension
						
						extend(id_read_for_extend,tab_premiers_kmers,tab_id_seq_pos,taille_kmer,stop,result) # recursion


'''chargement des données '''

tab_id_seq_pos = parserMultiFASTA('/Users/anissachibani/Desktop/M2/ALG/Projet/data/ecoli_2kb_perfect_forward_reads.fasta')
start, stop = parser_start_stop('/Users/anissachibani/Desktop/M2/ALG/Projet/data/start_stop_2kb.fa')

'''def des variables'''

taille_kmer=len(start) # taille de notre kmer
tab_premiers_kmers={} # matrice contant le 1kmer de tout les read de la forme {seq: id reads qui commencent par ça}  /!\ CA ?? 


# /!\ MIEUX EXPLIQUER LA PARTIE EN DESSOUS 
for id_read in tab_id_seq_pos: # remplissage de tab_premiers_kmers
	premier_kmer= tab_id_seq_pos[id_read][0][0:taille_kmer] # premier_kmer prend le premier kmer dans matrice
	if premier_kmer not in tab_premiers_kmers and premier_kmer!= '': # vérifie que premier_kmer n'est pas dans tab_premiers_kmers et qu'il n'est pas vide 
		tab_premiers_kmers[premier_kmer]=[id_read] #  
	else:
		if premier_kmer != '':
			tab_premiers_kmers[premier_kmer].append(id_read)
	try:
		pos_start = tab_id_seq_pos[id_read][0].find(start) # find(start) cherche s'il y a un start dans les read, et si oui, pos_start prend l'indice du start 
	except IndexError :
		pos_start=-1
	if pos_start != -1:
		tab_id_seq_pos[id_read].append(pos_start) #si il y a un start, la matrice prend la forme {id read : [sequence , position start]}


result ={'path': [], 'all_path': [],
          'seq': '', 'all_seq':[]} # creation du dictionnaire qui va nous permettre se stocker le "chemin" d'assemblage, et la séquence étendue

'''code principal : appel des fonction '''
for i in tab_id_seq_pos:
	if len(tab_id_seq_pos[i])>1 :
		pos=tab_id_seq_pos[i][1]
		result['seq']=tab_id_seq_pos[i][0][pos:] # initialisation de la séquence avec le premier read de la matrice contenant un start
		extend(i,tab_premiers_kmers,tab_id_seq_pos,taille_kmer,stop,result)
