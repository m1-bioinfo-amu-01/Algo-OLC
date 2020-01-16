import sys
from copy import deepcopy
import time
import resource


start_time = time.time() #initialisation du calcul du temps ecoulé

sys.setrecursionlimit(50000)  # augmente le maximum de récursion

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


''' cette fonction permet d'ouvrir le fichier contenant kmer stop et le kmer start :
    la premiere séquence est sauvegardée comme étant le stop et la deuxieme le start.
    Cela retourne un tuple de 2 variables. 
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

''' définition de la fonction qui permet de retourner le reverse complément de la sequence passée en argument '''

def rev_comp(sequence):
    reverse=''
    dico = {'A':'T','C':'G','G':'C','T':'A'}
    for nt in sequence:
        if nt in dico :
            reverse += dico[nt]
    return reverse

''' fonction qui permet de déterminer si le codon d'intérêt est présent dans le read dont l'identifiant est passé en argument.
Si cet identifiant commence par "-", cela signifie qu'il faut regarder dans le reverse complément.
Pour l'identifiant forward, s'il est présent dans le read que l'on test, pos = nombre positif, et s'il est absent, pos = -1.
Pour les read reverse, pos = nombre négatif si présent et None si absent '''

def try_stop_start(matrice, id_read, kmer_to_find):
    if id_read.startswith('-'):
        real_id_read = id_read[1:]
        rev= rev_comp(matrice[real_id_read][0])
        pos = rev.find(kmer_to_find)
        if pos != -1:
            return  pos*(-1)
        else :
            return None
    else:
        try:
            pos = matrice[id_read][0].find(kmer_to_find)
        except IndexError:
            pos = -1
        return pos


''' fonction récursive qui permet d'étendre la séquence du read initial
    - si le read contient un kmer stop, le programme retourne le dictionnaire resultat
    - sinon, il compare le read avec tab_premiers_kmers pour trouver d'autres reads qui partagent le même kmer
    - on vérifie ensuite que l'on est pas en train de comparer deux même kmers entre eux
    - on regarde ensuite si le reste de la partie chevauchante est la même
    - si oui, on relance, on rajoute le chemin et la séquence pour étendre la séquence
    - et enfin on relance la fonction avec le nouveau read grace à la récursion '''


def extend(id_read, tab_premiers_kmers, tab_id_seq_pos, taille_kmer, stop, result):
    res_copy = deepcopy(result)
    res_copy['path'].append(id_read)  # ajout du read dans le chemin de la copie profonde
    pos_stop = try_stop_start(tab_id_seq_pos,id_read,stop)  # on commence tout d'abord par vérifier s'il y a un stop pour pouvoir terminer le programme
    if pos_stop != -1:  # s'il y a un kmer stop, le path et la séquence contenus dans la copie est alors basculée dans le vrai dictionnaire resultat
        path_actuel = ','.join(res_copy['path']) 
        result['all_path'].append(path_actuel)
        seq_actuelle = ''.join(res_copy['seq'])
        result['all_seq'].append(seq_actuelle)
        return result  # fin programme + retourne dictionnaire resultat

    for pos_read_tab in range(0, len(tab_id_seq_pos[id_read][0]) - taille_kmer + 1):  # pos_read_tab parcours chaque
        # nucléotide des reads dans tab_id_seq_pos (cf parserMultiFASTA(nom_fichier))
        if id_read.startswith('-'): # si le read est dans le sens reverse
            id_read_2 = id_read[1:] # les reverse ne sont pas stockés directement, il faut récupérer l'id forward pour avoir la séquence forward
            kmer_test = rev_comp(tab_id_seq_pos[id_read_2][0][pos_read_tab:pos_read_tab + taille_kmer]) # ici on récupère le reverse complement du kmer de la séquence forward
        else:
            kmer_test = tab_id_seq_pos[id_read][0][pos_read_tab:pos_read_tab + taille_kmer]  # kmer_test récupère le premier kmer du read dans matrice = kmer à tester

        if kmer_test in tab_premiers_kmers:  # compare avec le dictionnaire des premiers kmers des reads

            for id_read_for_extend in tab_premiers_kmers[kmer_test]:  # id_read_for_extend parcours tab_premiers_kmers[test] pour tester tous les reads qui ont le même premier kmer

                if id_read_for_extend != id_read and id_read_for_extend not in res_copy['path'] and id_read_for_extend not in result['all_path']:  # on vérifie que id_read_for_extend ne retombe pas sur lui-même

                    if id_read.startswith('-'): # si reverse complement, on adapte chaque séquence en récuperant le forward et en le reverse complementant
                        seq_read_chevauchant = rev_comp(tab_id_seq_pos[id_read][0][pos_read_tab:])
                    else:
                        seq_read_chevauchant = tab_id_seq_pos[id_read][0][pos_read_tab:] # partie chevauchante dans le read que l'on veut étendre

                    if id_read_for_extend.startswith('-'):
                        id_read_for_extend_2 =id_read_for_extend[1:]
                        seq_read_extend_chevauchant = rev_comp(tab_id_seq_pos[id_read_for_extend_2][0][0:len(seq_read_chevauchant)])
                        seq_extension = rev_comp(tab_id_seq_pos[id_read_for_extend_2][0][len(seq_read_chevauchant):])
                    else:
                        seq_read_extend_chevauchant= tab_id_seq_pos[id_read_for_extend][0][0:len(seq_read_chevauchant)] # partie chevauchante dans le read qui pourrait faire l'extension
                        seq_extension = tab_id_seq_pos[id_read_for_extend][0][len(seq_read_extend_chevauchant):]# extention = partie non chevauchante du read utilisé pour l'extension

                    # on vérifie si le reste de la partie chevauchante est la même
                    if seq_read_chevauchant == seq_read_extend_chevauchant:

                        if len(seq_extension) != 0:  # vérifie que l'extension apportée par le read rallonge bien la séquence
                            res_copy['seq'] += seq_extension  # ajout de la partie non chevauchante dans result['seq']
                            res_tmp = extend(id_read_for_extend, tab_premiers_kmers, tab_id_seq_pos, taille_kmer, stop,res_copy)  # récursion
                            if res_tmp in not None:
                                return res_tmp
    return None

'''défninition de la fonction qui permet de remplir une matrice qui contient tous les premiers kmer de chaque read
  de la forme {seq: id reads qui commencent par cette sequence}.
  De plus, elle permet également de chercher dans chaque read la présence d'un kmer start.
  Si ce dernier est présent, on ajoute sa position à la matrice contenant tous les reads : 
  { id reads : sequence, position du start} 
  
'''
def fill_tab_premier_kmers(tab_id_seq,taille_kmer,start):

    tab_premiers_kmers = {}  # matrice contant le 1er kmer de tous les reads
# permet de parcourir les reads et de stocker le premier kmer dans un dictionnaire
    for id_read in tab_id_seq:  # remplissage de tab_premiers_kmers
        premier_kmer = tab_id_seq[id_read][0][0:taille_kmer]  # premier_kmer prend le premier kmer dans matrice
        rev_premier_kmer= rev_comp(premier_kmer)
        rev_id_read ='-'+id_read
        if  premier_kmer not in tab_premiers_kmers and premier_kmer != '':  # vérifie que premier_kmer n'est pas dans tab_premiers_kmers et qu'il n'est pas vide
            tab_premiers_kmers[premier_kmer] = [id_read]
        else:
            if premier_kmer != '':
                tab_premiers_kmers[premier_kmer].append(id_read)
        if rev_premier_kmer not in tab_premiers_kmers and rev_premier_kmer != '':  # vérifie que le premier_kmer du reverse complément n'est pas dans tab_premiers_kmers et qu'il n'est pas vide
            tab_premiers_kmers[rev_premier_kmer] = ['-' + id_read]  #
        else:
            if premier_kmer != '':
                tab_premiers_kmers[rev_premier_kmer].append('-'+id_read)
       #on profite du parcours de l'ensemble de tab_id_seq pour chercher dans les read, la présence du kmer start dans le sens forward et reverse
        pos_start = try_stop_start(tab_id_seq, id_read, start)
        if pos_start != -1 :
            tab_id_seq[id_read].append(pos_start)
        rev_pos_start = try_stop_start(tab_id_seq, rev_id_read, start)
        if rev_pos_start is not None:
            tab_id_seq[id_read].append(rev_pos_start)
    return tab_premiers_kmers

''' fonction qui permet de réaliser l'extension forward sans mismatch et d'écrire un fichier texte " result.txt "
en sortie contenant les paths et les séquences obtenues '''

def all_extentions(file_reads, file_start_stop,k):
    # chargement des donnees des fichiers / arguments en ligne de commande
    tab_id_seq_pos = parserMultiFASTA(file_reads)
    start, stop = parser_start_stop(file_start_stop)
    taille_kmer = k  # taille du kmer utilisée pour l'assemblage

    # création de la table contenant les premiers kmer
    tab_premiers_kmers = fill_tab_premier_kmers(tab_id_seq_pos,taille_kmer,start)

    # création du dictionnaire qui va nous permettre de stocker les chemins utilisés pour l'assemblage
    # et les séquence étendue pour chaque read
    result = {'path': [], 'all_path': [], 'seq': '', 'all_seq': []}
    all_result = [] # permet de stoker les resultats pour chaque nouveau start
    for i in tab_id_seq_pos: # parcours de la table contenant tous les reads
        if len(tab_id_seq_pos[i]) > 1:# s'il y a un read start
            result['path'] = [] # remise à zero pour recommencer à partir d'un nouveau start
            pos = tab_id_seq_pos[i][1]
            result['seq'] = tab_id_seq_pos[i][0][pos:]  # initialisation de la séquence avec le premier read de la
            # matrice contenant un start
            all_result.append(extend(i, tab_premiers_kmers, tab_id_seq_pos, taille_kmer, stop, result))
    #ecriture du fichier de sortie contenant les resultats
    exit_file = open("result.txt", "w")
    for i in all_result:
        if i != None:
            exit_file.write('path :'+ str(i['all_path'])+'\n')
            exit_file.write('extended sequence : '+ str(i['all_seq'])+'\n')
            exit_file.write('\n')
    exit_file.close()


''' récupération des arguments passés en ligne de commande dans l'ordre qui suit :
- chemin complet du fichier fasta contenant les read à assembler
- chemin complet du fichier fasta contenant les start et stop
- taille des kmers pour l'assemblage 
- option pour l'option a utiliser pour l'assemblage qui doit être parmis les suivantes : 
        -p , --perfect ''' # TODO ajouter l'option pour mismatch

path_file_read = sys.argv[1] # TODO warning if name contain spaces + full path
path_file_start_stop = sys.argv[2]

try:
    kmer_len = int(sys.argv[3])
except ValueError:
    print('WARNING kmer length is not a valid number')
option = sys.argv[4]

# assemblage parfait dans le sens forward
if option == '-p' or option == '--perfect':
    all_extentions(path_file_read,path_file_start_stop,kmer_len)
else :
    print('WARNING unknown option')


'''test temps et memoire'''
print("--- %s seconds ---" % (time.time() - start_time))
usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
print("--- %s kilobytes used ---" % usage)
