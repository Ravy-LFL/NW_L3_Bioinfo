#!/usr/local/bin/python3
# -*- coding: utf-8 -*-

import numpy as np
import sys 
import os

#mode = sys.argv[1]
mismatchX=-1
mismatch=1
match=2
gapouverture=-10
gapextensif=-1


#Ici la matrice de substitution blosum62 pour les alignement protéique
blosum62 = {
    ('W', 'F'): 1, ('L', 'R'): -2, ('S', 'P'): -1, ('V', 'T'): 0,('N','P'):-2,('N','K'):0,
    ('Q', 'Q'): 5, ('N', 'A'): -2, ('Z', 'Y'): -2, ('W', 'R'): -3,('K','N'):0,
    ('Q', 'A'): -1, ('S', 'D'): 0, ('H', 'H'): 8, ('S', 'H'): -1,
    ('H', 'D'): -1, ('L', 'N'): -3, ('W', 'A'): -3, ('Y', 'M'): -1,
    ('G', 'R'): -2, ('Y', 'I'): -1, ('Y', 'E'): -2, ('B', 'Y'): -3,
    ('Y', 'A'): -2, ('V', 'D'): -3, ('B', 'S'): 0, ('Y', 'Y'): 7,('K','P'):-1,
    ('G', 'N'): 0, ('E', 'C'): -4, ('Y', 'Q'): -1, ('Z', 'Z'): 4,
    ('V', 'A'): 0, ('C', 'C'): 9, ('M', 'R'): -1, ('V', 'E'): -2,
    ('T', 'N'): 0, ('P', 'P'): 7, ('V', 'I'): 3, ('V', 'S'): -2,
    ('Z', 'P'): -1, ('V', 'M'): 1, ('T', 'F'): -2, ('V', 'Q'): -2,
    ('K', 'K'): 5, ('P', 'D'): -1, ('I', 'H'): -3, ('I', 'D'): -3,
    ('T', 'R'): -1, ('P', 'L'): -3, ('K', 'G'): -2, ('M', 'N'): -2,
    ('P', 'H'): -2, ('F', 'Q'): -3, ('Z', 'G'): -2, ('X', 'L'): -1,
    ('T', 'M'): -1, ('Z', 'C'): -3, ('X', 'H'): -1, ('D', 'R'): -2,
    ('B', 'W'): -4, ('X', 'D'): -1, ('Z', 'K'): 1, ('F', 'A'): -2,
    ('Z', 'W'): -3, ('F', 'E'): -3, ('D', 'N'): 1, ('B', 'K'): 0,
    ('X', 'X'): -1, ('F', 'I'): 0, ('B', 'G'): -1, ('X', 'T'): 0,
    ('F', 'M'): 0, ('B', 'C'): -3, ('Z', 'I'): -3, ('Z', 'V'): -2,
    ('S', 'S'): 4, ('L', 'Q'): -2, ('W', 'E'): -3, ('Q', 'R'): 1,
    ('N', 'N'): 6, ('W', 'M'): -1, ('Q', 'C'): -3, ('W', 'I'): -3,
    ('S', 'C'): -1, ('L', 'A'): -1, ('S', 'G'): 0, ('L', 'E'): -3,
    ('W', 'Q'): -2, ('H', 'G'): -2, ('S', 'K'): 0, ('Q', 'N'): 0,
    ('N', 'R'): 0, ('H', 'C'): -3, ('Y', 'N'): -2, ('G', 'Q'): -2,
    ('Y', 'F'): 3, ('C', 'A'): 0, ('V', 'L'): 1, ('G', 'E'): -2,
    ('G', 'A'): 0, ('K', 'R'): 2, ('E', 'D'): 2, ('Y', 'R'): -2,('N','T'):0,('T','K'):-1,
    ('M', 'Q'): 0, ('T', 'I'): -1, ('C', 'D'): -3, ('V', 'F'): -1,('F','V'):-1,('F','K'):-3,
    ('T', 'A'): 0, ('T', 'P'): -1, ('B', 'P'): -2,('T', 'E'): -1,('B','B'):4,('F','G'):-3,('G','F'):-3,('F','D'):-3,('D','F'):-3}

 

#Matrice score en local
def MatriceScoreLocal(s1,s2):
	tailleligne2=len(s2)
	tailleligne1=len(s1)
	matrice=np.array([[0]*(tailleligne2+1)]*(tailleligne1+1))
	for i in range(tailleligne2+1) :  #Ces boucles initialises les scores associées aux gaps de la premiere lignes et de la premiere colonne
		matrice[0][i]=0
	for i in range(tailleligne1+1) :
		matrice[i][0]=0
	print(matrice)
	return(matrice)
	
#Lecture du fichier fasta ( nous ressort la taille nottamment)
def Lecture(sequence1, sequence2) :
        f1 = open(sequence1, "r")
        f2 = open(sequence2, "r")
        ligne=""
        s1=""
        s2=""
        taille1=0
        taille2=0
        for i in f1 :
                ligne=f1.readline()
                if ligne[1] != ">" and ligne!="" :
                        taille1+=len(ligne)
                        s1+=ligne.strip()
        for i in f2 :
                ligne=f2.readline()
                if ligne[1] != ">" and ligne!="" :
                        taille2+=len(ligne)
                        s2+=ligne.strip()
        f1.close()
        f2.close()
        return(taille1,taille2,s1,s2)

#initialisation de la Matrice de score ( sa taille est determiné par la taille des sequences et son remplissage est determiné par la fonction score )
def Matricescore(s1,s2) :
        tailleligne2=len(s2)
        tailleligne1=len(s1)
        matrice=np.array([[0]*(tailleligne2+1)]*(tailleligne1+1)) #On a la matrice de bonne taille
        matrice[0][0] = 0
        for i in range(tailleligne2+1) :  #Ces boucles initialises les scores associées aux gaps de la premiere lignes et de la premiere colonne
                matrice[0][i]= -i -9
        for i in range(tailleligne1+1) :
                matrice[i][0]= -i -9
        matrice[0][0]=0
        print(matrice)
        return(matrice)

#initialisation de la matrice TraceBack
def Matricetraceback(tailleseq1,tailleseq2,s1,s2) :
	matrice=np.array([[None]*(tailleseq2+1)]*(tailleseq1+1))
	for i in range(tailleseq2-1):
		matrice[0][i+2]=s2[i]      #Chaque case de la premiere ligne est remplie par la séquence 2
	for i in range(tailleseq1-1):
		matrice[i+2][0]=s1[i]      #Chaque case de la premire colonne est remplie par la premiere séquence
	for i in range(tailleseq2+1):
		matrice[1][i]="←"          #Affecte une fleche vers la gauche a chaque case de la seconde ligne
	for i in range(tailleseq1+1):
		matrice[i][1]="↑"          #Affecte une fleche vers le haut a chaque case de la seconde colonne
	matrice[1][1]="Done"
	matrice[0][1]=None
	matrice[1][0]=None
	print(matrice)
	return(matrice)


def DiagProt(character1,character2,MS): #Cette fonction calcul le score en diagonal en allant chercher le score de correspondance entre deux protéines dans la matrice de substitution 
	print("dans diagprot le character1 est :"+str(character1))
	print("dans diagprot le character2 est :"+str(character2))
	score=0
	match=blosum62[str(character1),str(character2)]
	score=match+MS[i-1][j-1]
	return(score)

def ScoreDIAG(charactere1,charactere2,MS):      #Trouve le score diag et le met dans la matrice
	score=0
	minuscule1=charactere1.lower()
	minuscule2=charactere2.lower()
	if minuscule1==minuscule2 :     #Test si c'est un match
		score=match+MS[i-1][j-1]
	else :
		if minuscule1=="a" and minuscule2=="t" or minuscule1=="t" and minuscule2=="a" or minuscule1=="c" and minuscule2=="g" or minuscule1=="g" and minuscule2=="c" :
			score=mismatch+MS[i-1][j-1]                 #Test si c'est un mismatch pyrimidine/pyrimidine ou un mismatch purine/purine et attribue le score en conséquen
		else :
			score = mismatchX+MS[i-1][j-1]                      #Ici c'est un mismatch pyrimidine/purine
	return(score)

def ScoreUp(MS,MT,i,j): #Trouve le score UP et remplie la matrice de score
        score=0
        if i==0 :       #Si i vaut 0 c'est qu'on est sur la premiere ligne alors c'est forcement un gap d'ouverture
                score=gapouverture+MS[i-1][j]
        else :          #On n'est pas sur la premiere ligne
                if MT[i-1][j]=="↑" :
                        score=MS[i-1][j]+gapextensif      #La case au dessus etait deja un gap
                else :
                        score=MS[i-1][j]+gapouverture     #La case au dessus n'etait pas un gap
        return(score)

def Scoreleft(MS,MT,i,j):       #Trouve le score left et le remplie dans la matrice SC
        score=0
        if j==0:
                score=gapouverture+MS[i][j-1]
        else :
                if MT[i][j-1]=="←" :
                        score=MS[i][j-1]+gapextensif      #La case a gauche etait deja un gap
                else :
                        score=MS[i][j-1]+gapouverture     #La case a gauche n'etait pas un gap
        return(score)





def bestscoreLocal(MS,MT,character1,character2,i,j):	#calcul les scores et remplie les matrices pour un alignement local en remplacant tout les scores négatif par des zéros
	scoreup=ScoreUp(MS,MT,i,j)
	if scoreup < 0 :
		scoreup=0
	print("up :"+str(scoreup))
	scoreleft=Scoreleft(MS,MT,i,j)
	if scoreleft<0 :
		scoreleft=0
	print("left :"+str(scoreleft))
	scorediag=ScoreDIAG(character1,character2,MS)
	if scorediag < 0 :
		scorediag=0
	print("diag :"+str(scorediag))
	symbole=""
	if scorediag>scoreleft and scorediag>scoreleft :
		MT[i+2][j+2]="↖"
		MS[i+1][j+1]=scorediag
	elif scoreleft>scorediag and scoreleft>scoreup :
		MT[i+2][j+2]="←"
		MS[i+1][j+1]=scoreleft
	elif scoreup>scorediag and scoreup>scoreleft :
		MT[i+2][j+2]="↑"
		MS[i+1][j+1]=scoreup
	elif scorediag==scoreleft:
		MT[i+2][j+2]="DL"
		MS[i+1][j+1]=scorediag
	elif scorediag==scoreup:
		MT[i+2][j+2]="DU"
		MS[i+1][j+1]=scorediag
	elif scoreleft==scoreup:
		MT[i+2][j+2]="LU"
		MS[i+1][j+1]=scoreleft
	elif scoreleft==scoreup and scoreup==scorediag :
		MT[i+2][j+2]="3"
		MS[i+1][j+1]=scoreup

def ChercheLeMax(MS): #Cette fonction permet de retourner les coordonées du plus grands score dans la matrice de score
	maxi=0
	for i in range(len(s1)):
		for j in range(len(s2)):
			if MS[i][j]>maxi :
				maxi=MS[i][j]
				coordoneI=i
				coordoneJ=j
	print("la coordone I :"+str(coordoneI))
	print("la coordone J :"+str(coordoneJ))
	print("le maximum :"+str(maxi))
	return(maxi,coordoneI,coordoneJ)

def bestscore(MS,MT,character1,character2,i,j): #Fonction qui determine le meilleur score
	print("dans bestscore le character1 est :"+str(character1))
	print("dans bestscore le character2 est :"+str(character2))
	scoreup=ScoreUp(MS,MT,i,j)      #Attribue le score en up
	print("up :"+str(scoreup))
	if mode=="-n":
		scorediag=ScoreDIAG(character1,character2,MS)   #Attribue le score en diag
		print("diag :"+str(scorediag))
	elif mode=="-p":
		scorediag=DiagProt(character1,character2,MS)
		print("diag :"+str(scorediag))
	scoreleft=Scoreleft(MS,MT,i,j)  #Attribue le score en left
	print("left :"+str(scoreleft))	#Les differents print precedant permettent de verifier si c'est les resultat sur la matrice de TraceBack son concordant
	scoreMAX=dict()
	if scorediag>scoreleft and scorediag>scoreup :  #Ici on verifie si le score diagonal est le plus grand
		scoreMAX["origine"]="↖"
		print("diagmax")
	elif scoreleft>scorediag and scoreleft>scoreup :        #On verifie si le score left est le plus grand 
		scoreMAX["origine"]="←"
		print("diagmax")
	elif scoreup>scoreleft and scoreup>scorediag :          #On verifie si le score up est le plus grand
		scoreMAX["origine"]="↑"
	elif scorediag==scoreleft :             #Dans ce cas la valeur de left et diag sont les memes
		scoreMAX["origine"]="DL"
	elif scorediag==scoreup :               #Dans ce cas les valeurs de up et diag sont les memes
		scoreMAX["origine"]="DU"
	elif scoreleft==scoreup :               #Dans ce cas les valeurs de up et left sont les memes
                scoreMAX["origine"]="LU"
	elif scorediag==scoreleft and scorediag==scoreup :      #Dans ce cas les 3 valeurs sont les memes
		scoreMAX["origine"]="3"
	else :
		print("erreur")
		print("score en diagonale :"+scorediag)
		print("score en up :"+scoreup)
		print("score en left :"+scoreleft)
	scoreMAX["Valeur"]=max(scorediag,scoreleft,scoreup)    #La valeur la plus haute est associé à la clef Valeur
	print(scoreMAX["Valeur"])		#cet ligne imprime la valeur choisis
	return(scoreMAX)


def AligneLocal(MS,MT,coordoneI,coordoneJ):	#Fonction qui permet de construire les séquences comparées et leur séquence qualité associées
	i=coordoneI
	j=coordoneJ
	alignement={}
	alignement["s1"]=""
	alignement["s2"]=""
	alignement["seqQuali"]=""
	while MS[i][j]!=0:	#A la difference de l'alignement global on s'arrete quand on croise un zero
		print(MT[i+1][j+1])
		if MT[i+1][j+1]=="↖" or MT[i+1][j+1]=="3" or MT[i+1][j+1]=="DL" or MT[i+1][j+1]=="DU" :
			alignement["s2"]+=MT[0][j+1]
			alignement["s1"]+=MT[i+1][0]
			if MT[0][j+1].lower()==MT[i+1][0].lower():
				alignement["seqQuali"]+="|"
			else:
				alignement["seqQuali"]+=":"
			j-=1
			i-=1
		elif MT[i+1][j+1]=="↑" or MT[i+1][j+1]=="LU":
			alignement["s1"]+=MT[i+1][0]
			alignement["s2"]+="_"
			alignement["seqQuali"]+=" "
			i-=1
		elif MT[i+1][j+1]=="←":
			alignement["s1"]+="_"
			alignement["s2"]+=MT[0][j+1]
			j-=1
	alignement["s1"]=alignement["s1"][::-1]
	alignement["s2"]=alignement["s2"][::-1]
	alignement["seqQuali"]=alignement["seqQuali"][::-1]
	return(alignement)

def Aligne(MS,MT):	#Effectue l'alignement
	i=MT.shape[0]-1	#ici les cases sont initialisé en bas a droite
	j=MT.shape[1]-1	#le dictionnaire contient la premiere sequence, la seconde sequence et la séquence qualité ( le score )
	alignement={}
	alignement["s1"]=""
	alignement["s2"]=""
	alignement["seqQuali"]=""
	while MT[i][j]!="Done":
		#print(MT[i][j])
		if MT[i][j]=="↖" or MT[i][j]=="3" or MT[i][j]=="DL" or MT[i][j]=="DU" :	#Si on croise une diagonale ( meme en DiagLeft ou DiagUp ) ca sera le chemin choisis
			alignement["s2"]+=MT[0][j]	#La s1 prend la lettre d'indice i de la sequence a gauche
			alignement["s1"]+=MT[i][0]	#La s2 prend la lettre d'indice j de la sequence en ligne
			if MT[0][j].lower()==MT[i][0].lower():	#Verifie si c'est un match, et si ca l'est on associe le symbole |
				alignement["seqQuali"]+="|"
			else:
				alignement["seqQuali"]+=":"
			j-=1				#On se déplace dans la case de diagonale haut/gauche
			i-=1
		elif MT[i][j]=="↑" or  MT[i][j]=="LU":	#Si on croise un Up c'est le chemin qu'on prendra
			alignement["s1"]+=MT[i][0]		
			alignement["s2"]+="_"
			alignement["seqQuali"]+=" "
			i-=1	#On se déplace vers la case de haut
		elif MT[i][j]=="←":      #La s1 prend la lettre d'indice i de la sequence a gauche
			alignement["s1"]+="_"
			alignement ["s2"]+=MT[0][j]
			alignement["seqQuali"]+=" "
			j-=1	#On pars vers la case a gauche
	alignement["s1"]=alignement["s1"][::-1]
	alignement["s2"]=alignement["s2"][::-1]
	alignement["seqQuali"]=alignement["seqQuali"][::-1]
	return(alignement)

def countmatch(symbAlign):	#On compte le nombre de match
	match=0
	for i in range(len(symbAlign)) :
		if symbAlign[i]=="|":
			match+=1
	return match

def countmissmatch(symbAlign):	#On compte le nombre de missmatch
	missmatch=0
	for i in range(len(symbAlign)):
		if symbAlign[i]==":":
			missmatch+=1
	return missmatch

def countgap(symbAlign):	#On compte le nombre de gap
	gap=0
	for i in range(len(symbAlign)):
		if symbAlign[i]==" ":
			gap+=1
	return gap


#Ce if permet de verifier si un argument a bien été choisis et permet de rappeler a l'utilisateur quels sont ses possibilités
if len(sys.argv) != 2 :
	sys.exit("ERREUR VEUILLEZ CHOISIR UN ARGUMENT '-p' POUR UN ALIGNEMENT PROTEIQUE OU '-n' POUR UN ALIGNEMENT NUCLEOTIDIQUE OU '-l' POUR UN ALIGNEMENT LOCAL")
else :
	mode = sys.argv[1]



if mode == "-n" :
#Programme principale pour nucléotide
	seq1=""
	seq2=""
	print("saisissez le nom du premier fichier fasta  :")
	seq1=input()
	extension1=os.path.splitext(seq1)
	if ".fa" not in extension1[1]:
		while ".fa" not in extension1[1] :
			print("le premier fichier n'est pas un fichier fasta, veuillez en saisir un nouveau :")
			seq1=input()
			extension1=os.path.splitext(seq1)
	print("saisissez le nom du second fichier fasta :")
	seq2=input()
	extension2=os.path.splitext(seq2)
	if ".fa" not in extension2[1]:
		while ".fa" not in extension2[1] :
			print("le second fichier n'est pas un fichier fasta, veuillez en saisir un nouveau :")
			seq2=input()
			extension2=os.path.splitext(seq2)
	taille1,taille2,s1,s2=Lecture(seq1,seq2)
	matriceSC=(Matricescore(s1,s2))
	matriceTB=Matricetraceback(taille1,taille2,s1,s2)
	for i in range(len(s1)):
        	for j in range(len(s2)):
                	matriceSC[i+1][j+1]=bestscore(matriceSC,matriceTB, s1[i],s2[j], i, j)["Valeur"]
                	matriceTB[i+2][j+2]=bestscore(matriceSC, matriceTB, s1[i], s2[j], i, j)["origine"]
	alignement=Aligne(matriceSC,matriceTB)	#On recupere l'alignement
	nbMatchGapMissmatch={"Nombre de match ":countmatch(alignement["seqQuali"]),"Nombre de missmatch ":countmissmatch(alignement["seqQuali"]),"Nombre de gap ":countgap(alignement["seqQuali"])}

	print(matriceSC)
	print(matriceTB)
	print("Alignement des séquences (avec alignement de la qualité):")
	print(alignement["s1"])
	print(alignement["seqQuali"])
	print(alignement["s2"])
	print("\n")
	print(nbMatchGapMissmatch)



elif mode == "-p" :
#Programme principale pour peptide
	seq1=""
	seq2=""
	print("saisissez le nom du premier fichier fasta  :")
	seq1=input()
	extension1=os.path.splitext(seq1)
	if ".fa" not in extension1[1]:
		while ".fa" not in extension1[1] :
			print("le premier fichier n'est pas un fichier fasta, veuillez en saisir un nouveau :")
			seq1=input()
			extension1=os.path.splitext(seq1)
	print("saisissez le nom du second fichier fasta :")
	seq2=input()
	extension2=os.path.splitext(seq2)
	if ".fa" not in extension2[1]:
		while ".fa" not in extension2[1] :
			print("le second fichier n'est pas un fichier fasta, veuillez en saisir un nouveau :")
			seq2=input()
			extension2=os.path.splitext(seq2)
	taille1,taille2,s1,s2=Lecture(seq1,seq2)
	matriceSC=(Matricescore(s1,s2))
	matriceTB=Matricetraceback(taille1,taille2,s1,s2)
	print(s1)
	print(s2)
	for i in range(len(s1)):
		for j in range(len(s2)):
			print("le carac s1 est :"+str(s1[i]))
			print("le carac s2 est :"+str(s2[j]))
			matriceSC[i+1][j+1]=bestscore(matriceSC,matriceTB,s1[i],s2[j],i,j)["Valeur"]
			matriceTB[i+2][j+2]=bestscore(matriceSC,matriceTB,s1[i],s2[j],i,j)["origine"]
			print(matriceSC)
			print(matriceTB)
	alignement=Aligne(matriceSC,matriceTB)
	nbMatchGapMissmatch={"Nombre de match ":countmatch(alignement["seqQuali"]),"Nombre de missmatch ":countmissmatch(alignement["seqQuali"]),"Nombre de gap ":countgap(alignement["seqQuali"])}

	print(matriceSC)
	print(matriceTB)
	print("Alignement des séquences proteiques (avec alignement de la qualité):")
	print(alignement["s1"])
	print(alignement["seqQuali"])
	print(alignement["s2"])
	print("\n")
	print(nbMatchGapMissmatch)

if mode == "-l" :
	#Programme principale local
	seq1=""
	seq2=""
	print("saisissez le nom du premier fichier fasta  :")
	seq1=input()
	extension1=os.path.splitext(seq1)
	if ".fa" not in extension1[1] :
		while ".fa" not in extension1[1] :
			print("le premier fichier n'est pas un fichier fasta, veuillez en saisir un nouveau :")
			seq1=input()
			extension1=os.path.splitext(seq1)
	print("saisissez le nom du second fichier fasta :")
	seq2=input()
	extension2=os.path.splitext(seq2)
	if ".fa" not in extension2[1]:
		while ".fa" not in extension2[1] :
			print("le second fichier n'est pas un fichier fasta, veuillez en saisir un nouveau :")
			seq2=input()
			extension2=os.path.splitext(seq2)
	taille1,taille2,s1,s2=Lecture(seq1,seq2)
	matriceSC=(MatriceScoreLocal(s1,s2))
	matriceTB=Matricetraceback(taille1,taille2,s1,s2)
	print(s1)
	print(s2)
	coordone1=0
	coordone2=0
	maxi=0
	for i in range(len(s1)):
		for j in range(len(s2)):
			bestscoreLocal(matriceSC,matriceTB,s1[i],s2[j],i,j)
		print(matriceSC)
		print(matriceTB)
	maxi,coordoneI,coordoneJ=ChercheLeMax(matriceSC)
	#print("VERIFIONS LE MAX")
	#print("voici le max en sortie de boucle"+str(maxi))
	#print("voici la coordoneI en sortie de boucle:"+str(coordoneI))
	#print("voici la coordoneJ en sortie de boucle:"+str(coordoneJ))
	#print("voici avec la coordooneI et J la sortie de la matriceSC: "+str(matriceSC[coordoneI][coordoneJ]))
	#print("voici avec la coordone I et J la sortie de la matriceTB :"+str(matriceTB[coordoneI][coordoneJ]))
	alignement=AligneLocal(matriceSC,matriceTB,coordoneI,coordoneJ)
	nbMatchGapMissmatch={"Nombre de match ":countmatch(alignement["seqQuali"]),"Nombre de missmatch ":countmissmatch(alignement["seqQuali"]),"Nombre de gap ":countgap(alignement["seqQuali"])}
	
	print(matriceSC)
	print(matriceTB)
	print("Alignement localement des séquences (avec alignement de la qualité):")
	print(alignement["s1"])
	print(alignement["seqQuali"])
	print(alignement["s2"])
	print("\n")
	print(nbMatchGapMissmatch)
