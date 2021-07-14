import sys, os
from Bio import SeqIO

blastn_file = sys.argv[1]
fasta_file = sys.argv[2]
tol_mismatch = int(sys.argv[3])

f = open(blastn_file, "r")

dico_spacer = {}
fasta_dico = SeqIO.index(sys.argv[2], "fasta")

# lecture du fichier new_spacers.blastn pour decider si bon hit ou non
for line in f:
        list_line = line.split("\t")
        spacer_ID = list_line[0]
        if spacer_ID not in dico_spacer:
		dico_spacer[spacer_ID] = []
	spacer_length = len(fasta_dico[spacer_ID].seq)
	ref_ID = list_line[1]
	al_length = int(list_line[3])
	al_mismatch = int(list_line[4])
	real_mismatch = spacer_length - al_length + al_mismatch
	if (real_mismatch <= tol_mismatch):
		dico_spacer[spacer_ID].append(ref_ID)

f.close()

dico_ref =  {"pNZ123":[], "DGCC7710":[], "AcrIIA5":[], "AcrIIA6":[], "Acr-phage123":[], "no_hit":[], "ambiguous":[]}

# ajout des spacers qui ne sont pas dans le blastn dans le dico avec pas de hit
for spacer in fasta_dico:
	if spacer not in dico_spacer:
		dico_ref["no_hit"].append(spacer)
#print(len(dico_ref["no_hit"]))
compteur = 0

# repartition des spacers selon l'endroit ou ils alignent
for spacer in dico_spacer:
	if len(dico_spacer[spacer]) == 0:
		dico_ref["no_hit"].append(spacer)
	elif len(dico_spacer[spacer]) > 1:
		dico_ref["ambiguous"].append(spacer)
	elif len(dico_spacer[spacer]) == 1:
		compteur += 1
		if dico_spacer[spacer][0].startswith("pNZ123"):
			dico_ref["pNZ123"].append(spacer)
                if dico_spacer[spacer][0].startswith("DGCC7710"):
                        dico_ref["DGCC7710"].append(spacer)
                if dico_spacer[spacer][0].startswith("AcrIIA5"):
                        dico_ref["AcrIIA5"].append(spacer)
                if dico_spacer[spacer][0].startswith("AcrIIA6"):
                        dico_ref["AcrIIA6"].append(spacer)
                if dico_spacer[spacer][0].startswith("Acr-phage123"):
                        dico_ref["Acr-phage123"].append(spacer)
print(compteur)
for key in dico_ref:
	print("{} : {} spacers".format(key, len(dico_ref[key])))

# extraction des spacers no hits
counter = 0
no_hit_spacer_count = {}

for element in dico_ref["no_hit"]:
	spacer = fasta_dico[element].seq
	#print(spacer)
	if spacer not in no_hit_spacer_count:
		no_hit_spacer_count[spacer] = 1
	else:
		no_hit_spacer_count[spacer] += 1
	counter += 1
print("Spacers with no hits: {}".format(counter))
ordered_list = sorted(no_hit_spacer_count.items(), key=lambda x: x[1], reverse=True)
#print(ordered_list)
print("Unique spacers with no hits: {}".format(len(ordered_list)))
for i in range(0, 5):
	print("{} ({})".format(ordered_list[i][0], ordered_list[i][1]))
