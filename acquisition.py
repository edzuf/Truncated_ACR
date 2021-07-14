import sys, os
from Bio import SeqIO

blastn_file = sys.argv[1]
sample = sys.argv[3]
tol_mismatch = int(sys.argv[4])

f = open(sys.argv[1], "r")

dico_counter = {}

# lecture du fichier blastn et comptage du nb de repeat par read
for line in f:
        list_line = line.split("\t")
        read_ID = list_line[0]
        if read_ID not in dico_counter:
                dico_counter[read_ID] = [0, []]
        al_mismatch = int(list_line[4])
        al_length = int(list_line[3])
        real_mismatch = 36 - al_length + al_mismatch
	# exclusion des reads avec un trop grand nombre de mismatchs
	if real_mismatch >= 4:
                dico_counter[read_ID][0] = -1
	# calcul des positions reelles start et end du repeat sur le read
	elif (real_mismatch <= tol_mismatch and dico_counter[read_ID][0] != -1):
		start_query = int(list_line[6])
        	end_query = int(list_line[7])
		start_repeat = int(list_line[8])
        	end_repeat = int(list_line[9])
		if start_repeat > end_repeat:
			intermed = end_repeat
			end_repeat = start_repeat
			start_repeat = intermed
		diff_pos_start = start_repeat - 1
		diff_pos_end = 36 - end_repeat
		real_start = start_query - diff_pos_start
		real_end = end_query + diff_pos_end
        	dico_counter[read_ID][0] += 1
		# ajout des positions start et end du repeat ds liste
                dico_counter[read_ID][1].append(real_start)
               	dico_counter[read_ID][1].append(real_end)

f.close()

# dictionnaires pour classement des reads selon nb de repeat
dico_rambiguous = {}
dico_r1 = {}
dico_r2 = {}
dico_r3 = {}
dico_r4 = {}

# enregistrement du contenu du fichier fasta dans un dictionnaire
record_dict = SeqIO.index(sys.argv[2], "fasta")

# liste de tous les nouveaux spacers extraits
spacer_list = []

# pour chaque read ac nouveau spacer, cle=read ID, value=spacer
new_spacer_dict = {}

# ajout des reads qui ne sont pas ds le fichier blastn dans le dict ac 0 repeat
for read in record_dict:
	if read not in dico_counter:
		dico_r1[read] = [0,[]]

# classement des reads dans le bon dictionnaire
for read in dico_counter:
	if dico_counter[read][0] == -1:
		dico_rambiguous[read] = dico_counter[read]
	elif (dico_counter[read][0] == 0 or dico_counter[read][0] == 1):
		dico_r1[read] = dico_counter[read]
	elif dico_counter[read][0] == 2:
		dico_r2[read] = dico_counter[read]
	elif dico_counter[read][0] == 3:
		dico_r3[read] = dico_counter[read]
		read_seq = record_dict[read]
		dico_counter[read][1].sort()
		# extraction du nouveau spacer
		spacer_start = dico_counter[read][1][1] + 1
                spacer_end = dico_counter[read][1][2] - 1
                spacer = read_seq.seq[spacer_start-1:spacer_end]
		spacer_list.append(str(spacer))
		new_spacer_dict[read] = str(spacer)
		#print(spacer)
	elif dico_counter[read][0] >= 4:
		dico_r4[read] = dico_counter[read]

spacer_set = set(spacer_list)

# affichage des nouveaux spacers sous format fasta
f2 = open("{}_new_spacers_5.fasta".format(sample), "w")

for cle in new_spacer_dict:
	f2.write(">{}\n".format(cle))
	f2.write(new_spacer_dict[cle])
	f2.write("\n")

f2.close() 

# calcul du pourcentage d'acquisition
#acquisition = float(len(dico_r3))/(len(dico_r2)+len(dico_r3))*100

print("Nombre de reads ambigus: {}".format(len(dico_rambiguous)))
print("Nombre de reads avec moins de 2 repeats: {}".format(len(dico_r1)))
print("Nombre de reads avec 2 repeats: {}".format(len(dico_r2)))
print("Nombre de reads avec 3 repeats: {}".format(len(dico_r3)))
print("Nombre de reads avec 4 repeats ou plus: {}".format(len(dico_r4)))

#print("Le pourcentage d'acquisition est de {}".format(acquisition))
#print("Il y a {} spacers differents".format(len(spacer_set)))
print("Il y a {} nouveaux spacers dans le dictionnaire".format(len(new_spacer_dict)))
