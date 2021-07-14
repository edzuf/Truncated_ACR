import sys, os
from Bio import SeqIO

blastn_file = sys.argv[1]
spacer_fasta_file = sys.argv[2]
ref_fasta_file = sys.argv[3]
sample = sys.argv[4]

f = open(blastn_file, "r")

spacers = SeqIO.index(spacer_fasta_file, "fasta")
references = SeqIO.index(ref_fasta_file, "fasta")
plasmid_seq = references["pNZ123.xdna"].seq

dico_info = {}

# lecture du fichier blast pour trouver les spacers qui alignent sur le plasmide ac 0 mismatch
for line in f:
        list_line = line.split("\t")
        spacer_ID = list_line[0]
        ref_ID = list_line[1]
        if (ref_ID == "pNZ123.xdna"):
		spacer_length = len(spacers[spacer_ID].seq)
		al_length = int(list_line[3])
        	al_mismatch = int(list_line[4])
        	real_mismatch = spacer_length - al_length + al_mismatch
        	if (real_mismatch <= 0):
                	# on extrait alors le PAM
			start_plasmid = int(list_line[8])
			end_plasmid = int(list_line[9])
			if (start_plasmid < end_plasmid):
				amont = plasmid_seq[start_plasmid-8:start_plasmid-1]
				PAM = plasmid_seq[end_plasmid:end_plasmid+7]
				start_PAM = end_plasmid+1
				end_PAM = end_plasmid+7
				list_to_add = [start_PAM, end_PAM, str(PAM), str(spacers[spacer_ID].seq), str(amont)]
				dico_info[spacer_ID] = list_to_add
			else:
				amont = plasmid_seq[start_plasmid:start_plasmid+7].reverse_complement()
				PAM = plasmid_seq[end_plasmid-8:end_plasmid-1].reverse_complement()
				start_PAM = end_plasmid-1
                                end_PAM = end_plasmid-7
                                list_to_add = [start_PAM, end_PAM, str(PAM), str(spacers[spacer_ID].seq), str(amont)]
                                dico_info[spacer_ID] = list_to_add
			#for element in dico_info.keys():
				#print(dico_info[element])

f.close()
#print(len(dico_info))

# ecriture des PAMs sous format fasta
f_PAM = open("{}_seq_amont.fasta".format(sample), "w")

for spacer in dico_info.keys():
	f_PAM.write(">{}\n".format(spacer))
	f_PAM.write(dico_info[spacer][4])
	f_PAM.write("\n")

f_PAM.close()
