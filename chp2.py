# ########################################
#              Basic Seq Object
##########################################

from Bio.Seq import Seq;

seq = Seq("AGTACACTGGT"); # Create a seqeunce since has methods attached to this obj
# print("Initial Sequence:   "+seq);
# # Sequence Methods
# print("Complement:         "+seq.complement());
# print("Reverse Complement: "+seq.reverse_complement());

# ########################################
#                   SeqIO
##########################################
from Bio import SeqIO
# for seq_record in SeqIO.parse("ls_orchid.fasta","fasta"):
#     print(seq_record.id);
#     print(repr(seq_record.seq)); # Print represent
#     print(len(seq_record)); # Prints length of the string

# ###########################
#      Parse into list
##############################
records = list(SeqIO.parse("ls_orchid.fasta", "fasta"))
print(records[0].id)  # first record
print(records[-1].id)  # last record

recGenBank = list(SeqIO.parse("ls_orchid.gbk", "genbank"))
print(records[0].id)  # first record
print(records[-1].id)  # last record
for i in recGenBank:
    print(i.id);