# ########################################
#              Basic Seq Object
##########################################

from Bio.Seq import Seq;

seq = Seq("AGTACACTGGT"); # Create a seqeunce since has methods attached to this obj

# ###########################
#    Similiar string methods
##############################

# Print sequence in list with index
# Enumerate associates index with the letter
for index, letter in enumerate(seq): 
    print("%i %s" % (index, letter)); 

# print length 
print("Length of sequence: ");
print("------------------------------")
print(len(seq));

# Access letters in string
print("Count letters in string:")
print("------------------------------")
print(seq[0]);
print("\n")

# Count substrings
print("Count substrings with seq.count(substring)):")
print("------------------------------")
print(seq.count("AC"));
print("\n")

# Percentage GC count
from Bio.SeqUtils import GC
gcCount = 100 * float(seq.count("G") + seq.count("C")) / len(seq);
gcCountMethod = GC(seq); # There is a built in function to do the above
print("Get the GC ratio with GC(seq):")
print("------------------------------")
print(gcCount,gcCountMethod);
print("\n")

# Slice sequence
slicedSeq = seq[4:12];
print("slice the sequence with seq[4:12]:")
print("------------------------------")
print(slicedSeq);
print("\n")

# Slice with a stride
print("Using stride, get the first and every 3rd codon, where seq[start:end:stride]")
print("------------------------------")
print(seq[0::3]);
print("\n")

print("Reverse the string using negative stride seq[::-1]");
print("------------------------------")
print(seq[::-1]);
print("\n")

print("Parse to str() or use directly in format string \">Name\\n%s\\n\" % seq");
print("------------------------------")
format_string = ">Name\n%s\n" % seq;
print(format_string);
print("\n")

# Concatenating sequences (+, for loop a list, .join)
seq1 = Seq("AGTCG");
seq2 = Seq("ACCCT");

print("Concatenation with +")
print("------------------------------")
print(seq1+seq2)
print("\n")
seq3 = Seq("ACCAG");
seqList = [seq1,seq2,seq3]

print("Loop through list and concatenate")
print("------------------------------")
concatSeq = Seq("");
for seq in seqList:
    concatSeq+=seq;

print(concatSeq)
print("\n")

print("Use .join() to concatenate strings and add spacer")
print("------------------------------")
spacer = Seq("N"*10)
print(spacer.join(seqList))
print("\n")

# Changing the case with .upper() and lower()
print("Changing the case with .upper() and lower()")
print("------------------------------")
print("\n")

# complement and reverse_complement of sequence
print("seq.complement() and seq.reverse_complement() of sequence")
print("------------------------------")
print("\n")

# The actual biological transcription process works from the template strand, doing a reverse complement
# (TCAG ! CUGA) to give the mRNA. However, in Biopython and bioinformatics in general, we typically
# work directly with the coding strand because this means we can get the mRNA sequence just by switching
# T ! U.
print("Transcribe the DNA sequence into RNA from the coding strand using transcribe()")
print("------------------------------")
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG");
template_dna = coding_dna.reverse_complement(); # Get the complement strand

mRna = coding_dna.transcribe();

print(mRna)
print("\n")

print("Translate the mRna into a protein using translate(). Look at other parameters it takes like start/stop codon, table type needs to be set")
print("------------------------------")

# The translation tables available in Biopython are based on those from the NCBI (see the next section of
# this tutorial). By default, translation will use the standard genetic code (NCBI table id 1). Suppose we are
# dealing with a mitochondrial sequence. We need to tell the translation function to use the relevant genetic
# code instead
protein = mRna.translate(table="Vertebrate Mitochondrial");

print(protein)
print("\n")

print("Translation tables")
print("------------------------------")
from Bio.Data import CodonTable
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
print(standard_table)