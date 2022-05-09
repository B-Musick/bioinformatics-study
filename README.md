# Chapter 2 - General Info About Biopython

## Create a Sequence Object
<code>
from Bio.Seq import Seq;

seq = Seq("AGTACACTGGT"); # Create a seqeunce since has seq methods attached to this obj
</code>

## Import Sequence from FASTA File

### File From NCBI
- Go to <a href="https://www.ncbi.nlm.nih.gov/nuccore/">NCBI Website</a> to locate nucleotide sequences.
- Info about <a href="https://biopython.org/wiki/SeqIO">SeqIO</a>
### Import sequence and parse into list
<code>
from Bio import SeqIO     

records = list(SeqIO.parse("ls_orchid.fasta", "fasta"));

print(records[0].id)  # first record  
print(records[-1].id)  # last record
</code>

## Import Sequence from GenBank File

### Import sequence and parse into list
<code>
from Bio import SeqIO     

records = list(SeqIO.parse("ls_orchid.gbk", "genbank"))
print(records[0].id)  # first record  
print(records[-1].id)  # last record
</code>




# Chapter 9 
<a href="https://www.ncbi.nlm.nih.gov/Web/Search/entrezfs.html">Entrez</a></br>
<a href="https://pubmed.ncbi.nlm.nih.gov/">PubMed</a>

# Chapter 20
<a href="https://biopython.org/wiki/Category:Cookbook">Project Cookbook</a>

## Project
<a href="https://biopython.org/wiki/Phylo_cookbook">Bio.Phylo</a></br>
<a href="https://biopython.org/wiki/Phylo">Phylogenetic Trees</a>



