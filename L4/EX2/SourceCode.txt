"""
Download from NCBI the FASTA files containing the COVID-19 genome and the Influenza genome.
Use AI to compare the codon frequencies between the two genomes.

a) Make a chart that shows the top 10 most frequent codons for COVID-19;
b) Same as previous but for influenza;
c) Compare the two results and show the most frequent codons shared by the two genomes;
d) Show in the output of the console the top three amino acids from each genome;
e) Formulate a prompt for the AI such that the three amino acids are used to ask
   which foods contain less of those amino acids.
"""

# --- Automatically ensure BioPython is available ---
import sys
import subprocess
import importlib

try:
    from Bio import SeqIO
except ImportError:
    print("BioPython not found — installing it now...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "biopython"])
    importlib.invalidate_caches()
    from Bio import SeqIO

from collections import Counter
import matplotlib.pyplot as plt
import os

base_path = r"E:\Facultate\BIOINF\EX2"
flu_path = os.path.join(base_path, "sequence1.fasta")       # Influenza
covid_path = os.path.join(base_path, "sequence2.fasta")     # SARS-CoV-2


for path in [flu_path, covid_path]:
    if not os.path.exists(path):
        raise FileNotFoundError(f"❌ File not found: {path}")

print("✅ FASTA files found, starting analysis...")


def get_codons(seq):
    seq = str(seq).upper().replace("U", "T")
    return [seq[i:i+3] for i in range(0, len(seq) - 2, 3) if len(seq[i:i+3]) == 3]

def codon_frequency(file):
    record = next(SeqIO.parse(file, "fasta"))
    return Counter(get_codons(record.seq))

def codon_to_amino_acid(codon):
    table = {
        'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
        'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
        'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
        'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
        'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
        'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
        'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
        'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
    }
    return table.get(codon, '?')

def top_amino_acids(codon_counts, n=3):
    aa_counts = Counter()
    for codon, freq in codon_counts.items():
        aa = codon_to_amino_acid(codon)
        if aa != '*':
            aa_counts[aa] += freq
    return aa_counts.most_common(n)

covid_codons = codon_frequency(covid_path)
flu_codons = codon_frequency(flu_path)


covid_top10 = covid_codons.most_common(10)
plt.figure(figsize=(8, 5))
plt.bar([c for c, _ in covid_top10], [f for _, f in covid_top10], color='orange')
plt.title("Top 10 Most Frequent Codons - SARS-CoV-2")
plt.xlabel("Codon")
plt.ylabel("Frequency")
plt.tight_layout()
plt.savefig(os.path.join(base_path, "covid_top10_codons.png"))
plt.show()

flu_top10 = flu_codons.most_common(10)
plt.figure(figsize=(8, 5))
plt.bar([c for c, _ in flu_top10], [f for _, f in flu_top10], color='skyblue')
plt.title("Top 10 Most Frequent Codons - Influenza")
plt.xlabel("Codon")
plt.ylabel("Frequency")
plt.tight_layout()
plt.savefig(os.path.join(base_path, "flu_top10_codons.png"))
plt.show()

common_codons = set(c for c, _ in covid_top10) & set(c for c, _ in flu_top10)
print("\nCommon Top Codons between SARS-CoV-2 and Influenza:", common_codons)


labels = sorted(set(c for c, _ in covid_top10 + flu_top10))
covid_vals = [covid_codons[c] for c in labels]
flu_vals = [flu_codons[c] for c in labels]

plt.figure(figsize=(10, 6))
x = range(len(labels))
plt.bar(x, covid_vals, width=0.4, label='SARS-CoV-2', color='orange', align='center')
plt.bar([i + 0.4 for i in x], flu_vals, width=0.4, label='Influenza', color='skyblue', align='center')
plt.xticks([i + 0.2 for i in x], labels)
plt.ylabel("Frequency")
plt.title("Codon Frequency Comparison (Top 10 Sets Combined)")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(base_path, "codon_comparison.png"))
plt.show()


covid_topAA = top_amino_acids(covid_codons)
flu_topAA = top_amino_acids(flu_codons)

print("\nTop 3 amino acids in SARS-CoV-2:")
for aa, freq in covid_topAA:
    print(f"  {aa}: {freq}")

print("\nTop 3 amino acids in Influenza:")
for aa, freq in flu_topAA:
    print(f"  {aa}: {freq}")

covid_names = ", ".join([aa for aa, _ in covid_topAA])
flu_names = ", ".join([aa for aa, _ in flu_topAA])



