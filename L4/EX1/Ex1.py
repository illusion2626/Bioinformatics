""" Implement an application that converts the coding region of a gene into an amino acid sequence. Use the genetic code from from below:(on moodle)
Take a sequence of DNA and convert it into an ARN sequence S_DNA= ...ATG......TAA... ; S_ARN= ...AUG......UAA..."""


genetic_code = {
    'UUU':'Phe', 'UUC':'Phe', 'UUA':'Leu', 'UUG':'Leu',
    'CUU':'Leu', 'CUC':'Leu', 'CUA':'Leu', 'CUG':'Leu',
    'AUU':'Ile', 'AUC':'Ile', 'AUA':'Ile', 'AUG':'Met',
    'GUU':'Val', 'GUC':'Val', 'GUA':'Val', 'GUG':'Val',
    'UCU':'Ser', 'UCC':'Ser', 'UCA':'Ser', 'UCG':'Ser',
    'CCU':'Pro', 'CCC':'Pro', 'CCA':'Pro', 'CCG':'Pro',
    'ACU':'Thr', 'ACC':'Thr', 'ACA':'Thr', 'ACG':'Thr',
    'GCU':'Ala', 'GCC':'Ala', 'GCA':'Ala', 'GCG':'Ala',
    'UAU':'Tyr', 'UAC':'Tyr', 'UAA':'Stop', 'UAG':'Stop',
    'CAU':'His', 'CAC':'His', 'CAA':'Gln', 'CAG':'Gln',
    'AAU':'Asn', 'AAC':'Asn', 'AAA':'Lys', 'AAG':'Lys',
    'GAU':'Asp', 'GAC':'Asp', 'GAA':'Glu', 'GAG':'Glu',
    'UGU':'Cys', 'UGC':'Cys', 'UGA':'Stop', 'UGG':'Trp',
    'CGU':'Arg', 'CGC':'Arg', 'CGA':'Arg', 'CGG':'Arg',
    'AGU':'Ser', 'AGC':'Ser', 'AGA':'Arg', 'AGG':'Arg',
    'GGU':'Gly', 'GGC':'Gly', 'GGA':'Gly', 'GGG':'Gly'
}

def dna_to_mrna(dna_seq):

    return dna_seq.upper().replace('T', 'U')

def translate_rna_to_protein(rna_seq):
   
    protein = []
    for i in range(0, len(rna_seq) - 2, 3):
        codon = rna_seq[i:i+3]
        amino_acid = genetic_code.get(codon, '')
        if amino_acid == 'Stop' or amino_acid == '':
            break
        protein.append(amino_acid)
    return '-'.join(protein)


S_DNA = "CATTGATGCCATATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
S_RNA = dna_to_mrna(S_DNA)
protein = translate_rna_to_protein(S_RNA)

print("DNA: ", S_DNA)
print("mRNA:", S_RNA)
print("Protein:", protein)
