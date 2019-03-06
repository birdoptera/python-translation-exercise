#! /usr/bin/env python3

import sys


#def translate_RNA_codon(genetic_code):
#    return RNA_codon_table[genetic_code]

def translate_sequence(rna_sequence, genetic_code):
    def translate_RNA_codon(codon):
        return genetic_code[codon]
#    print(rna_sequence)
    translation = ''
    if len(rna_sequence) > 2:
        for n in range(0, len(rna_sequence)-(len(rna_sequence) % 3), 3):
            if "*" in translate_RNA_codon(rna_sequence.upper()[n:n+3]):
                return translation
            else:
                translation += translate_RNA_codon(rna_sequence.upper()[n:n+3])
    else: return ''

    return translation


    """Translates a sequence of RNA into a sequence of amino acids.

    Translates `rna_sequence` into string of amino acids, according to the
    `genetic_code` given as a dict. Translation begins at the first position of
    the `rna_sequence` and continues until the first stop codon is encountered
    or the end of `rna_sequence` is reached.

    If `rna_sequence` is less than 3 bases long, or starts with a stop codon,
    an empty string is returned.
    """

    pass

def get_all_translations(rna_sequence, genetic_code):
    """Get a list of all amino acid sequences encoded by an RNA sequence.

    All three reading frames of `rna_sequence` are scanned from 'left' to
    'right', and the generation of a sequence of amino acids is started
    whenever the start codon 'AUG' is found. The `rna_sequence` is assumed to
    be in the correct orientation (i.e., no reverse and/or complement of the
    sequence is explored).

    The function returns a list of all possible amino acid sequences that
    are encoded by `rna_sequence`.

    If no amino acids can be translated from `rna_sequence`, an empty list is
    returned.
    """
    def translate_RNA_codon(codon):
        return genetic_code[codon]
    seq = rna_sequence.upper()
    print(seq)
    # Get an empty list to store our translations in
    list_of_peptides = []
    # Figure out the index (position) of the last possible codon
    seq_length = len(seq)
    index_of_last_codon = seq_length - 3
    # Loop over all possible positions that could be the start of a start codon
    for base_index in range(index_of_last_codon + 1):
        # get next codon by slicing next three nucleotides
        codon = seq[base_index: base_index + 3]
        # check if next codon is a start codon
        if codon == "AUG":
            # Now we are ready to translate and store the result in our list.
            # Just inserting a print statement below, but you can replace this
            # with code to do the translation and IF an amino acid sequence is
            # returned, append it to 'list_of_peptides'
            print(seq[base_index:])
            peptide = translate_sequence(rna_sequence[base_index:], genetic_code)
            if peptide != '':
                list_of_peptides += peptide
    return list_of_peptides

def get_reverse(sequence):
    """Reverse orientation of `sequence`.

    Returns a string with `sequence` in the reverse order.

    If `sequence` is empty, and empty string is returned.
    """
    if sequence != '':

        new = sequence.upper()
        new = new[::-1]
        return new
    else:
        return ''
    pass

#    sequence = sequence.upper()
#    rev_seq = ""
#    for c in sequence:
#        rev_seq = c + rev_seq
#    return rev_seq

#    pass


def get_complement(sequence):
    """Get the complement of `sequence`.

    Returns a string with the complementary sequence of `sequence`.

    If `sequence` is empty, and empty string is returned.
    """
    new = sequence.upper()
    if sequence != '':
        new = new.replace('A', 'W')
        new = new.replace('G', 'X')
        new = new.replace('C', 'Y')
        new = new.replace('U', 'Z')
        new = new.replace('W', 'U')
        new = new.replace('X', 'C')
        new = new.replace('Y', 'G')
        new = new.replace('Z', 'A')
        return new

    else:
        return ''
    pass

# or
#    sequence = sequence.upper()
#    comp_bases = {
#        'A' : 'U',
#        'U' : 'A',
#        'G' : 'T',
#        'T' : 'G',
#    comp_seq = ""
#    for c in sequence:
#        comp_seq += comp_bases[c]
#    return comp_seq

def reverse_and_complement(sequence):
    """Get the reversed and complemented form of `sequence`.

    Returns a string that is the reversed and complemented sequence
    of `sequence`.

    If `sequence` is empty, and empty string is returned.
    """
    if sequence != '':

        new = sequence.upper()
        new = new[::-1]
        new = new.replace('A', 'W')
        new = new.replace('G', 'X')
        new = new.replace('C', 'Y')
        new = new.replace('U', 'Z')
        new = new.replace('W', 'U')
        new = new.replace('X', 'C')
        new = new.replace('Y', 'G')
        new = new.replace('Z', 'A')
        return new

    else:
        return ''
#   or
#   return get_reverse(get_complement(sequence))

    pass

def get_longest_peptide(rna_sequence, genetic_code):
    """Get the longest peptide encoded by an RNA sequence.

    Explore six reading frames of `rna_sequence` (three reading frames of the
    current orientation, and the reversed and complemented form) and return (as
    a string) the longest sequence of amino acids that it encodes, according to
    the `genetic_code`.

    If no amino acids can be translated from `rna_sequence` nor its reverse and
    complement, an empty list is returned.
    """
    def translate_RNA_codon(codon):
        return genetic_code[codon]
   
    sequence1 = rna_sequence.upper()

    sequence2 = rna_sequence.upper()
    sequence2 = sequence2[::-1]

    sequence3  = rna_sequence.upper()
    sequence3 = sequence3.replace('A', 'W')
    sequence3 = sequence3.replace('G', 'X')
    sequence3 = sequence3.replace('C', 'Y')
    sequence3 = sequence3.replace('U', 'Z')
    sequence3 = sequence3.replace('W', 'U')
    sequence3 = sequence3.replace('X', 'C')
    sequence3 = sequence3.replace('Y', 'G')
    sequence3 = sequence3.replace('Z', 'A')

    sequence4 = rna_sequence.upper()
    sequence4 = sequence4[::-1]

    sequence4 = sequence4.replace('A', 'W')
    sequence4 = sequence4.replace('G', 'X')
    sequence4 = sequence4.replace('C', 'Y')
    sequence4 = sequence4.replace('U', 'Z')
    sequence4 = sequence4.replace('W', 'U')
    sequence4 = sequence4.replace('X', 'C')
    sequence4 = sequence4.replace('Y', 'G')
    sequence4 = sequence4.replace('Z', 'A')

    sequences = [sequence1, sequence2, sequence3, sequence4]
    for seq in sequences:
        list_of_peptides = []
        seq_length = len(seq)
        index_of_last_codon = seq_length - 3
        for base_index in range(index_of_last_codon + 1):
            codon = seq[base_index: base_index + 3]
            if codon == "AUG":
                print(codon)
                peptide = translate_sequence(rna_sequence[base_index:], genetic_code)
                if peptide != '':
                    list_of_peptides += peptide

    longest_peptide = max(list_of_peptides, key=len)
    return longest_peptide


    pass


if __name__ == '__main__':
    genetic_code = {'GUC': 'V', 'ACC': 'T', 'GUA': 'V', 'GUG': 'V', 'ACU': 'T', 'AAC': 'N', 'CCU': 'P', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'CAU': 'H', 'AAU': 'N', 'AGU': 'S', 'GUU': 'V', 'CAC': 'H',
 'ACG': 'T', 'CCG': 'P', 'CCA': 'P', 'ACA': 'T', 'CCC': 'P', 'UGU': 'C', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'UGC': 'C', 'CAG': 'Q', 'GAU': 'D', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'AGG': 'R',
 'GGG': 'G', 'UCC': 'S', 'UCA': 'S', 'UAA': '*', 'GGA': 'G', 'UAC': 'Y', 'GAC': 'D', 'UAG': '*', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'GGC': 'G', 'AUG': 'M', 'CUG': 'L', 'GAG': 'E', 'CUC': 'L',
 'AGA': 'R', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'GCU': 'A', 'GAA': 'E', 'AUU': 'I', 'UUG': 'L', 'UUA': 'L',
 'UGA': '*', 'UUC': 'F'}
    rna_seq = ("AUG"
            "UAC"
            "UGG"
            "CAC"
            "GCU"
            "ACU"
            "GCU"
            "CCA"
            "UAU"
            "ACU"
            "CAC"
            "CAG"
            "AAU"
            "AUC"
            "AGU"
            "ACA"
            "GCG")
    longest_peptide = get_longest_peptide(rna_sequence = rna_seq,
            genetic_code = genetic_code)
    assert isinstance(longest_peptide, str), "Oops: the longest peptide is {0}, not a string".format(longest_peptide)
    message = "The longest peptide encoded by\n\t'{0}'\nis\n\t'{1}'\n".format(
            rna_seq,
            longest_peptide)
    sys.stdout.write(message)
    if longest_peptide == "MYWHATAPYTHQNISTA":
        sys.stdout.write("Indeed.\n")
