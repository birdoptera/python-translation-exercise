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
    framenum = range(1,4)
    for i in framenum:
    
        def translate_with_open_reading_frames(seq, framenum):
            open = FALSE
            translation = ""
#        print(seq)
            seqlength = len(seq) - (framenum-1) 
            for n in range(framenum-1, seqlength- (seqlength % 3), 3):
                codon = translate_RNA_codon(seq[n:n+3])
                open = (open or codon == "M") and not (codon == "*")
                translation += codon if open else "*"
            return translation


#    rna_sequence.upper()
#    if "AUG" in rna_sequence:
#        start = rna_sequence.find("AUG")
#        return translate_sequence(rna_sequence[start:], genetic_code)
#    else:
#        return ''
#    pass

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
