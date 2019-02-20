#! /usr/bin/env python3

import unittest

import translate

class TestTranslateBaseClass(unittest.TestCase):
    def setUp(self):
        self.genetic_code = {'GUC': 'V', 'ACC': 'T', 'GUA': 'V', 'GUG': 'V', 'ACU': 'T', 'AAC': 'N', 'CCU': 'P', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'CAU': 'H', 'AAU': 'N', 'AGU': 'S', 'GUU': 'V', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'CCA': 'P', 'ACA': 'T', 'CCC': 'P', 'UGU': 'C', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'UGC': 'C', 'CAG': 'Q', 'GAU': 'D', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'AGG': 'R', 'GGG': 'G', 'UCC': 'S', 'UCA': 'S', 'UAA': '*', 'GGA': 'G', 'UAC': 'Y', 'GAC': 'D', 'UAG': '*', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'GGC': 'G', 'AUG': 'M', 'CUG': 'L', 'GAG': 'E', 'CUC': 'L', 'AGA': 'R', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'GCU': 'A', 'GAA': 'E', 'AUU': 'I', 'UUG': 'L', 'UUA': 'L', 'UGA': '*', 'UUC': 'F'}


    def run_translate_sequence(self, rna_seq,
            expected_amino_acid_seq,
            gen_code = None):
        if gen_code is None:
            gen_code = self.genetic_code
        amino_acid_seq = translate.translate_sequence(
                rna_sequence = rna_seq,
                genetic_code = gen_code)
        message = (
                "\n\n"
                "Calling `translate_sequence` with `rna_sequence`:\n"
                "\t{0!r}\n"
                "Expecting {1!r}, but {2!r} was returned".format(
                        rna_seq,
                        expected_amino_acid_seq,
                        amino_acid_seq))
        self.assertEqual(amino_acid_seq, expected_amino_acid_seq, message)

    def run_get_all_translations(self, rna_seq,
            expected_results,
            gen_code = None):
        if gen_code is None:
            gen_code = self.genetic_code
        amino_acid_seqs = translate.get_all_translations(
                rna_sequence = rna_seq,
                genetic_code = gen_code)
        message = (
                "\n\n"
                "Calling `run_get_all_translations` with `rna_sequence`:\n"
                "\t{0!r}\n"
                "Expecting {1!r}, but {2!r} was returned".format(
                        rna_seq,
                        expected_results,
                        amino_acid_seqs))
        if amino_acid_seqs:
            amino_acid_seqs = sorted(amino_acid_seqs)
        if expected_results:
            expected_results = sorted(expected_results)
        self.assertEqual(amino_acid_seqs, expected_results, message)

    def run_get_reverse(self, seq, expected_result):
        rev_seq = translate.get_reverse(seq)
        message = (
                "\n\n"
                "Calling `get_reverse` with argument:\n"
                "\t{0!r}\n"
                "Expecting\n"
                "\t{1!r}\n"
                "to be returned, but got\n"
                "\t{2!r}\n".format(
                        seq,
                        expected_result,
                        rev_seq))
        self.assertEqual(rev_seq, expected_result, message)








if __name__ == '__main__':
    unittest.main() 
