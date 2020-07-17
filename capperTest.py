import logging
import math
import unittest
from capper import Read

class CapperTest(unittest.TestCase):
    def test_minimizer_interval_iterator(self):
        reads = Read.parse_reads("minimizers_incorrect_1kg") #, maxReads=1000)

        for read in reads:
            #print(read)
            i, j, k = 0, -1, 0

            for left, right, bottom, top in read.minimizer_interval_iterator(read.minimizers, start_fn=lambda x: x.minimizer_start, length=29):
                #print(" Left {}, Right {}, Bottom {}, Top {}".format(left, right, bottom, top))
                assert left < right
                assert i <= left
                assert i < right
                i = right

                assert j == bottom or j+1 == bottom
                assert bottom < top
                j = bottom

                assert k == top or k + 1 == top
                k = top

                # Check each minimizer is part of the interval
                for minimizer in read.minimizers[bottom:top]:
                    assert minimizer.minimizer_start <= left
                    assert minimizer.minimizer_start + len(minimizer.minimizer) >= right

                # Check each preceding minimizer occurs before the interval
                for minimizer in read.minimizers[:bottom]:
                    assert minimizer.minimizer_start + len(minimizer.minimizer) <= left

                # Check each proceding minimizer occurs after the interval
                for minimizer in read.minimizers[top:]:
                    assert minimizer.minimizer_start >= right

            assert i <= len(read.read)
            assert k == len(read.minimizers)
            
    def test_recompute_adam_cap_nocap(self):
        
        # Here's a version where no minimizers are located (last fields all 0)
        parts = ["TTCTGCTGTTGTCAGAAGGCTACTTCTGATATTTATGCTTTCATTCTCCATTAATCATTTAATATTTTATTTTCCTTCATCAGCCAACACTGCAGGCAGTTGAGGCTTTTTCCATGATGAAGTAAATCCAACTATGATTTCTGAAGAATA",
                 "????????+??????????????????++?????????????????????????????????????+??????????????????&5??????5?????????????????????????????????????+5?????????5???????",
                 "1",
                 "AAAATAAAATATTAAATGATTAATGGAGA", "44", "41", "39", "1", "0",
                 "GAAATCATAGTTGGATTTACTTCATCATG", "112", "107", "43", "1", "0", 
                 "ATAGTTGGATTTACTTCATCATGGAAAAA", "106", "96", "49", "1", "0",
                 "TGGATTTACTTCATCATGGAAAAAGCCTC", "101", "94", "40", "1", "0",
                 "AGCCTCAACTGCCTGCAGTGTTGGCTGAT", "78", "68", "49", "1", "0",
                 "TCAACTGCCTGCAGTGTTGGCTGATGAAG", "74", "64", "42", "1", "0",
                 "AGGCAGTTGAGGCTTTTTCCATGATGAAG", "93", "83", "49", "1", "0",
                 "CAGGCAGTTGAGGCTTTTTCCATGATGAA", "92", "82", "39", "1", "0",
                 "AACACTGCAGGCAGTTGAGGCTTTTTCCA", "85", "79", "41", "1", "0",
                 "GATGAAGGAAAATAAAATATTAAATGATT", "52", "42", "49", "1", "0",
                 "TTGTCAGAAGGCTACTTCTGATATTTATG", "8", "0", "47", "1", "0",
                 "TTTTATTTTCCTTCATCAGCCAACACTGC", "64", "54", "48", "1", "0",
                 "CATTTAATATTTTATTTTCCTTCATCAGC", "55", "53", "39", "1", "0",
                 "ATGAAAGCATAAATATCAGAAGTAGCCTT", "15", "9", "44", "1", "0",
                 "TCATTCTCCATTAATCATTTAATATTTTA", "40", "35", "44", "1", "0",
                 "ATGCTTTCATTCTCCATTAATCATTTAAT", "34", "30", "43", "1", "0",
                 "TATTTATGCTTTCATTCTCCATTAATCAT", "29", "27", "41", "1", "0",
                 "TGATATTTATGCTTTCATTCTCCATTAAT", "26", "16", "49", "1", "0",
                 "CTGATATTTATGCTTTCATTCTCCATTAA", "25", "15", "39", "1", "0",
                 "1.07374e+09", "134.858", "0", "winner"]
    
        serialized = '\t'.join(parts)
        
        read = Read(serialized, True)
        
        # When no minimizers are explored we complain with inf
        assert read.recompute_adam_cap() == float('inf')
        
    def test_recompute_adam_cap_onemin(self):
    
        # + = q10, ? = q30
        # Here's a version where only one minimizer is located
        # Best way to hit it will be at the one q10 base.
        #         --------TTGTCAGAAGGCTACTTCTGATATTTATG----------
        parts = ["TTCTGCTGTTGTCAGAAGGCTACTTCTGATATTTATGCTTTCATTCTCCATTAATCATTTAATATTTTATTTTCCTTCATCAGCCAACACTGCAGGCAGTTGAGGCTTTTTCCATGATGAAGTAAATCCAACTATGATTTCTGAAGAATA",
                 "???????????????????????+??????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????",
                 "1",
                 "AAAATAAAATATTAAATGATTAATGGAGA", "44", "41", "39", "1", "0",
                 "GAAATCATAGTTGGATTTACTTCATCATG", "112", "107", "43", "1", "0", 
                 "ATAGTTGGATTTACTTCATCATGGAAAAA", "106", "96", "49", "1", "0",
                 "TGGATTTACTTCATCATGGAAAAAGCCTC", "101", "94", "40", "1", "0",
                 "AGCCTCAACTGCCTGCAGTGTTGGCTGAT", "78", "68", "49", "1", "0",
                 "TCAACTGCCTGCAGTGTTGGCTGATGAAG", "74", "64", "42", "1", "0",
                 "AGGCAGTTGAGGCTTTTTCCATGATGAAG", "93", "83", "49", "1", "0",
                 "CAGGCAGTTGAGGCTTTTTCCATGATGAA", "92", "82", "39", "1", "0",
                 "AACACTGCAGGCAGTTGAGGCTTTTTCCA", "85", "79", "41", "1", "0",
                 "GATGAAGGAAAATAAAATATTAAATGATT", "52", "42", "49", "1", "0",
                 "TTGTCAGAAGGCTACTTCTGATATTTATG", "8", "0", "47", "1", "1",
                 "TTTTATTTTCCTTCATCAGCCAACACTGC", "64", "54", "48", "1", "0",
                 "CATTTAATATTTTATTTTCCTTCATCAGC", "55", "53", "39", "1", "0",
                 "ATGAAAGCATAAATATCAGAAGTAGCCTT", "15", "9", "44", "1", "0",
                 "TCATTCTCCATTAATCATTTAATATTTTA", "40", "35", "44", "1", "0",
                 "ATGCTTTCATTCTCCATTAATCATTTAAT", "34", "30", "43", "1", "0",
                 "TATTTATGCTTTCATTCTCCATTAATCAT", "29", "27", "41", "1", "0",
                 "TGATATTTATGCTTTCATTCTCCATTAAT", "26", "16", "49", "1", "0",
                 "CTGATATTTATGCTTTCATTCTCCATTAA", "25", "15", "39", "1", "0",
                 "1.07374e+09", "134.858", "0", "winner"]
    
        serialized = '\t'.join(parts)
        
        read = Read(serialized, True)
        
        recomputed_cap = read.recompute_adam_cap()
        
        # We should be able to just pay the cost of the cheapest base.
        assert recomputed_cap == 10, "Cap is {} when it shouldn't be".format(recomputed_cap)
        
        
    def test_recompute_adam_cap_overlap(self):
    
        # + = q10, 0 = q15, ? = q30, 
        # Here's a version where only one minimizer is located
        # Best way to hit it will be at the shared q15 base.
        #         --------TTGTCAGAAGGCTACTTCTGATATTTATG----------
        #                        ----------CTGATATTTATGCTTTCATTCTCCATTAA
        parts = ["TTCTGCTGTTGTCAGAAGGCTACTTCTGATATTTATGCTTTCATTCTCCATTAATCATTTAATATTTTATTTTCCTTCATCAGCCAACACTGCAGGCAGTTGAGGCTTTTTCCATGATGAAGTAAATCCAACTATGATTTCTGAAGAATA",
                 "???????????+????????????????0????????????????????+????????????????????????????????????????????????????????????????????????????????????????????????????",
                 "1",
                 "AAAATAAAATATTAAATGATTAATGGAGA", "44", "41", "39", "1", "0",
                 "GAAATCATAGTTGGATTTACTTCATCATG", "112", "107", "43", "1", "0", 
                 "ATAGTTGGATTTACTTCATCATGGAAAAA", "106", "96", "49", "1", "0",
                 "TGGATTTACTTCATCATGGAAAAAGCCTC", "101", "94", "40", "1", "0",
                 "AGCCTCAACTGCCTGCAGTGTTGGCTGAT", "78", "68", "49", "1", "0",
                 "TCAACTGCCTGCAGTGTTGGCTGATGAAG", "74", "64", "42", "1", "0",
                 "AGGCAGTTGAGGCTTTTTCCATGATGAAG", "93", "83", "49", "1", "0",
                 "CAGGCAGTTGAGGCTTTTTCCATGATGAA", "92", "82", "39", "1", "0",
                 "AACACTGCAGGCAGTTGAGGCTTTTTCCA", "85", "79", "41", "1", "0",
                 "GATGAAGGAAAATAAAATATTAAATGATT", "52", "42", "49", "1", "0",
                 "TTGTCAGAAGGCTACTTCTGATATTTATG", "8", "0", "47", "1", "1",
                 "TTTTATTTTCCTTCATCAGCCAACACTGC", "64", "54", "48", "1", "0",
                 "CATTTAATATTTTATTTTCCTTCATCAGC", "55", "53", "39", "1", "0",
                 "ATGAAAGCATAAATATCAGAAGTAGCCTT", "15", "9", "44", "1", "0",
                 "TCATTCTCCATTAATCATTTAATATTTTA", "40", "35", "44", "1", "0",
                 "ATGCTTTCATTCTCCATTAATCATTTAAT", "34", "30", "43", "1", "0",
                 "TATTTATGCTTTCATTCTCCATTAATCAT", "29", "27", "41", "1", "0",
                 "TGATATTTATGCTTTCATTCTCCATTAAT", "26", "16", "49", "1", "0",
                 "CTGATATTTATGCTTTCATTCTCCATTAA", "25", "15", "39", "1", "1",
                 "1.07374e+09", "134.858", "0", "winner"]
    
        serialized = '\t'.join(parts)
        
        read = Read(serialized, True)
        
        recomputed_cap = read.recompute_adam_cap()
        
        assert recomputed_cap == 15, "Cap is {} when it shouldn't be".format(recomputed_cap)
        
    def test_recompute_adam_cap_flank(self):
    
        # + = q10, 0 = q15, ? = q30, 
        # Here's a version where only one minimizer is located
        # Best way to hit it will be at the shared q10 base and hope we break minimizer 2.
        # Minimizer 2 is beat with probability 0.02730720856393431 by any of 7 candidates in the agglomeration
        #         --------TTGTCAGAAGGCTACTTCTGATATTTATG----------                 
        #                        ----------CTGATATTTATGCTTTCATTCTCCATTAA             
        #         000000000011111111112222222222333333333344444444445555555555
        #         012345678901234567890123456789012345678901234567890123456789
        parts = ["TTCTGCTGTTGTCAGAAGGCTACTTCTGATATTTATGCTTTCATTCTCCATTAATCATTTAATATTTTATTTTCCTTCATCAGCCAACACTGCAGGCAGTTGAGGCTTTTTCCATGATGAAGTAAATCCAACTATGATTTCTGAAGAATA",
                 "?????????????????????+????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????",
                 "1",
                 "AAAATAAAATATTAAATGATTAATGGAGA", "44", "41", "39", "1", "0",
                 "GAAATCATAGTTGGATTTACTTCATCATG", "112", "107", "43", "1", "0", 
                 "ATAGTTGGATTTACTTCATCATGGAAAAA", "106", "96", "49", "1", "0",
                 "TGGATTTACTTCATCATGGAAAAAGCCTC", "101", "94", "40", "1", "0",
                 "AGCCTCAACTGCCTGCAGTGTTGGCTGAT", "78", "68", "49", "1", "0",
                 "TCAACTGCCTGCAGTGTTGGCTGATGAAG", "74", "64", "42", "1", "0",
                 "AGGCAGTTGAGGCTTTTTCCATGATGAAG", "93", "83", "49", "1", "0",
                 "CAGGCAGTTGAGGCTTTTTCCATGATGAA", "92", "82", "39", "1", "0",
                 "AACACTGCAGGCAGTTGAGGCTTTTTCCA", "85", "79", "41", "1", "0",
                 "GATGAAGGAAAATAAAATATTAAATGATT", "52", "42", "49", "1", "0",
                 "TTGTCAGAAGGCTACTTCTGATATTTATG", "8", "0", "47", "1", "1",
                 "TTTTATTTTCCTTCATCAGCCAACACTGC", "64", "54", "48", "1", "0",
                 "CATTTAATATTTTATTTTCCTTCATCAGC", "55", "53", "39", "1", "0",
                 "ATGAAAGCATAAATATCAGAAGTAGCCTT", "15", "9", "44", "1", "0",
                 "TCATTCTCCATTAATCATTTAATATTTTA", "40", "35", "44", "1", "0",
                 "ATGCTTTCATTCTCCATTAATCATTTAAT", "34", "30", "43", "1", "0",
                 "TATTTATGCTTTCATTCTCCATTAATCAT", "29", "27", "41", "1", "0",
                 "TGATATTTATGCTTTCATTCTCCATTAAT", "26", "16", "49", "1", "0",
                 "CTGATATTTATGCTTTCATTCTCCATTAA", "25", "15", "39", "1", "1",
                 "1.07374e+09", "134.858", "0", "winner"]
    
        serialized = '\t'.join(parts)
        
        read = Read(serialized, True)
        
        recomputed_cap = read.recompute_adam_cap()
        right_answer = 10 + -10 * math.log10(1 - ((1 - 0.02730720856393431)**7))
        
        assert recomputed_cap == right_answer, "Cap is {} when it should be {}".format(recomputed_cap, right_answer)
    
    
if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    unittest.main()
