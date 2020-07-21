import logging
import math
import unittest
from capper import Read

class CapperTest(unittest.TestCase):
    def test_overlap_interval_iterator(self):
        reads = Read.parse_reads("new/minimizers_incorrect_1kg") #, maxReads=1000)

        for read in reads:
            #print(read)
            i, j, k = 0, -1, 0

            for left, right, bottom, top in read.overlap_interval_iterator(read.minimizers,
                                            start_fn=lambda x: x.minimizer_start, length_fn=lambda x : 29):
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
        
    @unittest.skip("Cap is way off in this read due to minimizer probability quantization in vg")
    def test_recompute_adam_cap_hard(self):
    
        # This read gets the cap way off somehow
        # TODO: work it out

        #                                                                                                    11111111111111111111111111111111111111111111111111
        #          11111111112222222222333333333344444444445555555555666666666677777777778888888888999999999900000000001111111111222222222233333333334444444444
        #012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
        #TGCGGGCTGGCAATCGGCTCGCCCTCCTTCCAAGACAAGCTGCTGCTCGGCCCCTCTGAGATCCGTGTGCCTTCGGTGGTGCCTGTGCAGACGGTGCGATAGCTAGTCGCGCTCCTGTACAGCGGTTCCCTCGTTGTGGCACAGAGTGAA
        #?????????????????????????????????????+????????????????????????????????????????????????????????????'+5????????+???+++????????????+?????????????????????
        #                                  !                              !
        #
        #------CTGGCAATCGGCTCGCCCTCCTTCCAAGA----------                                                                                                         
        #       ------TCGGCTCGCCCTCCTTCCAAGACAAGCTG----------                                                                                                  
        #              ----------TCCTTCCAAGACAAGCTGCTGCTCGGCCC---------                                                                                        
        #                        ----------ACAAGCTGCTGCTCGGCCCCTCTGAGATC--                                                                                     
        #                           ----------AGCTGCTGCTCGGCCCCTCTGAGATCCGT----------                                                                          
        #                                      --------TCGGCCCCTCTGAGATCCGTGTGCCTTCG-------                                                                    
        #                                            ----------TCTGAGATCCGTGTGCCTTCGGTGGTGCC---------                                                          
        #                                                      ----------GTGTGCCTTCGGTGGTGCCTGTGCAGACG----------
        
        # These aren't disrupted:
        #                                                                 ---------GGTGGTGCCTGTGCAGACGGTGCGATAGC----------                                     
        #                                                                           --------TGTGCAGACGGTGCGATAGCTAGTCGCGC----------                            
        #                                                                                    -TGCAGACGGTGCGATAGCTAGTCGCGCTC----------                          
        #                                                                                      ----------CGATAGCTAGTCGCGCTCCTGTACAGCGG-----                    
        #                                                                                            ----------CTAGTCGCGCTCCTGTACAGCGGTTCCCT----------         
        #                                                                                                       ---TCGCGCTCCTGTACAGCGGTTCCCTCGTT----------     
        #                                                                                                           -----TCCTGTACAGCGGTTCCCTCGTTGTGGCA-----    
        #                                                                                                            ----------ACAGCGGTTCCCTCGTTGTGGCACAGAGT-  
        #                                                                                                              ----------AGCGGTTCCCTCGTTGTGGCACAGAGTGA-
        
        # Base 34 scores 39.209729
        
        # Base 37:
        #DEBUG:capper:At base 37
        #DEBUG:capper:Cost to break base: 10.000000
        #DEBUG:capper:Overlap agglomeration: window_start: 0 window_length: 7 minimizer_start: 6 minimizer_length: 29 minimizer: TCTTGGAAGGAGGGCGAGCCGATTGCCAG hits: 2 hits_in_extension: 2 hash: 000462e7b41b1caf beat_prob: 6.693035912203998e-05
        #DEBUG:capper:Out of minimizer: need to beat with probability 0.000067 and 8 chances for 32.713886 points
        #DEBUG:capper:Overlap agglomeration: window_start: 7 window_length: 7 minimizer_start: 13 minimizer_length: 29 minimizer: CAGCTTGTCTTGGAAGGAGGGCGAGCCGA hits: 2 hits_in_extension: 2 hash: 2b647149fc6d160d beat_prob: 0.16950138145732682
        #DEBUG:capper:In minimizer: no additional cost
        #DEBUG:capper:Overlap agglomeration: window_start: 14 window_length: 10 minimizer_start: 24 minimizer_length: 29 minimizer: GGGCCGAGCAGCAGCTTGTCTTGGAAGGA hits: 2 hits_in_extension: 2 hash: 118da2866b09a303 beat_prob: 0.06856742650692253
        #DEBUG:capper:In minimizer: no additional cost
        #DEBUG:capper:Overlap agglomeration: window_start: 24 window_length: 3 minimizer_start: 34 minimizer_length: 29 minimizer: GATCTCAGAGGGGCCGAGCAGCAGCTTGT hits: 2 hits_in_extension: 2 hash: 04cc30719faef05a beat_prob: 0.018740680446793116
        #DEBUG:capper:In minimizer: no additional cost
        #DEBUG:capper:Overlap agglomeration: window_start: 27 window_length: 11 minimizer_start: 37 minimizer_length: 29 minimizer: AGCTGCTGCTCGGCCCCTCTGAGATCCGT hits: 2 hits_in_extension: 2 hash: 040e7fa2a286d1f7 beat_prob: 0.015846230703142866
        #DEBUG:capper:In minimizer: no additional cost
        #DEBUG:capper:Total cost for best solution ending in this base: 42.713886
        
        # 0.000067 << 1/255 so our quantizing minimizer beat probabilities to 8 bits in the approximation is probably to blame for the mismatch in caps here.
                
        
        parts = ["TGCGGGCTGGCAATCGGCTCGCCCTCCTTCCAAGACAAGCTGCTGCTCGGCCCCTCTGAGATCCGTGTGCCTTCGGTGGTGCCTGTGCAGACGGTGCGATAGCTAGTCGCGCTCCTGTACAGCGGTTCCCTCGTTGTGGCACAGAGTGAA",
                 "?????????????????????????????????????+????????????????????????????????????????????????????????????'+5????????+???+++????????????+?????????????????????",
                 "2",
                 "AGCTGCTGCTCGGCCCCTCTGAGATCCGT", "37", "27", "49", "2", "2",
                 "GGCACCACCGAAGGCACACGGATCTCAGA", "54", "44", "48", "2", "2", 
                 "CGAAGGCACACGGATCTCAGAGGGGCCGA", "46", "38", "44", "2", "2", 
                 "GTGTGCCTTCGGTGGTGCCTGTGCAGACG", "64", "54", "49", "2", "2", 
                 "GATCTCAGAGGGGCCGAGCAGCAGCTTGT", "34", "24", "41", "2", "2", 
                 "GGGCCGAGCAGCAGCTTGTCTTGGAAGGA", "24", "14", "48", "2", "2", 
                 "CAGCTTGTCTTGGAAGGAGGGCGAGCCGA", "13", "7", "45", "2", "2", 
                 "TCTTGGAAGGAGGGCGAGCCGATTGCCAG", "6", "0", "45", "2", "2", 
                 "GGTGGTGCCTGTGCAGACGGTGCGATAGC", "74", "65", "48", "0", "0",
                 "CGATAGCTAGTCGCGCTCCTGTACAGCGG", "96", "86", "44", "0", "0", 
                 "TCGCGCTCCTGTACAGCGGTTCCCTCGTT", "106", "103", "42", "0", "0", 
                 "GCGCGACTAGCTATCGCACCGTCTGCACA", "83", "75", "47", "0", "0",
                 "TCCTGTACAGCGGTTCCCTCGTTGTGGCA", "112", "107", "39", "0", "0",
                 "GAGCGCGACTAGCTATCGCACCGTCTGCA", "85", "84", "40", "0", "0",
                 "ACAGCGGTTCCCTCGTTGTGGCACAGAGT", "118", "108", "40", "0", "0",
                 "AGGGAACCGCTGTACAGGAGCGCGACTAG", "102", "92", "49", "0", "0", 
                 "TCACTCTGTGCCACAACGAGGGAACCGCT", "120", "110", "40", "0", "0",
                 "1.50515", "58.0915", "0", "align"]
    
        serialized = '\t'.join(parts)
        
        read = Read(serialized, True)
        read.visualize()
        
        recomputed_cap = read.recompute_adam_cap()
        right_answer = 58.0915
        
        assert recomputed_cap == right_answer, "Cap is {} when it should be {}".format(recomputed_cap, right_answer)
    
if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    unittest.main()
