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

if __name__ == '__main__':
    unittest.main()
