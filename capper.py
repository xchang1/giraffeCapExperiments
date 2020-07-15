import math
import matplotlib.pyplot as plt
import numpy as np
import time

# negative infinity
n_inf = float("-inf")

def log_add(i, j):
    l = lambda x: math.log10(math.pow(10, x) + 1)
    if i < j:
        return l(j - i) + i if j - i <= 10 else j
    return l(i - j) + j if i - j <= 10 else i

class Minimizer:
    """
    Represents a minimizer within a read.
    """

    total_window_length = 29 + 11 - 1 # Currently hard-coded, and therefore assumed the same for all minimizers

    def __init__(self, minimizer, minimizer_start, window_start, window_length, hits, hits_in_extension, read_length):
        self.minimizer, self.minimizer_start, self.window_start, self.window_length, self.hits, self.hits_in_extension = \
            minimizer, minimizer_start, window_start, window_length, hits, hits_in_extension
        # Sanity checks
        assert read_length >= 0
        assert minimizer_start >= 0 and minimizer_start + len(minimizer) <= read_length
        assert window_start >= 0 and window_start + window_length <= read_length
        assert window_start <= minimizer_start
        assert window_start + window_length >= minimizer_start + len(minimizer)

    def __str__(self):
        return "window_start: {} window_length: {} minimizer_start: {} minimizer_length: {} minimizer: {} hits: {} hits_in_extension: {}".\
            format(self.window_start, self.window_length-Minimizer.total_window_length+1, self.minimizer_start, len(self.minimizer), self.minimizer, self.hits, self.hits_in_extension)

class Read:
    """
    Represents a read and its minimizers.
    """

    def __init__(self, line, correct):
        # Format is tab sep:
        #READ_STRING QUAL_STRING (MINIMIZER_STRING START WINDOW_START WINDOW_LENGTH HITS HITS_IN_EXTENSIONS?)xN MAP_Q STAGE
        tokens = line.split()
        self.read, self.qual_string = tokens[0], [ ord(i) - ord("!") for i in tokens[1] ]
        assert len(self.read) == len(self.qual_string)
        #raw mapq, adam's mapq_extend_cap, my probability cluster lost cap,  last correct stage
        self.map_q, self.adam_cap, self.xian_cap, self.stage = float(tokens[-4]), float(tokens[-3]), float(tokens[-2]), tokens[-1]
        self.correct = correct
        self.minimizers = sorted([ Minimizer(tokens[i+3], int(tokens[i+4]), int(tokens[i+5]), int(tokens[i+6]), int(tokens[i+7]), int(tokens[i+8]), len(self.read)) for i in range(0, len(tokens)-7, 6) ], key=lambda x : x.window_start) # Sort by start coordinate
        # Check they don't overlap and every window is accounted for
        p_window = 0
        for minimizer in self.minimizers:
            assert p_window == minimizer.window_start
            p_window += minimizer.window_length - Minimizer.total_window_length + 1
        assert p_window + Minimizer.total_window_length == len(self.read)+1

    def __str__(self):
        return "Read map_q:{} adam_cap:{} xian_cap:{} fast_cap:{} \n read_string: {}\n qual_string: {}\n ".format(self.map_q, self.adam_cap, self.xian_cap, self.fast_cap(), self.read, " ".join(map(str, self.qual_string))) + "\n\t".join(map(str, self.minimizers)) + "\n"

    def minimizer_interval_iterator(self, minimizers, start_fn, length):
        """
        Iterator over common intervals of overlapping minimizers.
        :param minimizers: the list of minimizers to iterate over
        :param start_fn: returns the start of the minimizer in the minimizer
        :param length: the length of the minimizer
        :return: yields sequence of (left, right, bottom, top) tuples, where left if the first base
        of the minimizer interval (inclusive), right is the last base of the minimizer interval (exclusive),
        bottom is the index of the first minimizer in the interval and top is the index of the last minimizer
        in the interval (exclusive).
        """
        # Handle no minimizer case
        if len(minimizers) == 0:
            return

        stack = [minimizers[0]]  # Minimizers currently being iterated over
        left = [start_fn(minimizers[0])]  # The left end of a minimizer interval
        bottom = [0]  # The index of the first minimizer in the interval in the sequence of minimizers

        def get_preceding_intervals(right):
            # Get all intervals that precede a given point "right"
            while left[0] < right:

                # Case where the left-most minimizer ends before the start of the new minimizer
                if start_fn(stack[0]) + length <= right:
                    yield left[0], start_fn(stack[0]) + length, bottom[0], bottom[0] + len(stack)

                    # If the stack contains only one minimizer there is a gap between the minimizer
                    # and the new minimizer, otherwise just shift to the end of the leftmost minimizer
                    left[0] = right if len(stack) == 1 else start_fn(stack[0]) + length

                    bottom[0] += 1
                    del(stack[0])

                # Case where the left-most minimizer ends at or after the beginning of the new new minimizer
                else:
                    yield left[0], right, bottom[0], bottom[0] + len(stack)
                    left[0] = right

        # For each minimizer in turn
        for minimizer in minimizers[1:]:
            assert len(stack) > 0

            # For each new minimizer we return all intervals that
            # precede start_fn(minimizer)
            for interval in get_preceding_intervals(start_fn(minimizer)):
                yield interval

            stack.append(minimizer)  # Add the new minimizer for the next loop

        # Intervals of the remaining intervals on the stack
        for interval in get_preceding_intervals(len(self.read)):
            yield interval

    def get_log_prob_of_base_error_in_interval(self, left, right, sum_fn=log_add):
        """
        Gives the log10 prob of a base error in the given interval of the read
        :param left: Leftmost base in interval of read, inclusive
        :param right: Rightmost base in interval of read, exclusive
        :param sum_fn: Function for adding/maxing log probs in interval
        :return: log10 prob as float
        """
        assert left < right  # We don't allow zero length intervals
        p = -(self.qual_string[left]) / 10
        for i in range(left+1, right):
            p = sum_fn(-(self.qual_string[i]) / 10, p)
        return p

    def get_log_prob_of_minimizer_skip(self, minimizer):
        """ Get probability of the minimizer not being part of a considered extension.
        :param minimizer:
        :return: log10 prob as float
        """
        assert minimizer.hits_in_extension <= minimizer.hits
        if minimizer.hits == 0:  # Case there are no hits in graph so must skip
            return 0.0
        if minimizer.hits_in_extension == 0:  # Case no hits are in extensions so must skip
            return 0.0
        if minimizer.hits > 1:
            return 0.0
        return n_inf

        #if minimizer.hits == minimizer.hits_in_extension:  # Case all hits in extensions so can't skip
        #    return n_inf
        #return math.log10(1.0 - minimizer.hits_in_extension/minimizer.hits)

    def fast_cap(self, start_fn=lambda x: x.minimizer_start, minimizer_length=29,
                 sum_fn=log_add, minimizer_filter_fn=lambda x : True,
                 include_skip_prob=True):
        """ Compute cap on the reads qual score by considering possibility that base errors and unlocated
        minimizer hits prevented us finding the true alignment.

        Algorithm uses a "sweep line" dynamic programming approach.
        For a read with minimizers aligned to it:

                     000000000011111111112222222222
                     012345678901234567890123456789
        Read:        ******************************
        Minimizer 1:    *****
        Minimizer 2:       *****
        Minimizer 3:                   *****
        Minimizer 4:                      *****

        For each distinct read interval of overlapping minimizers, e.g. in the example
        the intervals 3,4,5; 6,7; 8,9,10; 18,19,20; 21,22; and 23,24,25
        we consider base errors that would result in the minimizers in the interval being incorrect

        We use dynamic programming sweeping left-to-right over the intervals to compute the probability of the minimum number
        of base errors needed to disrupt all the minimizers, or for the minimizers to have not been located within extensions.

        :param start_fn: Function that returns start of the minimizer
        :param minimizer_length: Length of the minimizer
        :param sum_fn: method for adding up probability of one or more base errors in an interval of bases
        :return: Phred scaled log prob of read being wrongly aligned
        """
        # This step filters minimizers we will definitely skip - needs to be coordinated with get_log_prob_of_minimizer_skip
        # function above
        minimizer_filter_fn = lambda x : x.hits == 1 and x.hits_in_extension == 1
        filtered_minimizers = list(filter(minimizer_filter_fn, self.minimizers))

        c = np.full(len(filtered_minimizers) + 1, n_inf)  # Log10 prob of having mutated minimizers,
        # such that c[i+1] is log prob of mutating minimizers 0, 1, 2, ..., i
        c[0] = 0.0

        pBottom = 0
        for left, right, bottom, top in self.minimizer_interval_iterator(filtered_minimizers, start_fn=start_fn, length=minimizer_length):
            # If a new minimizer in the interval include probability of skipping the minimizer
            if include_skip_prob and pBottom == bottom:
                p = c[bottom] + self.get_log_prob_of_minimizer_skip(filtered_minimizers[bottom])
                if c[bottom+1] < p:
                    c[bottom+1] = p
                pBottom += 1

            # Calculate prob of all intervals up to top being disrupted
            p = c[bottom] + self.get_log_prob_of_base_error_in_interval(left, right, sum_fn=sum_fn)

            # Replace min-prob for minimizers in the interval
            for i in range(bottom+1, top+1):
                if c[i] < p:
                    c[i] = p

        # Checks
        if include_skip_prob:
            assert pBottom == len(filtered_minimizers)
        assert c[-1] != n_inf

        return -c[-1] * 10

    @staticmethod
    def parse_reads(reads_file, correct=True, max_reads=-1):
        reads = []
        with open(reads_file) as fh: # This masks a bug
            for line in fh:
                try:
                    # This is try/except loop is here because some minimizers contain overlapping windows and we ignore
                    # those minimizers right now
                    reads.append(Read(line, correct))
                except:
                    pass
                if max_reads != -1 and len(reads) > max_reads:
                    break
        return reads

class Reads:
    def __init__(self, correct_reads_file, incorrect_reads_file, max_correct_reads=-1, max_incorrect_reads=-1):
        self.reads = Read.parse_reads(incorrect_reads_file, False, max_reads=max_incorrect_reads) + \
                     Read.parse_reads(correct_reads_file, max_reads=max_correct_reads)

    def bin_reads(self, bin_fn):
        bins = {}
        for read in self.reads:
            i = bin_fn(read)
            if i not in bins:
                bins[i] = []
            bins[i].append(read)
        return bins

    def get_roc(self, map_q_fn=lambda read : round(read.map_q)):
        """
        Compute a pseudo ROC curve for the reads
        :param map_q_fn: function which determines which value is used as the map-q
        :return:
        """
        bins = self.bin_reads(bin_fn=map_q_fn)
        reads_seen, reads_correct = 0, 0
        roc_bins = []
        for map_q in sorted(bins.keys(), reverse=True):
            for read in bins[map_q]:
                reads_seen += 1
                reads_correct += 1 if read.correct else 0
            roc_bins.append((math.log10((reads_seen - reads_correct)/reads_seen) if reads_seen > reads_correct else -8, reads_correct/len(self.reads), len(bins[map_q]), map_q))
        return roc_bins

    @staticmethod
    def plot_rocs(rocs):
        colors = [ "r--", "bs", 'g^' ] # This needs to be extended if you want more than three lines on the plot
        for roc, color in zip(rocs, colors):
            plt.plot(list(map(lambda x : x[0], roc)), list(map(lambda x : x[1], roc)), color)
        plt.show()

def main():
    start_time = time.time()

    # Parse the reads
    reads = Reads("minimizers_correct_1kg", "minimizers_incorrect_1kg", max_correct_reads=10000)
    print("Got ", len(reads.reads), " in ", time.time()-start_time, " seconds")

    # Print some of the funky reads
    for i, read in enumerate(reads.reads):
        if not read.correct and read.map_q >= 60 and read.fast_cap() >= 60:
            print("Read {}".format(i), read)

    # Make ROC curves
    roc_unmodified = reads.get_roc()
    print("Roc unmodified", roc_unmodified)
    roc_adam_modified = reads.get_roc(map_q_fn=lambda read : round(min(read.adam_cap, read.map_q, 150)))
    print("Roc adam modified ", roc_adam_modified)
    roc_new_sum_modified = reads.get_roc(map_q_fn=lambda read : round(min(read.fast_cap()+30, read.map_q, 150)))
    print("Roc mode modified ", roc_new_sum_modified)
    Reads.plot_rocs([ roc_unmodified, roc_adam_modified, roc_new_sum_modified ])

if __name__ == '__main__':
    main()