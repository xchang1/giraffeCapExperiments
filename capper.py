import math
import matplotlib.pyplot as plt
import numpy as np
import time

total_window_length = 29+11-1
# negative infinity
n_inf = float("-inf")

def log_add(i, j):
    l = lambda x: math.log10(math.pow(10, x) + 1)
    if i < j:
        return l(j - i) + i if j - i <= 7.5 else j
    return l(i - j) + j if i - j <= 7.5 else i

class Minimizer:
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
            format(self.window_start, self.window_length-total_window_length+1, self.minimizer_start, len(self.minimizer), self.minimizer, self.hits, self.hits_in_extension)

class Read:
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
        #- bug here
        #print("Check read: \n", self)
        p_window = 0
        for minimizer in self.minimizers:
            #print(minimizer)
            assert p_window == minimizer.window_start
            p_window += minimizer.window_length - total_window_length + 1
        assert p_window + total_window_length == len(self.read)+1

    def __str__(self):
        return "Read map_q:{} adam_cap:{} xian_cap:{} \n".format(self.map_q, self.adam_cap, self.xian_cap) + "\n\t".join(map(str, self.minimizers)) + "\n"

    def minimizer_interval_iterator(self, start_fn, length):
        """
        Iterator over common intervals of overlapping minimizers.
        :param start_fn: returns the start of the minimizer in the minimizer
        :param length: the length of the minimizer
        :return: yields sequence of (left, right, bottom, top) tuples, where left if the first base
        of the minimizer interval (inclusive), right is the last base of the minimizer interval (exclusive),
        bottom is the index of the first minimizer in the interval and top is the index of the last minimizer
        in the interval (exclusive).
        """
        # Handle no minimizer case
        if len(self.minimizers) == 0:
            return

        stack = [ self.minimizers[0] ]  # Minimizers currently being iterated over
        left = [start_fn(self.minimizers[0])]  # The left end of a minimizer interval
        bottom = [ 0 ]  # The index of the first minimizer in the interval in the sequence of minimizers

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
        for minimizer in self.minimizers[1:]:
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
        p = -self.qual_string[left] / 10
        for i in range(left+1, right):
            p = sum_fn(-self.qual_string[i] / 10, p)
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

    def fast_cap(self, start_fn=lambda x: x.minimizer_start, minimizer_length=29, sum_fn=max):
        """ Compute cap on the reads qual score by considering possibility that base errors and unlocated
        minimizer hits prevented us finding the true alignment

        :param start_fn: Function that returns start of the minimizer
        :param minimizer_length: Length of the minimizer
        :param sum_fn:
        :return: Phred scaled log prob of read being wrongly aligned
        """
        c = np.full(len(self.minimizers) + 1, n_inf)  # Log10 prob of having mutated minimizers,
        # such that c[i+1] is log prob of mutating minimizers 0, 1, 2, ..., i
        c[0] = 0.0

        pBottom = 0
        for left, right, bottom, top in self.minimizer_interval_iterator(start_fn=start_fn, length=minimizer_length):
            # If a new minimizer in the interval include probability of skipping the minimizer
            if pBottom == bottom:
                p = c[bottom] + self.get_log_prob_of_minimizer_skip(self.minimizers[bottom])
                if c[bottom+1] < p:
                    c[bottom+1] = p
                pBottom += 1

            # Calculate prob of all intervals up to top being disrupted
            p = c[bottom] + self.get_log_prob_of_base_error_in_interval(left, right, sum_fn=sum_fn)

            # Replace min-prob for intervals
            for i in range(bottom+1, top+1):
                if c[i] < p:
                    c[i] = p

        assert pBottom == len(self.minimizers)
        assert c[-1] != n_inf

        #print("He he", c)
        return -c[-1] * 10

    def new_cap(self, sum_fn=log_add):
        # Layout states
        m = np.full(len(self.read) + 2, n_inf) # "Mismatch" bases plus start and end state
        c = np.full(len(self.read) + 2, n_inf) # "Conserved" windows states
        m[0] = 0.0 # Start state

        # the log10 prob of conserving a window starting at position i, indexing read from 1
        alpha = np.full(len(self.read)+2, n_inf)
        for minimizer in self.minimizers:
            p = n_inf if (minimizer.hits_in_extension == minimizer.hits and minimizer.hits_in_extension > 0) else 0.0 #math.log10(1.0 - minimizer.hits_in_clusters/minimizer.hits) if minimizer.hits_in_clusters < minimizer.hits else n_inf
            #p = math.log10(1.0 - 1/minimizer.hits) if minimizer.hits_in_extension and minimizer.hits > 1 else n_inf
            for k in range(minimizer.window_length-total_window_length + 1):
                alpha[minimizer.window_start+k] = p

        # the log10 prob of base error at position i, indexing read from 1
        beta = np.full(len(self.read)+2, 0.0)
        for i in range(len(self.read)):
            beta[i+1] = -self.qual_string[i]/10

        n_alpha = np.zeros(len(self.read)+2)  # log10 prob of NOT conserving a window
        n_beta = np.zeros(len(self.read)+2)  # log10 prob of no base error at position i
        for i in range(len(self.read)+2):
            n_alpha[i] = math.log10(1.0 - math.pow(10, alpha[i])) if alpha[i] != 0.0 else n_inf
            n_beta[i] = math.log10(1.0 - math.pow(10, beta[i])) if beta[i] != 0.0 else n_inf

        # For each column in read, indexing positions from 1 and treating end state as first
        # position beyond the end of the read
        for i in range(1, len(self.read)+2):

            # Calculate prob of conserved state
            c[i] = alpha[i] + sum_fn(c[i-1], m[i-1])

            # Calculate prob of mismatch state
            t = beta[i]  # sum of beta and n_beta values for successively longer transitions
            # Prob of transition from previous conserved states
            if i - total_window_length >= 0:
                m[i] = c[i-total_window_length] + n_alpha[i-total_window_length+1] + t
            # Now probs of transitions from previous mismatch states
            for j in range(i-1, max(0, i-total_window_length)-1, -1):
                m[i] = sum_fn(t + n_alpha[j+1] + m[j], m[i])
                t += n_beta[j]

        #print(" Alphas ", alpha)
        #print(" Beta ", beta)
        #print(" Ms ", m)
        #print(" Cs ", c)

        return -10 * m[-1] # Return the end state prob log prob in phred scale

    @staticmethod
    def parse_reads(reads_file, correct=True, maxReads=-1):
        reads = []
        with open(reads_file) as fh: # This masks a bug
            for line in fh:
                try:
                    reads.append(Read(line, correct))
                    #print("Hello ", reads[-1])
                    if reads[-1].map_q >= 60 and not correct and reads[-1].fast_cap() >= 60:
                        print("Hello ", reads[-1], reads[-1].qual_string)
                        print("Read ", "new cap:", reads[-1].new_cap(), " fast cap: ", reads[-1].fast_cap(), " original cap:", reads[-1].map_q, "adam cap:", reads[-1].adam_cap, "xian cap: ", reads[-1].xian_cap, "Correct: ", correct)
                except:
                    pass
                if maxReads != -1 and len(reads) > maxReads:
                    break
        #print("Got {} reads, correct: {}", len(reads), correct)
        return reads

class Reads:
    def __init__(self, correct_reads_file, incorrect_reads_file):
        self.reads = Read.parse_reads(incorrect_reads_file, False) + Read.parse_reads(correct_reads_file) #, maxReads=10000)

    def bin_reads(self, bin_fn):
        bins = {}
        for read in self.reads:
            i = bin_fn(read)
            if i not in bins:
                bins[i] = []
            bins[i].append(read)
        return bins

    def get_roc(self, map_q_fn=lambda read : round(read.map_q)):
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
        colors = [ "r--", "bs", 'g^' ]
        for roc, color in zip(rocs, colors):
            plt.plot(list(map(lambda x : x[0], roc)), list(map(lambda x : x[1], roc)), color)
        plt.show()

def main():
    start_time = time.time()
    reads = Reads("minimizers_correct_1kg", "minimizers_incorrect_1kg")
    print("Got ", len(reads.reads), " in ", time.time()-start_time, " seconds")
    roc_unmodified = reads.get_roc()
    print("Roc unmodified", roc_unmodified)
    roc_adam_modified = reads.get_roc(map_q_fn=lambda read : round(min(read.adam_cap, read.map_q)))
    print("Roc adam modified ", roc_adam_modified)
    roc_new_sum_modified = reads.get_roc(map_q_fn=lambda read : round(min(read.fast_cap(), read.map_q)))
    print("Roc mode modified ", roc_new_sum_modified)
    Reads.plot_rocs([ roc_unmodified, roc_adam_modified, roc_new_sum_modified ])

if __name__ == '__main__':
    main()