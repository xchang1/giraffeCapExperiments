import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import heapq
import logging
import sys
import time

logger = logging.getLogger(__name__)

RC_TABLE = str.maketrans("ACGT", "TGCA")


def reverse_complement(dna):
    return ''.join(reversed(dna.translate(RC_TABLE)))


# negative infinity
n_inf = float("-inf")


def log_add(i, j):
    def log_add_one(x):
        return math.log10(math.pow(10, x) + 1)

    if i < j:
        return log_add_one(j - i) + i if j - i <= 10 else j
    return log_add_one(i - j) + j if i - j <= 10 else i


def prob_to_log10prob(prob):
    """
    Convert a raw probability to a log10 probability.
    """

    return math.log10(prob)


def log10prob_to_phred(log10prob):
    """
    Convert a log10 probability to a Phred score.
    
    :return: the Phred score as a float
    """

    return -10 * log10prob


def log10prob_for_at_least_one(prob_each, count):
    """
    Return the log10 probability for at least one event happening, out of count
    chances each with independent raw probability log10_prob_each.
    """

    return log10prob_not(prob_to_log10prob(1.0 - prob_each) * count)


def log10prob_not(log10prob):
    """
    Return the log10 probability for the opposite of an event with the given log10 probability.
    """

    return math.log10(1.0 - math.pow(10, log10prob))


class Minimizer:
    """
    Represents a minimizer within a read.
    """

    total_window_length = 29 + 11 - 1  # Currently hard-coded, and therefore assumed the same for all minimizers

    def __init__(self, minimizer, minimizer_start, window_start, window_length, hits, hits_in_extension, source_read):
        self.minimizer, self.minimizer_start, self.window_start = minimizer, minimizer_start, window_start
        self.window_length, self.hits, self.hits_in_extension = window_length, hits, hits_in_extension

        if self.minimizer != source_read[minimizer_start:minimizer_start + len(self.minimizer)]:
            # It's a reverse strand minimizer
            self.is_reverse = True
            assert reverse_complement(self.minimizer) == source_read[
                                                         minimizer_start:minimizer_start + len(self.minimizer)]
        else:
            self.is_reverse = False

        # Work out how easy this minimizer is to beat.
        # At hash = 2^64 - 1, you are almost certain to beat it, unless you actually hit it,
        # in which case you won't beat it.
        # At hash = 0 you cannot beat it.
        # At hash = 1 only 1 of the 2^64 possibilities beats it.
        self.beat_prob = self.hash() / (2 ** 64)  # A float on 0-1 that represents the probability of beating this
        # minimizer's hash with a random other minimizer.

        # Consistency checks
        assert len(source_read) >= 0
        assert minimizer_start >= 0 and minimizer_start + len(minimizer) <= len(source_read)
        assert window_start >= 0 and window_start + window_length <= len(source_read)
        assert window_start <= minimizer_start
        assert window_start + window_length >= minimizer_start + len(minimizer)

    def __str__(self):
        return "window_start: {} window_length: {} minimizer_start: {} " \
               "minimizer_length: {} minimizer: {} hits: {} hits_in_extension: {} hash: {:016x} beat_prob: {}". \
            format(self.window_start, self.window_length - Minimizer.total_window_length + 1, self.minimizer_start,
                   len(self.minimizer), self.minimizer, self.hits, self.hits_in_extension, self.hash(), self.beat_prob)

    def encode(self):
        """
        Encode the minimizer as an integer.
        
        :return: the minimizer 2-bit encoded as A=0, C=1, G=2, T=3, last base
        in least significant position.
        """

        encoding = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        encoded = 0
        for base in self.minimizer:
            # Check the base for acceptability
            assert base in encoding
            # Encode and shift in each base.
            encoded = (encoded << 2 | encoding[base])
        return encoded

    def hash(self):
        """
        Get the integer hash value of the current minimizer.
        Hashes are over the 2^64 space.
        
        :return: the minimizer's hash value as a nonnegative int in 64 bits
        """

        # Numpy 64 bit ints warn on overflow and ctypes 64 bit ints don't implement inportant ops like xor.
        # So we use Python bigints and mask frequently.
        mask = 0xFFFFFFFFFFFFFFFF

        scratch = self.encode() & mask

        # Implement Thomas Wang's 64 bit integer hash function as in gbwt library
        scratch = ((scratch ^ mask) + ((scratch << 21) & mask)) & mask
        scratch = scratch ^ (scratch >> 24)
        scratch = (((scratch + (scratch << 3) & mask) & mask) + ((scratch << 8) & mask)) & mask
        scratch = (scratch ^ (scratch >> 14)) & mask
        scratch = (((scratch + ((scratch << 2) & mask)) & mask) + ((scratch << 4) & mask)) & mask
        scratch = (scratch ^ (scratch >> 28)) & mask
        scratch = (scratch + ((scratch << 31) & mask)) & mask
        return scratch


class UnacceptableReadError(RuntimeError):
    pass


class Read:
    """
    Represents a read and its minimizers.
    """

    def __init__(self, line, correct):

        # Format is tab sep (on one line):
        # READ_STR QUAL_STR IGNORED
        # (MINIMIZER_STR START WINDOW_START WINDOW_LENGTH HITS HITS_IN_EXTENSIONS?)xN
        # UNCAPPED_MAP_Q CAP CAP STAGE

        tokens = line.split()
        self.read, self.qual_string = tokens[0], [ord(i) - ord("!") for i in tokens[1]]
        assert len(self.read) == len(self.qual_string)
        # raw mapq, adam's mapq_extend_cap, my probability cluster lost cap,  last correct stage
        # Note that last correct stage may not be there
        if len(tokens) % 6 == 1:
            # Stage should be there
            self.uncapped_map_q, self.vg_computed_cap, self.xian_cap, self.stage = \
                float(tokens[-4]), float(tokens[-3]), float(tokens[-2]), tokens[-1]
        else:
            # No stage
            self.uncapped_map_q, self.vg_computed_cap, self.xian_cap = \
                float(tokens[-3]), float(tokens[-2]), float(tokens[-1])
             self.stage = "unknown"
        self.correct = correct
        self.minimizers = sorted([Minimizer(tokens[i + 3], int(tokens[i + 4]), int(tokens[i + 5]), int(tokens[i + 6]),
                                            int(tokens[i + 7]), int(tokens[i + 8]), self.read) for i in
                                  range(0, len(tokens) - 7, 6)],
                                 key=lambda x: x.window_start)  # Sort by start coordinate

        self.minimizers = sorted(self.minimizers, key=lambda x: x.minimizer_start)  # Fix broken minimizer ordering

        # Check they don't overlap and every window is accounted for
        p_window_start = -1
        p_window_end = -1
        for minimizer in self.minimizers:
            if p_window_start > minimizer.window_start:
                # Read contains minimizers that are out of order or overlap on start coordinates
                raise UnacceptableReadError(str(self))
            p_window_start = minimizer.window_start

            if p_window_end > minimizer.window_start + minimizer.window_length:
                # Read contains minimizers that are out of order or overlap on end coordinates
                raise UnacceptableReadError(str(self))

            p_window_end = minimizer.window_start + minimizer.window_length

            if minimizer.window_start + minimizer.window_length > len(self.read):
                # Read contains minimizer that extends beyond the read length
                raise UnacceptableReadError(str(self))

    def __str__(self):
        return "Read uncapped_map_q:{} vg_computed_cap:{} xian_cap:{} faster_cap:{} unique_cap:{} balanced_cap:{} stage: {}" \
            "\n\tread_string: {}\n\tqual_string: {} \n {}\n".format(self.uncapped_map_q, self.vg_computed_cap, self.xian_cap,
                                                                    self.faster_cap(), self.faster_unique_cap(),
                                                                    self.faster_balanced_cap(), self.stage,
                                                                    self.read, " ".join(map(str, self.qual_string)),
                                                                    "\n\t".join(map(str, self.minimizers)))

    def visualize(self, out=sys.stdout, minimizers=None):
        """
        Print out the read sequence with the minimizers aligned below it.
        """

        # Work out the set of minimizers we want
        if minimizers is None:
            minimizers = self.minimizers

        # First all the position numbers

        # Work out how namy digits we need for position
        digits = int(math.ceil(math.log10(len(self.read))))
        # Pad each number to that width
        texts = [('{:>' + str(digits) + '}').format(i) for i in range(len(self.read))]

        for row in range(digits):
            # For each digit we need
            for text in texts:
                # Print that digit of each padded number
                out.write(text[row])
            out.write('\n')

        # Then the sequence
        out.write(self.read)
        out.write('\n')

        # Then the minimizers
        for minimizer in minimizers:
            # Work out what orientation to print
            minimizer_text = minimizer.minimizer if not minimizer.is_reverse else reverse_complement(
                minimizer.minimizer)
            for i in range(len(self.read)):
                if i < minimizer.window_start:
                    # Before agglomeration
                    out.write(' ')
                elif i < minimizer.minimizer_start:
                    # In agglomeration before minimizer
                    out.write('-')
                elif i < minimizer.minimizer_start + len(minimizer.minimizer):
                    # In minimizer
                    out.write(minimizer_text[i - minimizer.minimizer_start])
                elif i < minimizer.window_start + minimizer.window_length:
                    # In agglomeration after minimizer (before required agglomeration length has elapsed)
                    out.write('-')
                else:
                    # After minimizer
                    out.write(' ')
            out.write('\n')

    def recompute_vg_computed_cap(self):
        """
        Recompute the "Adam Cap": probability of the most probable way to have
        errors at a set of bases, such that each error disrupts all minimizers
        whose agglomerations it hits, in order to disrupt all minimizers that
        were explored. Explored is defined as having participated in any
        extended cluster.
        
        :return: Phred scaled log prob of read being wrongly aligned
        """

        # Filter down to just located minimizers
        minimizers = [m for m in self.minimizers if m.hits_in_extension > 0]
        next_minimizer = 0

        if len(minimizers) == 0:
            # No minimizers left = can't cap
            # Should probably actually treat as 0.
            return float('inf')

        # Make a DP table of Phred score for cheapest solution up to here, with an error at this base
        table = np.full(len(self.read), n_inf)

        # Have a priority queue of agglomerations, with earliest-ending and then latest-starting being at the front
        # Agglomerations are (end, -start, minimizer)
        active_agglomerations = []

        # Have another priority queue of windows, with earliest-ending and then latest-starting being at the front
        # Windows are (end, -start)
        # These aren't all overlapped, but none have ended yet.
        active_windows = []

        # Track the latest-starting window ending before here, as (start, end), if any
        prev_window = None

        for i in range(len(self.read)):
            # Go through all the bases in the read.

            logger.debug("At base %d", i)

            while next_minimizer < len(minimizers) and minimizers[next_minimizer].window_start == i:

                # While the next agglomeration starts here
                agg = minimizers[next_minimizer]
                # Put it in the queue
                heapq.heappush(active_agglomerations, (agg.window_start + agg.window_length, -agg.window_start, agg))

                logger.debug("Agglomeration starts here: %s", agg)

                for window_offset in range(agg.window_length - agg.total_window_length + 1):
                    # And also all its windows
                    # Represent them as end, start for sort order
                    heapq.heappush(active_windows, (
                        agg.window_start + window_offset + agg.total_window_length,
                        -(agg.window_start + window_offset)))

                    logger.debug("Creates window %d to %d", agg.window_start + window_offset,
                                 agg.window_start + window_offset + agg.total_window_length)

                # Move on to the next minimizer
                next_minimizer += 1

            if len(active_agglomerations) == 0:
                # Nothing to do here
                logger.debug("Ignore base: no agglomerations")
            else:

                # Phred cost of an error here
                base_cost = self.qual_string[i]

                logger.debug("Cost to break base: %f", base_cost)
                for _, _, active_agg in active_agglomerations:
                    # For each minimizer this base is in an agglomeration of (i.e. that
                    # is active)

                    logger.debug("Overlap agglomeration: %s", active_agg)

                    # Work out the Phred cost of breaking the minimizer for
                    # the windows of it that we overlap by having an error here

                    if active_agg.minimizer_start <= i < active_agg.minimizer_start + len(active_agg.minimizer):
                        # We are in the minimizer itself.
                        # No cost to break
                        minimizer_break_score = 0
                        logger.debug("In minimizer: no additional cost")
                    else:
                        # We are in the flank.

                        # How many new possible minimizers would an error here create in this agglomeration,
                        # to compete with its minimizer?
                        # No more than one per position in a minimizer sequence.
                        # No more than 1 per base from the start of the agglomeration to here, inclusive.
                        # No more than 1 per base from here to the last base of the agglomeration, inclusive.
                        possible_minimizers = min(len(active_agg.minimizer),
                                                  min(i - active_agg.window_start + 1,
                                                      (active_agg.window_start + active_agg.window_length) - i))

                        # How probable is it to beat the current minimizer?
                        beat_prob = active_agg.beat_prob

                        # We have that many chances to beat this minimizer. What's that in Phred?
                        minimizer_break_score = log10prob_to_phred(
                            log10prob_for_at_least_one(beat_prob, possible_minimizers))

                        logger.debug("Out of minimizer: need to beat with probability %f and %d chances for %f points",
                                     beat_prob, possible_minimizers, minimizer_break_score)

                    # And AND them all up
                    base_cost += minimizer_break_score

                if prev_window is not None:
                    # AND that with the cheapest thing in the latest-starting window ending before here
                    prev_cost = np.min(table[prev_window[0]:prev_window[1]])
                    logger.debug("Cheapest previous base in %d to %d is %f", prev_window[0], prev_window[1], prev_cost)
                    base_cost += prev_cost

                logger.debug("Total cost for best solution ending in this base: %f", base_cost)

                # Save into the DP table
                table[i] = base_cost

            # Kick out agglomerations and windows we have passed

            while len(active_agglomerations) > 0 and active_agglomerations[0][0] == i + 1:
                # We are the last base in this agglomeration.
                # Pop it. It won't be active anymore.
                logger.debug("End agglomeration: %s", active_agglomerations[0][2])
                heapq.heappop(active_agglomerations)

            while len(active_windows) > 0 and active_windows[0][0] == i + 1:
                # We are the last base in this window

                logger.debug("End window %s to %s", -active_windows[0][1], active_windows[0][0])

                if prev_window is None or prev_window[0] < -active_windows[0][1]:
                    # We have no latest-starting window ending before here, or this one starts later.
                    # Take it (making sure to decode its weird encoding)
                    prev_window = (-active_windows[0][1], active_windows[0][0])
                    logger.debug("It is new latest-starting window ending before here")
                # Pop it
                heapq.heappop(active_windows)

        # Now the DP table is filled

        # Scan the latest-ending, latest-starting window to find the total cost.
        assert prev_window is not None
        winner = np.min(table[prev_window[0]:prev_window[1]])
        logger.debug("Cheapest answer in final window %d to %d is %f", prev_window[0], prev_window[1], winner)
        return winner

    def agglomeration_interval_iterator(self, minimizers, start_fn=lambda x: x.window_start,
                                        length_fn=lambda x: x.window_length):
        """
        Iterator over common intervals of overlapping minimizer agglomerations.
        
        :param minimizers: the list of minimizers to iterate over
        :param start_fn: given a minimizer, returns the start of its agglomeration.
        :param length_fn: given a minimizer, returns the length of its agglomeration.
        
        :return: yields sequence of (left, right, bottom, top) tuples, where
        left is the first base of the agglomeration interval (inclusive), right
        is the last base of the agglomeration interval (exclusive), bottom is
        the index of the first minimizer with an agglomeration in the interval
        and top is the index of the last minimizer with an agglomeration in the
        interval (exclusive).
        """

        return self.overlap_interval_iterator(minimizers, start_fn, length_fn)

    def overlap_interval_iterator(self, items, start_fn, length_fn):
        """
        Iterator over common intervals (rectangles) of overlapping items.
        Common intervals break when an item starts or ends.
        Assumes items are arranged such that they do not contain each other.
        Assumes items are sorted.
        
        :param items: the list of items to iterate over.
        :param start_fn: given an item, returns the start of an item.
        :param length_fn: given an item, returns the length of the item.
        
        :return: yields sequence of (left, right, bottom, top) tuples, where
        left is the first base of the interval (inclusive), right is the last
        base of the interval (exclusive), bottom is the index of the first item
        in the interval and top is the index of the last item in the interval
        (exclusive).
        """
        # Handle no item case
        if len(items) == 0:
            return

        stack = [items[0]]  # Items currently being iterated over
        left = [start_fn(items[0])]  # The left end of an item interval
        bottom = [0]  # The index of the first item in the interval in the sequence of items

        def get_preceding_intervals(right):
            # Get all intervals that precede a given point "right"
            while left[0] < right:

                # Case where the left-most item ends before the start of the new item
                if start_fn(stack[0]) + length_fn(stack[0]) <= right:
                    yield left[0], start_fn(stack[0]) + length_fn(stack[0]), bottom[0], bottom[0] + len(stack)

                    # If the stack contains only one item there is a gap between the item
                    # and the new item, otherwise just shift to the end of the leftmost item
                    left[0] = right if len(stack) == 1 else start_fn(stack[0]) + length_fn(stack[0])

                    bottom[0] += 1
                    del (stack[0])

                # Case where the left-most item ends at or after the beginning of the new new item
                else:
                    yield left[0], right, bottom[0], bottom[0] + len(stack)
                    left[0] = right

        # For each item in turn
        for item in items[1:]:
            assert len(stack) > 0

            # For each new item we return all intervals that
            # precede start_fn(item)
            for interval in get_preceding_intervals(start_fn(item)):
                yield interval

            stack.append(item)  # Add the new item for the next loop

        # Intervals of the remaining intervals on the stack
        for interval in get_preceding_intervals(len(self.read)):
            yield interval

    def get_log_prob_of_disruption_in_column(self, minimizers, index):
        """ Gives the log10 prob of a base error in the given column of the read that disrupts all the
        overlapping minimizers.
        :param minimizers: Minimizers to be disrupted
        :param index: Index of position in read
        :return: log10 prob of disruption as float
        """
        p = -(self.qual_string[index]) / 10
        for m in minimizers:
            if not (m.minimizer_start <= index < m.minimizer_start + len(m.minimizer)):
                # Index is out of range of the minimizer

                # How many new possible minimizers would an error here create in this agglomeration,
                # to compete with its minimizer?
                # No more than one per position in a minimizer sequence.
                # No more than 1 per base from the start of the agglomeration to here, inclusive.
                # No more than 1 per base from here to the last base of the agglomeration, inclusive.
                possible_minimizers = min(len(m.minimizer),
                                          min(index - m.window_start + 1,
                                              (m.window_start + m.window_length) - index))

                # We have that many chances to beat this minimizer.
                p += log10prob_for_at_least_one(m.beat_prob, possible_minimizers)
        logger.debug("Overall log10 prob for disrupting all minimizers via column %d: %f", index, p)
        return p

    def get_log_prob_of_disruption_in_interval(self, minimizers, left, right, sum_fn=log_add):
        """
        Gives the log10 prob of a base error in the given interval of the read,
        accounting for the disruption of minimizers
        :param minimizers: Minimizers to be disrupted
        :param left: Leftmost base in interval of read, inclusive
        :param right: Rightmost base in interval of read, exclusive
        :param sum_fn: Function for adding/maxing log probs in interval
        :return: log10 prob as float
        """
        if left == right:
            # 0-length intervals do not need to be disrupted.
            return 0
        p = self.get_log_prob_of_disruption_in_column(minimizers, left)
        for i in range(left + 1, right):
            p = sum_fn(self.get_log_prob_of_disruption_in_column(minimizers, i), p)
            
        # With really bad base qualities our assumption that OR ~= total of
        # probabilities which we might make will break down because the errors
        # aren't actually exclusive.
        
        # Cap at 0 (certainty)
        p = min(p, 0.0)
        
        logger.debug("Overall log10 prob for disrupting %d to %d: %f", left, right, p)
        return p

    def faster_balanced_cap(self, sum_fn=log_add):
        """ As faster_cap(), but scores minimizers
        :param sum_fn: Function for adding/maxing log probs in interval
        :return: Phred scaled log prob of read being wrongly aligned
        """
        x = self.faster_unique_cap(sum_fn=sum_fn)
        y = self.faster_cap(sum_fn=sum_fn)
        return min(x + 30, y)

    def faster_unique_cap(self, sum_fn=log_add):
        """ As faster_cap(), but only consider unique minimizers in the read.
        :param sum_fn: Function for adding/maxing log probs in interval
        :return: Phred scaled log prob of read being wrongly aligned
        """
        return self.faster_cap(sum_fn=sum_fn, minimizer_filter_fn=lambda x: 0 < x.hits_in_extension == x.hits)

    def faster_cap(self, sum_fn=log_add, minimizer_filter_fn=lambda x: x.hits_in_extension > 0):
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

        We use dynamic programming sweeping left-to-right over the intervals to compute the probability of
        the minimum number of base errors needed to disrupt all the minimizers.

        :param sum_fn: method for adding up probability of one or more base errors in an interval of bases
        :param minimizer_filter_fn: Function to filter which minimizers are considered
        :return: Phred scaled log prob of read being wrongly aligned
        """
        # This step filters minimizers we will definitely skip - needs to be coordinated with
        # get_log_prob_of_minimizer_skip function above
        filtered_minimizers = list(filter(minimizer_filter_fn, self.minimizers))
        
        self.visualize(minimizers=filtered_minimizers)

        c = np.full(len(filtered_minimizers) + 1, n_inf)  # Log10 prob of having mutated minimizers,
        # such that c[i+1] is log prob of mutating minimizers 0, 1, 2, ..., i
        c[0] = 0.0

        for left, right, bottom, top in self.agglomeration_interval_iterator(filtered_minimizers):
            # Calculate prob of all intervals up to top being disrupted
            p = c[bottom] + self.get_log_prob_of_disruption_in_interval(filtered_minimizers[bottom:top], left, right,
                                                                        sum_fn=sum_fn)

            # Replace min-prob for minimizers in the interval
            for i in range(bottom + 1, top + 1):
                if c[i] < p:
                    c[i] = p

        assert c[-1] != n_inf
        cap = -c[-1] * 10
        logger.debug("Overall cap: %s", cap)
        return cap

    @staticmethod
    def parse_reads(reads_file, correct=True, max_reads=-1):
        reads = []
        with open(reads_file) as fh:  # This masks a bug
            for line in fh:
                reads.append(Read(line, correct))
                print(reads[-1])
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

    def get_roc(self, map_q_fn=lambda read: round(read.uncapped_map_q)):
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
            roc_bins.append((
                math.log10((reads_seen - reads_correct) / reads_seen) if reads_seen > reads_correct else -8,
                reads_correct / len(self.reads), len(bins[map_q]), map_q))
        return roc_bins

    @staticmethod
    def plot_rocs(rocs):
        colors = ["r--", "bs", 'g^']  # This needs to be extended if you want more than three lines on the plot
        for roc, color in zip(rocs, colors):
            plt.plot(list(map(lambda x: x[0], roc)), list(map(lambda x: x[1], roc)), color)
        plt.savefig('roc.svg')
        #plt.show()


def main():
    logging.basicConfig(level=logging.DEBUG)
    start_time = time.time()

    # Parse the reads
    reads = Reads("minimizers_correct_1kg", "minimizers_incorrect_1kg", max_correct_reads=10000)
    print("Got ", len(reads.reads), " in ", time.time() - start_time, " seconds")

    # Print some of the funky reads
    for i, read in enumerate(reads.reads):
        if not read.correct and round(0.85 * min(2.0 * read.faster_cap(), read.uncapped_map_q)) >= 60:
            print("Read {} {}".format(i, read))

    # Make ROC curves
    roc_unmodified = reads.get_roc()
    print("Roc unmodified", roc_unmodified)
    roc_adam_modified = reads.get_roc(map_q_fn=lambda r: round(min(r.vg_computed_cap, r.uncapped_map_q, 60)))
    print("Roc adam modified ", roc_adam_modified)
    roc_new_sum_modified = reads.get_roc(map_q_fn=lambda r: round(0.85 * min(2.0 * r.faster_cap(), r.uncapped_map_q, 70)))
    print("Roc mode modified ", roc_new_sum_modified)
    Reads.plot_rocs([roc_unmodified, roc_adam_modified, roc_new_sum_modified])

    # plt.scatter([x.vg_computed_cap for x in reads.reads if not x.correct],
    # [x.faster_balanced_cap() for x in reads.reads if not x.correct])
    # plt.scatter([x.vg_computed_cap for x in reads.reads if x.correct],
    # [x.faster_cap() for x in reads.reads if x.correct], 'g^')
    # plt.show()
    
    plt.scatter([x.vg_computed_cap for x in reads.reads], [x.faster_cap() for x in reads.reads])
    plt.savefig('compare.svg')


if __name__ == '__main__':
    main()
