import math
import matplotlib.pyplot as plt
import time

total_window_length = 29+11-1

class Minimizer:
    def __init__(self, minimizer, minimizer_start, window_start, window_length, hits, hits_in_clusters, read_length):
        self.minimizer, self.minimizer_start, self.window_start, self.window_length, self.hits, self.hits_in_clusters = \
            minimizer, minimizer_start, window_start, window_length, hits, hits_in_clusters
        # Sanity checks
        assert read_length >= 0
        assert minimizer_start >= 0 and minimizer_start + len(minimizer) <= read_length
        assert window_start >= 0 and window_start + window_length <= read_length
        assert window_start <= minimizer_start
        assert window_start + window_length >= minimizer_start + len(minimizer)
        assert hits_in_clusters >= 0
        assert hits_in_clusters <= hits

    def __str__(self):
        return "window_start: {} window_length: {} minimizer_start: {} minimizer_length: {} minimizer: {} hits: {} hits_in_clusters: {}".\
            format(self.window_start, self.window_length-total_window_length+1, self.minimizer_start, len(self.minimizer), self.minimizer, self.hits, self.hits_in_clusters)

class Read:
    def __init__(self, line, correct):
        # Format is tab sep:
        #READ_STRING QUAL_STRING (MINIMIZER_STRING START WINDOW_START WINDOW_LENGTH HITS HITS_IN_CLUSTERS)xN MAP_Q STAGE
        tokens = line.split()
        self.read, self.qual_string = tokens[0], [ ord("?") - ord("!") for i in tokens[1] ]
        assert len(self.read) == len(self.qual_string)
        self.map_q, self.stage = float(tokens[-2]), tokens[-1]
        self.correct = correct
        self.minimizers = sorted([ Minimizer(tokens[i+3], int(tokens[i+4]), int(tokens[i+5]), int(tokens[i+6]), int(tokens[i+7]), int(tokens[i+8]), len(self.read)) for i in range(0, len(tokens)-5, 6) ], key=lambda x : x.window_start) # Sort by start coordinate
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
        return "Read\n" + "\n\t".join(map(str, self.minimizers)) + "\n"

    def new_cap(self):
        #for each minimizer
        pass


    @staticmethod
    def parse_reads(reads_file, correct):
        reads = []
        with open(reads_file) as fh: # This masks a bug
            for line in fh:
                try:
                    reads.append(Read(line, correct))
                except:
                    pass
        #print("Got {} reads, correct: {}", len(reads), correct)
        return reads


class Reads:
    def __init__(self, correct_reads_file, incorrect_reads_file):
        self.reads = Read.parse_reads(correct_reads_file, True) + Read.parse_reads(incorrect_reads_file, False)

    def bin_reads(self, bin_fn):
        bins = {}
        for read in self.reads:
            i = bin_fn(read)
            if i not in bins:
                bins[i] = []
            bins[i].append(read)
        return bins

    def get_roc(self, map_q_fn=lambda read : read.map_q):
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
    def plot_rocs(roc):
        plt.plot(list(map(lambda x : x[0], roc)), list(map(lambda x : x[1], roc)), 'r--')
        plt.show()

start_time = time.time()
reads = Reads("minimizers_correct_1kg_500000", "minimizers_incorrect_1kg")
print("Got ", len(reads.reads), " in ", time.time()-start_time, " seconds")
roc = reads.get_roc()
print("Roc", roc)
Reads.plot_rocs(roc)