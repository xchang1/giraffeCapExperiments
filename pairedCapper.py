from capper import Read, Reads, n_inf, p_inf
import time


class PairedRead(Read):
    """
    Represents one end of a paired read and its minimizers
    """

    def __init__(self, line, correct):
        # Format is tab sep (on one line):
        # READ_STR, QUAL_STR,  # fragment clusters,
        # bool read was found by rescue, bool read was used to rescue its pair,
        # fragment length,
        # (minimizer, offset of minimizer in read, agglomeration start, agglomeration length, # hits,
        # hits in a cluster that got aligned)xN,
        # uncapped mapq, fragment cluster cap, mapq of score group, mapq extend cap, last correct stage

        # Notes:
        # (1) The uncapped mapq is found from the scores of pairs of alignments, so it should be the same for both
        # reads in a pair.

        # (2) Fragment cluster cap is 1 /  # equivalent or better fragment clusters converted to phred and it should
        # also be the same for both reads.

        # (3) The score group mapq is the mapq calculated from the scores of all alignments that also paired
        # with this alignment's pair.
        # This isn't calculated when we found the alignments by rescue.

        # (4) Mapq extend cap is Adam's cap, the same as for single end

        tokens = line.split()
        self._parse_read(tokens, correct)

        self.fragment_clusters, self.rescued, self.rescuer, self.fragment_length = \
            int(tokens[2]), bool(tokens[3]), bool(tokens[4]), int(tokens[5])

        # raw mapq, adam's mapq_extend_cap, my probability cluster lost cap,  last correct stage
        self.map_q, self.xian_cap, self.score_group_map_q, self.adam_cap, self.stage = \
            float(tokens[-5]), float(tokens[-4]), float(tokens[-3]), float(tokens[-2]), tokens[-1]

        self._parse_minimizers(tokens[6:-5])

        self._check_minimizers()

        self.pair = None  # Pair is initially none

    @staticmethod
    def parse_reads(reads_file, correct=True, mixed=False, max_reads=-1):
        reads = []
        pair = None
        with open(reads_file) as fh:  # This masks a bug
            for line in fh:
                if mixed:  # Deal with read pairs where one is correctly
                    # mapped and the other is not
                    tokens = line.split()
                    correct = bool(tokens[-1])
                    line = "\t".join(tokens[:-1])

                reads.append(PairedRead(line, correct))
                if pair is not None:  # Has a pair
                    reads[-1].pair = pair
                    pair.pair = reads[-1]
                    pair = None
                else:
                    pair = reads[-1]
                if max_reads != -1 and pair is None and len(reads) > max_reads:
                    break
        assert pair is None  # Should not have any unpaired reads
        return reads


class PairedReads(Reads):
    def __init__(self, correct_reads_file, incorrect_reads_file, mixed_reads_file,
                 max_correct_reads=-1, max_incorrect_reads=-1, max_mixed_reads=-1):
        self.reads = PairedRead.parse_reads(incorrect_reads_file, False, max_reads=max_incorrect_reads) + \
                     PairedRead.parse_reads(correct_reads_file, max_reads=max_correct_reads) + \
                     PairedRead.parse_reads(mixed_reads_file, mixed=True, max_reads=max_mixed_reads)


def main():
    start_time = time.time()

    # Parse the reads
    reads = PairedReads("minimizers_paired/minimizers_correct_1kg_paired",
                        "minimizers_paired/minimizers_incorrect_1kg_paired",
                        "minimizers_paired/minimizers_mixed_1kg_paired",
                        )  # max_correct_reads=10000)

    print("Got ", len(reads.reads), " in ", time.time() - start_time, " seconds")

    def f(x):  # Remove 0 values from mapq calcs
        return x if x != 0 and x != -0 and x != n_inf else p_inf

    def proposed_cap(r):
        # The proposed map_q function
        return round(min(f(r.xian_cap), f(r.score_group_map_q/2.0),
                     r.faster_cap() + r.pair.faster_cap(), r.map_q, 60))

    def current_cap(r):
        """
        It should be the minimum of the
        (i) uncapped mapq (r.map_q),
        (ii) fragment cluster cap (r.xian_cap),
        (iii) score group mapq divided by 2 (r.score_group_map_q), and
        (iv) mapq extend cap (r.adam_cap)
        If a value is 0, -0, inf, or -inf then it is ignored
        :param r:
        :return:
        """
        return round(min(r.map_q, f(r.xian_cap), f(r.score_group_map_q/2.0), max(r.adam_cap, r.pair.adam_cap), 60))

    # Print some of the funky reads
    for i, read in enumerate(reads.reads):
        if not read.correct and proposed_cap(read) >= 60:
            print("Read {} {}".format(i, read))

    # Make ROC curves
    roc_unmodified = reads.get_roc()
    print("Roc unmodified", roc_unmodified)
    roc_adam_modified = reads.get_roc(map_q_fn=current_cap)
    print("Roc adam modified ", roc_adam_modified)
    start_time = time.time()
    roc_new_sum_modified = reads.get_roc(map_q_fn=proposed_cap)
    print("Roc mode modified (time:{}) ".format(time.time()-start_time), roc_new_sum_modified)
    Reads.plot_rocs([roc_unmodified, roc_adam_modified, roc_new_sum_modified])

    # plt.scatter([x.adam_cap for x in reads.reads if not x.correct],
    # [x.faster_balanced_cap() for x in reads.reads if not x.correct])
    # plt.scatter([x.adam_cap for x in reads.reads if x.correct],
    # [x.faster_cap() for x in reads.reads if x.correct], 'g^')
    # plt.show()


if __name__ == '__main__':
    main()
