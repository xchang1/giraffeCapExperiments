from capper import Read, Reads, n_inf, p_inf, UnacceptableReadError
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time
import math


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
        # read sequence, quality,
        self._parse_read(tokens, correct)

        # fragment clusters, bool read was found by rescue, bool read was used to rescue its pair, fragment length,
        self.fragment_clusters, self.rescued, self.rescuer, self.fragment_length = \
            int(tokens[2]), bool(int(tokens[3])), bool(int(tokens[4])), int(tokens[5])

        unmapped = ";" not in tokens[-2]

        if unmapped:
            #If this was unmapped and has no minimizers, we'd miss the alignment score token so everything would be off by token
            self.uncapped_map_q, self.xian_cap, self.score_group_map_q, self.vg_computed_cap, self.xian_new_cap, self.capped_map_q, self.stage = \
                float(tokens[-7]), float(tokens[-6]), float(tokens[-5]), float(tokens[-4]), float(tokens[-3]), float(tokens[-2]), tokens[-1]
            self.alignment_scores = []
            self.multiplicities = []

        else:
            # uncapped mapq, old xian cap, mapq of score group, adam cap, new xian cap, final mapq, (alignment score1, alignment score 2, fragment length log
            # likelihood, multiplicity, score of alignment pair)xN last correct stage
            self.uncapped_map_q, self.xian_cap, self.score_group_map_q, self.vg_computed_cap, self.xian_new_cap, self.capped_map_q, self.stage = \
                float(tokens[-8]), float(tokens[-7]), float(tokens[-6]), float(tokens[-5]), float(tokens[-4]), float(tokens[-3]), tokens[-1]

            #print("howdy", fragment_length_scores)
            all_aln_scores = [tuple(float(x) for x in token.split(",")) for token in tokens[-2].split(";")[:-1]] 
            self.alignment_scores = []
            self.multiplicities = []
            for (score1, score2, log_likelihood, multiplicity, pair_score) in all_aln_scores:
                self.alignment_scores.append(pair_score)
                self.multiplicities.append(multiplicity)


        #print("Alignment scores", self.alignment_scores)

        #[[float(j) for j in i.split(",") if j != ""] for i in self.alignment_scores.split(";") if i != ""]
        #print("hello", self.fragment_length_scores, self.alignment_pairs, self.alignment_scores)
        #assert(len(alignment_pairs) == len(self.alignment_scores))

        #self.alignment_scores = sorted([float(i.replace(";", "")) for i in self.alignment_scores.split(";") if i.replace(";", "") != ""])

        if self.uncapped_map_q < 0:
            self.uncapped_map_q = 0

        if not unmapped:

            # (minimizer, offset of minimizer in read, agglomeration start, agglomeration length,  # hits, # hits in a cluster that got aligned)xN,
            try:
                self._parse_minimizers(tokens[6:-8])
            except ValueError:
                # Probably a single-end read
                raise UnacceptableReadError("Probably single-ended")
            except IndexError:
                # Probably a single-end read
                raise UnacceptableReadError("Probably single-ended")
        else:
            self.minimizers = []

        self._check_minimizers()

        self.pair = None  # Pair is initially none

        self.faster_cap_precomputed = self.faster_cap()  # Precompute the faster cap

    def __str__(self):
        return "Read uncapped_map_q:{} rescued: {} rescuer: {} fragment_clusters: {} fragment_length: {} \n" \
            "vg_computed_cap:{} xian_cap:{} score_group_map_q: {} " \
            "faster_cap:{} unique_cap:{} balanced_cap:{} stage: {}" \
            "\n\tread_string: {}\n\tqual_string: {} \n {}\n".format(self.uncapped_map_q, self.rescued, self.rescuer,
                                                                    self.fragment_clusters, self.fragment_length,
                                                                    self.vg_computed_cap, self.xian_cap,
                                                                    self.score_group_map_q, self.faster_cap(),
                                                                    self.faster_unique_cap(),
                                                                    self.faster_balanced_cap(), self.stage,
                                                                    self.read, " ".join(map(str, self.qual_string)),
                                                                    "\n\t".join(map(str, self.minimizers)))


    @staticmethod
    def parse_reads(reads_file, correct=True, mixed=False, max_reads=-1):
        reads = []
        pair = None
        with open(reads_file) as fh:
            for line in fh:
                if mixed:  # Deal with read pairs where one is correctly
                    # mapped and the other is not
                    tokens = line.split()
                    correct = bool(int(tokens[-1]))
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
    reads = PairedReads("minimizers_1kg_paired/minimizers_correct",
                        "minimizers_1kg_paired/minimizers_incorrect",
                        "minimizers_1kg_paired/minimizers_mixed",
                        max_correct_reads=10000)

    print("Got ", len(reads.reads), " in ", time.time() - start_time, " seconds")

    def f(x):  # Remove 0 values from mapq calcs
        return x if x != n_inf else p_inf  # the str 0.0 is there because -0.0 is actually 0

    def proposed_cap(r):
        # The proposed map_q function

        # Note the f function is different to f2 in the current cap - I suspect a bug in the way that
        # the way we modify the score_group_map_q and xian_cap exists in Giraffe
        # that treats 0.0s weirdly.

        # r.score_group_map_q/6.0,

        #return int(min(r.xian_new_cap, (r.faster_cap() + r.pair.faster_cap()),
        #               r.uncapped_map_q, 60))

        return int(min(r.xian_cap, (r.faster_cap() + r.pair.faster_cap()),
                       r.uncapped_map_q, 60))

    def current_cap(r):
        """
        It should be the minimum of the
        (i) uncapped mapq (r.uncapped_map_q),
        (ii) fragment cluster cap (r.xian_cap),
        (iii) score group mapq divided by 2 (r.score_group_map_q), and
        (iv) mapq extend cap (r.vg_computed_cap)
        If a value is 0, inf, or -inf then it is ignored
        :param r:
        :return:
        """
        #map_q = int(min(r.uncapped_map_q, r.xian_cap, r.score_group_map_q/2.0, (r.vg_computed_cap+r.pair.vg_computed_cap)/2.0, 60))
        #assert map_q == r.capped_map_q
        return r.capped_map_q #map_q

    # Print some of the funky reads
    wrong = 0
    stages_wrong = {}
    rescued_wrong = 0
    rescuer_wrong = 0
    for i, read in enumerate(reads.reads):
        if not read.correct and proposed_cap(read) >= 30:
            print("Read {} {}".format(i, read))
            wrong += 1
            if read.stage not in stages_wrong:
                stages_wrong[read.stage] = 0
            stages_wrong[read.stage] += 1
            if read.rescued:
                rescued_wrong += 1
            if read.rescuer:
                rescuer_wrong += 1
    print("We got {} wrong, {} rescued-wrong, {} rescuer-wrong".format(wrong, rescued_wrong, rescuer_wrong))
    print("We got them wrong at these stages: ", stages_wrong)

    # Make ROC curves
    roc_unmodified = reads.get_roc()
    print("Roc unmodified", roc_unmodified)
    roc_adam_modified = reads.get_roc(map_q_fn=current_cap)
    print("Roc adam modified ", roc_adam_modified)
    start_time = time.time()
    roc_new_sum_modified = reads.get_roc(map_q_fn=proposed_cap)
    print("Roc mode modified (time:{}) ".format(time.time()-start_time), roc_new_sum_modified)
    Reads.plot_rocs([roc_unmodified, roc_adam_modified, roc_new_sum_modified])

    # plt.scatter([x.vg_computed_cap for x in reads.reads if not x.correct],
    # [x.faster_balanced_cap() for x in reads.reads if not x.correct])
    # plt.scatter([x.vg_computed_cap for x in reads.reads if x.correct],
    # [x.faster_cap() for x in reads.reads if x.correct], 'g^')
    # plt.show()
    
    #plt.clf()
    #plt.scatter([x.capped_map_q for x in reads.reads], [proposed_cap(x) for x in reads.reads])
    #plt.savefig('compare.png')
    
    #for read in reads.reads:
    #    if abs(read.capped_map_q - proposed_cap(read)) > 2:
    #        print("Got {} but vg got {} for {}".format(proposed_cap(read), read.capped_map_q, read))


if __name__ == '__main__':
    main()
