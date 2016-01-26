__author__ = 'letovesnoi'

from quast23.libs import N50


class BasicTranscriptsMetrics():
    """Class of basic assembled transcripts metrics without alignment"""

    def __init__(self):
        # total number of transcripts:
        self.number = 0

        # length of transcript and number of assembled transcripts having this length over all transcripts:
        self.len_distribution = {}
        # len of transcripts:
        self.min_len = 0
        self.max_len = 0
        self.avg_len = 0.0
        self.tot_len = 0

        # N50:
        self.n50 = 0
        self.lengths = []

        # number of transcripts longer than 500 bp:
        self.num_transcripts_500 = 0
        # number of transcripts longer than 1000 bp:
        self.num_transcripts_1000 = 0

        # INITIALIZE METRICS WITH ANNOTATION:
        # if args_annotation != None:
            # number of transcripts having length more then min isoform without introns length and less then max isoform without introns length:
            # self.num_consistent_len = 0
            # number of transcripts having length out isoform range:
            # self.num_outliers_len = 0
            # percent of transcripts having length out isoform range:
            # self.percent_outliers_len = 0.0


    def update_metrics(self, assembled_transcript_seq, TRANSCRIPT_LENS):
        self.number += 1
        len_assembled_transcript = len(assembled_transcript_seq)
        self.tot_len += len_assembled_transcript
        if len_assembled_transcript not in self.len_distribution:
            self.len_distribution[len_assembled_transcript] = 0
        self.len_distribution[len_assembled_transcript] += 1
        self.lengths.append(len_assembled_transcript)

        if len_assembled_transcript >= TRANSCRIPT_LENS[0]:
            self.num_transcripts_500 += 1

        if len_assembled_transcript >= TRANSCRIPT_LENS[1]:
            self.num_transcripts_1000 += 1

        # if args.annotation != None:
        #     self.update_metrics_w_annotation(assembled_transcript_seq, basic_isoforms_metrics)


    # def update_metrics_w_annotation(self, assembled_transcript, basic_isoforms_metrics):
    #     len_assembled_transcript = len(assembled_transcript)
    #     if len_assembled_transcript > basic_isoforms_metrics.min_len_wout_introns and len_assembled_transcript < basic_isoforms_metrics.max_len_wout_introns:
    #         self.num_consistent_len += 1


    def get_basic_metrics(self, args, transcripts_dict, logger, TRANSCRIPT_LENS):
        logger.info('  Getting BASIC TRANSCRIPTS metrics...')
        for id_transcript in transcripts_dict:
            self.update_metrics(transcripts_dict[id_transcript], TRANSCRIPT_LENS)
        # get averages:
        if self.number != 0:
            self.avg_len = self.tot_len * 1.0 / self.number
            self.max_len = max(self.len_distribution.keys())
            self.min_len = min(self.len_distribution.keys())

        self.n50 = N50.N50(self.lengths)

        # if args.annotation != None:
        #     self.get_metrics_w_annotation()

        logger.info('  Done.')


    # def get_metrics_w_annotation(self):
    #     if self.number != 0:
    #         self.num_outliers_len = self.number - self.num_consistent_len
    #         self.percent_outliers_len = self.num_outliers_len * 1.0 / self.number