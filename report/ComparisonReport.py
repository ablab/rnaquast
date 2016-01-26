__author__ = 'letovesnoi'

import os

from general import UtilsPipeline

from report import DistributionReport
from report import TXTMetricsReport


class ComparisonReport():
    """Class which generate distributions for all assemblies data sets at one from extended report"""

    def __init__(self, output_dir):

        # get folders for compairison report:
        self.output_dir = UtilsPipeline.create_empty_folder(os.path.join(output_dir, 'comparison_output'))

        self.txt_comparison_report = None

        self.distribution_report = None


    def get_comparison_report(self, args, transcripts_metrics, db_genes_metrics, reads_coverage, logger,
                              WELL_FULLY_COVERAGE_THRESHOLDS, PRECISION, TRANSCRIPT_LENS):

        logger.print_timestamp()
        logger.info('Getting COMPARISON report...')

        self.txt_comparison_report = \
            TXTMetricsReport.TXTMetricsReport(self.output_dir, transcripts_metrics, db_genes_metrics, reads_coverage, logger,
                                              WELL_FULLY_COVERAGE_THRESHOLDS, PRECISION, TRANSCRIPT_LENS)

        if not args.no_plots:
            self.distribution_report = \
                DistributionReport.DistributionReport(transcripts_metrics, db_genes_metrics, self.output_dir, logger,
                                                      PRECISION)

        logger.info('  saved to {}'.format(self.output_dir))