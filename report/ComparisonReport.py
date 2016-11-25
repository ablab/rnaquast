__author__ = 'letovesnoi'

import os

from general import UtilsPipeline

from report import DistributionReport
from report import TXTMetricsReport


class ComparisonReport():
    """Class which generate distributions for all assemblies data sets at one from extended report"""

    def __init__(self):

        # get folders for compairison report:
        self.output_dir = None

        self.txt_comparison_report = None

        self.path_well_expressed_list_by_reads = None
        self.path_fully_expressed_list_by_reads = None

        self.distribution_report = None


    def get_comparison_report(self, args, output_dir, labels, transcripts_metrics, db_genes_metrics, reads_coverage, logger,
                              WELL_FULLY_COVERAGE_THRESHOLDS, PRECISION, TRANSCRIPT_LENS):

        logger.print_timestamp()
        logger.info('Getting COMPARISON report...')

        if len(transcripts_metrics) != 0:
            self.output_dir = UtilsPipeline.create_empty_folder(os.path.join(output_dir, 'comparison_output'))
        else:
            self.output_dir = output_dir

        self.txt_comparison_report = \
            TXTMetricsReport.TXTMetricsReport(args.blast, self.output_dir, labels, transcripts_metrics, db_genes_metrics, reads_coverage, logger,
                                              WELL_FULLY_COVERAGE_THRESHOLDS, PRECISION, TRANSCRIPT_LENS)

        if reads_coverage is not None:
            self.path_well_expressed_list_by_reads = os.path.join(self.output_dir, 'reads.{}%-covered.list'.format(int(WELL_FULLY_COVERAGE_THRESHOLDS.well_isoform_threshold * 100)))
            reads_coverage.print_well_expressed_isoforms(self.path_well_expressed_list_by_reads, logger)

            self.path_fully_expressed_list_by_reads = os.path.join(self.output_dir, 'reads.{}%-covered.list'.format(int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_isoform_threshold * 100)))
            reads_coverage.print_fully_expressed_isoforms(self.path_fully_expressed_list_by_reads, logger)


        if not args.no_plots:
            self.distribution_report = \
                DistributionReport.DistributionReport(transcripts_metrics, db_genes_metrics, self.output_dir, logger,
                                                      PRECISION)

        logger.info('  saved to {}'.format(self.output_dir))