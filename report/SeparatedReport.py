__author__ = 'letovesnoi'

import os

from general import UtilsPipeline

from report import DistributionReport
from report import TXTMetricsReport


class SeparatedReport():
    """Class which generate extended report for each assembly separately"""

    def __init__(self, label, output_dir, transcripts_metrics, WELL_FULLY_COVERAGE_THRESHOLDS):
        # get folders for separated reports:
        self.output_dir = UtilsPipeline.create_empty_folder(os.path.join(output_dir, '{}_output'.format(label)))

        self.distribution_report = None

        self.txt_metrics_report = None

        # OTHER REPORTS:
        if transcripts_metrics.simple_metrics is not None:
            # UNALIGNED:
            self.path_fa_unaligned = os.path.join(self.output_dir, '{}.unaligned.fasta'.format(label))

            # MULTIPLE ALIGNED:
            self.path_fa_paralogous = os.path.join(self.output_dir, '{}.paralogs.fasta'.format(label))

            # MISASSEMBLED:
            self.path_fa_misassembled_together = os.path.join(self.output_dir, '{}.misassembled.fasta'.format(label))
            self.path_fa_misassembled_by_blat = os.path.join(self.output_dir, '{}.misassembled.blat.fasta'.format(label))
            self.path_fa_misassembled_by_blast = os.path.join(self.output_dir, '{}.misassembled.blast.fasta'.format(label))

            # UNIQUE ALIGNED:
            self.path_fa_unique_aligned = os.path.join(self.output_dir, '{}.correct.fasta'.format(label))

        if transcripts_metrics.assembly_completeness_metrics is not None:
            self.path_fully_assembled_genes = os.path.join(self.output_dir, '{}.{}%-assembled.genes'.format(label, int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_isoform_threshold * 100)))
            self.path_well_assembled_genes = os.path.join(self.output_dir, '{}.{}%-assembled.genes'.format(label, int(WELL_FULLY_COVERAGE_THRESHOLDS.well_isoform_threshold * 100)))
            self.path_fully_assembled_isoforms = os.path.join(self.output_dir, '{}.{}%-assembled.isoforms'.format(label, int(WELL_FULLY_COVERAGE_THRESHOLDS.fully_isoform_threshold * 100)))
            self.path_well_assembled_isoforms = os.path.join(self.output_dir, '{}.{}%-assembled.isoforms'.format(label, int(WELL_FULLY_COVERAGE_THRESHOLDS.well_isoform_threshold * 100)))


        if transcripts_metrics.assembly_correctness_metrics is not None:
            self.path_fa_unannotated = os.path.join(self.output_dir, '{}.unannotated.fasta'.format(label))


    def get_separated_report(self, args, label, transcripts_dict, transcripts_metrics, db_genes_metrics, reads_coverage, logger,
                             WELL_FULLY_COVERAGE_THRESHOLDS, PRECISION, TRANSCRIPT_LENS):
        logger.info()
        logger.info('Getting SEPARATED report for {}...'.format(label))

        self.txt_comparison_report = \
            TXTMetricsReport.TXTMetricsReport(args.blast, self.output_dir, [label], [transcripts_metrics], db_genes_metrics, reads_coverage, logger,
                                              WELL_FULLY_COVERAGE_THRESHOLDS, PRECISION, TRANSCRIPT_LENS)

        if not args.no_plots:
            self.distribution_report = \
                DistributionReport.DistributionReport([transcripts_metrics], db_genes_metrics, self.output_dir, logger, PRECISION)

        logger.print_timestamp('  ')
        logger.info('  Getting OTHER reports...')

        if transcripts_metrics.simple_metrics is not None:
            # UNALIGNED:
            transcripts_metrics.simple_metrics.print_fa_transcripts(transcripts_dict, transcripts_metrics.simple_metrics.ids_unaligned,
                                                                    self.path_fa_unaligned, logger, transcripts_name='Unaligned')

            # MULTIPLE ALIGNED:
            transcripts_metrics.simple_metrics.print_fa_transcripts(transcripts_dict, transcripts_metrics.simple_metrics.ids_mul_aligned,
                                                                    self.path_fa_paralogous, logger, transcripts_name='Paralogous')

            # MISASSEMBLED:
            # blat misassemblies intersection blast misassemblies:
            transcripts_metrics.simple_metrics.print_fa_transcripts(transcripts_dict, transcripts_metrics.simple_metrics.ids_misassembled_together,
                                                                    self.path_fa_misassembled_together, logger, transcripts_name='Misassembled')

            # blat misassemblies:
            transcripts_metrics.simple_metrics.print_fa_transcripts(transcripts_dict, transcripts_metrics.simple_metrics.ids_misassembled_by_blat,
                                                                    self.path_fa_misassembled_by_blat, logger, transcripts_name='Misassembled by BLAT')

            # blast misassemblies:
            transcripts_metrics.simple_metrics.print_fa_transcripts(transcripts_dict, transcripts_metrics.simple_metrics.ids_misassembled_by_blast,
                                                                    self.path_fa_misassembled_by_blast, logger, transcripts_name='Misassembled by BLASTN')

            # unique aligned:
            transcripts_metrics.simple_metrics.print_fa_transcripts(transcripts_dict, transcripts_metrics.simple_metrics.ids_unique_aligned,
                                                                    self.path_fa_unique_aligned, logger, transcripts_name='Unique aligned')

        if transcripts_metrics.assembly_completeness_metrics is not None:
            t_cov = transcripts_metrics.assembly_correctness_metrics.transcripts_coverage
            if t_cov:
                # UNANNOTATED:
                t_cov.print_unannotated_transcripts(transcripts_dict, self.path_fa_unannotated, logger)

            gene_db_cov = transcripts_metrics.assembly_completeness_metrics.isoforms_coverage
            if gene_db_cov:
                gene_db_cov.print_fully_assembled_genes(self.path_fully_assembled_genes, logger)
                gene_db_cov.print_well_assembled_genes(self.path_well_assembled_genes, logger)
                gene_db_cov.print_fully_assembled_isoforms(self.path_fully_assembled_isoforms, logger)
                gene_db_cov.print_well_assembled_isoforms(self.path_well_assembled_isoforms, logger)


        logger.info('  Done.')