__author__ = 'lenk'

from datetime import datetime

from general import UtilsAlignment
from general import UtilsCoverage

from objects import AlignedTranscript

from metrics import BasicTranscriptsMetrics
from metrics import SimpleTranscriptsMetrics
from metrics import AssemblyCompletenessMetrics
from metrics import InternalIsoformsCoverage
from metrics import AssemblyCorrectnessMetrics
from metrics import OneTranscriptCoverage


class TranscriptsMetrics():
    """Class of metrics of assembled transcripts with and without alignments and annotations"""

    def __init__(self, args, label):
        self.label = label

        # METRICS WITHOUT ALIGNMENT:
        self.basic_metrics = None

        # METRICS WITH ALIGNMENT:
        self.simple_metrics = None

        # METRCIS WITH ALIGNMENT AND ANNOTATION
        # coverages of aligned transcripts by annotated isoform:
        self.assembly_correctness_metrics = None
        # coverages of annotated isoforms by aligned transcripts:
        self.assembly_completeness_metrics = None
        # dictionary of transcripts coverages:
        self.transcripts_coverage_dict = None

        # INITIALIZE BASIC METRICS WITHOUT ALIGNMENT:
        if args.transcripts is not None:
            self.basic_metrics = BasicTranscriptsMetrics.BasicTranscriptsMetrics()

        # INITIALIZE SIMPLE METRICS WITH ALIGNMENT:
        if args.alignment is not None and args.reference is not None and args.transcripts is not None:
            self.simple_metrics = SimpleTranscriptsMetrics.SimpleTranscriptsMetrics()

            # temporary:
            args.fusion_misassemble_analyze = None
            # self.fusion_misassemble_metrics = FusionMisassembleMetrics.FusionMisassembleMetrics(args.fusion_misassemble_analyze, self, transcripts_file, sam_file, logger)

        # INITIALIZE METRICS WITH ALIGNMENT AND ANNOTATION:
        # ASSEMBLY CORRECTNESS METRICS
        self.assembly_correctness_metrics = AssemblyCorrectnessMetrics.AssemblyCorrectnessMetrics(args)

        # ASSEMBLY COMPLETENESS METRICS:
        # metrics of coverages of annotated isoforms by aligned transcripts:
        self.assembly_completeness_metrics = \
            AssemblyCompletenessMetrics.AssemblyCompletenessMetrics(args)


    def get_transcripts_metrics(self, args, type_organism, reference_dict, transcripts_path, transcripts_dict, label, threads,
                                sqlite3_db_genes, db_genes_metrics, reads_coverage, logger, tmp_dir, log_dir,
                                WELL_FULLY_COVERAGE_THRESHOLDS, TRANSCRIPT_LENS):
        logger.print_timestamp('  ')

        # GET BASIC TRANSCRIPTS METRICS:
        if self.basic_metrics is not None:
            self.basic_metrics.get_basic_metrics(args, transcripts_dict, logger, TRANSCRIPT_LENS)

        # GET ALIGNMENT METRICS:
        if self.simple_metrics is not None:
            self.simple_metrics.get_metrics(args.blast, reference_dict, transcripts_dict, self.basic_metrics.number,
                                            logger)

        # GET ASSEMBLY CORRECTNESS METRICS:
        if self.assembly_correctness_metrics is not None:
            self.assembly_correctness_metrics.get_assembly_correctness_metrics(self.simple_metrics, logger)

        if self.assembly_completeness_metrics is not None:
            self.assembly_completeness_metrics. \
                get_assembly_completeness_metrics(args, sqlite3_db_genes, db_genes_metrics, reads_coverage,
                                                  transcripts_path, type_organism, tmp_dir, label, threads,
                                                  WELL_FULLY_COVERAGE_THRESHOLDS, logger, log_dir)


    def processing_assembled_psl_file(self, assembled_psl_file, sorted_exons_attr, strand_specific, logger,
                                      sqlite3_db_genes, type_isoforms, WELL_FULLY_COVERAGE_THRESHOLDS):
        init_time = datetime.now()
        init_time -= init_time
        simple_time = init_time
        assembly_correctness_time = init_time
        assembly_completeness_time = init_time
        mapped_coverage_time = init_time
        best_mapped_time = init_time
        transcript_time = init_time

        logger.print_timestamp('  ')
        logger.info('  Processing assembled aligned transcripts...')

        with open(assembled_psl_file, 'r') as fin:
            line1 = fin.readline().strip()
            line2 = fin.readline().strip()
            line_count = 0
            while line1 != '':
                logger.print_timestamp('Line #' + str(line_count))
                prev_time_stamp = datetime.now()

                best_lines, best_alignments, line1, line2 = \
                    UtilsAlignment.get_curr_single_transcript_lines_alignments(0, line1, line2, fin, True)

                logger.info('get_curr_single_transcript_lines_alignments DONE in ' + str(datetime.now() - prev_time_stamp))
                prev_time_stamp = datetime.now()
                # GET BEST MAPPED ALIGNMENTS:
                # in case when we havn't annotation:
                best_mapped_lines, best_mapped_alignments, best_mapped_aligned_transcripts, \
                best_mapped_aligned_transcripts_coverages, best_mapped_internal_isoforms_coverages,\
                curr_best_mapped_time, curr_transcript_time = \
                    self.get_best_mapped_from_best_aligned(best_lines, best_alignments, sorted_exons_attr,
                                                           strand_specific, sqlite3_db_genes, type_isoforms,
                                                           WELL_FULLY_COVERAGE_THRESHOLDS, logger)

                logger.info('get_best_mapped_from_best_aligned DONE in ' + str(datetime.now() - prev_time_stamp))
                prev_time_stamp = datetime.now()

                best_mapped_time += curr_best_mapped_time
                transcript_time += curr_transcript_time

                # FILTERING OVER MAPPED ISOFORMS LENGTHS:
                # if isoforms_len_range != None:
                #     filtered_lines = []
                #     filtered_alignments = []
                #     filtered_aligned_transcripts = []
                #     filtered_aligned_transcripts_coverages = []
                #     for i_alignment in range(len(best_mapped_alignments)):
                #         id_chr = best_mapped_aligned_transcripts[i_alignment].alignment.target_fragment.name
                #         strand = best_mapped_aligned_transcripts[i_alignment].strand
                #         id_isoform = best_mapped_aligned_transcripts_coverages[i_alignment].id_mapped_isoform
                #         if id_isoform == None:
                #             continue
                #         elif annotated_isoforms[strand][id_chr].len_wout_introns_dict[id_isoform] >= isoforms_len_range[0] and \
                #                         annotated_isoforms[strand][id_chr].len_wout_introns_dict[id_isoform] <= isoforms_len_range[1]:
                #             filtered_lines.append(best_mapped_lines[i_alignment])
                #             filtered_alignments.append(best_mapped_alignments[i_alignment])
                #             filtered_aligned_transcripts.append(best_mapped_aligned_transcripts[i_alignment])
                #             filtered_aligned_transcripts_coverages.append(best_mapped_aligned_transcripts_coverages[i_alignment])
                #     best_mapped_lines = filtered_lines
                #     best_mapped_alignments = filtered_alignments
                #     best_mapped_aligned_transcripts = filtered_aligned_transcripts
                #     best_mapped_aligned_transcripts_coverages = filtered_aligned_transcripts_coverages

                if self.simple_metrics is not None:
                    simple_time += self.simple_metrics.update_metrics_by_best_mapped_alignments(best_mapped_alignments)

                for i_alignment in range(len(best_mapped_alignments)):
                    # UPDATE SIMPLE TRANSCRIPTS METRICS:

                    if self.simple_metrics is not None:
                        # update metrics of assembled transcripts with alignments by best single alignment:
                        simple_time += self.simple_metrics.update_metrics_by_best_mapped_transcript\
                            (best_mapped_aligned_transcripts[i_alignment])

                    # SET COVERAGES:
                    if self.assembly_correctness_metrics.transcripts_coverage is not None and self.assembly_completeness_metrics.isoforms_coverage is not None:
                        # update coverage of transcripts:
                        assembly_correctness_time += \
                            self.assembly_correctness_metrics.update_assembly_correctness_metrics\
                            (best_mapped_aligned_transcripts[i_alignment],
                             best_mapped_aligned_transcripts_coverages[i_alignment], WELL_FULLY_COVERAGE_THRESHOLDS)

                        # update coverage of annotations:
                        id_isoform = best_mapped_aligned_transcripts_coverages[i_alignment].id_mapped_isoform

                        if id_isoform is not None:
                            assembly_completeness_time += self.assembly_completeness_metrics.\
                                update_assembly_completeness_metrics(sqlite3_db_genes, best_mapped_internal_isoforms_coverages[i_alignment], id_isoform)

                    logger.info('> Processed alignment # ' + str(i_alignment) + ' in ' + str(datetime.now() - prev_time_stamp))
                    prev_time_stamp = datetime.now()

                logger.info('Processed ' + str(len(best_mapped_alignments)) + ' alignments in ' + str(datetime.now() - prev_time_stamp))
                prev_time_stamp = datetime.now()

                logger.print_timestamp('Line #' + str(line_count) + ' done')
                line_count += 1

        logger.info('  Done.')

        logger.debug('  ELAPSED TIME: ' + str(simple_time + assembly_correctness_time + assembly_completeness_time +
                                              mapped_coverage_time + best_mapped_time))
        logger.debug('    simple time: ' + str(simple_time))
        logger.debug('    assembly correctness time: ' + str(assembly_correctness_time))
        logger.debug('    assembly completeness time: ' + str(assembly_completeness_time))
        logger.debug('    mapped coverage time: ' + str(mapped_coverage_time))
        logger.debug('    best mapped time: ' + str(best_mapped_time))
        logger.debug('    transcript time: ' + str(transcript_time))


    def processing_misassembled_psl_file(self, misassembled_union_file, logger, is_psl_lines=True):
        logger.print_timestamp('  ')
        logger.info('  Processing misassembled aligned transcripts...')

        with open(misassembled_union_file, 'r') as fin:
            line1 = fin.readline().strip()
            line2 = fin.readline().strip()
            while line1 != '':
                best_lines, best_alignments, line1, line2 = \
                    UtilsAlignment.get_curr_single_transcript_lines_alignments(0, line1, line2, fin, is_psl_lines)

                self.simple_metrics.update_metrics_by_misassembled_alignments(best_alignments, is_psl_lines)

        logger.info('  Done.')


    def get_best_mapped_from_best_aligned(self, best_lines, best_alignments, sorted_exons_attr, strand_specific,
                                          sqlite3_db_genes, type_isoforms, WELL_FULLY_COVERAGE_THRESHOLDS, logger):
        start_time = datetime.now()
        prev_time_stamp = datetime.now()

        best_aligned_transcripts = []
        best_aligned_transcripts_coverages = []
        best_aligned_internal_isoforms_coverages = []
        for i_line_alignment in range(len(best_lines)):
            prev_for_time_stamp = datetime.now()
            curr_aligned_transcript, curr_aligned_transcript_coverage, curr_internal_isoforms_coverage, \
            elapsed_transcript_time = \
                self.get_aligned_transcript_and_coverages(best_alignments[i_line_alignment], sorted_exons_attr,
                                                          strand_specific, sqlite3_db_genes, type_isoforms,
                                                          WELL_FULLY_COVERAGE_THRESHOLDS, logger)

            logger.info('>> alignment #' + str(i_line_alignment) + ' get_aligned_transcript_and_coverages in ' + str(datetime.now() - prev_for_time_stamp))
            prev_for_time_stamp = datetime.now()

            best_aligned_transcripts.append(curr_aligned_transcript)

            best_aligned_transcripts_coverages.append(curr_aligned_transcript_coverage)

            best_aligned_internal_isoforms_coverages.append(curr_internal_isoforms_coverage)
    
            logger.info('alignment #' + str(i_line_alignment) + ' added in ' + str(datetime.now() - prev_for_time_stamp))
            prev_for_time_stamp = datetime.now()

        logger.info('> get_aligned_transcript_and_coverages DONE in ' + str(datetime.now() - prev_time_stamp))
        prev_time_stamp = datetime.now()

        # IN CASE WHEN WE HAVN'T ANNOTATION:
        if sqlite3_db_genes is None:
            elapsed_time = datetime.now() - start_time

            return best_lines, best_alignments, best_aligned_transcripts, \
                   best_aligned_transcripts_coverages, best_aligned_internal_isoforms_coverages, \
                   elapsed_time, elapsed_transcript_time

        # IN CASE WHEN WE HAVE ANNOTATION:
        else:
            # choose best annotated transcript alignments over all best alignments:
            transcripts_covered_bases = {}
            isoforms_covered_fraction = {}
            for i_alignment in range(len(best_aligned_transcripts_coverages)):
                curr_id_isoform = best_aligned_transcripts_coverages[i_alignment].id_mapped_isoform
                curr_isoforms_coverage = best_aligned_internal_isoforms_coverages[i_alignment]
                if curr_id_isoform is not None:
                    transcripts_covered_bases[i_alignment] = best_aligned_transcripts_coverages[i_alignment].covered_bases[curr_id_isoform]
                    isoforms_covered_fraction[i_alignment] = curr_isoforms_coverage.assembled_fraction[curr_id_isoform]

            curr_max_keys = UtilsCoverage.get_ids_best_mapped(transcripts_covered_bases, isoforms_covered_fraction)

            # for transcripts all aligned to unannotated regions:
            if curr_max_keys == []:
                elapsed_time = datetime.now() - start_time

                return best_lines, best_alignments, best_aligned_transcripts, \
                       best_aligned_transcripts_coverages, best_aligned_internal_isoforms_coverages,\
                       elapsed_time, elapsed_transcript_time

            # form lines, alignments, transcripts and coverages corresponded best mapping:
            best_mapped_lines = []
            best_mapped_alignments = []
            best_mapped_aligned_transcripts = []
            best_mapped_aligned_transcripts_coverages = []
            best_mapped_internal_isoforms_coverages = []
            for i_line_alignment in curr_max_keys:
                best_mapped_lines.append(best_lines[i_line_alignment])
                best_mapped_alignments.append(best_alignments[i_line_alignment])
                best_mapped_aligned_transcripts.append(best_aligned_transcripts[i_line_alignment])
                best_mapped_aligned_transcripts_coverages.append(best_aligned_transcripts_coverages[i_line_alignment])
                best_mapped_internal_isoforms_coverages.append(best_aligned_internal_isoforms_coverages[i_line_alignment])

            elapsed_time = datetime.now() - start_time

            logger.info('> choosen best annotated transcript in ' + str(datetime.now() - prev_time_stamp))
            prev_time_stamp = datetime.now()

            return best_mapped_lines, best_mapped_alignments, best_mapped_aligned_transcripts, \
                   best_mapped_aligned_transcripts_coverages, best_mapped_internal_isoforms_coverages,\
                   elapsed_time, elapsed_transcript_time


    # get aligned transcript, transcript coverage and internal isoforms coverage:
    def get_aligned_transcript_and_coverages(self, psl_alignment, sorted_exons_attr, strand_specific, sqlite3_db_genes,
                                             type_isoforms, WELL_FULLY_COVERAGE_THRESHOLDS, logger):

        prev_time_stamp = datetime.now()
        # CREATE ALIGNED TRANSCRIPT:
        start_time = datetime.now()

        # print 'aligned transcript: ', datetime.now()
        aligned_transcript = AlignedTranscript.AlignedTranscript(psl_alignment, sorted_exons_attr, strand_specific,
                                                                 sqlite3_db_genes, type_isoforms)

        logger.info('>>> got aligned transcript in ' + str(datetime.now() - prev_time_stamp))
        prev_time_stamp = datetime.now()
        # print 'done: ', datetime.now()

        elapsed_transcript_time = datetime.now() - start_time

        # GET COVERAGES:
        aligned_transcript_coverage = None
        internal_isoforms_coverage = None
        if self.assembly_correctness_metrics.transcripts_coverage is not None and \
                        self.assembly_completeness_metrics.isoforms_coverage is not None:
            # get coverages:
            aligned_transcript_coverage = \
                OneTranscriptCoverage.OneTranscriptCoverage(aligned_transcript.ids_internal_isoforms,
                                                            aligned_transcript.alignment.blocks_num)
            # print 'internal isoforms coverage: : ', datetime.now()
            internal_isoforms_coverage = \
                InternalIsoformsCoverage.InternalIsoformsCoverage(aligned_transcript.internal_isoforms)
            # print 'done: : ', datetime.now()

            logger.info('>>> calculated coverage in ' + str(datetime.now() - prev_time_stamp))
            prev_time_stamp = datetime.now()


            # print 'bases exons blocks covered: ', datetime.now()
            # set exons and blocks overlap (covered bases):
            for internal_isoform in aligned_transcript.internal_isoforms:
                # exon.start, exon.end: 1-based coordinates; start must be <= end
                exons_starts = [exon.start - 1 for exon in aligned_transcript.children_exons_dict[internal_isoform.id]]
                exons_ends = [exon.end - 1 for exon in aligned_transcript.children_exons_dict[internal_isoform.id]]
                exons_ids = [exon.id for exon in aligned_transcript.children_exons_dict[internal_isoform.id]]

                target_cov_pos, query_cov_pos = \
                    UtilsCoverage.get_coverage_positions(exons_ids, exons_starts, exons_ends, range(aligned_transcript.alignment.blocks_num),
                                                         aligned_transcript.alignment.target_fragment.starts, aligned_transcript.alignment.target_fragment.ends)

                internal_isoforms_coverage.update_internal_isoforms_coverage(sqlite3_db_genes, internal_isoform.id, target_cov_pos)

                aligned_transcript_coverage.update_transcript_coverage(internal_isoform.id, query_cov_pos)

            # print 'done: ', datetime.now()
            logger.info('>>> calculated coverage positions for ' + str(len(aligned_transcript.internal_isoforms)) + ' internal isoforms in ' + str(datetime.now() - prev_time_stamp))
            prev_time_stamp = datetime.now()

            internal_isoforms_coverage.get_internal_isoforms_coverage(aligned_transcript.internal_isoforms,
                                                                      aligned_transcript.children_exons_dict)

            aligned_transcript_coverage.get_transcript_coverage(aligned_transcript, internal_isoforms_coverage,
                                                                WELL_FULLY_COVERAGE_THRESHOLDS)
            logger.info('>>> finalized in ' + str(datetime.now() - prev_time_stamp))
            prev_time_stamp = datetime.now()
        return aligned_transcript, aligned_transcript_coverage, internal_isoforms_coverage, elapsed_transcript_time
