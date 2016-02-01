#!/usr/bin/env python

__author__ = 'lenk'

import os
import sys
import shutil

rquast_dirpath = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))

sys.path.append(os.path.join(rquast_dirpath, 'quast23', 'libs'))
sys.path.append(os.path.join(rquast_dirpath, 'quast23'))
sys.path.append(os.path.join(rquast_dirpath, 'libs'))
sys.path.append(os.path.join(rquast_dirpath, 'libs', 'joblib'))

import fastaparser

from general import rqconfig
from general import log

from general import UtilsGeneral
from general import UtilsPipeline
from general import UtilsTools
from general import UtilsAlignment
from general import UtilsAnnotations

from objects import SortedExonsAttributes

from metrics import TranscriptsMetrics
from metrics import GeneDatabaseMetrics
from metrics import ReadsCoverage

from report import ShortReport
from report import SeparatedReport
from report import ComparisonReport

logger = log.get_logger(rqconfig.LOGGER_DEFAULT_NAME)

from general.rqconfig import PRECISION


def main_utils():
    program_name = sys.argv[0][:sys.argv[0].rfind('.')]

    # parse running string of main program and get all arguments:
    args = UtilsPipeline.get_arguments()

    WELL_FULLY_COVERAGE_THRESHOLDS = rqconfig.well_fully_coverage_thresholds(args.lower_threshold, args.upper_threshold)

    ALIGNMENT_THRESHOLDS = rqconfig.alignment_thresholds()

    # run rnaQUAST on test_data:
    if args.test:
        UtilsPipeline.run_rnaQUAST_on_test_data(args, rquast_dirpath, program_name)
        # UtilsPipeline.run_rnaQUAST_on_debug_data(args, rquast_dirpath, program_name)
        sys.exit()

    # create output directory:
    args.output_dir = UtilsPipeline.create_output_folder(args.output_dir, program_name)
    # create temporary directory:
    tmp_dir = UtilsPipeline.create_empty_folder(os.path.join(args.output_dir, 'tmp'))

    # SET LOGGER:
    if args.debug:
        rqconfig.debug = True
        logger.set_up_console_handler(debug=True)
    else:
        logger.set_up_console_handler()
    logger.set_up_file_handler(args.output_dir)
    logger.print_command_line([os.path.realpath(__file__)] + sys.argv[1:], wrap_after=None)
    logger.start()

    # THREADING:
    args.threads = UtilsPipeline.get_num_threads(args.threads, logger)


    # GET SINGLE FASTA REFERENCE FILE:
    fa_flag = False
    args.database = None
    if args.reference is not None:
        list_ext_fasta = ['fa', 'fasta', 'fna', 'ffn', 'frn', 'fsa']
        # for file with list of pathes to scaffolds/patches/chromosomes:
        ext_database = 'txt'

        for ext in list_ext_fasta:
            if ext in args.reference:
                fa_flag = True
        if ext_database in args.reference:
            args.database = args.reference
            args.reference = UtilsGeneral.glue_scaffolds_together(args.reference, args.database, tmp_dir, logger)
        elif fa_flag == False:
            logger.warning('Strange FASTA extension.')

    # READ REFERENCE FROM MULTIFASTA:
    reference_dict = None
    ids_chrs = None
    if args.reference is not None:
        logger.print_timestamp()
        logger.info('Getting reference...')
        reference_dict = UtilsGeneral.list_to_dict(fastaparser.read_fasta(args.reference))
        logger.info('Done.')

        genome_len = UtilsGeneral.get_genome_len(reference_dict)

        ids_chrs = reference_dict.keys()

        # correction for fasta contained Y, W and etc:
        for id_chr in ids_chrs:
            reference_dict[id_chr] = UtilsGeneral.correct_nucl_seq(reference_dict[id_chr])


    # for strand specific data we store + and - keys in dictionaries and only + for non strand specific data:
    strands = UtilsGeneral.get_strands(args, logger)

    type_organism = 'eukaryotes'
    if args.prokaryote:
        type_organism = 'prokaryotes'

    # USE ANNOTATION:
    sqlite3_db_genes = None
    sorted_exons_attr = None
    db_genes_metrics = None
    type_genes, type_isoforms, type_exons = UtilsAnnotations.default_type_genes, UtilsAnnotations.default_type_isoforms, UtilsAnnotations.default_type_exons
    if args.gene_database is not None:
        annotation_name = os.path.split(args.gene_database)[1]
        annotation_label = annotation_name[:annotation_name.rfind('.g')]

        args.gene_database = UtilsAnnotations.clear_gtf_by_reference_chr(args.gene_database, ids_chrs, tmp_dir, annotation_label, logger)

        sqlite3_db_genes = UtilsAnnotations.create_sqlite3_db(args.gene_database, annotation_label, args.disable_infer_genes,
                                              args.disable_infer_transcripts, args.store_db, args.output_dir, tmp_dir,
                                              logger)

        type_genes, type_isoforms, type_exons = \
            UtilsAnnotations.get_type_features(sqlite3_db_genes, UtilsAnnotations.default_type_genes,
                                               UtilsAnnotations.default_type_isoforms,
                                               UtilsAnnotations.default_type_exons, logger)

        if UtilsAnnotations.default_type_exons == type_exons:
            type_organism = 'eukaryotes'
        else:
            type_organism = 'prokaryotes'

        db_genes_metrics = GeneDatabaseMetrics.GeneDatabaseMetrics(sqlite3_db_genes, type_genes, type_isoforms, logger)

        ALIGNMENT_THRESHOLDS.ERR_SPACE_TARGET_FAKE_BLAT = db_genes_metrics.max_intron_len + 100
        logger.info('\nSets maximum intron size equal {}. Default is 1500000 bp.\n'.format(ALIGNMENT_THRESHOLDS.ERR_SPACE_TARGET_FAKE_BLAT))

        # set exons starts / ends and ids for binning strategy:
        sorted_exons_attr = \
            SortedExonsAttributes.SortedExonsAttributes(sqlite3_db_genes, type_exons, strands, ids_chrs, reference_dict, logger)

    reads_coverage = None
    if args.reads_alignment is not None or \
            ((args.single_reads is not None or (args.left_reads is not None and args.right_reads is not None))
             and args.reference is not None and args.gene_database is not None):
        reads_coverage = \
            ReadsCoverage.ReadsCoverage(args.reads_alignment, args.tophat, args.reference, args.single_reads,
                                        args.left_reads, args.right_reads, reference_dict, sqlite3_db_genes, type_isoforms,
                                        sorted_exons_attr, args.strand_specific, db_genes_metrics.tot_isoforms_len,
                                        genome_len, tmp_dir, args.threads, WELL_FULLY_COVERAGE_THRESHOLDS, logger)


    if args.transcripts is not None:
        # GET TRANSCRIPTS:
        transcripts_dicts = []
        for i_transcripts in range(len(args.transcripts)):
            logger.print_timestamp('  ')
            logger.info('  Getting transcripts from {}...'.format(args.transcripts[i_transcripts]))
            transcripts_dicts.append(UtilsGeneral.list_to_dict(fastaparser.read_fasta(args.transcripts[i_transcripts])))
            logger.info('  Done.')

        # get labels for folders names and names of transcripts in reports:
        all_labels_from_dirs = False
        if args.labels is None:
            args.labels = UtilsPipeline.process_labels(args.transcripts, args.labels, all_labels_from_dirs)
    else:
        logger.error('Usage: --transcripts TRANSCRIPTS', exit_with_code=2, to_stderr=True)
        sys.exit(2)


    # GET PSL ALIGNMENT FILE:
    if args.alignment is None and args.reference is not None and args.transcripts is not None:
        if args.blat:
            args.alignment = UtilsTools.run_blat(args.database, args.reference, transcripts_dicts, args.labels, args.threads, tmp_dir, logger)
        else:
            args.alignment = UtilsTools.run_gmap(args.reference, genome_len, args.transcripts, args.labels, args.threads, tmp_dir, logger)

        #if args.fusion_misassemble_analyze:
        #    if not (args.left_reads is not None and args.right_reads is not None):
        #        logger.error('Usage: --left_reads LEFT_READS --right RIGHT_READS for analyse fusions and misassemblies',
        #                     exit_with_code=2, to_stderr=True)
        #        sys.exit(2)


    # FOR MISASSEMBLIES SEARCH:
    # GET DATABASE FOR FA ISOFORMS:
    args.blast = False
    if args.reference is not None and args.gene_database is not None and args.alignment is not None:
        blastn_run = os.path.join(rqconfig.rnaOUAST_LOCATION, '.', 'blastn')
        if not os.path.isfile(blastn_run):
            blastn_run = "blastn"

        if UtilsGeneral.which(blastn_run) is None:
            logger.warning('blastn not found! Please add blastn to PATH for better MISASSEMBLIES metrics.')
        else:
            args.blast = True

            isoforms_fa_path = os.path.join(tmp_dir, '{}.isoforms.fa'.format(annotation_label))
            isoforms_list = UtilsGeneral.dict_to_list(UtilsAnnotations.get_fa_isoforms(sqlite3_db_genes, type_isoforms, type_exons, reference_dict, logger))
            fastaparser.write_fasta(isoforms_fa_path, isoforms_list)

            isoforms_blast_db = UtilsTools.get_blast_db(isoforms_fa_path, annotation_label, tmp_dir, logger)


    # LOGGING INPUT DATA:
    logger.print_input_files(args)


    # INITIALIZATION TRANSCRIPTS METRICS AND REPORTS:
    transcripts_metrics = []
    separated_reports = []
    for i_transcripts in range(len(args.transcripts)):
        # INITIALIZE TRANSCRIPTS METRICS:
        #if args.sam_file is not None:
        #    sam_file_tmp = args.sam_file[i_transcripts]
        #else:
        transcripts_metrics.append(TranscriptsMetrics.TranscriptsMetrics(args, args.labels[i_transcripts]))

        # INITIALIZE SEPARATED REPORTS:
        separated_reports.append(SeparatedReport.SeparatedReport(args.labels[i_transcripts], args.output_dir, transcripts_metrics[i_transcripts], WELL_FULLY_COVERAGE_THRESHOLDS))

    if args.transcripts is not None:
        '''from joblib import Parallel, delayed

        n = len(args.transcripts)
        run_n = n / args.threads
        for i_run in range(run_n):
            tmp = Parallel(n_jobs=args.threads)(delayed(process_one_trascripts_file)(args, i_transcripts, reference_dict, annotation_dict,
                                                                                          annotated_exons, annotated_isoforms, strands, transcripts_metrics,
                                                                                          basic_isoforms_metrics, separated_reports)
                                                     for i_transcripts in range(i_run * args.threads, args.threads * (i_run + 1), 1))
            for i in range(args.threads):
                i_transcripts = i + i_run * args.threads
                transcripts_metrics[i_transcripts] = tmp[i][0]
                separated_reports[i_transcripts] = tmp[i][1]

        if n - run_n * args.threads != 0:
            tmp = Parallel(n_jobs=n - run_n * args.threads)(delayed(process_one_trascripts_file)(args, i_transcripts, reference_dict, annotation_dict,
                                                                                                 annotated_exons, annotated_isoforms, strands, transcripts_metrics,
                                                                                                 basic_isoforms_metrics, separated_reports)
                                                            for i_transcripts in range(run_n * args.threads, n, 1))
            for i in range(n - run_n * args.threads):
                i_transcripts = i + run_n * args.threads
                transcripts_metrics[i_transcripts] = tmp[i][0]
                separated_reports[i_transcripts] = tmp[i][1]'''

        alignments_reports = []
        blast_alignments = []
        for i_transcripts in range(len(args.transcripts)):
            logger.info()
            logger.info('Processing transcripts from {}:'.format(args.transcripts[i_transcripts]))

            if args.blast == True:
                blast_alignments.append\
                    (UtilsTools.align_transcripts_to_isoforms_by_blastn
                     (args.transcripts[i_transcripts], isoforms_blast_db, tmp_dir, args.labels[i_transcripts], logger))
            else:
                blast_alignments.append(None)

            # PROCESS TRANSCRIPTS ALIGNMENTS:
            if transcripts_metrics[i_transcripts].simple_metrics is not None:
                # GET FILES WITH ALIGNMENTS REPORTS:
                alignments_reports.append\
                    (UtilsAlignment.AlignmentsReport.get_alignments_report
                     (args.labels[i_transcripts], args.alignment[i_transcripts], blast_alignments[i_transcripts],
                      transcripts_dicts[i_transcripts], tmp_dir, args.min_alignment, logger, ALIGNMENT_THRESHOLDS))

                # UPDATE METRICS BY ASSEMBLED TRANSCRIPTS:
                transcripts_metrics[i_transcripts].processing_assembled_psl_file\
                    (alignments_reports[i_transcripts].blat_report.assembled_psl_file, sorted_exons_attr,
                     args.strand_specific, logger, sqlite3_db_genes, type_isoforms, WELL_FULLY_COVERAGE_THRESHOLDS)

                # UPDATE METRICS BY MISASSEMBLED TRANSCRIPTS:
                # by blat:
                transcripts_metrics[i_transcripts].processing_misassembled_psl_file\
                    (alignments_reports[i_transcripts].blat_report.misassembled_psl_union_file, logger, True)
                # by blast:
                if args.blast == True:
                    transcripts_metrics[i_transcripts].processing_misassembled_psl_file\
                        (alignments_reports[i_transcripts].blast6_report.misassembled_blast6_union_file, logger, False)

            # GET METRICS:
            tot_isoforms_len = None if db_genes_metrics is None else db_genes_metrics.tot_isoforms_len
            transcripts_metrics[i_transcripts].get_transcripts_metrics\
                (args, type_organism, reference_dict, args.transcripts[i_transcripts], transcripts_dicts[i_transcripts],
                 sqlite3_db_genes, tot_isoforms_len, reads_coverage, logger, tmp_dir, WELL_FULLY_COVERAGE_THRESHOLDS,
                 rqconfig.TRANSCRIPT_LENS)

            # GET SEPARATED REPORT:
            separated_reports[i_transcripts].get_separated_report\
                (args, args.labels[i_transcripts], transcripts_dicts[i_transcripts], transcripts_metrics[i_transcripts],
                 db_genes_metrics, reads_coverage, logger, WELL_FULLY_COVERAGE_THRESHOLDS, PRECISION, rqconfig.TRANSCRIPT_LENS)

    # GET COMPARISON REPORT:
    comparison_report = None
    if len(separated_reports) > 1:
        comparison_report = ComparisonReport.ComparisonReport(args.output_dir)
        comparison_report.get_comparison_report(args, transcripts_metrics, db_genes_metrics, reads_coverage, logger,
                                                WELL_FULLY_COVERAGE_THRESHOLDS, PRECISION, rqconfig.TRANSCRIPT_LENS)

    # GET SHORT REPORT:
    short_report = \
        ShortReport.ShortReport(args, db_genes_metrics, transcripts_metrics, args.output_dir, separated_reports,
                                comparison_report, logger, WELL_FULLY_COVERAGE_THRESHOLDS, PRECISION,
                                rqconfig.TRANSCRIPT_LENS)

    # REMOVE TEMPORARY DIRECTORY FROM OUTPUT DIRECTORY:
    if os.path.exists(tmp_dir) and not args.debug:
        logger.debug('Remove temporary directory {}'.format(tmp_dir))
        shutil.rmtree(tmp_dir)
        logger.debug('Done.')

    # LOGGING RESULTS PATHES:
    logger.print_path_results(args, separated_reports, comparison_report, short_report)

    if args.debug:
        UtilsGeneral.profile_memory(args, reference_dict, db_genes_metrics, transcripts_metrics,
                                    separated_reports, comparison_report, logger)

    # FINISH LOGGING:
    logger.finish_up()


if __name__ == '__main__':
    try:
        return_code = main_utils()
        exit(return_code)
    except Exception:
        _, exc_value, _ = sys.exc_info()
        logger.exception(exc_value)
        logger.error('Exception caught!', exit_with_code=1)