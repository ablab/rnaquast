#!/usr/bin/env python

__author__ = 'lenk'

import os
import sys
import shutil

rquast_dirpath = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))

sys.path.insert(0, os.path.join(rquast_dirpath, 'quast_libs'))

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
        UtilsPipeline.run_rnaQUAST_on_test_data(args, rquast_dirpath, program_name, logger)
        # UtilsPipeline.run_rnaQUAST_on_debug_data(args, rquast_dirpath, program_name, logger)
        sys.exit()

    UtilsPipeline.get_abspath_input_data(args)

    # create output directory:
    args.output_dir = UtilsPipeline.create_output_folder(args.output_dir, program_name)
    # create temporary directory:
    tmp_dir = UtilsPipeline.create_empty_folder(os.path.join(args.output_dir, 'tmp'))
    # create directory for log files:
    log_dir = UtilsPipeline.create_empty_folder(os.path.join(args.output_dir, 'logs'))

    # SET LOGGER:
    if args.debug:
        rqconfig.debug = True
        logger.set_up_console_handler(debug=True)
    else:
        logger.set_up_console_handler()
    logger.set_up_file_handler(log_dir)
    logger.print_command_line([os.path.realpath(__file__)] + sys.argv[1:], wrap_after=None)
    logger.start(args.blat, args.busco, args.gene_mark, tmp_dir)

    UtilsPipeline.get_input_data_exist_error(args, logger)

    # THREADING:
    args.threads = UtilsPipeline.get_num_threads(args.threads, logger)

    if args.meta:
        logger.info('\nYOU RUN QUALITY ASSESSMENT FOR METATRANSCRIPTOME ASSEMBLIES')

    # GET segregate FILES:
    if args.reference and args.gtf and len(args.reference) != len(args.gtf):
        logger.error('Numbers of references and gene databases are different', exit_with_code=1)

    args.reference = \
        UtilsPipeline.get_single_file(args.reference, tmp_dir, 'reference', rqconfig.list_ext_fa, args.meta, logger)

    args.gtf = \
        UtilsPipeline.get_single_file(args.gtf, tmp_dir, 'gene_database', rqconfig.list_ext_gtf, args.meta, logger)

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
        # for id_chr in ids_chrs:
        #     reference_dict[id_chr] = UtilsGeneral.correct_nucl_seq(reference_dict[id_chr])


    # for strand specific data we store + and - keys in dictionaries and only + for non strand specific data:
    strands = UtilsGeneral.get_strands(args, logger)

    if args.prokaryote:
        type_organism = 'prokaryotes'
    else:
        type_organism = 'eukaryotes'

    # USE ANNOTATION:
    sqlite3_db_genes = None
    sorted_exons_attr = None
    db_genes_metrics = None
    type_genes, type_isoforms, type_exons = \
        UtilsAnnotations.default_type_genes, \
        UtilsAnnotations.default_type_isoforms, \
        UtilsAnnotations.default_type_exons

    if args.gtf is not None or args.gene_db is not None:
        if args.gene_db is not None:
            gene_db_name = os.path.split(args.gene_db)[1]
            label_db = gene_db_name[:gene_db_name.rfind('.db')]
        else:
            gtf_name = os.path.split(args.gtf)[1]
            label_db = gtf_name[:gtf_name.rfind('.g')]

            if ids_chrs is not None:
                args.gtf = UtilsAnnotations.clear_gtf_by_reference_chr(args.gtf, ids_chrs, tmp_dir, label_db, logger)

        sqlite3_db_genes = \
            UtilsAnnotations.create_sqlite3_db(args.gene_db, args.gtf, label_db,
                                               args.disable_infer_genes, args.disable_infer_transcripts,
                                               args.output_dir, tmp_dir, logger)

        type_genes, type_isoforms, type_exons = \
            UtilsAnnotations.get_type_features(sqlite3_db_genes, UtilsAnnotations.default_type_genes,
                                               UtilsAnnotations.default_type_isoforms,
                                               UtilsAnnotations.default_type_exons, args.prokaryote, logger)

        # if UtilsAnnotations.default_type_exons == type_exons:
        #     type_organism = 'eukaryotes'
        # else:
        #     type_organism = 'prokaryotes'

        db_genes_metrics = GeneDatabaseMetrics.GeneDatabaseMetrics(sqlite3_db_genes, type_genes, type_isoforms, logger)

        ALIGNMENT_THRESHOLDS.ERR_SPACE_TARGET_FAKE_BLAT = db_genes_metrics.max_intron_len + 100
        logger.info('\nSets maximum intron size equal {}. Default is 1500000 bp.\n'.format(ALIGNMENT_THRESHOLDS.ERR_SPACE_TARGET_FAKE_BLAT))

        # set exons starts / ends and ids for binning strategy:
        if ids_chrs is not None:
            sorted_exons_attr = \
                SortedExonsAttributes.SortedExonsAttributes(sqlite3_db_genes, type_exons, strands, ids_chrs, reference_dict, logger)

    reads_coverage = None
    if args.reads_alignment is not None or \
            ((args.single_reads is not None or (args.left_reads is not None and args.right_reads is not None))
             and args.reference is not None and sqlite3_db_genes is not None):
        reads_coverage = \
            ReadsCoverage.ReadsCoverage(args.reads_alignment, args.reference, args.single_reads,
                                        args.left_reads, args.right_reads, reference_dict, sqlite3_db_genes, type_isoforms,
                                        sorted_exons_attr, args.strand_specific, db_genes_metrics.tot_isoforms_len,
                                        genome_len, tmp_dir, args.threads, WELL_FULLY_COVERAGE_THRESHOLDS, logger, log_dir)


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
        logger.warning('No transcripts. Use --transcripts option.')


    # GET PSL ALIGNMENT FILE:
    if args.alignment is None and args.reference is not None and args.transcripts is not None:
        if args.blat:
            args.alignment = UtilsTools.run_blat(None, args.reference, transcripts_dicts, args.labels,
                                                 args.threads, tmp_dir, logger, log_dir)
        else:
            args.alignment = UtilsTools.run_gmap(args.reference, genome_len, args.transcripts, args.labels,
                                                 args.threads, args.gmap_index, tmp_dir, logger, log_dir)

        #if args.fusion_misassemble_analyze:
        #    if not (args.left_reads is not None and args.right_reads is not None):
        #        logger.error('Usage: --left_reads LEFT_READS --right RIGHT_READS for analyse fusions and misassemblies',
        #                     exit_with_code=2, to_stderr=True)
        #        sys.exit(2)


    # FOR MISASSEMBLIES SEARCH:
    # GET DATABASE FOR FA ISOFORMS:
    args.blast = False
    if args.reference is not None and sqlite3_db_genes is not None and args.alignment is not None:
        blastn_run = os.path.join(rqconfig.rnaQUAST_LOCATION, '.', 'blastn')
        if not os.path.isfile(blastn_run):
            blastn_run = "blastn"

        if UtilsGeneral.which(blastn_run) is None:
            logger.warning('blastn not found! Please add blastn to PATH for better MISASSEMBLIES metrics.')
        else:
            args.blast = True

            isoforms_fa_path = os.path.join(tmp_dir, '{}.isoforms.fa'.format(label_db))
            isoforms_list = UtilsGeneral.dict_to_list(UtilsAnnotations.get_fa_isoforms(sqlite3_db_genes, type_isoforms, type_exons, reference_dict, logger))
            fastaparser.write_fasta(isoforms_fa_path, sorted(isoforms_list))

            isoforms_blast_db = UtilsTools.get_blast_db(isoforms_fa_path, label_db, tmp_dir, logger, log_dir)


    # LOGGING INPUT DATA:
    logger.print_input_files(args)


    # INITIALIZATION TRANSCRIPTS METRICS AND REPORTS:
    transcripts_metrics = []
    separated_reports = []
    if args.transcripts is not None:
        alignments_reports = []
        blast_alignments = []
        for i_transcripts in range(len(args.transcripts)):
            # INITIALIZE TRANSCRIPTS METRICS:
            #if args.sam_file is not None:
            #    sam_file_tmp = args.sam_file[i_transcripts]
            #else:
            transcripts_metrics.append(
                TranscriptsMetrics.TranscriptsMetrics(args, args.labels[i_transcripts]))

            # INITIALIZE SEPARATED REPORTS:
            separated_reports.append(SeparatedReport.SeparatedReport(args.labels[i_transcripts], args.output_dir, transcripts_metrics[i_transcripts], WELL_FULLY_COVERAGE_THRESHOLDS))

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

            logger.info()
            logger.info('Processing transcripts from {}:'.format(args.transcripts[i_transcripts]))

            if args.blast:
                blast_alignments.append\
                    (UtilsTools.align_transcripts_to_isoforms_by_blastn
                     (args.transcripts[i_transcripts], isoforms_blast_db, tmp_dir, args.labels[i_transcripts], logger, log_dir))
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
                if args.blast:
                    transcripts_metrics[i_transcripts].processing_misassembled_psl_file\
                        (alignments_reports[i_transcripts].blast6_report.misassembled_blast6_union_file, logger, False)

            # GET METRICS:
            transcripts_metrics[i_transcripts].get_transcripts_metrics\
                (args, type_organism, reference_dict, args.transcripts[i_transcripts], transcripts_dicts[i_transcripts],
                 args.labels[i_transcripts], args.threads, sqlite3_db_genes, db_genes_metrics, reads_coverage, logger,
                 tmp_dir, log_dir, WELL_FULLY_COVERAGE_THRESHOLDS, rqconfig.TRANSCRIPT_LENS)

            # GET SEPARATED REPORT:
            separated_reports[i_transcripts].get_separated_report\
                (args, args.labels[i_transcripts], transcripts_dicts[i_transcripts], transcripts_metrics[i_transcripts],
                 db_genes_metrics, reads_coverage, logger, WELL_FULLY_COVERAGE_THRESHOLDS, PRECISION, rqconfig.TRANSCRIPT_LENS)

    # GET COMPARISON REPORT:
    comparison_report = None
    if len(separated_reports) != 1:
        comparison_report = ComparisonReport.ComparisonReport()
        comparison_report.get_comparison_report(args, args.output_dir, args.labels, transcripts_metrics,
                                                db_genes_metrics, reads_coverage, logger,
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
