__author__ = 'letovesnoi'

import subprocess
import sys
import os
import shutil

this_location = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.join(this_location, 'quast_libs'))

from quast_libs import qutils
from quast_libs import fastaparser

if qutils.is_python2():
    from quast_libs.site_packages.joblib2 import Parallel, delayed
else:
    from quast_libs.site_packages.joblib3 import Parallel, delayed

from general import log
from general import rqconfig


logger = log.get_logger('parallel_blat_run')

def parallel_blat_run(transcripts_dict, reference_pathes, threads, tmp_dir, label, logger, log_dir):
    log_out_1 = os.path.join(log_dir, label + '.blat.out.log')

    logger.print_timestamp()
    logger.info('Getting psl files by BLAT for {}...'.format(label))

    output_dirs = []
    tmp_dirs = []

    us_threads = min(threads, len(transcripts_dict))

    # CREATE TEMPORARY DIRECTORIES FOR THREADS:
    for i_thread in range(us_threads):
        output_dirs.append(os.path.join(tmp_dir, '{}_thread'.format(i_thread)))
        if not os.path.exists(output_dirs[-1]):
            os.mkdir(output_dirs[-1])
        tmp_dirs.append(os.path.join(output_dirs[-1], 'tmp'))
        if not os.path.exists(tmp_dirs[-1]):
            os.mkdir(tmp_dirs[-1])

    # SPLIT CONTIGS TO FILES FOR EACH THREAD:
    transcripts_pathes_threads = split_file_with_transcripts(transcripts_dict, us_threads, output_dirs)

    # PARALLEL RUNS BLAT:
    tmp_names_psl = Parallel(n_jobs=us_threads)(delayed(align_fa_transcripts_to_psl_by_blat)
                                                (transcripts_pathes_threads[i_thread], reference_pathes,
                                                 output_dirs[i_thread], label, log_out_1)
                                                for i_thread in range(us_threads))

    # GLUING PSL FILES:
    f_pathes = ''
    for i_thread in range(us_threads):
        f_pathes += tmp_names_psl[i_thread] + ' '

    logger.print_timestamp()
    logger.info('Gluing psl files by cat for {}...'.format(label))
    out_name_psl = os.path.join(tmp_dir, '{}.psl'.format(label))
    command = 'cat {} >> {}'.format(f_pathes, out_name_psl)
    exit_code = subprocess.call(command, shell=True)
    if exit_code != 0:
        logger.error(message='cat failed!', exit_with_code=exit_code, to_stderr=True)
    logger.info('  saved to {}.'.format(out_name_psl))
    logger.info('  logs can be found in {}.'.format(log_out_1))

    # REMOVE TEMPORARY DIRECTORIES FOR THREADS:
    logger.debug('Remove temporary directories...')
    for i_thread in range(us_threads):
        shutil.rmtree(output_dirs[i_thread])
    logger.debug('Done.')

    return out_name_psl


def split_file_with_transcripts(transcripts_dict, us_threads, output_dirs):
    f_fa_pathes = []

    file_n = len(transcripts_dict)
    thread_n = file_n // us_threads

    id_transcripts = transcripts_dict.keys()[:]
    for i_thread in range(us_threads):
        fpath = os.path.join(output_dirs[i_thread], '{}.fasta'.format(i_thread))
        f_fa_pathes.append(fpath)
        for i_transcript in range(thread_n):
            fastaparser.write_fasta(fpath, [(id_transcripts[i_thread * thread_n + i_transcript], transcripts_dict[id_transcripts[i_thread * thread_n + i_transcript]])], mode='a')

    for i_thread in range(file_n - us_threads * thread_n):
        fastaparser.write_fasta(f_fa_pathes[i_thread], [(id_transcripts[us_threads * thread_n + i_thread], transcripts_dict[id_transcripts[us_threads * thread_n + i_thread]])], mode='a')

    return f_fa_pathes


def align_fa_transcripts_to_psl_by_blat(transcripts_path, reference_pathes, output_dir, label, log_out_1):
    tmp_out_names_psl = []

    # create folder for alignments one file with transcripts to several chromosomes/scaffolds/patches:
    if len(reference_pathes) > 1:
        alignment_dir_i = os.path.join(output_dir, 'tmp', 'alignment_dir_{}'.format(label))
        command = 'mkdir {}'.format(alignment_dir_i)
        exit_code = subprocess.call(command, shell=True)
        if exit_code != 0:
            #logger.error(message='mkdir {} failed!'.format(alignment_dir_i), exit_with_code=2, to_stderr=True)
            sys.exit(exit_code)
    else:
        alignment_dir_i = os.path.join(output_dir, 'tmp')

    for i_reference in range(len(reference_pathes)):

        if len(reference_pathes) > 1:
            tmp_out_names_psl.append(os.path.join(alignment_dir_i, '{}.{}.psl'.format(label, i_reference)))
        else:
            tmp_out_names_psl.append(os.path.join(alignment_dir_i, '{}.psl'.format(label)))

        #logger.print_timestamp()
        #logger.info('Getting psl file by blat for {} and {}...'.format(transcripts_pathes[i_transcripts], reference_pathes[i_reference]))

        blat_run = os.path.join(rqconfig.rnaQUAST_LOCATION, '.', 'blat')
        if not os.path.isfile(blat_run):
            blat_run = "blat"

        command = '{} {} {} {} -q=rna -trimHardA -trimT -noHead 1>> {}'.\
            format(blat_run, reference_pathes[i_reference], transcripts_path, tmp_out_names_psl[-1], log_out_1)
        exit_code = subprocess.call(command, shell=True)
        if exit_code != 0:
            #logger.error(message='blat failed!', exit_with_code=2, to_stderr=True)
            sys.exit(exit_code)

        #logger.info('  saved to {}'.format(tmp_out_names_psl[i_transcripts][-1]))

    if len(reference_pathes) > 1:
        # glue all files with alignments for all chromosomes/scaffolds/patches and one file with transcripts:
        #logger.info('Gluing psl files by pslSort for {}...'.format(alignment_dir_i))

        pslSort_run = os.path.join(rqconfig.rnaQUAST_LOCATION, '.', 'pslSort')

        tmp_alignment_dir_i = os.path.join(alignment_dir_i, 'tmp')

        OUTPSL = os.path.join(output_dir, 'tmp', '{}.psl'.format(label))

        if not os.path.isfile(pslSort_run):
            pslSort_run = "pslSort"
        command = '{} dirs {} {} {} -nohead 1>> {}'.\
            format(pslSort_run, OUTPSL, tmp_alignment_dir_i, alignment_dir_i, log_out_1)
        exit_code = subprocess.call(command, shell=True)
        if exit_code != 0:
            #logger.error(message='pslSort failed! File with reference more than 4G, please install pslSort '
            #                     '(http://hgwdev.cse.ucsc.edu/~kent/exe/linux/), add them to PATH and restart', exit_with_code=2, to_stderr=True)
            sys.exit(exit_code)

        # remove temporary directory with separately transcripts to chromosomes/scaffolds/patches alignments:
        if os.path.exists(alignment_dir_i):
            #logger.debug('Remove temporary directory {}'.format(alignment_dir_i))
            shutil.rmtree(alignment_dir_i)
            #logger.debug('Done.')
    else:
        OUTPSL = tmp_out_names_psl[0]

    return OUTPSL
