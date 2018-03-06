__author__ = 'letovesnoi'

import subprocess
import sys
import os
import shutil

from joblib import Parallel, delayed

this_location = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.join(this_location, 'quast_libs'))

from quast_libs import fastaparser

from general import log
import UtilsTools


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
    tmp_names_psl = Parallel(n_jobs=us_threads)(delayed(UtilsTools.align_fa_transcripts_to_psl_by_blat)
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
    thread_n = file_n / us_threads

    id_transcripts = transcripts_dict.keys()[:]
    for i_thread in range(us_threads):
        fpath = os.path.join(output_dirs[i_thread], '{}.fasta'.format(i_thread))
        f_fa_pathes.append(fpath)
        for i_transcript in range(thread_n):
            fastaparser.write_fasta(fpath, [(id_transcripts[i_thread * thread_n + i_transcript], transcripts_dict[id_transcripts[i_thread * thread_n + i_transcript]])], mode='a')

    for i_thread in range(file_n - us_threads * thread_n):
        fastaparser.write_fasta(f_fa_pathes[i_thread], [(id_transcripts[us_threads * thread_n + i_thread], transcripts_dict[id_transcripts[us_threads * thread_n + i_thread]])], mode='a')

    return f_fa_pathes