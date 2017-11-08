__author__ = 'letovesnoi'

import os
import sys
import subprocess
import shutil

import math

import datetime

import rqconfig

import UtilsGeneral
import UtilsPipeline

import parallel_blat_run


def run_blat(args_database, args_reference, transcripts_dicts, args_labels, args_threads, tmp_dir, logger, log_dir):

    blat_run = os.path.join(rqconfig.rnaOUAST_LOCATION, '.', 'blat')
    if not os.path.isfile(blat_run):
        blat_run = "blat"

    if UtilsGeneral.which(blat_run) is None:
        logger.error('BLAT not found! Please add BLAT to PATH for ALIGNMENT metrics.')
    else:
        # for single file with scaffolds/patches/chromosomes:
        # split big single file with reference to files with scaffolds/patches/chromosomes
        # this save us from blat segfault for reference more than 4g:
        if args_database is None and os.path.getsize(args_reference) >= 4294967296:
            args_database = get_database_split_chr(tmp_dir, args_reference, logger)

        if os.path.getsize(args_reference) < 4294967296:
            args_reference = UtilsGeneral.get_upper_case_fasta(args_reference, tmp_dir, logger)
            reference_pathes = [args_reference]
        else:
            # get upper case (don't mask repeats):
            args_database = get_upper_case_database_split_chr(args_database, tmp_dir, logger)
            reference_pathes = args_database

        # RUN BLAT:
        args_alignment = []
        for i_transcripts in range(len(transcripts_dicts)):
            start_time = datetime.datetime.now()

            args_alignment.append(parallel_blat_run.parallel_blat_run(
                transcripts_dicts[i_transcripts], reference_pathes, args_threads, tmp_dir,
                args_labels[i_transcripts], logger, log_dir))

            end_time = datetime.datetime.now()
            spent_time = end_time - start_time
            logger.info('\nBLAT TIME: {}\n\n'.format(spent_time))

            # args.alignment.append(UtilsPipeline.align_fa_transcripts_to_psl_by_blat
            #                       (args.transcripts[i_transcripts], reference_pathes, args.output_dir,
            #                        args.labels[i_transcripts]))

        return args_alignment


def get_database_split_chr(output_dir, reference_path, logger):
    database_dir = os.path.join(output_dir, 'database_dir')
    command = 'mkdir {}'.format(database_dir)
    subprocess.call(command, shell=True)

    # path to file with pathes to scaffolds/chromosomes/patches:
    chrs_database_path = os.path.join(output_dir, 'scaffolds.database')

    # list of pathes to scaffolds/chromosomes/patches:
    database_pathes = []

    logger.print_timestamp()
    logger.info('Getting split scaffolds database...')

    fout1 = open(chrs_database_path, 'w')
    with open(reference_path, 'r') as fin:
        file_chr_name = None
        fout2 = None
        for line in fin:
            if line[0] == '>':
                file_chr_name = line.strip().split(' ')[0][1:] + '.fa'
                if fout2 is not None:
                    fout2.close()

                chrs_file_path = os.path.join(database_dir, file_chr_name)
                fout2 = open(chrs_file_path, 'a')

                fout1.write(chrs_file_path + '\n')
                database_pathes.append(chrs_file_path)

            fout2.write(line.strip() + '\n')
        fout1.close()

    logger.info('  saved to {}'.format(chrs_database_path))

    return database_pathes


def get_upper_case_database_split_chr(database, tmp_dir, logger):
    database_dir = os.path.join(tmp_dir, 'database_upper_dir')
    command = 'mkdir {}'.format(database_dir)
    subprocess.call(command, shell=True)

    # path to file with pathes to scaffolds/chromosomes/patches:
    chrs_database_path = os.path.join(tmp_dir, 'scaffolds.upper.database')

    # list of pathes to scaffolds/chromosomes/patches:
    reference_pathes = []

    logger.print_timestamp()
    logger.info('Getting upper case split scaffolds database...')

    fout0 = open(chrs_database_path, 'w')

    for in_name_fa in database:
        (tmp_dir_name, tmp_file_name) = os.path.split(in_name_fa)

        out_name_fa = os.path.join(database_dir, tmp_file_name)
        reference_pathes.append(out_name_fa)

        fin1 = open(in_name_fa, 'r')
        fout1 = open(out_name_fa, 'w')
        fout0.write(out_name_fa + '\n')
        for line1 in fin1:
            if line1[0] != '>':
                line1 = line1.upper()
            fout1.write(line1)
        fin1.close()
        fout1.close()

    fout0.close()

    logger.info('  saved to {}'.format(chrs_database_path))

    return reference_pathes


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

        blat_run = os.path.join(rqconfig.rnaOUAST_LOCATION, '.', 'blat')
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

        pslSort_run = os.path.join(rqconfig.rnaOUAST_LOCATION, '.', 'pslSort')

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


def get_blast_db(isoforms_fa_path, gtf_label, tmp_dir, logger, log_dir):
    program_name = 'makeblastdb'

    log_out = os.path.join(log_dir, program_name + '.log')

    logger.print_timestamp()
    logger.info('Getting blast database for {}'.format(isoforms_fa_path))

    isoforms_blast_db = os.path.join(tmp_dir, '{}.isoforms'.format(gtf_label))

    command = '{} -in {} -dbtype nucl -out {} >> {}'.format(program_name, isoforms_fa_path, isoforms_blast_db, log_out)
    exit_code = subprocess.call(command, shell=True)
    if exit_code != 0:
        logger.error(message='{} failed!'.format(program_name), exit_with_code=exit_code, to_stderr=True)
        sys.exit(exit_code)

    logger.info('  saved to {}'.format(isoforms_blast_db))

    return isoforms_blast_db


def align_transcripts_to_isoforms_by_blastn(transcripts_path, isoforms_blast_db, tmp_dir, label, logger, log_dir):
    program_name = 'blastn'

    log_out = os.path.join(log_dir, label + '.' + program_name + '.log')

    logger.print_timestamp('  ')
    logger.info('  Aligning {} to {} by blastn...'.format(transcripts_path, isoforms_blast_db))

    alignment_isoforms_path = '{}.blast6'.format(os.path.join(tmp_dir, label))

    command = '{} -query {} -out {} -db {} -num_alignments 10 -evalue 0.01 -outfmt "6 qseqid sseqid pident length ' \
              'mismatch gapopen qstart qend sstart send evalue bitscore sstrand" 1>> {}'.\
        format(program_name, transcripts_path, alignment_isoforms_path, isoforms_blast_db, log_out)
    exit_code = subprocess.call(command, shell=True)
    if exit_code != 0:
        logger.error(message='{} failed for {}!'.format(program_name, label), exit_with_code=exit_code, to_stderr=True)
        sys.exit(exit_code)

    logger.info('    saved to {}'.format(alignment_isoforms_path))

    logger.info('    log can be found in {}.'.format(log_out))

    return alignment_isoforms_path


def run_gmap(args_reference, genome_len, args_transcripts, args_labels, args_threads, args_gmap_index, tmp_dir, logger, log_dir):
    args_alignment = []

    if genome_len < 2 ** 32:
        gmap_run = 'gmap'
    else:
        gmap_run = 'gmapl'

    gmap_build = 'gmap_build'

    gmap_build_logger_out_path = os.path.join(log_dir, gmap_build + '.out.log')
    gmap_build_logger_err_path = os.path.join(log_dir, gmap_build + '.err.log')

    args_reference = UtilsGeneral.get_upper_case_fasta(args_reference, tmp_dir, logger)

    ref_label = os.path.split(args_reference)[-1][:os.path.split(args_reference)[-1].rfind('.f')]

    # if UtilsGeneral.which(gmap_run) is None or UtilsGeneral.which(gmap_build):
    #     logger.warning('gmap or gmap_build not found! Please add GMAP to PATH or run with BLAT for ALIGNMENT metrics.')
    # else:
    # RUN GMAP:
    # create index (gmap_build):
    if args_gmap_index is None:
        logger.print_timestamp()
        logger.info('Creating genome index by {}...'.format(gmap_build))

        start_time = datetime.datetime.now()

        command = '{gmap_build} -D {tmp_dir} -d {ref_index_name} {reference} 1>> {log_out_1} 2>> {log_out_2}'.\
            format(gmap_build=gmap_build, tmp_dir=tmp_dir, ref_index_name=ref_label, reference=args_reference,
                   log_out_1=gmap_build_logger_out_path, log_out_2=gmap_build_logger_err_path)
        exit_code = subprocess.call(command, shell=True)

        logger.info('  logs can be found in {} and {}.'.format(gmap_build_logger_out_path, gmap_build_logger_err_path))

        if exit_code != 0:
            logger.error(message='{} failed!'.format(gmap_build), exit_with_code=exit_code, to_stderr=True)

        end_time = datetime.datetime.now()
        spent_time = end_time - start_time

        logger.info('  saved to {}'.format(os.path.join(tmp_dir, ref_label)))

        logger.info('\nGMAP_BUILD TIME: {}\n\n'.format(spent_time))
    else:
        command = 'ln -s {} {}'.format(args_gmap_index, os.path.join(tmp_dir, ref_label))
        subprocess.call(command, shell=True)

    # align (gmap):
    for i_transcripts in range(len(args_transcripts)):
        gmap_run_logger_err_path = os.path.join(log_dir, gmap_run + '.' + args_labels[i_transcripts] + '.err.log')

        logger.print_timestamp()
        logger.info('Aligning {} to {}...'.format(args_labels[i_transcripts], ref_label))

        alignment_psl_path = os.path.join(tmp_dir, args_labels[i_transcripts] + '.psl')

        start_time = datetime.datetime.now()

        command = '{gmap} -D {tmp_dir} -d {ref_index_name} {transcripts} --format=1 -t {threads} -O > {alignment_out}' \
                  ' 2>> {log_out_2}'.\
            format(gmap=gmap_run, tmp_dir=tmp_dir, ref_index_name=ref_label, transcripts=args_transcripts[i_transcripts],
                   threads=args_threads, alignment_out=alignment_psl_path, log_out_2=gmap_run_logger_err_path)
        exit_code = subprocess.call(command, shell=True)

        logger.info('  log can be found in {}.'.format(gmap_run_logger_err_path))

        if exit_code != 0:
            logger.error(message='{} failed for {}!'.format(gmap_run, args_labels[i_transcripts]),
                         exit_with_code=exit_code, to_stderr=True)

        end_time = datetime.datetime.now()
        spent_time = end_time - start_time

        args_alignment.append(alignment_psl_path)

        logger.info('  saved to {}'.format(alignment_psl_path))

        logger.info('\nGMAP TIME: {}\n\n'.format(spent_time))

    return args_alignment


# https://github.com/alexdobin/STAR/releases
# It is strongly recommended to include major chromosomes as well as un-placed and un-localized scaffolds.
def run_STAR(threads, reference_path, gtf_path, single_reads, left_reads, right_reads, output_dir,
             sjdbGTFtagExonParentTranscript, sjdbGTFtagExonParentGene, genome_len, logger, log_dir):
    # Basic STAR workflow consists of 2 steps:
    program_name = 'STAR'

    star_logger_out_path = os.path.join(log_dir, program_name + '.out.log')
    star_logger_err_path = os.path.join(log_dir, program_name + '.err.log')

    logger.print_timestamp()
    logger.info('Running {}...'.format(program_name))

    # create STAR output directory:
    star_outdir = UtilsPipeline.create_folder(os.path.join(output_dir, 'star_out'))
    # out_sorted_bam_path = os.path.join(star_outdir, 'Aligned.sortedByCoord.out.bam')

    # 1 Generating genome indexes files (supplied the reference genome sequences (FASTA files)
    # and annotations (GTF file))
    genome_dir = os.path.join(star_outdir, 'genome_dir')

    if not os.path.exists(genome_dir):
        mode = '--runMode'

        # create tmp output directory:
        tmp_dir = UtilsPipeline.create_empty_folder(os.path.join(star_outdir, 'tmp_dir'))

        # create tmp_genome_dir directory:
        tmp_genome_dir = UtilsPipeline.create_empty_folder(os.path.join(tmp_dir, 'genome_dir'))

        genomeSAindexNbases = min(14, math.log(genome_len, 2) / 2 - 1)

        command = '{program_name} {mode} genomeGenerate --runThreadN {threads} --genomeDir {tmp_genome_dir} ' \
                  '--genomeFastaFiles {reference} --genomeSAindexNbases {genomeSAindexNbases}'.\
            format(program_name=program_name, mode=mode, threads=threads, tmp_genome_dir=tmp_genome_dir,
                   reference=reference_path, genomeSAindexNbases=genomeSAindexNbases)

        if gtf_path is not None:
            command += ' --sjdbGTFfile {gtf} --sjdbGTFtagExonParentTranscript {parent_transcript} --sjdbGTFtagExonParentGene {parent_gene}'.\
                format(gtf=gtf_path, parent_transcript=sjdbGTFtagExonParentTranscript, parent_gene=sjdbGTFtagExonParentGene)
        command += ' 1>> {log_out_1} 2>> {log_out_2}'.format(log_out_1=star_logger_out_path, log_out_2=star_logger_err_path)

        logger.print_timestamp()
        logger.info('  ' + command)

        exit_code = subprocess.call(command, shell=True)

        logger.info('    logs can be found in {} and {}.'.format(star_logger_out_path, star_logger_err_path))

        if exit_code != 0:
            logger.error('{program_name_mode} failed!'. format(program_name_mode=program_name + ' ' + mode,
                                                               program_name=program_name))
        else:
            command = 'mv {} {}'.format(tmp_genome_dir, star_outdir)
            subprocess.call(command, shell=True)

    # 2 Mapping reads to the genome (supplied the genome files generated in the 1st step, as well as the RNA-seq reads
    # (sequences) in the form of FASTA or FASTQ files.)
    readFilesIn = ''
    if single_reads is not None:
        readFilesIn += single_reads + ' '
    if right_reads is not None and left_reads is not None:
        readFilesIn += left_reads + ' ' + right_reads
    command = '{program_name} --runThreadN {threads} --genomeDir {genome_dir} --readFilesIn {readFilesIn} ' \
              '--outFileNamePrefix {out_file_name_prefix} --outSAMtype SAM ' \
              '--limitBAMsortRAM 1000706316'.\
        format(program_name=program_name, threads=threads, genome_dir=genome_dir, readFilesIn=readFilesIn,
               out_file_name_prefix=star_outdir + '/', log_out_1=star_logger_out_path, log_out_2=star_logger_err_path)

    # for compressed read files:
    if '.gz' in single_reads and '.gz' in left_reads and '.gz' in right_reads:
        command += ' --readFilesCommand zcat'
    # if '.bz2' in single_reads and '.bz2' in left_reads and '.bz2' in right_reads:
    #     command += ' --readFilesCommand bzcat'

    command += ' 1>> {log_out_1} 2>> {log_out_2}'

    logger.print_timestamp()
    logger.info('  ' + command)
    exit_code = subprocess.call(command, shell=True)
    if exit_code != 0:
        star_outdir = None

        logger.error('{program_name} failed!'.format(program_name=program_name))
    else:
        logger.info('  saved to {}.'.format(star_outdir))

    logger.info('  logs can be found in {} and {}.'.format(star_logger_out_path, star_logger_err_path))

    return star_outdir


def get_sam_by_STAR(threads, reference_path, gtf_path, single_reads, left_reads, right_reads, output_dir,
                           sjdbGTFtagExonParentTranscript, sjdbGTFtagExonParentGene, genome_len, logger, log_dir):
    star_outdir = run_STAR(threads, reference_path, gtf_path, single_reads, left_reads, right_reads, output_dir,
                           sjdbGTFtagExonParentTranscript, sjdbGTFtagExonParentGene, genome_len, logger, log_dir)

    if star_outdir is not None:
        out_sam_path = os.path.join(star_outdir, 'Aligned.out.sam')
    else:
        out_sam_path = None

    return out_sam_path


# https://ccb.jhu.edu/software/tophat/manual.shtml
def run_tophat(bowtie2_index_path, reference_path, single_reads, reads_1_path, reads_2_path, output_dir, threads,
               logger, log_dir):
    program_name = 'tophat'

    # We are changing the FASTA file ending from .fna to .fa,
    # because tophat wants a file names .fa, its picky that way
    if os.path.splitext(reference_path)[1] != '.fa':
        new_ref_path = os.path.join(output_dir, os.path.basename(reference_path)[:os.path.basename(reference_path).rfind('.f')] + '.fa')
        command = 'ln -s {} {}'.format(reference_path, new_ref_path)
        subprocess.call(command, shell=True)
        reference_path = new_ref_path

    tophat_logger_err_path = os.path.join(log_dir, program_name + '.err.log')

    tophat_outdir = UtilsPipeline.create_folder(os.path.join(output_dir, program_name + '_out'))

    if bowtie2_index_path is None:
        bowtie2_index_path = get_genome_bowtie2_index(reference_path, logger, log_dir)

    logger.print_timestamp()
    logger.info('Running {}...'.format(program_name))

    reads = ''
    if reads_1_path is not None and reads_2_path is not None:
        reads += reads_1_path + ' ' + reads_2_path
    if single_reads is not None:
        if reads_1_path is not None and reads_2_path is not None:
            reads += ','
        reads += single_reads

    command = \
        '{program_name} -o {output_dir} {index} {reads} -p {threads} 2>> {log_out_2}'.\
            format(program_name=program_name, output_dir=tophat_outdir, index=bowtie2_index_path, reads=reads,
                   threads=threads, log_out_2=tophat_logger_err_path)
    exit_code = subprocess.call(command, shell=True)
    if exit_code != 0:
        tophat_outdir = None

        logger.error('{program_name} failed!'. format(program_name=program_name))
    else:
        logger.info('  saved to {}.'.format(tophat_outdir))

    logger.info('  log can be found in {}.'.format(tophat_logger_err_path))

    return tophat_outdir


def get_sam_by_tophat(bowtie2_index_path, reference_path, single_reads, reads_1_path, reads_2_path, output_dir,
                      threads, logger, log_dir):
    tophat_outdir = run_tophat(bowtie2_index_path, reference_path, single_reads, reads_1_path, reads_2_path,
                               output_dir, threads, logger, log_dir)

    if tophat_outdir is not None:
        out_bam_path = os.path.join(tophat_outdir, 'accepted_hits.bam')

        out_sam_path = bam2sam(out_bam_path, tophat_outdir, logger)
    else:
        out_sam_path = None

    # out_sorted_bam_path = get_sort_bam(out_bam_path, tophat_outdir, logger, type='position')

    return out_sam_path


def bam2sam(in_bam_path, output_dir, logger):
    program_name = 'samtools view'

    logger.print_timestamp()
    logger.info('Running {}...'.format(program_name))

    in_bam_name = os.path.split(in_bam_path)[-1]
    out_sam_path = os.path.join(output_dir, in_bam_name[:in_bam_name.rfind('.bam')] + '.sam')

    command = '{program_name} -h -o {sam} {bam}'.format(program_name=program_name, sam=out_sam_path, bam=in_bam_path)
    exit_code = subprocess.call(command, shell=True)
    if exit_code != 0:
        logger.error('{program_name} failed! Please add {program_name} in your PATH.'.format(program_name=program_name))
        sys.exit(2)

    logger.info('  saved to {}.'.format(out_sam_path))

    return out_sam_path


def chg_ref_names_for_tophat_GTF(in_ref_path, out_ref_path, logger):
    logger.info('Modify the names of the reference sequences for exactly matching GTF/GFF column which indicates '
                'the chromosome or contig on which the feature is located')

    fin = open(in_ref_path, 'r')

    fout = open(out_ref_path, 'w')

    for in_line in fin:
        if in_line[0] == '>':
            out_line = in_line.strip().split()[0] + '\n'
            logger.debug(out_line)
        else:
            out_line = in_line

        fout.write(out_line)

    fin.close()

    fout.close()

    logger.info('  saved to {}'.format(out_ref_path))

    return out_ref_path


# The default is name.
def get_sort_bam(in_bam_path, output_dir, logger, type='name'):
    in_bam_file_name = os.path.split(in_bam_path)[-1]
    out_bam_path = os.path.join(output_dir, in_bam_file_name[:in_bam_file_name.rfind('.')] + '.sorted.bam')

    program_name = 'samtools sort'

    logger.print_timestamp()
    logger.info('Sorting {in_bam} by {program_name}'.format(in_bam=in_bam_path, program_name=program_name))

    if type == 'name':
        command = '{program_name} -no {in_bam} {out_bam}'.format(program_name=program_name, in_bam=in_bam_path,
                                                                 out_bam=out_bam_path)
    else:
        command = '{program_name} {in_bam} {out_bam}'.format(program_name=program_name, in_bam=in_bam_path,
                                                             out_bam=out_bam_path[:-4])

    exit_code = subprocess.call(command, shell=True)
    if exit_code != 0:
        logger.error('{program_name} failed! Please add {program_name} in your PATH.'.format(program_name=program_name),
                     exit_with_code=exit_code)

    logger.info('  saved to {}.'.format(out_bam_path))

    return out_bam_path


# Please note that it is highly recommended that a FASTA file with the sequence(s) the genome being indexed
# be present in the same directory with the Bowtie index files and having the name <genome_index_base>.fa.
def get_genome_bowtie2_index(reference_path, logger, log_dir):
    program_name = 'bowtie2-build'

    bowtie_logger_out_path = os.path.join(log_dir, program_name + '.out.log')

    bowtie_logger_err_path = os.path.join(log_dir, program_name + '.err.log')

    logger.print_timestamp()
    logger.info('Indexing {reference} by {program_name}...'.format(reference=reference_path, program_name=program_name))

    out_bowtie2_index_path = reference_path[:reference_path.rfind('.f')]

    command = '{program_name} {reference} {index} 1>> {log_out_1} 2>> {log_out_2}'.\
        format(program_name=program_name, reference=reference_path, index=out_bowtie2_index_path,
               log_out_1=bowtie_logger_out_path, log_out_2=bowtie_logger_err_path)
    exit_code = subprocess.call(command, shell=True)
    if exit_code != 0:
        out_bowtie2_index_path = None

        logger.error('{program_name} failed!'.format(program_name=program_name))
    else:
        logger.info('  saved to {}.*.'.format(out_bowtie2_index_path))

    logger.info('  logs can be found in {} and {}.'.format(bowtie_logger_out_path, bowtie_logger_err_path))

    return out_bowtie2_index_path


# def simulate_fq_reads_by_flux(args, logger):
#     out_name_fq = os.path.join(args.output_dir, 'tmp', 'reads.flux.fq')
#     logger.print_timestamp()
#     logger.info('Getting fq file with reads by flux...')
#
#     command = 'flux-simulator -p {}'.format(args.par)
#     exit_code = subprocess.call(command, shell=True)
#     if exit_code != 0:
#         logger.error(message='Flux failed!', exit_with_code=exit_code, to_stderr=True)
#         sys.exit(exit_code)
#     command = 'mv {}.fastq {}'.format(args.par[:args.par.rfind('.')], out_name_fq)
#     subprocess.call(command, shell=True)
#
#     logger.info('  saved to {}'.format(out_name_fq))
#
#     return out_name_fq


# def assemble_fq_to_fa_by_trinity(args, logger):
#     command = None
#     out_dir_fa = os.path.join(args.output_dir, 'tmp')
#     out_name_fa = os.path.join(out_dir_fa, 'Trinity.fasta')
#     logger.print_timestamp()
#     logger.info('Getting fa file with transcripts by Trinity...')
#     if args.left_reads != None and args.right_reads != None:
#         command = 'Trinity --seqType fq --JM 10G --left {} --right {} --output {} --CPU {}'.\
#             format(args.left_reads, args.right_reads, out_dir_fa, args.threads)
#     elif args.paired_reads != None:
#         command = 'Trinity --seqType fq --JM 10G --single {} --run_as_paired --output {} --CPU {}'.format(args.paired_reads, out_dir_fa, args.threads)
#     elif args.single_reads != None:
#         command = 'Trinity --seqType fq --JM 10G --single {} --output {} --CPU {}'.format(args.single_reads, out_dir_fa, args.threads)
#     exit_code = subprocess.call(command, shell=True)
#     if exit_code != 0:
#         logger.error(message='Trinity failed!', exit_with_code=exit_code, to_stderr=True)
#         sys.exit(exit_code)
#
#     logger.info('  saved to {}'.format(out_name_fa))
#
#     return out_name_fa


# def assemble_fq_to_fa_by_spades(args, logger):
#     command = None
#     out_dir_fa = os.path.join(args.output_dir, 'tmp')
#     out_name_fa = os.path.join(out_dir_fa, 'contigs.fasta')
#     logger.print_timestamp()
#     logger.info('Getting fa file with transcripts by SPAdes...')
#     if args.left_reads != None and args.right_reads != None:
#         command = 'spades.py --sc -1 {} -2 {} -o {} --threads {}'.format(args.left_reads, args.right_reads, out_dir_fa, args.threads)
#     elif args.paired_reads != None:
#         command = 'spades.py --sc --12 {} -o {} --threads {}'.format(args.paired_reads, out_dir_fa, args.threads)
#     elif args.single_reads != None:
#         command = 'spades.py --sc -s {} -o {} --threads {}'.format(args.single_reads, out_dir_fa, args.threads)
#     exit_code = subprocess.call(command, shell=True)
#     if exit_code != 0:
#         logger.error(message='SPAdes failed!', exit_with_code=exit_code, to_stderr=True)
#         sys.exit(exit_code)
#
#     logger.info('  saved to {}'.format(out_name_fa))
#
#     return out_name_fa


# def assemble_fq_reads_to_fa_transcripts(args, logger):
#     out_names_fa = []
#     for assembler in args.assembler:
#         if assembler == 'Trinity':
#             out_names_fa.append(assemble_fq_to_fa_by_trinity(args, logger))
#         elif assembler == 'SPAdes':
#             out_names_fa.append(assemble_fq_to_fa_by_spades(args, logger))
#     return out_names_fa


# def get_file_with_misassemblies_by_reads(args, logger, transcripts_metrics, transcripts_file, in_sam_file):
#     tool_for_mis_by_reads_name = 'misfinder.py'
#     tool_for_mis_by_reads_path = os.path.join(rqconfig.rnaOUAST_LOCATION, tool_for_mis_by_reads_name)
#
#     mis_by_reads_file = os.path.join(args.output_dir, 'tmp', '{}.mis_by_reads'.format(transcripts_metrics.label))
#
#     logger.print_timestamp()
#     logger.info('Getting file with misassemblies for {} by misfinder...'.format(transcripts_file))
#
#     if in_sam_file != None:
#         command = '{} -o {} -c {} -r1 {} -r2 {} -s {}'.format(tool_for_mis_by_reads_path, os.path.join(args.output_dir, 'tmp'), transcripts_file, args.left_reads, args.right_reads, in_sam_file)
#     else:
#         command = '{} -o {} -c {} -r1 {} -r2 {}'.format(tool_for_mis_by_reads_path, os.path.join(args.output_dir, 'tmp'), transcripts_file, args.left_reads, args.right_reads)
#
#     exit_code = subprocess.call(command, shell=True)
#     if exit_code != 0:
#         logger.error(message='{} failed!'.format(tool_for_mis_by_reads_name), exit_with_code=2, to_stderr=True)
#         sys.exit(exit_code)
#
#     path_results_mis = os.path.join(args.output_dir, 'tmp', 'result.txt')
#     command = 'mv {} {}'.format(path_results_mis, mis_by_reads_file)
#     subprocess.call(command, shell=True)
#
#     if in_sam_file == None:
#         path_results_sam = os.path.join(args.output_dir, 'tmp', 'aligned.sam')
#         out_sam_file = os.path.join(args.output_dir, 'tmp', '{}.sam'.format(transcripts_metrics.label))
#         command = 'mv {} {}'.format(path_results_sam, out_sam_file)
#         subprocess.call(command, shell=True)
#
#     logger.info('  saved to {}'.format(mis_by_reads_file))
#
#     return mis_by_reads_file


# def get_isoforms_fa_from_gtf(annotation_label, reference_path, annotation_path, tmp_dir, logger):
#     logger.info('Getting fasta with isoforms for {}...'.format(annotation_path))
#
#     isoforms_fa_path = os.path.join(tmp_dir, '{}.isoforms.fa'.format(annotation_label))
#
    # create links for reference and annotation: ngsutils create temporary files in them folders:
    # reference_link = os.path.join(tmp_dir, os.path.split(reference_path)[1])
    # if not os.path.exists(reference_link):
    #     command = 'ln -s {} {}'.format(reference_path, tmp_dir)
    #     subprocess.call(command, shell=True)
    #
    # annotation_link = os.path.join(tmp_dir, os.path.split(annotation_path)[1])
    # if not os.path.exists(annotation_link):
    #     command = 'ln -s {} {}'.format(annotation_path, tmp_dir)
    #     subprocess.call(command, shell=True)
    #
    # isoforms_bed_path = os.path.join(tmp_dir, '{}.isoforms.bed'.format(annotation_label))
    #
    # create bed file with account of alternative splicing:
    # command = 'gtfutils tobed -regions {} > {}'.format(annotation_link, isoforms_bed_path)
    # exit_code = subprocess.call(command, shell=True)
    # if exit_code != 0:
    #     logger.error(message='gtfutils tobed failed! Please install and add to PATH gtfutils: Tools for next-generation sequencing analysis '
    #                          '(http://ngsutils.org/)', exit_with_code=2, to_stderr=True)
    #     sys.exit(2)
    #
    # create fa file from bed:
    # command = 'bedutils tofasta -name {} {} > {}'.format(isoforms_bed_path, reference_link, isoforms_fa_path)
    # command = 'bedtools getfasta -fi {} -bed {} -fo {} -s -name'.format(reference_link, annotation_path, isoforms_fa_path)
    # exit_code = subprocess.call(command, shell=True)
    # if exit_code != 0:
    #     logger.error(message='bedtools getfasta failed! Please install and add to PATH bedtools: a powerful toolset '
    #                          'for genome arithmetic (http://bedtools.readthedocs.org/en/latest/)', exit_with_code=2, to_stderr=True)
    #     sys.exit(2)
    #
    # logger.info('  saved to {}'.format(isoforms_fa_path))
    #
    # return isoforms_fa_path