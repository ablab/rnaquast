__author__ = 'lenk'

import os
import sys
import subprocess

import gffutils

import UtilsGeneral

default_type_genes = ['gene', 'miRNA_gene']

default_type_isoforms = ['transcript', 'RNA', 'mRNA', 'miRNA', 'ncRNA', 'tRNA', "rRNA", "tmRNA"]

default_type_exons = ['exon', 'coding_exon', 'noncoding_exon']


# general transform (included all transforms):
def transform(feature):
    feature = transform_fancy_id(feature)

    feature = transform_identical_gene_transcript_id(feature)

    # feature = transform_appropriate_parent(feature)

    return feature


def transform_fancy_id(feature):
    if feature.featuretype not in default_type_genes + default_type_isoforms:
        exon_location = '{}:{}:{}-{}:{}'.format(feature.featuretype, feature.seqid, feature.start, feature.stop, feature.strand)
        feature_id = exon_location
        if feature.featuretype == 'CDS':
            feature_id += ':' + feature.frame
        feature.attributes['fancy_id'] = [feature_id]
    return feature


def transform_identical_gene_transcript_id(feature):
    if 'transcript_id' in feature.attributes:
        for i_transcript in range(len(feature.attributes['transcript_id'])):
            feature.attributes['transcript_id'][i_transcript] += '_transcript'
    if 'gene_id' in feature.attributes:
        for i_gene in range(len(feature.attributes['gene_id'])):
            feature.attributes['gene_id'][i_gene] += '_gene'
    return feature


# set 1-level parent relation:
def transform_appropriate_parent(feature):
    if feature.featuretype in default_type_isoforms and 'gene_id' in feature.attributes:
        feature.attributes['Parent'] = feature.attributes['gene_id']
    elif feature.featuretype in default_type_exons and 'transcript_id' in feature.attributes:
        feature.attributes['Parent'] = feature.attributes['transcript_id']
    return feature


def get_id_spec():
    id_spec = {}
    for type in default_type_genes:
        id_spec[type] = ['ID', 'gene_id']

    for type in default_type_isoforms:
        id_spec[type] = ['ID', 'transcript_id']

    for type in default_type_exons:
        id_spec[type] = 'fancy_id'

    id_spec['CDS'] = 'fancy_id'

    id_spec['UTR'] = 'fancy_id'

    id_spec['start_codon'] = 'fancy_id'

    id_spec['stop_codon'] = 'fancy_id'

    return id_spec


def child_func(parent, child):
    if 'Parent' not in child.attributes:
        child.attributes['Parent'] = []
    child.attributes['Parent'].append(parent.id)


def create_sqlite2_db_in_file(in_gff_path, disable_infer_genes, disable_infer_transcripts,
                              sqlite3_db_path, tmp_sqlite3_db_path, logger):

    logger.print_timestamp()
    logger.info('Creating sqlite3 db by gffutils...')

    id_spec = get_id_spec()

    gffutils.create_db(in_gff_path, dbfn=tmp_sqlite3_db_path, force=True, verbose=True,
                                  merge_strategy="merge", transform=transform, force_merge_fields=['source'],
                                  id_spec=id_spec, disable_infer_genes=disable_infer_genes,
                                  disable_infer_transcripts=disable_infer_transcripts)

    command = 'mv {} {}'.format(tmp_sqlite3_db_path, sqlite3_db_path)
    subprocess.call(command, shell=True)

    logger.info('  saved to {}.'.format(sqlite3_db_path))

    return sqlite3_db_path


def load_sqlite3_db(sqlite3_db_path, logger):
    logger.print_timestamp()
    logger.info('Loading sqlite3 db by gffutils from {} to memory...'.format(sqlite3_db_path))

    sqlite3_db = gffutils.FeatureDB(sqlite3_db_path)

    logger.info('Done.')

    return sqlite3_db


# create database for gff/gtf file:
def create_sqlite3_db(in_gene_db, in_gff_path, label_db, disable_infer_genes, disable_infer_transcripts,
                      output_dir, tmp_dir, logger, prokaryote):
    tmp_sqlite3_db_path = os.path.join(tmp_dir, label_db + '.db')
    sqlite3_db_path = os.path.join(output_dir, label_db + '.db')

    if in_gene_db is not None:
        sqlite3_db_genes = load_sqlite3_db(in_gene_db, logger)
    elif os.path.exists(sqlite3_db_path):
        sqlite3_db_genes = load_sqlite3_db(sqlite3_db_path, logger)
    # elif store_sqlite3_db:
    else:
        sqlite3_db_path = create_sqlite2_db_in_file(in_gff_path, disable_infer_genes, disable_infer_transcripts,
                                                    sqlite3_db_path, tmp_sqlite3_db_path, logger)

        sqlite3_db_genes = load_sqlite3_db(sqlite3_db_path, logger)
    # else:
    #     logger.print_timestamp()
    #     logger.info('Creating sqlite3 db by gffutils...')
    #
    #     id_spec = get_id_spec()
    #
    #     sqlite3_db_genes = gffutils.create_db(in_gff_path, dbfn=':memory:', force=True, verbose=True, merge_strategy="merge",
    #                                   transform=transform, force_merge_fields=['source'], id_spec=id_spec,
    #                                   disable_infer_genes=disable_infer_genes,
    #                                   disable_infer_transcripts=disable_infer_transcripts)
    #
    #     logger.info('Done.')

    # add transcripts equal genes at prokaryotes:
    sqlite3_db_genes = add_transcripts_prokaryotes(sqlite3_db_genes, logger, prokaryote)

    # add exons equal transcripts at prokaryotes:
    sqlite3_db_genes = add_exons_prokaryotes(sqlite3_db_genes, logger, prokaryote)

    # for feature in sqlite3_db_genes.all_features():
    #     print(feature)

    # for gene in sqlite3_db_genes.features_of_type(type_genes):
    #     print(gene)
    #     for transcript in sqlite3_db_genes.children(gene.id, featuretype=type_isoforms):
    #         print(transcript)
    #         for exon in sqlite3_db_genes.children(transcript.id, featuretype=type_exons):
    #             print(exon)
    #     print '\n\n'

    return sqlite3_db_genes


# def get_type_features(sqlite3_db_genes, default_type_genes, default_type_isoforms, default_type_exons, logger):
#     type_genes = get_type_genes(sqlite3_db_genes, default_type_genes, logger)
#
#     type_isoforms = get_type_isoforms(sqlite3_db_genes, default_type_isoforms, type_genes, logger)
#
#     type_exons = get_type_exons(sqlite3_db_genes, default_type_exons, type_genes, type_isoforms, logger)
#
#     return type_genes, type_isoforms, type_exons
#
#
# def get_type_genes(sqlite3_db_genes, default_type_genes, logger):
#     if len(list(sqlite3_db_genes.features_of_type(default_type_genes))) != 0:
#         type_genes = default_type_genes
#     else:
#         type_genes = None
#         logger.error('Annotated genes not founded.', exit_with_code=2, to_stderr=True)
#     return type_genes
#
#
# def get_type_isoforms(sqlite3_db_genes, default_type_isoforms, type_genes, logger):
#     if len(list(sqlite3_db_genes.features_of_type(default_type_isoforms))) != 0:
#         type_isoforms = default_type_isoforms
#     elif len(list(sqlite3_db_genes.features_of_type(type_genes))) != 0:
#         type_isoforms = type_genes
#     else:
#         type_isoforms = None
#         logger.error('Annotated isoforms not founded.', exit_with_code=2, to_stderr=True)
#     return type_isoforms
#
#
# def get_type_exons(sqlite3_db_genes, default_type_exons, type_genes, type_isoforms, logger):
#     if len(list(sqlite3_db_genes.features_of_type(default_type_exons))) != 0:
#         type_exons = default_type_exons
#    # for prokaryotes or bad annotated species:
#     elif len(list(sqlite3_db_genes.features_of_type(type_isoforms))) != 0:
#         type_exons = type_isoforms
#     elif len(list(sqlite3_db_genes.features_of_type(type_genes))) != 0:
#         type_exons = type_genes
#     else:
#         type_exons = None
#         logger.error('Annotated exons not founded.', exit_with_code=2, to_stderr=True)
#     return type_exons


def add_transcripts_prokaryotes(genes_db, logger, prokaryote=False):
    logger.print_timestamp()
    logger.info('Add transcripts equal genes in genes database... ')
    missed_transcripts = []

    for gene in genes_db.features_of_type(default_type_genes):
        if len(list(genes_db.children(gene.id, featuretype=default_type_isoforms))) == 0 or prokaryote:
            transcript = \
                gffutils.Feature(seqid=gene.seqid, source='equal_gene', featuretype='transcript',
                                 start=gene.start, end=gene.end, score=gene.score,
                                 strand=gene.strand, frame=gene.frame,
                                 attributes={'ID': [gene.id + '_t'], 'gene_id': [gene.id]},
                                 id=gene.id + '_t')
            missed_transcripts.append(transcript)

    for transcript in missed_transcripts:
        # logger.debug(str(transcript))
        gene = genes_db[transcript.attributes['gene_id'][0]]
        genes_db.add_relation(gene, transcript, level=1, child_func=child_func(gene, transcript))
        genes_db.update(UtilsGeneral.get_iterator([transcript]), merge_strategy='create_unique')

    logger.info('Done.')

    return genes_db


def add_exons_prokaryotes(genes_db, logger, prokaryote=False):
    logger.print_timestamp()
    logger.info('Add exons equal transcripts in genes database... ')

    missed_exons = []

    for gene in genes_db.features_of_type(default_type_genes):
        for transcript in genes_db.children(gene.id, featuretype=default_type_isoforms):
            if len(list(genes_db.children(transcript.id, featuretype=default_type_exons))) == 0 or prokaryote:
                exon = \
                    gffutils.Feature(seqid=transcript.seqid, source='equal_transcript', featuretype='exon',
                                     start=transcript.start, end=transcript.end, score=transcript.score,
                                     strand=transcript.strand, frame=transcript.frame,
                                     attributes={'ID': [transcript.id + '_e'], 'transcript_id': [transcript.id], 'gene_id': [gene.id]},
                                     id=transcript.id + '_e')
                missed_exons.append(exon)

    for exon in missed_exons:
        # logger.debug(str(exon))
        transcript = genes_db[exon.attributes['transcript_id'][0]]
        gene = genes_db[exon.attributes['gene_id'][0]]
        genes_db.add_relation(transcript, exon, level=1, child_func=child_func(transcript, exon))
        genes_db.add_relation(gene, exon, level=2, child_func=child_func(gene, exon))
        genes_db.update(UtilsGeneral.get_iterator([exon]), merge_strategy='create_unique')

    logger.info('Done.')

    return genes_db


def get_fa_isoforms(sqlite3_db_genes, reference_dict, logger):
    logger.print_timestamp()
    logger.info("Extracting isoforms sequences...")

    inconsistent_ref_db = False

    isoforms_dict = {}

    isoforms = list(sqlite3_db_genes.features_of_type(default_type_isoforms))

    for transcript in isoforms:
        if transcript.seqid not in reference_dict:
            continue

        isoforms_dict[transcript.id] = ''

        exons = list(sqlite3_db_genes.children(transcript.id, featuretype=default_type_exons, order_by='start'))

        for exon in exons:
            # start, end: 1-based coordinates; start must be <= end
            isoforms_dict[transcript.id] += reference_dict[exon.seqid][exon.start - 1:exon.end]

        if transcript.strand == '-':
            isoforms_dict[transcript.id] = UtilsGeneral.rev_comp(isoforms_dict[transcript.id])

        if len(isoforms_dict[transcript.id]) == 0:
            inconsistent_ref_db = True

            logger.debug('Inconsistent length chromosome / scaffold and transcript start / end: {} skipped'.format(transcript.id))

            isoforms_dict.pop(transcript.id)

    logger.info('Done.')

    if inconsistent_ref_db:
        logger.warning('Inconsistent reference sequences and genes database')

    return isoforms_dict


def clear_gtf_by_reference_chr(in_gtf, ids_chrs, tmp_dir, gtf_label, logger):
    out_gtf = os.path.join(tmp_dir, gtf_label + '.cleared.gtf')

    out_handle = open(out_gtf, 'w')

    skipped_chrs = False

    with open(in_gtf, 'r') as in_handle:
        for line in in_handle:
            if len(line.strip()) == 0:
                continue
            gtf_chr = line.strip().split()[0]
            if gtf_chr[0] == '#':
                continue
            if gtf_chr in ids_chrs:
                out_handle.write(line)
            else:
                skipped_chrs = True
                logger.debug('{} skipped'.format(gtf_chr))

    out_handle.close()

    if skipped_chrs:
        logger.warning('Some chromosomes / scaffolds or patches in GTF / GFF are skipped because they are absent in FASTA file with reference')

    return out_gtf
