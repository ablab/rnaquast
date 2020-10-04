__author__ = 'letovesnoi'

import os

import glob

import subprocess

import shutil

from datetime import datetime

from metrics import IsoformsCoverage

from general import UtilsPipeline

class CegmaMetrics():

    def __init__(self):
        self.complete_completeness = None
        self.partial_completeness = None


    def get_metrics(self, args_threads, transcripts_path, tmp_dir, label, logger):
        logger.print_timestamp()
        logger.info('  Getting CEGMA (Core Eukaryotic Genes Mapping Approach) metrics...')

        # run CEGMA:
        cegma_completeness_report_path = self.get_cegma_completeness_report(args_threads, transcripts_path, tmp_dir, label, logger)

        # parse CEGMA completeness report:
        if cegma_completeness_report_path is not None:
            self.complete_completeness = self.get_complete_completeness(cegma_completeness_report_path)
            self.partial_completeness = self.get_partial_completeness(cegma_completeness_report_path)

            logger.info('  Done.')


    # get CEGMA (Core Eukaryotic Genes Mapping Approach) completeness report
    def get_cegma_completeness_report(self, args_threads, transcripts_path, tmp_dir, label, logger):
        cegma_completeness_report_path = None

        command = 'cegma -g {} -o {} -T {}'.format(transcripts_path, os.path.join(tmp_dir, label), args_threads)
        exit_code = subprocess.call(command, shell=True)
        if exit_code != 0:
            logger.warning(message='CEGMA failed! Please install and add to PATH CEGMA: '
                                   '(http://korflab.ucdavis.edu/Datasets/cegma/)')
        else:
            cegma_completeness_report_path = os.path.join(tmp_dir, '{}.completeness_report'.format(label))
            # command = 'mv {} {}'.format(os.path.join(tmp_dir, 'output.completeness_report'), cegma_completeness_report_path)
            # subprocess.call(command, shell=True)

        return cegma_completeness_report_path


    def get_complete_completeness(self, cegma_completeness_report_path):
        complete = []
        with open(cegma_completeness_report_path, 'r') as fin:
            for line in fin:
                if ' Complete ' in line:
                    tmp_list = line.strip().split(' ')
                    for element in tmp_list:
                        if element != '':
                            complete.append(element)
                    complete_completeness = complete[2]
                    return complete_completeness


    def get_partial_completeness(self, cegma_completeness_report_path):
        partial = []
        with open(cegma_completeness_report_path, 'r') as fin:
            for line in fin:
                if ' Partial ' in line:
                    tmp_list = line.strip().split(' ')
                    for element in tmp_list:
                        if element != '':
                            partial.append(element)
                    partial_completeness = partial[2]
                    return partial_completeness


    def print_metrics(self, path_txt_completeness, logger):
        logger.info('      Getting CEGMA transcripts metrics report...')

        with open(path_txt_completeness, 'a') as fout:
            fout.write('METRICS OF TRANSCRIPTS WITH CEGMA:\n')
            fout.write('{:<100}'.format('Complete %completeness') + str(self.complete_completeness) + '\n')
            fout.write('{:<100}'.format('Partial %completeness') + str(self.partial_completeness) + '\n\n')

        logger.info('        saved to {}'.format(path_txt_completeness))


class BuscoMetrics():

    def __init__(self, busco_completeness_report_path):
        self.complete_completeness = BuscoMetrics.get_complete_completeness(busco_completeness_report_path)
        self.partial_completeness = BuscoMetrics.get_partial_completeness(busco_completeness_report_path)


    # get BUSCO (Benchmarking Universal Single-Copy Orthologs) results
    @classmethod
    def get_busco_metrics(cls, args_busco, args_prokaryote, args_threads, transcripts_path, tmp_dir, label, logger, log_dir):
        busco_metrics = None

        busco_completeness_report_path = \
            BuscoMetrics.get_busco_completeness_report(args_busco, args_prokaryote, args_threads, transcripts_path, tmp_dir,
                                                       label, logger, log_dir)

        if busco_completeness_report_path is not None:
            busco_metrics = cls(busco_completeness_report_path)

        return busco_metrics


    @classmethod
    def get_busco_completeness_report(cls, args_busco, args_prokaryote, args_threads, transcripts_path, tmp_dir, label, logger, log_dir):
        # run BUSCO:
        logger.print_timestamp()
        logger.info('  Running BUSCO (Benchmarking Universal Single-Copy Orthologs)...')

        # since busco download lineage data and log to current directory
        initial_dir = os.getcwd()
        os.chdir(tmp_dir)

        out_name = label + '_BUSCO'
        out_dirpath = os.path.join(tmp_dir, out_name)
        log_out = os.path.join(log_dir, '{}.busco.out.log'.format(label))
        log_err = os.path.join(log_dir, '{}.busco.err.log'.format(label))
        busco_completeness_report_path = None

        program_name = 'busco'
        command = '{busco} -o {output_name} --out_path {out_path} -i {transcripts} -m transcriptome -f -c {threads}'.\
            format(busco=program_name, output_name=out_name, out_path=tmp_dir, transcripts=transcripts_path, threads=args_threads)
        if args_busco == 'auto-lineage':
            # type_str = ''
            if args_prokaryote:
                type_str = 'prok'
            else:
                type_str = 'euk'
            busco_completeness_report_mask = \
                os.path.join(out_dirpath, 'auto_lineage', 'run_{}aryota_*'.format(type_str), 'short_summary*.txt')
            command += ' --' + args_busco + '-' + type_str
        else:
            if os.path.exists(args_busco):
                if not os.path.isabs(args_busco):
                    args_busco = os.path.abspath(args_busco)
            busco_completeness_report_mask = \
                os.path.join(out_dirpath, 'short_summary.*.txt')
            command += ' -l ' + args_busco
        command += ' 1>> {log_out_1} 2>> {log_out_2}'.format(log_out_1=log_out, log_out_2=log_err)

        logger.debug('    ' + command)

        exit_code = subprocess.call(command, shell=True)

        os.chdir(initial_dir)

        for file in glob.glob(busco_completeness_report_mask):
            busco_completeness_report_path = file
        if exit_code != 0 or not busco_completeness_report_path:
            logger.error(message='{} failed for {}!'.format(program_name, label))
        else:
            logger.info('    saved to {}.'.format(busco_completeness_report_path))

        logger.info('    logs can be found in {} and {}.'.format(log_out, log_err))

        return busco_completeness_report_path


    @classmethod
    def get_complete_completeness(cls, busco_completeness_report_path):
        with open(busco_completeness_report_path, 'r') as fin_handle:
            for line in fin_handle:
                tmp = line.strip().split()
                if len(tmp) == 4 and tmp[1] == 'Complete' and tmp[2] == 'BUSCOs' and tmp[3] == '(C)':
                    complete = int(tmp[0])
                if len(tmp) == 5 and tmp[1] == 'Total' and tmp[2] == 'BUSCO' and tmp[3] == 'groups' and tmp[4] == 'searched':
                    total = int(tmp[0])

        complete_completeness = complete * 100.0 / total

        return complete_completeness

    @classmethod
    def get_partial_completeness(cls, busco_completeness_report_path):
        with open(busco_completeness_report_path, 'r') as fin_handle:
            for line in fin_handle:
                tmp = line.strip().split()
                if len(tmp) == 4 and tmp[1] == 'Fragmented' and tmp[2] == 'BUSCOs' and tmp[3] == '(F)':
                    part = int(tmp[0])
                if len(tmp) == 5 and tmp[1] == 'Total' and tmp[2] =='BUSCO' and tmp[3] == 'groups' and tmp[4] == 'searched':
                        total = int(tmp[0])

        partial_completeness = part * 100.0 / total

        return partial_completeness


    def print_metrics(self, path_txt_completeness, logger):
        logger.info('      Getting BUSCO transcripts metrics report...')

        with open(path_txt_completeness, 'a') as fout:
            fout.write('METRICS OF TRANSCRIPTS WITH BUSCO:\n')
            fout.write('{:<100}'.format('Complete %completeness') + str(self.complete_completeness) + '\n')
            fout.write('{:<100}'.format('Partial %completeness') + str(self.partial_completeness) + '\n\n')

        logger.info('        saved to {}'.format(path_txt_completeness))


class GeneMarkS_TMetrics():

    def __init__(self, GeneMarkS_T_report_path):
        self.genes = self.get_genes_from_report(GeneMarkS_T_report_path)

    # get GeneMarkS-T results
    @classmethod
    def get_GeneMarkS_T_metrics(cls, type_organism, args_threads, args_ss, transcripts_path, tmp_dir,
                                label, logger, log_dir):
        geneMarkS_T_metrics = None

        GeneMarkS_T_report_path = \
            GeneMarkS_TMetrics.get_GeneMarkS_T_report(type_organism, args_threads, args_ss, transcripts_path, tmp_dir,
                                                      label, logger, log_dir)

        if GeneMarkS_T_report_path is not None:
            geneMarkS_T_metrics = cls(GeneMarkS_T_report_path)

        return geneMarkS_T_metrics


    @classmethod
    def get_GeneMarkS_T_report(cls, type_organism, args_threads, args_ss, transcripts_path, tmp_dir, label, logger, log_dir):
        # run GeneMarkS-T:
        logger.print_timestamp()
        logger.info('  Running GeneMarkS-T (Gene Prediction in Transcripts)...')

        out_dir_path = UtilsPipeline.create_empty_folder(os.path.join(tmp_dir, label + '_GeneMarkS-T'))

        initial_dir = os.getcwd()

        os.chdir(out_dir_path)

        transcripts_name = os.path.split(transcripts_path)[-1]
        tmp_GeneMarkS_T_report_path = os.path.join(out_dir_path, transcripts_name + '.lst')
        GeneMarkS_T_report_path = None
        tmp_log_path = os.path.join(out_dir_path, 'gms.log')
        log_out_path = os.path.join(log_dir, label + '.GeneMarkS_T.out.log')
        log_err_path = os.path.join(log_dir, label + '.GeneMarkS_T.err.log')

        GeneMarkS_T_run = 'gmst.pl'
        command = '{} {} --output {} 2>> {}'.format(GeneMarkS_T_run, transcripts_path, tmp_GeneMarkS_T_report_path,
                                                    log_err_path)
        if type_organism == 'prokaryotes':
            command += ' --prok'

        if args_ss:
            command += ' --strand direct'

        logger.debug(command)

        exit_code = subprocess.call(command, shell=True)

        os.chdir(initial_dir)

        if exit_code != 0 or not os.path.exists(tmp_GeneMarkS_T_report_path):
            logger.error(message='GeneMarkS-T failed for {}!'.format(label))
        else:
            GeneMarkS_T_report_path = tmp_GeneMarkS_T_report_path

            logger.info('    saved to {}'.format(GeneMarkS_T_report_path))

        if os.path.exists(tmp_log_path):
            shutil.move(tmp_log_path, log_out_path)

        logger.info('    logs can be found in {} and {}.'.format(log_out_path, log_err_path))

        return GeneMarkS_T_report_path


    def get_genes_from_report(self, GeneMarkS_T_report_path):
        # GeneMarkS_T_all = []
        self.genes = 0

        with open(GeneMarkS_T_report_path, 'r') as in_handle:
            f_read_genes = False
            attributes = {}
            for line in in_handle:
                tmp_list = line.strip().split()

                if tmp_list == []:
                    f_read_genes = False
                    continue

                if tmp_list[0:2] == ['Model', 'information:']:
                    attributes = {}
                elif tmp_list[0:3] == ['FASTA', 'definition', 'line:']:
                    attributes['t_name'] = tmp_list[3:6]
                elif tmp_list == ['#', 'Length']:
                    attributes['genes'] = []
                    f_read_genes = True
                elif f_read_genes:
                    attributes['genes'].append({})
                    attributes['genes'][-1]['Gene'] = tmp_list[0]
                    attributes['genes'][-1]['Strand'] = tmp_list[1]
                    attributes['genes'][-1]['LeftEnd'] = tmp_list[2]
                    attributes['genes'][-1]['RightEnd'] = tmp_list[3]
                    attributes['genes'][-1]['Gene'] = tmp_list[4]
                    attributes['genes'][-1]['Class'] = tmp_list[5]

                    # GeneMarkS_T_all.append(attributes)

                    self.genes += 1

        return self.genes


class AssemblyCompletenessMetrics():
    """Class of annotation coverage metrics by aligned transcripts"""

    def __init__(self, args):
        self.isoforms_coverage = None
        if (args.gene_db is not None or args.gtf is not None) and args.alignment is not None and args.reference is not None and args.transcripts is not None:
            # INITIALIZE ASSEMBLY COMPLETENESS METRICS WITH ALIGNMENT AND ANNOTATION:
            self.isoforms_coverage = IsoformsCoverage.IsoformsCoverage()

        # INITIALIZE CEGMA METRICS:
        # self.cegma_metrics = None
        # if args.cegma:
        #     self.cegma_metrics = CegmaMetrics()

        # INITIALIZE BUSCO METRICS:
        self.busco_metrics = None

        # INITIALIZE GeneMarkS-T METRICS:
        self.geneMarkS_T_metrics = None


    # UPDATE COVERAGE OF ANNOTATION BY SPECIFIC ISOFORM:
    def update_assembly_completeness_metrics(self, sqlite3_db_genes, internal_isoforms_coverage, id_mapped_isoform):
        start_time = datetime.now()

        # update coverage of annotated isoforms:
        self.isoforms_coverage.update_isoforms_coverage_by_specific_isoform(sqlite3_db_genes, internal_isoforms_coverage, id_mapped_isoform)

        elapsed_time = datetime.now() - start_time

        return elapsed_time


    def get_assembly_completeness_metrics(self, args, sqlite3_db_genes, db_genes_metrics, reads_coverage, transcripts_path,
                                          type_organism, tmp_dir, label, threads, WELL_FULLY_COVERAGE_THRESHOLDS,
                                          logger, log_dir):
        # get average metrics of coverage of annotated isoforms (included exons coverages) by aligned transcripts:
        logger.info('  Getting SENSITIVITY metrics...')

        # get coverages of annotated isoforms:
        if self.isoforms_coverage is not None:
            self.isoforms_coverage.get_isoforms_coverage(sqlite3_db_genes, db_genes_metrics, reads_coverage, WELL_FULLY_COVERAGE_THRESHOLDS)

        # CEGMA:
        # if self.cegma_metrics is not None:
        #     self.cegma_metrics.get_metrics(args.threads, transcripts_path, tmp_dir, self.label, logger)

        if args.busco:
            self.busco_metrics = \
                BuscoMetrics.get_busco_metrics(args.busco, args.prokaryote, threads, transcripts_path, tmp_dir, label, logger, log_dir)

        if args.gene_mark or not ((args.gtf is not None or args.gene_db is not None) and args.alignment is not None and
                                          args.reference is not None and args.transcripts is not None):
            self.geneMarkS_T_metrics = \
                GeneMarkS_TMetrics.get_GeneMarkS_T_metrics(type_organism, threads, args.strand_specific,
                                                           transcripts_path, tmp_dir, label, logger, log_dir)

        logger.info('  Done.')