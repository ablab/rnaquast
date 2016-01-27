__author__ = 'letovesnoi'

import os

import subprocess

from datetime import datetime

import IsoformsCoverage

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

    def __init__(self):
        self.complete_completeness = None
        self.partial_completeness = None


    # get BUSCO (Benchmarking Universal Single-Copy Orthologs) results
    def get_metrics(self, args_clade, args_threads, transcripts_path, tmp_dir, label, logger):
        busco_completeness_report_path = \
            self.get_busco_completeness_report(args_clade, args_threads, transcripts_path, tmp_dir, label, logger)

        if busco_completeness_report_path is not None:
            self.complete_completeness = self.get_complete_completeness(busco_completeness_report_path)
            self.partial_completeness = self.get_partial_completeness(busco_completeness_report_path)

            logger.info('  Done.')


    def get_busco_completeness_report(self, args_clade, args_threads, transcripts_path, tmp_dir, label, logger):
        busco_completeness_report_path = None

        # run BUSCO:
        logger.print_timestamp()
        logger.info('  Running BUSCO (Benchmarking Universal Single-Copy Orthologs)...')

        initial_dir = os.getcwd()

        os.chdir(tmp_dir)

        out_name = label + '_BUSCO'
        out_dirpath = os.path.join(tmp_dir, 'run_' + out_name)
        log_out = 'busco.log'
        logger.info('    Logging to {}...'.format(os.path.join(out_dirpath, log_out)))

        busco_run = 'BUSCO_v1.1b1.py'
        command = \
            '{busco} -o {output_dir} -in {transcripts} -l {clade} -m trans -f -c {threads} > {log_out}'.\
                format(busco=busco_run, output_dir=out_name, transcripts=transcripts_path, clade=args_clade,
                       threads=args_threads, log_out=log_out)
        logger.debug(command)

        exit_code = subprocess.call(command, shell=True)

        os.chdir(initial_dir)

        if exit_code != 0:
            logger.warning(message='BUSCO failed! Please install and add to PATH BUSCO requirements.')
        else:

            busco_completeness_report_path = os.path.join(out_dirpath, 'short_summary_{}_BUSCO'.format(label))

        return busco_completeness_report_path


    def get_complete_completeness(self, busco_completeness_report_path):
        with open(busco_completeness_report_path, 'r') as fin_handle:
            for line in fin_handle:
                tmp = line.strip().split()
                if len(tmp) == 4 and tmp[1] == 'Complete' and tmp[2] == 'Single-copy' and tmp[3] == 'BUSCOs':
                    complete = int(tmp[0])
                if len(tmp) == 5 and tmp[1] == 'Total' and tmp[2] == 'BUSCO' and tmp[3] == 'groups' and tmp[4] == 'searched':
                    total = int(tmp[0])

        complete_completeness = complete * 100.0 / total

        return complete_completeness


    def get_partial_completeness(self, busco_completeness_report_path):
        with open(busco_completeness_report_path, 'r') as fin_handle:
            for line in fin_handle:
                tmp = line.strip().split()
                if len(tmp) == 3 and tmp[1] == 'Fragmented' and tmp[2] == 'BUSCOs':
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

    def __init__(self):
        self.genes = None


    # get GeneMarkS-T results
    def get_metrics(self, args_clade, args_threads, transcripts_path, tmp_dir, label, logger):
        GeneMarkS_T_report_path = \
            self.get_GeneMarkS_T_report(args_clade, args_threads, transcripts_path, tmp_dir, label, logger)

        if GeneMarkS_T_report_path is not None:
            self.genes = self.get_genes_from_report(GeneMarkS_T_report_path)

            logger.info('    saved to {}'.format(GeneMarkS_T_report_path))


    def get_GeneMarkS_T_report(self, type_organism, args_threads, transcripts_path, tmp_dir, label, logger):
        GeneMarkS_T_report_path = None

        # run GeneMarkS-T:
        logger.print_timestamp()
        logger.info('  Running GeneMarkS-T (Gene Prediction in Transcripts)...')

        out_dir_path = UtilsPipeline.create_empty_folder(os.path.join(tmp_dir, label + '_GeneMarkS-T'))
        transcripts_name = os.path.split(transcripts_path)[-1]
        GeneMarkS_T_report_path_tmp = os.path.join(out_dir_path, transcripts_name + '.lst')
        log_path = os.path.join(out_dir_path, 'gms.log')
        logger.info('    Logging to {}...'.format(log_path))

        GeneMarkS_T_run = 'gmst.pl'
        command = '{} {} --output {}'.format(GeneMarkS_T_run, transcripts_path, GeneMarkS_T_report_path_tmp)
        if type_organism == 'prokaryotes':
            command += ' --prok'
        logger.debug(command)

        exit_code = subprocess.call(command, shell=True)
        if exit_code != 0:
            logger.warning(message='GeneMarkS-T failed! Please install and add to PATH GeneMarkS-T requirements.')
        else:
            GeneMarkS_T_report_path = GeneMarkS_T_report_path_tmp

        return GeneMarkS_T_report_path


    def get_genes_from_report(self, GeneMarkS_T_report_path):

        return self.genes


class AssemblyCompletenessMetrics():
    """Class of annotation coverage metrics by aligned transcripts"""

    def __init__(self):
        # INITIALIZE ASSEMBLY COMPLETENESS METRICS WITH ALIGNMENT AND ANNOTATION:
        self.isoforms_coverage = IsoformsCoverage.IsoformsCoverage()


    # UPDATE COVERAGE OF ANNOTATION BY SPECIFIC ISOFORM:
    def update_assembly_completeness_metrics(self, sqlite3_db_genes, internal_isoforms_coverage, id_mapped_isoform):
        start_time = datetime.now()

        # update coverage of annotated isoforms:
        self.isoforms_coverage.update_isoforms_coverage_by_specific_isoform(sqlite3_db_genes, internal_isoforms_coverage, id_mapped_isoform)

        elapsed_time = datetime.now() - start_time

        return elapsed_time


    def get_assembly_completeness_metrics(self, sqlite3_db_genes, tot_isoforms_len, reads_coverage,
                                          WELL_FULLY_COVERAGE_THRESHOLDS, logger):
        # get average metrics of coverage of annotated isoforms (included exons coverages) by aligned transcripts:
        logger.info('  Getting SENSITIVITY metrics...')

        # get coverages of annotated isoforms:
        self.isoforms_coverage.get_isoforms_coverage(sqlite3_db_genes, tot_isoforms_len, reads_coverage, WELL_FULLY_COVERAGE_THRESHOLDS)

        logger.info('  Done.')


    def print_fully_assembled_isoforms(self, path_fully_assembled_list, logger):
        logger.info('    Getting Fully assembled isoforms list...')

        with open(path_fully_assembled_list, 'w') as fout:
            for id_isoform in self.isoforms_coverage.ids_fully_assembled_isoforms:
                fout.write(id_isoform + '\n')

        logger.info('      saved to {}'.format(path_fully_assembled_list))