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
        self.complete_completeness = 0.0
        self.partial_completeness = 0.0


    # get BUSCO (Benchmarking Universal Single-Copy Orthologs) results
    def get_busco_metrics(self, args_clade, args_threads, transcripts_path, tmp_dir, label, logger, log_dir):
        busco_completeness_report_path = \
            self.get_busco_completeness_report(args_clade, args_threads, transcripts_path, tmp_dir, label, logger, log_dir)

        if busco_completeness_report_path is not None:
            self.complete_completeness = self.get_complete_completeness(busco_completeness_report_path)
            self.partial_completeness = self.get_partial_completeness(busco_completeness_report_path)


    def get_busco_completeness_report(self, args_clade, args_threads, transcripts_path, tmp_dir, label, logger, log_dir):
        busco_completeness_report_path = None

        # run BUSCO:
        logger.print_timestamp()
        logger.info('  Running BUSCO (Benchmarking Universal Single-Copy Orthologs)...')

        initial_dir = os.getcwd()

        os.chdir(tmp_dir)

        out_name = label + '_BUSCO'
        out_dirpath = os.path.join(tmp_dir, 'run_' + out_name)
        log_out = os.path.join(log_dir, '{}.busco.out.log'.format(label))
        log_err = os.path.join(log_dir, '{}.busco.err.log'.format(label))

        program_name = 'BUSCO_v1.1b1.py'
        command = \
            '{busco} -o {output_dir} -in {transcripts} -l {clade} -m trans -f -c {threads} 1>> {log_out_1} 2>> {log_out_2}'.\
                format(busco=program_name, output_dir=out_name, transcripts=transcripts_path, clade=args_clade,
                       threads=args_threads, log_out_1=log_out, log_out_2=log_err)
        logger.debug('    ' + command)

        exit_code = subprocess.call(command, shell=True)

        os.chdir(initial_dir)

        if exit_code != 0:
            logger.error(message='{} failed!'.format(program_name))
        else:
            busco_completeness_report_path = os.path.join(out_dirpath, 'short_summary_{}_BUSCO'.format(label))

            logger.info('    saved to {}.'.format(busco_completeness_report_path))

        logger.info('    logs can be found in {} and {}.'.format(log_out, log_err))

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
        self.genes = 0


    # get GeneMarkS-T results
    def get_GeneMarkS_T_metrics(self, args_clade, args_threads, transcripts_path, tmp_dir, label, logger, log_dir):
        GeneMarkS_T_report_path = \
            self.get_GeneMarkS_T_report(args_clade, args_threads, transcripts_path, tmp_dir, label, logger, log_dir)

        if GeneMarkS_T_report_path is not None:
            self.genes = self.get_genes_from_report(GeneMarkS_T_report_path)


    def get_GeneMarkS_T_report(self, type_organism, args_threads, transcripts_path, tmp_dir, label, logger, log_dir):
        GeneMarkS_T_report_path = None

        # run GeneMarkS-T:
        logger.print_timestamp()
        logger.info('  Running GeneMarkS-T (Gene Prediction in Transcripts)...')

        out_dir_path = UtilsPipeline.create_empty_folder(os.path.join(tmp_dir, label + '_GeneMarkS-T'))
        transcripts_name = os.path.split(transcripts_path)[-1]
        GeneMarkS_T_report_path_tmp = os.path.join(out_dir_path, transcripts_name + '.lst')
        log_path = os.path.join(log_dir, label + '.GeneMarkS_T.log')

        GeneMarkS_T_run = 'gmst.pl'
        command = '{} {} --output {} 2>> {}'.format(GeneMarkS_T_run, transcripts_path, GeneMarkS_T_report_path_tmp,
                                                    log_path)
        if type_organism == 'prokaryotes':
            command += ' --prok'
        logger.debug(command)

        exit_code = subprocess.call(command, shell=True)
        if exit_code != 0:
            logger.error(message='GeneMarkS-T failed!')
        else:
            GeneMarkS_T_report_path = GeneMarkS_T_report_path_tmp

            logger.info('    saved to {}'.format(GeneMarkS_T_report_path))

        logger.info('  log can be found in {}.'.format(os.path.join(log_dir, log_path)))

        return GeneMarkS_T_report_path


    def get_genes_from_report(self, GeneMarkS_T_report_path):
        # GeneMarkS_T_all = []

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
        if args.gene_database is not None and args.alignment is not None and args.reference is not None and args.transcripts is not None:
            # INITIALIZE ASSEMBLY COMPLETENESS METRICS WITH ALIGNMENT AND ANNOTATION:
            self.isoforms_coverage = IsoformsCoverage.IsoformsCoverage()

        # INITIALIZE CEGMA METRICS:
        # self.cegma_metrics = None
        # if args.cegma:
        #     self.cegma_metrics = CegmaMetrics()

        # INITIALIZE BUSCO METRICS:
        self.busco_metrics = None
        if args.busco and args.clade is not None:
            self.busco_metrics = BuscoMetrics()

        self.gene_marks_t_metrics = None
        if args.gene_mark or not (args.gene_database is not None and args.alignment is not None and args.reference is not None and args.transcripts is not None):
            self.gene_marks_t_metrics = GeneMarkS_TMetrics()


    # UPDATE COVERAGE OF ANNOTATION BY SPECIFIC ISOFORM:
    def update_assembly_completeness_metrics(self, sqlite3_db_genes, internal_isoforms_coverage, id_mapped_isoform):
        start_time = datetime.now()

        # update coverage of annotated isoforms:
        self.isoforms_coverage.update_isoforms_coverage_by_specific_isoform(sqlite3_db_genes, internal_isoforms_coverage, id_mapped_isoform)

        elapsed_time = datetime.now() - start_time

        return elapsed_time


    def get_assembly_completeness_metrics(self, args_clade, threads, transcripts_path, tmp_dir, label, type_organism,
                                          sqlite3_db_genes, tot_isoforms_len, reads_coverage,
                                          WELL_FULLY_COVERAGE_THRESHOLDS, logger, log_dir):
        # get average metrics of coverage of annotated isoforms (included exons coverages) by aligned transcripts:
        logger.info('  Getting SENSITIVITY metrics...')

        # get coverages of annotated isoforms:
        if self.isoforms_coverage is not None:
            self.isoforms_coverage.get_isoforms_coverage(sqlite3_db_genes, tot_isoforms_len, reads_coverage, WELL_FULLY_COVERAGE_THRESHOLDS)

        # GET ASSEMBLY COMPLETENESS METRICS:
        # CEGMA:
        # if self.cegma_metrics is not None:
        #     self.cegma_metrics.get_metrics(args.threads, transcripts_path, tmp_dir, self.label, logger)

        # BUSCO:
        if self.busco_metrics is not None:
            self.busco_metrics.get_busco_metrics(args_clade, threads, transcripts_path, tmp_dir, label, logger, log_dir)

        if self.gene_marks_t_metrics is not None:
            self.gene_marks_t_metrics.get_GeneMarkS_T_metrics(type_organism, threads, transcripts_path, tmp_dir, label, logger, log_dir)

        logger.info('  Done.')