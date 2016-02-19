############################################################################
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

from __future__ import with_statement
import os
import sys
from datetime import datetime
from general import rqconfig
from general import UtilsGeneral

import logging

_loggers = {}

def get_logger(name):
    if name in _loggers.keys():
        return _loggers[name]
    else:
        _loggers[name] = QLogger(name)
        return _loggers[name]


class QLogger(object):
    _logger = None  # logging.getLogger('quast')
    _name = ''
    _log_fpath = ''
    _start_time = None
    _indent_val = 0
    _num_notices = 0
    _num_warnings = 0
    _num_nf_errors = 0

    _list_warnings = []
    _list_non_fatal_errors = []

    def __init__(self, name):
        self._name = name
        self._logger = logging.getLogger(name)
        self._logger.setLevel(logging.DEBUG)

    def set_up_console_handler(self, indent_val=0, debug=False):
        self._indent_val = indent_val

        for handler in self._logger.handlers:
            if isinstance(handler, logging.StreamHandler):
                self._logger.removeHandler(handler)

        console_handler = logging.StreamHandler(sys.stdout, )
        console_handler.setFormatter(logging.Formatter(indent_val * '  ' + '%(message)s'))
        console_handler.setLevel(logging.DEBUG if debug else logging.INFO)
        self._logger.addHandler(console_handler)

    def set_up_debug_level(self):
        for handler in self._logger.handlers:
            handler.setLevel(logging.DEBUG)

    def set_up_file_handler(self, output_dirpath):
        for handler in self._logger.handlers:
            if isinstance(handler, logging.FileHandler):
                self._logger.removeHandler(handler)

        self._log_fpath = os.path.join(output_dirpath, self._name + '.log')
        file_handler = logging.FileHandler(self._log_fpath, mode='w')
        file_handler.setLevel(logging.INFO)

        self._logger.addHandler(file_handler)

    def start(self):
        if self._indent_val == 0:
            self._logger.info('')
            self.print_version()
            self._logger.info('')
            self.print_system_info()

        self._start_time = self.print_timestamp('Started: ')
        self._logger.info('')
        self._logger.info('Logging to ' + self._log_fpath)

    def finish_up(self, numbers=None, check_test=False):
        self._logger.info('\nLog saved to ' + self._log_fpath)

        finish_time = self.print_timestamp('Finished: ')
        self._logger.info('Elapsed time: ' + str(finish_time - self._start_time))
        if numbers:
            self.print_numbers_of_notifications(prefix="Total ", numbers=numbers)
        else:
            self.print_numbers_of_notifications()

        if self._num_warnings != 0:
            self._logger.info('\nWARNINGS:')
            for warning_msg in self._list_warnings:
                self._logger.info('  WARNING: ' + str(warning_msg) if warning_msg else '')

        if self._num_nf_errors != 0:
            self._logger.info('\nnon-fatal ERRORs:')
            for non_fatal_error_msg in self._list_non_fatal_errors:
                self._logger.info('  non-fatal ERROR: ' + str(non_fatal_error_msg) if non_fatal_error_msg else '')

        self._logger.info('\nThank you for using {}!\n'.format(self._name))
        if check_test:
            if (numbers is not None and numbers[2] > 0) or self._num_nf_errors > 0:
                self._logger.info('\nTEST FAILED! Please find non-fatal errors in the log and try to fix them!')
            elif (numbers is not None and numbers[1] > 0) or self._num_warnings > 0:
                self._logger.info('\nTEST PASSED with WARNINGS!')
            else:
                self._logger.info('\nTEST PASSED!')

        for handler in self._logger.handlers:
            self._logger.removeHandler(handler)

        global _loggers
        del _loggers[self._name]

    def debug(self, message='', indent=''):
        self._logger.debug(indent + message)

    def info(self, message='', indent=''):
        self._logger.info(indent + message)

    def info_to_file(self, message='', indent=''):
        for handler in self._logger.handlers:
            if isinstance(handler, logging.FileHandler):
                self._logger.removeHandler(handler)

        with open(self._log_fpath, 'a') as f:
            f.write(indent + message + '\n')

        file_handler = logging.FileHandler(self._log_fpath, mode='a')
        file_handler.setLevel(logging.DEBUG)

        self._logger.addHandler(file_handler)

    def notice(self, message='', indent=''):
        self._num_notices += 1
        self._logger.info(indent + ('NOTICE: ' + str(message) if message else ''))

    def warning(self, message='', indent=''):
        self._num_warnings += 1
        self._list_warnings.append(message)
        self._logger.warning(indent + ('WARNING: ' + str(message) if message else ''))

    def error(self, message='', exit_with_code=0, to_stderr=False, indent='', fake_if_nested_run=False):
        if fake_if_nested_run and self._indent_val > 0:
            self.info('')
            self.notice(message)
            return

        if message:
            msg = indent + 'ERROR! ' + str(message)
            if exit_with_code:
                msg += "\n\nIn case you have troubles running rnaQUAST, you can write to rnaquast_support@ablab.spbau.ru\n"
                msg += "Please provide us with {}.log file from the output directory.".format(self._name)
        else:
            msg = ''

        if to_stderr or not self._logger.handlers:
            sys.stderr.write(msg + '\n')
        else:
            self._logger.error('')
            self._logger.error(msg)

        if exit_with_code:
            exit(exit_with_code)
        else:
            self._num_nf_errors += 1
            self._list_non_fatal_errors.append(message)

    def exception(self, e, exit_code=0):
        if self._logger.handlers:
            self._logger.error('')
            self._logger.exception(e)
        else:
            print >> sys.stderr, str(e)

        if exit_code:
            exit(exit_code)

    def print_command_line(self, args, indent='',
                           wrap_after=80, only_if_debug=False):
        text = ''
        line = indent

        for i, arg in enumerate(args):
            if ' ' in arg or '\t' in arg:
                args[i] = "'" + arg + "'"

            line += arg

            if i == len(args) - 1:
                text += line

            elif wrap_after is not None and len(line) > wrap_after:
                text += line + ' \\\n'
                line = ' ' * len(indent)

            else:
                line += ' '

        if only_if_debug:
            self.debug(text)
        else:
            self.info(text)

    def print_timestamp(self, message=''):
        now = datetime.now()
        current_time = now.strftime("%Y-%m-%d %H:%M:%S")
        self.info('')
        self.info(message + current_time)
        return now

    def print_version(self, to_stderr=False):
        version, build = UtilsGeneral.get_version(rqconfig.GENERAL_LOCATION)
        if to_stderr:
            print >> sys.stderr, "Version", str(version) + (", " + str(build) if build != "unknown" else "")
        else:
            self.info("Version " + str(version) + (", " + str(build) if build != "unknown" else ""))

    def print_system_info(self):
        self._logger.info("System information:")
        try:
            import platform
            self._logger.info("  OS: " + platform.platform() + " (%s)" % rqconfig.platform_name)
            self._logger.info("  Python version: " + str(sys.version_info[0]) + "." + str(sys.version_info[1]) + '.'\
                      + str(sys.version_info[2]))
            import multiprocessing
            self._logger.info("  CPUs number: " + str(multiprocessing.cpu_count()))
        except:
            self._logger.info("  Problem occurred when getting system information")

    def print_numbers_of_notifications(self, prefix="", numbers=None):
        if not numbers:
            numbers = (self._num_notices, self._num_warnings, self._num_nf_errors)
        self._logger.info(prefix + "NOTICEs: %d; WARNINGs: %d; non-fatal ERRORs: %d" %
                          numbers)

    def get_numbers_of_notifications(self):
        return (self._num_notices, self._num_warnings, self._num_nf_errors)

    def print_input_files(self, args):
        self.print_timestamp()
        # TRANSCRIPTS:
        if args.transcripts != None:
            self.info('Transcripts:')
            for i_transcripts in range(len(args.labels)):
                self.info('  %s -> %s' % (args.transcripts[i_transcripts], args.labels[i_transcripts]))
        # REFERENCE:
        if args.reference != None:
            self.info('Reference:')
            self.info('  %s' % (args.reference))
        # GENE DATABASE:
        if args.gene_database != None:
            self.info('Gene database:')
            self.info('  %s' % (args.gene_database))
        # ALIGNMENT:
        if args.alignment != None:
            self.info('Alignments:')
            for i_transcripts in range(len(args.labels)):
                self.info('  %s' % (args.alignment[i_transcripts]))
        self.info()

    def print_path_results(self, args, separated_reports, comparison_report, short_report):
        self.info('\nTHE QUALITY OF TRANSCRIPTOME ASSEMBLY DONE. RESULTS: {}\n'.format(args.output_dir))
        self.info('Separated assemblies reports are saved to')
        for iTranscripts in range(len(separated_reports)):
            self.info('  {}'.format(separated_reports[iTranscripts].output_dir))
            #self.info('  {}'.format(metricsOfTranscripts[iTranscripts].label))
            #self.info('    TXT version basic isoforms metrics (tables) saved to {}'.format(separated_reports[iTranscripts].isoforms_metrics_report.path_txt_basic))
            #self.info('    TXT version basic transcripts metrics (tables) saved to {}'.format(separated_reports[iTranscripts].transcripts_metrics_report.path_txt_basic))
            #self.info('    TXT version simple transcripts metrics (tables) saved to {}'.format(separated_reports[iTranscripts].transcripts_metrics_report.path_txt_simple))
            #self.info('    TXT version assembly correctness (tables) saved to {}'.format(separated_reports[iTranscripts].transcripts_metrics_report.path_txt_correctness))
            #self.info('    TXT version assembly completeness (tables) saved to {}'.format(separated_reports[iTranscripts].isoforms_metrics_report.path_txt_completeness))
            #self.info('    Fusions and misassemblies saved to {}'.format(separated_reports[iTranscripts].fusion_misassemble_dir))
            #self.info('    Mapped coverages saved to {}'.format(separated_reports[iTranscripts].mapped_coverage_dir))
        if comparison_report != None:
            self.info('Comparison assemblies report are saved to {}'.format(comparison_report.output_dir))
        self.info('  PDF version (tables and plots) saved to {}'.format(short_report.path_pdf))
        self.info('  TXT version (tables and plots) saved to {}'.format(short_report.path_txt))
