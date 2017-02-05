#!/usr/bin/env python

__author__ = 'letovesnoi'

import sys
import os
import subprocess

import datetime

import argparse

import logging


def get_path(output_dir, program_name):
    # if --output_dir not use, create default path for results:
    out_path = program_name + '_' + datetime.datetime.now().strftime('%Y_%m_%d_%H_%M_%S') + '.txt'
    if output_dir is not None:
        out_path = os.path.join(output_dir, out_path)
    return out_path

def get_arguments():
    # use --help for running without arguments:
    if len(sys.argv) == 1:
        command = 'python {} -h'.format(sys.argv[0])
        subprocess.call(command, shell=True)
        sys.exit(0)

    parser = \
        argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                description="Paste short reports together utils\n"
                                            "\nUsage:\npython %(prog)s --reports SHORT_REPORTS_TXT --output",
                                conflict_handler='resolve',
                                prog=sys.argv[0])

    # INPUT DATA:
    parser.add_argument('-r', '--reports', help='Two files with two short rnaQUAST reports', type=str, nargs=2)

    parser.add_argument('-o', '--output_dir', help='Directory to store result [default: current_directory/paste_together_results_<datetime>]', type=str)

    parser.add_argument('-d', '--debug', help='Report detailed information, typically used only for detecting problems.', action='store_true')

    args = parser.parse_args()

    return args


def get_metrics(path):
    names = []
    metrics = {}

    fin = open(path, 'r')

    for line in fin:
        tmp = [x for x in line.strip().split(2 * ' ') if x != '']
        if len(tmp) > 1:
            names.append(tmp[0])
            metrics[tmp[0]] = tmp[1:]
        else:
            names.append(line)

    num = len(metrics[metrics.keys()[0]])

    fin.close()

    logging.debug('names: ' + str(names))
    logging.debug('metrics: ' + str(metrics))
    logging.debug('num: ' + str(num))

    return names, metrics, num


def paste_2_metrics(path1, path2):
    new_metrics = {}

    names1, metrics1, num1 = get_metrics(path1)
    names2, metrics2, num2 = get_metrics(path2)

    for name in names1:
        if '\n' in name:
            continue
        if name in names2:
            new_metrics[name] = metrics1[name] + metrics2[name]
        else:
            new_metrics[name] = metrics1[name] + num2 * ['*']

    logging.debug('new metrics: ' + str(new_metrics))

    return names1, new_metrics, num1


def print_txt(names, metrics, path_txt):
    fout_txt = open(path_txt, 'w')

    for name in names:
        if '\n' in name:
            fout_txt.write(name)
        elif name in metrics:
            column_width_str = '{:<50}'
            fout_txt.write(column_width_str.format(name))
            for value in metrics[name]:
                column_width_str = '{:<25}'
                fout_txt.write(column_width_str.format(value))
            fout_txt.write('\n')
    fout_txt.close()

    logging.info('  saved to\n' + '    ' + '{}\n'.format(path_txt))


if __name__ == '__main__':
    try:
        program_name = 'paste_together'

        args = get_arguments()

        if args.debug:
            logging.basicConfig(level=logging.DEBUG)
        else:
            logging.basicConfig(level=logging.INFO)

        logging.info('Getting SHORT SUMMARY report...')

        out_path = get_path(args.output_dir, program_name)

        curr_report = args.reports[0]
        for report in args.reports[1:]:
            curr_names, curr_metrics, num = paste_2_metrics(curr_report, report)
            curr_report = report

        print_txt(curr_names, curr_metrics, out_path)

    except Exception:
        _, exc_value, _ = sys.exc_info()
        logging.exception(exc_value)
        logging.error('Exception caught!')
