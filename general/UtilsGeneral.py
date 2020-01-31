__author__ = 'lenk'

import os
import sys
import subprocess
try:
    from itertools import imap
except ImportError:
    imap = map

def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def check_python_version(SUPPORTED_PYTHON_VERSIONS):
    if sys.version[0:3] not in SUPPORTED_PYTHON_VERSIONS:
        sys.stderr.write("ERROR! Python version " + sys.version[0:3] + " is not supported!\n" +\
                         "Supported versions are " + ", ".join(SUPPORTED_PYTHON_VERSIONS) + "\n")
        sys.exit(1)


def get_version(location):
    version_fpath = os.path.join(location, 'VERSION')
    version = "unknown"
    build = "unknown"
    if os.path.isfile(version_fpath):
        version_file = open(version_fpath)
        version = version_file.readline()
        if version:
            version = version.strip()
        else:
            version = "unknown"
        build = version_file.readline()
        if build:
            build = build.strip().lower()
        else:
            build = "unknown"
    return version, build


def get_version_by_key(program_name, key, tmp_dir, v_ident=None, b_ident=None):
    version = 'unknown'
    build = 'unknown'

    version_path = os.path.join(tmp_dir, program_name + '.version.log')
    command = program_name + ' ' + key + ' 1>> ' + version_path + ' 2>> ' + version_path
    subprocess.call(command, shell=True)

    with open(version_path, 'r') as fin:
        for line in fin:
            if v_ident and v_ident in line:
                version = line.strip().partition(v_ident)[2].split()[0]
            if b_ident and b_ident in line:
                build = line.strip().partition(b_ident)[2].split()[0]

    command = 'rm ' + version_path
    subprocess.call(command, shell=True)

    return version, build


def read_fasta(fasta_path):
    seq = ''
    with open(fasta_path, 'r') as fin:
        fin.readline()
        for line in fin:
            seq += line.strip()
    return seq


def read_dict_from_multi_fasta(fasta_path):
    dict = {}
    with open(fasta_path, 'r') as fin:
        curr_seqname = ''
        curr_seq = ''
        for line in fin:
            strip_line = line.strip()
            if strip_line == '':
                continue
            if strip_line[0] == '>':
                curr_seqname = strip_line.split(' ')[0][1:]
                dict[curr_seqname] = ''
            else:
                curr_seq = strip_line
                dict[curr_seqname] += curr_seq
    return dict


def list_to_dict(l):
    res = {}
    for e in l:
        res[e[0].split()[0].strip()] = e[1].strip()
    return res


def dict_to_list(dict):
    list = []
    for key in dict:
        list.append((key, dict[key]))
    return list


def comp(letter):
    return {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
            'U': 'A', 'R': 'Y', 'Y': 'R', 'K': 'M', 'M': 'K',
            'S': 'S', 'W': 'W', 'B': 'V', 'D': 'H', 'H': 'D',
            'V': 'B', '-': '-'}[letter.upper()]


def rev_comp(seq):
    import itertools
    return ''.join(imap(comp, seq[::-1]))


def hamming_dist(s1, s2):
    if len(s1) != len(s2):
        return -1
    d = 0
    for i in range(0, len(s1)):
        if s1[i] != s2[i]:
            d += 1
    return d


def get_ids_chr_scf_from_multi_fasta(fasta_path):
    ids_chr_scf = set([])
    with open(fasta_path, 'r') as fin:
        ids_chr_scf.add(fin.readline().strip().split(' ')[0][1:])
        seq = fin.readline().strip()
        while seq != '':
            while seq[0] != '>':
                seq = fin.readline().strip()
                if seq == '':
                    break
            if seq != '':
                ids_chr_scf.add(seq.split(' ')[0][1:])
                seq = fin.readline().strip()
    return ids_chr_scf


def read_fasta_chr_from_database(database):
    chrs_dict = {}
    with open(database, 'r') as fin:
        for fasta_file in fin:
            tmp = read_dict_from_multi_fasta(fasta_file.strip())
            chrs_dict[tmp.keys()[0]] = tmp.values()[0]
    return chrs_dict


def print_fasta_seq(fasta_path, name_seq, seq):
    with open(fasta_path, 'a') as fout:
        i_current = 0
        fout.write('>{}\r\n'.format(name_seq))
        while(i_current in range(len(seq))):
            if i_current + 80 < len(seq):
                fout.write(seq[i_current:i_current + 80])
                fout.write('\r\n')
                i_current += 80
            else:
                fout.write(seq[i_current:])
                fout.write('\r\n')
                break


# CREATE MEMORY TRACKED OBJECTS:
def profile_memory(args, reference_dict, db_genes_metrics, transcripts_metrics, separated_reports, comparison_report, logger):
    try:
        from pympler import classtracker
        from pympler import asizeof

        logger.print_timestamp()
        logger.info('Memory profile:')

        tracker = classtracker.ClassTracker()

        if args.reference is not None:
            logger.info('reference_dict [bit]\t' + str(asizeof.asizeof(reference_dict)))


        if args.gtf is not None or args.gene_db:
            tracker.track_object(db_genes_metrics)

        for i_transcripts in range(len(transcripts_metrics)):
            if args.transcripts is not None:
                tracker.track_object(transcripts_metrics[i_transcripts].basic_metrics)

            if args.alignment is not None and args.reference is not None and args.transcripts is not None:
                tracker.track_object(transcripts_metrics[i_transcripts].simple_metrics)

            if (args.gtf is not None or args.gene_db is not None) \
                    and args.alignment is not None and args.reference is not None and args.transcripts is not None:
                tracker.track_object(transcripts_metrics[i_transcripts].assembly_correctness_metrics)

                tracker.track_object(transcripts_metrics[i_transcripts].assembly_completeness_metrics.exons_coverage)
                tracker.track_object(transcripts_metrics[i_transcripts].assembly_completeness_metrics.isoforms_coverage)
                tracker.track_object(transcripts_metrics[i_transcripts].assembly_completeness_metrics)

            tracker.track_object(separated_reports[i_transcripts])

        if comparison_report is not None:
            tracker.track_object(comparison_report)

        tracker.create_snapshot()
        tracker.stats.print_summary()
    except Exception:
        logger.warning('Can\'t profile memory: please install Pympler.')


# get strand for strand specific or non strand specific option:
def get_strands(args, logger):
    if args.strand_specific:
        logger.info('Using strand specific transcripts...')
        return ['+', '-', '.']
    else:
        logger.info('Using non strand specific transcripts...')
        return [None]


# Get sorted-out index array sort_index and sorted array, corresponding to increasing values of the elements:
def get_order_indexes_elements(elements):
    sort_array = elements[:]

    sort_index = list(range(len(sort_array)))

    return qsort_i(0, len(elements) - 1, sort_array, sort_index)


# Quicksort:
# we sort sort_array and get sortIndex - array of indices, corresponding to increasing values of the sort_array:
def qsort_i(low, high, sort_array, sort_index):
    i = low
    j = high
    # pick a pivot element from the list:
    middle = sort_array[int((i + j) / 2)]
    # Reorder the list so that all elements with values less than the pivot come before the pivot,
    # while all elements with values greater than the pivot come after it:
    while i <= j:
        while sort_array[i] < middle:
            i += 1
        while sort_array[j] > middle:
            j -= 1
        if i <= j:
            # reorder elements of array:
            temp = sort_array[i]
            sort_array[i] = sort_array[j]
            sort_array[j] = temp
            # reorder indices by the value of array[index]:
            temp = sort_index[i]
            sort_index[i] = sort_index[j]
            sort_index[j] = temp
            i += 1
            j -= 1
    # Recursively apply the above steps to the sub-list of elements with smaller values and
    # separately the sub-list of elements with greater values:
    if low < j:
        qsort_i(low, j, sort_array, sort_index)
    if i < high:
        qsort_i(i, high, sort_array, sort_index)
    # return value of indices, corresponding to increasing values of the sort_array and sort_array:
    return sort_index, sort_array


# get position of element x in sort_array so that array with added this element stay sorted
# if array contains the same elements as x, return the first position in sortArray of the same elements
def get_bin_search_position_of_element(sort_array, x):
    first = 0
    last = len(sort_array)
    while first < last:
        middle = first + (last - first) // 2
        if x <= sort_array[middle]:
            last = middle
        else:
            first = middle + 1
    return last


def get_iterator(objects):
    for obj in objects:
        yield obj


# correction for fasta contained Y, W and etc:
def correct_nucl_seq(nucl_seq):
    corrected_nucl_seq = ''

    for i_nucl in range(len(nucl_seq)):
        nucl = nucl_seq[i_nucl]
        if nucl not in ['A', 'C', 'G', 'T', 'N']:
            nucl = 'N'
        corrected_nucl_seq += nucl

    return corrected_nucl_seq


from posixpath import curdir, sep, pardir, join, abspath, commonprefix


def relpath(path, start=curdir):
    """Return a relative version of a path"""
    if not path:
        raise ValueError("No path specified")
    start_list = abspath(start).split(sep)
    path_list = abspath(path).split(sep)
    # Work out how much of the filepath is shared by start and path.
    i = len(commonprefix([start_list, path_list]))
    rel_list = [pardir] * (len(start_list) - i) + path_list[i:]
    if not rel_list:
        return curdir
    return join(*rel_list)


def call_subprocess(args, logger, stdin=None, stdout=None, stderr=None,
                    indent='', only_if_debug=True, env=None):
    import subprocess

    printed_args = args[:]
    if stdin:
        printed_args += ['<', stdin.name]
    if stdout:
        printed_args += ['>>' if stdout.mode == 'a' else '>', stdout.name]
    if stderr:
        printed_args += ['2>>' if stderr.mode == 'a' else '2>', stderr.name]

    for i, arg in enumerate(printed_args):
        if arg.startswith(os.getcwd()):
            printed_args[i] = relpath(arg)

    logger.print_command_line(printed_args, indent, only_if_debug=only_if_debug)

    return_code = subprocess.call(args, stdin=stdin, stdout=stdout, stderr=stderr, env=env)

    if return_code != 0:
        logger.debug(' ' * len(indent) + 'The tool returned non-zero.' +
                     (' See ' + relpath(stderr.name) + ' for stderr.' if stderr else ''))
        # raise SubprocessException(printed_args, return_code)

    return return_code


def get_upper_case_fasta(reference, tmp_dir, logger):
    (tmp_dir_name, tmp_file_name) = os.path.split(reference)
    tmp_name = tmp_file_name[:tmp_file_name.rfind('.')]
    tmp_extension = tmp_file_name[tmp_file_name.rfind('.'):]
    out_name_fa = os.path.join('{}'.format(tmp_dir), tmp_name + '.upper' + tmp_extension)

    logger.print_timestamp()
    logger.info('Getting upper case fasta...')

    fout = open(out_name_fa, 'w')
    with open(reference, 'r') as fin:
        line = fin.readline()
        while line != '':
            if line[0] != '>':
                line = line.upper()
            fout.write(line)
            line = fin.readline()
    fout.close()

    logger.info('  saved to {}'.format(out_name_fa))

    return out_name_fa


def glue_scaffolds_together(reference, database, tmp_dir, logger):
    out_name_fa = os.path.join(tmp_dir, os.path.split(database)[-1] + '.glue.fa')
    logger.print_timestamp()
    logger.info('Gluing scaffolds files by cat for {}...'.format(database))
    f_pathes = ''
    with open(reference, 'r') as fin0:
        for f_path in fin0:
            f_pathes += f_path.strip() + ' '
    command = 'cat {} >> {}'.format(f_pathes, out_name_fa)
    exit_code = subprocess.call(command, shell=True)
    if exit_code != 0:
        logger.error(message='cat failed!', exit_with_code=exit_code, to_stderr=True)
        sys.exit(exit_code)
    logger.info('  saved to {}'.format(out_name_fa))
    return out_name_fa


def get_genome_len(reference_dict):
    genome_len = 0
    for id_chr in reference_dict:
        genome_len += len(reference_dict[id_chr])

    return genome_len