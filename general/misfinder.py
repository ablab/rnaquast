#!/usr/bin/env python3
# -*- coding: utf-8 -*- 
import os
import sys
import argparse
import subprocess
import shutil
import copy
import numpy
import timeit
from time import gmtime, strftime
from Bio import SeqIO

NUMBER_FOR_ESTIMATION = 1000
FREQUENCY_FOR_PRINTING_UPDATES = 100000

MULT_COEFF_BEG_AND_END_READS = 15
MAX_DIST_BEG_END_READS_PARTS_OF_FRAGMENT = 5

DEVIATION_TLEN_IN_SIGMA = 3
MAX_DIST_TLEN_LEFT_PARTS_OF_FRAGMENT = 100
MAX_DIST_TLEN_RIGHT_PARTS_OF_FRAGMENT = 100
MIN_NUM_OF_READS_TLEN = 10

DEVIATION_ATYPICAL_COV_IN_SIGMA = 3
MIN_LEN_ATYPICAL_COV_FRAGMENT = 50
MAX_DIST_ATYPICAL_COV_PARTS_OF_FRAGMENT = 10

def is_tool(name):
    try:
        devnull = open(os.devnull)
        subprocess.Popen([name], stdout=devnull, stderr=devnull).communicate()
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            return False
    return True


def build_alignment_bowtie(bowtie_path, data_name, ref_file, reads1_file, reads2_file):
    if not os.path.exists(data_name):
        os.mkdir(data_name)
    if not os.path.exists(os.path.join(data_name, "index")):
        os.mkdir(os.path.join(data_name, "index"))
    subprocess.call([os.path.join(bowtie_path, "bowtie2-build"), ref_file, os.path.join(data_name, "index", "index")])
    subprocess.call([os.path.join(bowtie_path, "bowtie2"), "-q", "-x", os.path.join(data_name, "index", "index"), "-1", reads1_file, "-2", reads2_file, "-S", os.path.join(data_name, "aligned.sam")])
    shutil.rmtree(os.path.join(data_name, "index"))
    return os.path.join(data_name, "aligned.sam")


def build_alignment_bwa(bwa_path, data_name, ref_file, reads1_file, reads2_file):
    if not os.path.exists(data_name):
        os.mkdir(data_name)
    subprocess.call([os.path.join(bwa_path, "bwa"), "index", ref_file])
    with open(os.path.join(data_name, "aligned.sam"), "w") as sam_file:
        subprocess.call([os.path.join(bwa_path, "bwa"), "mem", ref_file, reads1_file, reads2_file], stdout = sam_file)
    os.remove(ref_file + ".amb")
    os.remove(ref_file + ".ann")
    os.remove(ref_file + ".bwt")
    os.remove(ref_file + ".pac")
    os.remove(ref_file + ".sa")
    return os.path.join(data_name, "aligned.sam")


def find_ref_len(ref_file):
    handle = open(ref_file, "rU")
    records = list(SeqIO.parse(handle, "fasta"))
    handle.close()
    coverage = {}
    for elem in records:
        name = elem.name
        seq = elem.seq
        empty_list = []
        for i in xrange(len(seq) + 1):
            empty_list.append([])
        coverage[name] = empty_list
    return coverage


def find_read_end(cigar): 
    if cigar == "*":
        return 0
    cigar_parsed = []
    i = 0
    count = 0
    while i < len(cigar):
        while (cigar[i + count]).isdigit():
            count += 1
        cigar_parsed.append((int(cigar[i:i+count]), cigar[i+count]))
        i += count + 1
        count = 0
    end = 0
    for (num, act) in cigar_parsed:
        if (act == "M") or (act == "D") or (act == "N")  or (act == "X")  or (act == "="):
            end += num
    return end


def parse_sam_record(line):
    if not line.startswith('$') and not line.startswith('@'):
        line = line.strip().split()
        record = {
        'QNAME' : line[0],
        'FLAG'  : int(line[1]),
        'RNAME' : line[2],
        'POS'   : int(line[3]),
#        'MAPQ'  : int(line[4]),
        'CIGAR' : line[5],
#        'RNEXT' : line[6],
#        'PNEXT' : int(line[7]),
        'TLEN'  : int(line[8]),
#        'SEQ'   : line[9],
#        'QUAL'  : line[10],
#        'optional' : []
        }
#        for optional in line[11:]:
#            record['optional'].append(optional.split(':'))
        return record


def one_of_mate_reads_unmapped(record):
    # 0x4 segment unmapped, 0x8 next segment in the template unmapped
    return (record ['FLAG'] & int(0x4)) or (record ['FLAG'] & int(0x8))


def estimate_tlen(sam_file): 
    num = 0
    tlen_list = []
    for line in open(sam_file):
        record = parse_sam_record(line)
        if record:
            cur_flag = record['FLAG']
            cur_tlen = record['TLEN']
            if one_of_mate_reads_unmapped(record) or (cur_flag & int(0x800)): # 0x800 supplementary alignment
                continue
            tlen_list.append(abs(cur_tlen))
            num += 1
        if num >= NUMBER_FOR_ESTIMATION:
            break
    average_tlen = numpy.mean(tlen_list)
    sigma_tlen = numpy.std(tlen_list)
    return (average_tlen, sigma_tlen)


def parse_file_with_reads(ref_file, sam_file):
    (average_tlen, sigma_tlen) = estimate_tlen(sam_file)
    max_tlen = average_tlen + DEVIATION_TLEN_IN_SIGMA * sigma_tlen
    large_tlen = {}
    coverage = find_ref_len(ref_file)
    number = 0
    failed_pairs = 0
    for line in open(sam_file):
        record = parse_sam_record(line)
        if record:
            number += 1
            if number % FREQUENCY_FOR_PRINTING_UPDATES == 0:
                print "Processing file... ", number, " lines"
            if one_of_mate_reads_unmapped(record):
                failed_pairs += 1
                continue
            begin = record['POS']
            end = record['POS'] + find_read_end(record['CIGAR'].strip())
            ref_name = record['RNAME']
            if abs(record['TLEN']) > max_tlen:
                if not ((record ['FLAG']) & int(0x800)): # 0x800 supplementary alignment
                    if not large_tlen.has_key(ref_name):
                        large_tlen[ref_name] = []
                    seq_name = record['QNAME']
                    if ((record ['FLAG']) & int(0x40)): # 0x40 the first segment in the template
                        seq_name += "/1"
                    elif ((record ['FLAG']) & int(0x80)): # 0x80 the last segment in the template
                        seq_name += "/2"
                    if (len(large_tlen[ref_name]) == 0) or (large_tlen[ref_name][-1][-1][2][:-2] != seq_name[:-2]):
                        large_tlen[ref_name].append([(begin, end, seq_name, record['TLEN'])])
                    else:
                        large_tlen[ref_name][-1].append((begin, end, seq_name, record['TLEN']))
            coverage[ref_name][begin].append(end)
    for ref_name in large_tlen:
        for i in xrange(len(large_tlen[ref_name])):
            large_tlen[ref_name][i] = sorted(large_tlen[ref_name][i])
        large_tlen[ref_name] = sorted(large_tlen[ref_name])
    print "Total number of reads is: ", number
    print "Failed pairs ", failed_pairs / 2
    return (coverage, large_tlen)


def estimate_coverage(coverage):
    coverage_list = []
    list_of_ending_reads = []
    list_of_begining_reads = []
    for (ref_name, ref_coverage) in coverage.items():
        cur_reads = {}
        number = 0
        for i in xrange(1, len(ref_coverage)):
            number += 1
            list_of_begining_reads.append(len(ref_coverage[i]))
            for read_end in ref_coverage[i]:
                if not cur_reads.has_key(read_end):
                    cur_reads[read_end] = 0
                cur_reads[read_end] += 1
            i_reads = 0
            for read_end in cur_reads:
                i_reads += cur_reads[read_end]
            coverage_list.append(i_reads)
            if cur_reads.has_key(i):
                list_of_ending_reads.append(cur_reads[i])
                del cur_reads[i]
            if number >= NUMBER_FOR_ESTIMATION:
                break
        if number >= NUMBER_FOR_ESTIMATION:
            break
    average_coverage = numpy.mean(coverage_list)
    sigma_coverage = numpy.std(coverage_list)
    number_of_begining_reads = numpy.mean(list_of_begining_reads)
    number_of_ending_reads = numpy.mean(list_of_ending_reads)
    return ((average_coverage, sigma_coverage), (number_of_begining_reads, number_of_ending_reads))


def get_atypical_coverage_and_number_of_beg_and_end_reads(coverage):
    ((average_coverage, sigma_coverage), (number_of_begining_reads, number_of_ending_reads)) = estimate_coverage(coverage)
    min_cov = average_coverage - DEVIATION_ATYPICAL_COV_IN_SIGMA * sigma_coverage
    max_cov = average_coverage + DEVIATION_ATYPICAL_COV_IN_SIGMA * sigma_coverage
    number_of_begining_reads *= MULT_COEFF_BEG_AND_END_READS
    number_of_ending_reads *= MULT_COEFF_BEG_AND_END_READS
    coverage_number = 0
    sum_of_len = 0
    atypical_coverage = {}
    beg_and_end_reads = {}
    for ref_name, ref_coverage in coverage.items():
        cur_reads = {}
        atypical_coverage[ref_name] = []
        beg_and_end_reads[ref_name] = []
        number = 0
        sum_of_len += len(ref_coverage)
        for i in xrange(1, len(ref_coverage)):
            number += 1
            for read_end in ref_coverage[i]:
                if not cur_reads.has_key(read_end):
                    cur_reads[read_end] = 0
                cur_reads[read_end] += 1
            if len(ref_coverage[i]) > number_of_begining_reads:
                beg_and_end_reads[ref_name].append((i, "beg", len(ref_coverage[i])))
            i_reads = 0
            for read_end in cur_reads:
                i_reads += cur_reads[read_end]
            coverage_number += i_reads
            if (i_reads > max_cov) or (i_reads < min_cov):
                atypical_coverage[ref_name].append((i, i_reads))
            if cur_reads.has_key(i):
                if cur_reads[i] > number_of_ending_reads:
                    beg_and_end_reads[ref_name].append((i, "end", cur_reads[i]))
                del cur_reads[i]
    coverage_in_average = round(float(coverage_number)/ sum_of_len, 1)
    return (atypical_coverage, beg_and_end_reads, coverage_in_average)


def write_tlen_to_file(file_name_tlen, large_tlen_list):
    print "writing tlen to file..."
    output_file_tlen = open(file_name_tlen,'w')
    for (ref_name, large_tlen_list_ref) in large_tlen_list.items():
        for pair in large_tlen_list_ref:
            output_file_tlen.write("ref: " + ref_name)
            if len(pair) != 2:
                for (begin, end, seq_name, tlen_val) in pair:
                    output_file_tlen.write("\tseq: " + seq_name + "\tbeg: " + str(begin) + "\tend: " + str(end) + "\ttlen: " + str(tlen_val))
            else:
                ((begin1, end1, seq_name1, tlen_val1), (begin2, end2, seq_name2, tlen_val2)) = pair
                output_file_tlen.write("\tbeg1: " + str(begin1) + "\tend1: " + str(end1) + "\tbeg2: " + str(begin2) + "\tend2: " + str(end2) + "\t" + seq_name1 + "\t" + seq_name2)
            output_file_tlen.write("\n")
    output_file_tlen.close()


def write_atypical_coverage_to_file(file_name_coverage, atypical_coverage):
    print "writing atypical coverage to file..."
    output_file_coverage = open(file_name_coverage,'w')
    for (ref_name, atypical_cov) in atypical_coverage.items():
        for (pos, num_of_reads) in atypical_cov:
            output_file_coverage.write("ref: " + ref_name + "\tpos: " + str(pos) + "\tcov: " + str(num_of_reads) + "\n")
    output_file_coverage.close()


def write_beg_and_end_reads_to_file(file_name_reads, beg_and_end_reads):
    print "writing number of begining and ending reads to file..."
    output_file_reads = open(file_name_reads, 'w')
    for ref_name, atypical_reads in beg_and_end_reads.items():
        for (pos, reads_type, num_of_reads) in atypical_reads:
            output_file_reads.write("ref: " + ref_name + "\tpos: " + str(pos) + "\t" + reads_type + "\tnum: " + str(num_of_reads) + "\n")
    output_file_reads.close()


def analize_atypical_cov(file_name, atypical_coverage, coverage_num):
    atypical_cov_fragments = {}
    for (ref_name, atypical_cov) in atypical_coverage.items():
        if len(atypical_cov) == 0:
            continue
        atypical_cov_fragments[ref_name] = []
        last_pos = atypical_cov[0][0]
        fragment_start = last_pos
        cur_sum = 0
        cur_len = 0
        for (pos, num_of_reads) in atypical_cov:
            if pos <= last_pos + MAX_DIST_ATYPICAL_COV_PARTS_OF_FRAGMENT:
                cur_sum += num_of_reads
                cur_len += 1
            else:
                if min((last_pos - fragment_start), cur_len) >= MIN_LEN_ATYPICAL_COV_FRAGMENT:
                    aver_coverage = round(float(cur_sum)/cur_len, 1)
                    atypical_cov_fragments[ref_name].append((fragment_start, last_pos - fragment_start + 1, aver_coverage))
                fragment_start = pos
                cur_sum = num_of_reads
                cur_len = 1
            last_pos = pos
        if min((last_pos - fragment_start), cur_len) >= MIN_LEN_ATYPICAL_COV_FRAGMENT:
            aver_coverage = round(float(cur_sum)/cur_len, 1)
            atypical_cov_fragments[ref_name].append((fragment_start, last_pos - fragment_start + 1, aver_coverage))
    output_file = open(file_name, 'a')
    output_file.write('Fragments with low or high coverage. avg coverage={}\n'.format(coverage_num))
    output_file.write('start\tend\tavg\n')
    for (ref_name, fragments) in atypical_cov_fragments.items():
        if len(fragments) > 0:
            output_file.write('query\t' + ref_name + '\n')
        for (start, length, aver_cov) in fragments:
            output_file.write(str(start) + '\t' + str(start + length - 1) + '\t' + str(aver_cov) + '\n')
    output_file.close()

def analize_large_tlen(file_name, large_tlen):
    tlen_fragments = {}
    for (ref_name, large_tlen_ref) in large_tlen.items():
        if len(large_tlen_ref) == 0:
            continue
        tlen_fragments[ref_name] = []

        def separate_fragments(fragment):
            def compare(pair1, pair2):
                (begin_right1, end_right1, seq_name_right1, tlen_val_right1) = pair1[-1]
                (begin_right2, end_right2, seq_name_right2, tlen_val_right2) = pair2[-1]
                if begin_right1 < begin_right2:
                    return -1
                elif begin_right1 > begin_right2:
                    return 1
                else:
                    return 0

            cur_fragment = sorted(fragment, cmp=compare)
            last_pair = cur_fragment[0]
            last_right_begin = last_pair[-1][0]
            left_fragment_beg = last_pair[0][0]
            left_fragment_end = last_pair[0][0]
            right_fragment_beg = last_pair[-1][0]
            right_fragment_end = last_pair[-1][0]
            num_of_reads = 0
            for pair in cur_fragment:
                (begin_left, end_left, seq_name_left, tlen_val_left) = pair[0]
                (begin_right, end_right, seq_name_right, tlen_val_right) = pair[-1]
                if begin_right <= last_right_begin + MAX_DIST_TLEN_RIGHT_PARTS_OF_FRAGMENT:
                    num_of_reads += 1
                    left_fragment_beg = min(left_fragment_beg, begin_left)
                    right_fragment_beg = min(right_fragment_beg, begin_right)
                    left_fragment_end = max(left_fragment_end, begin_left)
                    right_fragment_end = max(right_fragment_end, begin_right)
                else:
                    if num_of_reads >= MIN_NUM_OF_READS_TLEN:
                        tlen_fragments[ref_name].append((left_fragment_beg, left_fragment_end, right_fragment_beg, right_fragment_end, num_of_reads))
                    last_pair = pair
                    last_right_begin = last_pair[-1][0]
                    left_fragment_beg = last_pair[0][0]
                    left_fragment_end = last_pair[0][0]
                    right_fragment_beg = last_pair[-1][0]
                    right_fragment_end = last_pair[-1][0]
                    num_of_reads = 1
                last_right_begin = begin_right
            if num_of_reads >= MIN_NUM_OF_READS_TLEN:
                tlen_fragments[ref_name].append((left_fragment_beg, left_fragment_end, right_fragment_beg, right_fragment_end, num_of_reads))

        last_pair = large_tlen_ref[0]
        last_left_begin = last_pair[0][0]
        cur_fragment = []
        for pair in large_tlen_ref:
            (begin_left, end_left, seq_name_left, tlen_val_left) = pair[0]
            (begin_right, end_right, seq_name_right, tlen_val_right) = pair[-1]
            if abs(tlen_val_left) != abs(tlen_val_right):
                continue
            if begin_left <= last_left_begin + MAX_DIST_TLEN_LEFT_PARTS_OF_FRAGMENT:
                cur_fragment.append(pair)
            else:
                if len(cur_fragment) >= MIN_NUM_OF_READS_TLEN:
                    separate_fragments(cur_fragment)
                last_pair = pair
                last_left_begin = last_pair[0][0]
                cur_fragment = [pair]
            last_left_begin = begin_left
        if len(cur_fragment) >= MIN_NUM_OF_READS_TLEN:
            separate_fragments(cur_fragment)
    output_file = open(file_name, 'a')
    output_file.write('Fragments with large TLEN:\n')
    output_file.write('start\tend\tnumber\n')
    for (ref_name, fragments) in tlen_fragments.items():
        if len(fragments) > 0:
            output_file.write('query\t' + ref_name + '\n')
        for (left_beg, left_end, right_beg, right_end, num_of_reads) in fragments:
            output_file.write(str(left_beg) + '\t' + str(left_end) + '\t' + str(num_of_reads) + '\n'
                              + str(right_beg) + '\t' + str(right_end) + '\t' + str(num_of_reads) + '\n')
    output_file.close()


def analize_beg_and_end_reads(file_name, beg_and_end_reads):
    beg_and_end_pos = {}
    for (ref_name, atypical_reads) in beg_and_end_reads.items():
        if len(atypical_reads) == 0:
            continue
        beg_and_end_pos[ref_name] = []
        last_pos = atypical_reads[0][0]
        num = 0
        has_beg = False
        has_end = False
        fragment_start = last_pos
        for (pos, reads_type, num_of_reads) in atypical_reads:
            if pos <= last_pos + MAX_DIST_BEG_END_READS_PARTS_OF_FRAGMENT:
                num += num_of_reads
                if reads_type == "beg":
                    has_beg = True
                else:
                    has_end = True
            else:
                if (has_end and has_beg):
                    beg_and_end_pos[ref_name].append((fragment_start, last_pos, num))
                fragment_start = pos
                num = num_of_reads
                has_beg = False
                has_end = False
                if reads_type == "beg":
                    has_beg = True
                else:
                    has_end = True
            last_pos = pos
        if (has_end and has_beg):
            beg_and_end_pos[ref_name].append((fragment_start, last_pos, num))
    output_file = open(file_name, 'w')
    output_file.write("Fragments with large number of reads ends and begins:\n")
    output_file.write('start\tend\tnumber\n')
    for (ref_name, fragments) in beg_and_end_pos.items():
        if len(fragments) > 0:
            output_file.write('query\t' + ref_name + '\n')
        for (fragment_beg, fragment_end, num) in fragments:
            output_file.write(str(fragment_beg) + '\t' + str(fragment_end) + '\t' + str(num) + '\n')
    output_file.close()

if __name__ == "__main__":
    start = timeit.timeit()
    if len(sys.argv) == 1:
        print "Usage:", sys.argv[0], "-o <name of output directory> -c <file with contigs> -r1 <left reads> -r2 <right reads>"
        print "Please use the --help option to get more usage information."
        exit()

    parser = argparse.ArgumentParser(prog = sys.argv[0], description='Misassemblies detection without reference.')
    parser.add_argument("-v", "--version", action='version', version='%(prog)s 1.0')
    parser.add_argument("-o", "--out_dir", help="name of output directory", required=True)
    parser.add_argument("-c", "--contigs", help="file with contigs", required=True)
    parser.add_argument("-r1", "--reads_1", help="file with left reads", required=True)
    parser.add_argument("-r2", "--reads_2", help="file with right reads", required=True)
    parser.add_argument("--bowtie", help="path to bowtie-2 dir", default="")
    parser.add_argument("--bwa", help="path to bwa dir", default="")
    parser.add_argument("-a", "--aligner", help="method to align", choices=["bowtie", "bwa"], default="bwa")

    parser.add_argument('-s', '--sam_file', help='file with alignments in sam format', required=False)

    args = parser.parse_args()

    bowtie_path = args.bowtie
    bwa_path = args.bwa
    data_name = args.out_dir
    ref_file = args.contigs
    reads1_file = args.reads_1
    reads2_file = args.reads_2

    print strftime("Processing started... %Y-%m-%d %H:%M:%S", gmtime())

    if not os.path.isfile(ref_file):
        print >> sys.stderr, "Not a file\t" + ref_file
        exit(1)
    if not os.path.isfile(reads1_file):
        print >> sys.stderr, "Not a file\t" + reads1_file
        exit(1)
    if not os.path.isfile(reads2_file):
        print >> sys.stderr, "Not a file\t" + reads2_file
        exit(1)

    if args.sam_file:
        sam_file = args.sam_file
    else:
        if args.aligner == "bowtie":
            if (not is_tool(os.path.join(bowtie_path, "bowtie2"))) or (not is_tool(os.path.join(bowtie_path, "bowtie2-build"))):
                print >> sys.stderr, "Bowtie2 is not installed.\n Please, check http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#obtaining-bowtie-2"
                exit(1)
            sam_file = build_alignment_bowtie(bowtie_path, data_name, ref_file, reads1_file, reads2_file)
        elif args.aligner == "bwa":
            if not is_tool(os.path.join(bwa_path, "bwa")):
                print >> sys.stderr, "Bwa is not installed.\n Please, dowload from http://sourceforge.net/projects/bio-bwa/files/"
                exit(1)
            sam_file = build_alignment_bwa(bwa_path, data_name, ref_file, reads1_file, reads2_file)
        sam_file = os.path.join(data_name, "aligned.sam")

    (coverage, large_tlen) = parse_file_with_reads(ref_file, sam_file)
    (atypical_coverage, beg_and_end_reads, coverage_num) = get_atypical_coverage_and_number_of_beg_and_end_reads(coverage)
    write_tlen_to_file(os.path.join(data_name, "large_insert_size.txt"), large_tlen)
    write_atypical_coverage_to_file(os.path.join(data_name, "atypical_coverage.txt"), atypical_coverage)
    write_beg_and_end_reads_to_file(os.path.join(data_name, "beg_and_end_reads.txt"), beg_and_end_reads)
    analize_beg_and_end_reads(os.path.join(data_name, "result.txt"), beg_and_end_reads)
    analize_large_tlen(os.path.join(data_name, "result.txt"), large_tlen)
    analize_atypical_cov(os.path.join(data_name, "result.txt"), atypical_coverage, coverage_num)

    #os.remove(sam_file)

    print strftime("Processing finished! %Y-%m-%d %H:%M:%S", gmtime())
    print "Time elapsed: ", timeit.timeit() - start, "sec"


