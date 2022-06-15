#!/usr/bin/python

# Import modules
import vcf as pyvcf
import pysam
import argparse
import multiprocessing as mp
import queue
import time
import sys
import collections
import subprocess
import os
import glob
import math
from math import log10
import pandas as pd
from pathlib import Path
from statistics import mean
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import numpy as np
import re
from sklearn.mixture import GaussianMixture
from matplotlib.backends.backend_pdf import PdfPages
import warnings
import configparser

# Get version from git
__version__ = 'v2.2.0'

binsize = 10000000
base_phred_quality = 0
mapq = 0
read_length = 150

# Set arguments
parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Put here a description.')
parser.add_argument('-g', '--germline', type=str, help='Input germline indexed vcf.gz file', required=True)
parser.add_argument('-s', '--somatic', type=str, help='Input somatic SMuRF_filtered indexed vcf.gz file', required=True)
parser.add_argument('-b', '--bam', action='append', nargs="*", type=str, help='Input bam file', required=True)
parser.add_argument('-o', '--output', type=str, help='Output prefix')
parser.add_argument('-f', '--format', action='append', nargs="*", type=str, choices=['vcf','bed','txt'], help='Output format', required=False)
parser.add_argument('-t', '--threads', type=int, help='Number of threads',default=8, required=False)
parser.add_argument('-v', '--version', action='version', version=__version__)
args = parser.parse_args()

# Flatten input list of bam files
args.bam = [x for l in args.bam for x in l]
args.format = [x for l in args.format for x in l]

# Read the vcf, fix and add fields to the header
somatic_vcf_reader = pyvcf.Reader(filename=args.somatic, encoding='utf-8')
somatic_vcf_name = args.output
# somatic_vcf_name = os.path.basename(args.somatic)
# somatic_vcf_name = somatic_vcf_name.replace(".vcf.gz","")
somatic_samples = somatic_vcf_reader.samples
germline_vcf_reader = pyvcf.Reader(filename=args.germline, encoding='utf-8')

# Create tmp directory if it does not exists
try:
    os.stat('./walker_tmp')
except:
    os.mkdir('./walker_tmp')

# Define global variables
bam_sample_names = collections.defaultdict(dict)
done_list = []
def main():
    global somatic_vcf_reader, contig_list
    somatic_vcf_reader = fix_vcf_header(somatic_vcf_reader)
    somatic_vcf_reader = add_vcf_header(somatic_vcf_reader)


    # Create contig bins
    contig_list = create_bins(somatic_vcf_reader)

    # Create an input queue with the contigs and an empty output queue
    q = mp.Queue()
    q_out = mp.Queue()
    for contig in contig_list:
        q.put(contig)

    # Create number of processes to parse the vcf file
    processes = [mp.Process(target=parse_chr_vcf, args=(q, q_out, somatic_vcf_reader, germline_vcf_reader, args.bam)) for x in range(int(args.threads))]

    for p in processes:
        p.start()
    liveprocs = list(processes)
    while liveprocs:
        time.sleep(5)
        try:
            while 1:
                done = q_out.get(block=False, timeout=1)
                done_list.append(done)
        except queue.Empty:
            pass
    # Give tasks a chance to put more data in
        time.sleep(10)
        if not q.empty():
            continue
        liveprocs = [p for p in liveprocs if p.is_alive()]

    for p in processes:
        p.join()

def add_vcf_header( vcf_reader ):
    """
    Function to add a new field to the vcf header
    Input: A vcf reader object
    Return: The vcf reader object with new headers added
    """
    # Formats
    vcf_reader.formats['WS'] = pyvcf.parser._Format('WS',None,'Float','Walker score')

    return(vcf_reader )

def create_bins(vcf_reader):
    """
    Function to create bins.
    Input: VCF reader
    Output: list of bins
    """
    contig_list = []
    for contig in vcf_reader.contigs:
        contig_length = vcf_reader.contigs[contig][1]
        bin_start = 1
        bin_end = 1
        while( bin_end < contig_length ):
            bin_end = bin_start+binsize-1
            if bin_end > contig_length:
                bin_end = contig_length
            if not contig.startswith("Un"):
               contig_list.append(contig+"-"+str(bin_start)+"-"+str(bin_end))
            bin_start = bin_end+1
    return( contig_list )

def parse_chr_vcf(q, q_out, somatic_contig_vcf_reader, germline_contig_vcf_reader, bams):
    """
    Function to parse the vcf per contig.
    Write the new record to a vcf file.
    Input: Queue object
    Input: Queue out object
    Input: VCF reader object
    Input: List with the bam names
    """
    while True:
        try:
            # Get contig one by one from the queue
            contig = q.get(block=False,timeout=1)
            contig_chr, contig_start, contig_end = contig.split("-")
            contig_start = int(contig_start)
            contig_end = int(contig_end)
            if 'txt' in args.format:
                contig_walker_txt_writer = open('./walker_tmp/{}.walker.txt'.format(contig),'w', encoding='utf-8')
            if 'bed' in args.format:
                contig_walker_bed_writer = open('./walker_tmp/{}.walker.bed'.format(contig),'w', encoding='utf-8')
            if 'vcf' in args.format:
                contig_walker_vcf_writer = pyvcf.Writer(open('./walker_tmp/{}.walker.vcf'.format(contig),'w', encoding='utf-8'), somatic_contig_vcf_reader)
            try:
                # Try to parse the specific contig from the vcf
                somatic_contig_vcf_reader.fetch(contig_chr, contig_start, contig_end)
            except:
                # Skip contig if it is not present in the vcf file
                q_out.put( contig )
                continue

            for record in somatic_contig_vcf_reader.fetch(contig_chr, contig_start, contig_end):
                # if record.is_indel:
                #     continue
                sample_names = somatic_samples
                for sample_name in sample_names:
                    cis = {}
                    trans = {}
                    other = {}
                    if sample_name == '' or sample_name not in somatic_samples:
                        continue
                    if record.genotype(sample_name).is_het:
                        somatic_alignments = get_somatic_alignments( record, sample_name )
                        # print( somatic_alignments )
                        # if somatic_alignments == False:
                        #     continue
                        for somatic_alignment in somatic_alignments:
                            somatic_base_info = somatic_alignment[0]
                            somatic_align = somatic_alignment[1]
                            gl_records = get_linked_gl_records( germline_contig_vcf_reader, sample_name, somatic_align.reference_name, somatic_align.reference_start, somatic_align.reference_end )
                            gl_records.extend( get_linked_gl_records( germline_contig_vcf_reader, sample_name, somatic_align.next_reference_name, somatic_align.next_reference_start, somatic_align.next_reference_start+read_length) )
                            if len(gl_records) > 0:
                                for gl_record in gl_records:
                                    linked_gl_base_info = get_linked_gl_base( gl_record, sample_name, somatic_align )
                                    if linked_gl_base_info is None:
                                        continue
                                    gl_record_info = gl_record.CHROM+":"+str(gl_record.POS)
                                    if gl_record_info not in cis:
                                        cis[gl_record_info] = []
                                    if gl_record_info not in trans:
                                        trans[gl_record_info] = []
                                    if gl_record_info not in other:
                                        other[gl_record_info] = []

                                    if record.is_indel:
                                        # Deletion
                                        if (linked_gl_base_info[2] == gl_record.REF and somatic_base_info[2] == 0 and len(record.REF) == 1):
                                            cis[gl_record_info].append((linked_gl_base_info[3],somatic_base_info[3]))
                                        elif (linked_gl_base_info[2] == gl_record.ALT[0] and somatic_base_info[2] > 0 and len(record.ALT[0]) > 1):
                                            cis[gl_record_info].append((linked_gl_base_info[3],somatic_base_info[3]))
                                        elif (linked_gl_base_info[2] == gl_record.REF and somatic_base_info[2] > 0 and len(record.ALT[0]) > 1):
                                            trans[gl_record_info].append((linked_gl_base_info[3],somatic_base_info[3]))
                                        elif (linked_gl_base_info[2] == gl_record.ALT[0] and somatic_base_info[2] == len(record.REF) == 1):
                                            trans[gl_record_info].append((linked_gl_base_info[3],somatic_base_info[3]))
                                        # Insertion
                                        elif (linked_gl_base_info[2] == gl_record.REF and somatic_base_info[2] == 0 and len(record.REF) > 1):
                                            cis[gl_record_info].append((linked_gl_base_info[3],somatic_base_info[3]))
                                        elif (linked_gl_base_info[2] == gl_record.ALT[0] and somatic_base_info[2] < 0 and len(record.ALT[0]) == 1):
                                            cis[gl_record_info].append((linked_gl_base_info[3],somatic_base_info[3]))
                                        elif (linked_gl_base_info[2] == gl_record.REF and somatic_base_info[2] < 0 and len(record.ALT[0]) == 1):
                                            trans[gl_record_info].append((linked_gl_base_info[3],somatic_base_info[3]))
                                        elif (linked_gl_base_info[2] == gl_record.ALT[0] and somatic_base_info[2] == len(record.REF) > 1):
                                            trans[gl_record_info].append((linked_gl_base_info[3],somatic_base_info[3]))
                                        else:
                                            other[gl_record_info].append((linked_gl_base_info[3],somatic_base_info[3]))
                                    else:
                                        if (linked_gl_base_info[2] == gl_record.REF and somatic_base_info[2] == record.REF):
                                            cis[gl_record_info].append((linked_gl_base_info[3],somatic_base_info[3]))
                                        elif (linked_gl_base_info[2] == gl_record.ALT[0] and somatic_base_info[2] == record.ALT[0]):
                                            cis[gl_record_info].append((linked_gl_base_info[3],somatic_base_info[3]))
                                        elif (linked_gl_base_info[2] == gl_record.REF and somatic_base_info[2] == record.ALT[0]):
                                            trans[gl_record_info].append((linked_gl_base_info[3],somatic_base_info[3]))
                                        elif (linked_gl_base_info[2] == gl_record.ALT[0] and somatic_base_info[2] == record.REF):
                                            trans[gl_record_info].append((linked_gl_base_info[3],somatic_base_info[3]))
                                        else:
                                            other[gl_record_info].append((linked_gl_base_info[3],somatic_base_info[3]))
                    scores = []
                    for gl_info in cis:
                        cis_qual = sum(map(sum,cis[gl_info]))
                        trans_qual = sum(map(sum,trans[gl_info]))
                        other_qual = sum(map(sum,other[gl_info]))
                        lplist = bayes_gt( cis_qual, trans_qual, other_qual )
                        score = generate_score2( lplist )
                        scores.append(score)
                    if len(scores) > 0:
                        if 'txt' in args.format:
                            contig_walker_txt_writer.write( "\t".join([sample_name, str(record), str(mean(scores)), str(scores), str(cis), str(trans), str(other) ])+"\n" )
                        if 'bed' in args.format:
                            contig_walker_bed_writer.write( "\t".join([str(record.CHROM), str(record.POS-1), str(record.POS), sample_name, str(mean(scores)) ])+"\n" )
                        if 'vcf' in args.format:
                            for call in record.samples:
                                if call.sample == sample_name:
                                    update_call_data(call, ['WS'], [mean(scores)], somatic_vcf_reader)
                if 'vcf' in args.format:
                    format_list = list(somatic_vcf_reader.formats.keys())
                    format_list.remove('GT')
                    format_list.insert(0,'GT')
                    # Add VAF information to the format field of each sample
                    record.FORMAT = ":".join(format_list)
                    contig_walker_vcf_writer.write_record(record)
            if 'txt' in args.format:
                contig_walker_txt_writer.close()
            if 'bed' in args.format:
                contig_walker_bed_writer.close()
            if 'vcf' in args.format:
                contig_walker_vcf_writer.close()
            q_out.put( contig )

        # Break the loop if the queue is empty
        except queue.Empty:
            break

def update_call_data( call, edit_keys, edit_values, vcf_reader ):
    """
    Function to add or update a field to the format field.
    This will be automatically update in the call object
    Input: A call object
    Input: A list with format fields
    Input: A list with format values
    """
    f_keys = list(vcf_reader.formats.keys())
    d = dict(call.data._asdict())
    f_vals = []
    for key in f_keys:
        if key in edit_keys:
            f_vals.append(edit_values[edit_keys.index(key)] )
        elif key in d:
            f_vals.append(d[key] )
        else:
            f_vals.append(None)
    handy_dict = dict(zip(f_keys, f_vals))
    f_keys.remove('GT')
    f_keys.insert(0,'GT')
    call.data = collections.namedtuple('CallData',f_keys)(**handy_dict)

def generate_score2( gt_lplist ):
    best, second_best, worst = sorted([ (i, e) for i, e in enumerate(gt_lplist) ], key=lambda x: x[1], reverse=True)[0:3]
    gt_sum = 0
    qual = 0
    for gt in gt_lplist:
        try:
            gt_sum += 10**gt
        except OverflowError:
            gt_sum += 0
    if gt_sum > 0:
        gt_sum_log = math.log(gt_sum, 10)
        qual = abs(-10 * (gt_lplist[2] - gt_sum_log)) # phred-scaled probability site is non-reference in this sample

    return( qual )

def log_choose(n, k):
        """
        Returns log of given variable.
        :param n:
        :param k:
        """
        r = 0.0
        # swap for efficiency if k is more than half of n
        if k * 2 > n:
            k = n - k
        for d in range(1, k + 1):
            r += math.log(n, 10)
            r -= math.log(d, 10)
            n -= 1

        return (r)

def bayes_gt(cis, trans, other):
        """
        Returns genotype list
        :param cis:
        :param trans:
        :return list with probabilities of different types:
        """
        prob = [1e-3, 0.5, 0.9]

        total = cis + trans + other
        log_combo = log_choose(int(total), int(trans))
        lp_cis = log_combo + trans * math.log(prob[0], 10) + cis * math.log(1 - prob[0], 10)
        lp_fp = log_combo + (trans * math.log(prob[1], 10)) + (cis * math.log(1 - prob[1], 10))
        lp_trans = log_combo + trans * math.log(prob[2], 10) + cis * math.log(1 - prob[2], 10)

        return [lp_cis, lp_trans, lp_fp]
        #return [lp_cis, lp_trans]

def get_linked_gl_base( gl_record, sample_name, somatic_align):
    gl_base_info = None
    for bam in args.bam:
        F=pysam.AlignmentFile(bam,'rb')
        sample_name2 = get_sample_name(F)

        if sample_name != sample_name2:
            continue
        for gl_pileupcolumn in F.pileup(gl_record.CHROM, int(gl_record.POS)-1, int(gl_record.POS), truncate=True, stepper='nofilter',min_base_quality=int(base_phred_quality)):
            for gl_pileupread in gl_pileupcolumn.pileups:
                if gl_pileupread.alignment.query_name != somatic_align.query_name:
                    continue
                if (gl_pileupread.alignment.is_read1 and somatic_align.is_read2) or (gl_pileupread.alignment.is_read2 and somatic_align.is_read1):
                    continue
                if not isinstance(gl_pileupread.query_position, int):
                    continue
                if (len(gl_record.REF) == 1 and len(gl_record.ALT[0]) == 1):
                    gl_base = gl_pileupread.alignment.query_sequence[gl_pileupread.query_position]
                    gl_base_qual = gl_pileupread.alignment.query_qualities[gl_pileupread.query_position]
                    gl_base_info = (gl_record.CHROM, gl_record.POS, gl_base, gl_base_qual)
    return( gl_base_info )

def get_somatic_alignments( record, sample_name ):
    """
    Function to get reads overlapping the somatic variant.
    Input: Record object
    Input: Sample name
    Output: list of tuples ( base, read )
    """
    alignments = []
    for bam in args.bam:
        F=pysam.AlignmentFile(bam,'rb')
        sample_name2 = get_sample_name(F)

        if sample_name != sample_name2:
            continue

        for pileupcolumn in F.pileup(record.CHROM, int(record.POS)-1, int(record.POS), truncate=True, stepper='nofilter',min_base_quality=int(base_phred_quality)):
            for pileupread in pileupcolumn.pileups:
                # QC the read
                if ( check_pileupread( pileupread) ):
                    # If variant is SNV
                    if (len(record.REF) == 1 and len(record.ALT[0]) == 1):
                        base = pileupread.alignment.query_sequence[pileupread.query_position]
                        base_qual = pileupread.alignment.query_qualities[pileupread.query_position]
                        alignments.append(((record.CHROM,record.POS,base,base_qual),pileupread.alignment))
                    # If variant is INDEL
                    else:
                        base = pileupread.indel
                        base_qual = pileupread.alignment.query_qualities[pileupread.query_position]
                        alignments.append(((record.CHROM,record.POS,base,base_qual),pileupread.alignment))

    return( alignments )

def default_to_regular(d):
    if isinstance(d, collections.defaultdict):
        d = {k: default_to_regular(v) for k, v in d.items()}
    return d

def get_linked_gl_records( germline_contig_vcf_reader, sample_name, contig_chr, contig_start, contig_end ):
    records = []
    try:
        # Try to parse the specific contig from the vcf
        germline_contig_vcf_reader.fetch(contig_chr, contig_start, contig_end)
    except:
        # Skip contig if it is not present in the vcf file
        return( records )
    for record in germline_contig_vcf_reader.fetch(contig_chr, contig_start, contig_end):
        # if ((record.ID and "COSM" not in record.ID) or ("ControlEvidence" in record.FILTER) ) and not record.is_indel:
        # for call in record.samples:
            # if call.sample == sample_name and call.is_het:
        if record.genotype(sample_name).is_het:
            records.append(record)
    return( records )

def check_pileupread( pileupread ):
    """
    Function to check a pileup read.
    Returns True if the read needs to be kept and returns False if read can be skipped.
    Input: Pileupread object
    Return: True or False
    """
    check = True
    if pileupread.alignment.is_duplicate:
        check = False
    elif pileupread.is_del:
        check = False
    elif pileupread.is_refskip:
        check = False
    elif not pileupread.query_position:
        check = False
    elif pileupread.alignment.mapq < int(mapq):
        check = False
    elif pileupread.alignment.query_qualities[pileupread.query_position] < int(base_phred_quality):
        check = False

    return( check )

def get_sample_name( bamfile ):
    """
    Function to get the sample name from the bam file
    Input: An AlignmentFile object of the bam file
    Return: The sample or False if there is no SM tag in the bam header
    """
    header = bamfile.header
    sample_name = False
    if 'RG' in header:
        if type(header['RG']) is list:
            sample_name = header['RG'][0]['SM']
        else:
            sample_name = header['RG']['SM']

    return( sample_name )

def fix_vcf_header( vcf_reader ):
    """
    Function to fix fields in the vcf header
    Input: A vcf reader object
    Return: The vcf reader object with fixed headers
    """
    #dbNSFP_clinvar_clnsig has a Integer type but sometimes it is a String, e.g. 2|2
    vcf_reader.infos['dbNSFP_clinvar_clnsig'] = pyvcf.parser._Info("dbNSFP_clinvar_clnsig",1,"String","Field 'clinvar_clnsig' from dbNSFP", None, None)
    #dbNSFP_clinvar_golden_stars has a Integer type but sometimes it is a String, e.g. 0|1
    vcf_reader.infos['dbNSFP_clinvar_golden_stars'] = pyvcf.parser._Info("dbNSFP_clinvar_golden_stars",1,"String","Field 'clinvar_golden_stars' from dbNSFP", None, None)
    vcf_reader.infos['dbNSFP_hg18_chr'] = pyvcf.parser._Info("dbNSFP_hg18_chr",1,"String","Field 'hg18_chr' from dbNSFP", None, None)
    vcf_reader.infos['dbNSFP_hg19_chr'] = pyvcf.parser._Info("dbNSFP_hg19_chr",1,"String","Field 'hg19_chr' from dbNSFP", None, None)
    return( vcf_reader )

def merge_tmp_files():
    """
    Function to merge all the tmp contig files
    """
    start = time.time()
    header = False
    # Loop through all chromomsomes
    for contig in contig_list:
        if 'vcf' in args.format:
            if not header:
                os.system('cat walker_tmp/{}.walker.vcf > {}.walker.vcf'.format(contig, somatic_vcf_name))
                header = True
            else:
                os.system('grep -v \'^#\' walker_tmp/{}.walker.vcf >> {}.walker.vcf'.format(contig, somatic_vcf_name))
        if 'txt' in args.format:
            os.system('cat walker_tmp/{}.walker.txt >> {}.walker.txt'.format(contig, somatic_vcf_name))
        if 'bed' in args.format:
            os.system('cat walker_tmp/{}.walker.bed >> {}.walker.bed'.format(contig, somatic_vcf_name))

    time.sleep(5)
#    os.system("rm -rf walker_tmp")
if __name__ == "__main__":
    main()
    merge_tmp_files()
