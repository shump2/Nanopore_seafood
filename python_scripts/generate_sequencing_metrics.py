import argparse
import csv
import pysam
import glob
import os
import math

def phred_to_prob(phred_score):
    return 10 ** (-phred_score / 10)

def prob_to_phred(prob):
    return -10 * math.log10(prob)

def read_accuracy(quality_array):
    error_probs = [phred_to_prob(q) for q in quality_array]
    mean_error_prob = sum(error_probs) / len(error_probs)
    accuracy = 1 - mean_error_prob
    return accuracy

def ave_qual(quals):
    return -10 * math.log10(sum([10**(q/-10) for q in quals]) / len(quals))

parser = argparse.ArgumentParser(description='Generate input data file from multiple FASTQ files in a directory')
parser.add_argument('-d', '--dir', help='input directory containing FASTQ files', required=True)
parser.add_argument('-o', '--output', help='output CSV file', required=True)
args = parser.parse_args()

fastq_files = glob.glob(os.path.join(args.dir, '*.fastq'))

with open(args.output, 'w', newline='') as csvfile:
    fieldnames = ['sample_id', 'read_id', 'read_length', 'base_quality_score', 'read_accuracy']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()

    for fastq_file in fastq_files:
        sample_id = os.path.splitext(os.path.basename(fastq_file))[0] # Extract the basename and remove the extension
        with pysam.FastxFile(fastq_file) as f:
            for entry in f:
                read_id = entry.name
                read_length = len(entry.sequence)
                quality_array = entry.get_quality_array()
                base_quality_score = ave_qual(quality_array)
                accuracy = read_accuracy(quality_array)

                writer.writerow({
                    'sample_id': sample_id,
                    'read_id': read_id,
                    'read_length': read_length,
                    'base_quality_score': base_quality_score,
                    'read_accuracy': accuracy
                })

