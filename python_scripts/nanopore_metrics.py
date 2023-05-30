import argparse
import pysam
import math
from statistics import median

def phred_to_prob(phred_score):
    return 10 ** (-phred_score / 10)

def prob_to_phred(prob):
    return -10 * math.log10(prob)

parser = argparse.ArgumentParser(description='Calculate various metrics for a FASTQ file')
parser.add_argument('fastq_file', help='input fastq file')
parser.add_argument('--log_file', default='metrics.log', help='output log file')
args = parser.parse_args()

num_reads = 0
total_accuracy = 0
total_quality = 0
total_read_length = 0
shortest_read = float('inf')
longest_read = 0
total_gc_content = 0
total_n_content = 0
total_error_prob = 0
total_qscore = 0

with pysam.FastxFile(args.fastq_file) as f, open(args.log_file, 'w') as log_file:
    for entry in f:
        read_length = len(entry.sequence)
        total_read_length += read_length
        shortest_read = min(shortest_read, read_length)
        longest_read = max(longest_read, read_length)
        gc_content = entry.sequence.count('G') + entry.sequence.count('C')
        total_gc_content += gc_content
        n_content = entry.sequence.count('N')
        total_n_content += n_content
        qualities = entry.get_quality_array()
        total_error_prob += sum(phred_to_prob(q) for q in qualities)
        accuracy = sum(1 for base in entry.sequence if base != 'N') / read_length
        total_accuracy += accuracy
        total_quality += sum(qualities)
        num_reads += 1
        error_probs = [phred_to_prob(q) for q in qualities]
        mean_error_prob = sum(error_probs) / len(error_probs)
        mean_qscore = prob_to_phred(mean_error_prob)
        total_qscore += mean_qscore

    avg_qscore = total_qscore / num_reads
    avg_error_prob = phred_to_prob(avg_qscore)
    avg_accuracy_adjusted = 1 - avg_error_prob
    mean_read_length = total_read_length / num_reads
    avg_accuracy = total_accuracy / num_reads
    avg_quality = total_quality / (num_reads * mean_read_length)
    avg_gc_content = total_gc_content / total_read_length
    avg_n_content = total_n_content / total_read_length
    mean_error_prob = total_error_prob / (num_reads * mean_read_length)
    mean_quality_score_adjusted = prob_to_phred(mean_error_prob)

with open(args.log_file, 'w') as log_file:
    log_file.write(f"Number of reads: {num_reads}\n")
    log_file.write(f"Mean read length: {mean_read_length:.2f}\n")
    log_file.write(f"Shortest read: {shortest_read}\n")
    log_file.write(f"Longest read: {longest_read}\n")
    log_file.write(f"Average read accuracy: {avg_accuracy:.2%}\n")
    log_file.write(f"Adjusted average read accuracy: {avg_accuracy_adjusted:.2%}\n")
    log_file.write(f"Average read quality score: {avg_quality:.2f}\n")
    log_file.write(f"Adjusted average read quality score: {mean_quality_score_adjusted:.2f}\n")
    log_file.write(f"Average GC content: {avg_gc_content:.2%}\n")
    log_file.write(f"Average N content: {avg_n_content:.2%}\n")

