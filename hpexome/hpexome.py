import os
import os.path
import shutil
import subprocess
import tarfile
import tempfile
import urllib.request
from os import listdir
from os.path import join, basename, isdir, isfile, exists, abspath
from re import compile

import click
import pkg_resources


def check_files_exist(files):
    """Check a list of files printing those that not found"""
    files_not_found = [file for file in files if not isfile(file)]
    if any(files_not_found):
        for file in files_not_found:
            click.echo('File not found: ' + file, err=True)
        exit(1)

@click.command()
@click.option('-I', '--bam', 'bams', required=True, multiple=True, help='One or more BAM files or directories')
@click.option('-R', '--genome', 'genome_fasta_file', required=True, help='Reference genome in single FASTA file')
@click.option('--dbsnp', 'dbsnp_file', required=True, help='dbSNP file in VCF format')
@click.option('--sites', 'known_sites_files', required=True, multiple=True,
              help='VCF files containing known polymorphic sites to skip over in the recalibration algorithm')
@click.option('--indels', 'known_indels_files', multiple=True,
              help='Inputs the VCF file with known indels to be used')
@click.option('-L', '--intervals', 'intervals_files', multiple=True,
              help='One or more genomic intervals over which to operate')
@click.option('--unified_vcf', is_flag=True, default=False, help='Unify VCF files into a single one', show_default=True)
@click.option('-O', '--output_file_name', default='unified.vcf', help='Output file name for unified VCF',
              show_default=True)
@click.option('--min_prunning', default=2, help='Minimum support to not prune paths in the graph', show_default=True)
@click.option('--stand_call_conf', default=30, show_default=True,
              help='Minimum phred-scaled confidence threshold at which variants should be called')
@click.option('-nt', '--num_data_threads', type=click.INT, help='Number of data threads')
@click.option('-nct', '--num_threads_per_data_thread', type=click.INT, help='Number of threads per data thread')
@click.option('--scatter_count', type=click.INT)
@click.option('--job_runner', help='Use the specified job runner to dispatch command line jobs')
@click.option('--job_queue', help='Default queue for compute farm jobs')
@click.option('--job_native', multiple=True, help='Native arguments to pass to the job runner')
@click.option('--logging_level', help='Set the minimum level of logging')
@click.option('--dont_run', is_flag=True, default=False, show_default=True, help='Perform dry run')
@click.option('--java_path', default='java', help='Path to java', show_default=True)
@click.option('--java_mem', help='Maximum Java memory in GB.')
@click.argument('destination', default='.', type=click.Path())
def hpexome(bams, genome_fasta_file, dbsnp_file,
            known_sites_files, known_indels_files, intervals_files,
            unified_vcf, output_file_name, min_prunning, stand_call_conf,
            num_data_threads, num_threads_per_data_thread, scatter_count,
            job_runner, job_queue, job_native, logging_level, dont_run,
            java_path, java_mem, destination):
    """An automated workflow for processing whole-exome sequencing data"""
    # given a list of BAM files or directories containing them,
    # create a list of absolute path to BAM files
    m = compile('\\.bam$')
    bam_files = []
    for bam in bams:
        if isdir(bam):
            files = listdir(bam)
            bam_files.extend([abspath(join(bam, file)) for file in files if m.search(basename(file))])
        elif isfile(bam):
            bam_files.append(abspath(bam))
        else:
            click.echo('File or directory not found: ' + bam, err=True)
            exit(1)

    # get Queue jar file and script inside this Python package
    # start building command line
    queue_path = pkg_resources.resource_filename(__name__, 'Queue.jar')
    script_path = pkg_resources.resource_filename(__name__, 'HPexome.scala')
    command = [java_path, '-jar', queue_path, '-S', script_path]

    if not dont_run:
        command.append('-run')

    arguments = {
        '-I': bam_files, 
        '-R': abspath(genome_fasta_file), 
        '-dbsnp': abspath(dbsnp_file), 
        '-known': [abspath(f) for f in known_indels_files],
        '-knownSites': [abspath(f) for f in known_sites_files],
        '-L': [abspath(f) for f in intervals_files]
    }
    for argument, value in arguments.items():
        if value:
            files = [value] if isinstance(value, str) else value
            check_files_exist(files)
            for file in files:
                command.extend([argument, file])
    if unified_vcf:
        command.extend(['-unifiedVCF', '-filename', output_file_name])

    command.extend(['-stand_call_conf', str(stand_call_conf), '-minPruning', str(min_prunning)])

    if num_data_threads:
        command.extend(['-nt', str(num_data_threads)])
    if num_threads_per_data_thread:
        command.extend(['-nct', str(num_threads_per_data_thread)])
    if scatter_count:
        command.extend(['-scattercount', str(scatter_count)])
    if job_runner:
        command.extend(['-jobRunner', job_runner])
    if job_queue:
        command.extend(['-jobQueue', job_queue])
    for value in job_native:
        command.extend(['-jobNative', value])
    if logging_level:
        command.extend(['-l', logging_level])
    if java_mem:
        command.extend(['-Xmx', java_mem])

    if not exists(destination):
        os.mkdir(destination)

    click.echo('Executing command: ' + ' '.join(command), err=True)
    exit(subprocess.call(command, cwd=destination))
