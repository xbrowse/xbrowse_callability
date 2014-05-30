import argparse
import sausage
from pprint import pprint
import re
from pybedtools import BedTool
import pandas as pd


def display_gene_callability(
        gene,
        qmi_results,
        dt_results,
):
    """

    This function should do nothing but print things
    The goal is to see if we can print an ascii description of the callability in this gene,
    simplified to the expertise of a clinician
    So, conditions will be important. We don't want to print every QMI output, nobody is going to read that.
    We want to just summarize the relevant info - but "relevant" is for us to decide as a group

    Args:
        gene (dict): information about the gene, as a nested dictionary.
            There are keys for 'exons', 'transcripts', and 'cds'
        qmi_results (list): list of Qualify Missing Intervals outputs (each represented as a dict) that overlap this gene
            Note that "overlap" means it falls within the gene boundaries - it could be in an intron I suppose
        dt_results (list): like above, but output of DiagnoseTargets.
            We *think* that these entries correspond 1-1 to the exome capture targets, but this must be checked!
    Returns:
        Nothing! Just print a pretty summary

    """
    print '--------------------------------------'
    print 'Report for Sample:', args.sample_id
    print 'Gene Name: %s (%s)' % (gene['ID'], gene['external_name'])
    print '--------------------------------------'
    print ''

    missing_intervals_and_overlapping_features = get_features_with_missing_intervals(gene['cds'], qmi_results)

    not_baited_count = 0
    unknown_count = 0
    gc_content_count = 0

    unknown_bases_count = 0
    gc_content_bases_count = 0
    not_baited_count = 0

    total_gene = pd.DataFrame.from_records(gene['cds'])
    total_gene = total_gene.drop_duplicates()

    total_cds_bases = (total_gene['end'] - total_gene['start']).sum()

    missing_intervals = []
    for x in missing_intervals_and_overlapping_features:
        missing_intervals += x[1]

    df = pd.DataFrame.from_records(missing_intervals)
    df = df.drop_duplicates()
    total_bases = (df['MISSING_SIZE'].astype(int)).sum()

    not_baited = df[df['BAITED'] == 'false']
    not_baited_count = not_baited['MISSING_SIZE'].astype(int).sum()

    high_gc_content = df[(df['INTERPRETATION'] == 'GCCONTENT') & (df['BAITED'] == 'true')]
    high_gc_content_count = high_gc_content['MISSING_SIZE'].astype(int).sum()

    unknown = df[(df['INTERPRETATION'] == 'UNKNOWN') & (df['BAITED'] == 'true')]
    unknown_count = unknown['MISSING_SIZE'].astype(int).sum()

    print '%d bases within CDS features in %s were found to be missing' % (total_bases, gene['external_name'])
    print '---> Not baited: %d' % (not_baited_count)

    for interval in not_baited['INTERVAL'].tolist():
        print '    ---> %s' % (interval)
    print '---> High GC content: %d' % (high_gc_content_count)
    print '---> Unknown cause: %d' % (unknown_count)


def get_features_with_missing_intervals(features, qmi_results):
    """
    Return a list of tuples, each of which contains a feature (each in dictionary format) and a list of missing intervals
    (each in dictionary format) that overlap said feature.
    """
    interval_list = _get_missing_interval_list(qmi_results)
    qmi_bed_tool = BedTool(interval_list)

    results = []
    for feature in features:

        feature_bed_tool = BedTool([(u'chr' + feature['seq_region_name'], int(feature['start']), int(feature['end']))])
        overlapping_missing_intervals = list(qmi_bed_tool.intersect(feature_bed_tool))

        if len(overlapping_missing_intervals) > 0:
            missing_interval_list = []
            for x in overlapping_missing_intervals:
                # Get the relevant feature dict
                missing_interval = qmi_results[int(x[3])]

                # Adjust the missing interval to only the missing portion
                missing_interval['INTERVAL'] = '%s:%s-%s' % (x[0], x[1], x[2])
                missing_interval['MISSING_SIZE'] = int(x[2]) - int(x[1])
                missing_interval_list.append(missing_interval)

            results.append((feature, missing_interval_list))

    return results


def _get_missing_interval_list(qmi_results):
    """
    Return a list the missing intervals in bed format (chromosome, start, end, index in list of features)
    Index in list used pull full interval definition back out after overlap calculation.
    """
    interval_list = []
    for i, interval in enumerate(qmi_results):
        qmi_chromosome, qmi_start, qmi_end = _parse_interval(interval['INTERVAL'])
        interval_list.append(('chr' + qmi_chromosome, qmi_start, qmi_end, i))

    return interval_list


def _parse_interval(interval):
    """
    Return chromsome name, start coordinate, and end coordinate of a missing interval from qmi_results
    interval notation (chromosome:start-end)
    """

    match = re.search(r'([\dA-Za-z]+):(\d+)-?(\d+)?', interval)

    if not match:
        raise ValueError('Unexpected genomic interval format: ' + str(interval))

    chromosome_name = match.group(1)
    start_coordinate = match.group(2)
    end_coordinate = match.group(3)

    # One-base intervals will not have an end coordinate
    if not end_coordinate:
        end_coordinate = start_coordinate

    return chromosome_name, int(start_coordinate), int(end_coordinate)


def run(gene_id, sample_id):

    gene = sausage.get_gene(gene_id)
    qmi_results = sausage.get_qmi_output_for_region(
        sample_id,
        gene['seq_region_name'],
        gene['start'],
        gene['end']
    )
    dt_results = sausage.get_diagnose_targets_output(
        sample_id,
        gene['seq_region_name'],
        gene['start'],
        gene['end']
    )
    display_gene_callability(gene, qmi_results, dt_results)


if __name__ == "__main__":

    # feel free to edit command line interface - this was just quick
    parser = argparse.ArgumentParser('Display callability summary for one gene in one sample')
    parser.add_argument('gene_id', help='Ensembl Gene ID.')  # TODO: switch this to HGNC name
    parser.add_argument('sample_id')
    args = parser.parse_args()

    # no type checking or anything - should add some
    run(args.gene_id, args.sample_id)
