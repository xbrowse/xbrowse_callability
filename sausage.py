import argparse
import json
import requests
import os
import pandas as pd


DATA_DIR = os.path.join(os.path.dirname(__file__), 'data/')


def get_num_comment_lines(filepath):
    """
    How many header lines are there at the top of the file that we should ignore?
    There is probably a cleaner way to do this with file-like objects but escaping me right now
    """
    i = 0
    for line in open(filepath):
        if line.startswith('#') or line.strip() == '':
            i += 1
        else:
            return i


def get_gene(gene_id):
    """
    Get the gene structure from the ensembl REST API
    Includes gene, transcript, exons, and cds
    Stored as nested dictionaries
    """
    gene_request = requests.get(
        'http://beta.rest.ensembl.org/feature/id/' + gene_id,
        params={'feature': 'gene', 'content-type': 'application/json'}
    )
    gene = gene_request.json()[0]

    transcripts_request = requests.get(
        'http://beta.rest.ensembl.org/feature/id/' + gene_id,
        params={'feature': 'transcript', 'content-type': 'application/json'}
    )
    gene['transcripts'] = transcripts_request.json()

    exons_request = requests.get(
        'http://beta.rest.ensembl.org/feature/id/' + gene_id,
        params={'feature': 'exon', 'content-type': 'application/json'}
    )
    gene['exons'] = exons_request.json()

    cds_request = requests.get(
        'http://beta.rest.ensembl.org/feature/id/' + gene_id,
        params={'feature': 'cds', 'content-type': 'application/json'}
    )
    gene['cds'] = cds_request.json()

    return gene


def get_qmi_output_for_region(sample_id, chrom, start, stop):
    """
    QMI = Qualify Missing Intervals
    get the QMI output for one (contiguous) region in one individual

    **IMPORTANT**: everything is 1-based (vcf) instead of 0-based (bed)
    """
    qmi_file_path = DATA_DIR + sample_id + '.grp'
    headers = None
    ret = []
    for line in open(qmi_file_path):
        if line.startswith('#') or line.strip() == '':
            continue
        elif line.startswith('INTERVAL'):
            headers = line.strip('\n').split()
            continue

        # this is a data line
        _chrom, pos_str = line.split(None, 1)[0].split(':')
        if '-' in pos_str:
            _start, _stop = map(int, pos_str.split('-'))
        else:
            _start = _stop = int(pos_str)

        # out of region
        if _chrom != chrom or stop < _start or start > _stop:
            continue

        line = line.replace(', ', ',')
        fields = line.strip('\n').split()
        ret.append(dict(zip(headers, fields)))
    return ret


    # this proved useless because of awkward file parsing
    # num_comment_lines = get_num_comment_lines(qmi_file_path)
    # a = pd.read_csv(open(qmi_file_path), skiprows=num_comment_lines, sep='[^,]\s*', header=num_comment_lines-1)
    # return a


def get_diagnose_targets_output(sample_id, chrom, start, stop):
    """
    See above - get all the diagnosetargets outputs that overlap this region
    List of dictionaries
    """
    dt_file_path = DATA_DIR + sample_id + '.vcf'
    num_comment_lines = get_num_comment_lines(dt_file_path)
    targets = pd.read_csv(open(dt_file_path), skiprows=num_comment_lines-1, sep='\t')
    targets['#CHROM'] = targets['#CHROM'].astype(str)
    region_mask = (targets['#CHROM'] == chrom) & (start <= targets['POS']) & (stop >= targets['POS'])
    return targets[region_mask].T.to_dict().values()

