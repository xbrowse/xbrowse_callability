import argparse
import sausage


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

    # worst possible implementation
    print(gene)
    print(qmi_results)
    print(dt_results)

    return



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