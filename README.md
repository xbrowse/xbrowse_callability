xBrowse Callability
===================

This is Andrew and Brett's sandbox for exploring callability.
Briefly, it provides one method - `display_gene_callability` -
that prints out (ASCII) an intuitive summary for the callability results from one gene in one sample.

Some background: "Callability" broadly encompasses trying to answer the question
*If there were a variant in this gene, would we have seen it?*
We've had a couple takes at trying to set up an interface for this in xbrowse (eg. "Exome Coverage") but they were unsuccessful.

We want to be able to summarize a list of 100 genes for a clinician, intuitively:
*These 80 genes are fine. These 20 genes have potential issues. Here's why..*

The problem is that we've struggled with exactly what to present.
So, this repository was set up in order to support rapid iteration.
This one method has access to all of the necessary data.
Challenge for the analyst is to summarize everything in an intuitive way.

Once we settle on something, we'll implement in a pretty graphic in xbrowse

## Install

First install the requirements, preferably in a virtualenv:

    pip install -r requirements.txt

Then copy the callability output - all the output files from QualifyMissingIntervals and DiagnoseTargets -
into the `data/` directory

Then you can run with:

    python gene_callability -h

Ready set go!