=================================
Genotyping By Sequencing Analysis
=================================


Introduction
============

GBS is a fast, cheap and relatively easy way of assaying genetic variation in a
very high throughput manner. This workshop will show attendees
the way we analyse GBS data in the Borevitz lab.

I'll go over both the theory of GBS and how we analyse it.

.. contents::


Who am I?
---------

::

    Kevin Murray
    PhD Student,
    Borevitz Lab,
    ARC Centre of Excellence in Plant Energy Biology,
    ANU,
    Canberra, Australia.

Feel free to contact me at kevin.murray@anu.edu.au about any part of this
workshop. There is also my `github <https://github.com/kdmurray91>`_ and
`personal website <https://kdmurray.id.au/>`_.



What is GBS?
============

GBS is a reduced-representation sequencing library preparation method. The
general idea is to assay a small, reproducible subset of the genomes of all
samples. This is achieved using a restriction enzyme digest followed by a
PCR-based Illumina library preparation protocol. Elshire et al. proposed the
method and their paper explains the protocol in detail [ElshireGBS]_. The
important points of the protocol are reproduced below.


GBS wet lab overview
--------------------

DNA is extracted as usual, and digested with a restriction enzyme (in our case,
PstI). Sticky-ended adaptors are ligated to digested fragments, and amplified
using PCR. This also roughly size-selects for inserts between about 50 and 500
bp. A gel- or bead-based clean up is used to refine this size selection.
Individual libraries are quantified and pooled equimolar into a single pooled
library that is sequenced in one Illumina lane (we use a HiSeq 2500, and have
had issues with the two-dye chemistry). For a more detail description of the
protocol, please see the paper describing the protocol [ElshireGBS]_.

<++> DIAGRAM FROM ELSHIRE PAPER

An important and common modification to the original protocol is the use of
combinatorial adaptors. This involves using modified adaptors such that the
forward and reverse reads contain independent short nucleotide sequences that
identify the sample. The Borevitz lab (and others) use barcodes of differing
length, which avoids nucleotide imbalance that would occur if all reads
contained the restriction site at the same position. Nucleotide imbalance
causes the Illumina base-calling algorithms to fail [some students may need
this idea explained / expanded]


GBS data overview
-----------------

The most common form of GBS data you are likely to receive is a raw FASTQ file.
This FASTQ file will contain all reads from all samples that were sequenced in
an Illumina lane. You will also need some form of metadata table associating at
least each sample's DNA barcode with sample name or ID, and likely also
population and collection metadata where applicable. Various tools have
different format requirements for this metadata, so some interconversion may be
required (when is it not in bioinformatics!).

Forgive the banality, but I should take this opportunity to emphasise the
importance of backups, replication and versioning of both sequencing data *and*
associated metadata. We've lost many thousands of dollars worth of GBS
experiment due to failures in hardware, software or user.



Analysis of GBS data
====================

Workflow overview
-----------------

- Assess dataset quality with FastQC_
- Demultiplex reads with Axe_
- QC reads with gbsqc_
- Detect loci and call variants *de novo* using Stacks_
- Plot PCoA of SNP matrix

[MRC: insert sample plot to motivate the efforts]


Data
----

Our data is derived from a sequencing experiment in *Eucalyptus melliodora*
(Yellow Box Gum), a key species in the Box Gum Grassy Woodland ecosystem that
surrounds Canberra. This ecosystem is endangered due to agricultural
deforestation throughout Australia's eastern seaboard, so several landscape
genomic projects are ongoing. Key to these projects is the determination of
genetic relatedness between many hundred samples taken from remnant stands of
*E. melliodora*.

You have been given several data files:

- ``Emel-lb1234_R1.fastq.gz`` and ``Emel-lb1234_R2.fastq.gz``: Raw read files
  (forward and reverse)
- ``Emel-lb1234.axe``: The Axe keyfile, a mapping of DNA barcodes to sample
  names.

TODO:

 - Add papers about biology


Metadata
--------

The information which any sequencing experiment generates is useless without
well curated metadata. This sounds self-evident, however in our experience most
issues that arise during the analysis of GBS data are caused but incorrect or
missing metadata. The sample names associated with our samples are available
`here </samples.axe>`_  (pre-formatted in the format Axe requires).


Quality Control
---------------

As is customary for all NGS analyses, the first step in the analysis of GBS
data is to check the technical quality of the reads we have obtained. This is
done with FastQC:

.. code-block:: shell

  mkdir -p fastqc
  fastqc --extract -o fastqc Emel-lb1234_R[12].fastq.gz

Inspect the FastQC HTML output (files under ``./fastqc/``).


Demulitplexing
--------------

You may remember our samples come in one big FASTQ file. This is not what we
want, so we need to demultiplex the reads such that the samples are each in
their own file. We do this before quality trimming, so that reads are not
manipulated before being demultiplexed (as barcode sequences often have quite
low quality scores).

Demultiplexing is performed using Axe, as few other demultiplexers can handle
the rather eclectic needs that GBS has. Barcodes differ in length, and are
applied combinatorially (different of R1 and R2). The following incantation
should to the trick:

.. code-block:: shell

  mkdir -p demuxed
  axe-demux                         \
        -c                          \
        -z 6                        \
        -b Emel-lb1234.axe          \
        -t Emel-lb1234.stats        \
        -f Emel-lb1234_R1.fastq.gz  \
        -r Emel-lb1234_R2.fastq.gz  \
        -I demuxed/

Axe will have demultiplexed reads into individual interleaved files, under the
directory ``./demuxed``. Sample-wise read counts have been saved to the
``Emel-lb1234.stats`` file.  The following R snippet can be used to generate a
histogram of read counts across all samples. You can run it on the command
line, or locally after downloading the stats file if you want to play around
with other plots or stats.

.. code-block:: R

  axe <- read.delim("Emel-lb1234.stats", stringsAsFactors=F)
  # Remove count of reads without barcodes
  axe <- axe[axe$Sample != "No Barcode",]
  hist(axe$ReadCounts)


Quality and Adaptor Trimming
----------------------------

We need to remove both adaptor read-through and low-quality sections from our
reads. Additionally, due to the rather inane requirement of Stacks that all
reads be the same length, we need to enforce the truncation of long reads, and
remove shorter reads. We use a tool of our own named gbsqc, but Trimmomatic and
other similar tools will work just as well (albeit with more duct-tape). As we
have many files now, we need to loop over each of them. Since we have multiple
cores to use, we can utilise GNU parallel instead of a simple for loop [#]_.

<++>Pre-prepare samples file

.. code-block:: shell

  cut -f 3 < Emel-lb1234.axe >Emel-lb1234.samples
  mkdir -p qcd reports
  cat Emel-lb1234.samples | parallel -j 4 --verbose \
    gbsqc -q 25                                     \
          -l 64                                     \
          -y reports/{}.yml                         \
          -y reports/{}.yml                         \
      \| gzip \> qcd/{}-qc_il.fastq.gz


So now we have a directory containing a FASTQ file for each sample. In theory,
no contaminants are present in the reads.

.. [#] In case you've never seen GNU parallel before, I urge you to look it up
   and become familiar with its use. It sure comes in handy.


Variant calling
---------------

Stacks is used to assemble loci and call variants in a *de novo* fashion.
Stacks works by clustering reads into loci, then detecting variation between

.. code-block:: shell

    # This is a hack to prepare a list of -s samp1.fq -s samp2.fq ...
    samples=$(for s in $(cat Emel-lb1234.samples);  \
              do echo "-s qcd/${i}-qc_il.fastq.gz"; \
              done)
    denovo_map.pl                                   \
        -T 11                                       \
        -t                                          \
        -S                                          \
        -b 1                                        \
        -n 2                                        \
        -o stacks_output                            \
        -s $samples

This command will create a population file, an internal data format that stacks
uses to represent its state. To produce a VCF file for further analysis, we use
the `populations` command from `stacks`.

.. code-block:: shell

    populations                                     \
        -t 11                                       \
        -r 0.25                                     \
        -p 4                                        \
        -b 1                                        \
        -P stacks_output                            \
        -M emel_lball.map                           \
        -e pstI                                     \
        --write_single_snp                          \
        --vcf                                       \
        --fstats



Pitfalls of GBS
===============

No protocol or method produces perfect data, and GBS certainly produces it's
share of imperfections. Throughout this section, keep in mind that GBS is not
designed as an absolute method able to definitively determine relatedness.
Rather GBS is a cheap, reliable estimate of relatedness. For many, if not most,
applications in population genetics, this is more than sufficient.


Technical batch effects
-----------------------

One artifact we sometimes see is artifacts of the library preparation protocol.
In particular, we have seen cases where there is a strong lane effect on
genetic signal. This was traced to inconsistent size selection. Also keep in
mind that GBS relies o<+FINISH THIS SENTENCE+>


Input sample quality
--------------------

Input DNA quality can have a significant effect on the quality of results.
Partially degraded DNA will form libraries of poor quality or low complexity,
and can lead to systematic effects if sample quality is confounded with
biologically significant variables (which it often is).


Sample heterogeneity
--------------------

Frequently the concentration of DNA in individual libraries is too low to
reliably quantitiate. This can lead to quite variable coverage between samples,
that in turn can cause inaccuracies in the calculation of relatedness. The best
course of action in such a situation is simply to drop samples with too few
reads. The exact definition of "too few" is debatable, but we frequently use
500,000 reads as a hard cut off, and sometimes raise this to 1 million. Any
other choice is probably equally valid and equally arbitrary.

It is worth bearing this advice in mind as early as possible in the planning of
GBS experiments. GBS is a high throughput method, and samples fail at greater
frequency than other methods. If you have samples that are particularly
important, please consider sequencing them in at least duplicate. This is
especially true if your important samples are of lower quality (which they
often are).


Metadata mix-ups
----------------

This is not at all GBS specific, but as previously mentioned metadata is key to
the interpretation of any GBS dataset.  <++FINISH THIS+>


Contamination
-------------

As is the case for most *de novo* algorithms, there is an implicit assumption
that all reads come from the same individual. However biology can sometimes get
in the road of this reasonable assumption, particularly in plant species with
endo- or exophytic microorgansims. We have seen cases where up to 20% of reads
and a similar percentage of assembled loci come from fungal or bacterial
endosymbionts of *Eucalyptus*. This is not limited to plant species, there are
many organisms with similar microorgansimal communities.

If your samples are know or suspected to contain genetic material from other
species, it may be worth using taxonomic read classification tools such as
Kraken to partition reads into target and non-target species after QC, and
proceed with loci assembly and variant calling only with target species reads.
An alternative is to use BLAST or similar tools to taxonomically classify the
assembled loci, and exclude any non-target species' loci from the VCF file
before any post-analysis.



The Teacher's Pet Section
=========================

If you've managed to blaze through all the above, or are super-bored on the
way home, here are some extra things to try.

<+FILL IN OR REMOVE THIS SECTION+>

References
==========

Papers
------

.. [ElshireGBS]  Elshire RJ et al. (2011) **A Robust, Simple
    Genotyping-by-Sequencing (GBS) Approach for High Diversity Species.** *PLoS
    ONE* doi:`10.1371/journal.pone.0019379
    <https://dx.doi.org/10.1371/journal.pone.0019379>`_

Tools
-----

.. _FastQC: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
.. _Axe: https://github.com/kdmurray91/axe
.. _gbsqc: https://github.com/kdmurray91/libqcpp
.. _Stacks: http://catchenlab.life.illinois.edu/stacks/
