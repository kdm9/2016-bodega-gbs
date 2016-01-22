=================================
Genotyping by Sequencing Analysis
=================================



Introduction
============

GBS is a fast, cheap and relatively easy way of assaying genetic variation in a
very high throughput manner. This workshop will go hopefully show attendees
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
PstI). 


GBS data overview
-----------------

The most common form of GBS data you are likely to receive is a raw FASTQ file.
This FASTQ file will contain a



Analysis of GBS data
====================

Workflow overview
-----------------

- Assess dataset quality with FastQC_
- Demultiplex reads with Axe_
- QC reads with gbsqc_
- Detect loci and call variants *de novo* using Stacks_
- Plot PCoA of SNP matrix

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

As is customary, the first step in the analysis of GBS data is to check the
technical quality of the reads we have obtained. This is done with FastQC: ::

  fastqc -o Emel-lb1234 Emel-lb1234_R[12].fastq.gz

Inspect the FastQC HTML output (files under ``./Emel-lb1234/``).

Demulitplexing
--------------

You may remember our samples come in one big FASTQ file. This is not what we
want, so we need to demultiplex the reads such that the samples are each in
their own file. We do this before quality trimming, so that reads are not
manipulated before being demultiplexed (as barcode sequences often have quite
low read numbers).

Demultiplexing is performed using Axe, as few other demultiplexers can handle
the rather eclectic needs that GBS has. Barcodes differ in length, and are
applied combinatorially (different of R1 and R2). The following incantation
should to the trick: ::

  mkdir -p demuxed
  axe-demux                         \
        -c                          \
        -z                          \
        -b Emel-lb1234.axe          \
        -t Emel-lb1234.stats        \
        -f Emel-lb1234_R1.fastq.gz  \
        -r Emel-lb1234_R2.fastq.gz  \
        -I demuxed

Axe will have demultiplexed reads into individual interleaved files, under the
directory ``./demuxed``. Sample-wise read counts have been saved to the
``Emel-lb1234.stats`` file.

::

  axe <- read.delim("tab", stringsAsFactors=F)
  axe <- axe[axe$Sample != "No Barcode",]


Quality and Adaptor Trimming
----------------------------

We need to remove both adaptor read-through and low-quality sections from our
reads. Additionally, due to the rather inane requirement of Stacks that all
reads be the same length, we need to enforce the truncation of long reads, and
remove shorter reads. We use a tool of our own named gbsqc, but Trimmomatic and
other similar tools will work just as well. As we have many files now, we need
to loop over each of them.  Since we have multiple cores to use, we can utilise
GNU parallel instead of a simple for loop.

::

  cut -f 3 < Emel-lb1234.axe >Emel-lb1234.samples
  mkdir -p qcd report
  cat Emel-lb1234.samples | parallel -j 4 --verbose \
    gbsqc -q 25                                     \
          -l 64                                     \
          -y reports/{}.yml                         \
          -y reports/{}.yml                         \
      \| gzip \> qcd/{}-qc_il.fastq.gz


So now we have a directory containing a FASTQ file for each sample. In theory,
no contaminants are present in the reads.


TODO:

- Add GNU parallel footnote


Hocus Pocus
-----------





The Teacher's Pet Section
=========================

If you've managed to blaze through all the above, or are super-bored on the
way home, here are some extra things to try.



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
