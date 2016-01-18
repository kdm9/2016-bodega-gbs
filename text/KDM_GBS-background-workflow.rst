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

Kevin Murray
PhD Student, Borevitz Lab, CPEB, ANU, Canberra, Australia.

What is GBS?
============

GBS is a reduced-representation sequencing library preparation method. The
general idea is to assay a small, reproducible subset of the genomes of all
samples. This is achieved using a restriction enzyme digest followed by a
PCR-based illumina libary preparation protocol. Elshire et al. proposed the
method and their paper explains the protocol in detail [ElshireGBS]_. The
important points of the protocol are reproduced below.

GBS wet lab overview
--------------------



GBS data overview
-----------------



Analysis of GBS data
====================

Workflow overview
-----------------


- Assess dataset quality with FastQC_
- Demultiplex reads with Axe_
- QC reads with gbsqc_
- Detect loci and call variants *de novo* using Stacks_


Metadata
--------

The information which any sequencing experiment generates is useless without
well curated metadata associating  

Quality Control
---------------


The Teachers Pet Section
========================

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
