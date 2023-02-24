Q2_ITSxpress: A Qiime2 plugin to rapidly trim the Internally transcribed spacer (ITS) region of FASTQ files
===========================================================================================================
.. image:: https://travis-ci.org/USDA-ARS-GBRU/q2_itsxpress.svg?branch=master
  :target: https://travis-ci.org/USDA-ARS-GBRU/q2_itsxpress

.. image:: https://codecov.io/gh/USDA-ARS-GBRU/q2_itsxpress/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/USDA-ARS-GBRU/q2_itsxpress

.. image:: https://api.codacy.com/project/badge/Grade/4d00341b4abc4e04a77cf5ca6674cd3c
  :target: https://www.codacy.com/app/USDA-ARS-GBRU/q2_itsxpress?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=USDA-ARS-GBRU/q2_itsxpress&amp;utm_campaign=Badge_Grade

.. image:: https://zenodo.org/badge/138209572.svg
   :target: https://zenodo.org/badge/latestdoi/138209572



#####

**This is the end of life version 1 of q2_itsxpress and the command line version of ITSxpress.
The new version 2 of ITSxpress, has the Qiime2 plugin built in with command line version of ITSxpress. See 
master branch of `ITSxpress`_.**

#####

Authors
-------
* Adam R. Rivers, US Department of Agriculture, Agricultural Research Service
* Kyle C. Weber, US Department of Agriculture, Agricultural Research Service
* Sveinn V. Einarsson, US Department of Agriculture, Agricultural Research Service

Citation
--------
Rivers AR, Weber KC, Gardner TG et al. ITSxpress: Software to rapidly trim
internally transcribed spacer sequences with quality scores for marker gene
analysis. F1000Research 2018, 7:1418. doi: `10.12688/f1000research.15704.1`_

.. _`10.12688/f1000research.15704.1`: https://doi.org/10.12688/f1000research.15704.1

Introduction
------------

The internally transcribed spacer (ITS) is a region between the small subunit
and large subunit rRNA genes. In is a commonly used phylogenetic marker for
Fungi and other Eukaryotes. The ITS contains the 5.8s gene and two variable
length spacer regions. In amplicon sequencing studies it is common practice to
trim off the conserved (SSU, 5,8S or LSU) regions. Bengtsson-Palme et al. (2013)
published a software package ITSx_ to do this.

Q2-ITSxpress extends this work by rapidly trimming FASTQ sequences within
Qiime2.  Q2-ITSxpress is the Qiime2 plugin version of the stand alone command
line utility ITSxpress_. Q2_ITSxpress is designed to support the calling of
exact sequence variants rather than OTUs. This newer method of sequence
error-correction requires quality score data from each sequence, so each input
sequence must be trimmed. ITSxpress makes this possible by taking FASTQ data,
de-replicating the sequences then identifying the start and stop sites using
HMMSearch. Results are parsed and the trimmed files are returned. The ITS1,
ITS2 or the entire ITS region including the 5.8s rRNA gene can be selected.
ITSxpress uses the hmm models from ITSx so results are nearly identical.


Requirements/Dependencies
-------------------------

* Qiime2 is required to run Q2-itsxpress (for stand alone software see ITSxpress_)
* To install Qiime2 follow these instructions: https://docs.qiime2.org/2022.8/install/
* This end of life version 1 of q2-itsxpress and ITSxpress is **ONLY** compatible with Qiime2 version 2022.8. So make sure to follow the link above.

* We are using mamba because it resolves packages better and faster, but conda can be substituted.
	- Information on installing mamba or micromamba (either highly recommended) can be found here: `mamba installation guide`_

.. _`mamba installation guide`: https://mamba.readthedocs.io/en/latest/installation.html

Q2-itsxpress plugin installation
--------------------------------
1. Example on how to install and create new Qiime2-2022.8 environment.

.. code-block:: bash

  wget https://data.qiime2.org/distro/core/qiime2-2022.8-py38-osx-conda.yml
  mamba env create -n qiime2-2022.8 --file qiime2-2022.8-py38-osx-conda.yml

2. Activate the Qiime2 conda environment

.. code-block:: bash

  mamba activate qiime2-2022.8

3. Install Q2_itsxpress using BioConda_. Be sure to install itsxpress and q2_itsxpress in the Qiime2 environment using the following commands.

.. code-block:: bash

  mamba install -c bioconda itsxpress==1.8.1
  pip install q2-itsxpress

4. In your Qiime2 environment, refresh the plugins.

.. code-block:: bash

  qiime dev refresh-cache

5. Check to see if the ITSxpress plugin is installed. You should see an output similar to the image below.

.. code-block:: bash

  qiime itsxpress

.. image:: ../../screenshot.png

Usage
-----

Within Qiime2 you can trim paired-end or single-end reads using these commands

.. code-block:: bash

  qiime itsxpress trim-pair

  qiime itsxpress trim-pair-output-unmerged

  qiime itsxpress trim-single

1. qiime itsxpress trim-single

  This command takes single-end data and returns trimmed reads. The sequence may
  have been merged previously or have been generated from a long read technology
  like PacBio. Merged and long reads trimmed by this function can be used by
  Deblur but only long reads (not merged reads) trimmed by this function should
  be passed to Dada2. Its statistical model for estimating error rates was not
  designed for pre-merged reads.

+----------------------------------+---------------------------------------------------------------------------------------+
|    Command-requirement           | Description                                                                           |
+----------------------------------+---------------------------------------------------------------------------------------+
|   --i-per-sample-sequences       | - The artifact that contains the sequence file(s).                                    |
+ 			           + - Either Joined Paired or just a single fastq.                                        +
|                                  | - One file sequence in the qza data folder.                                           |
+----------------------------------+---------------------------------------------------------------------------------------+
|       --p-region                 | - The regions ITS2, ITS1, and ALL.                                                    |
+----------------------------------+---------------------------------------------------------------------------------------+
|				   | -	Select the taxonomic group sequenced: A, B, C, D, E, F, G, H, I, L, M, O, P,	   |
+	--p-taxa		   +	Q, R, S, T, U, V, ALL.								   +
| 				   |											   |
+----------------------------------+---------------------------------------------------------------------------------------+
|       --p-threads 	           | - The amount of threads to use.                                                       |
+----------------------------------+---------------------------------------------------------------------------------------+
|       --o-trimmed                | - The resulting trimmed sequences from ITSxpress in a qza format.                     |
+----------------------------------+---------------------------------------------------------------------------------------+
|      --cluster-id                | - The percent identity for clustering reads, set to 1 for exact dereplication.        |
+----------------------------------+---------------------------------------------------------------------------------------+


2. qiime itsxpress trim-pair

  This command takes paired-end data and returns merged, trimmed reads. The
  merged reads trimmed by this function can be used by Deblur but not
  Dada2. Its statistical model for estimating error rates was not
  designed for pre-merged reads, instead use `qiime itsxpress trim-pair-output-unmerged`.

+----------------------------------+---------------------------------------------------------------------------------------+
|    Command-requirement           | Description                                                                           |
+----------------------------------+---------------------------------------------------------------------------------------+
|   --i-per-sample-sequences       | - The artifact that contains the sequence file(s).                                    |
+ 			           + - Either Joined Paired or just a single fastq.                                        +
|                                  | - One file sequence in the qza data folder.                                           |
+----------------------------------+---------------------------------------------------------------------------------------+
|       --p-region                 | - The regions ITS2, ITS1, and ALL.                                                    |
+----------------------------------+---------------------------------------------------------------------------------------+
|				   | -	Select the taxonomic group sequenced: A, B, C, D, E, F, G, H, I, L, M, O, P,	   |
+	--p-taxa		   +	Q, R, S, T, U, V, ALL.								   +
| 				   |											   |
+----------------------------------+---------------------------------------------------------------------------------------+
|       --p-threads 	           | - The amount of threads to use.                                                       |
+----------------------------------+---------------------------------------------------------------------------------------+
|       --o-trimmed                | - The resulting trimmed sequences from ITSxpress in a qza format.                     |
+----------------------------------+---------------------------------------------------------------------------------------+
|      --cluster-id                | - The percent identity for clustering reads, set to 1 for exact dereplication.        |
+----------------------------------+---------------------------------------------------------------------------------------+

3. qiime itsxpress trim-pair-output-unmerged

  This command takes paired-end data and returns unmerged, trimmed reads. The
  merged reads trimmed by this function can be used by Dada2 but not Deblur.
  For Deblur use `qiime itsxpress trim-pair`.

+----------------------------------+---------------------------------------------------------------------------------------+
|    Command-requirement           | Description                                                                           |
+----------------------------------+---------------------------------------------------------------------------------------+
|   --i-per-sample-sequences       | - The artifact that contains the sequence file.                                       |
+ 			           + - Only paired will work.                                                              +
|                                  | - Two file sequences in the qza data folder.                                          |
+----------------------------------+---------------------------------------------------------------------------------------+
|       --p-region                 | - The regions ITS2, ITS1, and ALL.                                                    |
+----------------------------------+---------------------------------------------------------------------------------------+
|				   | -	Select the taxonomic group sequenced: A, B, C, D, E, F, G, H, I, L, M, O, P,	   |
+	--p-taxa		   +	Q, R, S, T, U, V, ALL.								   +
| 				   |											   |
+----------------------------------+---------------------------------------------------------------------------------------+
|       --p-threads 	           | - The amount of threads to use.                                                       |
+----------------------------------+---------------------------------------------------------------------------------------+
|       --o-trimmed                | - The resulting trimmed sequences from ITSxpress in a qza format.                     |
+----------------------------------+---------------------------------------------------------------------------------------+
|      --cluster-id                | - The percent identity for clustering reads, set to 1 for exact dereplication.        |
+----------------------------------+---------------------------------------------------------------------------------------+

Taxa Key
--------

+-+-------------------------------------+
|A| Alveolata				|
+-+-------------------------------------+
|B| Bryophyta				|
+-+-------------------------------------+
|C| Bacillariophyta			|
+-+-------------------------------------+
|D| Amoebozoa				|
+-+-------------------------------------+
|E| Euglenozoa				|
+-+-------------------------------------+
|F| Fungi				|
+-+-------------------------------------+
|G| Chlorophyta (green algae)		|
+-+-------------------------------------+
|H| Rhodophyta (red algae)		|
+-+-------------------------------------+
|I| Phaeophyceae (brown algae)		|
+-+-------------------------------------+
|L| Marchantiophyta (liverworts)	|
+-+-------------------------------------+
|M| Metazoa				|
+-+-------------------------------------+
|O| Oomycota				|
+-+-------------------------------------+
|P| Haptophyceae (prymnesiophytes)	|
+-+-------------------------------------+
|Q| Raphidophyceae			|
+-+-------------------------------------+
|R| Rhizaria				|
+-+-------------------------------------+
|S| Synurophyceae			|
+-+-------------------------------------+
|T| Tracheophyta (higher plants)	|
+-+-------------------------------------+
|U| Eustigmatophyceae			|
+-+-+-----------------------------------+
|ALL| All				|
+---+-----------------------------------+



Example
-------

Use case: Trimming the ITS2 region from a fungal amplicon
sequencing dataset with a PairedSequencesWithQuailty qza using two cpu threads.
The example file used is in the Tests folder under paired.qza.

.. code:: bash

  qiime itsxpress trim-pair --i-per-sample-sequences ~/parired.qza --p-region ITS2 \
  --p-taxa F --p-threads 2 --o-trimmed ~/Desktop/out.qza

License information
-------------------

This software is a work of the United States Department of Agriculture,
Agricultural Research Service and is released under a Creative Commons CC0
public domain attribution.

.. _ITSxpress: https://github.com/USDA-ARS-GBRU/itsxpress
.. _ITSx: http://microbiology.se/software/itsx/
.. _BioConda: https://bioconda.github.io/
