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


Authors
-------
* Adam R. Rivers, US Department of Agriculture, Agricultural Research Service
* Kyle C. Weber, US Department of Agriculture, Agricultural Research Service

Introduction
------------

The internally transcribed spacer (ITS) is a region between the small subunit
and large subunit rRNA genes. In is a commonly used phylogenetic marker for
Fungi and other Eukaryotes. The ITS contains the 5.8s gene and two variable
length spacer regions. In amplicon sequencing studies it is common practice to
trim off the conserved (SSU, 5,8S or LSU) regions. Bengtsson-Palme et al. (2013)
published a software package ITSx_ to do this.

Q2_ITSxpress extends this work by rapidly trimming FASTQ sequences within
Qiime2.  Q2_ITSxpress is the Qiime2 plugin version of the stand alone command
line utility ITSxpress_. Q2_ITSxpress is designed to support the calling of
exact sequence variants rather than OTUs. This newer method of sequence
error-correction requires quality score data from each sequence, so each input
sequence must be trimmed. ITSXpress makes this possible by taking FASTQ data,
de-replicating the sequences then identifying the start and stop sites using
HMMSearch. Results are parsed and the trimmed files are returned. The ITS 1,
ITS2 or the entire ITS region including the 5.8s rRNA gene can be selected.
ITSxpress uses the hmm models from ITSx so results are nearly identical.


Requirements/Dependencies
-------------------------

* Qiime2 is required to run Q2-itsxpress (for stand alone software see ITSxpress_)
* To install Qiime2 follow these instructions: https://docs.qiime2.org/2018.6/install/

Q2_itsxpress Installation
-------------------------

1. Activate the Qiime2 conda environment

.. code-block:: bash

  source activate qiime2-2018.6

2. Install Q2_itsxpress using BioConda_. Be sure to install Q2_itsxpres in the Qiime2 environment.

.. code-block:: bash

  conda install -c bioconda bioconda itsxpress
  pip install q2-itsxpress

3. In your Qiime2 environment, refresh the plugins.

.. code-block:: bash

	qiime dev refresh-cache

4. Check to see if the ITSxpress plugin is installed. You should see an output similar to the image below.

.. code-block:: bash

	qiime itsxpress

.. image:: https://i.gyazo.com/2216236a43c75a92174185b4d81a2eb5.png

Usage
-----

Within Qiime2 you can trim paired end or single end reads either of these commands

.. code-block:: bash

	qiime itsxpress trim-single

  qiime itsxpress trim-paired

1. qiime itsxpress trim-single

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

	qiime itsxpress trim-pair --i-per-sample-sequences ~/parired.qza --p-region
	ITS2 \ --p-taxa F --p-threads 2 --o-trimmed ~/Desktop/out.qza

License information
-------------------

This software is a work of the United States Department of Agriculture,
Agricultural Research Service and is released under a Creative Commons CC0
public domain attribution.

.. _ITSxpress: https://github.com/USDA-ARS-GBRU/itsxpress
.. _ITSx: http://microbiology.se/software/itsx/
.. _BioConda: https://bioconda.github.io/
