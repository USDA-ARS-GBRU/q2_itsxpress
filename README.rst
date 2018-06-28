ITSxpress-qiime2: A qiime2 plugin using ITSxpress to rapidly trim the Internally transcribed spacer (ITS) region of FASTQ files
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

.. image:: https://travis-ci.org/kweber1/q2_itsxpress.svg?branch=master
  :target: https://travis-ci.org/kweber1/q2_itsxpress
  
.. image:: https://codecov.io/gh/kweber1/q2_itsxpress/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/kweber1/q2_itsxpress
  
.. image:: https://api.codacy.com/project/badge/Grade/4d00341b4abc4e04a77cf5ca6674cd3c
  :target: https://www.codacy.com/app/kweber1/q2_itsxpress?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=kweber1/q2_itsxpress&amp;utm_campaign=Badge_Grade
  
Authors
-------
* Adam R. Rivers, US Department of Agriculture, Agricultural Research Service
  
* Kyle Weber

Introduction
------------

The internally transcribed spacer region is a region between highly conserved the small subunit (SSU) of rRNA and the large subunit (LSU) of the rRNA. In Eukaryotes it contains the 5.8s genes and two variable length spacer regions. In amplicon sequening studies it is common practice to trim off the conserved (SSU, 5,8S or LSU) regions. Bengtsson-Palme et al. (2013) published software the software package ITSx to do this.

ITSxpress is designed to support the calling of exact sequence variants rather than OTUs. This newer method of sequence error-correction requires quality score data from each sequence, so each input sequence must be trimmed. ITSXpress makes this possible by taking FASTQ data, de-replicating the sequences then identifying the start and stop sites using HMMSearch. Results are parsed and the trimmed files are returned. The ITS 1, ITS2 or the entire ITS region including the 5.8s rRNA gene can be selected. ITSxpress uses the hmm model from ITSx so results are comprable.

Installation:
-------------

Requirments/Dependencies
________________________

In order to consider installing ITSxpress-qiime2, it is **required to install qiime2.**

Here is the installation page to install qiime at <https://docs.qiime2.org/2018.6/install/>

Once qiime2 has been installed, you can now install the dependices for ITSxpress.

**Make sure you are in your qiime2 environment when installing the dependices and ITSxpress-qiime2.**

1. Add the bioconda channel to qiime2.

.. code-block:: bash

	conda config --add channels bioconda
			 
2. Install hmmer.

.. code-block:: bash

	conda install hmmer
		
3. Install bbmap.

.. code-block:: bash

	conda install bbmap
	
4. Install ITSxpress.

.. code-block:: bash

	conda install itsxpress 
	
or
	
.. code-block:: bash
	
	pip install itsxpress
		
ITSxpress-qiime2
________________

To install the ITSxpress plugin for qiime, a few steps are needed.
		
1. Install the project using pip. Make sure you install it in qiime2.

.. code-block:: bash

	pip install itsxpressqiime2
		
2. In your qiime2 environment, refresh the plugins.
	
.. code-block:: bash

	qiime dev refresh-cache
		
3. Check to see if the ITSxpress plugin is installed.

.. code-block:: bash

	qiime itsxpress
	
.. image:: https://i.gyazo.com/2216236a43c75a92174185b4d81a2eb5.png

Usage
-----

The main command being:

.. code-block:: bash

	qiime itsxpress

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
|				   | -	Select the taxonomic group sequenced: A, B, C, D, E, F, G, H, I, L, M, N, O, P,	   |
+	--p-taxa		   +	Q, R, S, T, U, V, X, Y.								   +
| 				   |											   |
+----------------------------------+---------------------------------------------------------------------------------------+
|       --p-threads 	           | - The amount of threads to use.                                                       | 
+----------------------------------+---------------------------------------------------------------------------------------+
|       --o-trimmed                | - The resulting trimmed sequences from ITSxpress in a qza format.                     |
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
|				   | -	Select the taxonomic group sequenced: A, B, C, D, E, F, G, H, I, L, M, N, O, P,	   |
+	--p-taxa		   +	Q, R, S, T, U, V, X, Y.								   +
| 				   |											   |
+----------------------------------+---------------------------------------------------------------------------------------+
|       --p-threads 	           | - The amount of threads to use.                                                       | 
+----------------------------------+---------------------------------------------------------------------------------------+
|       --o-trimmed                | - The resulting trimmed sequences from ITSxpress in a qza format.                     |
+----------------------------------+---------------------------------------------------------------------------------------+

Taxa Key
________

+-+-------------------------------------+
|A| Alveolata				|
+-+-------------------------------------+
|B| Bryophyta				|
+-+-------------------------------------+
|C| Bacillariophyta			|
+-+-------------------------------------+
|D| Amoebozoa				|
+-+-------------------------------------+
|E+ Euglenozoa				|
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
|N| Microsporidia			|
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
+-+-------------------------------------+
|X| Apusozoa				|
+-+-------------------------------------+
|Y| Parabasalia				|
+-+-------------------------------------+

Example
-------
Use case: Trimming the ITS1 region from a fungal amplicon sequencing dataset with a PairedSequencesWithQuailty qza using two cpu threads.The example file used is in the Tests folder.

.. code:: bash

	qiime itsxpress trim-pair --i-per-sample-sequences ~/in.qza --p-region ITS2 \
	--p-taxa F --p-threads 2 --o-trimmed ~/Desktop/out.qza

Known errors
------------
1. If an error is with hmmer:

.. image:: https://i.gyazo.com/fc4c32c86971af29347d035d809f3fb0.png

To resolve this issue:

1. Change your directory to itsxpress database in your qiime environment.

.. code-block:: bash

	cd ~/miniconda3/envs/qiime*/lib/python3.5/site-packages/itsxpress/ITSx_db

or
 
.. code-block:: bash
 
	cd ~/anaconda3/envs/qiime*/lib/python3.5/site-packages/itsxpress/ITSx_db
	
2. List the current files in the directory.

.. code-block:: bash

	ls -d .!(|.)
and

.. code-block:: bash

	ls
	
If you see this output:

.. image:: https://i.gyazo.com/d22a6fb16c6f08430e2c78dd81acd7a9.png

or

.. image:: https://i.gyazo.com/eedbdf635f2d8449ecb48048dce17420.png

You may need to delete the .DS_store file and place the true ITSx_dp file in its place.
You can download the ITSx_dp file from here <https://github.com/USDA-ARS-GBRU/itsxpress/tree/master/itsxpress/ITSx_db>

If you see this output:

.. image:: https://i.gyazo.com/fbe68023f85325d301b817575b793baa.png

This may mean another error such as hmmer may not be installed.

License information
-------------------

This software is a work of the United States Department of Agriculture, Agricultural Research Service. 17 U.S.C. 	Section 105 states that "Copyright protection under this title is not available for any work of the United States 	Government". While I anticipate that this work will be released under a CC0 public domain attribution, only the USDA 	ARS Office of Technology transfer has the authority to make that determination.
