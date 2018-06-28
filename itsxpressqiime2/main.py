"""ITSxpress-qiime2: A qiime2 plugin to rapidly trim ITS amplicon sequences from Fastq files
Author: Adam Rivers, USDA Agricultural Reseach Service and Kyle Weber
The internally transcribed spacer region is a region between highly conserved the small
subunit (SSU) of rRNA and the large subunit (LSU) of the rRNA. In Eukaryotes it contains
the 5.8s genes and two variable length spacer regions. In amplicon sequening studies it is
common practice to trim off the conserved (SSU, 5,8S or LSU) regions. Bengtsson-Palme
et al. (2013) published software the software package ITSx to do this.
ITSxpress is a high-speed implementation of the methods in	ITSx. It can process a typical
ITS amplicon sample with 100,000 read pairs in about 5 minutes, aproxamatly 100x faster.
It also trims fastq files rather than just fasta files.
Process:
	* Merges and error corrects reads using bbduk if reade are paired-end
	* Deduplicates reads using Vmatch to eliminate redundant hmm searches
	* Searches for conserved regions using the ITSx hmms, useing HMMsearch:
	  https://cryptogenomicon.org/2011/05/27/hmmscan-vs-hmmsearch-speed-the-numerology/
	* Parses everyting in python returning (optionally gzipped) fastq files.
Refernce:
	Johan Bengtsson-Palme, Vilmar Veldre, Martin Ryberg, Martin Hartmann, Sara Branco,
	Zheng Wang, Anna Godhe, Yann Bertrand, Pierre De Wit, Marisol Sanchez,
	Ingo Ebersberger, Kemal Sanli, Filipe de Souza, Erik Kristiansson, Kessy Abarenkov,
	K. Martin Eriksson, R. Henrik Nilsson. (2013). ITSx: Improved software detection
	and extraction of ITS1 and ITS2 from ribosomal ITS sequences of fungi and other
	eukaryotes for use in environmental sequencing. Methods in Ecology and Evolution,
	4: 914-919, 2013 (DOI: 10.1111/2041-210X.12073)
"""

import os

import yaml

from itsxpress import main as itsx

from itsxpress.definitions import (ROOT_DIR,
                                   taxa_dict)

from q2_types.per_sample_sequences import (SingleLanePerSamplePairedEndFastqDirFmt,
                                           SingleLanePerSampleSingleEndFastqDirFmt,
                                           FastqManifestFormat,
                                           YamlFormat)
import tempfile

def _view_artifact_type():

    """Opens the metadata file and looks for the 'type'.

    Returns:

        (str): The artifact type in the metadata file.

    Raises:

    ValueError: If the metadata file is missing or the 'type' is missing in the metadata file.

    """

    try:
        os.chdir("../")
        path = os.path.join(os.getcwd(), "metadata.yaml")
        fnOpen = open(path, "r")

        for line in fnOpen:
            if 'type:' in line:
                artifactType = line[6:]
                break

        fnOpen.close()
        return artifactType
    except:
        raise ValueError("The metadata file of the qza you entered is missing or the 'type:' in the file is missing.")

def _set_fastqs_and_check(per_sample_sequences, artifactType, sequence, singleEnd, threads):

    """Checks and writes the fastqs as well as if there paired end, interleaved and single end.

        Args:

            per_sample_sequences (SingleLanePerSampleSingleEndFastqDirFmt): The SingleLanePerSampleSingleEndFastqDirFmt type
            of the input.
            artifactType (str): The artifact type in the metadata file.
            sequence (list): The list of sequences and their IDs
            singleEnd (bool): If the sequences are singled ended or not
            threads (int): The amount of threads to use

        Returns:

            (zip lists): The path/name of the sequences.
            (bool): If single end is true or false

        Raises:

            ValueError1: BBTools error or fastq format issue.
            ValueError2: BBmerge error.

        """
    # Setting a temp folder
    dirt = tempfile.tempdir
    # Setting the fastq files and if singleEnd is used.
    fastq = os.path.join(str(per_sample_sequences.path), str(sequence[0]))
    if "SampleData[PairedEndSequencesWithQuality]" in artifactType:
        fastq2 = os.path.join(str(per_sample_sequences.path), str(sequence[1]))
    else:
        fastq2 = sequence[1]
    sequenceID = sequence[2]
    # checking fastqs
    try:
        itsx._check_fastqs(fastq, fastq2)
        # Parse input types
        paired_end, interleaved = itsx._is_paired(fastq, fastq2, singleEnd)
    except:
        raise ValueError("There is a problem with the fastq file(s) you selected or\n"
                         "BBtools was not found. check that the BBtools reformat.sh package is executable.")
    # Create SeqSample objects and merge if needed.
    try:
        if paired_end and interleaved:
            sobj = itsx.SeqSamplePairedInterleaved(fastq=fastq, tempdir=dirt)
            sobj._merge_reads(threads=threads)

        elif paired_end and not interleaved:
            sobj = itsx.SeqSamplePairedNotInterleaved(fastq=fastq, fastq2=fastq2, tempdir=dirt)
            sobj._merge_reads(threads=threads)

        elif not paired_end and not interleaved:
            sobj = itsx.SeqSampleNotPaired(fastq=fastq, tempdir=dirt)

    except:
        raise ValueError("BBmerge was not found. check that the BBmerge reformat.sh package is executible")

    return sequenceID, sobj

def _write_metadata(results):

    """Writes the metadata for the output qza as phred-offset33

    Args:

        results (SingleLanePerSampleSingleEndFastqDirFmt): The SingleLanePerSampleSingleEndFastqDirFmt type
        of the output.
    """

    metadata = YamlFormat()
    metadata.path.write_text(yaml.dump({'phred-offset': 33}))
    results.metadata.write_data(metadata, YamlFormat)

def _fastq_id_maker(per_sample_sequences, artifactType):

    """Iterates among the manifest to get the file path/name.

    Args:

        per_sample_sequences (SingleLanePerSampleSingleEndFastqDirFmt): The SingleLanePerSampleSingleEndFastqDirFmt type
        of the input.
        artifactType (str): The artifact type in the metadata file.

    Returns:

        (zip lists): The path/name of the sequences.
        (bool): If single end is true or false

    """

    path = os.path.join(str(per_sample_sequences.path), "MANIFEST")
    fn = open(path, "r")
    sampleForward = []
    sampleReverse = []
    outputNames = []
    singleEnd = False
    for line in fn:

        if "#" in line:
            continue

        elif "sample-id,filename,direction" in line:
            continue

        elif ",forward" in line:
            line = line.split(",")
            sampleForward.append(line[1].replace(",",""))
            outputNames.append(line[0].replace(",",""))

        elif ",reverse" in line:

            if "SampleData[PairedEndSequencesWithQuality]" in artifactType:
                line = line.split(",")
                sampleReverse.append(line[1].replace(",",""))

            else:
                line = line.split(",")
                sampleForward.append(line[1].replace(",",""))
                sampleReverse.append(None)

    if (len(sampleForward) != len(sampleReverse)
    and ("SampleData[PairedEndSequencesWithQuality]" in artifactType)):

        raise ValueError("The number of forward and reverse samples do not match.")

    else:
        sampleForward.sort()
        sampleReverse.sort()
        outputNames.sort()
        sampleIds = zip(sampleForward, sampleReverse, outputNames)
        if "SampleData[SequencesWithQuality]" in artifactType:
            singleEnd = True
    return sampleIds, singleEnd

def _taxa_prefix_to_taxa(taxa_prefix):

    """Turns the taxa prefix letter into the taxa

        Args:
            taxa_prefix (str): The taxa prefix that will be converted to taxa.

        Returns:

            (str): The Taxa
    """
    taxa_dict = {"A": "Alveolata", "B": "Bryophyta", "C": "Bacillariophyta", "D": "Amoebozoa", "E": "Euglenozoa",
                 "F": "Fungi","G": "Chlorophyta", "H": "Rhodophyta", "I": "Phaeophyceae", "L": "Marchantiophyta",
                 "M": "Metazoa","N": "Microsporidia","O": "Oomycota", "P": "Haptophyceae", "Q": "Raphidophyceae",
                 "R": "Rhizaria", "S": "Synurophyceae","T": "Tracheophyta", "U": "Eustigmatophyceae", "X": "Apusozoa",
                 "Y": "Parabasalia"}
    taxa_choice = taxa_dict[taxa_prefix]
    return taxa_choice

# The ITSxpress handling
def main(per_sample_sequences, threads, taxa, region):
    """The main communtion between the pluin and the ITSxpress program.

    Args:

        per_sample_sequences (SingleLanePerSampleSingleEndFastqDirFmt): The SingleLanePerSampleSingleEndFastqDirFmt type
        of the input.
        threads (int) : The number of threads to use.
        taxa (str): The taxa to be used for the search.
        region (str) : The region to be used for the search.

    Returns:

        (SingleLanePerSampleSingleEndFastqDirFmt): The SingleLanePerSampleSingleEndFastqDirFmt type
        of the output.

    Raises:

        ValueError1: hmmsearch error.

    """
    # Setting the current dir
    os.chdir(str(per_sample_sequences.path))
    # Finding the artifact type.
    artifactType = _view_artifact_type()
    # Setting the taxa
    taxa = _taxa_prefix_to_taxa(taxa)
    # Writing the manifest for the output qza
    manifest = FastqManifestFormat()
    manifest_fn = manifest.open()
    manifest_fn.write('sample-id,filename,direction\n')
    # Getting the sequences from the manifest
    sequences,singleEnd = _fastq_id_maker(per_sample_sequences, artifactType)
    sequenceList = set(sequences)
    barcode = 0
    # Creating result dir
    results = SingleLanePerSampleSingleEndFastqDirFmt()
    # Running the for loop for each sample
    for sequence in sequenceList:
        # writing fastqs and there attributes and checking the files
        sequenceID,sobj=_set_fastqs_and_check(per_sample_sequences,artifactType,sequence,singleEnd,threads)

        # Deduplicate
        sobj._deduplicate(threads=threads)
        try:
            # HMMSearch for ITS regions
            hmmfile = os.path.join(ROOT_DIR, "ITSx_db", "HMMs", taxa_dict[taxa])
            sobj._search(hmmfile=hmmfile, threads=threads)
        except:
            raise ValueError("hmmsearch was not found, make sure HMMER3 is installed and executible")

        # Parse HMMseach output.
        its_pos = itsx.ItsPosition(domtable=sobj.dom_file, region=region)
        # Create deduplication object.
        dedup_obj = itsx.Dedup(uc_file=sobj.uc_file, rep_file=sobj.rep_file, seq_file=sobj.seq_file)
        pathForward = results.sequences.path_maker(sample_id=sequenceID,
                                                   barcode_id=barcode,
                                                   lane_number=1,
                                                   read_number=1)

        manifest_fn.write("{},{},forward\n".format(sequenceID,pathForward.name))
        # Create trimmed sequences.
        dedup_obj.create_trimmed_seqs(str(pathForward), gzipped=True, itspos=its_pos)
        # Deleting the temp files.
        itsx.shutil.rmtree(sobj.tempdir)
        # Adding one to the barcode
        barcode += 1
    #Writing out the results.
    manifest_fn.close()
    _write_metadata(results)
    results.manifest.write_data(manifest, FastqManifestFormat)
    return results

# Separating the functions from the commands

# First command Trim for SingleLanePerSampleSingleEndFastqDirFmt
def trim_single(per_sample_sequences: SingleLanePerSampleSingleEndFastqDirFmt,
         region: str,
         taxa : str="F",
         threads: int=1)-> SingleLanePerSampleSingleEndFastqDirFmt:

    results = main(per_sample_sequences, threads, taxa, region)
    return results

# Second command Trim for SingleLanePerSamplePairedEndFastqDirFmt
def trim_pair(per_sample_sequences: SingleLanePerSamplePairedEndFastqDirFmt,
               region: str,
               taxa: str="F",
               threads: int=1) -> SingleLanePerSampleSingleEndFastqDirFmt:

    results = main(per_sample_sequences, threads, taxa, region)
    return results
