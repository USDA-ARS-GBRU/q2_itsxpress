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
import shutil

import yaml
from itsxpress import main as itsxpress
from itsxpress.definitions import taxa_dict, ROOT_DIR
from q2_types.per_sample_sequences import (SingleLanePerSamplePairedEndFastqDirFmt,
                                           SingleLanePerSampleSingleEndFastqDirFmt,
                                           FastqManifestFormat,
                                           YamlFormat)
from q2_types.per_sample_sequences._format import _SingleLanePerSampleFastqDirFmt
from yaml import YAMLError, load


def _view_artifact_type(per_sample_sequence: _SingleLanePerSampleFastqDirFmt) -> str:
    """Opens the metadata file and looks for the 'type'.
    Args:
        qzaPath: The path of the qza

    Returns:
        (str): The artifact type in the metadata file.

    Raises:
        ValueError: If the metadata file is missing or the 'type' is missing in the metadata file.

    """
    try:
        per_sample_sequence_str = str(per_sample_sequence.path)
        head = os.path.split(per_sample_sequence_str)
        path = os.path.join(str(head[0]), "metadata.yaml")
        with open(path, "r") as fn_open:
            return load(fn_open)["type"]

    except (NotADirectoryError,
            FileNotFoundError,
            YAMLError):

        raise ValueError("The metadata file of the qza you entered is missing or the 'type:' in the file is missing.")


def _set_fastqs_and_check(per_sample_sequences: _SingleLanePerSampleFastqDirFmt,
                          artifact_type: str,
                          sequence: dict,
                          single_end: bool,
                          threads: int) -> (str,
                                            object):
    """Checks and writes the fastqs as well as if there paired end, interleaved and single end.

        Args:
            per_sample_sequences (SingleLanePerSampleSingleEndFastqDirFmt): The SingleLanePerSampleSingleEndFastqDirFmt type
            of the input.
            artifact_type (str): The artifact type in the metadata file.
            sequence (dict): The dictionary of sequences and their IDs
            single_end (bool): If the sequences are singled ended or not
            threads (int): The amount of threads to use

        Returns:
            (lists): The sequenceIDs
            (object): Ihe sobj object

        Raises:
            ValueError1: BBTools error or fastq format issue.
            ValueError2: BBmerge error.

        """
    # Setting the fastq files and if singleEnd is used.
    fastq = os.path.join(str(per_sample_sequences.path), str(sequence["paths"][0]))
    if "SampleData[PairedEndSequencesWithQuality]" in artifact_type:
        fastq2 = os.path.join(str(per_sample_sequences.path), str(sequence["paths"][1]))
    else:
        fastq2 = None
    sequence_id = sequence["id"]
    # checking fastqs
    try:
        itsxpress._check_fastqs(fastq=fastq, fastq2=fastq2)
        # Parse input types
        paired_end, interleaved = itsxpress._is_paired(fastq=fastq,
                                                       fastq2=fastq2,
                                                       single_end=single_end)
    except (NotADirectoryError,
            FileNotFoundError,
            ModuleNotFoundError):

        raise ValueError("There is a problem with the fastq file(s) you selected or\n"
                         "BBtools was not found. check that the BBtools reformat.sh package is executable.")
    # Create SeqSample objects and merge if needed.
    try:
        if paired_end and interleaved:
            sobj = itsxpress.SeqSamplePairedInterleaved(fastq=fastq,
                                                        tempdir=None)
            sobj._merge_reads(threads=threads)
            return sequence_id, sobj

        elif paired_end and not interleaved:
            sobj = itsxpress.SeqSamplePairedNotInterleaved(fastq=fastq,
                                                           fastq2=fastq2,
                                                           tempdir=None)
            sobj._merge_reads(threads=threads)
            return sequence_id, sobj

        elif not paired_end and not interleaved:
            sobj = itsxpress.SeqSampleNotPaired(fastq=fastq,
                                                tempdir=None)
            return sequence_id, sobj

    except (ModuleNotFoundError,
            FileNotFoundError):

        raise ValueError("BBmerge was not found. check that the BBmerge reformat.sh package is executible")


def _write_metadata(results: SingleLanePerSampleSingleEndFastqDirFmt):
    """Writes the metadata for the output qza as phred-offset33

    Args:
        results (SingleLanePerSampleSingleEndFastqDirFmt): The SingleLanePerSampleSingleEndFastqDirFmt type
        of the output.

    """

    metadata = YamlFormat()
    metadata.path.write_text(yaml.dump({'phred-offset': 33}))
    results.metadata.write_data(metadata, YamlFormat)


def _fastq_id_maker(per_sample_sequences: _SingleLanePerSampleFastqDirFmt,
                    artifact_type: str) -> (tuple,
                                            bool):
    """Iterates among the manifest to get the file path/name.

    Args:
        per_sample_sequences (SingleLanePerSampleSingleEndFastqDirFmt): The SingleLanePerSampleSingleEndFastqDirFmt type
        of the input.
        artifact_type (str): The artifact type in the metadata file.

    Returns:
        (tuple): The path/name of the sequences.
        (bool): If single end is true or false

    """

    path = os.path.join(str(per_sample_sequences.path), "MANIFEST")
    sequences = []
    with open(path, "r") as fn:
        # {"id": "1", "paths": ["path/to/forward", "path/to/backward"]}
        single_end = False
        for line in fn:
            parts = line.split(",")

            if "#" in line or "sample-id,filename,direction" in line:
                continue

            if not any([sequence['id'] == parts[0] for sequence in sequences]):
                sequence = {"id": parts[0], "paths": [parts[1]]}
                sequences.append(sequence)
                continue

            for sequence in sequences:
                if sequence['id'] == parts[0] and "SampleData[PairedEndSequencesWithQuality]" in artifact_type:
                    sequence["paths"].append(parts[1])

    if (any([len(sequence['paths']) != 2 for sequence in sequences])
            and "SampleData[PairedEndSequencesWithQuality]" in artifact_type):
        raise ValueError("The number of forward and reverse samples do not match.")

    if "SampleData[SequencesWithQuality]" in artifact_type:
        single_end = True
    return sequences, single_end


def _taxa_prefix_to_taxa(taxa_prefix: str) -> str:
    """Turns the taxa prefix letter into the taxa

        Args:
            taxa_prefix (str): The taxa prefix that will be converted to taxa.

        Returns:
            (str): The Taxa

    """
    taxa_dic = {"A": "Alveolata", "B": "Bryophyta", "C": "Bacillariophyta", "D": "Amoebozoa", "E": "Euglenozoa",
                "F": "Fungi", "G": "Chlorophyta", "H": "Rhodophyta", "I": "Phaeophyceae", "L": "Marchantiophyta",
                "M": "Metazoa", "N": "Microsporidia", "O": "Oomycota", "P": "Haptophyceae", "Q": "Raphidophyceae",
                "R": "Rhizaria", "S": "Synurophyceae", "T": "Tracheophyta", "U": "Eustigmatophyceae", "X": "Apusozoa",
                "Y": "Parabasalia"}
    taxa_choice = taxa_dic[taxa_prefix]
    return taxa_choice


# First command Trim for SingleLanePerSampleSingleEndFastqDirFmt
def trim_single(per_sample_sequences: SingleLanePerSampleSingleEndFastqDirFmt,
                region: str,
                taxa: str = "F",
                threads: int = 1) -> SingleLanePerSampleSingleEndFastqDirFmt:
    results = main(per_sample_sequences=per_sample_sequences,
                   threads=threads,
                   taxa=taxa,
                   region=region)
    return results


# Second command Trim for SingleLanePerSamplePairedEndFastqDirFmt
def trim_pair(per_sample_sequences: SingleLanePerSamplePairedEndFastqDirFmt,
              region: str,
              taxa: str = "F",
              threads: int = 1) -> SingleLanePerSampleSingleEndFastqDirFmt:
    results = main(per_sample_sequences=per_sample_sequences,
                   threads=threads,
                   taxa=taxa,
                   region=region)
    return results


# The ITSxpress handling
def main(per_sample_sequences: _SingleLanePerSampleFastqDirFmt,
         threads: int,
         taxa: str,
         region: str) -> SingleLanePerSampleSingleEndFastqDirFmt:
    """The main communication between the plugin and the ITSxpress program.

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
    # Finding the artifact type.
    artifact_type = _view_artifact_type(per_sample_sequence=per_sample_sequences)
    # Setting the taxa
    taxa = _taxa_prefix_to_taxa(taxa)
    # Writing the manifest for the output qza
    manifest = FastqManifestFormat()
    manifest_fn = manifest.open()
    manifest_fn.write('sample-id,filename,direction\n')
    # Getting the sequences from the manifest
    sequences, single_end = _fastq_id_maker(per_sample_sequences=per_sample_sequences,
                                            artifact_type=artifact_type)
    barcode = 0
    # Creating result dir
    results = SingleLanePerSampleSingleEndFastqDirFmt()
    # Running the for loop for each sample
    for sequence in sequences:
        # writing fastqs and there attributes and checking the files
        sequence_id, sobj = _set_fastqs_and_check(per_sample_sequences=per_sample_sequences,
                                                  artifact_type=artifact_type,
                                                  sequence=sequence,
                                                  single_end=single_end,
                                                  threads=threads)

        # Deduplicate
        sobj._deduplicate(threads=threads)
        try:
            # HMMSearch for ITS regions
            hmmfile = os.path.join(ROOT_DIR, "ITSx_db", "HMMs", taxa_dict[taxa])
            sobj._search(hmmfile=hmmfile, threads=threads)
        except (ModuleNotFoundError,
                FileNotFoundError,
                NotADirectoryError):

            raise ValueError("hmmsearch was not found, make sure HMMER3 is installed and executable")

        # Parse HMMseach output.
        its_pos = itsxpress.ItsPosition(domtable=sobj.dom_file,
                                        region=region)
        # Create deduplication object.
        dedup_obj = itsxpress.Dedup(uc_file=sobj.uc_file,
                                    rep_file=sobj.rep_file,
                                    seq_file=sobj.seq_file)

        path_forward = results.sequences.path_maker(sample_id=sequence_id,
                                                    barcode_id=barcode,
                                                    lane_number=1,
                                                    read_number=1)

        manifest_fn.write("{},{},forward\n".format(sequence_id, path_forward.name))
        # Create trimmed sequences.
        dedup_obj.create_trimmed_seqs(str(path_forward),
                                      gzipped=True,
                                      itspos=its_pos)
        # Deleting the temp files.
        shutil.rmtree(sobj.tempdir)
        # Adding one to the barcode
        barcode += 1
    # Writing out the results.
    manifest_fn.close()
    _write_metadata(results=results)
    results.manifest.write_data(manifest,
                                FastqManifestFormat)
    return results