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
from q2_types.per_sample_sequences import SingleLanePerSamplePairedEndFastqDirFmt, \
                                          SingleLanePerSampleSingleEndFastqDirFmt, \
                                          FastqGzFormat, \
                                          FastqManifestFormat, \
                                          YamlFormat
import os
import yaml
from itsxpress import main as itsx
from itsxpress.definitions import ROOT_DIR,\
                                  taxa_dict

def _view_artifact_type(qzaPath):

    """Opens the metadata file and looks for the 'type'.

    Args:

        qzaPath (str): The path of the per_sequence_sample

    Returns:

        (str): The artifact type in the metadata file.

    Raises:

    ValueError: If the metadata file is missing or the 'type' is missing in the metadata file.

    """

    try:
        path = str(qzaPath[:-4]) + "/metadata.yaml"
        fnOpen = open(path, "r")

        for line in fnOpen:
            if 'type:' in line:
                artifactType = line[6:]
                break

        fnOpen.close()
        return artifactType
    except:
        raise ValueError("The metadata file of the qza you entered is missing or the 'type:' in the file is missing.")


def _amount_of_files_in_data(qzaPath):

    """Finds the number of file in the data folder of the qza.

    Args:

        qzaPath (int): The path of the per_sequence_sample.

    Returns:

        (int): The number of files in the data folder.

    """

    fnlist = list(os.listdir(qzaPath))
    fnNumber = len(fnlist)
    return fnNumber

def _write_metadata(results):

    """Writes the metadata for the output qza as phred-offset33

    Args:

        results (SingleLanePerSampleSingleEndFastqDirFmt): The SingleLanePerSampleSingleEndFastqDirFmt type
        of the output.
    """

    metadata = YamlFormat()
    metadata.path.write_text(yaml.dump({'phred-offset': 33}))
    results.metadata.write_data(metadata, YamlFormat)

def _fastq_id_maker(per_sample_sequences):

    """Iterates thru the sequences files to get the file path/name.

    Args:

        per_sample_sequences (SingleLanePerSampleSingleEndFastqDirFmt): The SingleLanePerSampleSingleEndFastqDirFmt type
        of the input.

    Returns:

        (list): The path/name of the sequences.

    """

    iterator = per_sample_sequences.sequences.iter_views(FastqGzFormat)
    sampleIds = []

    for ids, (fname, fp) in enumerate(iterator):
        sampleIds.append(str(fname))
        # useless code, just to call the ids, fp and holder.
        holder = ids, fp
        holder = holder

    return sampleIds

def _set_fastq_files(artifactType, qzapath, per_sample_sequences):

    """Gives the one or two fastq files that th ITSxpress program will be handed.

    Args:

        artifactType (str) : The artifact type in the metadata file.
        qzaPath (int): The path of the per_sequence_sample.
        per_sample_sequences (SingleLanePerSampleSingleEndFastqDirFmt): The SingleLanePerSampleSingleEndFastqDirFmt type
        of the input.

    Returns:

        (str): The first fastq location.
        (str): The second fastq location.
        (bool): The number of files in the data folder.
        (str): The number of files in the data folder.

    Raises:

        ValueError1: The sequences files in the qza is an invalid amount.
        ValueError2: The qza is a invalid type or the number of files is invalid

    """

    seqAmount = int(_amount_of_files_in_data(qzapath)) - 2
    sequences = _fastq_id_maker(per_sample_sequences)
    singleEnd = False

    if (("SampleData[PairedEndSequencesWithQuality]" in artifactType)
        and (seqAmount == 2)):
        try:
            fastq = qzapath + "/" + sequences[0]
            fastq2 = qzapath + "/" + sequences[1]
        except:
            raise ValueError("Invalid Sequence amount for type.")

    elif (("SampleData[SequencesWithQuality]" in artifactType)
        and (seqAmount == 1)):
        fastq = qzapath + "/" + sequences[0]
        fastq2 = None
        singleEnd = True

    elif (("SampleData[JoinedSequencesWithQuality]" in artifactType)
        and (seqAmount == 1)):
        fastq = qzapath + "/" + sequences[0]
        fastq2 = None

    else:
        raise ValueError("The qza and its type '{}' you entered contains more than\n"
                         "the amount of sequence files for the ITSxpress program. ".format(artifactType))
    return fastq, fastq2, singleEnd

# The ITSxpress handling
def _taxa_prefix_to_taxa(taxa_prefix):

    """Trun the taxa prefix letter into the taxa

        Args:
            taxa_prefix (str): The taxa prefix that will be converted to taxa.

        Returns:

            (str): The Taxa
    """
    taxa_dict = {"A": "Alveolata", "B": "Bryophyta", "C": "Bacillariophyta", "D": "Amoebozoa", "E": "Euglenozoa",
                 "F": "Fungi",
                 "G": "Chlorophyta", "H": "Rhodophyta", "I": "Phaeophyceae", "L": "Marchantiophyta", "M": "Metazoa",
                 "N": "Microsporidia",
                 "O": "Oomycota", "P": "Haptophyceae", "Q": "Raphidophyceae", "R": "Rhizaria", "S": "Synurophyceae",
                 "T": "Tracheophyta", "U": "Eustigmatophyceae", "X": "Apusozoa", "Y": "Parabasalia"}
    taxa_choice = taxa_dict[taxa_prefix]
    return taxa_choice

def main(fastq, fastq2, singleEnd, threads, taxa, region):
    """The main communtion between the pluin and the ITSxpress program.

    Args:

        fastq (str) : The first fastq location.
        fastq2 (str) : The second fastq location.
        singleEnd (bool) : boolean for if singleEnd is used or not.
        threads (int) : The number of threads to use.
        taxa (str): The taxa to be used for the search.
        region (str) : The region to be used for the search.

    Returns:

        (SingleLanePerSampleSingleEndFastqDirFmt): The SingleLanePerSampleSingleEndFastqDirFmt type
        of the output.

    Raises:

        ValueError1: BBTools error or fastq format issue.
        ValueError2: BBmerge error.
        ValueError3: hmmsearch error.

    """
    dirt = "/tmp"
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
            sobj = itsx.SeqSampleNotPaired(fastq=fastq, tempdir=dirt)
    except:
        raise ValueError("BBmerge was not found. check that the BBmerge reformat.sh package is executible")

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
    results = SingleLanePerSampleSingleEndFastqDirFmt()
    path = results.sequences.path_maker(sample_id="seq",
                                        barcode_id=1,
                                        lane_number=1,
                                        read_number=1)
    # Writing the manifest for the output qza
    manifest = FastqManifestFormat()
    manifest_fn = manifest.open()
    manifest_fn.write('sample-id,filename,direction\n')
    manifest_fn.write("seq,{},reverse".format(path))
    manifest_fn.close()
    # Create trimmed sequences.
    dedup_obj.create_trimmed_seqs(str(path), gzipped=True, itspos=its_pos)
    # Writing out the results.
    _write_metadata(results)
    results.manifest.write_data(manifest, FastqManifestFormat)
    # Deleting the temp files.
    itsx.shutil.rmtree(sobj.tempdir)
    return results
# Separating the functions from the commands
 # First command Trim for SingleLanePerSampleSingleEndFastqDirFmt

def trim_single(per_sample_sequences: SingleLanePerSampleSingleEndFastqDirFmt,
         region: str,
         taxa: str,
         threads: int)-> SingleLanePerSampleSingleEndFastqDirFmt:

    # Finding the main path of the qza.
    qzaPath = str(per_sample_sequences.path)
    # Finding the artifact type.
    artifactType = _view_artifact_type(qzaPath)
    # Setting the fastq files and if singleEnd is used.
    fastq, fastq2, singleEnd =_set_fastq_files(artifactType, qzaPath, per_sample_sequences)
    # setting the taxa
    taxa = _taxa_prefix_to_taxa(taxa)
    # Running the main ITSxpress program.
    results = main(fastq, fastq2, singleEnd, threads, taxa, region)
    return results

# Second command Trim for SingleLanePerSamplePairedEndFastqDirFmt
def trim_pair(per_sample_sequences: SingleLanePerSamplePairedEndFastqDirFmt,
               region: str,
               taxa: str,
               threads: int) -> SingleLanePerSampleSingleEndFastqDirFmt:

    # Finding the main path of the qza.
    qzaPath = str(per_sample_sequences.path)
    # Finding the artifact type.
    artifactType = _view_artifact_type(qzaPath)
    # Setting the fastq files and if singleEnd is used.
    fastq, fastq2, singleEnd = _set_fastq_files(artifactType, qzaPath, per_sample_sequences)
    # setting the taxa
    taxa = _taxa_prefix_to_taxa(taxa)
    # Running the main ITSxpress program.
    results = main(fastq, fastq2, singleEnd, threads, taxa, region)
    return results
