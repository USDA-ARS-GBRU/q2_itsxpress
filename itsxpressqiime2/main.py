from q2_types.per_sample_sequences import SingleLanePerSamplePairedEndFastqDirFmt, SingleLanePerSampleSingleEndFastqDirFmt, FastqGzFormat, FastqManifestFormat, YamlFormat
import os, yaml
from itsxpress import main as itsx
from itsxpress.definitions import ROOT_DIR, taxa_dict


def _view_artifact_type(qzaPath):

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

        raise ValueError("The metadata file of the QZA you entered is missing or the 'type' in the file is missing.")


def _amount_of_files_in_data(qzaPath):

    try:
        fnlist = list(os.listdir(qzaPath))

        fnNumber = len(fnlist)

        return fnNumber

    except:

        raise ValueError("Internal error with file sequence file counting.")


def _write_metadata(results):

    metadata = YamlFormat()
    metadata.path.write_text(yaml.dump({'phred-offset': 33}))
    results.metadata.write_data(metadata, YamlFormat)

def _fastq_id_maker(per_sample_sequences):


        iterator = per_sample_sequences.sequences.iter_views(FastqGzFormat)

        sampleIds = []

        for id, (fname, fp) in enumerate(iterator):

            sampleIds.append(str(fname))

        return sampleIds


def _set_fastq_files(artifactType, qzapath, per_sample_sequences):


    seqAmount = int(_amount_of_files_in_data(qzapath)) - 2

    sequences = _fastq_id_maker(per_sample_sequences)

    singleEnd = False

    sample_id = sequences[0]

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

        raise ValueError("The QZA file and its type '{}' you entered contains more than the amount of sequence files for the ITSxpress program. ".format(artifactType),artifactType)

    return fastq, fastq2, singleEnd, sample_id


def trim(per_sample_sequences: SingleLanePerSampleSingleEndFastqDirFmt,
         region: str,
         taxa: str,
         threads: int)-> SingleLanePerSampleSingleEndFastqDirFmt:




    qzaPath = str(per_sample_sequences.path)

    artifactType = _view_artifact_type(qzaPath)

    fastq, fastq2, singleEnd, sample_id = _set_fastq_files(artifactType, qzaPath, per_sample_sequences)

    #Ruuning ITSxpress

    try:
        itsx._check_fastqs(fastq, fastq2)
        paired_end, interleaved = itsx._is_paired(fastq, fastq2, singleEnd)

    except:
        raise ValueError("There is a problem with the fastq file you selected or\n"
                         "BBTools is not installed on this system.\n"
                         "Please check them before running the program again.")
    try:
        if paired_end and interleaved:

            sobj = itsx.SeqSamplePairedInterleaved(fastq=fastq, tempdir="/tmp")

            sobj._merge_reads(threads=threads)

        elif paired_end and not interleaved:

            sobj = itsx.SeqSamplePairedNotInterleaved(fastq=fastq, fastq2=fastq2,tempdir= "/tmp")

            sobj._merge_reads(threads=threads)

        elif not paired_end and not interleaved:

            sobj = itsx.SeqSampleNotPaired(fastq=fastq, tempdir="/tmp")
    except:

        raise ValueError("BBmerge is not installed on this system.")

    sobj._deduplicate(threads=threads)

    try:

        hmmfile = os.path.join(ROOT_DIR, "ITSx_db", "HMMs", taxa_dict[taxa])

        sobj._search(hmmfile=hmmfile, threads=threads)

    except:

        raise ValueError("hmmsearch was not found, make sure HMMER3 is installed and executible")

    its_pos = itsx.ItsPosition(domtable=sobj.dom_file, region=region)

    dedup_obj = itsx.Dedup(uc_file=sobj.uc_file, rep_file=sobj.rep_file, seq_file=sobj.seq_file)

    results = SingleLanePerSampleSingleEndFastqDirFmt()

    path = results.sequences.path_maker(sample_id="seq",
                                        barcode_id=1,
                                        lane_number=1,
                                        read_number=1)

    manifest = FastqManifestFormat()

    manifest_fn = manifest.open()

    manifest_fn.write('sample-id,filename,direction\n')

    manifest_fn.write("seq,{},reverse".format(path))

    manifest_fn.close()


    dedup_obj.create_trimmed_seqs(str(path), gzipped=True, itspos=its_pos)

    _write_metadata(results)

    results.manifest.write_data(manifest, FastqManifestFormat)

    itsx.shutil.rmtree(sobj.tempdir)


    return results