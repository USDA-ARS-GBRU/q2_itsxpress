from itsxpress import main as itsx
from q2_types.per_sample_sequences import SingleLanePerSamplePairedEndFastqDirFmt, FastqGzFormat, FastqManifestFormat, YamlFormat, SingleLanePerSampleSingleEndFastqDirFmt
import os
import gzip
import yaml
from itsxpress.definitions import ROOT_DIR, taxa_choices, taxa_dict
from q2_types.feature_data import FeatureData

def _read_sequences(ID, mainPath):

    fn = gzip.open('{}/{}'.format(mainPath,ID))

    return fn



def _artifact_type_finder(mainPath):



    PEWSO = False
    SWQ = False
    JSWO = False

    fn = (mainPath[:-4] + "/metadata.yaml")

    fnOpened= open(fn, "r")

    for lines in fnOpened:

        if 'type' in lines:

            artifactType = lines[6:]

            break

    if '[PairedEndSequencesWithQuality]' in artifactType:

        PEWSO = True

    elif '[SequencesWithQuality]' in artifactType:

        SWQ = True

    elif '[JoinedSequencesWithQuality]' in artifactType:

        JSWO = True

    else:

        raise ValueError("Incorrect artifact type given!")

    fnOpened.close()

    return PEWSO, SWQ, JSWO




def _amount_of_files(mainPath):

    fnlist = list(os.listdir(mainPath))

    fnNumber = len(fnlist)

    return fnNumber

def _write_metadata(results):

    metadata = YamlFormat()
    metadata.path.write_text(yaml.dump({'phred-offset': 33}))
    results.metadata.write_data(metadata, YamlFormat)



#The main function to set up and run itsxpress

def trim2(per_sample_sequences: SingleLanePerSamplePairedEndFastqDirFmt,
         region: str,
         taxa: str,
         threads: int)-> FeatureData:

    results = FeatureData()

    return results



def trim(per_sample_sequences: SingleLanePerSamplePairedEndFastqDirFmt,
         region: str,
         taxa: str,
         threads: int)-> SingleLanePerSampleSingleEndFastqDirFmt:
    try:

        #setup for ITSxpress

        single_end = False

        mainPath = str(per_sample_sequences.path)

        PEWSOType, SWQType, JSWQType = _artifact_type_finder(mainPath)

        fnNum = _amount_of_files(mainPath)

        iterator = per_sample_sequences.sequences.iter_views(FastqGzFormat)


        sampleIds = []
        for id, (fname,fp) in enumerate(iterator):

            sampleIds.append(str(fname))


        if (PEWSOType
            and  fnNum == 4):

            fastq = mainPath + "/" + str(sampleIds[0])

            fastq2 = mainPath + "/" + str(sampleIds[1])


        elif (SWQType

            and fnNum == 3):

            fastq = mainPath + "/" + str(sampleIds[0])

            fastq2 = None

            single_end = True

        elif (JSWQType
            and fnNum == 3):

            fastq = mainPath + "/" + str(sampleIds[0])

            fastq2 = None

        else:

            raise ValueError("The artifact type and number of sequence files inside do not match.")


        #running ITSxpress


        itsx._check_fastqs(fastq, fastq2)

        paired_end, interleaved = itsx._is_paired(fastq, fastq2, single_end)


        if paired_end and interleaved:

            sobj = itsx.SeqSamplePairedInterleaved(fastq=fastq, tempdir="./tmp")

            sobj.itsx_merge_reads(threads=threads)

        elif paired_end and not interleaved:

            sobj = itsx.SeqSamplePairedNotInterleaved(fastq=fastq, fastq2=fastq2 ,tempdir= "./tmp")

            sobj._merge_reads(threads=threads)

        elif not paired_end and not interleaved:

            sobj = itsx.SeqSampleNotPaired(fastq=fastq, tempdir="./tmp")



        sobj._deduplicate(threads=threads)

        hmmfile = os.path.join(ROOT_DIR, "ITSx_db", "HMMs", taxa_dict[taxa])

        sobj._search(hmmfile=hmmfile, threads=threads)


        its_pos = itsx.ItsPosition(domtable=sobj.dom_file, region=region)

        dedup_obj = itsx.Dedup(uc_file=sobj.uc_file, rep_file=sobj.rep_file, seq_file=sobj.seq_file)

        results = SingleLanePerSampleSingleEndFastqDirFmt()


        path = results.path



        manifest = FastqManifestFormat()
        manifest_fn = manifest.open()
        manifest_fn.write('sample-id,filename,direction\n')
        manifest_fn.write("{},{},reverse".format(sampleIds[0],path.name))
        manifest_fn.close()
        dedup_obj.create_trimmed_seqs((str(path)+"/"+sampleIds[0]),gzipped=True, itspos=its_pos)
        results.manifest.write_data(manifest, FastqManifestFormat)
        _write_metadata(results)

    except:
        raise ValueError("error")

    return results





