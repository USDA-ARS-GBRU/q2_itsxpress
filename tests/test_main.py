import gzip
import operator
import os

from nose.tools import assert_raises
from q2_types.per_sample_sequences import (SingleLanePerSamplePairedEndFastqDirFmt,
                                           SingleLanePerSampleSingleEndFastqDirFmt)

import q2_itsxpress._itsxpress as itsx

TEST_DIR = os.path.dirname(os.path.abspath(__file__))


def test_view_artifcat_type():
    test_file = os.path.join(TEST_DIR, "test_data", "paired", "445cf54a-bf06-4852-8010-13a60fa1598c", "data")
    test_data = SingleLanePerSamplePairedEndFastqDirFmt(test_file, "r")
    exp1 = itsx._view_artifact_type(per_sample_sequence=test_data)
    if not ("SampleData[PairedEndSequencesWithQuality]" in exp1):
        raise AssertionError()
    test_file2 = os.path.join(TEST_DIR, "test_data", "pairedBrokenMissingData", "50d5f31a-a761-4c04-990c-e7668fe6bf00",
                              "data")
    test_data2 = SingleLanePerSamplePairedEndFastqDirFmt(test_file2, "r")
    os.chdir(str(test_data2))
    assert_raises(ValueError, exp2=itsx._view_artifact_type(per_sample_sequence=test_data))


def test_write_metadata():
    results = SingleLanePerSampleSingleEndFastqDirFmt()
    itsx._write_metadata(results)
    path = results.path
    metadata = os.path.join(str(path), "metadata.yml")
    fn = open(metadata, "r")
    passed = False
    for line in fn:
        if "33" in line:
            passed = True
    fn.close()
    if not passed:
        raise AssertionError()


def test_fastq_id_maker():
    test_file = os.path.join(TEST_DIR, "test_data", "paired", "445cf54a-bf06-4852-8010-13a60fa1598c", "data")
    test_data = SingleLanePerSamplePairedEndFastqDirFmt(test_file, "r")
    artifact_type = "SampleData[PairedEndSequencesWithQuality]"
    exp1, exp2 = itsx._fastq_id_maker(per_sample_sequences=test_data, artifact_type=artifact_type)
    exp_list = []
    exp1_set = set(exp1)
    for sequence in exp1_set:
        exp_list.append(str(sequence[0]))
        exp_list.append(str(sequence[1]))
    if exp_list != ['4774-1-MSITS3_0_L001_R1_001.fastq.gz', '4774-1-MSITS3_1_L001_R2_001.fastq.gz']:
        raise AssertionError()
    if exp2 is not False:
        raise AssertionError()
    test_file2 = os.path.join(TEST_DIR, "test_data", "pairedAllForward", "445cf54a-bf06-4852-8010-13a60fa1598c", "data")
    test_data2 = SingleLanePerSamplePairedEndFastqDirFmt(test_file2, "r")
    try:
        itsx._fastq_id_maker(per_sample_sequences=test_data2, artifact_type=artifact_type)
        passed = True
    except (ValueError,
            NotADirectoryError,
            FileNotFoundError):

        passed = False
    if passed:
        raise AssertionError()
    test_file3 = os.path.join(TEST_DIR, "test_data", "paired", "445cf54a-bf06-4852-8010-13a60fa1598c", "data")
    test_data3 = SingleLanePerSamplePairedEndFastqDirFmt(test_file3, "r")
    artifact_type = "SampleData[SequencesWithQuality]"
    itsx._fastq_id_maker(per_sample_sequences=test_data3, artifact_type=artifact_type)


def test_taxa_prefix_to_taxa():
    exp1 = itsx._taxa_prefix_to_taxa(taxa_prefix="A")
    if not (exp1 == "Alveolata"):
        raise AssertionError()


def test_set_fastqs_and_check():
    test_file = os.path.join(TEST_DIR, "test_data", "paired", "445cf54a-bf06-4852-8010-13a60fa1598c", "data")
    test_data = SingleLanePerSamplePairedEndFastqDirFmt(test_file, "r")
    artifact_type = "SampleData[PairedEndSequencesWithQuality]"
    sequences, single_end = itsx._fastq_id_maker(per_sample_sequences=test_data, artifact_type=artifact_type)
    threads = 1
    sequence_set = set(sequences)
    fn_test_data = (os.path.join(TEST_DIR, "test_data", "seq.fq.gz"))
    fn_test_data2 = (os.path.join(TEST_DIR, "test_data", "seq2.fq.gz"))
    fn_test_data3 = (os.path.join(TEST_DIR, "test_data", "seq3.fq.gz"))
    for sequence in sequence_set:
        exp1, exp2 = itsx._set_fastqs_and_check(per_sample_sequences=test_data, artifact_type=artifact_type,
                                                sequence=sequence, single_end=single_end, threads=threads)
        if exp1 != '4774-1-MSITS3':
            raise AssertionError()
        fn_output_data = gzip.open(exp2.seq_file, "rt")
        out_data = fn_output_data.read()
        fn_output_data.close()
        fn_test_data_unzip = gzip.open(fn_test_data, "rt")
        test_data = fn_test_data_unzip.read()
        fn_test_data_unzip.close()
        seq_filebool = operator.eq(test_data, out_data)
        if seq_filebool is False:
            raise AssertionError()
        fn_test_data_unzip2 = gzip.open(fn_test_data2, "rt")
        test_data = fn_test_data_unzip2.read()
        fn_test_data_unzip.close()
        seq_filebool = operator.eq(test_data, out_data)
        if seq_filebool is True:
            raise AssertionError()
        fn_test_data_unzip3 = gzip.open(fn_test_data3, "rt")
        test_data = fn_test_data_unzip3.read()
        fn_test_data_unzip.close()
        seq_filebool = operator.eq(test_data, out_data)
        if seq_filebool is False:
            raise AssertionError()

    test_file2 = os.path.join(TEST_DIR, "test_data", "pairedBrokenWithMANIFEST",
                              "50d5f31a-a761-4c04-990c-e7668fe6bf00", "data")
    test_data2 = SingleLanePerSamplePairedEndFastqDirFmt(test_file2, "r")
    sequences, single_end = itsx._fastq_id_maker(per_sample_sequences=test_data2, artifact_type=artifact_type)
    sequence_set = set(sequences)
    passed = False
    try:
        for sequence in sequence_set:
            itsx._set_fastqs_and_check(per_sample_sequences=test_data2, artifact_type=artifact_type,
                                       sequence=sequence, single_end=single_end, threads=threads)
        passed = True
    except:
        passed = False
    if passed:
        raise AssertionError()


def test_main():
    test_file = os.path.join(TEST_DIR, "test_data", "paired", "445cf54a-bf06-4852-8010-13a60fa1598c", "data")
    test_data = SingleLanePerSamplePairedEndFastqDirFmt(test_file, "r")
    threads = 1
    taxa = "A"
    region = "ITS1"
    try:
        itsx.main(per_sample_sequences=test_data, threads=threads, taxa=taxa, region=region)
    except ValueError:
        raise AssertionError()
