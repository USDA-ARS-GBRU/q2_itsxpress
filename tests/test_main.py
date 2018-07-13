import os
from sys import getsizeof

from nose.tools import eq_, raises
from q2_types.per_sample_sequences import (SingleLanePerSampleSingleEndFastqDirFmt,
                                           SingleLanePerSamplePairedEndFastqDirFmt)

import q2_itsxpress._itsxpress as _itsxpress

# The test data dir
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
# Test info 1
TEST_FILE = os.path.join(TEST_DIR,
                         "test_data",
                         "paired",
                         "445cf54a-bf06-4852-8010-13a60fa1598c",
                         "data")

TEST_DATA = SingleLanePerSamplePairedEndFastqDirFmt(TEST_FILE, "r")
# Test info 2
TEST_FILE_PBMD = os.path.join(TEST_DIR,
                              "test_data",
                              "pairedBrokenMissingData",
                              "50d5f31a-a761-4c04-990c-e7668fe6bf00",
                              "data")

TEST_DATA_PBMD = SingleLanePerSamplePairedEndFastqDirFmt(TEST_FILE_PBMD, "r")
# Test info 3
TEST_FILE_PAF = os.path.join(TEST_DIR,
                             "test_data",
                             "pairedAllForward",
                             "445cf54a-bf06-4852-8010-13a60fa1598c",
                             "data")
TEST_DATA_PAF = SingleLanePerSamplePairedEndFastqDirFmt(TEST_FILE_PAF, "r")
# Test info 4
TEST_FILE_OUT = os.path.join(TEST_DIR,
                             "test_data",
                             "out",
                             "d9955749-00d5-44ae-a628-4b2da43000e1",
                             "data")
TEST_DATA_OUT = SingleLanePerSamplePairedEndFastqDirFmt(TEST_FILE_OUT, "r")
# Test info 5
TEST_FILE_SINGLEOUT = os.path.join(TEST_DIR,
                                   "test_data",
                                   "singleOut",
                                   "75aea4f5-f10e-421e-91d2-feda9fe7b2e1",
                                   "data")
TEST_DATA_SINGLEOUT = SingleLanePerSamplePairedEndFastqDirFmt(TEST_FILE_SINGLEOUT, "r")
# Test info 6
TEST_FILE_SINGLEIN = os.path.join(TEST_DIR,
                                  "test_data",
                                  "singleIn",
                                  "cfd0e65b-05fb-4329-9618-15ecd0aec9b3",
                                  "data")
TEST_DATA_SINGLEIN = SingleLanePerSampleSingleEndFastqDirFmt(TEST_FILE_SINGLEIN, "r")
# Test artifact1
ARTIFACT_TYPE_P = "SampleData[PairedEndSequencesWithQuality]"
# Test artifact2
ARTIFACT_TYPE_S = "SampleData[SequencesWithQuality]"


def test_view_artifcat_type():
    exp1 = _itsxpress._view_artifact_type(per_sample_sequence=TEST_DATA)
    eq_("SampleData[PairedEndSequencesWithQuality]", exp1)
    raises(ValueError, lambda: _itsxpress._view_artifact_type(per_sample_sequence=TEST_DATA_PBMD))


def test_write_metadata():
    results = SingleLanePerSampleSingleEndFastqDirFmt()
    _itsxpress._write_metadata(results)
    path = results.path
    metadata = os.path.join(str(path), "metadata.yml")
    with open(metadata, "rt") as fn:
        eq_("{phred-offset: 33}", fn.readline().replace("\n", ""))


def test_fastq_id_maker():
    exp1, exp2 = _itsxpress._fastq_id_maker(per_sample_sequences=TEST_DATA,
                                            artifact_type=ARTIFACT_TYPE_P)
    for sequence in exp1:
        eq_(sequence["paths"][0], '4774-1-MSITS3_0_L001_R1_001.fastq.gz')
        eq_(sequence["paths"][1], '4774-1-MSITS3_1_L001_R2_001.fastq.gz')
        eq_(exp2, False)
    raises(ValueError, lambda: _itsxpress._fastq_id_maker(per_sample_sequences=TEST_DATA_PAF,
                                                          artifact_type=ARTIFACT_TYPE_P))
    exp3 = _itsxpress._fastq_id_maker(per_sample_sequences=TEST_DATA,
                                      artifact_type=ARTIFACT_TYPE_S)
    eq_(exp3[1], True)


def test_taxa_prefix_to_taxa():
    exp1 = _itsxpress._taxa_prefix_to_taxa(taxa_prefix="A")
    eq_(exp1, "Alveolata")


def test_set_fastqs_and_check():
    sequences, single_end = _itsxpress._fastq_id_maker(per_sample_sequences=TEST_DATA,
                                                       artifact_type=ARTIFACT_TYPE_P)
    for sequence in sequences:
        exp1 = _itsxpress._set_fastqs_and_check(per_sample_sequences=TEST_DATA,
                                                artifact_type=ARTIFACT_TYPE_P,
                                                sequence=sequence,
                                                single_end=single_end,
                                                threads=1)
        eq_(exp1[0], "4774-1-MSITS3")

    sequences2, single_end2 = _itsxpress._fastq_id_maker(per_sample_sequences=TEST_DATA_PAF,
                                                         artifact_type=ARTIFACT_TYPE_S)
    for sequence2 in sequences2:
        exp2 = _itsxpress._set_fastqs_and_check(per_sample_sequences=TEST_DATA,
                                                artifact_type=ARTIFACT_TYPE_S,
                                                sequence=sequence2,
                                                single_end=single_end2,
                                                threads=1)
        eq_(exp2[0], "4774-1-MSITS3")

    sequences3= _itsxpress._fastq_id_maker(per_sample_sequences=TEST_DATA_PAF,
                                           artifact_type=ARTIFACT_TYPE_S)
    for sequence3 in sequences3[0]:
        exp3 = _itsxpress._set_fastqs_and_check(per_sample_sequences=TEST_DATA,
                                                artifact_type=ARTIFACT_TYPE_S,
                                                sequence=sequence3,
                                                single_end=True,
                                                threads=1)
        eq_(exp3[0], "4774-1-MSITS3")


def test_trim_pair():
    threads = 1
    taxa = "F"
    region = "ITS2"

    exp1 = _itsxpress.trim_pair(per_sample_sequences=TEST_DATA,
                                threads=threads,
                                taxa=taxa,
                                region=region)
    exp2 = getsizeof(exp1)
    exp3 = getsizeof(TEST_DATA_OUT)
    eq_(exp2, exp3)


def test_trim_single():
    threads = 1
    taxa = "F"
    region = "ITS2"

    exp1 = _itsxpress.trim_single(per_sample_sequences=TEST_DATA_SINGLEIN,
                                threads=threads,
                                taxa=taxa,
                                region=region)
    exp2 = getsizeof(exp1)
    exp3 = getsizeof(TEST_DATA_SINGLEOUT)
    eq_(exp2, exp3)

def test_trim_single_no_cluster():
    threads = 1
    taxa = "F"
    region = "ITS2"
    cluster_id = 1

    exp1 = _itsxpress.trim_single(per_sample_sequences=TEST_DATA_SINGLEIN,
                                threads=threads,
                                taxa=taxa,
                                region=region,
                                cluster_id=cluster_id)
    exp2 = getsizeof(exp1)
    exp3 = getsizeof(TEST_DATA_SINGLEOUT)
    eq_(exp2, exp3)

def test_trim_pair_no_hmmer():
    threads = 1
    taxa = "F"
    region = "ITS2"

    raises(ValueError, lambda: _itsxpress.trim_pair(per_sample_sequences=TEST_DATA,
                                                    threads=threads,
                                                    taxa=taxa,
                                                    region=region))


def test_trim_pair_no_bb():
    sequences, single_end = _itsxpress._fastq_id_maker(per_sample_sequences=TEST_DATA,
                                                       artifact_type=ARTIFACT_TYPE_P)
    for sequence in sequences:
        raises(ValueError, lambda: _itsxpress._set_fastqs_and_check(per_sample_sequences=TEST_DATA,
                                                                    artifact_type=ARTIFACT_TYPE_P,
                                                                    sequence=sequence,
                                                                    single_end=single_end,
                                                                    threads=1))
        raises(ValueError, lambda: _itsxpress._set_fastqs_and_check(per_sample_sequences=TEST_DATA,
                                                                    artifact_type=ARTIFACT_TYPE_S,
                                                                    sequence=sequence,
                                                                    single_end=single_end,
                                                                    threads=1))
