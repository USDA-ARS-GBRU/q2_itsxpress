import os

import subprocess

from nose.tools import assert_raises

from q2_types.per_sample_sequences import SingleLanePerSamplePairedEndFastqDirFmt

import itsxpressqiime2.main as itsxq

TEST_DIR = os.path.dirname(os.path.abspath(__file__))


def test_view_artifcat_type():
    testfile = os.path.join(TEST_DIR, "test_data","paired.qza","45cf54a-bf06-4852-8010-13a60fa1598c","data")
    os.chdir(testfile)
    exp1 = itsxq._view_artifact_type()
    assert exp1 == "SampleData[PairedEndSequencesWithQuality]"
    testfile2 = os.path.join(TEST_DIR, "test_data","pairedbroken.qza","45cf54a-bf06-4852-8010-13a60fa1598c","data")
    os.chdir(testfile2)
    assert_raises(subprocess.CalledProcessError, itsxq._view_artifact_type)

def test_fastq_id_maker():
    testfile = os.path.join(TEST_DIR, "test_data", "paired.qza")
    test_data = SingleLanePerSamplePairedEndFastqDirFmt(testfile, "r")
    artifactType = "SampleData[PairedEndSequencesWithQuality]"
    exp1 = itsxq._fastq_id_maker(test_data, artifactType)
    expList = []
    for sequence in exp1:
        expList.append(sequence[0],sequence[1])
    assert expList == ["4474-1MSITS3_0_L001_R1_001.fastq.gz","4474-1MSITS3_0_L001_R2_001.fastq.gz"]

def test_taxa_prefix_to_taxa():
    exp1 = itsxq._taxa_prefix_to_taxa("A")
    assert exp1 == "Alveolata"