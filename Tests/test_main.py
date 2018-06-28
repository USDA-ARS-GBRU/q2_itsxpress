import os

import itsxpressqiime2.main as itsxq

from nose.tools import assert_raises

from q2_types.per_sample_sequences import SingleLanePerSamplePairedEndFastqDirFmt

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
def test_view_artifcat_type():
    testFile = os.path.join(TEST_DIR, "test_data","paired","445cf54a-bf06-4852-8010-13a60fa1598c","data")
    testData = SingleLanePerSamplePairedEndFastqDirFmt(testFile,"r")
    os.chdir(str(testData))
    exp1 = itsxq._view_artifact_type()
    if not ("SampleData[PairedEndSequencesWithQuality]" in exp1):
        raise AssertionError()
    testFile2 = os.path.join(TEST_DIR, "test_data","pairedbroken","50d5f31a-a761-4c04-990c-e7668fe6bf00","data")
    testData2 = SingleLanePerSamplePairedEndFastqDirFmt(testFile2,"r")
    os.chdir(str(testData2))
    assert_raises(ValueError,exp2=itsxq._view_artifact_type())

def test_fastq_id_maker():
    testFile = os.path.join(TEST_DIR, "test_data","paired","445cf54a-bf06-4852-8010-13a60fa1598c","data")
    testData = SingleLanePerSamplePairedEndFastqDirFmt(testFile,"r")
    artifactType = "SampleData[PairedEndSequencesWithQuality]"
    exp1,exp2 = itsxq._fastq_id_maker(testData, artifactType)
    expList = []
    exp1Set = set(exp1)
    for sequence in exp1Set:
        expList.append(sequence[0])
        expList.append(sequence[1])
    if not expList == ['4774-1-MSITS3_0_L001_R1_001.fastq.gz', '4774-1-MSITS3_1_L001_R2_001.fastq.gz']:
        raise AssertionError()
    if exp2 != False:
        raise AssertionError()

def test_taxa_prefix_to_taxa():
    exp1 = itsxq._taxa_prefix_to_taxa("A")
    if not (exp1 == "Alveolata"):
        raise AssertionError()



