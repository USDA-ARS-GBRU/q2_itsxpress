import os

import itsxpressqiime2.main as itsxq

import operator

import gzip

from nose.tools import assert_raises

from q2_types.per_sample_sequences import (SingleLanePerSamplePairedEndFastqDirFmt,
                                           SingleLanePerSampleSingleEndFastqDirFmt)

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
def test_view_artifcat_type():
    testFile = os.path.join(TEST_DIR, "test_data","paired","445cf54a-bf06-4852-8010-13a60fa1598c","data")
    testData = SingleLanePerSamplePairedEndFastqDirFmt(testFile,"r")
    os.chdir(str(testData))
    exp1 = itsxq._view_artifact_type()
    if not ("SampleData[PairedEndSequencesWithQuality]" in exp1):
        raise AssertionError()
    testFile2 = os.path.join(TEST_DIR, "test_data","pairedBrokenMissingData","50d5f31a-a761-4c04-990c-e7668fe6bf00","data")
    testData2 = SingleLanePerSamplePairedEndFastqDirFmt(testFile2,"r")
    os.chdir(str(testData2))
    assert_raises(ValueError,exp2=itsxq._view_artifact_type())

def test_write_metadata():
    results = SingleLanePerSampleSingleEndFastqDirFmt()
    itsxq._write_metadata(results)
    path = results.path
    metadata = os.path.join(str(path),"metadata.yml")
    fn = open(metadata,"r")
    passed=False
    for line in fn:
        if "33" in line:
            passed=True
    fn.close()
    if not passed:
        raise AssertionError()


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
    testFile2 = os.path.join(TEST_DIR, "test_data","pairedAllForward","445cf54a-bf06-4852-8010-13a60fa1598c","data")
    testData2 = SingleLanePerSamplePairedEndFastqDirFmt(testFile2,"r")
    try:
        itsxq._fastq_id_maker(testData2, artifactType)
        passed = True
    except:
        passed = False
    if passed:
        raise AssertionError()
    testFile3 = os.path.join(TEST_DIR, "test_data","paired","445cf54a-bf06-4852-8010-13a60fa1598c","data")
    testData3 = SingleLanePerSamplePairedEndFastqDirFmt(testFile3,"r")
    artifactType = "SampleData[SequencesWithQuality]"
    itsxq._fastq_id_maker(testData3, artifactType)

def test_taxa_prefix_to_taxa():
    exp1 = itsxq._taxa_prefix_to_taxa("A")
    if not (exp1 == "Alveolata"):
        raise AssertionError()

def test_set_fastqs_and_check():
    testFile = os.path.join(TEST_DIR, "test_data","paired","445cf54a-bf06-4852-8010-13a60fa1598c","data")
    testData = SingleLanePerSamplePairedEndFastqDirFmt(testFile,"r")
    artifactType = "SampleData[PairedEndSequencesWithQuality]"
    sequences,singleEnd = itsxq._fastq_id_maker(testData, artifactType)
    threads = 1
    sequenceSet = set(sequences)
    fnTestData = (os.path.join(TEST_DIR, "test_data", "seq.fq.gz"))
    fnTestData2 = (os.path.join(TEST_DIR, "test_data", "seq2.fq.gz"))
    fnTestData3 = (os.path.join(TEST_DIR, "test_data", "seq3.fq.gz"))
    for sequence in sequenceSet:
        exp1,exp2= itsxq._set_fastqs_and_check(testData, artifactType, sequence, singleEnd, threads)
        if exp1 != '4774-1-MSITS3':
            raise AssertionError()
        fnOutputData = gzip.open(exp2.seq_file,"rt")
        outData = fnOutputData.read()
        fnOutputData.close()
        fnTestDataUnzip = gzip.open(fnTestData,"rt")
        testData = fnTestDataUnzip.read()
        fnTestDataUnzip.close()
        seqFilebool = operator.eq(testData,outData)
        if seqFilebool == False:
            raise AssertionError()
        fnTestDataUnzip2 = gzip.open(fnTestData2,"rt")
        testData = fnTestDataUnzip2.read()
        fnTestDataUnzip.close()
        seqFilebool = operator.eq(testData,outData)
        if seqFilebool == True:
            raise AssertionError()
        fnTestDataUnzip3 = gzip.open(fnTestData3,"rt")
        testData = fnTestDataUnzip3.read()
        fnTestDataUnzip.close()
        seqFilebool = operator.eq(testData,outData)
        if seqFilebool == False:
            raise AssertionError()

    testFile2 = os.path.join(TEST_DIR, "test_data","pairedBrokenWithMANIFEST","50d5f31a-a761-4c04-990c-e7668fe6bf00","data")
    testData2 = SingleLanePerSamplePairedEndFastqDirFmt(testFile2,"r")
    sequences,singleEnd = itsxq._fastq_id_maker(testData2, artifactType)
    sequenceSet = set(sequences)
    passed = False
    try:
        for sequence in sequenceSet:
            itsxq._set_fastqs_and_check(testData2, artifactType, sequence, singleEnd, threads)

        passed = True
    except:
        passed = False
    if passed:
        raise AssertionError()

def test_main():
    testFile = os.path.join(TEST_DIR, "test_data", "paired", "445cf54a-bf06-4852-8010-13a60fa1598c", "data")
    testData = SingleLanePerSamplePairedEndFastqDirFmt(testFile, "r")
    threads = 1
    taxa = "A"
    region = "ITS1"
    try:
        itsxq.main(testData, threads, taxa, region)
    except:
        raise AssertionError()
