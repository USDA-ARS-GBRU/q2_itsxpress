from q2_types.per_sample_sequences import (SequencesWithQuality,
                                           PairedEndSequencesWithQuality,
                                           JoinedSequencesWithQuality)
from q2_types.sample_data import SampleData
from qiime2.plugin import (Plugin,
                           Str,
                           Choices,
                           Int,
                           Bool,
                           Float,
                           Range)

from q2_itsxpress._itsxpress import (trim_single,
                                     trim_pair)

plugin = Plugin(
    name='itsxpress',
    version='1.6.1',
    package='q2_itsxpress',
    website='https://github.com/kweber1/q2_itsxpress             '
            'ITSxpress: https://github.com/USDA-ARS-GBRU/itsxpress',
    description='ITSxpress is designed to support the calling of exact sequence variants'
                'rather than OTUs. This newer method of sequence error-correction requires'
                'quality score data from each sequence, so each input sequence must be trimmed.'
                'ITSXpress makes this possible by taking FASTQ data, de-replicating the'
                'sequences then identifying the start and stop sites using HMMSearch.'
                'Results are parsed and the trimmed files are returned.'
                'The ITS 1, ITS2 or the entire ITS region including the 5.8s rRNA gene can be selected.'
                'ITSxpress uses the hmm model from ITSx so results are comprable.',
    short_description='Plugin for using ITSxpress to rapidly trim the\n'
                      'internally transcribed spacer (ITS) region of FASTQ files.'
)

taxaList = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'L', 'M', 'ALL', 'O', 'P', 'Q', 'R', 'S', 'T', 'U']

plugin.methods.register_function(
    function=trim_single,
    inputs={'per_sample_sequences': SampleData[SequencesWithQuality |
                                               JoinedSequencesWithQuality]},
    parameters={'region': Str % Choices(['ITS2', 'ITS1', 'ALL']),
                'taxa': Str % Choices(taxaList),
                'threads': Int,
                'slow': Bool,
                'cluster_id': Float % Range(0.97, 1.0)},
    outputs=[('trimmed', SampleData[SequencesWithQuality])],
    input_descriptions={'per_sample_sequences': 'The artifact that contains the sequence file(s).'
                                                ' Either Joined Paired or just a single fastq.'
                                                ' One file sequences in the qza data folder.'},
    parameter_descriptions={
        'region': ('\nThe regions ITS2, ITS1, and ALL that can be selected from.'),
        'taxa': ('\nThe selected taxonomic group sequenced that can be selected from.'),
        'threads': ('\nThe number of processor threads to use in the run.'),
        'slow': ('\nIf True, dereplication will be used instead of clustering at high identity, default is False'),
        'cluster_id': ('\n The clustering float vaule that will be used if clustering is set to True.')
    },
    output_descriptions={'trimmed': 'The trimmed sequences from ITSxpress.'},
    name='TrimSingle',
    description='ITSxpress trimSingle is used for qza types with\n'
                'SquencesWithQuality or JoinedSequencesWithQuality.'
                ' This means the qza must be in the\n'
                'SingleLanePerSampleSingleEndFastqDirFmt, and only contain\n'
                'one fastq file.\n'
                '\nThe taxa list and variable for it:\n'
                '\nA = Alveolata\n'
                '\nB = Bryophyta\n'
                '\nC = Bacillariophyta\n'
                '\nD = Amoebozoa\n'
                '\nE = Euglenozoa\n'
                '\nF = Fungi\n'
                '\nG = Chlorophyta (green algae)\n'
                '\nH = Rhodophyta (red algae)\n'
                '\nI = Phaeophyceae (brown algae)\n'
                '\nL = Marchantiophyta (liverworts)\n'
                '\nM = Metazoa\n'
                '\nO = Oomycota\n'
                '\nP = Haptophyceae (prymnesiophytes)\n'
                '\nQ = Raphidophyceae\n'
                '\nR = Rhizaria\n'
                '\nS = Synurophyceae\n'
                '\nT = Tracheophyta (higher plants)\n'
                '\nU = Eustigmatophyceae\n'
                '\nALL = All'
)

plugin.methods.register_function(
    function=trim_pair,
    inputs={'per_sample_sequences': SampleData[PairedEndSequencesWithQuality]},
    parameters={'region': Str % Choices(['ITS2', 'ITS1', 'ALL']),
                'taxa': Str % Choices(taxaList),
                'threads': Int,
                'slow': Bool,
                'cluster_id': Float % Range(0.97, 1.0)},
    outputs=[('trimmed', SampleData[SequencesWithQuality])],
    input_descriptions={'per_sample_sequences': 'The artifact that contains the sequence file(s). '
                                                'Only Paired can be used. '
                                                'Two files sequences in the qza data folder'},
    parameter_descriptions={
        'region': ('\nThe regions ITS2, ITS1, and ALL that can be selected from.'),
        'taxa': ('\nThe selected taxonomic group sequenced that can be selected from.'),
        'threads': ('\nThe number of processor threads to use in the run.'),
        'slow': ('\nIf True, dereplication will be used instead of clustering at high identity, default is False'),
        'cluster_id': ('\n The clustering float vaule that will be used if clustering is set to True.')
    },
    output_descriptions={'trimmed': 'The resulting trimmed sequences from ITSxpress'},
    name='TrimPaired',
    description='ITSxpress trimPair is used for the qza type with\n'
                'PairedEndSquencesWithQuality. This means the qza\n'
                'must be in the SingleLanePerSamplePairedEndFastqDirFmt\n'
                'and only contain two fastq files.\n'
                '\nThe taxa list and variable for it:\n'
                '\nA = Alveolata\n'
                '\nB = Bryophyta\n'
                '\nC = Bacillariophyta\n'
                '\nD = Amoebozoa\n'
                '\nE = Euglenozoa\n'
                '\nF = Fungi\n'
                '\nG = Chlorophyta (green algae)\n'
                '\nH = Rhodophyta (red algae)\n'
                '\nI = Phaeophyceae (brown algae)\n'
                '\nL = Marchantiophyta (liverworts)\n'
                '\nM = Metazoa\n'
                '\nO = Oomycota\n'
                '\nP = Haptophyceae (prymnesiophytes)\n'
                '\nQ = Raphidophyceae\n'
                '\nR = Rhizaria\n'
                '\nS = Synurophyceae\n'
                '\nT = Tracheophyta (higher plants)\n'
                '\nU = Eustigmatophyceae\n'
                '\nALL = All'

)
