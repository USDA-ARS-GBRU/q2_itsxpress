
from q2_types.sample_data import SampleData

from q2_types.per_sample_sequences import SequencesWithQuality, \
                                          PairedEndSequencesWithQuality, \
                                          JoinedSequencesWithQuality

from itsxpress import definitions as taxa

from qiime2.plugin import Plugin,\
                          Str, \
                          Choices, \
                          Int

from itsxpressqiime2.main import trim_single,\
                                 trim_pair

plugin = Plugin(
    name='itsxpress',
    version='1.5.3',
    package='Itsxpress-qiime2',
    website='https://github.com/kweber1/ITSxpress-qiime2              '
            'ITSxpress: https://github.com/USDA-ARS-GBRU/itsxpress',
    description='ITSxpress is designed to support the calling of exact sequence variants'
                'rather than OTUs. This newer method of sequence error-correction requires'
                'quality score data from each sequence, so each input sequence must be trimmed.'
                'ITSXpress makes this possible by taking FASTQ data, de-replicating the'
                'sequences then identifying the start and stop sites using HMMSearch.'
                'Results are parsed and the trimmed files are returned.'
                'The ITS 1, ITS2 or the entire ITS region including the 5.8s rRNA gene can be selected.'
                'ITSxpress uses the hmm model from ITSx so results are comprable.',
    short_description='A qiime2 plugin using ITSxpress to rapidly trim the\n'
                      'Internally transcribed spacer (ITS) region of FASTQ files.'
)

plugin.methods.register_function(
    function=trim_single,
    inputs={'per_sample_sequences': SampleData[SequencesWithQuality |
                                               JoinedSequencesWithQuality]},
    parameters={'region': Str %Choices(['ITS2','ITS1','ALL']),
                'taxa': Str %Choices(['A','B','C','D','E','F','G','H','I','L','M','N','O','P',
                                     'Q', 'R', 'S', 'T','U','X','Y']),
                'threads': Int},
    outputs=[('trimmed', SampleData[SequencesWithQuality])],
    input_descriptions={'per_sample_sequences': 'The artifact that contains the sequence file(s).'
                                                ' Either Joined Paired or just a single fastq.'
                                                ' One file sequences in the qza data folder.'},
    parameter_descriptions={
        'region': ('\nThe regions ITS2, ITS1, and ALL that can be selected from.'),
        'taxa': ('\nThe selected taxonomic group sequenced that can be selected from.'),
        'threads': ('\nThe number of processor threads to use in the run.')
    },
    output_descriptions={'trimmed': 'The trimmed sequences from ITSxpress.'},
    name='TrimSingle',
    description='ITSxpress trimSingle is used for qza types with\n'
                'SquencesWithQuality or JoinedSequencesWithQuality.'
                ' This means the qza must be in the\n'
                'SingleLanePerSampleSingleEndFastqDirFmt, and only contain\n'
                'one fastq file.\n'
                '\nThe taxa list and variable for it:\n'
                '\n# A = Alveolata\n'
                '\n# B = Bryophyta\n'
                '\n# C = Bacillariophyta\n'
                '\n# D = Amoebozoa\n'
                '\n# E = Euglenozoa\n'
                '\n# F = Fungi\n'
                '\n# G = Chlorophyta (green algae)\n'
                '\n# H = Rhodophyta (red algae)\n'
                '\n# I = Phaeophyceae (brown algae)\n'
                '\n# L = Marchantiophyta (liverworts)\n'
                '\n# M = Metazoa\n'
                '\n# N = Microsporidia\n'
                '\n# O = Oomycota\n'
                '\n# P = Haptophyceae (prymnesiophytes)\n'
                '\n# Q = Raphidophyceae\n'
                '\n# R = Rhizaria\n'
                '\n# S = Synurophyceae\n'
                '\n# T = Tracheophyta (higher plants)\n'
                '\n# U = Eustigmatophyceae\n'
                '\n# X = Apusozoa\n'
                '\n# Y = Parabasalia'
)

plugin.methods.register_function(
    function=trim_pair,
    inputs={'per_sample_sequences': SampleData[PairedEndSequencesWithQuality]},
    parameters={'region': Str %Choices(['ITS2','ITS1','ALL']),
                'taxa': Str %Choices(['A','B','C','D','E','F','G','H','I','L','M','N','O','P','Q', 'R', 'S', 'T','U','X','Y']),
                'threads': Int},
    outputs=[('trimmed', SampleData[SequencesWithQuality])],
    input_descriptions={'per_sample_sequences': 'The artifact that contains the sequence file(s). '
                                                'Only Paired can be used. '
                                                'Two files sequences in the qza data folder'},
    parameter_descriptions={
        'region': ('\nThe regions ITS2, ITS1, and ALL that can be selected from.'),
        'taxa': ('\nThe selected taxonomic group sequenced that can be selected from.'),
        'threads': ('\nThe number of processor threads to use in the run.')
    },
    output_descriptions={'trimmed': 'The resulting trimmed sequences from ITSxpress'},
    name='TrimPaired',
    description='ITSxpress trimPair is used for the qza type with\n'
                'PairedEndSquencesWithQuality. This means the qza\n'
                'must be in the SingleLanePerSamplePairedEndFastqDirFmt\n'
                'and only contain two fastq files.\n'
                '\nThe taxa list and variable for it:\n'
                '\n# A = Alveolata\n'
                '\n# B = Bryophyta\n'
                '\n# C = Bacillariophyta\n'
                '\n# D = Amoebozoa\n'
                '\n# E = Euglenozoa\n'
                '\n# F = Fungi\n'
                '\n# G = Chlorophyta (green algae)\n'
                '\n# H = Rhodophyta (red algae)\n'
                '\n# I = Phaeophyceae (brown algae)\n'
                '\n# L = Marchantiophyta (liverworts)\n'
                '\n# M = Metazoa\n'
                '\n# N = Microsporidia\n'
                '\n# O = Oomycota\n'
                '\n# P = Haptophyceae (prymnesiophytes)\n'
                '\n# Q = Raphidophyceae\n'
                '\n# R = Rhizaria\n'
                '\n# S = Synurophyceae\n'
                '\n# T = Tracheophyta (higher plants)\n'
                '\n# U = Eustigmatophyceae\n'
                '\n# X = Apusozoa\n'
                '\n# Y = Parabasalia'
)