
from q2_types.sample_data import SampleData

from q2_types.per_sample_sequences import SequencesWithQuality, \
                                          PairedEndSequencesWithQuality, \
                                          JoinedSequencesWithQuality

from itsxpress import definitions as taxa

from qiime2.plugin import Plugin,\
                          Str, \
                          Choices, \
                          Int

from itsxpressqiime2.main import trimSingle,\
                                 trimPair

plugin = Plugin(
    name='itsxpress',
    version='1.5.3',
    package='Itsxpress-qiime2',
    website='ITSxpress: https://github.com/USDA-ARS-GBRU/itsxpress\n'
            'ITSxpress-plugin: https://github.com/kweber1/ITSxpress-qiime2',
    description='Trimmer',
    short_description='Plugin for fast trimming'
)

plugin.methods.register_function(
    function=trimSingle,
    inputs={'per_sample_sequences': SampleData[SequencesWithQuality |
                                               JoinedSequencesWithQuality]},
    parameters={'region': Str %Choices(['ITS2',
                                        'ITS1',
                                        'ALL']),
                'taxa': Str %Choices(taxa.taxa_choices),
                'threads': Int},
    outputs=[('trimmed', SampleData[SequencesWithQuality])],
    input_descriptions={'per_sample_sequences': 'The artifact that contains the sequence file(s).\n'
                                                'Either Joined Paired or just a single fastq.\n'
                                                'One file sequences in the qza data folder'},
    parameter_descriptions={
        'region': ('The regions ITS2, ITS1, and ALL.\n'),
        'taxa': ('The selected taxonomic group sequenced.\n'),
        'threads': ('The number of processor threads to use.\n')
    },
    output_descriptions={'trimmed': 'The resulting trimmed sequences from ITSxpress.\nqi'},
    name='TrimSingle',
    description='ITSxpress is deigned to.....'
)


plugin.methods.register_function(
    function=trimPair,
    inputs={'per_sample_sequences': SampleData[PairedEndSequencesWithQuality]},
    parameters={'region': Str %Choices(['ITS2',
                                        'ITS1',
                                        'ALL']),
                'taxa': Str %Choices(taxa.taxa_choices),
                'threads': Int},
    outputs=[('trimmed', SampleData[SequencesWithQuality])],
    input_descriptions={'per_sample_sequences': 'The artifact that contains the sequence file(s).\n'
                                                'Only Paired can be used.'
                                                'Two files sequences in the qza data folder'},
    parameter_descriptions={
        'region': ('The regions ITS2, ITS1, and ALL.\n'),
        'taxa': ('The selected taxonomic group sequenced.\n'),
        'threads': ('The number of processor threads to use.\n')
    },
    output_descriptions={'trimmed': 'The resulting trimmed sequences from ITSxpress.\nqi'},
    name='TrimPaired',
    description='ITSxpress is deigned to.....'
)