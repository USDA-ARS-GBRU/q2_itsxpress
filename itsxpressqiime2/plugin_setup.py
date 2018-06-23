
from q2_types.sample_data import SampleData
from q2_types.per_sample_sequences import SequencesWithQuality, PairedEndSequencesWithQuality, JoinedSequencesWithQuality
from itsxpress import definitions as taxa
from qiime2.plugin import Plugin, Str, Choices, Int
from itsxpressqiime2.main import trim
from q2_types.feature_data import FeatureData, AlignedSequence

plugin = Plugin(
    name='itsxpress',
    version='1.5.1',
    package='itsxpress-qiime2',
    website='url',
    description='Trimmer',
    short_description='Plugin for fast trimming'
)


plugin.methods.register_function(
    function=trim,
    inputs={'per_sample_sequences': SampleData[PairedEndSequencesWithQuality| SequencesWithQuality | JoinedSequencesWithQuality]},
    parameters={'region': Str %Choices(['ITS2','ITS1','ALL']),
                'taxa': Str %Choices(taxa.taxa_choices),
                'threads': Int},
    outputs=[('sample', SampleData[SequencesWithQuality])],
    input_descriptions={'per_sample_sequences': 'The artifact that contains the sequence files.'},
    parameter_descriptions={
        'region': ('The regions ITS2, ITS1, and ALL.'),
        'taxa': ('The selected taxonomic group sequenced.'),
        'threads': ('The number of processor threads to use.')
    },
    output_descriptions={'sample': 'The resulting trimmed sequences.'},
    name='Trim',
    description='ITSxpress is deigned to.....'
)
