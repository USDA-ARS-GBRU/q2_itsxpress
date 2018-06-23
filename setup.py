from setuptools import setup

setup(
    name='ITSxpress-qiime2',
    version='1.00',
    packages=['ITSxpress-qiime2'],
    author=['Adam R. Rivers','Kyle Weber'],
    author_email=['adam.rivers@ars.usda.gov','kweber1@ufl.edu'],
    description="itsxpress qiime2 plugin",
    url='https://github.com/kweber1/ITSxpress-qiime2',
    install_requires=[
        'itsxpress',
        'hmmer',
        'bbmap'
    ],
    entry_points={
        'qiime2.plugins':['Itsxpress-qiime2=itsxpressqiime2.plugin_setup:plugin']
    },
    zip_safe=False
)

