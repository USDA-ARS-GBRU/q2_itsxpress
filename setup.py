from setuptools import setup

setup(
    name='itsxpressqiime2',
    version='1.07',
    packages=['itsxpressqiime2'],
    author='Kyle Weber',
    author_email='kweber1@ufl.edu',
    description="itsxpress qiime2 plugin",
    url='https://github.com/kweber1/ITSxpress-qiime2',
    install_requires=[
        'itsxpress'

    ],
    entry_points={
        'qiime2.plugins':['itsxpressqiime2=itsxpressqiime2.plugin_setup:plugin']
    },
    zip_safe=False
)

