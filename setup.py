from setuptools import setup

setup(
    name='q2_itsxpress',
    version='1.4',
    packages=['q2_itsxpress'],
    author='Kyle Weber',
    author_email='kweber1@ufl.edu',
    description="itsxpress qiime2 plugin",
    url='https://github.com/kweber1/q2_itsxpress',
    test_suite ='nose.collector',

    install_requires=[
        'itsxpress'

    ],
    entry_points={
        'qiime2.plugins':['q2_itsxpress=q2_itsxpress.plugin_setup:plugin']
    },
    zip_safe=False
)

