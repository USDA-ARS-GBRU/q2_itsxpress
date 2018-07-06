from setuptools import setup

setup(
    name='q2_itsxpress',
    version='1.4.2',
    packages=['q2_itsxpress'],
    author='Adam R. Rivers, Kyle C. Weber',
    description="A QIIME2 plugin to trim ITS regions using ITSxpress",
    long_description=open('README.rst').read(),
    url='https://github.com/kweber1/q2_itsxpress',
    test_suite='nose.collector',
    license="License :: CC0 1.0 Universal (CC0 1.0) Public Domain Dedication",
    classifiers=["Topic :: Scientific/Engineering :: Bio-Informatics",
                 "Programming Language :: Python :: 3.6",
                 "Programming Language :: Python :: 3.5",
                 "Development Status :: 3 - Alpha"],
    keywords="Aplicon sequencing fungal ITS QIIME2",
    python_requires=">3.5",
    includer_package_data=True,
    install_requires=[
        'itsxpress'
    ],
    entry_points={
        'qiime2.plugins': ['q2_itsxpress=q2_itsxpress.plugin_setup:plugin']
    },
    zip_safe=False
)
