from setuptools import setup

setup(
    name='q2_itsxpress',
    version='1.8.1',
    packages=['q2_itsxpress'],
    author='Adam R. Rivers, Kyle C. Weber, Sveinn V. Einarsson',
    author_email='adam.rivers@ars.usda.gov, kweber1@ufl.edu, seinarsson@ufl.edu',
    description="A QIIME2 plugin to trim ITS regions using ITSxpress",
    long_description=open('README.rst').read(),
    url='https://github.com/usda-ars-gbru/q2_itsxpress',
    test_suite='pytest',
    license="License :: CC0 1.0 Universal (CC0 1.0) Public Domain Dedication",
    classifiers=["Topic :: Scientific/Engineering :: Bio-Informatics",
                 "Programming Language :: Python :: 3.6",
                 "Programming Language :: Python :: 3.5",
                 "Development Status :: 3 - Alpha"],
    keywords="Amplicon sequencing fungal ITS QIIME2",
    python_requires=">3.5",
    include_package_data=True,
    install_requires=[
        'qiime'
        'itsxpress'
        'pandas',
    ],
    entry_points={
        'qiime2.plugins': ['q2_itsxpress=q2_itsxpress.plugin_setup:plugin']
    },
    zip_safe=False
)
