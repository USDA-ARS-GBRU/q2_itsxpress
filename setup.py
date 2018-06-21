from setuptools import setup

setup(
    name='itsxpressqiime2',
    version='1.00',
    packages=['itsxpressqiime2'],
    author=['Adam R. Rivers','Kyle Weber'],
    author_email=['adam.rivers@ars.usda.gov','kweber1@ufl.edu'],
    description="itsxpress qiime2 plugin",
    url='http://github.com/usda-ars-gbru/itsxpress',
    install_requires=[
        'itsxpress',
    ],
    entry_points={
        'qiime2.plugins':['itsxpress-qiime2=itsxpressqiime2.plugin_setup:plugin']
    },
    zip_safe=False
)

