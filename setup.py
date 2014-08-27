from setuptools import setup

setup(name='nucChIP',
	version='0.0',
	description='Tools for Mnase digested ChIP-seq',
	url='http://mrivas.hasdocs.com/nucChIP/',
	author='Marcelo Rivas-Astroza',
	license='GPL',
	packages=['nucChIP'],
	scripts=['bin/bam2bed','bin/checkMatches','bin/getAvrCov','bin/getCounts','bin/getFigures','bin/getProfile','bin/merged','bin/getBedTSS'],
	install_requires=['numpy','sys','argparse','HTSeq','scipy','pybedtools','scipy','pysam',],
	zip_safe=False)