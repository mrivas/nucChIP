from setuptools import setup

setup(name='mnaseChIP',
	version='0.0',
	description='Tools for Mnase ChIP-seq',
	url='http://mrivas.hasdocs.com/mnaseChIP/',
	author='Marcelo Rivas-Astroza',
	license='GPL',
	packages=['mnaseChIP'],
	scripts=['bin/bam2bed.py','bin/checkMatches.py','bin/getAvrCov.py','bin/getCounts.py','bin/getFigures.py','bin/getProfile.py','bin/merged.py'],
	install_requires=['numpy','sys','argparse','HTSeq','scipy','pybedtools','scipy','pysam',],
	zip_safe=False)
