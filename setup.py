from setuptools import setup, find_packages

with open('README.rst', mode='r') as f:
    l_description = f.read()

setup(name='micomplete',
        version='1.1.1',
        description='Quality control of assembled genomes',
        long_description=l_description,
        url='https://bitbucket.org/evolegiolab/micomplete',
        author='Eric Hugoson',
        author_email='eric@hugoson.org',
        license='GNU General Public License v3 or later (G  PLv3+)',
        
        classifiers=[
            'Development Status :: 5 - Production/Stable',
            'Intended Audience :: Science/Research',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Environment :: Console',
            'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7',
        ],
        keywords='bioinformatics completeness metagenomics',
        install_requires=[
            'biopython>=1.53',
            'numpy>=1.13',
            'matplotlib>=2.0.2',
            ],
        packages=find_packages(),
        include_package_data=True,
        entry_points={
            'console_scripts': [
                'miComplete = micomplete.micomplete:main',
                ]
            },
        scripts=['micomplete/utils/miCompletelist.sh'],
        )

