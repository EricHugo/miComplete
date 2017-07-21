from setuptools import setup, find_packages

with open('README.md', mode='r') as f:
    l_description = f.read()

setup(name='micomplete',
        version='0.1',
        description='Quality control for metagenome assembled genomes',
        long_description=l_description,
        url='https://bitbucket.org/evolegiolab/micomplete',
        author='Eric Hugoson',
        author_email='eric@hugoson.org',
        license='GNU General Public License v3 or later (G  PLv3+)',
        
        classifiers=[
            'Development Status :: 3 - Alpha',
            'Intended Audience :: Science/Research',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Environment :: Console',
            'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
            'Programming Language :: Python :: 2',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.4',
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',
        ],
        keywords='bioinformatics completeness metagenomics',
        install_requires=[
            'biopython',
            'numpy',
            'matplotlib',
            ],
        packages=['micomplete'],
        include_package_data=True,
        entry_points={
            'console_scripts': [
                'micomplete = micomplete.miComplete:main',
                ]
            },
        )

