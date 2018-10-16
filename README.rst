==============
**miComplete**
==============

- Eric Hugoson (eric@hugoson.org / eric.hugoson@imbim.uu.se / @EricHugo)
- Lionel Guy (guy.lionel@gmail.com / lionel.guy@imbim.uu.se / @LionelGuy)
 

Introduction
----------------
With the increasing rate of the production of genomics data, particularly metagenomic data, there is need for more and faster quality control of the resulting genomes. An increasing usage on 
metagenome assembled genomes (MAGs) means often working with incomplete genomes, which can be acceptable providing the researcher is cognizant of this. Therefore there is a clear use for software 
able to rapidly provide the stats that describe the quality and completeness of the genomes or genomic bins of interest.

miComplete allows a user to provide a list of genomes or genomic bins to retrieve some basic statistics regarding the given genomes (size, GC-content, N- and L50, N- and L90). Further a set of marker genes 
in HMM format can be provided to also retrieve completeness and redundance of those markers in each genome. Additionally, a set of weights for the marker genes can be provided to also retrieve the
weighted versions of completeness and redundance which can inform the user a bit more of the actual state completeness (see description). Alternatively, the user can calculate new weights for any given set
of marker genes provided.

miComplete is still in a relatively early state of development, there are a few missing features and bugs are very much expected. Feedback, bug reports, and feature requests are welcome through the 
issue system.

Description
---------------
miComplete is a compact software aimed at the rapid determining of the quality of assembled genomes, often metagenome assembled. miComplete also aims at providing a more reliable completeness and redundance 
metric via a system of weighting the impact of different marker genes presence or absence differently.

Completeness
^^^^^^^^^^^^^^^
In miComplete completeness is calculated based on the presence/absence a set of marker genes provided as a set of HMMs. The presence or absence of the marker genes is determined by HMMER3 (see dependencies) 
if the hit reported is below the cutoff e-value provided. Duplicated marker genes are also gathered, and found duplications are reported as redundance, but only if the e-value of the reported duplicated 
hit is at least equal to or less than the square root of the the best hit (though this can easily be altered as desired by the user).

Weights
^^^^^^^^^^^
Not all marker genes are equal in determining the completeness of a genome. Some genes associate closely together within a genome (not in the least genes in operons), and thus viewing two genes that typically 
associate together within an operon as providing the same completeness information as two unrelated genes would be misleading. miComplete is able to provide a weighted version of completeness and redundance 
that attempts to factor in how closely the provided marker genes typically associate with other provided marker genes.

Linkage
"""""""""""""""""
Weights can be calculated for any given set of marker genes in miComplete. This is can be done by a user by providing a set of reference genomes (note that these need to be single contig chromosomes). 
The reference genomes can be any set that the user wishes, but as general rule the larger and more diverse number of genomes the better weights. At the end of the run a boxplot of the distribution of 
weights for all markers is produced.

Dependencies
--------------
Python (>=3.5)


External software
^^^^^^^^^^^^^^^^^^^
Executables should be available in the user's ``$PATH``.

HMMER3
"""""""""""""""""
HMMER: biosequence analysis using profile hidden Markov models, by Sean Eddy and coworkers. Tested with v. 3.1b2. Available from <http://hmmer.org/>.

prodigal
""""""""""""""""
A gene prediction software by Doug Hyatt. Tested with v. 2.6.3. Download at: 
<https://github.com/hyattpd/Prodigal>

Python libraries
^^^^^^^^^^^^^^^^^^^
If built from the package these will be installed automatically, otherwise can easily be installed using ``pip`` (`Install pip <https://pip.pypa.io/en/stable/installing/>`_).

Required
""""""""""""""""""

- Biopython (>= 1.70) (``$ pip install biopython``)
- Numpy (>= 1.13.1) (``$ pip install numpy``)
- Matplotlib (>= 2.0.2) (``$ pip install matplotlib``)
- Termcolor (>= 1.1.0) (``$ pip install termcolor``)


Installation
--------------

Python package
^^^^^^^^^^^^^^^^^^^

miComplete can easily be installed along with all python dependencies::

   $ pip install micomplete

Assuming that the python bin is in your ``$PATH``, can then be run as::

   $ miComplete

Git
^^^^^^^^^^^^^^^^^^^
1. Choose an appropriate location, e.g. your home::

   $ cd $HOME
   
2. Clone the latest version of the repository::
   
   $ git clone http://bitbucket.org/evolegiolab/micomplete.git

3. Create symlink to some directory in your ``$PATH`` (in this example ``$HOME/bin``)::

   $ cd micomplete
   $ ls micomplete
   $ ln -s $(realpath micomplete/micomplete.py $HOME/bin/miComplete)
   
3. Optionally, add the folder ``micomplete`` in your ``PATH``. The scripts should be kept at their original location.

Usage
--------------

Positional arguments
^^^^^^^^^^^^^^^^^^^^^^^

   A file of sequence(s) along with type (fna, faa, gbk) provided in a tabular format

The file has to contain per line both a path (relative or absolute) to a genomic file as well as the type (separated by a tab)::

   /seq/genomic_sources/e_coli.fna   fna
   /seq/genomic_sources/l_pneumophila.gbk   gbk
   (...)

Optional arguments
^^^^^^^^^^^^^^^^^^^^^^^^

   -h, --help          show help message and exit
   -c, --completeness  Do completeness check (also requires a set of HMMs to have been provided)
   --hlist             Write list of Present, Absent and Duplicated markers for each organism to file
   --hmms HMMS         Specifies a set of HMMs to be used for completeness check or linkage analysis
   --weights WEIGHTS   Specify a set of weights for the HMMs specified, (optional)
   --linkage           Specifies that the provided sequences should be used to calculate the weights of the provided HMMs
   --evalue EVALUE     Specify e-value cutoff to be used for completeness check, default=1e-10
   --cutoff CUTOFF     Specify cutoff percentage of markers required to be present in genome for it be included in linkage calculat. Default = 0.9
   --threads THREADS   Specify number of threads to be used in parallel
   --log LOG           Log name (default=miComplete.log)
   -v, --verbose       Enable verbose logging
   --debug             Debug mode
   
Examples
^^^^^^^^^^^^^^^^^^^^^^^^

Sequence tab file, test_set.tab::

   test_set_common_fna/klebsiella_pneumoniae.fna   fna
   test_set_common_fna/pseudonomonas_aeruginosa.fna        fna
   test_set_common_fna/escherichia_coli.fna        fna
   test_set_common_fna/salmonella_enterica.fna     fna
   
Example 1 - Basic stats
""""""""""""""""""""""""

This example merely produces basic information about the given sequences::

   $ miComplete test_set.tab
   Name	Length	GC-content	N50	L50	N90	L90
   klebsiella_pneumoniae	5682322	57.12	5333942	1	5333942	1
   pseudonomonas_aeruginosa	6264404	66.56	6264404	1	6264404	1
   escherichia_coli	4641652	50.79	4641652	1	4641652	1
   salmonella_enterica	5133713	51.87	4809037	1	4809037	1
   
miComplete prints result to stdout in tabular format, this can favourably be redirected towards a file with a pipe and examined with spreadsheet reader. ::

   $ miComplete test_set.tab > results.tab

Example 2 - Completeness
""""""""""""""""""""""""

This example will produce the same basic statistics, but also completeness and redundance::

   $ miComplete test_set.tab -c --hmms share/Bact139.hmm
   Name	Length	GC-content	Present Markers	Completeness	Redundance	N50	L50	N90	L90
   escherichia_coli	4641652	50.79	139	1.000	1.000	4641652	1	4641652	1
   salmonella_enterica	5133713	51.87	138	0.993	1.000	4809037	1	4809037	1
   klebsiella_pneumoniae	5682322	57.12	136	0.978	1.000	5333942	1	5333942	1
   pseudonomonas_aeruginosa	6264404	66.56	135	0.971	1.000	6264404   1	6264404	1

That is great, but the run time is starting to increase significantly since we have to translate four genomes to proteomes. 
We can speed up the process by running all four parallel with ``--threads``::

   $ miComplete test_set.tab -c --hmms share/Bact139.hmm --threads 4 > results.tab
   
Example 3 - Weighted completeness
""""""""""""""""""""""""""""""""""

This example will also produce the weighted completeness::

   $ miComplete test_set.tab -c --hmms share/Bact139.hmm --weights share/Bact139.weights --threads 4
   Name	Length	GC-content	Present Markers	Completeness	Redundance	CompletenessW	RedundanceW	N50	L50	N90	L90
   escherichia_coli	4641652	50.79	139	1.000	1.000	1.000	1.000	4641652	1	4641652	1
   salmonella_enterica	5133713	51.87	138	0.993	1.000	0.991	1.000	4809037	1	4809037	1
   klebsiella_pneumoniae	5682322	57.12	136	0.978	1.000	0.982	1.000	5333942	1	5333942	1
   pseudonomonas_aeruginosa	6264404	66.56	135	0.971	1.000	0.965	1.000	6264404	1	6264404	1

Example 4 - Creating weights
""""""""""""""""""""""""""""

Finally we will create our own set of weights given a set of marker genes for which we do not already have weights::

   $ miComplete test_set.tab -c --hmms share/Bact109.hmm --linkage --threads 4 > Bact109.weights

Also produces a box plot of the distribution of weights for each marker gene.

