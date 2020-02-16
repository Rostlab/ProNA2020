# ProNA2020
ProNA2020: System predicting protein-DNA, protein-RNA and protein-protein binding sites from sequence
ProNA2020 can predict protein-DNA, protein-RNA and protein-protein binding sites with only sequence information. 
It is a two level prediction. In first level, ProNA2020 will predict whether the input protein is a binding protein. 
If the input protein is a binding protein, then ProNA2020 will give the residue level prediction.

## How to install
After downloading, you need install some necessary python package and sofware before you use ProNA2019.

#### Python3 package:

pybrain (>=0.3.3)

scikit-learn (=0.19.1-2)

botocore (>=1.10.33)

smart-open (>=1.5.7)

#### Software:

ncbi-blast+

fastprofkernel (>=1.0.24)

## How to run

Usage: prona2020 [options]

Options:

* -h:show this help message and exit
  
* -p:PATH,Directory containing the PredictProtein output files with the
               suffixes .chk .in .fasta .blastPsiMat .profbval .mdisorder and
               .profRdb.
               
* -o:FILENAME,Output file. If not specified, the output is written to STDOUT.
  
* -l:LABEL,Turn off protein level prediction by inputting binding label,
               e.g. "-l  Protein_DNA", which means the input protein is already known as a Protein- and DNA-binding protein.
               
* -d:DATABASE,Use your own local database for PSI-BLAST (homology based
               inference), default is using the profile (.chk) from big_80
               database(rostlab) which is a comprehensive blast database at
               80% sequence identity redundancy level
               
* -v:Print verbose or not (True/False), default is False

