<pre>
miRMaster is an open source python implementation made available under the GNU General Public License. It requires Python3.5 or above for smooth running. The software screens novel precursor sequences and mature miRNA from EST datasets<br> of plants. There is a huge inflow of data in terms of genomics and molecular biology. There are several softwares for screening miRNA from small RNA-seq data. However, there is a scarcity of working softwares for screening miRNA from EST data. Moreover the ones available are mostly based on homology. Here, we present an ML-based system with Random Forest algorithm named miRMaster for the computational prediction of miRNA from EST datasets.

System requirements:
-Ubuntu (>=16.04)

Dependencies required to run miRMaster:
-Blast executables
-ViennaRNA Package 
-Biopython
-scikit-learn


Python 2.7.9 onwards, and Python 3.4 onwards, include the package management
system "pip" which should allow you to install Biopython (and its dependency
NumPy if needed), upgrade or uninstall with just one terminal command::

    pip install biopython
    pip install --upgrade biopython
For installing scikit-learn modules type 
    pip install -U scikit-learn

Building local server
git clone https://github.com/Garima268/miRMaster.git
cd miRMaster

Once the dependencies are installed please enter the absolute paths of all the above in the config.py 
For example
  #Enter full path to blast executables
  path_blast = "/usr/bin/ncbi-blast-2.10.0+/bin"
  #Enter full path to ViennaRNA package directory
  path_vienna = "/usr/bin/ViennaRNA-2.4.14/"
  #Enter full path to Uniprot/nr database directory
  path_db = "/home/User/miRMaster/"
  #Enter name of Database fasta file
  db_name = "uniprot_sprot.fasta"
  #Enter path to miRMaster directory
  path_tool = "/home/User/miRMaster/"


Usage:
python miRMaster.py <EST.fasta>

</pre>

