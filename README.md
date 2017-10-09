# GOFindBias: Analysis tool for finding bias in the GAF files.
GOFindBias is developed to provide the user with some insightful statistics about the [GAF](http://www.geneontology.org/page/go-annotation-file-formats) file to determine if the conclusions on the gene ontology studies can be biased because of abstract terms or the high throughput experiments([1](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003063)). This tool selects and gives statistics about only entries having 'EXP','IDA','IPI','IMP','IGI','IEP' evidence code since they are found to be most reliable one[2](http://www.geneontology.org/book/export/html/799).


Statistics Provided by the tool are as follows:
1. Shannon's equitability.
2. Top 'n' PubMed and GO terms.
3. Statistical test to compare two different [GAF](http://www.geneontology.org/page/go-annotation-file-formats) files.


### Prerequisites:
#### Required modules. 

Modules are available in most GNU/Linux distributions, or from their respective websites.

* [Matplotib](https://matplotlib.org/)

* [Biopython](http://biopython.org/)

### Installation

Installing from source
```
git clone https://github.com/Pranavkhade/GOFindBias
cd GoFindBias
python setup.py install
```

Installing with pip
```
pip install gofindbias
```
OR
```
pip install git+git://github.com/Pranavkhade/GOFindBias
```

### Files and instructions

1. Collect the .gaf file you wish to analyse. For reference .gaf files you can visit ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/
2. You can also use the .gaf files obtained as an output from the [debias](https://github.com/Rinoahu/debias) program.
3. Instructions and help for the parameter is as follows:

```
positional arguments:
  NameOfTheGAFFile  Please provide the name of the .gaf file
  Log10[1 or 0]     Parse 1 if you want the scale to be Log to the base 10
  TopStatistics     How many top entires you want to include?

optional arguments:
  -h, --help        show this help message and exit
```
### Examples

1. `GOFindBias test/2014.gaf 1 10`
This command will parse 2014.gaf file for the analysis and all the GO Term counts will be represented on the Natural Log scale for better comparitive visualisation of the data. The last argument is the number of top 'n' entries with highest count in the GAF file. The output of graphs will be posted in the `/graph_output` folder with names corrosponding to GO/PubMed ID count and the ontology level (F/C/P). File named `Shannon's Statistics.txt` will have the information about the diversity of a given .gaf file.

2. `GOFindBias test/2016.gaf 0 50`
This will give the exact same output files but for a file named 2016.gaf. All the graph visualisation will not be in the log scale and the graphs will have information about the top 50 entires. 

### NOTE

1. If you are using Anaconda environment, make sure that the Python reads libraries from "~anaconda2/lib/python2.7/site-packages/lib". You can also simply copy the files in that location to an appropriate path where other python libraries are readable (importable) by Python.
2. You can find few test .gaf files in the /lib/test/ folders.
