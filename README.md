# GOFindBias: Analysis tool for finding bias in the GAF files.
GOFindBias is developed to provide the user with some insightful statistics about the [GAF](http://www.geneontology.org/page/go-annotation-file-formats) file to determine if the conclusions on the gene ontology studies can be biased because of abstract terms or the high throughput experiments([1](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003063)).

Statistics Provided by the tool are as follows:
1. Shannon's equitability.
2. Top 'n' PubMed and GO terms.
3. t test to compare two different [GAF](http://www.geneontology.org/page/go-annotation-file-formats) files.


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
### NOTE

If you are using Anaconda environment, make sure that the Python reads libraries from "~anaconda2/lib/python2.7/site-packages/lib". You can also simply copy the files in that location to an appropriate path where other python libraries are readable (importable) by Python.