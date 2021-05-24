# zk-SNARKs argument size comparison

The repository provides tools to evaluate the argument size of two (at the moment) IOP-based zk-SNARKs, Ligero [[AHIV17]](https://acmccs.github.io/papers/p2087-amesA.pdf) and Aurora [[BCRSVW19]](https://eprint.iacr.org/2018/828), when applied to a binary R1CS, comparing in each case the standard protocol with the optimised version obtained through a RMFE.

### Installation and usage
The only requirement is any version of Python2 or Python3. To compute a full comparison (takes several minutes) execute
```
    $ python optimiser.py -a -v
    $ python optimiser.py -l -v
```
The first one will produce a comparison using Aurora over an extension field of F2 of dimension 198 and, in the optimised case, an RMFE of parameters (48, 198). The second one instead will compare Ligero over arbitrary (small) fields with interactive repetitions using in the optimised case an RMFE of parameters (48, 160).
For a quick test try the following commands
```
    $ python optimiser.py -a -ld 8 -hd 8 -v
    $ python optimiser.py -l -ld 8 -hd 8 -v
```
which instruct the program to perform a comparison with respect to Aurora/Ligero with 2^8 constraints printing some information on screen. By replacing `-a` with `-ta` instead of producing a `.csv` file it will print on screen the parameters that minimises costs in each case. For more informations on the optional flags
```
    $ python optimiser.py -help
    $ python optimiser.py -help | less -r
```

### Other utilities
In orther to display the data produced and quickly test the improvements, two additional utilities are provided. 

* `csv_plotter.py` plots the data contained in the `.csv` files it receives. It is not meant to be general purpose and can handle well only files produced by optimiser.py. Moreover passing non-csv files will cause the program to raise an input error.
If the data is stored in the same folder of the program, you can run
```
    $ python csv_plotter.py *.csv
```
* `csv_analiser.py` prints  on screen some basic information on the improvement factor for the data in the .csv files provided. If the data contains comparisons with more than one RMFE, it will print the ratios between the standard version and each the optimised ones in the order given by the .csv file.
If the data is stored in the same folder of the program, you can run
```
    $ python csv_analyser.py *.csv
```


