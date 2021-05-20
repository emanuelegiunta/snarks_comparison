# zk-SNARKs Argument size comparison

The repository provides tool to evaluate the argument size of two (at the moment) IOP-based zkSNARKs Ligero [[AHIV17]](https://acmccs.github.io/papers/p2087-amesA.pdf) and Aurora [[BCRSVW19]](https://eprint.iacr.org/2018/828) when applied to a R1CS of prescribed size over the binary field, and compare the result with the same scheme aplying several optimisation.

### Installation and usage
The only requirement is any version of Python2 or Python3. To compute a full comparison execute
```
    $ python optimiser.py -a -v
    $ python optimiser.py -l -v
```
The first one will produce a comparison using Aurora over an extension field of F2 of dimension 198 and, in the optimised case, an RMFE of parameters (48, 198). The second one instead will compare ligero over arbitrary fields with interactive repetitions and, in the optimised case, an RMFE of parameters (48, 160).
For a quick test try the following commands
```
    $ python optimiser.py -a -ld 6 -hd 8 -v
    $ python optimiser.py -l -ld 6 -hd 8 -v
```
Which instructs the program to perform a comparison with respect to Aurora/Ligero for the number of constranints from 2^6 to 2^8 (inclusive) and turn the verbose flag on. By `-a` with `-ta` instead of producing a .csv file will print on screen the parameters that minimises costs in each case. For more informations on the optional flags
```
    $ python optimiser.py -help
```

### Other utilities
In orther to display the data produced and quickly test the improvement two additional utilities are provided. 

* `csv_plotter.py` plot the data contained in the .csv files it receives. It is not meant to be general purpose and can handle well only files produced by optimiser.py. Moreover passing non .csv files will cause the program to raise an input error.
If the data is stored in the same folder of the program, you can run
```
    $ python csv_plotter.py *.csv
```
* `csv_analiser.py` print on screen some basic information on the improvement factor for the data in the .csv files provided. If the data contains comparisons with more than one RMFE, it will print the ratios between the standard version and each the optimised RMFE in the order given by the .csv file.
If the data is stored in the same folder of the program, you can run
```
    $ python csv_analyser.py *.csv
```


