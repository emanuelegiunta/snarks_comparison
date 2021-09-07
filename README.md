# zk-SNARKs argument size comparison

The repository provides tools to evaluate the argument size of three (at the moment) IOP-based zkSNARKs, Ligero [[AHIV17]](https://acmccs.github.io/papers/p2087-amesA.pdf), Aurora [[BCRSVW19]](https://eprint.iacr.org/2018/828) and Ligero++ [[BFHVXZ20]](https://dl.acm.org/doi/abs/10.1145/3372297.3417893), when applied to a binary R1CS, comparing in each case the standard protocol with the optimised version obtained through a RMFE.
In the case of Ligero, further comparisons are provided with BooLigero [[GSV21]](https://eprint.iacr.org/2021/121), a different way to optimise the same scheme for binary constraints systems.


### Installation and usage
The requirements are any version of Python3 with matplotlib and numpy, which can be installed throught the commands
```
    $ pip3 install numpy
    $ pip3 install matplotlib
```
To compute a full comparison execute
```
    $ python3 optimiser.py -a -v
    $ python3 optimiser.py -l -v
    $ python3 optimiser.py -lpp -v
```
The first one will produce a comparison using Aurora over an extension field of F2 of dimension 198 and, in the optimised case, an RMFE of parameters (48, 198). The second one instead will compare Ligero over arbitrary (small) fields with interactive repetitions using in the optimised case an RMFE of parameters (48, 160).
For a quick test try the following commands
```
    $ python3 optimiser.py -a -ld 6 -hd 8 -v
    $ python3 optimiser.py -l -ld 6 -hd 8 -v
    $ python3 optimiser.py -lpp -ld 6 -hd 8 -v
```
which instruct the program to perform a comparison with respect to Aurora/Ligero/Ligero++ with 2^8 constraints printing some information on screen. By replacing `-a` with `-ta` instead of producing a `.csv` file it will print on screen the parameters that minimises costs in each case. For more informations on the optional flags
```
    $ python3 optimiser.py -help
    $ python3 optimiser.py -help | less -r
```

### Other utilities
In orther to display the data produced and quickly test the improvements, two additional utilities are provided. 

* `csv_plotter.py` plots the data contained in the `.csv` files it receives. It is not meant to be general purpose and can handle well only files produced by optimiser.py. Moreover passing non-csv files will cause the program to raise an input error.
If the data is stored in the same folder of the program, you can run
```
    $ python3 csv_plotter.py *.csv
```
* `csv_analiser.py` prints on screen some basic information on the improvement factor for the data in the .csv files provided. If the data contains comparisons with more than one RMFE, it will print the ratios between the standard version and each of the optimised ones in the order given by the .csv file.
If the data is stored in the same folder of the program, you can run
```
    $ python3 csv_analyser.py *.csv
```


