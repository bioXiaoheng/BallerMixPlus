# BalLeRMixPlus
## ---*Bal*ancing selection *L*ik*e*lihood *R*atio *Mix*ture models

This repository hosts the software package for BalLeRMixPlus, which is an extension of [BalLeRMix](https://github.com/bioXiaoheng/BalLeRMix) that can jointly detect recent positive selection and long-term balancing selection.  

------

In BalLeRMixPlus, we introduce the optional `--findPos` and `--findBal`, as well as `--fixAlpha <abeta>`, arguments to help specify the disperson parameter in beta-binomial distributions, with a value between 0 and 1 indicating positive selection, and a value larger than 1 indicating balancing selection. 

Note that BalLeRMixPlus does not consider multi-allelic balancing selection.

------

## Quick Guide 

```
usage: BalLeRMixPlus.py [-h]  -i INFILE [-o OUTFILE] --spect SPECTFILE [--getSpect] [--getConfig] 
                        [--nofreq] [--nosub] [--MAF] [--findBal] [--findPos] 
                        [--physPos] [--rec RRATE] [--fixSize] [-w W] [--noCenter] [-s STEP] 
                        [--fixX X] [--fixAlpha ABETA] [--rangeA SEQA] [--listA LISTA]
                       
```
You can use `python BalLeRMixPlus.py -h` to see the more detailed help page.

### 1. Input format
For B<sub>0</sub> and B<sub>2</sub> statistics, the user should first generate the __*tab-delimited*__ site frequency spectrum file, __*without header*__, *e.g.*:
> \<k\>|\<sample size n\>|\<proportion in the genome\>    
> :-----:|:-----:|-----
> 1|50|0.03572
> 2|50|0.02024
> ...
  
or the configuration file with polymorphism/substitution ratio, __*without header*__, *e.g.*:

> \<sample size n\> | \<\% of substitutions\> | \<\% of polymorphisms\>  
> :-----:|:-----:|:-----:   
> 50  |  0.7346  |  0.2654    

The input files should have four columns, presenting physical positions, genetic positions, number of derived (or minor) alleles observed, and total number of alleles observed (*i.e.* sample size). This file should be tab-delimited and should have a header, *e.g.*:

> physPos|genPos|x|n    
> :-----:|:-----:|:-----:|:-----:    
> 16|0.000016|50|50    
> 35|0.000035|12|50   
> ...
  
### 2. Running the *B* statistics
To perform B<sub>2</sub> scans on your input data, use

    python BalLeRMixPlus.py -i <input> --spect <derived allele frequency spectrum> -o <output>

To perform B<sub>2,MAF</sub> scans on your input data, use

    python BalLeRMixPlus.py -i <input> --spect <minor allele frequency spectrum> -o <output> --MAF

To perform B<sub>1</sub> scans on your input data, use

    python BalLeRMixPlus.py -i <input> --config <sub/poly configuration file> -o <output> --nofreq

To perform B<sub>0</sub> scans on your input data, use

    python BalLeRMixPlus.py -i <input> --config <derived allele frequency spectrum> -o <output> --nosub 

To perform B<sub>0,MAF</sub> scans on your input data, use

    python BalLeRMixPlus.py -i <input> --config <minor allele frequency spectrum> -o <output> --nosub --MAF

### 3. Generate helper files
To generate spectrum file for B<sub>2</sub>:

    python BalLeRMix.py -i <concatenated input> --getSpect --spect <spectrum file name>

To generate spectrum file for B<sub>2,MAF</sub>:

    python BalLeRMix.py -i <concatenated input> --getSpect --MAF --spect <spectrum file name>

To generate spectrum file for B<sub>1</sub>:

    python BalLeRMix.py -i <concatenated input> --getConfig --spect <config file name>

To generate spectrum file for B<sub>0</sub>:

    python BalLeRMix.py -i <concatenated input> --getSpect --nosbub --spect <spectrum file name>

To generate spectrum file for B<sub>0,MAF</sub>:

    python BalLeRMix.py -i <concatenated input> --getSpect --nosub --MAF --spect <spectrum file name>

### 4. Customizing the scan
All arguments besides the aforementioned ones are for customizing the scan. They are not necessary for the scan to complete.

- `[--physPos] [--rec RRATE] `:

   Because `BalLeRMix` uses genetic distances (in cM) to compute likelihood, to direct the software to use physical positions instead, you should use `--physPos`, and indicate the uniform recombination rate (cM/nt) in your species of interest with `--rec`. The default value is 10<sup>-6</sup> cM/nt.
   
   This argument will be automatically incurred if you choose to fix the window size (*e.g.*, 1000bp, 5kb, *etc.* ), in which case yuou want to make sure the software is correctly informed of the recombination rate. Using physical positions will also change how you define window sizes and step sizes, if you were to customize the scanning window.

- ` [--fixX X] [--rangeA SEQA] [--listA LISTA]`:

   These areguments allow you to specify the parameter space that the software optimizes over. The presumed equilibrium frequency is *x*, and the rate of decay in linkage disequilibrium is *A*. If you choose to look for multi-allelic balancing selection where more than two alleles are being balanced, *x* should be a vector of descending equilibrium frequencies, and should match the number of balanced alleles you chose (via `-m`) to scan for.

- ` [--fixSize] [-w R] [--noCenter] [-s STEP] [--physPos]`:

   These areguments are for customizing the scanning window. You probably won't need them because `BalLeRMix` is robust to window sizes. For more details on how these arguments work, check the v1 software manual.
