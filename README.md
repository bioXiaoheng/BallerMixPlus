# BalLeRMix+
## ---*Bal*ancing selection *L*ik*e*lihood *R*atio *Mix*ture models Plus

This repository hosts the software package for BalLeRMix+, an extension of [BalLeRMix](https://github.com/bioXiaoheng/BalLeRMix) that can jointly detect recent positive selection and long-term balancing selection.  

------

In BalLeRMix+, we introduce the optional `--findPos` and `--findBal`, as well as `--fixAlpha <abeta>`, arguments to help specify the disperson parameter in beta-binomial distributions, with a value between 0 and 1 indicating positive selection, and a value larger than 1 indicating balancing selection. Based on feedback from BalLeRMix users, we also added a `--minCount` argument to indicate the smallest number of allele counts in the input data, in case the user chooses to remove rare variants from the data.

To help visualize the output, we also added an R script in the `test/` folder in case you need. 

Note that BalLeRMix+ does not consider multi-allelic balancing selection as does BalLeRMix_v2.2.

------

## User Guide 

```
usage: BalLeRMix+_v1.py [-h]  -i INFILE [-o OUTFILE] --spect SPECTFILE [--minCount MINCOUNT] [--getSpect] [--getConfig] 
                        [--noFreq] [--MAF] [--findBal] [--findPos] [--usePhysPos] [--rec RRATE]
                        [--fixWinSize] [-w W] [--noCenter] [-s STEP] [--fixX X] [--fixAlpha ABETA] [--rangeA SEQA] [--listA LISTA]
                       
```
You can use `python BalLeRMix+_v1.py -h` to see a more detailed help page.

### 0. Installation and dependency
To install `BalLeRMix+`, simply navigate to your working directory and clone this repository:

     git clone https://github.com/bioXiaoheng/BallerMixPlus.git

The current release supports Python version 3.6 and above, and needs to have numpy (>=1.19.1) and scipy (>=1.5.0) installed. 
<details open>
 <summary>Commands to install requirements:</summary>
 
 ```
 # navigate to the cloned repository
 cd BallerMixPlus/
  
 # check for required packages
 pip install -r requirements.txt
 ```
  
</details>

### 1. Quick start and examples
To run the program, the user can build commands based on `python BalLeRMix+_v1.py -i <input file> [other arguments]`, with no particular order required for each argument. In particular, useful argument combinations are

| *B* variant | Input data type | To generate helper file | To run scan | 
|:-----------:|:---------------:|:-----------------------:|:-----------:|
| *B<sub>2</sub>* | Derived allele frequency (DAF) plus substitutions (x==n)| `--getSpect` `--spect <spect file>` | `-o <output file>` |
| *B<sub>2,MAF</sub>* | Minor allele frequency (MAF) plus substitutions (x==0)| `--getSpect` `--MAF``--spect <spect file>` | `--MAF`  `-o <output file>` |
| *B<sub>1</sub>* | Polymorphism (0<x<n) or substitution (x==n) calls | `--getConfig` `--spect <config file>` | `--noFreq` `-o <output file>` |
| *B<sub>0</sub>* | Derived allele frequency (DAF)| `--getSpect` `--noSub`  `--spect <spect file>` | `--noSub` `-o <output file>` |
| *B<sub>0,MAF</sub>* | Minor allele frequency (MAF)| `--getSpect` `--noSub` `--MAF`  `--spect <spect file>` | `--noSub` `--MAF`  `-o <output file>` |

Note that when generating helper files, the input must be concatenated from all sequences in the genome. One command-line option to concatenate them is (suppose your input files are named `input_chr<n>.txt`):
```
cat input_chr*.txt | awk '(NR == 1 || $1 != "position")' > Concatenated_input_for_helper_file.txt
```

There are two examples in `test/` folder, including their input files as well as some helper files. To obtain helper files from the concatenated input (`test/HC_CEU_Neut_Concatenated-DAF.txt`), you can run:
<details open>
  <summary>Click to check sample commands</summary>
  
```
# generate helper file for B1
python BalLeRMix+_v1.py -i test/HC_CEU_Neut_Concatenated-DAF.txt --getConfig --spect test_config.txt

# generate helper file for B2
python BalLeRMix+_v1.py -i test/HC_CEU_Neut_Concatenated-DAF.txt --getSpect --spect test_DFS.txt

# generate helper file for B2maf
python BalLeRMix+_v1.py -i test/HC_CEU_Neut_Concatenated-DAF.txt --getSpect --MAF --spect test_MFS.txt

# generate helper file for B0
python BalLeRMix+_v1.py -i test/HC_CEU_Neut_Concatenated-DAF.txt --getSpect --noSub --spect test_DFS-noSub.txt

# generate helper file for B0maf
python BalLeRMix+_v1.py -i test/HC_CEU_Neut_Concatenated-DAF.txt --getSpect --noSub --MAF --spect test_MFS-noSub.txt
```
  
</details>

The resulting files should be the same as `test/HC_CEU_Neut_config_for_B1.txt`, `test/HC_CEU_Neut_DAF_spect_for_B2.txt`, `test/HC_CEU_Neut_MAF_spect_for_B2maf.txt`, `test/HC_CEU_Neut_DAF-noSub_spect_for_B0.txt`, and `test/HC_CEU_Neut_MAF-noSub_spect_for_B0maf.txt`, respectively.

To run B<sub>1</sub>, B<sub>2</sub>, and B<sub>2,MAF</sub>, respectively, on the first example:
```
#run B1
python BalLeRMix+_v1.py -i test/Example1_fullSweep_200kya_DAF.txt -o testout_ex1_B1.txt --noFreq --spect test/HC_CEU_Neut_config_for_B1.txt

# run B2
python BalLeRMix+_v1.py -i test/Example1_fullSweep_200kya_DAF.txt -o testout_ex1_B2.txt --spect test/HC_CEU_Neut_DAF_spect_for_B2.txt

# run B2maf
python BalLeRMix+_v1.py -i test/Example1_fullSweep_200kya_MAF.txt -o testout_ex1_B2maf.txt --MAF --spect test/HC_CEU_Neut_MAF_spect_for_B2maf.txt

```
You will find that the output should be close, if not identical, to the files `test/output/Example1_B1.txt`, `test/output/Example1_B2.txt`, and `test/output/Example1_B2maf.txt`, respectively. Example1 is simulated to have undergone a selective sweep of *s*=0.01. Because *B<sub>0</sub>* and *B<sub>0,MAF</sub>* have little power for positive selection, we do not recommend the user to use them for inferring positive selection. 

To visualize the *B<sub>1</sub>*, *B<sub>2</sub>*, and *B<sub>2,MAF</sub>* scores, you can try:
```
Rscript test/plotScores.r <output from BalLeRMix+_v1.py> <image name with suffix>
```
You should be able to generate images resembling `test/ScorePlot_example1_B2.png`. Note that this R script requires the `ggplot2` package.

The same operations can be repeated on example 2, and you can check your output with `test/output/Example2_B1.txt`, `test/output/Example2_B2.txt`, and `test/output/Example2_B2maf.txt` accordingly.

By default, the `BalLeRMix+_v1.py` script uses all informative sites provided in the input. With additional arguments, users can customize the size and centers of their sliding windows for the analyses. Details on these areguments are listed in section 4 As an example, one can run
```
# compute B0maf (--noSub --MAF) while using physical positions as coordinates (--usePhysPos), scan with sliding windows of a fixed physical size (--fixWinSize) of 1000 nt (-w 1000), centered on every other (--step 2) informative sites
python BalLeRMix+_v1.py --noSub --MAF --spect test/HC_CEU_Neut_MAF-nosub_spect_for_B0maf.txt -i test/Example2_balancing_10MYA_MAF_nosub.txt -o testout_B0maf_1kb-500b.txt --usePhysPos --fixWinSize -w 1000 --step 2  
```
The output should be close, if not identical, to `test/output/Example2_B0maf_1kb-2sites.txt`.

### 2. Input and helper files
<details open>
<summary></summary>
    The input file should have four columns, presenting physical positions, genetic positions, number of derived (or minor) alleles observed, and total number of alleles observed (*i.e.* sample size). We ask that a sample size be provided for each locus so that loci with missing data (and hence a difference sample size) can be better accounted for. Meanwhile, all (x,n) combinations in your input must also exist in the helper file. 
  
  The input file should be tab-delimited and should have a header, *e.g.*: 

  > physPos|genPos|x|n    
  > :-----:|:-----:|:-----:|:-----:    
  > 16|0.000016|50|50    
  > 35|0.000035|12|50   
  > ...

  Note that for B<sub>2,MAF</sub>, it is recommended to only have minor allele counts in the data. If, however, the user choose to run it on derived allele counts, as long as `--MAF` argument is given, the script will automatically fold the allele count in your data. For B<sub>1</sub>, the script will read all loci with allele counts different from their sample sizes (*i.e.*, x != n ) as polymorphisms, and all that has identical values to the samples sizes as substitution (also known as fixed differences).

  In addition to the allele count input file, the user also need to provide one helper file that records the genome-wide variation level. In particular:

  For B<sub>2</sub> and B<sub>2,MAF</sub> statistics, a site frequency spectrum file is needed, and should be __*tab-delimited*__ and __*without header*__, *e.g.*:
  > \<k\>|\<sample size n\>|\<proportion in the genome\>    
  > :-----:|:-----:|-----
  > 1|50|0.03572
  > 2|50|0.02024
  > ...

  For B<sub>1</sub>, the helper file (also tab-delimited) records the genome-wide polymorphism/substitution ratio, __*without header*__, *e.g.*:

  > \<sample size n\> | \<\% of substitutions\> | \<\% of polymorphisms\>  
  > :-----:|:-----:|:-----:   
  > 50  |  0.7346  |  0.2654    

  The `BalLeRMix+_v1.py` program can help generate helper files from the concatenated input files. 
  For site frequency spectrum: 

      python BalLeRMix+_v1.py -i <concatenated input file> --getSpect --spect <spectrum file name>

  For minor allele frequency spectrum:

      python BalLeRMix+_v1.py -i <concatenated input file> --getSpect --MAF --spect <spectrum file name>

  For polymorphism-substitution configuration file:

      python BalLeRMix+_v1.py -i <concatenated input file> --getConfig --spect <config file name>

</details>

### 3. Running the *B* statistics 
<details open>
<summary></summary>
To perform B<sub>2</sub> scans on your input data, use

    python BalLeRMix+_v1.py -i <input> --spect <derived allele frequency spectrum> -o <output>

To perform B<sub>2,MAF</sub> scans on your input data, use

    python BalLeRMix+_v1.py -i <input> --spect <minor allele frequency spectrum> -o <output> --MAF

To perform B<sub>1</sub> scans on your input data, use

    python BalLeRMix+_v1.py -i <input> --config <sub/poly configuration file> -o <output> --noFreq

</details>
  
### 4. Customizing the sliding window for your scan
  
All arguments besides the aforementioned ones are for customizing the scan. Their usage is largely based on the previous `BalLeRMix` (find out more in its detailed manual [here](https://github.com/bioXiaoheng/BalLeRMix/blob/master/software/BalLeRMix_manual.pdf) ). Note that these arguments are not necessary for the scan to complete. 
  <details open>
    <summary>In general, they can be grouped into the following three categories:</summary>

- `[--usePhysPos] [--rec RRATE] `:

   Because `BalLeRMix` uses genetic distances (in cM) to compute likelihood, to direct the software to use physical positions instead, you should use `--physPos`, and indicate the uniform recombination rate (cM/nt) in your species of interest with `--rec`. The default value is 10<sup>-6</sup> cM/nt.
   
   This argument will be automatically incurred if you choose to fix the window size (*e.g.*, 1000bp, 5kb, *etc.* ), in which case yuou want to make sure the software is correctly informed of the recombination rate. Using physical positions will also change how you define window sizes and step sizes, if you were to customize the scanning window.

- ` [--fixX X] [--fixAlpha ABETA] [--rangeA SEQA] [--listA LISTA]`:

   These areguments allow you to customize the parameter space that the software optimizes over. The presumed equilibrium frequency is *x*, the extent of dispersion is *α*, and the rate of decay in linkage disequilibrium is *A*. 

- ` [--fixWinSize] [-w R] [--noCenter] [-s STEP] [--usePhysPos]`:

   These areguments are for customizing the scanning window. You probably won't need them because `BalLeRMix` and `BalLeRMix+` are robust to window sizes. The software help page (`python BalLeRMix+_v1.py --help`) provides more details on how these arguments work. The table below lists combinations the users can go with to customize their scans:
  
 |   Window Span             |  Position tested   |          Required Commands          |    Optional Commands  | 
 |:-------------------------:|:------------------:|:-----------------------------------:|:---------------------:|
 |Fixed sequence length (nt) | Informative sites  |`--fixWinSize` `-w <nt>` `[-s <sites>]` `[--usePhysPos]`| `[--rec <cM/nt>]`     |
 |Fixed sequence length (nt) |Evenly-spaced arbitrary positions|`--fixWinSize` `--noCenter``-w <nt>` `-s <nt>` `[--usePhysPos]`| `[--rec <cM/nt>]`     |
 |Fixed number of informative sites| Informative sites  |`-w <sites in radius>` `-s <sites per step>`| `[--usePhysPos]` `[--rec <cM/nt>]`|

</details>  
