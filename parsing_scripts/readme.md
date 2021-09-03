### Parsing empirical data for [BalLeRMix](https://github.com/bioXiaoheng/BalLeRMix) or BalLeRMix+ application

This folder contains example scripts for parsing empirical data into formats fit for BalLeRMix and BalLeRMix+. Note that all scripts in this directory are mostly meant to examplify typical workflows and have <u>only been tested with human and great ape data</u> so far, so *please excercise caution* before using them.

#### *What you need:* VCFs; *Bonus:* AXT and/or recombination map
 
 The must-have for parsing an input for BalLeRMix and BalLeRMix+ is variant calls from your samples, ideally covering the whole genome. One most commonly-used filetype for this is the [VCF](https://faculty.washington.edu/browning/beagle/intro-to-vcf.html) file. 

 To apply B<sub>1</sub>, B<sub>2</sub>, and B<sub>2,MAF</sub> on your data, you must include substitutions (also known as \"fixed differences\") in your input. To find out these sites, you can compare your samples with one genome sequence from another species closely-related to your species-of-interest, and designate the outgroup sequence as the ancestral sequence. To this end, pairwise alignment [AXT](https://genomebrowser.wustl.edu/goldenPath/help/axt.html) files can be very useful. [UCSC Genome Browser](https://hgdownload.soe.ucsc.edu/downloads.html) hosts such data for many organisms sequenced so far. 

 For species blessed with recombination maps, it is highly recommend to provide BalLeRMix and BalLeRMix+ the corresponding genetic position (in cM) of your SNPs to properly account for uneven recombination landscape. These maps comes in many formats, among which a relatively commonly-used one is a tab-delimited four-column plain text file that shows chromosome, physical position, recombination rate (often in cM/Mb), and cumulative map position in cM. You can find some examples [here](https://github.com/cflerin/dog_recombination) and [here](https://github.com/popgenmethods/pyrho#make_table).

#### Using the `parse_ballermix_input.py` script

 For the convenience of potential users, we provide `parse_ballermix_input.py` script to (jointly) parse empirical data from common formats to BalLeRMix-ready format. You need Python3.8 or above to run it. Again, note that this script has <u>only been tested on human and great ape data</u> so far, so use caution before applying this on your data.

 Suppose you already cloned the BallerMixPlus repository, check out the help page of this script by using

```Bash
 cd parsing_scripts/

 #check for requirements
 pip install -r requirements.txt

 #see help page
 python parse_ballermix_input.py -h

```
 You will see

```Bash
usage: parse_ballermix_input.py [-h] --vcf VCFFILE --ID_list POP_LIST -c CH [--axt AXTFILE] [--rec REC_RATE] [--rec_map REC_MAP] -o OUTFILE

optional arguments:
  -h, --help            show this help message and exit
  --vcf VCFFILE         Path and name of the vcf file.
  --ID_list POP_LIST    Path and name to the file containing list of sample IDs (identical to their column names in vcf) to be counted,
                        separated by comma. If not provided, all samples in the vcf will be counted.
  -c CH, --chr CH       ID of the chromosome. E.g. 2a for chr2a, 12 for chr12, etc.
  --axt AXTFILE         Path and name of the sequence alignment (in .axt or .axt.gz) for calling substitution and polarizing the allele
                        frequency. If not provided, then the output will only be applicable to B_0maf.
  --rec REC_RATE        Recombination rate in cM/nt. Default value is 1e-6 cM/nt.
  --rec_map REC_MAP     Path and name of the recombination map (hapmap format) of the same sequence. If not provided, a uniform recombination
                        rate will be applied with a default rate of 1e-6 cM/nt. Use "--rec" to specify another rate.
  -o OUTFILE, --output OUTFILE
                        Path and name of the output file.
```


 In the example shown here, we are using the chimpanzee genome panTro6 as a reference, and use