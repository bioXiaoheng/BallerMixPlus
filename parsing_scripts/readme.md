## Parsing empirical data for [BalLeRMix](https://github.com/bioXiaoheng/BalLeRMix) or BalLeRMix+ application

This folder contains example scripts for parsing empirical data into formats fit for BalLeRMix and BalLeRMix+. Note that all scripts in this directory are mostly meant to examplify typical workflows and have *only been tested with human and great ape data* so far, so *please excercise caution* before using them.

### 0. *What you need:* VCFs; *Bonus:* AXT and/or recombination map
 
- **VCF**: The must-have for parsing an input for BalLeRMix and BalLeRMix+ is variant calls from your samples, ideally covering the whole genome. One most commonly-used filetype for this is the [VCF](https://faculty.washington.edu/browning/beagle/intro-to-vcf.html) file. 

- **AXT**: To apply B<sub>1</sub>, B<sub>2</sub>, and B<sub>2,MAF</sub> on your data, you must include substitutions (also known as \"fixed differences\") in your input. To find out these sites, you can compare your samples with one genome sequence from another species closely-related to your species-of-interest, and designate the outgroup sequence as the ancestral sequence. To this end, pairwise alignment [AXT](https://genomebrowser.wustl.edu/goldenPath/help/axt.html) files can be very useful. [UCSC Genome Browser](https://hgdownload.soe.ucsc.edu/downloads.html) hosts such data for many organisms sequenced so far. 

- **Recombination map**: For species blessed with recombination maps, it is highly recommend to provide BalLeRMix and BalLeRMix+ the corresponding genetic position (in cM) of your SNPs to properly account for uneven recombination landscape. These maps comes in many formats, among which a relatively commonly-used one is a tab-delimited four-column plain text file that shows chromosome, physical position, recombination rate (often in cM/Mb), and cumulative map position in cM. You can find some examples [here](https://github.com/cflerin/dog_recombination) and [here](https://github.com/popgenmethods/pyrho).

### 1. Using the `parse_ballermix_input.py` script

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

```
usage: parse_ballermix_input.py [-h] --vcf VCFFILE -c CH -o OUTFILE [--ID_list POP_LIST] [--axt AXTFILE] [--rec_rate REC_RATE] [--rec_map REC_MAP]

optional arguments:
  -h, --help            show this help message and exit
  --vcf VCFFILE         Path and name of the vcf file. Format can be either .vcf or .vcf.gz.
  -c CH, --chr CH       ID of the chromosome. E.g. 2a for chr2a, 12 for chr12, etc.
  -o OUTFILE, --output OUTFILE
                        Path and name of the output file.
  --ID_list POP_LIST    Path and name to the file containing list of sample IDs (identical to their column names in vcf) to be counted,
                        separated by comma. If not provided, all samples in the vcf will be counted.
  --axt AXTFILE         Path and name of the sequence alignment (in .axt or .axt.gz) for calling substitution and polarizing the allele
                        frequency. If not provided, then the output will only be applicable to B_0maf.
  --rec_rate REC_RATE        Recombination rate in cM/nt. Default value is 1e-6 cM/nt.
  --rec_map REC_MAP     Path and name of the recombination map (hapmap format) of the same sequence. If not provided, a uniform recombination
                        rate will be applied with a default rate of 1e-6 cM/nt. Use "--rec" to specify another rate.
```
Notes:
* This script is only applicable to *diploid organisms with individual genotype calls*. While we understand that pooled-sequencing data is common for many non-model organisms, the algorithms behind polarizing alleles and matching recombination/genetic position are nonetheless the same. We hope this script can lend insight into the workflow or some of the algoritms, and the users are welcomed to copy or modify this script to best serve your data.
* We understand that sometimes it is procedurally easier to merge samples from multiple populations for variant calling. Therefore, the VCF file here does not have to only contain samples from the population you want to examine. Instead, we invite users to provide a list of samples (matching their column names in the VCF) in your population along with the merged VCF. See `Example3_YRI_samples_1KG-v3.20130502.txt` for its format.
* Note that some VCF files (but not all) contain information about the ancestral alleles ("AA"). This information can also be used to polarize allele frequency, but the current `parse_ballermix_input.py` script does not account for it yet.
* When the AXT file is provided, only sites within sequences mapped in both genomes should be counted. For example, matches in `0 chr1 10918 11034 chrUn_NW_019934090v1 15192 15299 - 6550` or `15208 chr10_gl383543_fix 1 57120 chr10 18044152 18101342 + 5201976` will not be counted.
* The table below shows what *B* statistics can be applied on the parsed data. Note that **input suitable for *B*<sub>2</sub> can be used to compute all other *B* variants too**.
 
 VCF | AXT | Recombination Map | Applicable *B* variant
 :---:|:---:|:---:|:---:
  :heavy_check_mark: |  |  | *B*<sub>0,MAF</sub>
 :heavy_check_mark: |  | :heavy_check_mark: | *B*<sub>0,MAF</sub>
 :heavy_check_mark: | :heavy_check_mark: | | *B*<sub>2</sub>
 :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | *B*<sub>2</sub>

#### 1.1 Walking through Example 3

 In the Example 3 shown here, we will parse the 1000 Genome Project data (downloaded [here](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/)) for YRI on a truncation of chromosome 22, using the chimpanzee genome panTro6 as a reference (pairwise alignment downloaded [here](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/vsPanTro6/)). For recombination map, we use the ones inferred by Spence & Song (2019) (downloaded from [here](https://github.com/popgenmethods/pyrho)).
 
 <details open>
 <summary>Here is how I prepared the downloaded files into these example files.</summary>

```Bash
 # Paths to the files are ommitted for brevity.
 # In reality, for the commands below to work, all files need to be in the same directory
 
 # truncate vcf, keep the first 2000 variant calls
 zcat ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz | head -2252 | gzip > Example3_first2000var.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
 
 # obtain the list of samples in YRI population
 cat integrated_call_samples_v3.20130502.ALL.panel | awk 'BEGIN{ORS = ","}{if ($2 == "YRI") print $1}' | sed 's/,$/\r/'> Example3_YRI_samples_1KG-v3.20130502.txt

 # extract only relevant entries from the alignment file (which is too huge to post on github)
 ## I wrote another short script 'getChrAxt.sh' just for convenience. Please doublecheck whether this script suits your data before using it
 sh getChrAxt.sh hg19.panTro6.net.axt.gz 22 | gzip > Example3_chr22_good.hg19.panTro6.net.axt.gz
 
 # truncate the recombination map accordingly
 ## unzip the downloaded file
 tar -xzf hg19_maps.tar.gz
 ## recall where the vcf ends:
 endpos=$( zcat ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz | cut -f 2 | head -2252 | tail -1 )
 ## truncate
 cat hg19/YRI/YRI_recombination_map_hapmap_format_hg19_chr_22.txt | awk -v e=$endpos '(NR == 1 || $2 <= e)' > Example3_YRI_0-1622e4_rec_map_hapmap_format_hg19_chr22.txt
```
 
 </details>
  
  To prepare these data into BalLeRMix-ready format, use the following command:
  
```Bash
  # if you only have vcf available
  ## indicate the chromosome id with -c/--chr
  ## the vcf contain variant calls for all 2504 samples in 1KG project; the ID_list provides the sample IDs that should be included
  ## assume a uniform recombination rate of 1.25e-6
  python parse_ballermix_input.py --vcf Example3_first2000var.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz -c 22 --ID_list Example3_YRI_samples_1KG-v3.20130502.txt --rec_rate 1.25e-6 -o Example3_vcf-only_rec1.25e-6_b0maf-ready.txt
  
  # if you also have recombination map
  python parse_ballermix_input.py --vcf Example3_first2000var.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz -c 22 --ID_list Example3_YRI_samples_1KG-v3.20130502.txt --rec_map Example3_YRI_0-1622e4_rec_map_hapmap_format_hg19_chr22.txt -o Example3_vcf-plus-recmap_b0maf-ready.txt
  
  # if you don't have rec map, but do have pairwise alignment with a close outgroup (chimpanzee here)
  ## assume a uniform recombination rate of 1.25e-6 cM/nt
  python parse_ballermix_input.py --vcf Example3_first2000var.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz -c 22 --ID_list Example3_YRI_samples_1KG-v3.20130502.txt --axt Example3_chr22_good.hg19.panTro6.net.axt.gz --rec_rate 1.25e-6 -o Example3_vcf-plus-axt_rec1.25e-6_ballermix-ready.txt

  # if you have both rec map and pairwise alignment
  python parse_ballermix_input.py --vcf Example3_first2000var.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz -c 22 --ID_list Example3_YRI_samples_1KG-v3.20130502.txt --axt Example3_chr22_good.hg19.panTro6.net.axt.gz --rec_map Example3_YRI_0-1622e4_rec_map_hapmap_format_hg19_chr22.txt -o Example3_vcf-axt-recmap_ballermix-ready.txt
                                                                                                         
  # note that "--rec_map" and "--rec_rate" are not mutually exclusive. When a map is provided, genomic regions not covered by the map are assumed to have a recombination rate specified by "--rec_rate". If it's not specified, the default is 1e-6
  python parse_ballermix_input.py --vcf Example3_first2000var.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz -c 22 --ID_list Example3_YRI_samples_1KG-v3.20130502.txt --axt Example3_chr22_good.hg19.panTro6.net.axt.gz --rec_map Example3_YRI_0-1622e4_rec_map_hapmap_format_hg19_chr22.txt --rec_rate 1.25e-6 -o Example3_vcf-axt-recmap_rec1.25e-6_ballermix-ready.txt
```
You can compare your output files from these commands with the corresponding files in `test_output/` folder.
