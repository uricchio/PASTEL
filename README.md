# PASTEL

You can install our software by running:

python setup.py install

If you don't have admin access on the machine you can do:

python setup.py install --user

The software assumes you have numpy installed.

I've also attached a sample script and a couple sample data files.  The script is called "runPerms.py", and can be run with:

python runPerms.py <input_summary_stat_file> <lower frequency threshold> <upper frequency threshold>

For example, to run it on the UK Biobank test file I included and use all alleles, it would be:

python runPerms.py sumstat_examples/body_HEIGHTz.sumstats.txt 0 1

The input file formats look like:
 
```
rsnum chromosome position effect_allele alt_allele ancestral_allele effect_allele_frequency beta INFO block_number pval 
rs28863004 1 526736 C G G 0.993967 -0.00127255 0.396081 0 7.7E-01
rs78497331 1 533198 C T C 0.998983 0.0294703 0.532402 0 5.3E-01
rs576404767 1 544584 C T C 0.99811 0.0214212 0.492517 0 4.4E-01
```

The last colum (pvalue) is optional.  You can condition on pvalue or INFO (i.e. imputation quality).  Note that in the GIANT example file all the INFO values are set to 1 since we didn't have those data immediately available.

Ancestral allele calls were mapped from the Thousand Genomes Project. The "block numbers" (which correspond to approximately independent blocks of the genome) from this paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4731402/.   

The code will output the observed value of S_{\beta} from our manuscript on the first line (after a #), and then the permuted values will follow.  
