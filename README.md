# ONT-process
based on the R10 chemistry deep sequencing of ONT sequencing data, for genome assembly, QC, and genomic features calling
The code is sequential labeled as the process running. 
Notes:
1, 6 Alignment is processed with alignment to the hg38 and T2T-CHM13 when keeping MM/ML tags (methylation). 
2. Assembly is processed directly from fastq/ fastq.gz files. single node, 16 cpus. multi node won't help.
3. sniffles used the T2T-CHM13 aligned files to call out structural variants, output is vcf files. paired with sv_circos.py to generate the circos plot. (mod: the circos plot moves the intrachromosomal SVs to the outer boundaries to better visualize, given the rarity)
4. ecDNA calling using CoRAL used hg38 alignment BAM files.
5, 5hmC/5mC calling match with EM-seq cov.gz files.
6, as mentioned earlier in note 1
7, determine the coverage of dark region in the T2T compared to hg38, used 05 python file
8, annotate genes in the dark regions, used 02a, 02b, 06a, 06b, 07a/b/c files to process.
