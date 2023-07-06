# VEME 2023 NGS Variant Calling Tutorial

## Taylor K. Paisie
### `https://taylorpaisie.github.io/VEME_NGS_Variant_Calling/`

### 1. Background and Metadata
#### What is Variant (or SNP) Calling?
#### Variant calling is the process of identifying and cataloging the differences between the observed sequencing reads and a reference genome
#### Variants are usually determined from alignments or BAM files Typical variant calling process:
1. Align reads to the reference genome
2. Correct and refine alignments
3. Determine variants from the alignments
4. Filter the resulting variants for the desired characteristics

#### In the variant calling process you will run into terminologies that are not always used consistently
#### We will attempt to define a few of these terms while making a note that these definitions are not “textbook” definitions designed to accurately capture the proper terminology
### 2. Assessing Read Quality
<figure>
    <img src="variant_calling_steps.png" width="230" height="300">
    <figcaption>Variant Calling Workflow</figcaption>
</figure>
	
1. Downloading SRA files:  
    * Make a directory to download fastq files:  
    `$ mkdir -p data/untrimmed_fastq`  
    `$ cd data/untrimmed_fastq`
    * Download fastq files:  
	`$ curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR197/007/SRR1972917/SRR1972917_1.fastq.gz`   
	`$ curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR197/007/SRR1972917/SRR1972917_2.fastq.gz`  
	`$ curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR197/008/SRR1972918/SRR1972918_1.fastq.gz`  
	`$ curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR197/008/SRR1972918/SRR1972918_2.fastq.gz`  
	`$ curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR197/009/SRR1972919/SRR1972919_1.fastq.gz`  
	`$ curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR197/009/SRR1972919/SRR1972919_2.fastq.gz`  
	`$ curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR197/000/SRR1972920/SRR1972920_1.fastq.gz`  
	`$ curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR197/000/SRR1972920/SRR1972920_2.fastq.gz`  
	
2. Lets take a look at one of our fastq files:
	* In order to view our the fastq file, we must decompress it:  
		`$ gunzip SRR1972917_1.fastq.gz`
	* We can view the first complete read in one of the files our dataset by using head to look at the first four lines:  
		`$ head -n 4 SRR1972917_1.fastq`  
    * Output:  
        `@SRR1972917.1 1/1
        TCCGTGGGGCTGGTACGACAGTATCGATGAGGGTGGACGCTTCAAGGTCAAGCGTATACAGGTCAACCCCAAAGCTAGCCTGAGCCTTCAGAAACACCACC  
        +  
        @CCFDDFFHGHHHGIJJIIJJJJGJJJJGIJIJJFIJJJIIJJJJHHFHHFFFFAADEFEDDDDDDDDDDDD??CCCDDDDDDCCCDDCCCDD:?CABDDB`

    * Each quality score represents the probability that the corresponding nucleotide call is incorrect  
    * This quality score is logarithmically based, so a quality score of 10 reflects a base call accuracy of 90%, but a quality score of 20 reflects a base call accuracy of 99%
    * These probability values are the results from the base calling algorithm and depend on how much signal was captured for the base incorporation  
    `Quality encoding: !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ`  
    `               |         |         |         |         |  `  
    `Quality score:    01........11........21........31........41`

3. Assessing read quality using FastQC
    * FastQC has a number of features which can give you a quick impression of any problems your data may have, so you can take these issues into consideration before moving forward with your analyses   
    * Rather than looking at quality scores for each individual read, FastQC looks at quality collectively across all reads within a sample  
    * The x-axis displays the base position in the read, and the y-axis shows quality scores  
    * For each position, there is a box-and-whisker plot showing the distribution of quality scores for all reads at that position  
    * The horizontal red line indicates the median quality score and the yellow box shows the 1st to 3rd quartile range  
    * This means that 50% of reads have a quality score that falls within the range of the yellow box at that position  
    * The whiskers show the absolute range, which covers the lowest (0th quartile) to highest (4th quartile) values  
    *The plot background is also color-coded to identify good (green), acceptable (yellow), and bad (red) quality scores  


    `$ fastqc -h`  
    `$ fastqc *.fastq*`  

    * Output should look like this:  
    `application/gzip   
    Started analysis of SRR1972917_1.fastq.gz  
    Approx 5% complete for SRR1972917_1.fastq.gz  
    Approx 10% complete for SRR1972917_1.fastq.gz  
    Approx 15% complete for SRR1972917_1.fastq.gz  
    Approx 20% complete for SRR1972917_1.fastq.gz  
    Approx 25% complete for SRR1972917_1.fastq.gz  
    Approx 30% complete for SRR1972917_1.fastq.gz  
    Approx 35% complete for SRR1972917_1.fastq.gz  
    Approx 40% complete for SRR1972917_1.fastq.gz  
    Approx 45% complete for SRR1972917_1.fastq.gz  
    Approx 50% complete for SRR1972917_1.fastq.gz  
    Approx 55% complete for SRR1972917_1.fastq.gz  
    Approx 60% complete for SRR1972917_1.fastq.gz  
    Approx 65% complete for SRR1972917_1.fastq.gz  
    Approx 70% complete for SRR1972917_1.fastq.gz  
    Approx 75% complete for SRR1972917_1.fastq.gz  
    Approx 80% complete for SRR1972917_1.fastq.gz  
    Approx 85% complete for SRR1972917_1.fastq.gz  
    Approx 90% complete for SRR1972917_1.fastq.gz  
    Approx 95% complete for SRR1972917_1.fastq.gz  
    Analysis complete for SRR1972917_1.fastq.gz`  

    * Lets now look at the files created by FastQC:  
    `$ ls` 

    * For each input FASTQ file, FastQC has created a .zip file and a .html file  
    * The .zip file extension indicates that this is actually a compressed set of multiple output files   
    * The .html file is a stable webpage displaying the summary report for each of our samples  
    * We want to keep our data files and our results files separate, so we will move these output files into a new directory within our results/ directory  
    `$ mkdir results/fastqc_untrimmed_reads`  
    `$ mv *.zip results/fastqc_untrimmed_reads`  
    `$ mv *.html results/fastqc_untrimmed_reads`  
    `$ cd results/fastqc_untrimmed_reads`
    * We can now open the .html file to view the FastQC results:  

    <figure>
    <img src="untrimmed_fastqc_picture.png" width="500" height="300">
    <figcaption>FastQC output of SRR1972917_1.fastq.gz</figcaption>
    </figure>



### 3. Trimming and Filtering

1. Trimming the bad quality reads from our fastq files
    * In the previous episode, we took a high-level look at the quality of each of our samples using FastQC  
    * We visualized per-base quality graphs showing the distribution of read quality at each base across all reads in a sample and extracted information about which samples fail which quality checks  
    * It is very common to have some quality metrics fail, and this may or may not be a problem for your downstream application  
    * For our variant calling workflow, we will be removing some of the low quality sequences to reduce our false positive rate due to sequencing error  
    * We will use a program called Trimmomatic to filter poor quality reads and trim poor quality bases from our samples  
    
    `$ trimmomatic`  
    ``Usage:  
       PE [-version] [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-summary <statsSummaryFile>] [-quiet] [-validatePairs] [-basein <inputBase> | <inputFile1> <inputFile2>] [-baseout <outputBase> | <outputFile1P> <outputFile1U> <outputFile2P> <outputFile2U>] <trimmer1>...  
   or:  
       SE [-version] [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-summary <statsSummaryFile>] [-quiet] <inputFile> <outputFile> <trimmer1>...  
   or:  
       -version``  

    * This output shows us that we must first specify whether we have paired end (PE) or single end (SE) reads  
    * Next, we specify what flag we would like to run. For example, you can specify threads to indicate the number of processors on your computer that you want Trimmomatic to use  
    * In most cases using multiple threads (processors) can help to run the trimming faster  
    * These flags are not necessary, but they can give you more control over the command  
    * The flags are followed by positional arguments, meaning the order in which you specify them is important  
    * In paired end mode, Trimmomatic expects the two input files, and then the names of the output files  
    * While, in single end mode, Trimmomatic will expect 1 file as input, after which you can enter the optional settings and lastly the name of the output file  

   * Trimmomatic command options:  
     * ILLUMINACLIP: Cut adapter and other Illumina-specific sequences from the read  
     * SLIDINGWINDOW: Perform a sliding window trimming, cutting once the average
  quality within the window falls below a threshold  
     * LEADING: Cut bases off the start of a read, if below a threshold quality  
     * TRAILING: Cut bases off the end of a read, if below a threshold quality  
     * CROP: Cut the read to a specified length  
     * HEADCROP: Cut the specified number of bases from the start of the read 4
     * MINLEN: Drop the read if it is below a specified length  
     * TOPHRED33: Convert quality scores to Phred-33  
     * TOPHRED64: Convert quality scores to Phred-64  

2. Running Trimmomatic:  
    `$ trimmomatic PE SRR1972917_1.fastq.gz  SRR1972917_2.fastq.gz \`  
    `SRR1972917_1.trim.fastq.gz SRR1972917_1un.trim.fastq.gz \`  
    `SRR1972917_2.trim.fastq.gz SRR1972917_2un.trim.fastq.gz \`  
    `SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15`  

    * Output:  
    ``TrimmomaticPE: Started with arguments:
 SRR1972917_1.fastq.gz SRR1972917_2.fastq.gz SRR1972917_1.trim.fastq.gz SRR1972917_1un.trim.fastq.gz SRR1972917_2.trim.fastq.gz SRR1972917_2un.trim.fastq.gz SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15
Multiple cores found: Using 4 threads
Using PrefixPair: 'AGATGTGTATAAGAGACAG' and 'AGATGTGTATAAGAGACAG'
Using Long Clipping Sequence: 'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG'
Using Long Clipping Sequence: 'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG'
Using Long Clipping Sequence: 'CTGTCTCTTATACACATCTCCGAGCCCACGAGAC'
Using Long Clipping Sequence: 'CTGTCTCTTATACACATCTGACGCTGCCGACGA'
ILLUMINACLIP: Using 1 prefix pairs, 4 forward/reverse sequences, 0 forward only sequences, 0 reverse only sequences
Quality encoding detected as phred33
Input Read Pairs: 4377867 Both Surviving: 1241328 (28.35%) Forward Only Surviving: 2133670 (48.74%) Reverse Only Surviving: 18680 (0.43%) Dropped: 984189 (22.48%)
TrimmomaticPE: Completed successfully``  

    * List out files created by Trimmomatic:  
    `$ ls SRR1972917*`
    * Trimmed files should be smaller in size than our untrimmed fastq files
  
3. Running a for loop on all fastq files  

    ``$ for infile in *_1.fastq.gz
        do
            base=$(basename ${infile} _1.fastq.gz)
            trimmomatic PE ${infile} ${base}_2.fastq.gz \  
            ${base}_1.trim.fastq.gz ${base}_1un.trim.fastq.gz \  
            ${base}_2.trim.fastq.gz ${base}_2un.trim.fastq.gz \  
            SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15   
        done``  

4. Moving trimmed fastq files to a new directory:
    * We have now completed the trimming and filtering steps of our quality control process! Before we move on, let’s move our trimmed FASTQ files to a new subdirectory within our data/ directory  
  
    `$ cd data/untrimmed_fastq`  
    `$ mkdir -p ../trimmed_fastq`  
    `$ mv *trim* ../trimmed_fastq`  
    `$ cd ../trimmed_fastq`  
    `$ ls -al`  






### 4. Reference Based Mapping  




### 5. Variant Calling
### 6. Visualizing the Results
### 7. Optional (if we have time): Automating Variant Calling Workflow