# Intro to genomic variant calling with Genome Analysis Toolkit (GATK) #

This QuickBytes outlines how to run a pipeline based on Genome Analysis Toolkit v4.xx (GATK4) best practices, a common pipeline for processing genomic data from Illumina platforms. Major modifications from “true” best practices are done to facilitate using this pipeline for both model and non-model organisms. This pipeline is designed for genome-wide data but is similar to what one would use for target capture data. Here we to summarize the steps and get you set up, but much more extensive documentation (including other Best Practices pipelines) can be found [here](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows). Specifically, the Best Practices informing this pipeline are the [data pre-processing workflow](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery) and the [germline short variant discovery workflow](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-). GATK also has workflows developed already, but those are meant to be run as a “black box”. Here we aim to give you sample commands to emulate these scripts, which will also allow you easily modify the pipeline.

The goal of this is to output Single Nucleotide Polymorphisms (SNPs) and optionally indels for a given dataset. This same pipeline can be used for humans, model organisms, and non-model organisms. Spots that can leverage information from model organisms are noted, but those steps can also be bypassed for the sake of generality. Because sample size and depth of coverage are often lower in non-model organisms, filtering recommendations and memory requirements will vary. Note that this assumes you are using paired-end data and will differ slightly if you use unpaired.

The basic steps are aligning and processing raw reads into binary alignment map (BAM) files, optionally getting descriptive metrics about the samples’ sequencing and alignment, calling variants to produce genomic variant call format (GVCF) files, and finally filtering those variants for analysis.
Please note that you must cite any program you use in a paper. At the end of this, we have provided citations you would include for the programs we ran here.

## Preliminary stuff ##

### Module and directories ###

We will be using conda to make an environment to load within our PBS script. First, if you haven’t already, set up conda as follows:

	module load miniconda3-4.7.12.1-gcc-4.8.5-lmtvtik
	# can also use other conda modules
	conda init bash
	<exit and re-log or restart shell>

This following line will create an environment and install the most recent versions of what we need. We assume you run this before starting up your job.

	conda create -n gatk-env -c bioconda -c conda-forge gatk4 bwa samtools picard

Alternatively, you can load these as modules if you are on Wheeler (Xena only has Samtools now), but they may not be the most recent versions:

	module load bwa-0.7.17-intel-18.0.2-7jvpfu2
	module load samtools-1.9-gcc-7.3.0-tnzvvzt
	module load picard-2.20.8-gcc-4.8.5-3yh2dzv
	module load gatk-4.1.4.1-gcc-4.8.5-python3-fqiavji

If you are parallelizing (see “Calling variants with HaplotypeCaller”), you need this:

	module load parallel-20170322-gcc-4.8.5-2ycpx7e
	source $(which env_parallel.bash)

The directories we will need (other than the home directory) are a raw_reads directory for the demultiplexed reads and the following for various intermediate files to go into. Alternatively, if you don’t want to move around all your reads, just replace the path in the BWA call with that path.

	mkdir alignment
	# next three are only if you get optional metrics
	mkdir alignments/alignment_summary
	mkdir alignments/insert_metrics
	mkdir alignments/depth
	mkdir alignments/dedup_temp
	mkdir bams
	mkdir raw_gvcfs
	mkdir combined_vcfs
	mkdir analysis_vcfs

We will also be using a few variables throughout this that we can set now. Most notable is the number of threads and the working directory (that we call “src”).

	src=$PBS_O_WORKDIR
	reference=$src/reference
	threads=[# of threads to use]

### Sample Names ###

To keep our script short, and outputs easy to understand, we will use consistent sample names for each step, and keep the sample names in a file. We assume this file is named “sample_list”. The file should have one sample name per line, as below:

	Sample1
	Sample2
	Sample3
	
A sample PBS script at the end of the document will include an example loop for each section. The way we implement it is:

	while read sample; do
		<tasks>
	done < $src/sample_list

### Demultiplexing ###

Because it is not covered by best practices, and is often done by the sequencing center, we will not go into the details of demultiplexing here. We recommend you use Illumina’s software [bcl2fastq](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html) if you have the data in .bcl format, and [saber](https://github.com/najoshi/sabre) if it has already been converted to fastq format and it does not have duel combinatorial barcodes.

## The Pipeline ##

### Alignment and Pre-processing ###

This section prepares BAM files for variant calling. First, we need to index our reference. We'll do this two ways, one for bwa and one for GATK:

	bwa index -p $reference ${reference}.fa
	samtools faidx ${reference}.fa -o ${reference}.fa.fai

Then, we need to align demultiplexed reads to a reference. For this step, we will use the Burrough-Wheeler Aligner’s (BWA) mem algorithm. Another common option is [Bowtie](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml). One important flag here is the -R flag, which is the read group and sample ID for a given sample. We assume that these samples are in the same read group. The command looks like this:

	bwa mem \
		-t $threads -M \
		-R "@RG\tID:${sample}\tPL:ILLUMINA\tLB:${sample}\tSM:${sample}" \
		$reference \
		$src/raw_reads/${sample}_R1.fastq.gz \
		$src/raw_reads/${sample}_R2.fastq.gz \
		> $src/alignment/${sample}.sam

The next step is to mark PCR duplicates to remove bias, sort the file, and convert it to the smaller BAM format for downstream use. GATK’s new MarkDuplicatesSpark performs all these tasks, but needs a temporary directory to store intermediate files:

	gatk MarkDuplicatesSpark \
		-I $src/alignment/”$sample”.sam \
		-M $src/bams/”$sample”_dedup_metrics.txt \
		--tmp-dir $src/alignments/dedup_temp \
		-O $src/bams/”$sample”_dedup.bam

The final step is to recalibrate base call scores. This applies machine learning to find where quality scores are over or underestimated based on things like read group and cycle number of a given base. This is strongly recommended, but is rarely possible for non-model organisms, as a file of known polymorphisms is needed. Note, however, that it can take a strongly filtered VCF from the end of the pipeline, before running the entire pipeline again (but others haven’t found much success with this method).

	gatk BaseRecalibrator \
		-I $src/bams/${sample}_dedup.bam \
		-R ${reference}.fa \
		--known-sites $src/[name-of-known-sites].vcf \
		-O $src/bams/${sample}_recal_data.table
		
	gatk ApplyBQSR \
		-I $src/bams/${sample}_dedup.bam \
		-R ${reference}.fa \
		--bqsr-recal-file $src/bams/${sample}_recal_data.table \
		-O $src/bams/${sample}_recal.bam

We recommend combining these steps per sample for efficiency and smoother troubleshooting. One issue is that we do not want large SAM files piling up. This can either be done by piping BWA output directly to MarkDuplicatesSpark or removing the SAM file after each loop. In case you want to save the SAM files, we did the latter (this isn’t a bad idea if you have the space, in case there is a problem with generating BAM files). If you are doing base recalibration, you can also add “rm ${sample}\_debup.bam” to get rid of needless BAM files. Later in the pipeline, we assume you did base recalibration, so will use the [sample]\_recal.bam file. If you did not use base recalibration, use [sample]\_dedup.bam file in its place.

### Collect alignment and summary statistics (optional)

This step is optional, and is not part of GATK’s best practices, but we recommend it. It will output important stats for assessing sample quality. Picard’s “CollectAlignmentSummaryMetrics” gives several helpful statistics about the alignment for a given sample. Picard’s “CollectInsertSizeMetrics” gives information about the distribution of insert sizes. Samtools’s “depth” gives information about the read depth of the sample.

	picard CollectAlignmentSummaryMetrics \
		R=${reference}.fa \
		I=$src/bams/${sample}_recal.bam \
		O=$src/alignments/alignment_summary/${sample}_alignment_summary.txt
		
	picard CollectInsertSizeMetrics \
		INPUT=$src/bams/${sample}_recal.bam \
		OUTPUT=$src/alignments/insert_metrics/${sample}_insert_size.txt \
		HISTOGRAM_FILE=$src/${sample}_insert_hist.pdf
		
	samtools depth \
		-a $src/bams/${sample}_recal.bam \
		> $src/alignments/depth/${sample}_depth.txt

### Calling variants with HaplotypeCaller

The first step is to go through each BAM files and call SNPs on it using HaplotypeCaller in GVCF mode (“-ERC GVCF” flag), resulting in GVCFs as output.

	gatk HaplotypeCaller \
		-R ${reference}.fa \
		-I $src/bams/${sample}_recal.bam \
		-O $src/gvcfs/${sample}_raw.g.vcf.gz \
		-ERC GVCF

One issue with HaplotypeCaller is that it takes a long time but is not programmed to be parallelized by default. However, we can use GNU parallel to easily solve that problem! One note is that we need env_parallel if we are using conda. We can do this two ways. If you have small inputs, it should be fine to parallelize by sample like this. Note that we restrict the memory such that each job can only max out the core it's on:

	cat sample_list | env_parallel --sshloginfline $PBS_NODEFILE \
		“gatk --java-options "-Xmx6g" HaplotypeCaller \
		-R ${reference}.fa \
		-I $src/bams/{}_recal.bam \
		-O $src/gvcfs/${}_raw.g.vcf.gz \
		-ERC GVCF”

However, if you are dealing with large files, HaplotypeCaller may take longer than your walltime! We'll do this by breaking the job into multiple intervals, which will be every contig in our reference. We'll then combine every GVCF for a given sample.:

	cut -f 1 ${reference}.faa.fai > $src/chromosomes.list
	
	while read sample; do
	      cat $src/chromosomes.list | env_parallel --sshloginfile $PBS_NODEFILE \
		  "gatk --java-options "-Xmx6g" HaplotypeCaller \
		  -R ${reference}.fa \
		  -I $src/bams/${sample}_recal.bam \
		  -O $src/gvcfs/${sample}_{}_raw.g.vcf.gz \
		  -ERC GVCF"
	done < $src/sample_list

### Consolidating and Genotyping ###

This next step has two options, GenomicsDBImport and CombineGVCFs. GATK recommends GenomicsDBImport, as it is more efficient for large datasets, but it performs poorly on references with many contigs. CombineGVCFs takes a lot of memory for datasets with many samples. Overall, it seems that GenomicsDBImport is more suited for model organisms, while CombineGVCFs is more convenient for non-model organisms (which often have smaller projects and less contiguous references). GenomicsDBImport must have intervals (generally corresponding to contigs or chromosomes) specified. GenomicsDBImport can take a file specifying GVCFs, but because CombineGVCFs cannot we just make a list of samples to combine programmatically and plug it in. Here is how we generate that command:

	gvcf_names=""
	while read sample; do
		gvcf_names="${src}/${gvcf_names}-V ${src}/gvcfs/${sample}_raw.g.vcf.gz "
	done < $src/sample_list

If you do use GenomicsDBImport, or want to genotype contigs/chromosomes independently, we'll need intervals for it to work with. This is the same used for parallelizing above:

	cut -f 1 ${reference}.faa.fai > $src/chromosomes.list	

An example GenomicsDBImport command is:

	gatk GenomicsDBImport \
		${gvcf_names} \
		--genomicsdb-workspace-path $src/genomic_database \
		--tmp-dir $src/temp
		-L chromosomes.list

And an example for CombineGVCFs is:

	gatk CombineGVCFs \
		-R ${reference}.fa \
		${gvcf_names} \
		-O $src/combined_vcfs/combined_gvcf.g.vcf.gz

The next step is to genotype the combined (cohort) GVCF file. Here’s a sample command for GenomicsDBImport:

	gatk GenotypeGVCFs \
		-R ${reference}.fa \
		-V gendb://genomic_database \
		-O $src/combined_vcfs/combined_vcf.vcf.gz
	
And one for CombineGVCFs:

	gatk GenotypeGVCFs \
		-R ${reference}.fa \
		-V $src/combined_vcfs/combined_gvcf.g.vcf.gz \
		-O $src/combined_vcfs/combined_vcf.vcf.gz

### Selecting and filtering variants

This first step is optional, but here we separate out indels and SNPs. Note that we don’t use indels down the line, but similar filters can be applied.

	gatk SelectVariants \
		-R ${reference}.fa \
		-V $src/combined_vcfs/combined_vcf.vcf.gz \
		-selectType SNP \
		-o $src/combined_vcfs/raw_snps.vcf.gz

	gatk SelectVariants \
		-R ${reference}.fa \
		-V $src/combined_vcfs/combined_vcf.vcf.gz \
		-selectType INDEL \
		-o $src/combined_vcfs/raw_indel.vcf.gz

Here are some good sample filters. The “DP_filter” is depth of coverage (you will probably want to change this), “Q_filter” is quality score, “QD_filter” is quality by depth (avoids artificial inflation of calls), “MQ_filter” is a mapping quality filter, and “FS_filter” is a strand bias filter (higher value means higher bias). Note that DP is better for low depth samples, while QD is better for high depth. More info can be found on [GATK’s website](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants).

	gatk VariantFiltration \
		-R ${reference}.fa \
		-V $src/combined_vcfs/raw_snps.vcf.gz \
		-O $src/analysis_vcfs/filtered_snps.vcf  \
		-filter-name “DP_filter” -filter “DP < 4” \
		-filter-name “Q_filter” -filter “QUAL < 30.0” \
		-filter-name “QD_filter” -filter “QD < 2.0” \
		-filter-name “MQ_filter” -filter “MQ < 40.0” \
		-filter-name “FS_filter” -filter “FS > 60.0”

This will give us our final VCF! Note that the filtered SNPs are still included, just with a filter tag. You can use something like [VCFtools’](http://vcftools.sourceforge.net/) “--remove-filtered-all” flag to get rid of them. It is often advisable to make this and other adjustments using VCFtools before moving on to other analyses.

## Sample PBS Script ##

Here is a sample PBS script combining everything we have above. One reason to break up steps like we did is for improved checkpointing (without having to write code checking if files are already present). Once you are finished running a block of code, you can just comment it out. Similarly, if you can only get part way through your sample list, you can copy is and remove samples that have already run.

	module load miniconda3-4.7.12.1-gcc-4.8.5-lmtvtik
	source activate gatk-env
	
	src=$PBS_O_WORKDIR
	reference=${src}/reference.fa
	threads=[# of threads to use]
	
	# Section for alignment, marking duplicates, and base recalibration.

	while read sample; do
		bwa mem \
			-t $threads -M \
			-R ‘@RG\tID:${sample}\tPL:ILLUMINA\tLB:”$sample”\tSML”$sample” \
			$reference \
			$src/raw_reads/${sample}_R1.fastq.gz \
			$src/raw_reads/${sample}_R2.fastq.gz \
			> $src/sams/${sample}.sam
		gatk MarkDuplicatesSpark \
			-I $src/sams/${sample}.sam \
			-M $src/bams/${sample}_dedup_metrics.txt \
			--tmp-dir $src/alignments/dedup_temp \
			-O $src/bams/${sample}_dedup.bam
		gatk BaseRecalibrator \
			-I $src/bams/${sample}_dedup.bam \
			-R $reference \
			--known-sites $src/[name-of-known-sites].vcf \
			-O $src/bams/${sample}_recal_data.table
		gatk ApplyBQSR \
			-I $src/bams/${sample}_dedup.bam \
			-R $reference \
			--bqsr-recal-file bams/${sample}_recal_data.table \
			-O $src/bams/${sample}_recal.bam
		rm $src/sams/${sample}.sam
	done < $src/sample_list

	# Loop for collecting metrics.
	# Remember to change from _recal to _dedup if you can’t do base recalibration.

	while read sample; do
		picard CollectAlignmentSummaryMetrics \
			R=${reference}.fa \
			I=$src/bams/${sample}_recal.bam \
			O=$src/alignments/alignment_summary/${sample}_alignment_summary.txt
		picard CollectInsertSizeMetrics \
			INPUT=$src/bams/${sample}_recal.bam \
			OUTPUT=$src/alignments/insert_metrics/${sample}_insert_size.txt \
			HISTOGRAM_FILE=$src/${sample}_insert_hist.pdf
		samtools depth \
			-a $src/bams/${sample}_recal.bam \
			> $src/alignments/depth/${sample}_depth.txt
	done < $src/sample_list

	# Loop for HaploType Caller, probably the most likely to need checkpoints.
	# This will look very different if you parallelize! No loop will be needed.
	# Go to the HaplotypeCaller section for more info. 
	
	while read sample; do
		gatk HaplotypeCaller \
			-R ${reference}.fa \
			-I $src/bams/${sample}_sorted_bqsr.bam \
			-O $src/gvcfs/${sample}_raw.g.vcf.gz \
			-ERC GVCF
	done < $src/sample_list

	# The rest, assuming you do CombineGVCFs

	gvcf_names=""
	while read sample; do
		gvcf_names = "${gvcf_names}-V ${src}/gvcfs/${sample}_raw.g.vcf.gz "
	done < $src/sample_list

	gatk CombineGVCFs \
		-R ${reference}.fa \
		${gvcf_names} \
		-O $src/combined_vcfs/combined_gvcf.g.vcf.gz

	gatk GenotypeGVCFs \
		-R ${reference}.fa \
		-V $src/combined_vcfs/combined_gvcf.g.vcf.gz \
		-O $src/combined_vcfs/combined_vcf.vcf.gz

	gatk SelectVariants \
		-R ${reference}.fa \
		-V $src/combined_vcfs/combined_vcf.vcf.gz \
		-selectType SNP \
		-o $src/combined_vcfs/raw_snps.vcf.gz

	gatk SelectVariants \
		-R ${reference}.fa \
		-V $src/combined_vcfs/combined_vcf.vcf.gz \
		-selectType INDEL \
		-o $src/combined_vcfs/raw_indel.vcf.gz

	gatk VariantFiltration \
		-R ${reference}.fa \
		-V $src/combined_vcfs/raw_snps.vcf.gz \
		-O $src/analysis_vcfs/filtered_snps.vcf  \
		-filter-name “DP_filter” -filter “DP < 4” \
		-filter-name “Q_filter” -filter “QUAL < 30.0” \
		-filter-name “QD_filter” -filter “QD < 2.0” \
		-filter-name “MQ_filter” -filter “MQ < 40.0” \
		-filter-name “FS_filter” -filter “FS > 60.0”

## Troubleshooting ##

If you need help troubleshooting an error, make sure to let us know the size of your dataset (number of individuals and approximate number of reads should suffice, unless coverage varies between individuals), GATK version, node details, and any error messages output.

## Citations ##

Li, H., & Durbin, R. (2009). Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics, 25(14), 1754–1760. https://doi.org/10.1093/bioinformatics/btp324

Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., … Durbin, R. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16), 2078–2079. https://doi.org/10.1093/bioinformatics/btp352

McKenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., … DePristo, M. A. (2010). The genome analysis toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data. Genome Research, 20(9), 1297–1303. https://doi.org/10.1101/gr.107524.110

Picard toolkit. (2019). Broad Institute, GitHub Repository. https://doi.org/http://broadinstitute.github.io/picard/
