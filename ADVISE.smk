#!/usr/bin/env Snakemake

# Used to analyse all ADVISE samples and cohorts in a consistent manner
# This includes both samples run on the MiSeq and HiSeq

# Environment configurations
workDir = "/mnt/thanos_lv/crushton/ADVISE/analysis_final/"
refGenome = "/morinlab/reference/igenomes/Homo_sapiens/GSC/GRCh38/Sequence/WholeGenomeFasta/GRCh38_no_alt.fa"
sampleDir = "/mnt/thanos_lv/crushton/ADVISE/analysis_final/01-BAMs/"
logDir = workDir + "logs/"
capSpace = "/mnt/thanos_lv/crushton/ADVISE/analysis_final/00-helperfiles/Targets.Validated.sort.bed6.gz"
capSpaceGATK = "/mnt/thanos_lv/crushton/ADVISE/analysis_final/00-helperfiles/Targets.Validated.sort.bed"
capSpacePicard = "/mnt/thanos_lv/crushton/ADVISE/analysis_final/00-helperfiles/Targets.Validated.list"
sampleFile = workDir + "00-helperfiles/samplelist.hiseq.txt"  # Contains two columns: sample name and path
ffpolish = "/home/crushton/Software/FFPolish/src/cli.py"

# Input files for running SAGE
sagePath = "/mnt/thanos_lv/crushton/ADVISE/analysis_final/00-helperfiles/sage-2.7.jar"
hotspotVCF = "/mnt/thanos_lv/crushton/ADVISE/analysis_final/00-helperfiles/KnownHotspots.hg38.vcf.gz"
highConfBED = "/mnt/thanos_lv/crushton/ADVISE/analysis_final/00-helperfiles/Targets.Validated.merge.bed3"

# Config files
mantaConfig = workDir + "00-helperfiles/configManta.py.ini"
strelkaConfig = workDir + "00-helperfiles/configureStrelkaSomaticWorkflow.py.ini"
dbSNPVCF = "/morinlab/reference/igenomes/Homo_sapiens/GSC/GRCh38/Annotation/Variation/dbsnp.all.151.chr.vcf.gz"
strelka1Script = "/software/strelka-workflow/1.0.14/bin/configureStrelkaWorkflow.pl"
strelka1Config = workDir + "00-helperfiles/ADVISE.strelka_config.ini"
augmentMAFStrelka = "/mnt/thanos_lv/crushton/ADVISE/Analysis/SNVs/Strelka_Aggregated/00-HelperFiles/augment_maf_strelka.py"
augmentMAF = "/mnt/thanos_lv/crushton/Software/lab_scripts/augment_maf/augment_maf.py"
postFiltScript = "/mnt/thanos_lv/crushton/ADVISE/analysis_final/00-helperfiles/filter_maf.py"

# VEP data paths
vepCache = "/home/crushton/.vep/"
ncbiBuild = "GRCh38"

# Load sample information
allSamples = []
tSamples = []
nSamples = {}
tBiopsySamples = {}

# Output folders for each step of the analysis
qualimapOut = workDir + "01-qualimap/"
hsMetOut = workDir + "02-hsmet/"
readGroupOut = workDir + "03-withrg/"
clipOut = workDir + "04-clipo/"
haplotypeCallOut = workDir + "05-haplotypecaller/"
mantaOut = workDir + "06-manta/"
strelka2Out = workDir + "07-strelka2/"
strelka2Merged = workDir + "08-strelka2VCFs/"
strelka2MAFOut = workDir + "09-strelka2MAFs/"
strelka1Out = workDir + "10-strelka1/"
strelka1Merged = workDir + "11-strelka1VCFs/"
strelka1MAFOut = workDir + "12-strelka1MAFs/"
strelka1AugOut = workDir + "13-strelka1AugMAF/"
mutect2Out = workDir + "14-mutect2/"
mutect2Filt =  workDir + "15-mutect2VCFs/"
mutect2MAFOut = workDir + "16-mutect2MAFs/"
haplotypeMAFOut = workDir + "17-haplotypeMAFs/"
ffpolishOut = workDir + "19-strelka2FFPolish/"
ffpolishMAFOut = workDir + "20-strelka2FFPolishMAF/"
sageOut = workDir + "22-sage/"
sagePassedOut = workDir + "23-sagePassed/"
sageMAFOut = workDir + "24-sageMAFs/"
sageFiltOut = workDir + "25-sagePostFilt/"
biopsyAugOut = workDir + "26-tumorsens/"

with open(sampleFile) as f:
	for line in f:
		line = line.rstrip("\n").rstrip("\r")
		allSamples.append(line)  # Includes both tumour and normal samples
		if "-N" in line:
			continue
		tSamples.append(line)
		sID, runID = line.split("-")
		runID = runID.split(".")[1]
		nID = sID + "-N." + runID
		nSamples[line] = nID

		# Store the associated tumour tissue biopsy sample
		# This only applies for patients which have a tumour
		if "EN" in sID or "OV" in sID:
			tBiopsyID = sID + "-T." + runID
			tBiopsySamples[line] = tBiopsyID



### START PIPELINES ###

rule all:
	input:
		expand(qualimapOut + "{sampleName}.qualimap.pdf", sampleName=allSamples),
		expand(hsMetOut + "{sampleName}.hsmet.tsv", sampleName=allSamples),
		expand(haplotypeMAFOut + "{sampleName}.hapcaller.maf", sampleName=allSamples),
		expand(strelka2MAFOut + "{tSample}.strelka2.maf", tSample=tSamples),
		expand(strelka1AugOut + "{tSample}.strelka1.augment.maf", tSample=tSamples),
		expand(mutect2MAFOut + "{tSample}.mutect2.maf", tSample=tSamples),
#		expand(ffpolishMAFOut + "{tSample}.strelka2.ffpolish.maf", tSample=tSamples),
		expand(sageFiltOut + "{tSample}.sage.filt.maf", tSample=tSamples),
		expand(biopsyAugOut + "{tSample}.biopsyVar.maf", tSample=tBiopsySamples.keys())

### STEP 1: QC ###
rule Qualimap:
	input:
		inBAM = sampleDir + "{allSamples}.bam"
	output:
		qualiPDF = "{qualimapOut}/{allSamples}.qualimap.pdf"
	threads:
		4
	log:
		logFile = logDir + "{allSamples}.qualimap.log"
	message:
		"Running Qualimap on {wildcards.allSamples}"
	conda:
		"envs/gatk4.yaml"
	shell:
		"qualimap bamqc -bam {input.inBAM} -outdir {qualimapOut} -outfile {wildcards.allSamples}.qualimap.pdf -sd -gff {capSpaceGATK} -nt {threads} 2> {log.logFile} > {log.logFile}"

rule CollectHsMetrics:
	input:
		inBAM = sampleDir + "{allSamples}.bam"
	output:
		hsMet = hsMetOut + "{allSamples}.hsmet.tsv",
		tarCov = hsMetOut + "{allSamples}.tarcov.tsv"
	log:
		logFile = logDir + "{allSamples}.CollectHsMetrics.log"
	conda:
		"envs/gatk4.yaml"
	message:
		"Running Picard CollectHsMetrics on {wildcards.allSamples}"
	shell:
		"picard CollectHsMetrics I={input.inBAM} O={output.hsMet} COVMAX=10000 PER_TARGET_COVERAGE={output.tarCov} BI={capSpacePicard} N=ADVISE TI={capSpacePicard} R={refGenome} 2> {log.logFile}"
		
### STEP 2: CALL GERMLINE VARIANTS ###
rule AddReadGroups:
	input:
		inBAM = sampleDir + "{allSamples}.bam"
	output:
		rgBAM = readGroupOut + "{allSamples}.withRG.bam"
	message:
		"Adding read groups to {wildcards.allSamples}"
	threads:
		3
	conda:
		"envs/samtools.yaml"
	shell:
		"samtools addreplacerg -@ {threads} -r \"ID:{wildcards.allSamples}\tPL:ILLUMINA\tCN:BCCRC\tPM:Nexterra\tSM:{wildcards.allSamples}\" {input.inBAM} | samtools view -b > {output.rgBAM} && samtools index {output.rgBAM}"

rule HaplotypeCaller:
	input:
		rgBAM = readGroupOut + "{allSamples}.withRG.bam"
	output:
		hapVCF = haplotypeCallOut + "{allSamples}.hapcaller.vcf"
	message:
		"Running HaplotypeCaller on {wildcards.allSamples}"
	log:
		logFile = logDir + "{allSamples}.haplotypecaller.log"
	conda:
		"envs/gatk4.yaml"
	shell:
		"gatk HaplotypeCaller -R {refGenome} -I {input.rgBAM} -O {output.hapVCF} --dbsnp {dbSNPVCF} -L {capSpaceGATK} 2> {log.logFile}"

### STEP 3: CALL SOMATIC VARIANTS
# Lets use Strelka, Strelka2, and MuTect2
rule MantaConfig:
	input:
		tBAM = sampleDir + "{tSamples}.bam",
		nBAM = lambda wildcards: sampleDir + nSamples[wildcards.tSamples] + ".bam"
	output:
		mantaWorkflow = mantaOut + "{tSamples}/runWorkflow.py"
	message:
		"Configurating Manta for {wildcards.tSamples}"
	log:
		logFile = logDir + "{tSamples}.configmanta.log"
	conda:
		"envs/strelka2manta.yaml"
	shell:
		"configManta.py --tumorBam {input.tBAM} --bam {input.nBAM} --exome --runDir {mantaOut}/{wildcards.tSamples} --callRegions {capSpace} --config {mantaConfig} 2> {log.logFile} > {log.logFile}"

rule MantaRun:
	input:
		mantaWorkflow = mantaOut + "{tSamples}/runWorkflow.py"
	output:
		candSmallIndels = mantaOut + "{tSamples}/results/variants/candidateSmallIndels.vcf.gz"
	log:
		logFile = logDir + "{tSamples}.mantarun.log"
	threads: 1
	message:
		"Running Manta on {wildcards.tSamples}"
	conda:
		"envs/strelka2manta.yaml"
	shell:
		"{input.mantaWorkflow} -m local -j {threads} 2> {log.logFile}"

rule Strelka2Config:
	input:
		tBAM = sampleDir + "{tSamples}.bam",
		nBAM = lambda wildcards: sampleDir + nSamples[wildcards.tSamples] + ".bam",
		candSmallIndels = mantaOut + "{tSamples}/results/variants/candidateSmallIndels.vcf.gz"
	output:
		strelka2Workflow = strelka2Out + "{tSamples}/runWorkflow.py"
	message:
		"Configure Strelka2 for {wildcards.tSamples}"
	log:
		logFile = logDir + "{tSamples}.configstrelka2.log"
	conda:
		"envs/strelka2manta.yaml"
	shell:
		"rmdir {strelka2Out}{wildcards.tSamples}/ ; configureStrelkaSomaticWorkflow.py --tumorBam {input.tBAM} --normalBam {input.nBAM} --outputCallableRegions --referenceFasta {refGenome} --indelCandidates {input.candSmallIndels} "
		"--exome --callRegions {capSpace} --runDir {strelka2Out}{wildcards.tSamples}/ --config {strelkaConfig} 2> {log.logFile} > {log.logFile}"

rule Strelka2Run:
	input:
		strelka2Workflow = strelka2Out + "{tSamples}/runWorkflow.py"
	output:
		strelka2SNVs = strelka2Out + "{tSamples}/results/variants/somatic.snvs.vcf.gz",
		strelka2indels = strelka2Out + "{tSamples}/results/variants/somatic.indels.vcf.gz"
	threads: 4
	message:
		"Running Strelka2 on {wildcards.tSamples}"
	log:
		logFile = logDir + "{tSamples}.strelka2run.log"
	conda:
		"envs/strelka2manta.yaml"
	shell:
		"{input.strelka2Workflow} -m local -j {threads} 2> {log.logFile}"

# Run Strelka1
rule ClipOverlap:
	input:
		inBAMs = sampleDir + "{allSamples}.bam"
	output:
		clipBAM = clipOut + "{allSamples}.clipO.bam"
	message:
		"Running ClipOverlap on {wildcards.allSamples}"
	log:
		logFile = logDir + "{allSamples}.clipo.log"
	run:
		if ".M" in input.inBAMs:
			# Symlink MiSeq BAMs since they have already been clipped
			shell("ln -s {input.inBAMs} {output.clipBAM}")
			shell("ln -s {input.inBAMs}.bai {output.clipBAM}.bai")
		else:
			shell("dellingr clip -i {input.inBAMs} -o - 2> {log.logFile} | samtools sort | samtools calmd -b - {refGenome} 2>> {log.logFile} > {output.clipBAM} && samtools index {output.clipBAM}")

rule Strelka1Config:
	input:
		tBAM = clipOut + "{tSamples}.clipO.bam",
		nBAM = lambda wildcards: clipOut + nSamples[wildcards.tSamples] + ".clipO.bam",
	output:
		strelka1Conf = strelka1Out + "{tSamples}/Makefile"
	message:
		"Configuring Strelka1 for {wildcards.tSamples}"
	log:
		logFile = logDir + "{tSamples}.configstrelka1.log"
	shell:
		"rmdir {strelka1Out}/{wildcards.tSamples}; {strelka1Script} --tumor {input.tBAM} --normal {input.nBAM} --ref {refGenome} --config {strelka1Config} --output-dir {strelka1Out}/{wildcards.tSamples} > {log.logFile} 2> {log.logFile}"

rule Strelka1Run:
	input:
		strelka1Conf = strelka1Out + "{tSamples}/Makefile"
	output:
		strelka1SNVs = strelka1Out + "{tSamples}/results/passed.somatic.snvs.vcf",
		strelka1Indels = strelka1Out + "{tSamples}/results/passed.somatic.indels.vcf"
	message:
		"Running Strelka1 on {wildcards.tSamples}"
	log:
		logFile = logDir + "{tSamples}.strelka1run.log"
	threads:
		4
	shell:
		"make -C {strelka1Out}/{wildcards.tSamples} -j {threads} 2> {log.logFile} > {log.logFile}"

# Run SAGE
rule SageRun:
	input:
		tBAM = readGroupOut + "{tSamples}.withRG.bam",
		nBAM = lambda wildcards: readGroupOut + nSamples[wildcards.tSamples] + ".withRG.bam"
	output:
		sageVcf = sageOut + "{tSamples}.sage.vcf"
	params:
		normName = lambda wildcards: nSamples[wildcards.tSamples]
	message:
		"Running SAGE on {wildcards.tSamples}"
	threads:
		1
	log:
		logFile = logDir + "{tSamples}.sage.log"
	conda:
		"envs/sage.yaml"
	shell:
		"java -jar {sagePath} -tumor_bam {input.tBAM} -tumor {wildcards.tSamples} -out {output.sageVcf} -ref_genome {refGenome} -hotspots {hotspotVCF} -panel_bed {highConfBED} -high_confidence_bed {highConfBED} -assembly hg38 -threads {threads} -reference_bam {input.nBAM} -validation_stringency SILENT -reference {params.normName} -max_read_depth 100000000 2> {log.logFile}"
		
rule SagePassed:
	input:
		sageVcf = sageOut + "{tSamples}.sage.vcf"
	output:
		sagePassed = sagePassedOut +  "{tSamples}.sage.passed.vcf"
	message:
		"Filtering SAGE calls for {wildcards.tSamples}"
	conda:
		"envs/samtools.yaml"
	shell:
		"egrep '^#|PASS' {input.sageVcf} | bedtools intersect -header -a - -b {highConfBED} | uniq > {output.sagePassed}"


# Run Mutect2
rule Mutect2Run:
	input:
		tBAM = readGroupOut + "{tSamples}.withRG.bam",
		nBAM = lambda wildcards: readGroupOut + nSamples[wildcards.tSamples] + ".withRG.bam"
	output:
		mutect2AllMut = mutect2Out + "{tSamples}.allVar.vcf"
	params:
		normRG = lambda wildcards: nSamples[wildcards.tSamples]
	message:
		"Running Mutect2 on {wildcards.tSamples}"
	conda:
		"envs/gatk4.yaml"
	log:
		logFile = logDir + "{tSamples}.mutect2.log"
	shell:
		"gatk Mutect2 -R {refGenome} -I {input.tBAM} -I {input.nBAM} -normal {params.normRG} -O {output.mutect2AllMut} -L {capSpaceGATK} 2> {log.logFile}"

rule FilterMutect2:
	input:
		mutect2AllMut = mutect2Out + "{tSamples}.allVar.vcf"
	output:
		mutect2Filt = mutect2Filt + "{tSamples}.mutect2.filt.vcf",
		mutect2VCF = mutect2Filt + "{tSamples}.mutect2.passed.vcf"
	message:
		"Filtering Mutect2 calls for {wildcards.tSamples}"
	log:
		logFile = logDir + "{tSamples}.mutect2filt.log"
	conda:
		"envs/gatk4.yaml"
	shell:
		"gatk FilterMutectCalls -V {input.mutect2AllMut} -O {output.mutect2Filt} -R {refGenome} 2> {log.logFile} &&"
		"egrep '^#|PASS' {output.mutect2Filt} > {output.mutect2VCF}"

### STEP 4: AGGREGATE AND ANNOTATE VARIANT CALLS ###
rule MergeStrelka2VCFs:
	input:
		strelka2SNVs = strelka2Out + "{tSamples}/results/variants/somatic.snvs.vcf.gz",
		strelka2indels = strelka2Out + "{tSamples}/results/variants/somatic.indels.vcf.gz"
	output:
		strelka2VCF = strelka2Merged + "{tSamples}.strelka2.passed.vcf"
	message:
		"Merging Strelka2 calls for {wildcards.tSamples}"
	log:
		logFile = logDir + "{tSamples}.mergestrelka2.log"
	conda:
		"envs/samtools.yaml"
	shell:
		"bcftools concat -a {input.strelka2SNVs} {input.strelka2indels} -O v | bcftools view -f PASS -o {output.strelka2VCF} 2> {log.logFile}"

rule VCF2MAFStrelka2:
	input:
		strelka2VCF = strelka2Merged + "{tSamples}.strelka2.passed.vcf"
	output:
		strelka2MAF = strelka2MAFOut + "{tSamples}.strelka2.maf"
	params:
		normID = lambda wildcards: nSamples[wildcards.tSamples]
	message:
		"Running VCF2MAF on Strelka2 calls for {wildcards.tSamples}"
	log:
		logFile = logDir + "{tSamples}.vcf2mafstrelka2.log"
	conda:
		"envs/vcf2maf.yaml"
	shell:
		"vcf2maf.pl --vep-path `dirname $(which variant_effect_predictor.pl)` --input-vcf {input.strelka2VCF} --output-maf {output.strelka2MAF} --vep-data {vepCache} --ref-fasta {refGenome} "
		"--species homo_sapiens --ncbi-build {ncbiBuild} --tumor-id {wildcards.tSamples} --vcf-tumor-id TUMOR --normal-id {params.normID} --vcf-normal-id NORMAL 2> {log.logFile} > {log.logFile}"

rule MergeStrelka1VCFs:
	input:
		strelka1SNVs = strelka1Out + "{tSamples}/results/passed.somatic.snvs.vcf",
		strelka1Indels = strelka1Out + "{tSamples}/results/passed.somatic.indels.vcf"
	output:
		strelka1VCF = strelka1Merged + "{tSamples}.strelka1.passed.vcf"
	log:
		logFile = logDir + "{tSamples}.mergestrelka1.log"
	conda:
		"envs/samtools.yaml"
	shell:
		"cat <(grep ^# {input.strelka1Indels}) <(cat {input.strelka1Indels} {input.strelka1SNVs} | grep -v ^# | sort -V) | bedtools intersect -header -a - -b {capSpaceGATK} | uniq > {output.strelka1VCF} 2> {log.logFile}"

rule VCF2MAFStrelka1:
	input:
		strelka1VCF = strelka1Merged + "{tSamples}.strelka1.passed.vcf"
	output:
		strelka1MAF = strelka1MAFOut + "{tSamples}.strelka1.maf"
	params:
		normID = lambda wildcards: nSamples[wildcards.tSamples]
	message:
		"Running VCF2MAF on Strelka1 calls for {wildcards.tSamples}"
	log:
		logFile = logDir + "{tSamples}.vcf2mafstrelka1.log"
	conda:
		"envs/vcf2maf.yaml"
	shell:
		"vcf2maf.pl --vep-path `dirname $(which variant_effect_predictor.pl)` --input-vcf {input.strelka1VCF} --output-maf {output.strelka1MAF} --vep-data {vepCache} --ref-fasta {refGenome} "
		"--species homo_sapiens --ncbi-build {ncbiBuild} --tumor-id {wildcards.tSamples} --vcf-tumor-id TUMOR --normal-id {params.normID} --vcf-normal-id NORMAL 2> {log.logFile} > {log.logFile}"

rule AugmentMAF:
	input:
		strelka1MAF = strelka1MAFOut + "{tSamples}.strelka1.maf",
		strelka1VCF = strelka1Merged + "{tSamples}.strelka1.passed.vcf"
	output:
		strelka1Aug = strelka1AugOut + "{tSamples}.strelka1.augment.maf"
	message:
		"Running Augment MAF on {input.strelka1MAF}"
	log:
		logFile = logDir + "{tSamples}.augmentmaf.log"
	shell:
		"python {augmentMAFStrelka} --replace -o {output.strelka1Aug} {input.strelka1MAF} {input.strelka1VCF} 2> {log.logFile} > {log.logFile}"

rule vcf2MAFMutect2:
	input:
		mutect2VCF = mutect2Filt + "{tSamples}.mutect2.passed.vcf"
	output:
		mutect2MAF = mutect2MAFOut + "{tSamples}.mutect2.maf"
	params:
		normID = lambda wildcards: nSamples[wildcards.tSamples]
	message:
		"Running VCF2MAF on Mutect2 calls for {wildcards.tSamples}"
	log:
		logFile = logDir + "{tSamples}.vcf2mafmutect2.log"
	conda:
		"envs/vcf2maf.yaml"
	shell:
		"vcf2maf.pl --vep-path `dirname $(which variant_effect_predictor.pl)` --input-vcf {input.mutect2VCF} --output-maf {output.mutect2MAF} --vep-data {vepCache} --ref-fasta {refGenome} "
		"--species homo_sapiens --ncbi-build {ncbiBuild} --tumor-id {wildcards.tSamples} --normal-id {params.normID} 2> {log.logFile} > {log.logFile}"

rule vcf2MAFHaplotypeCaller:
	input:
		hapVCF = haplotypeCallOut + "{allSamples}.hapcaller.vcf"
	output:
		hapMAF = haplotypeMAFOut + "{allSamples}.hapcaller.maf"
	message:
		"Running VCF2MAF on HaplotypeCaller calls for {wildcards.allSamples}"
	log:
		logFile = logDir + "{allSamples}.vcf2mafhaplotypecaller.log"
	conda:
		"envs/vcf2maf.yaml"
	shell:
		"vcf2maf.pl --vep-path `dirname $(which variant_effect_predictor.pl)` --input-vcf {input.hapVCF} --output-maf {output.hapMAF} --vep-data {vepCache} --ref-fasta {refGenome} "
		"--species homo_sapiens --ncbi-build {ncbiBuild} --normal-id {wildcards.allSamples} 2> {log.logFile} > {log.logFile}"

rule vcf2MAFSage:
	input:
		sagePassed = sagePassedOut +  "{tSamples}.sage.passed.vcf" 
	output:
		sageMAF = sageMAFOut + "{tSamples}.sage.maf"
	params:
		normID = lambda wildcards: nSamples[wildcards.tSamples]
	message:
		"Running VCF2MAF on HaplotypeCaller calls for {wildcards.tSamples}"
	log:
		logFile = logDir + "{tSamples}.vcf2mafsage.log"
	conda:
		"envs/vcf2maf.yaml"
	shell:
		"vcf2maf.pl --vep-path `dirname $(which variant_effect_predictor.pl)` --input-vcf {input.sagePassed} --output-maf {output.sageMAF} --vep-data {vepCache} --ref-fasta {refGenome} "
		"--species homo_sapiens --ncbi-build {ncbiBuild} --tumor-id {wildcards.tSamples} --normal-id {params.normID} 2> {log.logFile} > {log.logFile}"

### STEP 6: FILTER STRELKA2 CALLS ####
rule FFPolish:
	input:
		tBAM = sampleDir + "{tSamples}.bam",
		strelka2VCF = strelka2Merged + "{tSamples}.strelka2.passed.vcf"
	output:
		strelka2ffpolish = ffpolishOut + "{tSamples}_filtered.vcf"
	message:
		"Running FFPolish on {wildcards.tSamples}"
	log:
		logFile = logDir + "{tSamples}.strelka2ffpolish.log"
	shell:
		"python {ffpolish} filter -f --cleanup {refGenome} {input.strelka2VCF} {input.tBAM} -o {ffpolishOut} -p {wildcards.tSamples} 2> {log.logFile} 2> {log.logFile}"

rule vcf2MAFFFPolish:
	input:
		strelka2ffpolish = ffpolishOut + "{tSamples}_filtered.vcf"
	output:
		ffpolishMAF = ffpolishMAFOut + "{tSamples}.strelka2.ffpolish.maf"
	params:
		normID = lambda wildcards: nSamples[wildcards.tSamples]
	message:
		"Running VCF2MAF on Strelka2 + FFPolish calls for {wildcards.tSamples}"
	log:
		logFile = logDir + "{tSamples}.vcf2mafffpolish.log"
	conda:
		"envs/vcf2maf.yaml"
	shell:
		"vcf2maf.pl --vep-path `dirname $(which variant_effect_predictor.pl)` --input-vcf {input.strelka2ffpolish} --output-maf {output.ffpolishMAF} --vep-data {vepCache} --ref-fasta {refGenome} "
		"--species homo_sapiens --ncbi-build {ncbiBuild} --tumor-id {wildcards.tSamples} --vcf-tumor-id TUMOR --normal-id {params.normID} --vcf-normal-id NORMAL 2> {log.logFile} > {log.logFile}"


### STEP 7: FILTER SAGE CALLS  ###
rule filterSAGE:
	input:
		tBAM = sampleDir + "{tSamples}.bam",
		sageMAF = sageMAFOut + "{tSamples}.sage.maf"
	output:
		sageFilt = sageFiltOut + "{tSamples}.sage.filt.maf"
	message:
		" Post-filtering SAGE calls for {wildcards.tSamples}"
	log:
		logFile = logDir + "{tSamples}.sage.filt.log"
	shell:
		"python {postFiltScript} -r {refGenome} -f {log.logFile} -b {input.tBAM} -m {input.sageMAF} -o {output.sageFilt} --disable_vaf_filter "



### LAST STEP ###
# Calculate theoretical sensitivity of tumour biopsy mutations
# i.e. are there any reads at all supporting these variants in the various liquid biopsies

rule augmentMAFTumor:
	input:
		tBAM = sampleDir + "{tSamples}.bam",
		nBAM = lambda wildcards: sampleDir + nSamples[wildcards.tSamples] + ".bam",
		strelka1Aug = lambda wildcards: strelka1AugOut + tBiopsySamples[wildcards.tSamples] + ".strelka1.augment.maf"
	output:
		biopsyAugMaf = biopsyAugOut + "{tSamples}.biopsyVar.maf"
	params:
		sampleName = "{tSamples}",
		tumorName = lambda wildcards: tBiopsySamples[wildcards.tSamples]
	message:
		"Running Augment MAF using the tumour MAF on {wildcards.tSamples}"
	log:
		logFile = logDir + "{tSamples}.augmentMAFBiopsy.log"
	conda:
		"envs/augmentmaf.yaml"
	shell:
		"python {augmentMAF} -r --normal-bam {input.nBAM} --tumour-bam {input.tBAM} --maf {input.strelka1Aug} --ignore_overlap {refGenome} {output.biopsyAugMaf} > {log.logFile} 2> {log.logFile};"
		"sed -i 's/{params.tumorName}/{params.sampleName}/' {output.biopsyAugMaf}"
		
