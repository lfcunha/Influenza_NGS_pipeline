#!/usr/bin/make -Rf

# 16.10.2012 14:00:54 EDT
# Harm van Bakel <hvbakel@gmail.com>

###############
# DEFINE ARGS #
###############

type        ?=A
np          ?=12
kmersize    ?=27
kmercov     ?=2
varthr      ?=0.05
awkt        :=awk -F '\t' -v OFS='\t'
sortt       :=sort -t $$'\t'
currentdir   := $(notdir $(shell pwd))

##############
# CHECK ARGS #
##############

ifeq ($(type),A)
	inflblatdb =/sc/orga/projects/vanbah01a/reference-databases/influenza/IAVreference_cdhit.2bit
	inflblatbw =/sc/orga/projects/vanbah01a/reference-databases/influenza/IAVreference_cdhit_chimera_800bp_noN
else
	inflblatdb =/sc/orga/projects/vanbah01a/reference-databases/influenza/IBVreference_cdhit_noN.2bit
	inflblatbw =/sc/orga/projects/vanbah01a/reference-databases/influenza/IBVreference_cdhit_noN
endif

ifndef in
$(error missing argument in)
else
name?=$(in:.fastq.gz=)
endif

ifndef PICARD_HOME
$(error PICARD_HOME environment variable not set, please run the pipeline through the provided 'run-flu-assembly' wrapper script)
endif
ifndef GATK_JAR
$(error GATK_JAR environment variable not set, please run the pipeline through the provided 'run-flu-assembly' wrapper script)
endif

#########
# RULES #
#########

assemble: cutadapt chimeradetect inchworm cap3 complete mapback coverage

cutadapt: $(name)_cutadapt.fastq.gz $(name)_cutadapt_nochreads.fastq.gz

chimeradetect: $(name)_chimeradetect.fa $(name)_chimeradetect.chimerajunct $(name)_cutadapt_nochimera.fastq.gz

inchworm: $(name)_assembly.fa $(name)_assembly.chimerajunct

cap3: $(name)_assembly_incomplete_cap3.fa $(name)_assembly_cap3.fa $(name)_assembly_cap3_cdhit.fa

complete: $(name)_complete.fa $(name)_complete.info

mapback: $(name)_assembly_cap3_cdhit_bwa.bam $(name)_assembly_cap3_cdhit_bwa.bam.bai

coverage: $(name)_assembly_cap3_cdhit.info.cov $(name)_assembly_cap3_cdhit.plotvar $(name)_assembly_cap3_cdhit.report.pdf

report: $(name)_curated_bwa.bam $(name)_curated_bwa.bam.bai $(name)_curated.info.cov $(name)_curated.plotvar \
        $(name)_curated.plotcov $(name)_curated.plotdip $(name)_curated.report.pdf $(name)_curated.variants.calls.txt \
        $(name)_curated_star.bam $(name)_curated_star.bam.bai

lofreq: $(name)_curated_bwa.bam $(name)_curated_bwa.bam.bai $(name)_curated_bwa_loseq_knownsites_filt.vcf \
        $(name)_curated_bwa_rg_recal.table $(name)_curated_bwa_rg_bqsr.bam $(name)_curated_bwa_rg_bqsr.bam.bai \
        $(name)_curated_bwa_loseq_snvindel.vcf

annotate:
	iavsat_annotate.sh $(name) 
	iavsat_translate.py $(name)_cds.fa

log:
ifeq ("$(currentdir)", "curated")
	process_assembly.py --anno $(name).anno --var $(name)_curated.variants.calls.txt
else
	process_assembly.py --cov $(name)_assembly_cap3_cdhit.info.cov
endif


varpatch: $(name)_curated_bwa_loseq_snvindel_noconsensus.vcf.gz $(name)_curated_varpatch.fa

upload: annotate
	cripdb-submit.php $(name)

.PHONY: clean compact

compact:
	rm -rf $(name)_chimeradetect $(name)_assembly $(name)_cutadapt.fastq.gz $(name)_cutadapt_nochimera.fastq.gz \
	$(name)_cutadapt.fastq.stat $(name)_assembly_cap3.info $(name)_assembly_cap3.bam.err $(name)_assembly.fa.err \
	$(name)_assembly.chimerajunct $(name)_assembly.fa

clean: compact
	rm -rf $(name)_assembly* $(name)_chimeradetect* $(name)_complete* 

############################
# REMOVE ADAPTER SEQUENCES #
############################

%_cutadapt.tmp: %.fastq.gz
	cutadapt -f fastq --match-read-wildcards -e 0.1 -O 6 -g ^GAGCTAGTCTG -g ^GTCGAGCTCG -g ^GTTACGCGCC -g ^GGGGGG -g ^CGGGTTATT -g ^GGTAACGCGTGATC -a GATCGGAAGAGCACACGTCT -b ACACTCTTTCCCTACACGACGCTCTTCCGATCT -b GCCAGAGCCGTAAGGACGACTTGGCGAGAAGGCTAGA $< > $@
%_cutadapt.fastq.gz: %_cutadapt.tmp
	fastq-quality-trimmer.pl -min_len 40 -min_qual 4 -max_fail_qual 20 -max_fail 15 -output $(patsubst %.tmp,%,$<) $<
%_cutadapt.chreadids: %_cutadapt.fastq.gz $(inflblatbw).1.ebwt
	#influenza_detect-chimeric-reads.pl -p $(np) -f $< -b $(inflblatbw) > $@
	rm -rf ./$(name)_chreaddetect/
	mkdir ./$(name)_chreaddetect
	STAR \
        --chimSegmentMin        8 \
        --chimJunctionOverhangMin 8 \
        --alignSJoverhangMin 8 \
        --outSJfilterReads Unique \
        --outSJfilterOverhangMin 8 8 8 8 \
        --outSJfilterCountUniqueMin 2 2 2 2 \
        --outSJfilterCountTotalMin 2 2 2 2 \
        --outSJfilterDistToOtherSJmin 0 0 0 0 \
        --genomeDir 	$(inflblatbw) \
        --runThreadN $(np) \
        --outReadsUnmapped unmapped \
        --outStd SAM \
        --outSAMmode Full \
        --outFileNamePrefix ./$(name)_chreaddetect/refmap_star \
        --readFilesCommand zcat \
        --readFilesIn $< | samtools view -bS - \
        > ./$(name)_chreaddetect/refmap_star.bam
	samtools view ./$(name)_chreaddetect/refmap_star.bam | $(awkt) '$$6~/[0-9][0-9][0-9]+N/{print $$1}' | sort | uniq > $@
	$(awkt) '$$1!~/^@/{print $$1}' ./$(name)_chreaddetect/refmap_starChimeric.out.sam | uniq >> $@
	rm -rf ./$(name)_chreaddetect
%_cutadapt_nochreads.fastq.gz: %_cutadapt.fastq.gz %_cutadapt.chreadids
	fastq-filter-by-id.pl -f $(word 1, $+) -m $(word 2, $+) -v | pigz -b 2048 -p $(np) > $@

####################################################
# STAGE 1: INCHWORM ASSEMBLY FOR CHIMERA DETECTION #
####################################################

%_chimeradetect: %_cutadapt_nochreads.fastq.gz
	Trinity --no_run_chrysalis --no_run_butterfly --CPU $(np) --KMER_SIZE $(kmersize) --min_kmer_cov $(kmercov) --no_cleanup --seqType fq --JM 5G --output $@ --single $<
%_chimeradetect.assembly.fa: %_chimeradetect
	filter-inchworm-contigs.pl -minlength 100 -mincount 250 $</inchworm.K$(kmersize).L25.DS.fa > $@
%_chimeradetect.fa %_chimeradetect.chimerajunct: %_chimeradetect.assembly.fa $(inflblatdb)
	influenza_annotate-contigs.pl -f $(word 1, $+) -r $(word 2, $+) -t $(type) -o $(patsubst %.assembly.fa,%,$<)
%_cutadapt_nochimera.fastq.gz: %_cutadapt_nochreads.fastq.gz %_chimeradetect.chimerajunct
	fastq-filter-by-seq.pl -f $(word 1, $+) -m $(word 2, $+) -r -e 1 -v | pigz -b 2048 -p $(np) > $@

############################################################
# STAGE 2: INCHWORM ASSEMBLY AFTER REMOVING CHIMERIC READS #
############################################################

%_assembly: %_cutadapt_nochimera.fastq.gz
	Trinity --no_run_chrysalis --no_run_butterfly --CPU $(np) --KMER_SIZE $(kmersize) --min_kmer_cov $(kmercov) --no_cleanup --seqType fq --JM 5G --output $@ --single $<
%_assembly.assembly.fa: %_assembly
	filter-inchworm-contigs.pl -minlength 100 -mincount 250 $</inchworm.K$(kmersize).L25.DS.fa > $@
%_assembly.fa %_assembly.chimerajunct: %_assembly.assembly.fa $(inflblatdb)
	influenza_annotate-contigs.pl -f $(word 1, $+) -r $(word 2, $+) -t $(type) -o $(patsubst %.assembly.fa,%,$<)

#############################################
# STAGE3: MERGE SEGMENT FRAGMENTS WITH CAP3 #
#############################################

# Try to further merge incomplete viral fragments with cap3
%_assembly_complete.ids: %_assembly.fa
	grep -P '(complete\|contiguous)|(complete\|chimera)' $< | perl -pi -e 's/^>//' > $@
%_assembly_complete.fa: %_assembly.fa %_assembly_complete.ids
	fasta-filter-by-id.pl -f $(word 1, $+) -m $(word 2, $+) > $@
%_assembly_incomplete.ids: %_assembly.fa
	grep 'fragment|contiguous' $< | perl -pi -e 's/^>//' > $@
%_assembly_incomplete.fa: %_assembly.fa %_assembly_incomplete.ids
	fasta-filter-by-id.pl -f $(word 1, $+) -m $(word 2, $+) > $@
%_assembly_incomplete.fa.cap.contigs %_assembly_incomplete.fa.cap.singlets: %_assembly_incomplete.fa
	-cap3 $< -i 30  -j 41  -o 20  -s 400 -r > $<.err; \
	touch $(patsubst %.fa,%.fa.cap.contigs,$<); \
	touch $(patsubst %.fa,%.fa.cap.singlets,$<);
%_assembly_incomplete.fa.cap.contigs.fa: %_assembly_incomplete.fa.cap.contigs $(inflblatdb)
	influenza_annotate-contigs.pl -f $(word 1, $+) -r $(word 2, $+) -t $(type) -o $<
%_assembly_incomplete_cap3.fa: %_assembly_incomplete.fa.cap.singlets %_assembly_incomplete.fa.cap.contigs.fa
	cat $(word 1, $+) $(word 2, $+) > $@; \
	rm -f $(patsubst %.singlets,%.*,$<)
%_assembly_cap3.fa: %_assembly_complete.fa %_assembly_incomplete_cap3.fa
	cat $(word 1, $+) $(word 2, $+) | fasta2tabbed.pl - | $(sortt) -k2,2 -k1,1 | $(awkt) 'out!=$$2{print $$1,$$2; out=$$2}' \
	   | $(sortt) -k1,1 | $(awkt) '{print ">" $$1 "\n" $$2}' | fasta-reflow.pl - > $@
%_assembly_cap3_cdhit.fa: %_assembly_cap3.fa
	cd-hit-est -c 0.98 -i $< -o $@

####################################
# GATHER FULLY ASSEMBLED FRAGMENTS #
####################################

%_complete.ids: %_assembly_cap3_cdhit.fa
	grep -P '(complete\|contiguous)|(complete\|chimera)' $< | perl -pi -e 's/^>//' > $@
%_complete.fa: %_assembly_cap3_cdhit.fa %_complete.ids
	fasta-filter-by-id.pl -f $(word 1, $+) -m $(word 2, $+) > $@
%_complete.info: %_complete.fa
	fasta-get-info.pl $< > $@

##############################
# MAP READS BACK TO ASSEMBLY #
##############################

# BAM indexing
%.bam.bai: %.bam
	samtools index $<

# BWA mapping (for read coverage and variant analysis)
%_bwa_index.amb: %.fa
	bwa index -p $(patsubst %.amb,%,$@) $< 
%_bwa_tmp.bam: %_bwa_index.amb $(name)_cutadapt.fastq.gz
	bwa mem -t $(np) $(patsubst %.amb,%,$<) $(word 2, $+) | samtools view -S -bq 1 - > $@
	rm -f $(patsubst %.amb,%*,$<)
%_bwa.bam: %_bwa_tmp.bam
	samtools sort $< $(patsubst %.bam,%,$@)

# STAR mapping (for DIP analysis)
%_star_tmp.bam: $(name)_cutadapt.fastq.gz %.fa
	rm -rf ./$(name)_index/
	mkdir ./$(name)_index
	STAR --runMode genomeGenerate --genomeDir ./$(name)_index --genomeFastaFiles $(word 2, $+) --runThreadN 8 --genomeSAindexNbases 6
	rm -f Log.out
	STAR \
        --chimSegmentMin        8 \
        --chimJunctionOverhangMin 8 \
        --alignSJoverhangMin 8 \
        --outSJfilterReads Unique \
        --outSJfilterOverhangMin 8 8 8 8 \
        --outSJfilterCountUniqueMin 2 2 2 2 \
        --outSJfilterCountTotalMin 2 2 2 2 \
        --outSJfilterDistToOtherSJmin 0 0 0 0 \
        --genomeDir ./$(name)_index \
        --runThreadN $(np) \
        --outReadsUnmapped unmapped \
        --outStd SAM \
        --outSAMmode Full \
        --outFileNamePrefix $(patsubst %_tmp.bam,%_,$@) \
        --readFilesCommand zcat \
        --readFilesIn $< | samtools view -bS - \
        > $@
	rm -rf ./$(name)_index
%_star.bam:	%_star_tmp.bam
	samtools sort $< $(patsubst %_tmp.bam,%,$<)
	rm -f $(name)_curated_tmp.bam
%.plotdip: %_star.bam %.plotcov
	get-DIP-breakpoints-star.py $(patsubst %.bam,%_SJ.out.tab,$<) $(word 2, $+) $@

############################################################
# GET COVERAGE AND VARIANT REPORT FOR ALL TRIMMED SEGMENTS #
############################################################

# Coverage table
%.info: %.fa
	fasta-get-info.pl $< > $@
%.fulllength.bed: %.info
	$(awkt) '$$1!~/^#/{print $$1,0,$$2,$$1}' $< > $@
%.fulllength.tmp: %_bwa.bam %.fulllength.bed
	coverageBed -abam $< -b $(word 2, $+) | cut -f 1,5,6,8 > $@
%.fulllength.cov: %.fulllength.tmp
	add-header $< '#id' full_count full_bp_covered full_fraction_covered; mv -f $< $@
%.info.cov: %.info %.fulllength.cov
	join-by-ids -a 1 -1 1 -2 1 $< $(word 2, $+) > $@;

# Report coverage and variants
%.plotcov: %_bwa.bam %.fulllength.bed
	coverageBed -d -abam $< -b $(word 2, $+) | cut -f 4-6 > $@
%.plotvar: %_bwa.bam %.fa
	influenza_count-variant-bases.pl -b $< -f $(word 2, $+) > $@
%.variants.calls.txt: %.plotvar
	${awkt} '($$3>=0.05 && $$2!=4) || $$1~/^#/' $< > $@
%.report.pdf: %.plotcov %.plotvar %.plotdip
	influenza_finalreport.R -l 100 -i $< -v $(word 2, $+) -d $(word 3, $+) -o $(patsubst %.pdf,%,$@) -t $(varthr)

############################################
# LOSEQ VARIANT CALLING (INCLUDING INDELS) #
############################################

# Call and filter first-pass variant sites
%_curated_bwa_loseq_knownsites.vcf: %_curated.fa %_curated_bwa.bam %_curated_bwa.bam.bai
	lofreq call-parallel --pp-threads $(np) -f $< -o $@ $(word 2, $+)
%_curated_bwa_loseq_knownsites_filt.vcf: %_curated_bwa_loseq_knownsites.vcf
	bcftools filter -e "AF<0.005" $< 2> /dev/null > $@

# Recalibrate base qualities using picard and GATK
%_curated.dict: %_curated.fa
	java -Xmx2g -jar $(PICARD_HOME)/CreateSequenceDictionary.jar R=$< O=$@
%_curated_bwa_rg.bam: %_curated_bwa.bam
	java -Xmx2g -jar $(PICARD_HOME)/AddOrReplaceReadGroups.jar I=$< O=$@ RGLB=A RGPL=illumina RGPU=run RGSM=$(patsubst %_curated_bwa.bam,%,$<)
%_curated_bwa_rg_recal.table: %_curated.fa %_curated_bwa_rg.bam %_curated_bwa_rg.bam.bai %_curated_bwa_loseq_knownsites_filt.vcf %_curated.dict
	java -Xmx4g -jar $(GATK_JAR) -T BaseRecalibrator -R $< -I $(word 2, $+) --filter_reads_with_N_cigar -knownSites $(word 4, $+) -o $@
%_curated_bwa_rg_bqsr.bam: %_curated.fa %_curated_bwa_rg.bam %_curated_bwa_rg_recal.table %_curated.dict
	java -jar $(GATK_JAR) -T PrintReads -R $< -I $(word 2, $+) --filter_reads_with_N_cigar -BQSR $(word 3, $+) -o $@

# Do final calling including indels
%_curated_bwa_loseq_snvindel.vcf: %_curated.fa %_curated_bwa_rg_bqsr.bam %_curated_bwa_rg_bqsr.bam.bai
	lofreq call-parallel --call-indels --pp-threads $(np) -f $< -o $@ $(word 2, $+)

# Patch consensus with variants scoring above average
%_noconsensus.vcf: %.vcf
	bcftools filter -e "AF<0.5" $< 2> /dev/null > $@
%.vcf.gz: %.vcf
	bgzip -c $< > $@
%.vcf.gz.tbi: %.vcf.gz
	tabix -p vcf $<
%_curated_varpatch.fa: %_curated.fa %_curated_bwa_loseq_snvindel_noconsensus.vcf.gz %_curated_bwa_loseq_snvindel_noconsensus.vcf.gz.tbi
	cat $< | vcf-consensus $(word 2, $+) | fasta-reflow.pl - > $@
