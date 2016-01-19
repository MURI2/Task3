#!/bin/bash

#$ -l walltime=2:00:00,vmem=20gb

# Clones from 113
date
#Change path to where Assembly needs to start
cd /N/dc2/projects/muri2/Task3/<path>
AR=(<Sample_Numbers>)

for i in "${AR[@]}"
do
	cd Sample${i}
	gunzip *
	cat *R1_001.fastq > Sample_${i}_R1.fastq
	cat *R2_001.fastq > Sample_${i}_R2.fastq
        	
	# clean and trim
	echo "#!/bin/bash" > Task3${i}_qc_clean.sh
	echo "" >> Task3${i}_qc_clean.sh
	echo "#$ -l vmem=50gb walltime=48:00:00" >> Task3${i}_qc_clean.sh
	echo "" >> Task3${i}_qc_clean.sh
	echo "cd /N/dc2/projects/muri2/Task3/<path>/Sample${i}" >> Task3${i}_qc_clean.sh
	echo "" >> Task3${i}_qc_clean.sh
	# Removes residual Adaptor Sequence #
	echo "time cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC  -o Sample_${i}_R1_rmadapter.fastq -p Sample_${i}_R2_rmadapter.fastq Sample_${i}_R1.fastq Sample_${i}_R2.fastq" >> Task3${i}_qc_clean.sh
	echo "" >> Task3${i}_qc_clean.sh
	# Filters Reads Based On Quality Scores #
	echo "time cutadapt -q 15,10  -o Sample_${i}_R1_filtered.fastq -p Sample_${i}_R2_filtered.fastq Sample_${i}_R1_rmadapter.fastq Sample_${i}_R2_rmadapter.fastq" >> Task3${i}_qc_clean.sh
	echo "" >> Task3${i}_qc_clean.sh
	# Trims garbage from 5' end of read #
	echo "time cutadapt -u 15 -o Sample_${i}_R1_trimmed.fastq -p Sample_${i}_R2_trimmed.fastq Sample_${i}_R1_filtered.fastq Sample_${i}_R2_filtered.fastq" >> Task3${i}_qc_clean.sh
	echo "" >> Task3${i}_qc_clean.sh
	#Submits assembly script for E. coli Reads#
	echo "qsub -l walltime=20:00:00,vmem=64gb,nodes=1:ppn=4 Task3${i}_EcoliAssemble.sh" >> Task3${i}_qc_clean.sh
	echo "" >> Task3${i}_qc_clean.sh
	echo "qsub -l walltime=20:00:00,vmem=64gb,nodes=1:ppn=4 Task3${i}_RpalAssemble.sh" >> Task3${i}_qc_clean.sh
	echo "" >> Task3${i}_qc_clean.sh
	echo "mkdir fastqc" >> Task3${i}_qc_clean.sh
	echo "" >> Task3${i}_qc_clean.sh
	echo fastqc -o fastqc/ Sample_${i}_R1_trimmed.fastq" >> Task3${i}_qc_clean.sh
	echo "" >> Task3${i}_qc_clean.sh
	echo fastqc -o fastqc/ Sample_${i}_R2_trimmed.fastq" >> Task3${i}_qc_clean.sh
	echo "" >> Task3${i}_qc_clean.sh
	echo "exit" >> Task3${i}_qc_clean.sh

#### Assembly script for E. coli Reads####

	echo "#!/bin/bash" > Task3${i}_EcoliAssemble.sh
	echo "" >> Task3${i}_EcoliAssemble.sh
	echo "#$ -l vmem=50gb walltime=48:00:00 " >> Task3${i}_EcoliAssemble.sh
	echo "" >> Task3${i}_EcoliAssemble.sh
	echo "cd /N/dc2/projects/muri2/Task3/<path>/Sample${i}/" >> Task3${i}_EcoliAssemble.sh
	echo "" >> Task3${i}_EcoliAssemble.sh
	echo "echo \"BWA\" >&2" >> Task3${i}_EcoliAssemble.sh
	echo "time bwa mem /N/dc2/projects/muri2/Task3/RefGenome/E_coli/Ecoli_K12_MG1655.fna Sample_${i}_R1_trimmed.fastq Sample_${i}_R2_trimmed.fastq > Sample_${i}.Ecoli.sam" >> Task3${i}_EcoliAssemble.sh
   	echo "" >> Task3${i}_EcoliAssemble.sh
	echo "qsub -l walltime=2:00:00,vmem=20gb Task3${i}_compress.fastq.sh" >> Task3${i}_EcoliAssemble.sh
	echo "" >> Task3${i}_EcoliAssemble.sh
	echo "echo \"SAMTOOLS\" >&2" >> Task3${i}_EcoliAssemble.sh
	echo "samtools view -bS -T /N/dc2/projects/muri2/Task3/RefGenome/E_coli/Ecoli_K12_MG1655.fna Sample_${i}.Ecoli.sam > Sample_${i}.Ecoli.bam" >> Task3${i}_EcoliAssemble.sh
	echo "" >> Task3${i}_EcoliAssemble.sh
	echo "samtools sort Sample_${i}.Ecoli.bam Sample_${i}.Ecoli.sorted" >> Task3${i}_EcoliAssemble.sh
	echo "" >> Task3${i}_EcoliAssemble.sh
	echo "samtools index Sample_${i}.Ecoli.sorted.bam" >> Task3${i}_EcoliAssemble.sh
	echo "" >> Task3${i}_EcoliAssemble.sh
	echo "" >> Task3${i}_EcoliAssemble.sh
	echo "qsub -l walltime=2:00:00,vmem=20gb Task3${i}_EcoliRetroSeq.sh" >> Task3${i}_EcoliAssemble.sh
	
	        ### GATK for Ecoli###	
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/picard/picard-tools-1.107/AddOrReplaceReadGroups.jar I=Sample_${i}.Ecoli.sorted.bam O=Sample_${i}.Ecoli.sorted.fixed.bam SORT_ORDER=coordinate RGID=Ecoli RGLB=bar RGPL=illumina RGSM=Sample_${i} RGPU=6 CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT" >> Task3${i}_EcoliAssemble.sh
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/picard/picard-tools-1.107/MarkDuplicates.jar I=Sample_${i}.Ecoli.sorted.fixed.bam O=Sample_${i}.Ecoli.sorted.fixed.marked.bam M=Sample_${i}.metrics CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT" >> Task3${i}_EcoliAssemble.sh
	echo "GATK" >&2
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar  -R /N/dc2/projects/muri2/Task3/RefGenome/E_coli/Ecoli_K12_MG1655.fna -I Sample_${i}.Ecoli.sorted.fixed.marked.bam -T RealignerTargetCreator -o Sample_${i}.Ecoli.intervals" >> Task3${i}_EcoliAssemble.sh
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar -R /N/dc2/projects/muri2/Task3/RefGenome/E_coli/Ecoli_K12_MG1655.fna -I Sample_${i}.Ecoli.sorted.fixed.marked.bam -T IndelRealigner -targetIntervals Sample_${i}.Ecoli.intervals --filter_bases_not_stored -o Sample_${i}.Ecoli.sorted.fixed.marked.realigned.bam" >> Task3${i}_EcoliAssemble.sh
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar -nt 4 -T UnifiedGenotyper -R /N/dc2/projects/muri2/Task3/RefGenome/E_coli/Ecoli_K12_MG1655.fna -I Sample_${i}.Ecoli.sorted.fixed.marked.realigned.bam -glm BOTH -rf BadCigar -o Sample_${i}.Ecoli.sorted.fixed.marked.realigned.vcf" >> Task3${i}_EcoliAssemble.sh
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar  -T BaseRecalibrator -I Sample_${i}.Ecoli.sorted.fixed.marked.realigned.bam -R /N/dc2/projects/muri2/Task3/RefGenome/E_coli/Ecoli_K12_MG1655.fna -rf BadCigar --filter_bases_not_stored -knownSites Sample_${i}.Ecoli.sorted.fixed.marked.realigned.vcf -o Sample_${i}.Ecoli.recal_data.grp" >> Task3${i}_EcoliAssemble.sh
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar  -R /N/dc2/projects/muri2/Task3/RefGenome/E_coli/Ecoli_K12_MG1655.fna -I Sample_${i}.Ecoli.sorted.fixed.marked.realigned.bam -T PrintReads -rf BadCigar  -o Sample_${i}.Ecoli.mapped.bam -BQSR Sample_${i}.Ecoli.recal_data.grp" >> Task3${i}_EcoliAssemble.sh
	#java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar -nt 4 -T UnifiedGenotyper -R ../../Ecoli_RefGenome/Ecoli_K12_MG1655.fna -I Sample_${i}.sorted.fixed.marked.realigned.bam -rf BadCigar -o Sample_${i}.sorted.fixed.marked.realigned.vcf" >> Task3${i}_EcoliAssemble.sh
	echo "" >> Task3${i}_EcoliAssemble.sh
	echo "exit" >> Task3${i}_EcoliAssemble.sh
	
	
	#### Assembly script for R. palustris Reads####

	echo "#!/bin/bash" > Task3${i}_RpalAssemble.sh
	echo "" >> Task3${i}_RpalAssemble.sh
	echo "#$ -l vmem=50gb walltime=48:00:00 " >> Task3${i}_RpalAssemble.sh
	echo "" >> Task3${i}_RpalAssemble.sh
	echo "cd /N/dc2/projects/muri2/Task3/<path>/Sample${i}/" >> Task3${i}_RpalAssemble.sh
	echo "" >> Task3${i}_RpalAssemble.sh
	echo "echo \"BWA\" >&2" >> Task3${i}_RpalAssemble.sh
	echo "time bwa mem /N/dc2/projects/muri2/Task3/RefGenome/R_palustris/<R_Pal_Reference> Sample_${i}_R1_trimmed.fastq Sample_${i}_R2_trimmed.fastq > Sample_${i}.Rpal.sam" >> Task3${i}_RpalAssemble.sh
   	echo "" >> Task3${i}_RpalAssemble.sh
	echo "echo \"SAMTOOLS\" >&2" >> Task3${i}_RpalAssemble.sh
	echo "samtools view -bS -T /N/dc2/projects/muri2/Task3/RefGenome/R_palustris/<R_Pal_Reference> Sample_${i}.Rpal.sam > Sample_${i}.Rpal.bam" >> Task3${i}_RpalAssemble.sh
	echo "" >> Task3${i}_RpalAssemble.sh
	echo "samtools sort Sample_${i}.Rpal.bam Sample_${i}.Rpal.sorted" >> Task3${i}_RpalAssemble.sh
	echo "" >> Task3${i}_RpalAssemble.sh
	echo "samtools index Sample_${i}.Rpal.sorted.bam" >> Task3${i}_RpalAssemble.sh
	echo "" >> Task3${i}_RpalAssemble.sh
	echo "" >> Task3${i}_RpalAssemble.sh
	echo "qsub -l walltime=2:00:00,vmem=20gb Task3${i}_RpalRetroSeq.sh" >> Task3${i}_RpalAssemble.sh
	
	        ### GATK for Rpal###	
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/picard/picard-tools-1.107/AddOrReplaceReadGroups.jar I=Sample_${i}.Rpal.sorted.bam O=Sample_${i}.Rpal.sorted.fixed.bam SORT_ORDER=coordinate RGID=Rpal RGLB=bar RGPL=illumina RGSM=Sample_${i} RGPU=6 CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT" >> Task3${i}_RpalAssemble.sh
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/picard/picard-tools-1.107/MarkDuplicates.jar I=Sample_${i}.Rpal.sorted.fixed.bam O=Sample_${i}.Rpal.sorted.fixed.marked.bam M=Sample_${i}.Rpal.metrics CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT" >> Task3${i}_RpalAssemble.sh
	echo "GATK" >&2
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar  -R /N/dc2/projects/muri2/Task3/RefGenome/R_palustris/<R_Pal_Reference> -I Sample_${i}.Rpal.sorted.fixed.marked.bam -T RealignerTargetCreator -o Sample_${i}.Rpal.intervals" >> Task3${i}_RpalAssemble.sh
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar -R /N/dc2/projects/muri2/Task3/RefGenome/R_palustris/<R_Pal_Reference> -I Sample_${i}.Rpal.sorted.fixed.marked.bam -T IndelRealigner -targetIntervals Sample_${i}.Rpal.intervals --filter_bases_not_stored -o Sample_${i}.Ecoli.sorted.fixed.marked.realigned.bam" >> Task3${i}_RpalAssemble.sh
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar -nt 4 -T UnifiedGenotyper -R /N/dc2/projects/muri2/Task3/RefGenome/R_palustris/<R_Pal_Reference> -I Sample_${i}.Rpal.sorted.fixed.marked.realigned.bam -glm BOTH -rf BadCigar -o Sample_${i}.Rpal.sorted.fixed.marked.realigned.vcf" >> Task3${i}_RpalAssemble.sh
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar  -T BaseRecalibrator -I Sample_${i}.Rpal.sorted.fixed.marked.realigned.bam -R /N/dc2/projects/muri2/Task3/RefGenome/R_palustris/<R_Pal_Reference> -rf BadCigar --filter_bases_not_stored -knownSites Sample_${i}.Rpal.sorted.fixed.marked.realigned.vcf -o Sample_${i}.Rpal.recal_data.grp" >> Task3${i}_RpalAssemble.sh
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar  -R /N/dc2/projects/muri2/Task3/RefGenome/R_palustris/<R_Pal_Reference> -I Sample_${i}.Rpal.sorted.fixed.marked.realigned.bam -T PrintReads -rf BadCigar  -o Sample_${i}.Ecoli.mapped.bam -BQSR Sample_${i}.Rpal.recal_data.grp" >> Task3${i}_RpalAssemble.sh
	#java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar -nt 4 -T UnifiedGenotyper -R ../../Ecoli_RefGenome/Ecoli_K12_MG1655.fna -I Sample_${i}.sorted.fixed.marked.realigned.bam -rf BadCigar -o Sample_${i}.sorted.fixed.marked.realigned.vcf" >> Task3${i}_EcoliAssemble.sh
	echo "" >> Task3${i}_RpalAssemble.sh
	echo "exit" >> Task3${i}_RpalAssemble.sh
	
 	# compress.fastq again
	echo "#!/bin/bash" > Task3${i}_compress.fastq.sh
	echo "" >> Task3${i}_compress.fastq.sh
	echo "#$ -l vmem=50gb walltime=48:00:00 " >> Task3${i}_compress.fastq.sh
	echo "" >> Task3${i}_compress.fastq.sh
	echo "cd /N/dc2/projects/muri2/Task3/<path>/Sample${i}/" >> Task3${i}_compress.fastq.sh
	#echo "bzip2 Sample_${i}_R*.fastq" >> Task3${i}_compress.fastq.sh
	echo "rm Sample_${i}_R*_filtered.fastq" >> Task3${i}_compress.fastq.sh
	echo "rm Sample_${i}_R*_rmadapter.fastq" >> Task3${i}_compress.fastq.sh
	#echo "bzip2 Sample_${i}_R*_trimmed.fastq" >> Task3${i}_compress.fastq.sh
	echo "" >> Task3${i}_compress.fastq.sh
	echo "exit" >> Task3${i}_compress.fastq.sh
	
	##### call novel E.coli insertion sequences######
	echo "#!/bin/bash" > Task3${i}_EcoliRetroSeq.sh
	echo "" >> Task3${i}_EcoliRetroSeq.sh
	echo "#$ -l vmem=50gb walltime=2:00:00 " >> Task3${i}_EcoliRetroSeq.sh
	echo "" >> Task3${i}_EcoliRetroSeq.sh
	echo "cd /N/dc2/projects/muri2/Task3/<path>/Sample${i}/"  >>Task3${i}_EcoliRetroSeq.sh
	echo "mkdir E_coli_InsSeq" >> Task3${i}_EcoliRetroSeq.sh
	echo "cd E_coli_InsSeq" >>Task3${i}_EcoliRetroSeq.sh
	echo "time bwa mem /N/dc2/projects/muri2/Task3/RefGenome/E_coli/InsSeq/RefGenome/Ecoli_K12_MG1655.fna ../Sample_${i}_R1_trimmed.fastq ../Sample_${i}_R2_trimmed.fastq > Sample_${i}.sam" >> Task3${i}_EcoliRetroSeq.sh
	echo "samtools view -bS -T /N/dc2/projects/muri2/Task3/RefGenome/E_coli/InsSeq/RefGenome/Ecoli_K12_MG1655.fna Sample_${i}.sam > Sample_${i}.bam" >> Task3${i}_EcoliRetroSeq.sh
	echo "samtools sort Sample_${i}.bam Sample_${i}.sorted" >> Task3${i}_EcoliRetroSeq.sh
	echo "samtools index Sample_${i}.sorted.bam" >> Task3${i}_EcoliRetroSeq.sh
	echo "perl /N/dc2/scratch/megbehri/SAM_MURI/Tools/RetroSeq/bin/retroseq.pl -discover -eref /N/dc2/projects/muri2/Task3/RefGenome/E_coli/InsSeq/Ecoli_IS_RefFile.txt -bam Sample_${i}.sorted.bam -output Sample_${i}.IS.Reads -align" >> Task3${i}_EcoliRetroSeq.sh
	echo "perl /N/dc2/scratch/megbehri/SAM_MURI/Tools/RetroSeq/bin/retroseq.pl -call -bam Sample_${i}.sorted.bam -ref /N/dc2/projects/muri2/Task3/RefGenome/E_coli/InsSeq/RefGenome/Ecoli_K12_MG1655.fna -output Sample_${i}.Ecoli.IS -input Sample_${i}.IS.Reads -hets" >>Task3${i}_EcoliRetroSeq.sh
	echo "" >> Task3${i}_EcoliRetroSeq.sh
	echo "exit" >> Task3${i}_EcoliRetroSeq.sh
	
	##### call novel R. palustris insertion sequences######
	echo "#!/bin/bash" > Task3${i}_RpalRetroSeq.sh
	echo "" >> Task3${i}_RpalRetroSeq.sh
	echo "#$ -l vmem=50gb walltime=2:00:00 " >> Task3${i}_RpalRetroSeq.sh
	echo "" >> Task3${i}_RpalRetroSeq.sh
	echo "cd /N/dc2/projects/muri2/Task3/<path>/Sample${i}/"  >>Task3${i}_RpalRetroSeq.sh
	echo "mkdir R_pal_InsSeq" >> Task3${i}_RpalRetroSeq.sh
	echo "cd R_pal_InsSeq" >>Task3${i}_RpalRetroSeq.sh
	echo "time bwa mem /N/dc2/projects/muri2/Task3/RefGenome/R_palustris/InsSeq/RefGenome/<R_Pal_Reference> ../Sample_${i}_R1_trimmed.fastq ../Sample_${i}_R2_trimmed.fastq > Sample_${i}.sam" >> Task3${i}_RpalRetroSeq.sh
	echo "samtools view -bS -T /N/dc2/projects/muri2/Task3/RefGenome/R_palustris/InsSeq/RefGenome/<R_Pal_Reference> Sample_${i}.sam > Sample_${i}.bam" >> Task3${i}_RpalRetroSeq.sh
	echo "samtools sort Sample_${i}.bam Sample_${i}.sorted" >> Task3${i}_RpalRetroSeq.sh
	echo "samtools index Sample_${i}.sorted.bam" >> Task3${i}_RpalRetroSeq.sh
	echo "perl /N/dc2/scratch/megbehri/SAM_MURI/Tools/RetroSeq/bin/retroseq.pl -discover -eref /N/dc2/projects/muri2/Task3/RefGenome/R_palustris/InsSeq/Rpal_IS_RefFile.txt -bam Sample_${i}.sorted.bam -output Sample_${i}.IS.Reads -align" >> Task3${i}_RpalRetroSeq.sh
	echo "perl /N/dc2/scratch/megbehri/SAM_MURI/Tools/RetroSeq/bin/retroseq.pl -call -bam Sample_${i}.sorted.bam -ref /N/dc2/projects/muri2/Task3/RefGenome/R_palustris/InsSeq/RefGenome/<R_Pal_Reference> -output Sample_${i}.Rpal.IS -input Sample_${i}.IS.Reads -hets" >>Task3${i}_RpalRetroSeq.sh
	echo "" >> Task3${i}_RpalRetroSeq.sh
	echo "exit" >> Task3${i}_RpalRetroSeq.sh
		
	chmod u+x Task3${i}_qc_clean.sh
	chmod u+x Task3${i}_EcoliAssemble.sh
	chmod u+x Task3${i}_RpalAssemble.sh
	chmod u+x Task3${i}_compress.fastq.sh
	chmod u+x Task3${i}_EcoliRetroSeq.sh
	chmod u+x Task3${i}_RpalRetroSeq.sh
	
	qsub -l walltime=2:00:00,vmem=20gb Task3${i}_qc_clean.sh
	cd ..
	
	done

