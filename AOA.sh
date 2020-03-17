cd /media/dave/storage/macronutrient/BFI_paper/data/sequences/PKG_ENQ_723_data_transfer_RIVERS_FUNCTIONAL_GENES_NERCMACRO

# cut primer sequences
# remove any sequences that didn't contain both the forward and reverse primer
for f in *_b_*R1.fastq.gz; do
cutadapt -g "ATGGTCTGGCTWAGACG;min_overlap=16" -G "GCCATCCATCTGTATGTCCA;min_overlap=19" -j 0 --discard-untrimmed -o /media/dave/storage/macronutrient/BFI_paper/data/sequences/AOA/${f%L001_R1.fastq.gz}R1_primer_trimmed.fastq -p /media/dave/storage/macronutrient/BFI_paper/data/sequences/AOA/${f%L001_R1.fastq.gz}R2_primer_trimmed.fastq $f ${f%R1.fastq.gz}R2.fastq.gz
done

cd /media/dave/storage/macronutrient/BFI_paper/data/sequences/AOA

# make directories for merged sequence files, html output and discarded sequences
mkdir {qualTrimmedSeqs,htmlReports}

# quality filter
# remove seqs if > 20% is under q20
# below 200 bases after trimming
# paired reads won't merge - too long an amplicon
for f in *R1_primer_trimmed.fastq; do
  fastp -i $f -q 15 -u 20 -l 200 -r -y -A -5 -n 0 -o qualTrimmedSeqs/${f%R1_primer_trimmed.fastq}qualTrimmed.fastq -h htmlReports/${f%R1_primer_trimmed.fastq}quality_report.html
done

cd qualTrimmedSeqs

### Error correct all seqs ###
for f in *.fastq;
	do ~/bioinformatics/SPAdes-3.13.0-Linux/bin/spades.py -o ${f%q_trimmed_merged.fastq}_errorCorrected --only-error-correction -s $f -t 8 -m 32 --disable-gzip-output;
done

# copy all error corrected sequences and convert to fastq and length filter
mkdir allErrorCorrected
mv *errorCorrected/corrected/*.fastq allErrorCorrected

cd allErrorCorrected

# convert to fasta format
for f in *.fastq; do
  seqtk seq -a $f > ${f%_qualTrimmed.00.0_0.cor.fastq}.fna
done

### NEED TO LOOK AT LENGTH FILTERS
# get read lengths for one file
R
library(data.table)
dat <- fread("grep '^[^>]' FEB-AS2_2_b_GGACTCCT-GTAAGGAG.fna | awk '{print length}'")#"
summary(dat$V1)
q()

# take all reads between > 240 bp
for f in *.fna; do
awk '!/^>/ { next } { getline seq } length(seq) >= 280 {print $0 "\n" seq}' $f > ${f%.fna}_len_filtered.fna
done

# add sample labels to fasta headers
for f in *_len_filtered.fna; do
	samp=">$(echo $f | awk -F "_b_" '{print $1}' | sed "s/-/x/g; s/_/x/g")-"
	sed -i "s/>/$samp/" $f
done

# calculate final lib sizes
R
library(data.table)
dat <- fread("grep -c '^>' *_len_filtered.fna", sep = ":")
dat[order(V2, decreasing = F), ]
paste(dat[V2 <= 1800, V1], collapse = ",")
q()

# discard samples with fewer than 1800 seqs
rm -f {AUG-AS3_1_b_TCCTGAGC-AAGGAGTA_len_filtered.fna,AUG-CE1_1_b_TCCTGAGC-CTAAGCCT_len_filtered.fna,AUG-CE1_2_b_GGACTCCT-CTAAGCCT_len_filtered.fna,AUG-CE1_3_b_TAGGCATG-CTAAGCCT_len_filtered.fna,FEB-CW2_2_b_GGACTCCT-ACTGCATA_len_filtered.fna,FEB-GN1_1_b_CTCTCTAC-GTAAGGAG_len_filtered.fna}

# concatenate all length filtered fasta files into one
cat *_len_filtered.fna > allSequences.fna

# get total final lib sizes
grep -c "^>" allSequences.fna # 242715

# dereplicate reads by sample to preserve abundance info
~/bioinformatics/usearch11.0.667_i86linux32 -fastx_uniques_persample allSequences.fna -fastaout uniq_DNA_seqs.fna -minuniquesize 1 -sizeout

# use subsampled AOA db from PRINCe work
mkdir framebotOut

java -Xmx16g -jar ~/bioinformatics/RDPTools/FrameBot.jar framebot -o framebotOut/AOA_amoA -N -i 0.5 -l 50 /media/dave/storage/PRINCE/data/sequenceData/fungeneDBs/AOA_amoA/subSampled_AOA.fna uniq_DNA_seqs.fna

cd framebot

awk '{ if ($0 !~ />/) {print toupper($0)} else {print $0} }' AOA_amoA_corr_prot.fasta > uc_AOA_amoA_corr_prot.fasta

# convert fs corrected aa sequences to upper case for usearch
~/bioinformatics/usearch11.0.667_i86linux32 -fastx_uniques uc_AOA_amoA_corr_prot.fasta -fastaout AAVs.fasta -minuniquesize 1 -relabel AAV

# convert to amino acid variant OTU table
# might have to split into chunks and recombine afterwards
~/bioinformatics/usearch11.0.667_i86linux32 -otutab uc_AOA_amoA_corr_prot.fasta -otus AAVs.fasta -id 1 -minsize 0 -otutabout AOA_AAV.txt -qmask none -dbmask none

# use vsearch to dereplicate frameshift corrected nucleotide sequences
vsearch --derep_fulllength AOA_amoA_corr_nucl.fasta --output derep_nucl_seqs.fna --minuniquesize 2 --relabel OTU --fasta_width 0 --sizein --sizeout

# sort by cluster size
vsearch --sortbysize derep_nucl_seqs.fna --output derepSorted.fna

# pick otu centroids
vsearch --cluster_size derepSorted.fna --id 0.97 --centroids AOA_amoACentroids.fna

# create map reads to OTUs
vsearch --usearch_global AOA_amoA_corr_nucl.fasta --db AOA_amoACentroids.fna --id 0.97 --log AOA_amoAOtus.log --otutabout AOA_amoAOtuTab.txt --sizein
