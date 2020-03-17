cd /media/dave/storage/macronutrient/BFI_paper/data/sequences/PKG_ENQ_723_data_transfer_RIVERS_FUNCTIONAL_GENES_NERCMACRO

mkdir ../AOB

# target samples shared with AOA first
SAMPS=$(ls *R1.fastq.gz | grep -v "_b_*")
for f in ${SAMPS}; do
cutadapt -g "GGGGTTTCTACTGGTGGT;min_overlap=17" -G "CCCCTCKGSAAAGCCTTCTTC;min_overlap=20" -j 0 --discard-untrimmed -o /media/dave/storage/macronutrient/BFI_paper/data/sequences/AOB/${f%L001_R1.fastq.gz}R1_primer_trimmed.fastq -p /media/dave/storage/macronutrient/BFI_paper/data/sequences/AOB/${f%L001_R1.fastq.gz}R2_primer_trimmed.fastq $f ${f%R1.fastq.gz}R2.fastq.gz
done

cd ../AOB

# make directories for merged sequence files, html output and discarded sequences
mkdir {mergedSeqs,discardedSeqs,htmlReports}

# try pair end alignment. seems to work for ~60% of seqs
for f in *R1_primer_trimmed.fastq; do
  fastp -i $f -I ${f%R1_primer_trimmed.fastq}R2_primer_trimmed.fastq -q 20 -u 20 -l 150 -y -A -5 -c --overlap_len_require 10 --merge --merged_out mergedSeqs/${f%R1_primer_trimmed.fastq}q_trimmed_merged.fastq -o discardSeqs/${f%primer_trimmed.fastq}discarded.fastq -O discardSeqs/${f%R1_primer_trimmed.fastq}R2_discarded.fastq -h htmlReports/${f%R1_primer_trimmed.fastq}quality_report.html
done

cd mergedSeqs

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
  seqtk seq -a $f > ${f%_q_trimmed_merged.00.0_0.cor.fastq}.fna
done

# get read lengths for one file
### NEED TO LOOK AT LENGTH FILTERS
# get read lengths for one file
R
library(data.table)
dat <- fread("grep '^[^>]' FEB-AS1_2_CGTACTAG-GCGTAAGA.fna | awk '{print length}'")#"
summary(dat$V1)
q()

# take all reads between > 448
for f in *.fna; do
awk '!/^>/ { next } { getline seq } length(seq) >= 448 { print $0 "\n" seq }' $f > ${f%.fna}_len_filtered.fna
done


for f in *len_filtered.fna; do
	samp=">$(echo $f | cut -d "_" -f 1-2 | sed "s/-/x/g; s/_/x/g")-"
	sed -i "s/>/$samp/" $f
done

# calculate final lib sizes
R
library(data.table)
dat <- fread("grep -c '^>' *_len_filtered.fna", sep = ":")
dat[order(V2, decreasing = F), ]
paste(dat[V2 <= 2600, V1], collapse = ",")
q()

rm -f {AUG-CE1_1_TCCTGAGC-AGAGTAGA_len_filtered.fna,FEB-CW2_2_GGACTCCT-CTCTCTAT_len_filtered.fna,FEB-GA2_2_AAGAGGCA-GCGTAAGA_len_filtered.fna}

# concatenate all length filtered fasta files into one
cat *_len_filtered.fna > allSequences.fna

# get total final lib sizes
grep -c "^>" allSequences.fna # 836868

# dereplicate reads by sample to preserve abundance info
~/bioinformatics/usearch11.0.667_i86linux32 -fastx_uniques_persample allSequences.fna -fastaout uniq_DNA_seqs.fna -minuniquesize 1 -sizeout

# use subsampled AOA db from PRINCe work
mkdir framebotOut

java -Xmx16g -jar ~/bioinformatics/RDPTools/FrameBot.jar framebot -o framebotOut/AOB_amoA -N -i 0.5 -l 50 /media/dave/storage/PRINCE/data/sequenceData/fungeneDBs/AOB_amoA/derep_AOB_amoA.fna uniq_DNA_seqs.fna

cd framebot

awk '{ if ($0 !~ />/) {print toupper($0)} else {print $0} }' AOB_amoA_corr_prot.fasta > uc_AOB_amoA_corr_prot.fasta

# convert fs corrected aa sequences to upper case for usearch
~/bioinformatics/usearch11.0.667_i86linux32 -fastx_uniques uc_AOB_amoA_corr_prot.fasta -fastaout AAVs.fasta -minuniquesize 1 -relabel AAV

# convert to amino acid variant OTU table
# might have to split into chunks and recombine afterwards
~/bioinformatics/usearch11.0.667_i86linux32 -otutab uc_AOB_amoA_corr_prot.fasta -otus AAVs.fasta -id 1 -minsize 0 -otutabout AOB_AAV.txt -qmask none -dbmask none

# use vsearch to dereplicate frameshift corrected nucleotide sequences
vsearch --derep_fulllength AOB_amoA_corr_nucl.fasta --output derep_nucl_seqs.fna --minuniquesize 2 --relabel OTU --fasta_width 0 --sizein --sizeout

# sort by cluster size
vsearch --sortbysize derep_nucl_seqs.fna --output derepSorted.fna

# pick otu centroids
vsearch --cluster_size derepSorted.fna --id 0.97 --centroids AOB_amoACentroids.fna

# create map reads to OTUs
vsearch --usearch_global AOB_amoA_corr_nucl.fasta --db AOB_amoACentroids.fna --id 0.97 --log AOB_amoAOtus.log --otutabout AOB_amoAOtuTab.txt --sizein
