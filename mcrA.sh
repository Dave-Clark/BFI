cd /media/dave/storage/macronutrient/BFI_paper/data/sequences/PKG_ENQ_723_data_transfer_RIVERS_FUNCTIONAL_GENES_NERCMACRO

mkdir ../mcrA

# target samples shared with mcrA first
SAMPS=$(ls *R1.fastq.gz | grep -v "_b_*")
for f in ${SAMPS}; do
cutadapt -g "GGTGGTGTMGGDTTCACMCARTA;min_overlap=21" -G "CGTTCATBGCGTAGTTVGGRTAGT;min_overlap=22" -j 0 --discard-untrimmed -o /media/dave/storage/macronutrient/BFI_paper/data/sequences/mcrA/${f%L001_R1.fastq.gz}R1_primer_trimmed.fastq -p /media/dave/storage/macronutrient/BFI_paper/data/sequences/mcrA/${f%L001_R1.fastq.gz}R2_primer_trimmed.fastq $f ${f%R1.fastq.gz}R2.fastq.gz
done

cd ../mcrA

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

# take all reads between > 415
for f in *.fna; do
awk '!/^>/ { next } { getline seq } length(seq) >= 415 { print $0 "\n" seq }' $f > ${f%.fna}_len_filtered.fna
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
paste(dat[V2 <= 4900, V1], collapse = ",")
q()

rm -f {AUG-AS1_3_GTAGAGGA-CTCTCTAT_len_filtered.fna,AUG-CW2_3_GCTACGCT-AGAGTAGA_len_filtered.fna,FEB-AS1_2_CGTACTAG-GCGTAAGA_len_filtered.fna,FEB-CA3_1_CTCTCTAC-CTCTCTAT_len_filtered.fna,FEB-GA2_2_AAGAGGCA-GCGTAAGA_len_filtered.fna}


# concatenate all length filtered fasta files into one
cat *_len_filtered.fna > allSequences.fna

# get total final lib sizes
grep -c "^>" allSequences.fna # 773580

# dereplicate reads by sample to preserve abundance info
~/bioinformatics/usearch11.0.667_i86linux32 -fastx_uniques_persample allSequences.fna -fastaout uniq_DNA_seqs.fna -minuniquesize 1 -sizeout

# use subsampled AOA db from PRINCe work
mkdir framebotOut

java -Xmx16g -jar ~/bioinformatics/RDPTools/FrameBot.jar framebot -o framebotOut/mcrA -N -i 0.5 -l 50 /media/dave/storage/PRINCE/data/sequenceData/fungeneDBs/mcrA/derep_mcrA.fa uniq_DNA_seqs.fna

cd framebot

awk '{ if ($0 !~ />/) {print toupper($0)} else {print $0} }' mcrA_corr_prot.fasta > uc_mcrA_corr_prot.fasta

# convert fs corrected aa sequences to upper case for usearch
~/bioinformatics/usearch11.0.667_i86linux32 -fastx_uniques uc_mcrA_corr_prot.fasta -fastaout AAVs.fasta -minuniquesize 1 -relabel AAV

# convert to amino acid variant OTU table
# might have to split into chunks and recombine afterwards
~/bioinformatics/usearch11.0.667_i86linux32 -otutab uc_mcrA_corr_prot.fasta -otus AAVs.fasta -id 1 -minsize 0 -otutabout mcrA_AAV.txt -qmask none -dbmask none

# use vsearch to dereplicate frameshift corrected nucleotide sequences
vsearch --derep_fulllength mcrA_corr_nucl.fasta --output derep_nucl_seqs.fna --minuniquesize 2 --relabel OTU --fasta_width 0 --sizein --sizeout

# sort by cluster size
vsearch --sortbysize derep_nucl_seqs.fna --output derepSorted.fna

# pick otu centroids
vsearch --cluster_size derepSorted.fna --id 0.97 --centroids mcrACentroids.fna

# create map reads to OTUs
vsearch --usearch_global mcrA_corr_nucl.fasta --db mcrACentroids.fna --id 0.97 --log mcrAOtus.log --otutabout mcrAOtuTab.txt --sizein
