cd /media/dave/storage/macronutrient/BFI_paper/data/sequences

# cut primer sequences
# remove any sequences that didn't contain both the forward and reverse primer
# target samples shared with AOA first
SAMPS=$(ls *R1* | grep -v "_b_")

for f in *_b_*R1.fastq.gz; do
cutadapt -g "GGMATGGTKCCSTGGCA ;min_overlap=16" -G "GCCTCGATCAGRTTRTGG;min_overlap=18" -j 0 --discard-untrimmed -o /media/dave/storage/macronutrient/BFI_paper/data/sequences/nirK/${f%L001_R1.fastq.gz}primer_trimmed.fastq -p /media/dave/storage/macronutrient/BFI_paper/data/sequences/nirK/${f%L001_R1.fastq.gz}2_primer_trimmed.fastq $f ${f%L001_R1.fastq.gz}L001_R2.fastq.gz
done

# not returning any seqs
