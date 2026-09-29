export BAM_FILE='possorted_bam.bam'
source $10X_PATH/sourceme.bash 
# Save the header lines
samtools view -H $BAM_FILE > SAM_header

# Filter alignments using filter.txt. Use LC_ALL=C to set C locale instead of UTF-8
samtools view $BAM_FILE | LC_ALL=C grep -F -f filter.txt > filtered_SAM_body

# Combine header and body
cat SAM_header filtered_SAM_body > filtered.sam

# Convert filtered.sam to BAM format
samtools view -b filtered.sam > filtered.bam

4942
4942