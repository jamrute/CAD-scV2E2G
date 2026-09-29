parent_dir=$(pwd)

# Iterate over subfolders
for subfolder in "$parent_dir"/*/; do
  subfolder_name=$(basename "$subfolder")
  
  # Get the BAM files within subfolder
  bam_files="$subfolder"*.bam.bai
  
  # Rename each BAM file with the prefix
  for bam_file in $bam_files; do
    file_name=$(basename "$bam_file")
    new_bam_file="$subfolder${subfolder_name}_${file_name}"
    mv "$bam_file" "$new_bam_file"
  done
done