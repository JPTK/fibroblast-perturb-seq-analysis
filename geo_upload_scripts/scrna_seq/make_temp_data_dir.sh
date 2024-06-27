mkdir /data/torsten/lara-haematopoesis-mouse/tmp_space/scrna_seq
# mkdir /data/scrna_seq


for link in /data/torsten/lara-haematopoesis-mouse/geo_upload_space/scrna_seq/*; do
    if [ -L "$link" ]; then # Check if it's a symbolic link
        # Copy the file it points to into temp_dir with the link name
        # cp --dereference "$link" "/scratch/scrna_seq/$(basename "$link")"
        cp --dereference "$link" "/data/torsten/lara-haematopoesis-mouse/tmp_space/scrna_seq/$(basename "$link")"
    fi
done