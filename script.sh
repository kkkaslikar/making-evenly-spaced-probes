# echo "running blockParse.py"

# blockParse.py --file Tgfbr1.fa --minLength 19 --maxLength 25 --min_Tm 55 --max_Tm 60 --salt 50 --formamide 0

# echo "running bowtie2"

# bowtie2 -x /mnt/f/Bioinformatics_Data/Indexes/Mus_musculus_GRCm38_unmasked_bt2_index/MM.GRCm38.unmasked -U Tgfbr1.fastq --no-hd -t -k 100 --very-sensitive-local -S Tgfbr1_probes.sam


# echo "running outputClean.py"

# outputClean.py --file Tgfbr1_probes.sam --unique --salt 50 --formamide 0 --output Tgfbr1_cleaned_probes


echo "running structureCheck.py"

structureCheck.py --file Tgfbr1_cleaned_probes_chr.bed --formamide 0 --salt 50 --hybTemp 51


