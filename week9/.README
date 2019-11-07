# To copy the ER4 ChIP peaks to this week's folder
cp ../week6/ER4_peaks_std.bed ./

# To choose the top 100 peaks and save as a new file
sort -r -n -k 5 ER4_peaks_std.bed | head -100 > ER4_100_peak.bed

# To copy the mouse chr19 sequence to this week's folder and pull out sequences of top 100 ER4 peaks (get_seqs.py was written to accomplish this)
cp ../week6/mus_musculus_chr19.fa ./
python get_seqs.py ER4_100_peak.bed mus_musculus_chr19.fa ER4_100_peak_seqs.fa

# To identify conserved motifs within the selected ER4 binding sites and compare to JASPAR database
meme-chip -meme-maxw 20 -db JASPAR_CORE_2016.meme -o meme-chip_stuff ER4_100_peak_seqs.fa

# To get the sequences of ALL ER4 binding sites and find motifs within these sequences (uses top motif id'ed by meme)
python get_seqs.py ER4_peaks_std.bed mus_musculus_chr19.fa ER4_peak_seqs.fa
fimo meme-chip_stuff/combined.meme ER4_peak_seqs.fa

# To generate motif location density plot from FIMO output (the .gff file)
python motif_location_density.py fimo_out/fimo.gff ER4_peaks_std.bed motif_density_over_binding_site.png