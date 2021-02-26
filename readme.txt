You need to download four supplemental files from GEO:
	1) scK27_4lanes_rc_at_25951_bipeaks.txt
	2) scK27_4lanes_rc_at_79100peaks_width10k.txt
	3) scK4_rc_at_25951_bipeaks.txt
	4) scH3K4me3_7798_cell_counts_at_WBCpeaks.txt
These files can be genrated by the R code "iscChiC_seq_generate_rc_files_for_sc.r"

files in  predata/:

	bivalent_domain25951_annot_gene.txt  ## 25951 bivalent domains with gene annotation
	bulk21_H3K27me3_at_bipeaks_counts.txt ## H3K27me3 ChIP-seq readcount at	25951 bivalent domains
	bulk72_H3K4me3_at_25951_bipeaks_counts.txt ## H3K4me3 ChIP-seq readcount at 25951 bivalent domains
	bulk_wbc_file_name.txt ## H3K4me3 ChIP-seq filename
	bulk_wbc_rc_52798_72_mat.txt ## H3K4me3 ChIP-seq readcount at 52798 H3K4me3 peaks
	combined_sicer_E100_peaks.sort.merge.txt ## 79100 H3K27me3 peaks
	Encode_wbcK27_rc_at_79100peaks_width10k.txt ## H3K27me3 ChIP-seq readcount at 79100 H3K27me3 peaks
	H3K4me3_peakname3_2.txt ## 52798 H3K4me3 peaks with filtering outrange peaks
	H3K4me3_peakname3.txt ## 52798 H3K4me3 peaks
	K4_K27_bivalent_domain.txt ##25951 bivalent domains
	scH3K4me3_7798_cell_counts_at_WBCpeaks.txt ##H3K4me3 single cell readcount at 52798 H3K4me3 peaks
	scK27_4lanes_rc_at_25951_bipeaks.txt ##H3K27me3 single cell readcount at 25951 bivalent domains
	scK27_4lanes_rc_at_79100peaks_width10k.txt ## H3K27me3 single cell readcount at 79100 H3K27me3 peaks
	scK4_rc_at_25951_bipeaks.txt ##H3K4me3 single cell readcount at 25951 bivalent domains
	wc_bulk21.bed.txt  ## H3K27me3 ChIP-seq depth
	wc_encode_uniq.txt  ## H3K4me3 ChIP-seq depth
	wc_K27_uniq.txt ## H3K27me3 single cell depth 
	wc_K4_uniq.txt ## H3K4me3 single cell depth
        barcode_96.txt ## 96 barcodes

script and matlab files:
	run script_get_bivalent to generate the file of K4_K27_bivalent_domain.txt
	Follow guidance in iscChIC_mapping.m abd run it for single cell mapping
	Run isChIC_seq_K4_and_K27_clustering_comparison.m for obtaining the plots of
		1) K4 single cell clustering and annotation
		2) K27 single cell clustering and annotation
		3) Matching between K4 and K27 clusters
		4) heatmap showing correlation between K4 CV and K27 CV

Acknowledgment:

AdvancedColormap.m is downloaded from https://www.mathworks.com/matlabcentral/fileexchange/41583-advancedcolormap created by Andriy Nych 
