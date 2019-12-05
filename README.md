A. Download files
--------------------------------------

1. Download the files from github and change the directory to CRK. You need to login github, and download the testing123.git

<pre>
[wlku@matrix] mv testing123 CRK
[wlku@matrix] cd CRK
[wlku@matrix CRK] ls
[wlku@matrix CRK] AdvancedColormap.m  <b>data</b>  <b>Figures</b>  <b>GSE139857</b>  README.md  <b>src</b> violin.m  violinplot.m
</pre>

2. Download GSE139857_RAW.tar from GEO website and save it to the folder  <b>GSE139857</b>
<pre>
[wlku@matrix CRK] cd GSE139857
[wlku@matrix GSE105012] ls
[wlku@matrix GSE105012] GSE139857_RAW.tar
[wlku@matrix GSE105012] tar -xvf GSE139857_RAW.tar
[wlku@matrix GSE105012] gunzip *.gz
[wlku@matrix GSE105012] cd ..
</pre>



B. Demultiplexing and mapping single cell H3K4me3 data
--------------------------------------

1. 
<pre>
[wlku@matrix CRK] sh ./src/figure1_code/script_sicer_low_cell
</pre>


C. Quality analysis and clustering for single cell H3K4me3 data
--------------------------------------

D. Demultiplexing and mapping single cell H3K27me3 data
--------------------------------------

E. Quality analysis and clustering for single cell H3K27me3 data
--------------------------------------

F. Matching clusters between single cell H3K4me3 and H3K27me3 data
--------------------------------------

G. Computation for spatial and density cell-to-cell variation
--------------------------------------


--------------------------------------
