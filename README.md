A. Download files
--------------------------------------

1. Download the files from github and change the directory to CRK. You need to login github, and download the testing123.git

<pre>
[wlku@matrix] mv testing123 CRK
[wlku@matrix] cd CRK
[wlku@matrix CRK] ls
[wlku@matrix CRK] AdvancedColormap.m  <b>data</b>  <b>Figures</b>  <b>GSE105012</b>  README.md  <b>src</b> violin.m  violinplot.m
</pre>

2. Download GSE105012_RAW.tar from GEO website and save it to the folder  <b>GSE105012</b>
<pre>
[wlku@matrix CRK] cd GSE105012
[wlku@matrix GSE105012] ls
[wlku@matrix GSE105012] GSE105012_RAW.tar
[wlku@matrix GSE105012] tar -xvf GSE105012_RAW.tar
[wlku@matrix GSE105012] gunzip *.gz
[wlku@matrix GSE105012] cd ..
</pre>
