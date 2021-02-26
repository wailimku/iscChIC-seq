%==================================================
%This code is to show an example for mapping the H3K4me3 single cell data
%You will need to set the path. For example, '/data/kuw/biocore/wlku/pipeline/ChIC_test/';
%Download all fastq files for one histone modification marks (For example:H3K4me3: GSM4147373-GSM4147852).
%run script_mkdir
%run script_mv1
%Do the mapping in each folder using Bowtie2
%For example "bowtie2  -p 18 -q -5 0 -3 0 -x /data/kuw/biocore/Basic_data/Bowtie2Index/hg18/genome -U SRR10386888_1.fastq |samtools view -bS - > SRR10386888_hg18.bam"
%==================================================
% In matlab, run the following code to genreate the file,"script_mapping", in each folder.
%***********************************
% Start (for generating "script_mapping")
%***********************************
path = '/data/kuw/biocore/wlku/pipeline/ChIC_test/'; %% need to change
outputname = '_readname.txt.R1.bed.hg18'; %% need to change

common_primer1 = 'AGAACCATGTCGTCAGTGTCCCCCCCCC';
common_primer2 = 'AGAACCATGTCGTCAGTGTCCCCCCCC';
common_primer3 = 'AGAACCATGTCGTCAGTGTCCCCCCC';
common_primer4 = 'AGAACCATGTCGTCAGTGTCCCCCC';
common_primer5 = 'AGAACCATGTCGTCAGTGTCCCCC';
common_primer6 = 'AGAACCATGTCGTCAGTGTCCCC';
common_primer7 = 'AGAACCATGTCGTCAGTGTCCC';
common_primer8 = 'AGAACCATGTCGTCAGTGTCC';
common_primer9 = 'AGAACCATGTCGTCAGTGTC';
common_primer10 ='AGAACCATGTCGTCAGTGT';

      
listing = dir(strcat(char(path),'/SRR*/SRR*_2.fastq'));
msize = max(size(listing));

filename  = cell([msize, 1]);
for i = 1:msize
    filename(i) = cellstr(listing(i).name);
end;

foldername = filename;
for i = 1:msize
    c = strsplit(char(filename(i)),'_2');
    foldername(i) = cellstr(c(1));
end;    

for i = 1:msize
    fp = fopen(strcat(char(path),char(foldername(i)),'/','script_mapping'),'w');
    fprintf(fp,'%s%s%s%s\n', 'sed -i -e ''s/',char(common_primer1),'/\t/g'' ',strcat(char(path),char(foldername(i)),'/',char(filename(i))));
    fprintf(fp,'%s%s%s%s\n', 'sed -i -e ''s/',char(common_primer2),'/\t/g'' ',strcat(char(path),char(foldername(i)),'/',char(filename(i))));
    fprintf(fp,'%s%s%s%s\n', 'sed -i -e ''s/',char(common_primer3),'/\t/g'' ',strcat(char(path),char(foldername(i)),'/',char(filename(i))));
    fprintf(fp,'%s%s%s%s\n', 'sed -i -e ''s/',char(common_primer4),'/\t/g'' ',strcat(char(path),char(foldername(i)),'/',char(filename(i))));
    fprintf(fp,'%s%s%s%s\n', 'sed -i -e ''s/',char(common_primer5),'/\t/g'' ',strcat(char(path),char(foldername(i)),'/',char(filename(i))));
    fprintf(fp,'%s%s%s%s\n', 'sed -i -e ''s/',char(common_primer6),'/\t/g'' ',strcat(char(path),char(foldername(i)),'/',char(filename(i))));
    fprintf(fp,'%s%s%s%s\n', 'sed -i -e ''s/',char(common_primer7),'/\t/g'' ',strcat(char(path),char(foldername(i)),'/',char(filename(i))));
    fprintf(fp,'%s%s%s%s\n', 'sed -i -e ''s/',char(common_primer8),'/\t/g'' ',strcat(char(path),char(foldername(i)),'/',char(filename(i))));
    fprintf(fp,'%s%s%s%s\n', 'sed -i -e ''s/',char(common_primer9),'/\t/g'' ',strcat(char(path),char(foldername(i)),'/',char(filename(i))));
    fprintf(fp,'%s%s%s%s\n', 'sed -i -e ''s/',char(common_primer10),'/\t/g'' ',strcat(char(path),char(foldername(i)),'/',char(filename(i))));
    fprintf(fp,'\n');
    fclose(fp);
end;    

%+========================================================================
%+========================================================================
command2 = '|awk ''{{if(NR%4==1) {key1=$1;getline;} if(NR%4==2) {if(NF==2){key0=$2; key2=$1;} k=NF;getline;};if (NR%4==3) {key3=$1;};getline} if(k==2) print key1"\t"key2"\t"key0}''>';
for i = 1:msize
    fp = fopen(strcat(char(path),char(foldername(i)),'/','script_mapping'),'a');

    infile = strcat(char(path),char(foldername(i)),'/',filename(i));
    outfile = strcat(char(path),char(foldername(i)),'/barcode_readname.txt');
    fprintf(fp,'%s\t%s%s%s\n','less',char(infile),char(command2),char(outfile));
    fprintf(fp,'\n');
    fclose(fp);
end;    
 
%+========================================================================
%+========================================================================
barcode_com = cell([96,1]);

a = readtable(strcat(char(path),'/predata/barcode_96.txt'),'ReadVariableNames',0);
barcode = table2array(a);

barcode_com2 = barcode_com;
for i = 1:96
    barcode_com(i) = cellstr('barcode_readname.txt');
    barcode_com2(i) = cellstr(strcat(barcode(i),'_','barcode_readname.txt'));
end;    


for i = 1:msize
    fp = fopen(strcat(char(path),char(foldername(i)),'/','script_mapping'),'a');
    for j = 1:96
        filein = strcat(char(path),char(foldername(i)),'/',char(barcode_com(j)));
        fileout = strcat(char(path),char(foldername(i)),'/',char(barcode_com2(j)));
        fprintf(fp,'%s %s %s%s%s\n','grep',char(barcode(j)),char(filein),'>',char(fileout));
    end;
    fprintf(fp,'\n');
    fclose(fp);
end;    

commands1 = '|awk ''{if(substr($2,1,8)==';
commands2 = ') print}''';
for i = 1:msize
    fp = fopen(strcat(char(path),char(foldername(i)),'/','script_mapping'),'a');
    for j = 1:96
        filein = strcat(char(path),char(foldername(i)),'/',char(barcode_com2(j)));
        fileout = strcat(char(path),char(foldername(i)),'/',char(barcode_com2(j)),'.v2.txt');
        fprintf(fp,'%s %s %s%s%s%s%s%s%s\n','less',char(filein),char(commands1),'"',char(barcode(j)),'"',char(commands2),'>',char(fileout));
    end;
    fprintf(fp,'\n');
    fclose(fp);
end;    
    

for i = 1:msize
    fp = fopen(strcat(char(path),char(foldername(i)),'/','script_mapping'),'a');
    for j = 1:96
        filein = strcat(char(path),char(foldername(i)),'/',char(barcode_com2(j)),'.v2.txt');
        fprintf(fp,'%s %s\n','sed -i -e ''s/@//g''',char(filein));
    end;
    fprintf(fp,'\n');
    fclose(fp);
end;    


command3 = '|awk ''{print $1}''>';
for i = 1:msize
    fp = fopen(strcat(char(path),char(foldername(i)),'/','script_mapping'),'a');
    for j = 1:96
        filein = strcat(char(path),char(foldername(i)),'/',char(barcode_com2(j)),'.v2.txt');
        fileout = strcat(char(path),char(foldername(i)),'/',char(barcode(j)),'_readname.txt');
        
        fprintf(fp,'%s %s%s%s\n','less',char(filein),char(command3),char(fileout));
    end;
    fprintf(fp,'\n');
    fclose(fp);
end;    


command4 = 'bedtools bamtobed -i';
command4b = '|awk ''{if($5>=10) print}''';
for i = 1:msize
    a = dir(strcat(char(path),char(foldername(i)),'/','*hg18.bam'));
    bedin = strcat(char(path),char(foldername(i)),'/',cellstr(a.name));
    bedout = strcat(char(path),char(foldername(i)),'/','test_result_hg18.bed');
    fp = fopen(strcat(char(path),char(foldername(i)),'/','script_mapping'),'a');
    fprintf(fp,'%s %s%s%s%s\n',char(command4),char(bedin),char(command4b),'>',char(bedout));
    fprintf(fp,'\n');
    fclose(fp);
end;


command5 = 'grep -Fwf';
for i = 1:msize
    bedin = strcat(char(path),char(foldername(i)),'/','test_result_hg18.bed');
    fp = fopen(strcat(char(path),char(foldername(i)),'/','script_mapping'),'a');
    for j = 1:96
        filein = strcat(char(path),char(foldername(i)),'/',char(barcode(j)),'_readname.txt');
        fileout = strcat(char(path),char(foldername(i)),'/',char(barcode(j)),char(outputname)); 
        fprintf(fp,'%s %s %s%s%s\n', char(command5),char(filein),char(bedin),'>',char(fileout));
    end;
    fprintf(fp,'\n');
    fclose(fp);
end;

%***********************************
% End(for generating "script_mapping")
%***********************************
%=========================================
%run "sh script_mapping" in each folder
%==================================================

command6 = 'ls ./*/*';
command7 = '|awk ''{print "sort -u -k1,1 -k2,2 -k3,3 "$1" > "$1".uniq" }''>script_sort_uniq';
command8 = strcat(char(command6),char(outputname),char(command7));
fp = fopen(strcat(char(path),'script_generate_sort_script'),'w');
fprintf(fp,'%s\n', char(command8));
fclose(fp);


%============================================
%run "sh script_generate_sort_script" in path
%run "sh script_sort_uniq"
%==================================================
command9 = 'wc -l ./*/*.uniq|awk ''{print $1"\t"$2}''>wc_uniq2.txt';
command10 = 'sed -e ''$ d''  wc_uniq2.txt>wc_uniq.txt';
command11 = 'rm wc_uniq2.txt';
fp = fopen(strcat(char(path),'script_get_readnumber'),'w');
fprintf(fp,'%s\n', char(command9));
fprintf(fp,'%s\n', char(command10));
fprintf(fp,'%s\n', char(command11));

fclose(fp);
%============================================
%run "sh script_get_readnumber" in path
%==================================================
nw = 480; %number of well, need to change;
kcheck = zeros(96,nw);
a = readtable(strcat(char(path),'wc_uniq.txt'),'ReadVariableNames',0);
b = table2array(a(:,2));
d = table2array(a(:,1));
wc_size = reshape(table2array(a(:,1)),[96,nw]);
check = zeros(max(size(d)),1);
for i =1:nw
    [ia,ib] = sort(wc_size(:,i),'descend');
    ic = ib(1:30); % top 30 cells are considered
    id = ia(1:30); % top 30 cells are considered
    q = find(id<3000); %read cut off=3000
    ic(q)=[];
    if(min(size(ic))>0)
        m = (i-1)*96+ic;
        check(m) = 1;
        kcheck(m,i) = 1;
    end;
end;
q = find(check==1);

sel_file = b(q);
T1 = table(sel_file);

% T1 is the files for selected cells.
