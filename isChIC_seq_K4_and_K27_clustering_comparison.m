

path0 = '/data/kuw/biocore/wlku/pipeline/iscChIC_seq_test'; %% you will need to change this path
path = strcat(char(path),'/predata');

a= readtable(strcat(char(path),'/bulk21_H3K27me3_at_bipeaks_counts.txt'));
bulk_rc = table2array(a(:,1:21));

a = readtable(strcat(char(path),'/wc_bulk21.bed.txt'));
bulk_depth = table2array(a(:,1));
bulkfile = table2array(a(:,2));

ss = sum(bulk_rc);
bulk_cpm = bulk_rc;
for i = 1:21
    bulk_cpm(:,i) = bulk_rc(:,i)*1000000./sum(bulk_rc(:,i));
end;

wbc_type_name={'B_cell_ENCFF529PGJ_noDup.bed                                                 '
    'B_cell_ENCFF705DLT_noDup.bed                                                 '
    'CD14_positive_monocyte_ENCFF001FYR_noDup.bed                                 '
    'CD14_positive_monocyte_ENCFF127BPU_noDup.bed                                 '
    'CD4_positive_alpha_beta_memory_T_cell_ENCFF576BBU_noDup.bed                  '
    'CD4_positive_alpha_beta_memory_T_cell_ENCFF991QCC_noDup.bed                  '
    'CD4_positive_alpha_beta_T_cell_ENCFF691MOI_noDup.bed                         '
    'CD4_positive_alpha_beta_T_cell_ENCFF815XHM_noDup.bed                         '
    'CD4_positive_CD25_positive_alpha_beta_regulatory_T_cell_ENCFF010QBJ_noDup.bed'
    'CD4_positive_CD25_positive_alpha_beta_regulatory_T_cell_ENCFF549FJX_noDup.bed'
    'CD8_positive_alpha_beta_T_cell_ENCFF784OKP_noDup.bed                         '
    'CD8_positive_alpha_beta_T_cell_ENCFF921LVS_noDup.bed                         '
    'effector_memory_CD4_positive_alpha_beta_T_cell_ENCFF085HVZ_noDup.bed         '
    'effector_memory_CD4_positive_alpha_beta_T_cell_ENCFF173ACS_noDup.bed         '
    'naive_thymus_derived_CD4_positive_alpha_beta_T_cell_ENCFF196VTT_noDup.bed    '
    'naive_thymus_derived_CD4_positive_alpha_beta_T_cell_ENCFF509MTV_noDup.bed    '
    'Natural_killer_H3K27me3_SRR609646_7_0_0_mapq10_noDup.bed                     '
    'NKcells_H3K27me3_SRR5006015_0_0_mapq10_noDup.bed'                             
    'NKcells_H3K27me3_SRR5006016_0_0_mapq10_noDup.bed'                             
    'T_cell_ENCFF503UGV_noDup.bed                                                 '
    'T_cell_ENCFF687JRX_noDup.bed                                                 '};
wbc_cpm = log2(bulk_cpm+1);
import bioma.data.DataMatrix     
DM_B = DataMatrix(wbc_cpm(:,1:2));
DM_Mono = DataMatrix(wbc_cpm(:,3:4))
DM_T = DataMatrix(wbc_cpm(:,7:8));
DM_NK = DataMatrix(wbc_cpm(:,17:19));

wbc_cpm_type=[mean(DM_B')',mean(DM_Mono')',mean(DM_T')',mean(DM_NK')'];
A27wbc_type = wbc_cpm_type;
%=========================================================
a= readtable(strcat(char(path),'/bulk72_H3K4me3_at_25951_bipeaks_counts.txt'));
bulk_rc = table2array(a(:,1:72));
bulk_cpm = bulk_rc;
for i = 1:72
    bulk_cpm(:,i) = bulk_rc(:,i)*1000000/sum(bulk_rc(:,i));
end;
a= readtable(strcat(char(path),'/bulk_wbc_file_name.txt'));

bulkfile = table2array(a(:,2));
wbc_type_name={'B_cell',
'CD14-positive_monocyte',
'CD4-positive_CD25-positive_alpha-beta_regulatory_T_cell',
'CD4-positive_alpha-beta_memory_T_cell',
'CD4-positive_helper_T_cell',
'CD8-positive_alpha-beta_T_cell',
'T-cell',
'common_myeloid_progenitor_CD34-positive',
'mononuclear_cell',
'naiveT',
'natural_killer_cell',
'neutrophil'};
wbc_cpm=log2(bulk_cpm+1);  
bulk_cpm = wbc_cpm;


import bioma.data.DataMatrix     
DM_B = DataMatrix(bulk_cpm(:,3:7));
DM_Mono = DataMatrix(bulk_cpm(:,15:16));
DM_T = DataMatrix(bulk_cpm(:,71:72));
DM_NK = DataMatrix(bulk_cpm(:,61:62));

wbc_cpm_type=[mean(DM_B')',mean(DM_Mono')',mean(DM_T')',mean(DM_NK')'];
A4wbc_type = wbc_cpm_type;

%=====================
a= readtable(strcat(char(path),'/Encode_wbcK27_rc_at_79100peaks_width10k.txt'));
bulk_rc = table2array(a(:,1:21));

a = readtable(strcat(char(path),'/wc_bulk21.bed.txt'));

bulk_depth = table2array(a(:,1));
bulkfile = table2array(a(:,2));

ss = sum(bulk_rc);
bulk_cpm = bulk_rc;
for i = 1:21
    bulk_cpm(:,i) = bulk_rc(:,i)*1000000./sum(bulk_rc(:,i));
end;

wbc_type_name={'B_cell_ENCFF529PGJ_noDup.bed                                                 '
    'B_cell_ENCFF705DLT_noDup.bed                                                 '
    'CD14_positive_monocyte_ENCFF001FYR_noDup.bed                                 '
    'CD14_positive_monocyte_ENCFF127BPU_noDup.bed                                 '
    'CD4_positive_alpha_beta_memory_T_cell_ENCFF576BBU_noDup.bed                  '
    'CD4_positive_alpha_beta_memory_T_cell_ENCFF991QCC_noDup.bed                  '
    'CD4_positive_alpha_beta_T_cell_ENCFF691MOI_noDup.bed                         '
    'CD4_positive_alpha_beta_T_cell_ENCFF815XHM_noDup.bed                         '
    'CD4_positive_CD25_positive_alpha_beta_regulatory_T_cell_ENCFF010QBJ_noDup.bed'
    'CD4_positive_CD25_positive_alpha_beta_regulatory_T_cell_ENCFF549FJX_noDup.bed'
    'CD8_positive_alpha_beta_T_cell_ENCFF784OKP_noDup.bed                         '
    'CD8_positive_alpha_beta_T_cell_ENCFF921LVS_noDup.bed                         '
    'effector_memory_CD4_positive_alpha_beta_T_cell_ENCFF085HVZ_noDup.bed         '
    'effector_memory_CD4_positive_alpha_beta_T_cell_ENCFF173ACS_noDup.bed         '
    'naive_thymus_derived_CD4_positive_alpha_beta_T_cell_ENCFF196VTT_noDup.bed    '
    'naive_thymus_derived_CD4_positive_alpha_beta_T_cell_ENCFF509MTV_noDup.bed    '
    'Natural_killer_H3K27me3_SRR609646_7_0_0_mapq10_noDup.bed                     '
    'NKcells_H3K27me3_SRR5006015_0_0_mapq10_noDup.bed'                             
    'NKcells_H3K27me3_SRR5006016_0_0_mapq10_noDup.bed'                             
    'T_cell_ENCFF503UGV_noDup.bed                                                 '
    'T_cell_ENCFF687JRX_noDup.bed                                                 '};

wbc_cpm = log2(bulk_cpm+1);

import bioma.data.DataMatrix     
DM_B = DataMatrix(wbc_cpm(:,1:2));
DM_Mono = DataMatrix(wbc_cpm(:,3:4));
DM_T = DataMatrix(wbc_cpm(:,7:8));
DM_NK = DataMatrix(wbc_cpm(:,17:19));

wbc_cpm_type=[mean(DM_B')',mean(DM_Mono')',mean(DM_T')',mean(DM_NK')'];

[PValues_T_B, TScores1] = mattest(DM_T, DM_B);
[PValues_T_NK, TScores1] = mattest(DM_T, DM_NK);
[PValues_T_Mono, TScores1] = mattest(DM_T, DM_Mono);

[PValues_B_T, TScores1] = mattest(DM_B, DM_T);
[PValues_B_NK, TScores1] = mattest(DM_B, DM_NK);
[PValues_B_Mono, TScores1] = mattest(DM_B, DM_Mono);

[PValues_NK_B, TScores1] = mattest(DM_NK, DM_B);
[PValues_NK_T, TScores1] = mattest(DM_NK, DM_T);
[PValues_NK_Mono, TScores1] = mattest(DM_NK, DM_Mono);

[PValues_Mono_B, TScores1] = mattest(DM_Mono, DM_B);
[PValues_Mono_T, TScores1] = mattest(DM_Mono, DM_T);
[PValues_Mono_NK, TScores1] = mattest(DM_Mono, DM_NK);


pcut=0.05;
fcut = 0.4;
fcut2 = 0.4;

ncut = log2(1+1);
n2cut = log2(100000000+1);
qp_B = find(PValues_B_T<pcut & PValues_B_Mono<pcut & PValues_B_NK<pcut & wbc_cpm_type(:,1)-wbc_cpm_type(:,2)>fcut & wbc_cpm_type(:,1)-wbc_cpm_type(:,3)>fcut & wbc_cpm_type(:,1)-wbc_cpm_type(:,4)>fcut & wbc_cpm_type(:,1)>ncut &  wbc_cpm_type(:,2)<n2cut &  wbc_cpm_type(:,3)<n2cut &  wbc_cpm_type(:,4)<n2cut);
qp_Mono = find(PValues_Mono_T<pcut & PValues_Mono_B<pcut & PValues_Mono_NK<pcut & wbc_cpm_type(:,2)-wbc_cpm_type(:,1)>fcut2 & wbc_cpm_type(:,2)-wbc_cpm_type(:,3)>fcut2 & wbc_cpm_type(:,2)-wbc_cpm_type(:,4)>fcut2 & wbc_cpm_type(:,2)>ncut &  wbc_cpm_type(:,1)<n2cut &  wbc_cpm_type(:,3)<n2cut &  wbc_cpm_type(:,4)<n2cut );
qp_T = find(PValues_T_B<pcut & PValues_T_Mono<pcut & PValues_T_NK<pcut & wbc_cpm_type(:,3)-wbc_cpm_type(:,1)>fcut & wbc_cpm_type(:,3)-wbc_cpm_type(:,2)>fcut & wbc_cpm_type(:,3)-wbc_cpm_type(:,4)>fcut & wbc_cpm_type(:,3)>ncut &  wbc_cpm_type(:,2)<n2cut &  wbc_cpm_type(:,1)<n2cut &  wbc_cpm_type(:,4)<n2cut );
qp_NK = find(PValues_NK_T<pcut & PValues_NK_Mono<pcut & PValues_NK_B<pcut & wbc_cpm_type(:,4)-wbc_cpm_type(:,2)>fcut & wbc_cpm_type(:,4)-wbc_cpm_type(:,3)>fcut & wbc_cpm_type(:,4)-wbc_cpm_type(:,1)>fcut & wbc_cpm_type(:,4)>ncut &  wbc_cpm_type(:,2)<n2cut &  wbc_cpm_type(:,3)<n2cut &  wbc_cpm_type(:,1)<n2cut );

%=========================================================

a= readtable(strcat(char(path),'/wc_encode_uniq.txt'),'delimiter','\t');
bulk_depth = table2array(a(:,1));
a= readtable(strcat(char(path),'/bulk_wbc_rc_52798_72_mat.txt'));

bulk_peakname= table2array(a(:,1));
bulk_rc = table2array(a(:,2:73));

bulk_cpm = bulk_rc;
for i = 1:72
    bulk_cpm(:,i) = bulk_rc(:,i)*1000000/bulk_depth(i);
end;
a= readtable(strcat(char(path),'/bulk_wbc_file_name.txt'));
bulkfile = table2array(a(:,2));
wbc_type_name={'B_cell',
'CD14-positive_monocyte',
'CD4-positive_CD25-positive_alpha-beta_regulatory_T_cell',
'CD4-positive_alpha-beta_memory_T_cell',
'CD4-positive_helper_T_cell',
'CD8-positive_alpha-beta_T_cell',
'T-cell',
'common_myeloid_progenitor_CD34-positive',
'mononuclear_cell',
'naiveT',
'natural_killer_cell',
'neutrophil'};
wbc_cpm=log2(bulk_cpm+1);  
bulk_cpm = wbc_cpm;

import bioma.data.DataMatrix     
DM_B = DataMatrix(bulk_cpm(:,3:7));
DM_Mono = DataMatrix(bulk_cpm(:,15:16));
DM_T = DataMatrix(bulk_cpm(:,71:72));
DM_NK = DataMatrix(bulk_cpm(:,61:62));

wbc_cpm_type=[mean(DM_B')',mean(DM_Mono')',mean(DM_T')',mean(DM_NK')'];

[PValues_T_B, TScores1] = mattest(DM_T, DM_B);
[PValues_T_NK, TScores1] = mattest(DM_T, DM_NK);
[PValues_T_Mono, TScores1] = mattest(DM_T, DM_Mono);

[PValues_B_T, TScores1] = mattest(DM_B, DM_T);
[PValues_B_NK, TScores1] = mattest(DM_B, DM_NK);
[PValues_B_Mono, TScores1] = mattest(DM_B, DM_Mono);

[PValues_NK_B, TScores1] = mattest(DM_NK, DM_B);
[PValues_NK_T, TScores1] = mattest(DM_NK, DM_T);
[PValues_NK_Mono, TScores1] = mattest(DM_NK, DM_Mono);

[PValues_Mono_B, TScores1] = mattest(DM_Mono, DM_B);
[PValues_Mono_T, TScores1] = mattest(DM_Mono, DM_T);
[PValues_Mono_NK, TScores1] = mattest(DM_Mono, DM_NK);


pcut=0.05;
fcut = 0.2;
fcut2 = 0.2;

ncut = log2(1+1);
n2cut = log2(100000000+1);
q_B = find(PValues_B_T<pcut & PValues_B_Mono<pcut & PValues_B_NK<pcut & wbc_cpm_type(:,1)-wbc_cpm_type(:,2)>fcut & wbc_cpm_type(:,1)-wbc_cpm_type(:,3)>fcut & wbc_cpm_type(:,1)-wbc_cpm_type(:,4)>fcut & wbc_cpm_type(:,1)>ncut &  wbc_cpm_type(:,2)<n2cut &  wbc_cpm_type(:,3)<n2cut &  wbc_cpm_type(:,4)<n2cut);
q_Mono = find(PValues_Mono_T<pcut & PValues_Mono_B<pcut & PValues_Mono_NK<pcut & wbc_cpm_type(:,2)-wbc_cpm_type(:,1)>fcut2 & wbc_cpm_type(:,2)-wbc_cpm_type(:,3)>fcut2 & wbc_cpm_type(:,2)-wbc_cpm_type(:,4)>fcut2 & wbc_cpm_type(:,2)>ncut &  wbc_cpm_type(:,1)<n2cut &  wbc_cpm_type(:,3)<n2cut &  wbc_cpm_type(:,4)<n2cut );
q_T = find(PValues_T_B<pcut & PValues_T_Mono<pcut & PValues_T_NK<pcut & wbc_cpm_type(:,3)-wbc_cpm_type(:,1)>fcut & wbc_cpm_type(:,3)-wbc_cpm_type(:,2)>fcut & wbc_cpm_type(:,3)-wbc_cpm_type(:,4)>fcut & wbc_cpm_type(:,3)>ncut &  wbc_cpm_type(:,2)<n2cut &  wbc_cpm_type(:,1)<n2cut &  wbc_cpm_type(:,4)<n2cut );
q_NK = find(PValues_NK_T<pcut & PValues_NK_Mono<pcut & PValues_NK_B<pcut & wbc_cpm_type(:,4)-wbc_cpm_type(:,2)>fcut & wbc_cpm_type(:,4)-wbc_cpm_type(:,3)>fcut & wbc_cpm_type(:,4)-wbc_cpm_type(:,1)>fcut & wbc_cpm_type(:,4)>ncut &  wbc_cpm_type(:,2)<n2cut &  wbc_cpm_type(:,3)<n2cut &  wbc_cpm_type(:,1)<n2cut );
%===================================================================
a= readtable(strcat(char(path),'/scK27_4lanes_rc_at_25951_bipeaks.txt'),'ReadVariableNames',0);
K27sc_count = table2array(a(:,1:9207));
%save(char(strcat(char(path),'/scK27_4lanes_rc_at_25951_bipeaks.mat')),'K27sc_count','-v7.3');
%a =load(strcat(char(path),'/scK27_4lanes_rc_at_25951_bipeaks.mat'));
%K27sc_count = a.K27sc_count;
a= readtable(strcat(char(path),'/scK27_4lanes_rc_at_79100peaks_width10k.txt'),'ReadVariableNames',0);
sc_count = table2array(a(:,1:9207));
%save(char(strcat(char(path),'/scK27_4lanes_rc_at_79100peaks_width10k.mat')),'sc_count','-v7.3');
%a =load(strcat(char(path),'/scK27_4lanes_rc_at_79100peaks_width10k.mat'));
%sc_count = a.sc_count;
a = readtable(strcat(char(path),'/wc_K27_uniq.txt'),'ReadVariableNames',0);
depth2 = table2array(a(:,1));
filename = table2array(a(:,2));
sc_cpm = sc_count;
K27sc_cpm = K27sc_count;
for i = 1:9207
    sc_cpm(:,i) = sc_count(:,i)*1000000./(depth2(i));
    K27sc_cpm(:,i) = K27sc_count(:,i)*1000000./(depth2(i));
end;
yy = sc_count./depth2';
qdel = find(sum(yy)<=0.15 |sum(sc_count)<3000);
sc_count(:,qdel)=[];
K27sc_count(:,qdel)=[];
sc_cpm(:,qdel)=[];
K27sc_cpm(:,qdel)=[];
yy(:,qdel)=[];
depth2(qdel)=[];
filename(qdel)=[];
sclog_cpm = log2(sc_cpm+1);
K27sclog_cpm = log2(K27sc_cpm+1);
B = sclog_cpm(:,:);
A = B;
q = find(A<0);
A(q)=0;
q = find(A>0);
A(q)=1;
qdel = find(sum(A')<10 );
B(qdel,:)=[];
cc = corr(B);
dd =1-cc;
lcc1 = exp(-dd/max(max(dd)));
for i = 1:max(size(cc))
    lcc1(i,i) = 1;
end;    
lcc2 = diag(1./((sum(lcc1)-1).^0.5));
lcc3 = diag(ones(max(size(cc)),1))-lcc2*lcc1*lcc2;
[v1,w1] = eigs(lcc3,max(size(cc)));
rng('default')
tic
ee = randn(max(size(cc)), max(size(cc)));
ee(:)=0;
m=1;
for i = 1:10
    for k = 1:10
        [ia,ib] = kmeans(v1(:,max(size(cc))-k:max(size(cc))),6);
        for j = 1:6
            q = find(ia==j);
            ee(q,q) = ee(q,q) +1;      
        end;
        m = m+1;
        
    end;
    toc
end;

rng('default')
nclus = 6;
[k1,k2] = kmeans(ee,nclus);

clus_cpmlog={};
for kkkk= 1:10
    qmin = [];
    for i=1:nclus
        q = find(k1==i);
        qmin = [qmin,max(size(q))];
    end;
    npeak2 = size(sc_count,1);
    clus_count = zeros(npeak2,nclus);
    clus_depth = zeros(nclus,1);
    clus_cpm = zeros(npeak2,1);
    for i=1:nclus
        q = find(k1==i);
        rr = randperm(max(size(q)),floor(max(size(q))/3));
        q= q(rr);    
        clus_count(:,i) = sum(sc_count(:,q)')';
        clus_depth(i) = sum(sum(sc_count(:,q)'));
        clus_cpm(:,i)= clus_count(:,i)*1000000./(clus_depth(i));
    end;
    clus_cpmlog{kkkk} = log2(clus_cpm+1);
end;

import bioma.data.DataMatrix
Y{:,1}  = DataMatrix([clus_cpmlog{1}(:,1),clus_cpmlog{2}(:,1),clus_cpmlog{3}(:,1),clus_cpmlog{4}(:,1),clus_cpmlog{5}(:,1),clus_cpmlog{6}(:,1)]);%%40WT
Y{:,2}  = DataMatrix([clus_cpmlog{1}(:,2),clus_cpmlog{2}(:,2),clus_cpmlog{3}(:,2),clus_cpmlog{4}(:,2),clus_cpmlog{5}(:,2),clus_cpmlog{6}(:,2)]);%%40WT
Y{:,3}  = DataMatrix([clus_cpmlog{1}(:,3),clus_cpmlog{2}(:,3),clus_cpmlog{3}(:,3),clus_cpmlog{4}(:,3),clus_cpmlog{5}(:,3),clus_cpmlog{6}(:,3)]);%%40WT
Y{:,4}  = DataMatrix([clus_cpmlog{1}(:,4),clus_cpmlog{2}(:,4),clus_cpmlog{3}(:,4),clus_cpmlog{4}(:,4),clus_cpmlog{5}(:,4),clus_cpmlog{6}(:,4)]);%%40WT
Y{:,5}  = DataMatrix([clus_cpmlog{1}(:,5),clus_cpmlog{2}(:,5),clus_cpmlog{3}(:,5),clus_cpmlog{4}(:,5),clus_cpmlog{5}(:,5),clus_cpmlog{6}(:,5)]);%%40WT
Y{:,6}  = DataMatrix([clus_cpmlog{1}(:,6),clus_cpmlog{2}(:,6),clus_cpmlog{3}(:,6),clus_cpmlog{4}(:,6),clus_cpmlog{5}(:,6),clus_cpmlog{6}(:,6)]);%%40WT

K27_ccorr = zeros(nclus,4);
tic
for i= 1:nclus
    xq = setdiff(1:1:nclus,i);
    PP= ones(79100,nclus-1);
    qt = 1:1:79100;
    for j= 1:nclus-1
       [PValues, TScores] = mattest(Y{:,i}, Y{:,xq(j)});
       PP(:,j) = double(PValues);
       q = find(PP(:,j)<0.05 & mean(Y{:,i}')'-mean(Y{:,xq(j)}')'>0.4);
       qt = intersect(qt,q);
    end;    
    for k = 1:4
        if(k==1)
            qxx = [qp_B'];
        elseif(k==2) 
            qxx = [qp_Mono'];
        elseif(k==3) 
            qxx = [qp_T'];
        elseif(k==4) 
            qxx = [qp_NK'];
        end;  
        q1 = qt;
        q = intersect(q1,qxx);
        a = real(-log2((1-sum(hygepdf(0:max(size(q))*min(size(q)),size(sc_count,1),max(size(qxx)),max(size(q1))))+0.00000000000000001)));
        K27_ccorr(i,k)=a
    end; 
    toc
end;  
K27_k1 = k1;
K27_qmin = qmin;

imagesc(K27_ccorr);
AdvancedColormap('wr')
text(1,1,num2str(round(K27_ccorr(1,1),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(2,1,num2str(round(K27_ccorr(1,2),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(3,1,num2str(round(K27_ccorr(1,3),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(4,1,num2str(round(K27_ccorr(1,4),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(1,2,num2str(round(K27_ccorr(2,1),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(2,2,num2str(round(K27_ccorr(2,2),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(3,2,num2str(round(K27_ccorr(2,3),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(4,2,num2str(round(K27_ccorr(2,4),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(1,3,num2str(round(K27_ccorr(3,1),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(2,3,num2str(round(K27_ccorr(3,2),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(3,3,num2str(round(K27_ccorr(3,3),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(4,3,num2str(round(K27_ccorr(3,4),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(1,4,num2str(round(K27_ccorr(4,1),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(2,4,num2str(round(K27_ccorr(4,2),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(3,4,num2str(round(K27_ccorr(4,3),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(4,4,num2str(round(K27_ccorr(4,4),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(1,5,num2str(round(K27_ccorr(5,1),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(2,5,num2str(round(K27_ccorr(5,2),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(3,5,num2str(round(K27_ccorr(5,3),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(4,5,num2str(round(K27_ccorr(5,4),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(1,6,num2str(round(K27_ccorr(6,1),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(2,6,num2str(round(K27_ccorr(6,2),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(3,6,num2str(round(K27_ccorr(6,3),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(4,6,num2str(round(K27_ccorr(6,4),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
set(gca,'Xtick',[1,2,3,4],'XtickLabel',{'B','Mono','T','NK'});
set(gca,'Ytick',[1,2,3,4,5,6],'YtickLabel',{'clus 1','clus 2','clus 3','clus 4','clus 5','clus 6'});
saveas(gcf,strcat(char(path0),'/figure/K27_clus_matching.fig'));

%=========================================================================
a= readtable(strcat(char(path),'/scK4_rc_at_25951_bipeaks.txt'),'ReadVariableNames',0);
K4sc_count = table2array(a);
save(char(strcat(char(path),'/scK4_rc_at_25951_bipeaks.mat')),'K4sc_count','-v7.3');
%a =load(strcat(char(path),'/scK4_rc_at_25951_bipeaks.mat'));
%K4sc_count = a.K4sc_count;
a= readtable(strcat(char(path),'/scH3K4me3_7798_cell_counts_at_WBCpeaks.txt'),'ReadVariableNames',1);
sc_count = table2array(a);
save(char(strcat(char(path),'/scH3K4me3_7798_cell_counts_at_WBCpeaks.mat')),'sc_count','-v7.3');
%a =load(strcat(char(path),'/scH3K4me3_7798_cell_counts_at_WBCpeaks.mat'));
%sc_count = a.sc_count;
a = readtable(strcat(char(path),'/wc_K4_uniq.txt'),'ReadVariableNames',0);
depth2 = table2array(a(:,1));
filename = table2array(a(:,2));
sc_cpm = sc_count*1000000./depth2';
K4sc_cpm = K4sc_count*1000000./depth2';
yy = sc_count./depth2';
q = find(sum(yy)<=0.15| sum(sc_count)<1000);
sc_count(:,q)=[];
K4sc_count(:,q)=[];
sc_cpm(:,q)=[];
K4sc_cpm(:,q)=[];
yy(:,q)=[];
depth2(q)=[];
filename(q)=[];
sclog_cpm = log2(sc_cpm+1);
K4sclog_cpm = log2(K4sc_cpm+1);

B = sclog_cpm(:,:);
A = B;
q = find(A<0);
A(q)=0;
q = find(A>0);
A(q)=1;
qdel = find(sum(A')<100 );
B(qdel,:)=[];
cc = corr(B);
dd = max(max(cc)) - cc;
lcc1 = exp(-dd/max(max(dd)));
for i = 1:max(size(cc))
    lcc1(i,i) = 1;
end;    
lcc2 = diag(1./((sum(lcc1)-1).^0.5));
lcc3 = diag(ones(max(size(cc)),1))-lcc2*lcc1*lcc2;
[v1,w1] = eigs(lcc3,max(size(cc)));
tic
ee = randn(max(size(cc)), max(size(cc)));
ee(:)=0;
m=1;
for i = 1:10
    for k = 1:10
        [ia,ib] = kmeans(v1(:,max(size(cc))-k:max(size(cc))),6);
        for j = 1:6
            q = find(ia==j);
            ee(q,q) = ee(q,q) +1;
            
        end;
        m = m+1;
        
    end;
    toc
end;

rng('default')
nclus = 6;
[k1,k2] = kmeans(ee,nclus);


clus_cpmlog={};

for kkkk= 1:10
    qmin = [];
    for i=1:nclus
        q = find(k1==i);
        qmin = [qmin,max(size(q))];
    end;
    npeak2 = size(sc_count,1);
    clus_count = zeros(npeak2,nclus);
    clus_depth = zeros(nclus,1);
    clus_cpm = zeros(npeak2,1);
    for i=1:nclus
        q = find(k1==i);
        rr = randperm(max(size(q)),floor(max(size(q))/3));
        q= q(rr);    
        clus_count(:,i) = sum(sc_count(:,q)')';

        clus_depth(i) = sum(sum(sc_count(:,q)'));
        clus_cpm(:,i)= clus_count(:,i)*1000000./(clus_depth(i));
    end;
    clus_cpmlog{kkkk} = log2(clus_cpm+1);
end;

import bioma.data.DataMatrix

Y{:,1}  = DataMatrix([clus_cpmlog{1}(:,1),clus_cpmlog{2}(:,1),clus_cpmlog{3}(:,1),clus_cpmlog{4}(:,1),clus_cpmlog{5}(:,1),clus_cpmlog{6}(:,1)]);%%40WT
Y{:,2}  = DataMatrix([clus_cpmlog{1}(:,2),clus_cpmlog{2}(:,2),clus_cpmlog{3}(:,2),clus_cpmlog{4}(:,2),clus_cpmlog{5}(:,2),clus_cpmlog{6}(:,2)]);%%40WT
Y{:,3}  = DataMatrix([clus_cpmlog{1}(:,3),clus_cpmlog{2}(:,3),clus_cpmlog{3}(:,3),clus_cpmlog{4}(:,3),clus_cpmlog{5}(:,3),clus_cpmlog{6}(:,3)]);%%40WT
Y{:,4}  = DataMatrix([clus_cpmlog{1}(:,4),clus_cpmlog{2}(:,4),clus_cpmlog{3}(:,4),clus_cpmlog{4}(:,4),clus_cpmlog{5}(:,4),clus_cpmlog{6}(:,4)]);%%40WT
Y{:,5}  = DataMatrix([clus_cpmlog{1}(:,5),clus_cpmlog{2}(:,5),clus_cpmlog{3}(:,5),clus_cpmlog{4}(:,5),clus_cpmlog{5}(:,5),clus_cpmlog{6}(:,5)]);%%40WT
Y{:,6}  = DataMatrix([clus_cpmlog{1}(:,6),clus_cpmlog{2}(:,6),clus_cpmlog{3}(:,6),clus_cpmlog{4}(:,6),clus_cpmlog{5}(:,6),clus_cpmlog{6}(:,6)]);%%40WT

K4_ccorr = zeros(nclus,4);

tic
for i= 1:nclus
    xq = setdiff(1:1:nclus,i);
    PP= ones(52798,nclus-1);
    qt = 1:1:52798;
    for j= 1:nclus-1
       [PValues, TScores] = mattest(Y{:,i}, Y{:,xq(j)});
       PP(:,j) = double(PValues);
       q = find(PP(:,j)<0.05 & mean(Y{:,i}')'-mean(Y{:,xq(j)}')'>0.2);
       qt = intersect(qt,q);
    end;    
    for k = 1:4
        if(k==1)
            qxx = [q_B'];
        elseif(k==2) 
            qxx = [q_Mono'];
        elseif(k==3) 
            qxx = [q_T'];
        elseif(k==4) 
            qxx = [q_NK'];
        end;  
        q1 = qt;
        q = intersect(q1,qxx);
        a = real(-log2((1-sum(hygepdf(0:max(size(q))*min(size(q)),size(sc_count,1),max(size(qxx)),max(size(q1))))+0.00000000000000001)));
        K4_ccorr(i,k)=a
    end; 
    toc
end;   

K4_k1 = k1;
K4_qmin = qmin;


imagesc(K4_ccorr);
AdvancedColormap('wr')
text(1,1,num2str(round(K4_ccorr(1,1),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(2,1,num2str(round(K4_ccorr(1,2),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(3,1,num2str(round(K4_ccorr(1,3),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(4,1,num2str(round(K4_ccorr(1,4),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(1,2,num2str(round(K4_ccorr(2,1),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(2,2,num2str(round(K4_ccorr(2,2),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(3,2,num2str(round(K4_ccorr(2,3),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(4,2,num2str(round(K4_ccorr(2,4),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(1,3,num2str(round(K4_ccorr(3,1),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(2,3,num2str(round(K4_ccorr(3,2),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(3,3,num2str(round(K4_ccorr(3,3),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(4,3,num2str(round(K4_ccorr(3,4),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(1,4,num2str(round(K4_ccorr(4,1),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(2,4,num2str(round(K4_ccorr(4,2),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(3,4,num2str(round(K4_ccorr(4,3),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(4,4,num2str(round(K4_ccorr(4,4),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(1,5,num2str(round(K4_ccorr(5,1),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(2,5,num2str(round(K4_ccorr(5,2),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(3,5,num2str(round(K4_ccorr(5,3),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(4,5,num2str(round(K4_ccorr(5,4),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(1,6,num2str(round(K4_ccorr(6,1),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(2,6,num2str(round(K4_ccorr(6,2),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(3,6,num2str(round(K4_ccorr(6,3),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(4,6,num2str(round(K4_ccorr(6,4),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
set(gca,'Xtick',[1,2,3,4],'XtickLabel',{'B','Mono','T','NK'});
set(gca,'Ytick',[1,2,3,4,5,6],'YtickLabel',{'clus 1','clus 2','clus 3','clus 4','clus 5','clus 6'});
saveas(gcf,strcat(char(path0),'/figure/K4_clus_matching.fig'));
%====================================================================
K4a_count = K4sc_count;
q = find(K4a_count>0);
K4a_count(q)=1;
for i = 1:6
    q1 = find(K4_k1==i);
    K4_acount(:,i) = sum(K4a_count(:,q1)')';
    K4_bulkcount(:,i) = sum(K4sc_count(:,q1)')';
    K4_bulkcpm(:,i) = log2(K4_bulkcount(:,i)*1000000/sum(K4_bulkcount(:,i))+1);
    K4_bulkcv(:,i) = std(K4sclog_cpm(:,q1)')'./(mean(K4sclog_cpm(:,q1)')'+0.0001);
    toc
end;  

K27a_count = K27sc_count;
q = find(K27a_count>0);
K27a_count(q)=1;
for i = 1:6
    q1 = find(K27_k1==i);
    K27_acount(:,i) =sum(K27a_count(:,q1)')';
    K27_bulkcount(:,i) = sum(K27sc_count(:,q1)')';
    K27_bulkcpm(:,i) = log2(K27_bulkcount(:,i)*1000000/sum(K27_bulkcount(:,i))+1);
    K27_bulkcv(:,i) = std(K27sclog_cpm(:,q1)')'./(mean(K27sclog_cpm(:,q1)')'+0.0001);
    toc
end; 

[q4_selx,q4_sely] = find(K4_ccorr>18);
q4_del = setdiff(1:1:nclus,q4_selx);
[q27_selx,q27_sely] = find(K27_ccorr>18);
q27_del = setdiff(1:1:nclus,q27_selx);

K4_bulkcv(:,q4_del)=[];
K27_bulkcv(:,q27_del)=[];
K4_bulkcpm(:,q4_del)=[];
K27_bulkcpm(:,q27_del)=[];
K4_bulkcount(:,q4_del)=[];
K27_bulkcount(:,q27_del)=[];
K4_acount(:,q4_del)=[];
K27_acount(:,q27_del)=[];

rng('default')

K4clus_cpmlog={};
K4nclus = 6;
for kkkk= 1:10
    npeak2 = size(K4sc_count,1);
    clus_count = zeros(npeak2,K4nclus);
    clus_depth = zeros(K4nclus,1);
    clus_cpm = zeros(npeak2,1);
    for i=1:K4nclus
        q = find(K4_k1==i);
        rr = randperm(max(size(q)),floor(max(size(q))/1.5));
        q= q(rr);    
        clus_count(:,i) = sum(K4sc_count(:,q)')';
        clus_depth(i) = sum(sum(K4sc_count(:,q)'));
        clus_cpm(:,i)= clus_count(:,i)*1000000./(clus_depth(i));
    end;
    K4clus_cpmlog{kkkk} = clus_cpm;
end;



K27clus_cpmlog={};
K27nclus = 6;
for kkkk= 1:10
    npeak2 = size(K27sc_count,1);
    clus_count = zeros(npeak2,K27nclus);
    clus_depth = zeros(K27nclus,1);
    clus_cpm = zeros(npeak2,1);
    for i=1:K27nclus
        q = find(K27_k1==i);
        rr = randperm(max(size(q)),floor(max(size(q))/1.5));
        q= q(rr);    
        clus_count(:,i) = sum(K27sc_count(:,q)')';
        clus_depth(i) = sum(sum(K27sc_count(:,q)'));
        clus_cpm(:,i)= clus_count(:,i)*1000000./(clus_depth(i));
    end;
    K27clus_cpmlog{kkkk} = clus_cpm;
end;


import bioma.data.DataMatrix

Y{:,1}  = DataMatrix([K4clus_cpmlog{1}(:,1),K4clus_cpmlog{2}(:,1),K4clus_cpmlog{3}(:,1),K4clus_cpmlog{4}(:,1),K4clus_cpmlog{5}(:,1),K4clus_cpmlog{6}(:,1)]);%%40WT
Y{:,2}  = DataMatrix([K4clus_cpmlog{1}(:,2),K4clus_cpmlog{2}(:,2),K4clus_cpmlog{3}(:,2),K4clus_cpmlog{4}(:,2),K4clus_cpmlog{5}(:,2),K4clus_cpmlog{6}(:,2)]);%%40WT
Y{:,3}  = DataMatrix([K4clus_cpmlog{1}(:,3),K4clus_cpmlog{2}(:,3),K4clus_cpmlog{3}(:,3),K4clus_cpmlog{4}(:,3),K4clus_cpmlog{5}(:,3),K4clus_cpmlog{6}(:,3)]);%%40WT
Y{:,4}  = DataMatrix([K4clus_cpmlog{1}(:,4),K4clus_cpmlog{2}(:,4),K4clus_cpmlog{3}(:,4),K4clus_cpmlog{4}(:,4),K4clus_cpmlog{5}(:,4),K4clus_cpmlog{6}(:,4)]);%%40WT
Y{:,5}  = DataMatrix([K4clus_cpmlog{1}(:,5),K4clus_cpmlog{2}(:,5),K4clus_cpmlog{3}(:,5),K4clus_cpmlog{4}(:,5),K4clus_cpmlog{5}(:,5),K4clus_cpmlog{6}(:,5)]);%%40WT
Y{:,6}  = DataMatrix([K4clus_cpmlog{1}(:,6),K4clus_cpmlog{2}(:,6),K4clus_cpmlog{3}(:,6),K4clus_cpmlog{4}(:,6),K4clus_cpmlog{5}(:,6),K4clus_cpmlog{6}(:,6)]);%%40WT


Z{:,1}  = DataMatrix([K27clus_cpmlog{1}(:,1),K27clus_cpmlog{2}(:,1),K27clus_cpmlog{3}(:,1),K27clus_cpmlog{4}(:,1),K27clus_cpmlog{5}(:,1),K27clus_cpmlog{6}(:,1)]);%%40WT
Z{:,2}  = DataMatrix([K27clus_cpmlog{1}(:,2),K27clus_cpmlog{2}(:,2),K27clus_cpmlog{3}(:,2),K27clus_cpmlog{4}(:,2),K27clus_cpmlog{5}(:,2),K27clus_cpmlog{6}(:,2)]);%%40WT
Z{:,3}  = DataMatrix([K27clus_cpmlog{1}(:,3),K27clus_cpmlog{2}(:,3),K27clus_cpmlog{3}(:,3),K27clus_cpmlog{4}(:,3),K27clus_cpmlog{5}(:,3),K27clus_cpmlog{6}(:,3)]);%%40WT
Z{:,4}  = DataMatrix([K27clus_cpmlog{1}(:,4),K27clus_cpmlog{2}(:,4),K27clus_cpmlog{3}(:,4),K27clus_cpmlog{4}(:,4),K27clus_cpmlog{5}(:,4),K27clus_cpmlog{6}(:,4)]);%%40WT
Z{:,5}  = DataMatrix([K27clus_cpmlog{1}(:,5),K27clus_cpmlog{2}(:,5),K27clus_cpmlog{3}(:,5),K27clus_cpmlog{4}(:,5),K27clus_cpmlog{5}(:,5),K27clus_cpmlog{6}(:,5)]);%%40WT
Z{:,6}  = DataMatrix([K27clus_cpmlog{1}(:,6),K27clus_cpmlog{2}(:,6),K27clus_cpmlog{3}(:,6),K27clus_cpmlog{4}(:,6),K27clus_cpmlog{5}(:,6),K27clus_cpmlog{6}(:,6)]);%%40WT

Y(q4_del)=[];
Z(q27_del)=[];

xxx4 = quantilenorm(K4_bulkcv);
xxx4 = quantilenorm(log2([double(Y{:,1}')',double(Y{:,2}')',double(Y{:,3}')',double(Y{:,4}')']+1));

Y1{:,1} = DataMatrix(xxx4(:,1:6));
Y1{:,2} = DataMatrix(xxx4(:,7:12));
Y1{:,3} = DataMatrix(xxx4(:,13:18));
Y1{:,4} = DataMatrix(xxx4(:,19:24));

K4_xt1={};
K4_xt2={};

for i= 1:4
    K4_xt1{i} = 1:1:max(size(K4_bulkcv));
    K4_xt2{i} = 1:1:max(size(K4_bulkcv));
end;    

for i= 1:4
    xq = setdiff(1:1:4,i);
    for j= 1:3
        [PValues, TScores] = mattest(Y1{:,i}, Y1{:,xq(j)});
        
        fdr1 = mafdr(PValues);
        
        mv1 = mean(double(Y1{:,i})')';
        mv2 = mean(double(Y1{:,xq(j)})')';
        q = find(double(fdr1)<0.05 & mv1-mv2>0.2 & mv1>0);
        K4_xt1{i} = intersect(K4_xt1{i},q);
        q = find(double(fdr1)<0.05 & mv1-mv2<-0.2 & mv1>=0);
        K4_xt2{i} = intersect(K4_xt2{i},q);        
    end;
end;    



xxx27 = quantilenorm(K27_bulkcv);
xxx27 = quantilenorm(log2([double(Z{:,1}')',double(Z{:,2}')',double(Z{:,3}')',double(Z{:,4}')']+1));

Z1{:,1} = DataMatrix(xxx27(:,1:6));
Z1{:,2} = DataMatrix(xxx27(:,7:12));
Z1{:,3} = DataMatrix(xxx27(:,13:18));
Z1{:,4} = DataMatrix(xxx27(:,19:24));

K27_xt1={};
K27_xt2={};

for i= 1:4
    K27_xt1{i} = 1:1:max(size(K27_bulkcv));
    K27_xt2{i} = 1:1:max(size(K27_bulkcv));
end;    
for i= 1:4
    xq = setdiff(1:1:4,i);
    for j= 1:3
       	[PValues, TScores] = mattest(Z1{:,i}, Z1{:,xq(j)});
        fdr1 = mafdr(PValues);
        mv1 = mean(double(Z1{:,i})')';
        mv2 = mean(double(Z1{:,xq(j)})')';
        q = find(double(fdr1)<0.05 & mv1-mv2<-0.2 & mv1>=0);
        K27_xt1{i} = intersect(K27_xt1{i},q);
        q = find(double(fdr1)<0.05  & mv1-mv2>0.2 & mv1>0);
        K27_xt2{i} = intersect(K27_xt2{i},q);        
    end;
end;    


a  = readtable(strcat(char(path),'/bivalent_domain25951_annot_gene.txt'),'ReadVariableNames',0);
b = table2array(a(:,4));
bigene = table2array(a(:,5));

cc_pv1 = zeros(4,4);
for i=1:4
    for j= 1:4
        q = find(abs(b)<=100); %% need to closely overelap to TSS
        q1a = intersect(q,K4_xt1{i});
        q2a = intersect(q,K27_xt1{j});
        q3a = intersect(q1a,q2a);
        q1b = intersect(q,K4_xt2{i});
        q2b = intersect(q,K27_xt2{j}); 
        q3b = intersect(q1b,q2b);
        cc_pv1(i,j) = real(-log2(1-sum(hygepdf(0:max(size(q3a))*min(size(q3a)), max(size(q)), max(size(q1a)), max(size(q2a))))));
    end;
end;    
imagesc(zscore(cc_pv1')');
AdvancedColormap('wr')
text(1,1,num2str(round(cc_pv1(1,1),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(2,1,num2str(round(cc_pv1(1,2),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(3,1,num2str(round(cc_pv1(1,3),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(4,1,num2str(round(cc_pv1(1,4),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(1,2,num2str(round(cc_pv1(2,1),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(2,2,num2str(round(cc_pv1(2,2),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(3,2,num2str(round(cc_pv1(2,3),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(4,2,num2str(round(cc_pv1(2,4),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(1,3,num2str(round(cc_pv1(3,1),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(2,3,num2str(round(cc_pv1(3,2),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(3,3,num2str(round(cc_pv1(3,3),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(4,3,num2str(round(cc_pv1(3,4),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(1,4,num2str(round(cc_pv1(4,1),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(2,4,num2str(round(cc_pv1(4,2),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(3,4,num2str(round(cc_pv1(4,3),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(4,4,num2str(round(cc_pv1(4,4),2,'significant')),'fontsize',16,'HorizontalAlignment','center')

x27clus = {'K27clus1','K27clus2','K27clus3','K27clus4','K27clus5','K27clus6'};
x4clus = {'K4clus1','K4clus2','K4clus3','K4clus4','K4clus5','K4clus6'};
x27clus(q27_del)=[];
x4clus(q4_del)=[];
set(gca,'Xtick',[1,2,3,4],'XtickLabel',x27clus);
set(gca,'Ytick',[1,2,3,4],'YtickLabel',x4clus);
saveas(gcf,strcat(char(path0),'/figure/K4clus_to_K27clus_matching_heatmap.fig'));




axk4=mean(K4_acount)*0.05;
axk27=mean(K27_acount)*0.05;
xxx4 = quantilenorm(log2(K4_bulkcv+1));
K4_xt1={};
for i= 1:4
    K4_xt1{i} = 1:1:max(size(K4_bulkcv));
end;    

for i= 1:4
    xq = setdiff(1:1:4,i);
    for j= 1:3
        q = find( xxx4(:,i)-xxx4(:,xq(j))>0.19 & xxx4(:,i)>0 & K4_acount(:,i)>=axk4(i)); % 0.19 or 0.2
        K4_xt1{i} = intersect(K4_xt1{i},q);
     
    end;
end;    

xxx27 = quantilenorm(log2(K27_bulkcv+1));
K27_xt2={};
for i= 1:4
    K27_xt2{i} = 1:1:max(size(K27_bulkcv));
end;    
for i= 1:4
    xq = setdiff(1:1:4,i);
    for j= 1:3
        q = find(xxx27(:,i)-xxx27(:,xq(j))>0.19 & xxx27(:,i)>0 & K27_acount(:,i)>=axk27(i));% 0.19 or 0.2
        K27_xt2{i} = intersect(K27_xt2{i},q);        
    end;
end;    
xxxx= quantilenorm([K4_bulkcv,K27_bulkcv]);
cc_pv1 = zeros(4,4);
cc_pv2 = zeros(4,4);
cc_pv3 = zeros(4,4);
for i=1:4
    for j= 1:4
        q = find(abs(b)>=0 );
        q1a = intersect(q,K4_xt1{i});
        q2a = intersect(q,K27_xt2{j});
        q3a = intersect(q1a,q2a);
        %scatter(xxxx(q3a,i),xxxx(q3a,j+4));
        cc_pv(i,j) = corr(xxxx(q3a,i),xxxx(q3a,j+4),'type','Spearman');
    end;
end;    
imagesc(zscore(cc_pv')');
AdvancedColormap('wr')
text(1,1,num2str(round(cc_pv(1,1),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(2,1,num2str(round(cc_pv(1,2),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(3,1,num2str(round(cc_pv(1,3),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(4,1,num2str(round(cc_pv(1,4),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(1,2,num2str(round(cc_pv(2,1),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(2,2,num2str(round(cc_pv(2,2),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(3,2,num2str(round(cc_pv(2,3),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(4,2,num2str(round(cc_pv(2,4),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(1,3,num2str(round(cc_pv(3,1),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(2,3,num2str(round(cc_pv(3,2),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(3,3,num2str(round(cc_pv(3,3),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(4,3,num2str(round(cc_pv(3,4),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(1,4,num2str(round(cc_pv(4,1),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(2,4,num2str(round(cc_pv(4,2),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(3,4,num2str(round(cc_pv(4,3),2,'significant')),'fontsize',16,'HorizontalAlignment','center')
text(4,4,num2str(round(cc_pv(4,4),2,'significant')),'fontsize',16,'HorizontalAlignment','center')


x27clus = {'K27clus1','K27clus2','K27clus3','K27clus4','K27clus5','K27clus6'};
x4clus = {'K4clus1','K4clus2','K4clus3','K4clus4','K4clus5','K4clus6'};
x27clus(q27_del)=[];
x4clus(q4_del)=[];
set(gca,'Xtick',[1,2,3,4],'XtickLabel',x27clus);
set(gca,'Ytick',[1,2,3,4],'YtickLabel',x4clus);
xlabel('K27 Cell-to-cell variation');
ylabel('K4 Cell-to-cell variation');
saveas(gcf,strcat(char(path0),'/figure/K4CV_to_K27CV_corr_heatmap.fig'));




