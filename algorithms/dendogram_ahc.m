
close all
clear all

load fisheriris;

data.name='fisher';
data.X=meas;
data.Labels=[ones(50,1);2*ones(50,1);3*ones(50,1)];
data.n=150;
data.d=4;
data.c=3;
data.cluster_samples=50*[1 1 1];


path='C:\Users\tamir\Desktop\תמיר\matlab scrifts\datasets\';
name='AHC_data';
DATASET=load([path name]);
DATASET=DATASET.data;
X=DATASET.X;
Labels=DATASET.Labels';
n=DATASET.n;
d=DATASET.d;
c=DATASET.c;

        % #1 - minimum of the euclidian distance between 2 samples in same
        %         groups = 'signle'
        % #2 - maximum of the euclidian distance between 2 samples from diffrent
        %         groups = 'complete'
        % #3 - Avereage distance between the clusters = 'average'
        % #4 - Distance between the means of the group = 'centroid'
        % #5 - Sum of squared error criterion = 'ward'
        % #6 - The majority from all 5 criterias
        
Z=linkage(X,'single');
T=cluster(Z,'maxclust',c);

figure();
cutoff=median([Z(end-2,3) Z(end-1,3)]);
dendrogram(Z,'ColorThreshold',cutoff);
title('Dendogram for AHC with X dataset');
grid on;

