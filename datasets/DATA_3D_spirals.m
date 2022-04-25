
close all
clear all
clf;

d=3;       % Dim of samples
g=3;         % Number of gaos

start=0*[0 0 0;-5 0 0;30 0 0;];
r=[0 50;15 85;30 100];
z=[-5 5;-10 10;-20 20];
theta=[0 2*pi;0 2*pi;0 1.5*pi];
res=100;

eps=1.5;

X=[];
Labels=[];
for k=1:g
    R=linspace(r(k,1),r(k,2),res);
    THETA=linspace(theta(k,1),theta(k,2),res);
    Z=linspace(z(k,1),z(k,2),res);
    for t=1:res
        RANDOM=(rand(1,d)*2-1)*eps;
        X=cat(1,X,start(k,:)+RANDOM+...
            [R(t)*cos(THETA(t)) R(t)*sin(THETA(t)) Z(t)]);
        Labels=cat(1,Labels,k);
    end
end


colors={{'r'},{'k'},{[0.9290 0.6940 0.1250]},{'magenta'},{'b'},{'g'},{'yellow'}};
figure();
for i=1:g
    ADD=X(Labels==i,:);
    cluster_samples(i)=length(find(Labels==i));
    plot3(ADD(:,1),ADD(:,2),ADD(:,3),'o','color',cell2mat(colors{i}),...
        'LineWidth',1,'MarkerSize',4);
    hold on;
end
grid on;

xlabel('x');
ylabel('y');
zlabel('z');


data.name='spi';
data.X=X;
data.Labels=Labels;
data.n=size(X,1);
data.d=d;
data.c=length(cluster_samples);
data.cluster_samples=cluster_samples;


clear i s theta mu sigma ADD cluster_colors cluster_samples
clear iter dens noise noise_p p point norms_points