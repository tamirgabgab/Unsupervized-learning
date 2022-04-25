

clear all
clf;

d=3;       % Dim of samples
c=3;       % Number of clusters
center=[0 0 0;0 1 2;4 2 0];
r=[6;3;2];
res=[200;40;20]; 
eps=0.2;

diff=max(r)./r;

cluster_samples=randi([50,150],1,2);
%n=sum(cluster_samples);   % Number of samples
cluster_colors=rand(c,3);

dens=10;
noise=0;
noise_range=[-5,15];
noise_thr=5;

X=[];
Labels=[];
figure(1);
for j=1:c
    for theta=0:2*pi/res(j):diff(j)*2*pi*(res(j)-1+1)*4/res(j)
    dev=eps*(2*rand(1,3)-1);
    ADD=r(j,1)*[cos(theta) sin(theta) theta/4]+center(j,:)+dev;
    X=cat(1,X,ADD);
    Labels=cat(1,Labels,j);
    plot3(ADD(:,1),ADD(:,2),ADD(:,3),'o','color',cluster_colors(j,:),...
        'LineWidth',1,'MarkerSize',4);
    hold on;
    end
end
grid on;
n=length(X);


iter=1;
noise_p=[];
norms_points=zeros(1,c);
while length(noise_p)<noise && iter<500 && noise>0
    point=rand(1,d)*(noise_range(1)-noise_range(2))+noise_range(2);
    for p=1:c
        norms_points(p)=norm(point-center(p,:));
    end
    if min(norms_points)>noise_thr
       noise_p=cat(1,noise_p,point);
    end
    iter=iter+1;
end
if noise>0
    plot3(noise_p(:,1),noise_p(:,2),noise_p(:,3),'hr',...
            'LineWidth',1,'MarkerSize',2);
    X=cat(1,X,noise_p);
    Labels=cat(1,Labels,-ones(noise,1));
    n=n+length(noise_p);
end

zlim([0 20]);
%ax=gca;
%ax.ZTick=-10:2:10;

clear s theta mu ADD

data.name='Gaos3 noise';
data.X=X;
data.Labels=Labels;
data.n=n;
data.d=d;
data.c=c;
data.cluster_samples=cluster_samples;

%clear i s theta mu_G mu sigma ADD cluster_colors cluster_samples
%clear iter dens noise noise_p p point norms_points