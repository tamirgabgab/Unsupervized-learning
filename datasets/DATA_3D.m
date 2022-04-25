
close all
clear all
clf;

d=4;       % Dim of samples
g=4;         % Number of gaos
mu_G=[0 0 0 1 ;...
      10 0 1 0;....
      0 10 0 1;
      10 10 0 1];
      %0 5 5];
sigma_G(:,:,:,1)=1.5*eye(d);
sigma_G(:,:,:,2)=1.5*eye(d);
sigma_G(:,:,:,3)=1.5*eye(d);
sigma_G(:,:,:,4)=1.5*eye(d);
%sigma_G(:,:,:,5)=[1 2 0;0 4 2; 0 2 2];


cluster_samples=[60,80,90,50];
n=sum(cluster_samples);   % Number of samples
cluster_colors=rand(g,3);

colors={{'r'},{'k'},{[0.9290 0.6940 0.1250]},{'magenta'},{'b'},{'g'},{'yellow'}};

dens=30;
noise=0;
noise_range=[-2,10];
noise_thr=3;

X=[];
Labels=[];
figure(1);
for i=1:g
    mu=mu_G(i,:);
    %sigma=0.5*rand(d,d)*0.5+diag((0+rand(d,1))/dens);
    %sigma=sqrt(sigma_G(:,:,:,i)*sigma_G(:,:,:,i)');
    sigma=sigma_G(:,:,:,i)*sigma_G(:,:,:,i)';
    sigma_G(:,:,:,i)=sigma;
    ADD=mvnrnd(mu,sigma,cluster_samples(i));
    X=cat(1,X,ADD);
    Labels=cat(1,Labels,i*ones(cluster_samples(i),1));
    plot3(ADD(:,1),ADD(:,2),ADD(:,3),'o','color',cell2mat(colors{i}),...
        'LineWidth',1,'MarkerSize',4);
    hold on;
end
grid on;

xlabel('x');
ylabel('y');
zlabel('z');

iter=1;
noise_p=[];
norms_points=zeros(1,g);
while length(noise_p)<noise && iter<500 && noise>0
    point=rand(1,d)*(noise_range(1)-noise_range(2))+noise_range(2);
    for p=1:g
        norms_points(p)=norm(point-mu_G(p,:));
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


clear s theta mu ADD

data.name='5_Gaos_3d';
data.X=X;
data.Labels=Labels;
data.n=n;
data.d=d;
data.c=length(cluster_samples);
data.cluster_samples=cluster_samples;
data.mu=mu_G;
data.sigma=sigma_G;


clear i s theta mu sigma ADD cluster_colors cluster_samples
clear iter dens noise noise_p p point norms_points