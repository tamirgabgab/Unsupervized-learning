

clear all
clf;

d=3;       % Dim of samples
c=4;       % Number of clusters
center=[0 0 0;0 0 0;12 7 -7;7 0 -12];
r=[8 8 8;2 2 2;6 6 6;3 3 3];
rotation=[0 0 0;6 1 3;6 3 2;2 2 1];
k=1*[1 1 1 1];
%edge=[-10 10;0 10;10 40];
eps=0.1;
res_theta=50;   % res balls
res_phi=50;
cluster_colors=rand(c,3);
colors={{'r'},{'g'},{'k'},{'magenta'},{'b'},{'cyan'},{'yellow'}};

noise=0;
noise_range=[-10,10];
noise_thr=4;

X=[];
Labels=[];


theta=linspace(-pi/2,pi/2,res_theta);
phi=linspace(0,2*pi,res_phi);

for j=1:c
    R_A=[1 0 0;...
        0 cos(rotation(j,1)) sin(rotation(j,1));...
        0 -sin(rotation(j,1)) cos(rotation(j,1))];
    R_B=[cos(rotation(j,2)) 0 -sin(rotation(j,2));...
         0 1 0;...
         sin(rotation(j,2)) 0 cos(rotation(j,2))];
    R_C=[cos(rotation(j,3)) sin(rotation(j,3)) 0;...
         -sin(rotation(j,3)) cos(rotation(j,3)) 0;...
         0 0 1];
    for t=1:length(theta)
        for p=1:length(phi)
            ADD=[r(j,1)*cos(theta(t))*cos(phi(p))...
                 r(j,2)*cos(theta(t))*sin(phi(p))...
                 r(j,3)*sin(theta(t))]+center(j,:);
            ADD=(R_A*R_B*R_C*ADD')'+center(j,:);
            X=cat(1,X,ADD);
            Labels=cat(1,Labels,j);
        end
    end
end

%{
a=1;
for i=1:length(X)
    for j=1:length(X)
        if i~=j
            norms(a)=norm(X(i,:)-X(j,:));
            a=a+1;
        end
    end
end
%}

Y=[];
Y_Labels=[];
cluster_samples=zeros(1,c);
for j=1:c
    indexes=find(Labels==j);
    cluster_samples(j)=length(indexes(1):floor(1/k(j)):indexes(end));
    Y=cat(1,Y,X(indexes(1):floor(1/k(j)):indexes(end),:));
    Y_Labels=cat(1,Y_Labels,j*ones(cluster_samples(j),1));
end

Y=Y+eps*(2*rand(length(Y),d)-1);

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
    %pause(5/noise);
    %plot3(noise_p(:,1),noise_p(:,2),noise_p(:,3),'hr',...
    %        'LineWidth',1,'MarkerSize',2);
    Y=cat(1,Y,noise_p);
    Y_Labels=cat(1,Y_Labels,-ones(noise,1));
end


%m=ones(length(Y),1)*(edge(:,2)'-edge(:,1)')./(max(Y)-min(Y));     %Normalization
%Y=ones(length(Y),1)*edge(:,1)'+m.*(Y-ones(length(Y),1)*min(Y));


figure(1);
for i=1:length(Y)
    %pause(1/length(Y));
    if Y_Labels(i)>0    
        plot3(Y(i,1),Y(i,2),Y(i,3),'o','color',cell2mat(colors{Y_Labels(i)}),...
            'LineWidth',1,'MarkerSize',4);
    else
        plot3(noise_p(:,1),noise_p(:,2),noise_p(:,3),'hr',...
        'LineWidth',1,'MarkerSize',2);
    end
    hold on;
end
grid on;
xlabel('x');
ylabel('y');
zlabel('z');


clear s theta mu ADD

data.name='Elipse';
data.X=Y;
data.Labels=Y_Labels;
data.n=length(Y);
data.d=d;
data.c=c;
data.cluster_samples=cluster_samples;

clear a cluster colors dr dz i j k Labels n noise_p
clear norms_points point p iter indexes

