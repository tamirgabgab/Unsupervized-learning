

clear all
clf;

d=3;       % Dim of samples
c=3;       % Number of clusters
center=[0 0 0;8 5 -2;-2 0 -1];
r=[10 2 2];
h=[10 5 3];
k=[1 0.5 0.4];
%edge=[-10 10;0 10;10 40];
eps=0.001;
res=10;   % res rings
res_space=linspace(4,20,res-1);

cluster_colors=rand(c,3);

noise=100;
noise_range=[-10,10];
noise_thr=4;

X=[];
Labels=[];

for j=1:c
    a=1;
    for dr=0:r(j)/res:r(j)*(res-1)/res
        if dr==0
            X=cat(1,X,center(j,:)-[0 0 h(j)/2]);
            Labels=cat(1,Labels,j);
            continue;
        end
        for theta=0:2*pi/res_space(a):2*pi*(res_space(a)-1)/res_space(a)
            X=cat(1,X,[dr*cos(theta) dr*sin(theta) -h(j)/2]+center(j,:));
            Labels=cat(1,Labels,j);
        end
        a=a+1;
    end
    a=res-1;
    for dz=-h(j)/2:h(j)/(res_space(a)-1):h(j)/2
        for theta=0:2*pi/res_space(a):2*pi*(res_space(a)-1)/res_space(a)
            X=cat(1,X,[r(j)*cos(theta) r(j)*sin(theta) dz]+center(j,:));
            Labels=cat(1,Labels,j);
        end
    end
    
    
    for dr=r(j)*(res-1)/res:-r(j)/res:0
        if dr<=0
            X=cat(1,X,center(j,:)+[0 0 h(j)/2]);
            Labels=cat(1,Labels,j);
            a=a+1;
            continue;
        end
        for theta=0:2*pi/res_space(a):2*pi*(res_space(a)-1)/res_space(a)
            X=cat(1,X,[dr*cos(theta) dr*sin(theta) +h(j)/2]+center(j,:));
            Labels=cat(1,Labels,j);
        end
        a=a-1;
    end
end

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
    pause(1/length(Y));
    if Y_Labels(i)>0    
        plot3(Y(i,1),Y(i,2),Y(i,3),'o','color',cluster_colors(Y_Labels(i),:),...
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

data.name='Galil';
data.X=Y;
data.Labels=Y_Labels;
data.n=length(Y);
data.d=d;
data.c=c;
data.cluster_samples=cluster_samples;

clear a cluster colors dr dz i j k Labels n noise_p
clear norms_points point p X iter indexes

