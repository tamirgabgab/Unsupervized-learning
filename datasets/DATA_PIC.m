
clear all
close all


COLORS={{[0 0.4470 0.7410]},{[0.8500 0.3250 0.0980]},{[0.4940 0.1840 0.5560]},...
        {[0.4660 0.6740 0.1880]},{[0.3010 0.7450 0.9330]},{[0.6350 0.0780 0.1840]},...
        {'k'},{'r'},{'b'},('m')};
%{
path='C:\Users\tamir\Desktop\תמיר\matlab scrifts\datasets\report\temp\';
names={'D','a1','t1','a2','s','e','t2'};
X=[];
Labels=[];
vec_add=2*[0 0 0;30 0 0;60 0 0;90 0 0;...
         200 0 0;230 0 0;260 0 0];

r=linspace(35,35,length(names));
theta=linspace(0,2*pi*(1-1/length(names)),length(names)); 
z=linspace(0,10,length(names));

for i=1:length(names)
    vec_add(i,:)=[r(i)*cos(theta(i)) r(i)*sin(theta(i)) z(i)];
end

         
vec_add=vec_add([3 2 1 7 6 5 4],:);

         
for i=1:length(names)    
    DATASET=load([path names{i}]);
    DATASET=DATASET.data;
    X=cat(1,X,vec_add(i,:)+DATASET.X);
    %X=X*(i^2/4);
    Labels=cat(1,Labels,i*DATASET.Labels);
end
X(:,2)=10*X(:,2);
figure(1);
for k=1:max(Labels)
    pp=find(Labels==k);
    plot3(X(pp,1),X(pp,2),X(pp,3),'o','color',cell2mat(colors{k}),...
          'LineWidth',1,'MarkerSize',4);
    hold on;  
end
title('Real clusters');
grid on;
xlabel('x');
ylabel('y');
zlabel('z');

data.name='dataset-letters';
data.X=X;
data.Labels=Labels;
data.n=length(Labels);
data.d=size(X,2);
data.c=length(names);
%}


A=imread('clusters.jpg');
A=imresize(A,1,'nearest');

res=10;
c=2;
d=3;
thresh=multithresh(A,c);
V=imquantize(A,thresh);
n=size(A);
eps=0.1;
k=[1 1]*0.025;
edge=[-10 10;5 15;0 5];


cluster_colors=rand(c,3);

index=1;
X=[];
for i=1:2:n(1)
    for j=1:2:n(2)
        
        if V(i,j,1)<length(thresh)+1 ||...
           V(i,j,2)<length(thresh)+1 ||...
           V(i,j,3)<length(thresh)+1
       
           state=[V(i,j,1) V(i,j,2) V(i,j,3)]-1;
           RGB(index,:)=state(1)*c^2+state(2)*c+state(3);
           X(index,:)=[j/res (n(1)-i)/res];
           index=index+1;
        end      
    end
end


for i=1:length(X)
    m=1;
    for j=1:length(X)
        if m>=16
            RGB(i)=mode(colors);
            m=1;
            break;
        end
        if norm(X(i,:)-X(j,:))<1
            colors(m)=RGB(j);
            m=m+1;
            continue;
        end

    end
end

%{
figure(1);
for i=1:length(X)
    Labels(i,1)=find(unique(RGB)==RGB(i));
    plot(X(i,1),X(i,2),'o','color',cluster_colors(Labels(i,1),:));
    hold on;
end
%}


numbers=unique(RGB);
Y=[];
a=1;
for i=1:c
    index=find(RGB==numbers(i));
    for j=1:length(index)
        Y=cat(1,Y,X(index(j),:)+eps*(2*rand(1,2)-1));
        Labels(a,1)=i;
        a=a+1;
    end
end

Z=[];
Z_Labels=[];
cluster_samples=zeros(1,c);
for j=1:c
    indexes=find(Labels==j);
    cluster_samples(j)=length(indexes(1):floor(1/k(j)):indexes(end));
    Z=cat(1,Z,Y(indexes(1):floor(1/k(j)):indexes(end),:));
    Z_Labels=cat(1,Z_Labels,j*ones(cluster_samples(j),1));
end

t=0*(linspace(0,4*pi,length(Z)))';
Y=cat(2,Z,t+0.01*rand(length(t),1));
Labels=Z_Labels;

m=ones(length(Y),1)*(edge(:,2)'-edge(:,1)')./(max(Y)-min(Y));   %Normalization
Y=ones(length(Y),1)*edge(:,1)'+m.*(Y-ones(length(Y),1)*min(Y));

figure(2);
grid on;
for i=1:length(Y)
    %pause(1/length(Y));
    %{
    if Labels(i)==1
        Y(i,1)=Y(i,1)+5;
        Y(i,2)=Y(i,2)-1;
        Y(i,3)=Y(i,3)+0.5+sin(2*i);
    end
    if Labels(i)==2
        Y(i,1)=Y(i,1)-2;
        Y(i,2)=Y(i,2)-1;
        Y(i,3)=Y(i,3)+1;
    end
    %}
    plot3(Y(i,1),Y(i,2),Y(i,3),'o','color',cluster_colors(Labels(i),:),...
        'LineWidth',1,'MarkerSize',4);
    hold on;
end
grid on;
xlabel('x');
ylabel('y');
zlabel('z');
xlim([min(Y(:,1))-1 max(Y(:,1))+1]);
ylim([min(Y(:,2))-1 max(Y(:,2))+1]);
zlim([min(Y(:,3))-1 max(Y(:,3))+1]);


data.name='intersections';
data.X=Y;
data.Labels=Labels;
data.n=length(Y);
data.d=d;
data.c=c;


clear i j k t m a A colors index numbers state RGB X Z Z_Labels V
%}


