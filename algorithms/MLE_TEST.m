

clear all
close all

DATASET=load('elipse.mat');
DATASET=DATASET.data;
X=DATASET.X;
X=X(1:2900,:);
Labels=DATASET.Labels';
n=DATASET.n;
d=DATASET.d;
c=DATASET.c;
%mu_correct=DATASET.mu_correct;
%sigma_correct=DATASET.sigma_correct;
colors={{'r'},{'k'},{'b'},{'magenta'},{'cyan'},{'g'},{'yellow'},{'#A2142F'}};

u_mle_real=zeros(c,length(X));
for i=1:length(X)
    u_mle_real(Labels(i),i)=1;
end

figure(1);
for k=1:length(X)
    plot3(X(k,1),X(k,2),X(k,3),'o','color',cell2mat(colors{Labels(k)}),...
            'LineWidth',1,'MarkerSize',4);
    hold on;  
end
title('Real clusters');
grid on;


%% ------------- MLE algorithm performence - Home Work - excersize 1 %%%%%%%
%stage 1 initalization 
for k=1:c %for each cluster
     P_w_guess(k)=1/c; %set initial guess
end
for k=1:c %for each cluster
     sigma_guess(:,:,k)=eye(size(X,2));%set initial sigma[if x(m,n),creat I(n,n)]
end

mu_guess=datasample(X,c,1); %Random sample choice[y = datasample(data,k,dim)- k observations sampled] 

mu_guess=k_means(X,mu_guess,c); %Intialization using k-means


%MLE estimation function:
[ u_mle,mu_est,sigma_est,P_w_est ] = mle_est( X,c,mu_guess,sigma_guess,P_w_guess ); %MLE
%Graphic display - MLE friendships

%error_mu=norm((mu_est-mu_correct)./mu_correct,'fro');
%error_sigma=0;
%for k=1:c
%    error_sigma=error_sigma+norm((sigma_est(:,:,k)-sigma_correct(:,:,k))./sigma_correct(:,:,k),'fro');
%end


for i=1:length(u_mle)
    est_Label(i)=find(u_mle(:,i)==max(u_mle(:,i)));
end

mone=0;
mehane=n*(n-1)/2;
for i=1:length(X)
    for j=i+1:length(X)
        mone=mone+double(not(xor(Labels(i)==Labels(j),est_Label(i)==est_Label(j))));
    end
end
RI=mone/mehane;

figure(2);
for k=1:length(X)
    plot3(X(k,1),X(k,2),X(k,3),'o','color',cell2mat(colors{est_Label(k)}),...
            'LineWidth',1,'MarkerSize',4);
    hold on;  
end
grid on;

figure(3);
t=1:size(X,1);
for k=1:c  %for each cluster
    subplot(c,1,k);
    bar(t,u_mle_real(k,:),cell2mat(colors{k}));
    xlabel('j - Samples');
    ylabel(['U',num2str(k),'j']);
    title('Friendships values');
end
title('Real clusters');

figure(4);
for k=1:c  %for each cluster
    color_ind=mod(k-1,length(colors))+1;
    subplot(c,1,k);
    bar(t,u_mle(k,:),cell2mat(colors{color_ind}));
    xlabel('j - Samples');
    ylabel(['U',num2str(k),'j']);
    title(['Friendships values for group dataset 6 random init ',num2str(k),'- MLE']);
end