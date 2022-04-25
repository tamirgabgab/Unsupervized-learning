clc;
close all;
clear;
% create data from differnt densities
Flag_plot=1;
clusterN=3;
N=1000;
%creats 3-d gaussian 
[X,lable] = Creat_dataset(Flag_plot,N,clusterN);
a=1;


%% ------------- MLE algorithm performence - Home Work - excersize 1 %%%%%%%
%stage 1 initalization 
for k=1:clusterN %for each cluster
     P_w_guess(k)=1/clusterN; %set initial guess
end
for k=1:clusterN %for each cluster
     sigma_guess(:,:,k)=eye(size(X,2));%set initial sigma[if x(m,n),creat I(n,n)]
end
mu_guess=datasample(X,clusterN,1); %Random sample choice[y = datasample(data,k,dim)- k observations sampled] 

% mu_guess=k_means(X,mu_guess,clusterN); %Intialization using k-means


%MLE estimation function:
[ u_mle,mu_est,sigma_est,P_w_est ] = mle_est( X,clusterN,mu_guess,sigma_guess,P_w_guess ); %MLE
%Graphic display - MLE friendships


figure(2);
t=1:size(X,1);
colors={{'r'},{'g'},{'cyan'},{'magenta'},{'b'},{'k'},{'yellow'}};
for k=1:clusterN  %for each cluster
    color_ind=mod(k-1,length(colors))+1;
    subplot(clusterN,1,k);
    bar(t,u_mle(k,:),cell2mat(colors{color_ind}));
    xlabel('j - Samples');
    ylabel(['U',num2str(k),'j']);
    title(['Friendships values for group dataset 6 random init ',num2str(k),'- MLE']);
end