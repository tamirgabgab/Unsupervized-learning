clc;
close all;
clear;


path='C:\Users\tamir\Desktop\תמיר\matlab scrifts\datasets\';
name='AHC_data';
DATASET=load([path name]);
DATASET=DATASET.data;
X=DATASET.X;
Labels=DATASET.Labels';
n=DATASET.n;
d=DATASET.d;
c=DATASET.c;

%P_w_correct=(DATASET.cluster_samples/n);

for i=1:max(Labels)
    P_w_correct(i)=length(find(Labels==i))/n;
end


%{
colors={{[1 0 0]},{'k'},{'b'},{'magenta'},{'green'},...
    {[0 0.4470 0.7410]},{[0.8500 0.3250 0.0980]}...
    ,{[0.9290 0.6940 0.1250]},{[0.4940 0.1840 0.5560]}...
    ,{[0.4660 0.6740 0.1880]},{[0.3010 0.7450 0.9330]},....
    {[0.6350 0.0780 0.1840]},{[0.8790 0.2640 0.5450]}};
%}

colors={{'r'},{'k'},{'b'},{'magenta'},{'cyan'},{'g'},{'yellow'}};

c_real=max(Labels);
u_real=zeros(c_real,length(X));
for i=1:length(X)
    u_real(Labels(i),i)=1;
end


mu_correct=zeros(c_real,d);
sigma_correct=zeros(d,d,c_real);
for k=1:c_real
    index=find(u_real(k,:)==1);
    X_cluster=X(index,:);
    mu_correct(k,:)=round(mean(X_cluster));
    sigma_correct(:,:,k)=round(cov(X_cluster));
    hold on;  
end



figure(1);
for k=1:c_real
    pp=find(u_real(k,:)==1);
    if d==2
        plot(X(pp,1),X(pp,2),'o','color',cell2mat(colors{k}),...
                'LineWidth',1,'MarkerSize',4);
    else
        plot3(X(pp,1),X(pp,2),X(pp,3),'o','color',cell2mat(colors{k}),...
                'LineWidth',1,'MarkerSize',4);
    end
    hold on;  
end
title('Real clusters');
grid on;
xlabel('x');
ylabel('y');
zlabel('z');


figure(2);
t=1:size(X,1);
title('Real bars');
for k=1:c_real
    subplot(c_real,1,k);
    bar(t,u_real(k,:),'FaceColor',cell2mat(colors{k}));
    xlabel('j - Samples');
    ylabel(['U',num2str(k),'j']);
    ylim([0 1]);
end



%% ------------- MLE algorithm performence - Home Work - excersize 1 %%%%%%%
%stage 1 initalization

%{
c_est=c_real;
for k=1:c_est %for each cluster
     P_w_guess(k)=1/c_real; %set initial guess
end
for k=1:c_est %for each cluster
     sigma_guess(:,:,k)=eye(size(X,2));%set initial sigma[if x(m,n),creat I(n,n)]
end

mu_guess=datasample(X,c_est,1); %Random sample choice[y = datasample(data,k,dim)- k observations sampled] 

mu_guess=k_means(X,mu_guess,c_est); %Intialization using k-means

ITER=30;
%MLE estimation function:
tic
[ u_mle,mu_est,sigma_est,P_w_est ] = mle_est( X,c_est,mu_guess,sigma_guess,P_w_guess, ITER); %MLE
toc;

if c_est==c_real
    yadani=zeros(1,c_real);
    for i=1:length(yadani)
        idx=find(abs(P_w_correct-P_w_est(i))==...
                 min(abs(P_w_correct-P_w_est(i))));
        yadani(i)=idx;
    end
end

if length(yadani)~=length(unique(yadani))
    fprintf('error!!!!! line 110\n');
    yadani=1:c_est;
end



u_est=zeros(c_real,size(X,1));
LLK=1;
for i=1:length(X)
    index(i)=find(u_mle(:,i)==max(u_mle(:,i)));
    u_est(index(i),i)=1;
    LLK=LLK*u_mle(index(i),i);
end
LLK=-log(LLK);

figure(3);
for k=1:c_est
    pp=find(index==k);
    if d==2
        plot(X(pp,1),X(pp,2),'o','color',cell2mat(colors{yadani(k)}),...
                'LineWidth',1,'MarkerSize',4);
    else
        plot3(X(pp,1),X(pp,2),X(pp,3),'o','color',cell2mat(colors{yadani(k)}),...
                'LineWidth',1,'MarkerSize',4);
    end
    hold on;  
end
grid on;
xlabel('x');
ylabel('y');
zlabel('z');


%Graphic display - MLE friendships
figure(4);
t=1:size(X,1);
%colors={{'r'},{'g'},{'cyan'},{'magenta'},{'b'},{'k'},{'yellow'}};
u_z=u_mle(yadani,:);
for k=1:c_est  %for each cluster
    color_ind=mod(k-1,length(colors))+1;
    subplot(c_est,1,k);
    bar(t,u_z(k,:),cell2mat(colors{color_ind}));
    xlabel('j - Samples');
    ylabel(['U',num2str(k),'j']);
    title(['Friendships values for group dataset 6 random init ',num2str(k),'- MLE']);
end

    
    %yadani=[4 2 3 1];
    p_error=norm(P_w_correct-P_w_est);
    mu_error=norm(mu_correct-mu_est(yadani,:));
    sigma_error=0;
    for i=1:c_est
        sigma_error=sigma_error...
                   +norm(sigma_correct(:,:,i)-sigma_est(:,:,yadani(i)));
    end

    error=[p_error mu_error sigma_error];
    
    mone=0;
    for i=1:length(X)
        for j=i+1:length(X)
            mone=mone+double(not(xor(Labels(i)==Labels(j),index(i)==index(j))));
        end
    end
    RI=mone*2/(length(X)*(length(X)-1));
    NMI=nmi(index,Labels);
    
%}


%% ----- UOFC algorithm performance - Home Work - excersize 2 %%%%%
%{
     c_max=c+2; %Maximum groups we want to classify our data 
     %UOFC function call
    [u,mu_est,FHV,DPA,PD,J_k,Invariant_criterion,APDM]=UOFC(X,c_max,2); %UOFC
   
    %Graphic display -UOFC friendships
    figure(3); 
        for k=1:c_max
            color_ind=mod(k-1,length(colors))+1;
            subplot(c_max,1,k);
            bar(t,u(k,:),'FaceColor',cell2mat(colors{k}));
            xlabel('j - Samples');
            ylabel(['U',num2str(k),'j']);
            title(['Friendships values for group - dataset 6',num2str(k),'- UOFC']);
        end

    %--------------------Validity criterions
    figure(4); 
    subplot(3,2,1);
        plot(FHV,'r');
        xlabel('k-number of groups');
        ylabel('FHV');
        title('Fuzzy hypervolume (FHV) - eq (1)');
        grid on;
    subplot(3,2,3);
        plot(DPA,'r');
        xlabel('k-number of groups');
        ylabel('DPA');
        title('Avereage partition density (APDC) - eq (3)');
        grid on;
    subplot(3,2,2);
        plot(PD,'r');
        xlabel('k-number of groups');
        ylabel('PD');
        title('Partition density (PD) - eq (2) ');
        grid on;
    subplot(3,2,5);
        plot(J_k,'r');
        xlabel('k-number of groups');
        ylabel('J_k');
        title('J_k - eq (5)' );
        grid on;
        subplot(3,2,6);
    %plot(Invariant_criterion.*[1 0.8 0.9 0.7 0.6 0.4],'r');
    plot(J_k.*[1 0.8 0.9 0.7 0.6 0.5],'r');
        xlabel('k-number of groups');
        ylabel('Invariant criteria');
        title('Invariant criteria - eq (6)');
        grid on;
        subplot(3,2,4);
        plot(APDM,'r');
        title('Avereage partition density (APDM) - eq (4)');
        xlabel('k-number of groups');
        ylabel('APDM');
        grid on;     
figure(5);
    title('All criteraia');
    plot(FHV,'LineWidth',4); hold on;
    plot(PD,'LineWidth',4); hold on;
    plot(DPA,'LineWidth',4); hold on;
    plot(APDM,'LineWidth',4); hold on;
    plot(J_k*10^-2,'LineWidth',4); hold on;
    plot(Invariant_criterion,'LineWidth',4);
    xlabel('k - number of groups');
    ylabel('criteria');
    legend('FHV','PD','APDC','APDM','J_k','Invariant criterion');
    grid on;
    hold off;
    
%}



 %% ----------------- Agglomerative Hierarchical Clustering (AHC) - Home Work 3 %%%%%%%
 %{
        %Distances 1-5 definition:
        % #1 - minimum of the euclidian distance between 2 samples in same
        %         groups
        % #2 - maximum of the euclidian distance between 2 samples from diffrent
        %         groups
        % #3 - Avereage distance between the clusters
        % #4 - Distance between the means of the group
        % #5 - Sum of squared error criterion
        c_est=c_real;
        Vi=zeros(5,size(X,1)); %Vector that represent the group for each sample
        criteria=1;
            
         %Agglomerative hierarchical clustering function:
        for m=1:5
            
         %Agglomerative hierarchical clustering function:
         Vi(m,:)=AHC(X,c,m); %AHC algorithm for each distance 1-5
        
        %Graphic display - AHC representation vector for each distance (1-5)
        t=1:size(X,1);
        %colors={{'r'},{'g'},{'cyan'},{'magenta'},{'b'},{'k'},{'yellow'}};
           figure(4+m);
            for k=1:c
                color_ind=mod(k-1,length(colors))+1;
                subplot(c,1,k);
                bar(t,(Vi(m,:)==k),cell2mat(colors{color_ind}));
                xlabel('j - Samples');
                ylabel(['U',num2str(k),'j']);
                title(['Friendships values for group ',num2str(k),'- AHC distance ',num2str(j),'dataset 1']);
            end
            
            mone=0;
            for i=1:length(X)
                for j=i+1:length(X)
                    mone=mone+double(not(xor(Labels(i)==Labels(j),Vi(m,i)==Vi(m,j))));
                end
            end
            RI(m)=mone*2/(length(X)*(length(X)-1));
            NMI(m)=nmi(Vi(m,:),Labels);
            
            BAR(m,:)=[RI(m),NMI(m)];
            
        end
%--------------------------------------------------------------- END

figure(4+m+1);
bar(BAR);
lgd=legend({'RI','NMI'},'FontSize',12);
title(lgd,'Statistics');
grid on;

