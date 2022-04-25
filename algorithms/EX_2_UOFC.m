

clc;
close all;
clear;
clear ylabel

path='C:\Users\tamir\Desktop\תמיר\matlab scrifts\datasets\';
name='gaos_d3_c4_n10000';
DATASET=load([path name]);
DATASET=DATASET.data;
X=DATASET.X;
Labels=DATASET.Labels';
n=DATASET.n;
d=DATASET.d;
c=DATASET.c;
c_real=c;

bins=100;
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

colors={{[0 0.4470 0.7410]},{[0.8500 0.3250 0.0980]},{[0.4940 0.1840 0.5560]},...
        {[0.4660 0.6740 0.1880]},{[0.3010 0.7450 0.9330]},{[0.6350 0.0780 0.1840]},...
        {'k'},{'r'},{'b'},('m')};

c_real=max(Labels);
u_real=zeros(c_real,length(X));
for i=1:length(X)
    u_real(Labels(i),i)=1;
end


mu_correct=zeros(c_real,d);
sigma_correct=zeros(d,d,c_real);
for w=1:c_real
    index=find(u_real(w,:)==1);
    X_cluster=X(index,:);
    mu_correct(w,:)=round(mean(X_cluster));
    sigma_correct(:,:,w)=round(cov(X_cluster));
    hold on;  
end



figure(1);
for w=1:c_real
    pp=find(u_real(w,:)==1);
    if d==2
        plot(X(pp,1),X(pp,2),'o','color',cell2mat(colors{w}),...
                'LineWidth',1,'MarkerSize',4);
    else
        plot3(X(pp,1),X(pp,2),X(pp,3),'o','color',cell2mat(colors{w}),...
                'LineWidth',1,'MarkerSize',4);
    end
    hold on;  
end
title('Real clusters');
grid on;
xlabel('x');
ylabel('y');
zlabel('z');


figure();
t=1:size(X,1);
title('Real bars');
for w=1:c_real
    subplot(c_real,1,w);
    bar(t,u_real(w,:),'FaceColor',cell2mat(colors{w}));
    xlabel('j - Samples');
    ylabel(['U',num2str(w),'j']);
    ylim([0 1]);
end


    c_max=c_real+2; 
    q=2;
    [u,mu_est,FHV,DPA,PD,J_k,Invariant_criterion,APDM,y_est]=UOFC(X,c_max,q); 
  
    for w=1:c_max
        %{
        figure();
        for k=1:w
            pp=find(y_est(1,:,w)==k);
            if d==2
                plot(X(pp,1),X(pp,2),'o','color',cell2mat(colors{k}),...
                        'LineWidth',1,'MarkerSize',4); hold on;
                plot(mu_est(k,1),mu_est(k,2),'xk',...
                    'LineWidth',5,'MarkerSize',15);
            else
                plot3(X(pp,1),X(pp,2),X(pp,3),'o','color',cell2mat(colors{k}),...
                        'LineWidth',1,'MarkerSize',4); hold on;
                plot3(mu_est(k,1),mu_est(k,2),mu_est(k,3),'xk',...
                    'LineWidth',5,'MarkerSize',15);
            end
            hold on;
        end
        grid on;
        title(['Estimated Clusters (UOFC) for c=',num2str(w)]);
        xlabel('x');
        ylabel('y');
        zlabel('z');
        %}
        figure();
        u_est=zeros(w,size(X,1));
        for m=1:size(X,1)
            u_est(y_est(1,m,w),m)=1;
        end 
        for k=1:w                                          
            color_ind=mod(k-1,length(colors))+1;
            subplot(w,1,k);
            bar(t,u_est(k,:),'FaceColor',cell2mat(colors{k}));
            xlabel('j - Samples');
            ylabel(['U',num2str(k),'j']);
            title(['Friendships values for group - dataset',num2str(k),'- UOFC']);
        end
        
        sill(:,w)=silhouette(X,y_est(1,:,w));
        BAR_sill(w,:)=[median(sill(:,w)),mean(sill(:,w)),...
                       var(sill(:,w))];        
        RI(w)=rand_index(y_est(1,:,w),Labels);
        NMI(w)=nmi(y_est(1,:,w),Labels);
        
        
    end   
        

    figure();
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
        h8=plot(J_k,'r');
        xlabel('k-number of groups');
        ylabel('J_k');
        title('J_k - eq (5)' );
        grid on;
        cla(h8);
        subplot(3,2,6);
    %Invariant_criterion=Invariant_criterion(end:-1:1);
    plot(Invariant_criterion,'r');
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
lg={};
figure();
sill(:,1)=zeros(size(X,1),1);
for w=1:c_max
    subplot(c_max,1,w);
        histogram(sill(:,w),bins,'FaceColor',cell2mat(colors{w}));
        title(['Silhouette cofficiants for c= ',num2str(w)]);
        xlabel('Cofficiant value');
        ylabel('N_{s}');
        xlim([-1 1]);
        grid on;
        lg{w}=['c=',num2str(w)];
end       
figure();
bar(BAR_sill);
lgd=legend({'med','\mu','\sigma^{2}'},'FontSize',12);
title(lgd,'Statistics');
title('Silhouette statistics : S_{i}');
xticklabels(lg)
grid on;

