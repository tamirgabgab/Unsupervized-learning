



clc;
close all;
clear;


path='C:\Users\tamir\Desktop\תמיר\matlab scrifts\datasets\report\';
name='very_dense_gaos_ahc';
DATASET=load([path name]);
DATASET=DATASET.data;
X=DATASET.X;
Labels=DATASET.Labels';
n=DATASET.n;
d=DATASET.d;
c=DATASET.c;
bins=200;

res=1;
X=X(1:res:end,:);
Labels=Labels(:,1:res:end);

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
sill_R=silhouette(X,Labels);


mu_correct=zeros(c_real,d);
sigma_correct=zeros(d,d,c_real);
for k=1:c_real
    index=find(u_real(k,:)==1);
    X_cluster=X(index,:);
    mu_correct(k,:)=round(mean(X_cluster));
    sigma_correct(:,:,k)=round(cov(X_cluster));
end



figure();
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
title('\color{blue}Real clusters');
xlabel('x');
ylabel('y');
zlabel('z');
grid on;


figure();
t=1:size(X,1);
for k=1:c_real
    subplot(c_real,1,k);
    bar(t,u_real(k,:),'FaceColor',cell2mat(colors{k}));
    xlabel('j - Samples');
    ylabel(['U',num2str(k),'j']);
    ylim([0 1]);
    title('Real Friendships value u_R');
    grid on;
end

        %Distances 1-5 definition:
        % #1 - minimum of the euclidian distance between 2 samples in same
        %         groups
        % #2 - maximum of the euclidian distance between 2 samples from diffrent
        %         groups
        % #3 - Avereage distance between the clusters
        % #4 - Distance between the means of the group
        % #5 - Sum of squared error criterion
        % #6 - The majority from all 5 criterias
        
        c_est=c_real;
        Vi=zeros(6,size(X,1)); 
            
        for m=1:size(Vi,1)
            if  m<=5
                Vi(m,:)=AHC(X,c_est,m);
            else
                for i=1:size(X,1)
                    Vi(m,i)=mode(Vi(1:m-1,i));
                end
            end
               
        t=1:size(X,1);            
            figure(4);
            subplot(2,3,m);
            for k=1:c_est
                pp=find(Vi(m,:)==k);
                if d==2
                    plot(X(pp,1),X(pp,2),'o','color',cell2mat(colors{k}),...
                            'LineWidth',1,'MarkerSize',4);
                else
                    plot3(X(pp,1),X(pp,2),X(pp,3),'o','color',cell2mat(colors{k}),...
                            'LineWidth',1,'MarkerSize',4);
                end
                hold on;  
            end            
            title({['Criteria  ',num2str(m)];['RI=',num2str(rand_index(Vi(m,:),Labels))];...
                               ['NMI=',num2str(nmi(Vi(m,:),Labels))]});
            xlabel('x');
            ylabel('y');
            zlabel('z');
            grid on;
            
            %{
            figure();
            for k=1:c_est
                color_ind=mod(k-1,length(colors))+1;
                subplot(c_est,1,k);
                bar(t,(Vi(m,:)==k),'FaceColor',cell2mat(colors{color_ind}));
                xlabel('i^{th} Samples');
                ylabel(['U_{',num2str(k),'j}']);
                title(['Friendships values U_{AHC} for group ',num2str(k),'- AHC with criteria ',num2str(m)]);
                grid on;
            end
            
            RI(m)=rand_index(Labels,Vi(m,:));           
            NMI(m)=nmi(Vi(m,:),Labels);           
            BAR(m,:)=[RI(m),NMI(m)];           
            %}
        end

figure();
methods=[{'Single'},{'Complete'},{'Average'},{'Centroid'},{'Ward'},{'Weighted'}];
for i=1:length(methods)
    subplot(2,3,i);
    Z=linkage(X,methods{i});
    T(:,i)=cluster(Z,'maxclust',c_est);
    %dendrogram(Z,'ColorThreshold',median([Z(end-c_est,3) Z(end-1,3)]));
    %dendrogram(Z,'ColorThreshold',Z(end-c_est+2,3));
    dendrogram(Z,0,'ColorThreshold',Z(end-c_est+2,3));
    title({[methods{i},' Linkage'];['RI=',num2str(rand_index(T(:,i),Labels))];...
                               ['NMI=',num2str(nmi(T(:,i),Labels))]});
    ylabel('Heigths');
    xticks('');
    grid on;
end           
sill=zeros(size(X,1),1); 
L=zeros(1,5);
for k=1:size(Vi,1)
    sill(:,k)=silhouette(X,Vi(k,:))-sill_R;
    %z=histogram(sill(:,k),bins);
    %L(k)=max(z.Values);
    BAR_sill(k,:)=[median(sill(:,k)),mean(sill(:,k)),...
                   var(sill(:,k)),my_entropy(sill(:,k))];
end
%L=max(L);        
%{       
figure();
for k=1:size(Vi,1)            
        subplot(size(Vi,1),1,k);
        plot(sill(:,k),'color',cell2mat(colors{k}));
        title(['Rlative Silhouette cofficiants for criteria ',num2str(k)]);
        xlabel('i^{th} smaple');
        ylabel=('cofficiant value');
        xlim([1 size(Vi,2)]);
        ylim([-1 1]);
        grid on;
end
%}

figure();
x_res=4;
for k=1:size(Vi,1)             
        subplot(size(Vi,1),1,k);
        his=histogram(sill(:,k),bins,'FaceColor',cell2mat(colors{k}));
        title(['Relative Silhouette Coefficients for criteria ',num2str(k)]);
        xlabel('cofficiant value');        
        xlim([-1 1]);
        %ylim([0 L+mod(L,x_res)]);
        %yticks(0:round(L/x_res):L+round(L/x_res))
        grid on;      
end

%{
lim=-2:0.05:2;
figure();
for k=1:size(Vi,1)
    title('Fit the cofficiants into normal disribution');
    pd=fitdist(sill(:,k)+0.001*rand(size(sill,1),1),'normal');
    pdf_est=pdf(pd,lim);
    plot(lim,pdf_est,'color',cell2mat(colors{k}),...
                'LineWidth',2);
    hold on;
    plot([pd.mu pd.mu],[0 max(pdf_est)],'--','color',cell2mat(colors{k}),...
                'LineWidth',1);
    hold on;
    xlabel('x');
    
end
xlim([-1 1]);
ylim([0 3]);
lgd=legend({'1','\mu_1','2','\mu_2','3','\mu_3',...
            '4','\mu_4','5','\mu_5','6','\mu_6'},'FontSize',12);
lgd.NumColumns=6;
grid on;
 %}

%{
figure();
bar(BAR);
lgd=legend({'RI','NMI'},'FontSize',12);
title(lgd,'Statistics');
grid on;
%}

figure();
bar(BAR_sill);
lgd=legend({'med_S','\mu_S','\sigma^{2}_S','H_S'},'FontSize',12);
title(lgd,'Statistics');
yticks(-1:0.1:1);
grid on;



