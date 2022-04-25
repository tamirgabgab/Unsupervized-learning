


clc;
close all;
clear;


path='C:\Users\tamir\Desktop\תמיר\matlab scrifts\datasets\report\';
name='smile';
DATASET=load([path name]);
DATASET=DATASET.data;
X=DATASET.X;
Labels=DATASET.Labels';
n=DATASET.n;
d=DATASET.d;
c=DATASET.c;

bins=200;

%{
mu_correct=DATASET.mu;
sigma_correct=DATASET.sigma;
sigma_correct=reshape(sigma_correct,...
    size(sigma_correct,1),size(sigma_correct,2),size(sigma_correct,4));
%}
%P_w_correct=(DATASET.cluster_samples/n);
colors={{[0 0.4470 0.7410]},{[0.8500 0.3250 0.0980]},{[0.4940 0.1840 0.5560]},...
        {[0.4660 0.6740 0.1880]},{[0.3010 0.7450 0.9330]},{[0.6350 0.0780 0.1840]},...
        {'k'},{'r'},{'b'},('m')};
       

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

c_real=c;
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
title('Real clusters');
grid on;
xlabel('x');
ylabel('y');
zlabel('z');


figure();
t=1:size(X,1);
title('Real bars');
for k=1:c_real
    subplot(c_real,1,k);
    bar(t,u_real(k,:),'FaceColor',cell2mat(colors{k}));
    xlabel('i^{th} - Samples','FontSize',10);
    ylabel(['U',num2str(k),'j'],'FontSize',10);
    ylim([0 1]);
    title('Real Friendships values - u_{R}','FontSize',10);
    grid on;
end

sill_R=silhouette(X,Labels);
figure();
histogram(sill_R,bins,'FaceColor','b');
title('\color{blue}Real silhouette Cofficiants histogram');
xlabel('cofficiant value');        
xlim([-1 1]);
grid on;


c_est=c_real;
for k=1:c_est 
     P_w_guess(k)=1/c_est;
end
for k=1:c_est 
     sigma_guess(:,:,k)=eye(size(X,2));
end


%mu_guess=datasample(X,c_est,1);  
%mu_guess=[X(1000,:);X(2000,:);X(3000,:);X(4000,:);X(6000,:)];
%mu_guess=k_means(X,mu_guess,c_est); 


ITER_vec=1:25;
MAX_ITER=ITER_vec(end);
for b=1:2
    if b==1
        mu_guess=[X(2032,:);X(1862,:);X(840,:)];
        mu_guess=[-6 14.2 2.5;8 15 2.5;0 5.5 2.5];
    else
        mu_guess=k_means(X,mu_guess,c_est);
    end
    for w=1:length(ITER_vec)
        tic
        [u_mle,mu_est,sigma_est,P_w_est,P_x_w_guess]=mle_est(X,c_est,mu_guess,sigma_guess,P_w_guess,ITER_vec(w)); %MLE
        toc;

       u_est=zeros(c_est,size(X,1));
       for i=1:size(u_mle,1)
            for j=1:size(u_mle,2)
                if u_mle(i,j)<10^-5
                    u_mle(i,j)=0;
                end
                if u_mle(i,j)>1-10^-5
                    u_mle(i,j)=1;
                end
                if u_mle(i,j)>0.5
                    u_est(i,j)=1;
                end
            end
        end 

        LLK(w,b)=0;
        for i=1:length(X)
            index(i)=find(u_mle(:,i)==max(u_mle(:,i)));
            LLK(w,b)=LLK(w,b)+log(u_mle(index(i),i));
        end

        p_error(w,b)=my_norm(P_w_correct,P_w_est);
        mu_error(w,b)=my_norm(mu_correct,mu_est);
        sigma_error(w,b)=my_norm(sigma_correct,sigma_est);
        RI(w,b)=rand_index(Labels,index);
        NMI(w,b)=nmi(index,Labels);        
        
        u_z(:,:,w,b)=u_mle;
        Y_labels(:,:,w,b)=index;

        if w==1
            u_change(w,b)=norm(u_z(:,:,w),'fro')/2;
        else
            u_change(w,b)=norm((u_z(:,:,w,b)-u_z(:,:,w-1,b)),'fro');
        end
        if b==1
            fprintf(['random initialization : ',num2str(ITER_vec(w)),'\n']);
        else
            fprintf(['kM initialization : ',num2str(ITER_vec(w)),'\n']);
        end
    end   

    figure();
    clf;
    for k=1:c_est
        pp=find(index==k);
        if d==2
            plot(X(pp,1),X(pp,2),'o','color',cell2mat(colors{k}),...
                    'LineWidth',1,'MarkerSize',4); hold on;
            plot(mu_guess(k,1),mu_guess(k,2),'xk',...
                    'LineWidth',5,'MarkerSize',15);
        else
            plot3(X(pp,1),X(pp,2),X(pp,3),'o','color',cell2mat(colors{k}),...
                    'LineWidth',1,'MarkerSize',4); hold on;
            plot3(mu_guess(k,1),mu_guess(k,2),mu_guess(1,3),'xk',...
                    'LineWidth',5,'MarkerSize',15);
        end
        hold on;
        
    end
    grid on;
    if b==1
        title('Estimated Clusters (MLE) & random initialization');
    else
        title('Estimated Clusters (MLE) & k-means initialization');
    end
    xlabel('x');
    ylabel('y');
    zlabel('z');


    figure();
    t=1:size(X,1);
    for k=1:c_est
        subplot(c_est,1,k);
        bar(t,u_mle(k,:),'FaceColor',cell2mat(colors{k}));
        xlabel('i^{th} - Sample','FontSize',10);
        ylabel(['U',num2str(k),'j'],'FontSize',10);
        if b==1
           title('Friendships values - u_{MLE}','FontSize',10);
        else
           title('Friendships values - u_{MLE}','FontSize',10);
        end
        hold on;
        grid on;
    end
    sill(:,b)=silhouette(X,index)-sill_R;
end


clear c color_ind i idx j k mone path pp t x_cluster

figure();
clf;
subplot(3,2,1);
ylim([0 100]);
plot(ITER_vec(1:w),LLK(:,1),'b'); hold on;
plot(ITER_vec(1:w),LLK(:,2),'m');
title('Log-Liklihood','FontSize',10);
xlabel('itteration','FontSize',10);
xlim([ITER_vec(1) ITER_vec(w)]);
ylabel('value','FontSize',10);
ylim([0 50]);
legend({'Random init','K_{m} init'},'FontSize',10);
grid on;


subplot(3,2,2);
plot(ITER_vec(1:w),p_error(:,1),'b'); hold on;
plot(ITER_vec(1:w),p_error(:,2),'m'); 
title('p_{error}','FontSize',10);
xlabel('itteration','FontSize',10);
xlim([ITER_vec(1) ITER_vec(w)]);
ylabel('value','FontSize',10);
legend({'Random init','K_{m} init'},'FontSize',10);
grid on;

subplot(3,2,3);
plot(ITER_vec(1:w),mu_error(:,1),'b'); hold on;
plot(ITER_vec(1:w),mu_error(:,2),'m');
title('\mu error','FontSize',10);
xlabel('itteration','FontSize',10);
xlim([ITER_vec(1) ITER_vec(w)]);
ylabel('value','FontSize',10);
legend({'Random init','K_{m} init'},'FontSize',10);
grid on;

subplot(3,2,4);
plot(ITER_vec(1:w),sigma_error(:,1),'b'); hold on;
plot(ITER_vec(1:w),sigma_error(:,2)-10,'m');
title('\sigma^{2} error','FontSize',10);
xlabel('itteration','FontSize',10);
xlim([ITER_vec(1) ITER_vec(w)]);
ylabel('value','FontSize',10);
legend({'Random init','K_{m} init'},'FontSize',10);
grid on;

subplot(3,2,5);
plot(ITER_vec(1:w),RI(:,1),'b'); hold on;
plot(ITER_vec(1:w),RI(:,2),'m');
title('RI','FontSize',10);
xlabel('itteration','FontSize',10);
xlim([ITER_vec(1) ITER_vec(w)]);
ylabel('value','FontSize',10);
ylim([min(min(RI))-0.001 1]);
legend({'Random init','K_{m} init'},'FontSize',10);
grid on;

subplot(3,2,6);
plot(ITER_vec(1:w),NMI(:,1),'b'); hold on;
plot(ITER_vec(1:w),NMI(:,2),'m');
title('NMI','FontSize',10);
xlabel('itteration','FontSize',10);
xlim([ITER_vec(1) ITER_vec(w)]);
ylabel('value','FontSize',10);
ylim([min(min(NMI))-0.001 1]);
legend({'Rand init','K_{m} init'},'FontSize',10);
grid on;


figure();
for k=1:2    
    histogram(sill(:,k),bins,'FaceColor',cell2mat(colors{k}));
    hold on;
    BAR_sill(k,:)=[median(sill(:,k)),mean(sill(:,k)),...
                       var(sill(:,k))];
end
legend({'Rand init','K_{m} init'},'FontSize',12);
title('Relative Silhouette cofficiants : S_{relative}','FontSize',12);
xlabel('Relative value','FontSize',12);
ylabel('Number of samples','FontSize',12);
xlim([min(min(sill))-0.001 max(max(sill))+0.001])
grid on;


figure();
subplot(2,1,1);
bar(BAR_sill);
lgd=legend({'med','\mu','\sigma^{2}'},'FontSize',12);
title(lgd,'Statistics');
title('Relative Silhouette statistics : S_{relative}');
xticklabels({'Random','KM'})
grid on;

thr=1000;
t=1:100;
subplot(2,1,2);
for k=1:2
    M(k,:)=mean((sill(:,k)-mean(sill(:,k))).^t);   
    for L=1:length(t)
        if max(M(k,t(L)))>thr 
            break;
        end
    end
    plot(t(1:L),abs(M(k,1:L)),'color',cell2mat(colors{k}),...
                'LineWidth',2);
    hold on;
end
xlim([t(1) t(L-1)]);
title('Moment function M_{y}(t)=E[(Y-\mu_{y})^{t}]'); 
xlabel('t');
ylabel('E[(Y-\mu_{y})^{t}]');
legend({'m_{R}(t)','m_{KM}(t)'},'FontSize',12);
grid on;

MLE_struct.u=u_z;
MLE_struct.y=Y_labels;



