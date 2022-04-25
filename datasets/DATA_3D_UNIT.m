




close all
clear all

%Gaos shape or cube
type='cube';   
d=3;
c=6;
c_real=c;
c_est=c;

cluster_samples=1000*ones(1,c);
rotation=[0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0];
colors={{'r'},{'k'},{[0.9290 0.6940 0.1250]},{'magenta'},{'b'},{'g'},{'yellow'},...
        {[0.3010 0.7450 0.9330]},{[0.6350 0.0780 0.1840]}};
sigma=[1 2 4];
alpha=0.99919224588;
r=0:0.5:5;
ITTER=20;
c_est=c_real;
bins=200;
t=1:500;
thr=500;

for k=1:c_est 
     P_w_guess(k)=1/c_est;
end
for k=1:c_est 
    sigma_guess(:,:,k)=eye(d);
end

for s=1:length(sigma)
    for j=1:length(r)
        %r=30;                                           %distance from center
        centers=r(j)*[1 0 0;-1 0 0;0 1 0;0 -1 0;0 0 1;0 0 -1];
        mu_correct=centers;
        sigma_correct=zeros(d,d,c);
        for i=1:c
            sigma_correct(:,:,i)=eye(d)*sigma(s)^2;
        end
                                                       %ball in radius 3*sigma cover aplha from cube in len = L
        X=[];
        Labels=[];       
        for k=1:c
            R_A=[1 0 0;...
                0 cos(rotation(k,1)) sin(rotation(k,1));...
                0 -sin(rotation(k,1)) cos(rotation(k,1))];
            R_B=[cos(rotation(k,2)) 0 -sin(rotation(k,2));...
                 0 1 0;...
                 sin(rotation(k,2)) 0 cos(rotation(k,2))];
            R_C=[cos(rotation(k,3)) sin(rotation(k,3)) 0;...
                 -sin(rotation(k,3)) cos(rotation(k,3)) 0;...
                 0 0 1];
            if strcmp(type,'cube')
                L=sigma(s)*(36*pi/alpha)^(1/3);
                len=L*ones(c,d);
                ADD=centers(k,:)+len(k,:).*(rand(cluster_samples(k),d)-0.5);
            else       
                ADD=mvnrnd(mu_correct(k,:),sigma_correct(:,:,k),cluster_samples(k));
            end
            X=cat(1,X,ADD*R_A*R_B*R_C);
            Labels=cat(1,Labels,k*ones(cluster_samples(k),1));            
            plot3(ADD(:,1),ADD(:,2),ADD(:,3),'o','color',cell2mat(colors{k}),...
                        'LineWidth',1,'MarkerSize',4);
            title(['Real cluster (',type,') and r=',...
                  num2str(r(j)),' \sigma^{2}=',num2str(sigma(s))]);
            xlabel('x'); xlim(r(end)*[-1 1]);
            ylabel('y'); ylim(r(end)*[-1 1]);
            zlabel('z'); zlim(r(end)*[-1 1]);            
            hold on;
            
        end
        grid on;
              
        sill_R=silhouette(X,Labels);
        %{
        figure();
        histogram(sill_R,bins,'FaceColor','b');
        title('\color{blue}Real silhouette Cofficiants histogram');
        xlabel('cofficiant value');        
        xlim([-1 1]);
        grid on;
        %}
        
        for b=1:2
            if b==1
                ran=randi([1 size(X,1)],1,c_est);
                ran=[1 5 7 18 24 33];
                mu_guess=X(ran,:);
            else
                mu_guess=k_means(X,mu_guess,c_est);
            end
            
            [u_mle,mu_est,sigma_est,P_w_est,P_x_w_guess]=...
            mle_est(X,c_est,mu_guess,sigma_guess,P_w_guess,ITTER);
            disp(mu_est);
            if norm(double(isnan(u_mle)))>0 || ...
               norm(double(isnan(mu_est)))>0 || ...
               norm(double(isnan(reshape(sigma_est,d*d*c,1,1))))>0
               disp('gfdsgfdsgs');
            end
                                                
            u_est=zeros(c_est,size(X,1));
            for i=1:size(u_mle,1)
                 for v=1:size(u_mle,2)
                     if u_mle(i,v)<10^-5
                         u_mle(i,v)=0;
                     end
                     if u_mle(i,v)>1-10^-5
                         u_mle(i,v)=1;
                     end
                     if u_mle(i,v)>0.5
                         u_est(i,v)=1;
                     end
                 end
            end
            
            for i=1:size(X,1)
                index(i)=find(u_mle(:,i)==max(u_mle(:,i)));
            end
            mu_error(j,s,b)=my_norm(mu_correct,mu_est);
        	sigma_error(j,s,b)=my_norm(sigma_correct,sigma_est);
            RI(j,s,b)=rand_index(Labels,index);
            NMI(j,s,b)=nmi(index,Labels);
            S=silhouette(X,index)-sill_R;
            BAR_sill(j,s,b,1)=median(S);
            BAR_sill(j,s,b,2)=mean(S);
            BAR_sill(j,s,b,3)=var(S);
            M(j,s,b,:)=mean((S-mean(S)).^t);
            %REL_sill(j,s,b,:)=silhouette(X,index)-sill_R;
            fprintf(['r=',num2str(r(j)),...
                '  sigma=',num2str(sigma(s)),'  b=',num2str(b),'\n']);            
        end 
        
    end
end

RES(:,:,:,1)=mu_error;
RES(:,:,:,2)=sigma_error;
RES(:,:,:,3)=RI;
RES(:,:,:,4)=NMI;

tit={'\mu error','\sigma^{2} error','RI','NMI'};
for b=1:2
    figure(b+2);
    if b==1
        suptitle('Random initialization')
    else
        suptitle('k-mean initialization')
    end
    for PLOT=1:4
        subplot(2,2,PLOT)
        for s=1:length(sigma)
            plot(r,RES(:,s,b,PLOT),'color',cell2mat(colors{s}),...
                'LineWidth',2);
            hold on;
            lg{s}=['\sigma^{2}=',num2str(sigma(s))];
        end
        title(tit{PLOT},'FontSize',10);
        xlabel('r','FontSize',10);
        ylabel('Value','FontSize',10);
        legend(lg,'FontSize',10);
        grid on;
    end
end


figure();
LG={};
i=1;
for b=1:2
    for s=1:length(sigma)
        if b==1
            LG{i}=['R - \sigma^{2}=',num2str(sigma(s))];
        else
            LG{i}=['KM - \sigma^{2}=',num2str(sigma(s))];
        end        
        i=i+1;
    end
end

for PLOT=3:4
    subplot(2,1,PLOT-2);
    for s=1:length(sigma)
        plot(r,smooth(RES(:,s,1,PLOT)),'--','color',cell2mat(colors{s}),...
        'LineWidth',1);
        xlabel('r','FontSize',10);
        if PLOT==3
            ylabel('RI','FontSize',10); 
        else
            ylabel('NMI','FontSize',10); 
        end              
        hold on;
    end    
    for s=1:length(sigma)
        plot(r,RES(:,s,2,PLOT),'color',cell2mat(colors{s}),...
        'LineWidth',1.5);
        xlabel('r','FontSize',10);
        if PLOT==3
            ylabel('RI','FontSize',10); 
        else
            ylabel('NMI','FontSize',10); 
        end       
        hold on;
    end
    title([type,' disribution']);
    legend(LG,'FontSize',10);
    grid on;
end



figure();
rad=[3 length(r)];
u=0;
for z=1:length(rad)
    u=u+1;
    subplot(2,2,u);
    A=BAR_sill(rad(z),:,:,:);
    A=reshape(A,length(sigma)*2,3);
    bar(A);
    lgd=legend({'med','\mu','\sigma^{2}'},'FontSize',12);
    title(lgd,'Statistics');
    title(['Relative Silhouette statistics : S_{relative} ans R=',num2str(r(rad(z)))]);
    xticklabels({'R - \sigma^{2}=1 ','KM - \sigma^{2}=1 ',...
                 'R - \sigma^{2}=2 ','KM - \sigma^{2}=2 ',...
                 'R - \sigma^{2}=4','KM - \sigma^{2}=4 '...
                 });
    grid on;
    u=u+1;
    subplot(2,2,u)
    for s=1:length(sigma)
        for b=1:2
            for L=1:length(t)
                if max(M(rad(z),s,b,t(L)))>thr 
                    break;
                end
            end
            U=reshape(M(rad(z),s,b,1:L),1,L);
            if b==1
                plot(t(1:L),abs(U),'--','color',cell2mat(colors{s}),...
                    'LineWidth',2);
            else
                plot(t(1:L),abs(U),'color',cell2mat(colors{s}),...
                    'LineWidth',2);
            end
            hold on;
        end
    end
    xlim([t(1) t(L-1)]);
    title('Moment function M_{y}(t)=E[(Y-\mu_{y})^{t}]'); 
    xlabel('t');
    ylabel('E[(Y-\mu_{y})^{t}]');
    legend('R - \sigma^{2}=1 ','KM - \sigma^{2}=1 ',...
           'R - \sigma^{2}=2 ','KM - \sigma^{2}=2 ',...
           'R - \sigma^{2}=4 ','KM - \sigma^{2}=4 '...
           );
    grid on;
end