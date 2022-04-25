function [u,mu_est,FHV,DPA,PD,J_k,Invariant_criterion,APDM,y_est] = UOFC( X,c,q )
%target - recive X and c%   Detailed explanation goes here
%UOFC - Unsupervised Optimal Fuzzy Clustering algorithm

%---------------inputs:
%c- number of groups
%X - dataset as NxD matrix

%--------------outputs :
%u-friendships values
%mu_est - estimated means
%FHV, DPA, PD, J_k, Invariant_criterion - clustering validities

k=1; %Numer of groups
N=size(X,1);
D=size(X,2);
mu_est=zeros(c,D); 
%dist=zeros(c,N); %Distance between centers and data points
%dist_euc=zeros(c,N);

m=mean(X); %Total mean vector
mu_est(k,:)=m; %First guessing of center

%Finding the distance between points to the center
for n=1:N
        dist(k,n)=norm(X(n,:)-mu_est(k,:))+10^(-9);
end

%Clustering Vaidity parameters
FHV=zeros(1,c); %Fuzzy hypervolume
DPA=zeros(1,c); %Avereage partition density
APDM=zeros(1,c);
PD=zeros(1,c);  %Partition density
J_k=zeros(1,c); %Normalized by k (number of groups) patition indexes criterion
Invariant_criterion=zeros(1,c); %Invariant criterion

X_centered=X-repmat(mean(X),N,1); %Centering the data

%cacluate the first validity data
FHV(1)=det((X_centered'*X_centered)./N).^0.5;%The Fuzzy hypervolume criterion:
DPA(1)=N./FHV(1);%The avarage partition Denstiy
PD(1)=DPA(1);
J_k(1)=sum(dist(k,:).^2);
ni=zeros(1,c); %Size of group i (1<=i<=c)

u_est(1,:,k)=ones(1,size(X,1));
y_est(1,:,k)=ones(1,size(X,1));

while(k<c) %Updating data condition
 Sw=zeros(D,D); %Within scatter matrix
 k=k+1;
 u_group=zeros(k,N);
 for l=1:k
     u_group(l,1)=10^-15;
 end
 times=zeros(1,k)+10^-8;
 
 %Calculating center far away from all data, as 5 Variance
 dist(k,:)=(5*norm(var(X))).^2*ones(1,N); 
 
 %FKM
 [u,mu_est,dist]=fuzzy_k_means(X,mu_est,q,dist,k,30);

 %Fuzzy MLE - Fix by the Fuzzy Cov matrix and MLE
 [u,mu_est,Fi,~] = fmle( X,mu_est,u,k );
 
 for m=1:size(X,1)
     y_est(1,m,k)=find(u(:,m)==max(u(:,m)));
 end
 
 %Calculating the cluster validity parameters
 for l=1:N %for each data
     for ii=1:k %for each cluster
         %if distance inside elepsoide - eq 14 in the article 
        if (((X(l,:)-mu_est(ii,:))*(Fi(:,:,ii)\((X(l,:)-mu_est(ii,:)))'))<1) 
            u_group(ii,l)=u(ii,l);
            times(ii)=times(ii)+1; %counter
        end       
     end 
 end
 
 %Calculating the m_k and V_APDM
 m_k=zeros(k,1);
 for i=1:N
     MAX=find(u(:,i)==max(u(:,i)));
     m_k(MAX)=m_k(MAX)+max(u(:,i));
 end
 
 for n=1:k
    FHV(k)=FHV(k)+det(Fi(:,:,n)).^0.5; 
    DPA(k)=DPA(k)+sum(nonzeros(u_group(n,:)))./(det(Fi(:,:,n)).^0.5);
    APDM(k)=DPA(k)+sum(m_k(n,1)./(det(Fi(:,:,n)).^0.5));
 end
 DPA(k)=DPA(k)./k;
 APDM(k)=APDM(k)./k;
 PD(k)=sum(sum(nonzeros(u_group)))./FHV(k);
 
 %Normalized by k (number of groups) patition indexes criterion
  for n=1:k
     J_k(k)=J_k(k)+(u(n,:).^q)*(dist(n,:).^2)';
  end
  J_k(k)=J_k(k)*k;
  
  %%% Invariant criterion %%%
  %Calculating size of every group
  for l=1:k
      ni(l)=sum(u(l,:));
  end
  
  %Calculating withing clustering
  for l=1:k
      Sw=Sw+Fi(:,:,l)*ni(l);
  end
  
    %Calculating between clustering
    SB=(mu_est(1:k,1:D)-repmat(m,k,1))'*diag(ni(1:k))*(mu_est(1:k,1:D)-repmat(m,k,1));
    Invariant_criterion(k)=trace(Sw\SB); %Invariant criterion calculation
    fprintf('groups assume : %d \n',k);
end
    

end
