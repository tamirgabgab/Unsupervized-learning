function [ P_w_x,P_x_w ] = gauss_prob_w_x( X,c,mu,sigma,P_w )
%traget - given the mu and sigma the function will give probability of
%group W depand on data X

% input:
%c- number of groups
%X- data
%Mean, sigma and Probability of group

%output:
%P_w_x  - dependent probability

%Conditional probability of x given wi
sample_size=size(X,1);
P_x_w=zeros(c,sample_size);
P_w_x=zeros(c,sample_size);

%for each cluster
    for k=1:c
    %create distribution function from mu and sigma,for the X dataset
        if all(isnan(sigma(:,:,k))>0)
            sigma(:,:,k)=eye(size(X,2));
        end
    gm = gmdistribution(mu(k,:),sigma(:,:,k));
    P_x_w(k,:)=pdf(gm,X);
    end

    
%Calculating conditional probability for wi given xi 
%for range of X 
%P(w|x)=eq 22 - page 19 class book
for l=1:sample_size
          P=0;
          %mone calculation
          for s=1:c %sun for all cluster
             %summerize the Pxw*Pw
             P=P+P_x_w(s,l)*P_w(s);
          end
         
          %mechane/mone
          for k=1:c
            P_w_x(k,l)=P_x_w(k,l)*P_w(k)./P; %Calculating probability for each group
          end 
end
      

end

