function  [ P_w_x_guess,mu_est,sigma_est,P_w_est,P_x_w_guess ] = mle_est( X,c,mu_guess,sigma_guess,P_w_guess,ITER  )
%target - recive c,x and guessed mu,sigma,Pwx , case 2 in the book
%MLE Estimation of mean, variance and group probability

%------------input:
%c-groups
%X-dataset
%mu_est,sigma_est,pw_est - first guess
%N_iterations - number of iterations
s=size(sigma_guess,1)*size(sigma_guess,2)*size(sigma_guess,3);
%-----------output:
%P_w_x_guess - friendships
%MLE estimation of mu,sigma,Pw
    if norm(double(isnan(X)))>0 || ...
       norm(double(isnan(mu_guess)))>0 || ...
       norm(double(isnan(reshape(sigma_guess,s,1,1))))>0
           disp('gfdsgfdsgs');
    end

%At first , estimating probability of group W depend on data X
    %create probability function for each cluster, base on mu,sig
    [P_w_x_guess,P_x_w_guess]=gauss_prob_w_x(X,c,mu_guess,sigma_guess,P_w_guess);
    
    
    %------------Initializing parameters
    tol=10^-10;
    Max_iterations=ITER;
    P_w_next=zeros(size(P_w_guess)); %Next calculation for groups probability
    mu_next=mu_guess;
    sigma_next=sigma_guess;
    mu_prev=zeros(size(mu_guess)); %Next calculation for mean
    sigma_prev=zeros(size(sigma_guess)); %Next calculation for variance
    eps_sigma=(10^-8)*eye(size(X,2));
    eror_tot=10;
    error=ones(1,c);
    N=size(X,1); %Number of samples
    iterations=0;
    
while (iterations<Max_iterations && norm(error)>tol)
     %for each cluster
     for k=1:c
           %save the previous mu,sig
            mu_prev(k,:)=mu_next(k,:);
            sigma_prev(:,:,k)=sigma_next(:,:,k);
            
            %calculate next mu,sig
            %eq 19 - page 19
            P_w_next(k)=mean(P_w_x_guess(k,:));  %Next iteration for measuring the group probability
            %eq 20
            mu_next(k,:)=P_w_x_guess(k,:)*X/sum(P_w_x_guess(k,:)); %Next iteration for measuring the mean

            %repat - duplicate data for eq 21
            X_centered=X-repmat(mu_next(k,:),N,1); %Centering the data ==  (X-mu)(X-Mu)'
            %eq 21 - create the sigma
            sigma_next(:,:,k)=eps_sigma+X_centered'*diag(P_w_x_guess(k,:))*X_centered/sum(P_w_x_guess(k,:)); %Next iteration for measuring the variance
            
            error(k)=norm(mu_prev(k,:)-mu_next(k,:))+norm(sigma_prev(:,:,k)-sigma_next(:,:,k),'fro');
     end %End for k=1:c loop
       
         %Generate new probability fucnction by the new data
         P_w_x_guess=gauss_prob_w_x( X,c,mu_next,sigma_next,P_w_next );
         %fprintf('progress : %d \n',iterations+1);
         iterations=iterations+1;         
end %While \ covergence loop
    

%generate data as function output
        mu_est=mu_next;
        sigma_est=sigma_next;
        P_w_est=P_w_next;
        if(iterations==Max_iterations) %check if stoped by max iteration
        %    disp('max iterations');
        end
        %disp(iterations)
        
end
    
    




