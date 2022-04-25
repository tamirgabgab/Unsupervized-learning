function [u,mu_next,dist]=fuzzy_k_means(X,mu_guess,q,dist,C,MAX)
%target - recive X,mu_guess,q and c
%Fuzzy K-Means algorithm for estimating the centers of the data in a soft
%clustering method

%------------------inputs:
%c- number of groups
%q-fuzziness coeficient
%dist - distance between centers and the data points
%X - dataset
%mu_guess - first estimations of the mean


%-------------------output:
%mu_next - estimations of the means
%u - friendships
%dist - euclidian distances between data points and centers 

%Initialize parameters:
N=size(X,1);%number of samples
mu_next=mu_guess; %Next iteration for finding center
% % dist_mu=ones(1,C); %distance between the center in the present and previous iteration
tol=10^-10;
%%kld=ones(1,N); %Kullback Divergence Leibler 
error=1;
iterations=0; 
max_iterations=MAX; %Maximum iterations allowed

u=ones(C,N)*(1/C); %Friendships - like p_w_x 
%while(error>tol && iterations<max_iterations)
while(iterations<max_iterations) % Error = Maximum KLD + maximum distance of prob
% %     mu_prev=mu_next;
    u_prev=u;

    %Calculating probabilities
    for s=1:C
        for k=1:N
            %set friendship for each data eq 28
          u(s,k)=(dist(s,k).^(-2/(q-1)))./(sum(dist(:,k).^(-2/(q-1))));
        end
    end
    
    %Calculating new centers
    for s=1:C
        %eq 27
        %(Friendship(c).* data )/(Friendship(c))
       mu_next(s,:)=(u(s,:).^q)*X/(sum(u(s,:).^q));
    end
    iterations=iterations+1;
    fprintf('fuzzy itterations : %d \n',iterations);
    
    %Measuring the distances between every point to a center
    for s=1:C
        for k=1:N
           dist(s,k)=norm(X(k,:)-mu_next(s,:))+10^(-9);
        end
    end
    error=max(max(abs(u-u_prev)));
end

if(iterations==max_iterations)
    disp('exceed number of iterations');
end

end