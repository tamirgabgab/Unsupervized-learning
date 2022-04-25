function [ Vi ] = AHC( X,c,choice )
%Agglomerative hierarchical clustering

%   -------inputs: 
%   X- dataset size NxD
%   c- number of groups (c<N)
%   choice- choice of the distance (max,min,avg,mean,squared error)

%   -------outputs:
%   Vi-clustering vector  -size 1xN - each index i connected to data sample i
% and the value in the vector tell us to which group the sample belong

N=size(X,1); %number of data samples
Vi=1:N; %Clusteing vector 
k=N; %Number of groups


while(k>c)
min_dist=inf; %Initializing minimum distance
        for ii=1:k-1 %for each 
            for jj=(ii+1):k
            %calculate distance by requested choise
            dist=dist_meas(X,Vi,ii,jj,choice); %Measuring the distance between the groups
                if (dist<min_dist) %Finding the minimum distance
                min_dist=dist; 
                min_i=ii; %Number of group i
                min_j=jj; %Number of group j
                end
            end 
        end
    %Merge between groups
    Vi(Vi==min_j)=min_i;
    
    if(min_j<k)
       Vi(Vi>min_j)=Vi(Vi>min_j)-1; 
    end
    
    k=k-1;
    fprintf('progress : %d \n',k);
end



end

