
function plot_dataset(DATASET)

    DATASET=DATASET.data;
    X=DATASET.X;
    Labels=DATASET.Labels;
    n=DATASET.n;
    d=DATASET.d;
    c=DATASET.c;
    bins=100;

    colors={{[1 0 0]},{'k'},{'b'},{'magenta'},{'green'},...
        {[0 0.4470 0.7410]},{[0.8500 0.3250 0.0980]}...
        ,{[0.9290 0.6940 0.1250]},{[0.4940 0.1840 0.5560]}...
        ,{[0.4660 0.6740 0.1880]},{[0.3010 0.7450 0.9330]},....
        {[0.6350 0.0780 0.1840]},{[0.8790 0.2640 0.5450]}};

    figure();    
    for k=1:c
        pp=find(Labels==k);
        if d==2
            plot(X(pp,1),X(pp,2),'o','color',cell2mat(colors{k}),...
                    'LineWidth',1,'MarkerSize',4);
        else
            plot3(X(pp,1),X(pp,2),X(pp,3),'o','color',cell2mat(colors{k}),...
                    'LineWidth',1,'MarkerSize',4);
        end
        hold on;  
    end
    title('\color{blue}Clustering plot');
    grid on;

    u=zeros(c,n);
    for i=1:n
        u(Labels(i),i)=1;
    end
    
    figure();   
    t=1:n;
    for k=1:c
        subplot(c,1,k);
        bar(t,u(k,:),'FaceColor',cell2mat(colors{k}));
        xlabel('j^{th} - Samples');
        ylabel(['U',num2str(k),'j']);
        ylim([0 1]);
        grid on;
    end
    suptitle('\color{blue}friendship values');


    figure();
    histogram(silhouette(X,Labels),bins,'FaceColor','b');
    title('\color{blue}Silhouette Cofficiants histogram');
    xlabel('cofficiant value');        
    xlim([-1 1]);
    grid on;

end