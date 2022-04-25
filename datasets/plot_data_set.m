
function v=plot_data_set(string)
    
    DATASET=load(string);
    DATASET=DATASET.data;
    name=DATASET.name;
    X=DATASET.X;
    Labels=DATASET.Labels;
    n=DATASET.n;
    d=DATASET.d;
    c=DATASET.c;

    cluster_colors=rand(c,3);
    if d==2
        plot(min(X(:,1)),min(X(:,2)),'LineWidth',0.1); hold on;
        plot(max(X(:,1)),max(X(:,2)),'LineWidth',0.1); hold on;
    elseif d>=3
        plot3(min(X(:,1)),min(X(:,2)),min(X(:,3)),'LineWidth',0.1); hold on;
        plot3(max(X(:,1)),max(X(:,2)),max(X(:,3)),'LineWidth',0.1); hold on;
    end
    grid on;
    pause(1);
    for i=1:n
        pause(0.5/n);
        if Labels(i,1)<0 && d==2
            pause(1/n);
            plot(X(i,1),X(i,2),'*r',...
            'LineWidth',1,'MarkerSize',3);
        elseif Labels(i,1)<0 && d>=3
            pause(1/n);
            plot3(X(i,1),X(i,2),X(i,3),'*r',...
            'LineWidth',1,'MarkerSize',3); 
        elseif Labels(i,1)>=0 && d==2
             plot(X(i,1),X(i,2),'o','color',cluster_colors(Labels(i),:),...
            'LineWidth',1,'MarkerSize',4);            
        elseif Labels(i,1)>=0 && d>=3
             plot3(X(i,1),X(i,2),X(i,3),'o','color',cluster_colors(Labels(i),:),...
            'LineWidth',1,'MarkerSize',4);
        end
        hold on;
    end
    grid on;
    str=strcat('Data set : ',name);
    for i=1:length(str)
        pause(0.1);
        title(strcat(str(1:i)));
    end
    pause(1);
    v='Complete!!!';
end