
function [idx,C]=fuzzy_c_means(X,u_0,m,MAXitter,tol)
    if MAXitter==0
        return
    end
    u_prev=u_0;
    itter=0;
    [c,n]=size(u_prev);
    u_next=999*ones(c,n);
    while itter<=MAXitter || max(max(abs(u_next-u_prev)))>tol
        if itter>0
            u_prev=u_next;
        end
        for j=1:c
            C(j,:)=(u_prev(j,:).^m)*X./sum((u_prev(j,:)).^m);
        end
        
        for i=1:n
            s=zeros(1,c);
            for j=1:c
                for k=1:c
                    s(j)=(norm(X(i,:)-C(j,:))/norm(X(i,:)-C(k,:)))^(2/(m-1));
                end
                u_next(j,i)=sum(s)^-1;
            end
        end
        itter=itter+1;
    end
    idx=u_next;
end