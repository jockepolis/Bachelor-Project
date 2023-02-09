function [val,err] = QMCint(f,D,N,M)
    
    V   = cumprod(D(:,2)-D(:,1));
    dim = numel(V);
    V   = V(end);
    p = sobolset(dim);
    net(p,N);
    r=zeros(dim,N);
    
    for j=1:M
        r(:,:) = repmat(D(:,1),1,N)+net(p,N).*repmat(D(:,2)-D(:,1),1,N);
        I(j)=V*mean(f(r));
    end
    
    val = mean(I);
    err = std(I);
end
