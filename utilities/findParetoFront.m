function idx = findParetoFront(X)

    d = length(X(1,:));
    N = length(X(:,1));
    
    idx = false(N,1);
    temp = false(N,d);
    for i=1:N
        
        temp(:) = false;
        for j=1:d
            temp(:,j) = X(i,j) < X(:,j);
        end
                
        if max(sum(temp,2)) < d
            idx(i) = true;
        end
        
    end