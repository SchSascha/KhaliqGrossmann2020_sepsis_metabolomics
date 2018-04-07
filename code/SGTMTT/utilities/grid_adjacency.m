function A=grid_adjacency(n,m)

A=zeros(n,m);

for i=1:n
    for j=1:m
%         u=[i-1,j];
%         d=[i+1,j];
%         l=[i,j-1];
%         r=[i,j+1];
        
        if i > 1
            A(n*(i-1)+j, n*(i-2)+j)=1;
        end
        
        if i < n-1
            A(n*(i-1)+j, n*i+j)=1;
        end
        
        if j > 1
            A(n*(i-1)+j, n*(i-1)+j-1)=1;
        end
        
        if j < m-1
            A(n*(i-1)+j, n*(i-1)+j+1)=1;
        end
    end
end