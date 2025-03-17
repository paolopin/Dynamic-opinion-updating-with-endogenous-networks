function F = net_sel(a,i,n,f,V)

% clear
% a=[0,0.0344827586206897,0.0689655172413793,0.103448275862069,0.137931034482759,0.172413793103448,0.206896551724138,0.241379310344828,0.275862068965517,0.310344827586207,0.344827586206897,0.379310344827586,0.413793103448276,0.448275862068966,0.482758620689655,0.517241379310345,0.551724137931035,0.586206896551724,0.620689655172414,0.655172413793103,0.689655172413793,0.724137931034483,0.758620689655172,0.793103448275862,0.827586206896552,0.862068965517241,0.896551724137931,0.931034482758621,0.965517241379310,1];
% i=16;
% n=30;
% f=.1;
% V=.002;

% given a vector <a> of types (which is supposed to be increasing), 
% and its i'th element, this algorithm returns the best connections for i, 
% so to maximize the objective function:
% k [ V - f(1 - f)(mu_i - ti)^2 - f(sigma_i)2 ]

% let us start with a naive approach: the node considers all possible
% connected subsets, which are (n-1)*(n-2)

A=zeros(n,n);
    
for j=1:n
    if j~=i
        for k=j:n
            neigh=zeros(1,n);
            neigh(j:k)=neigh(j:k)+1;
            neigh(i)=0;
            degree=sum(neigh);
            m=mean(a(logical(neigh)));
            v=var(a(logical(neigh)));
            
            A(j,k)=degree*(V-f*(1-f)*(m-a(i))^2-f*v);
        end
    end
end

%find the largest element, if it is greater than 0
    
neigh=zeros(1,n);
    
if max(max(A))>0    
    [M,I] = max(A(:));
    [I_row, I_col] = ind2sub(size(A),I);

    neigh(I_row:I_col)=neigh(I_row:I_col)+1;
    neigh(i)=0;
end

F=neigh;





            

