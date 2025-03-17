% Krause model
n=101; %number of nodes
k=.22; %distance of interaction (threshold around 0.21-0.211)
T=80; % iterations
in=0.5; %inertia

S=zeros(n,T); %size of peer group

M=zeros(n,T);
for i=1:n
    M(i,1)=(i-1)/(n-1);
end
A=zeros(n,n);
for i=1:n
    for j=1:n
        if i~=j && abs(M(j,1)-M(i,1))<=k
            A(i,j)=1;
        end
    end
end

G=digraph(A);
figure
%plot(G,'Layout','force')
deg_out=outdegree(G);
deg_in=indegree(G);
axisx=linspace(0,n,n);
figure
plot(axisx,deg_in)
diam=diameter(A);

for t=2:T
    for i=1:n
        a=M(:,t-1);
        x=a(i);
        ref=(a>x-k & a<x+k);
        h=sum(ref);
        S(i,t)=h;
        if h>1
            mean=(sum(a(ref))-x ) / (h-1);
        else
            mean=x;
        end
        M(i,t)=in*x+(1-in)*mean;
    end
end

figure
for i=1:n
    plot(M(i,:));
    hold on
end
