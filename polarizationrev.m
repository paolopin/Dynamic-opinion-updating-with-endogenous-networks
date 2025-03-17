function p=polarizationrev(M,alpha,t,K,n,s)

low=linspace(0,1,s+1);
v=M(t,:);
for i=1:s
    b(i,1)=numel(v(v>=low(i)&v<low(i+1)))/n;
    d(i,1)=mean(v(v>=low(i)&v<low(i+1)));
end
d=fillmissing(d,'constant',0);
for i=1:s
    for j=1:s
        term(i,j)=b(i)^(1+alpha)*(b(j)*abs(d(i)-d(j)));
    end
end
for i=1:s
    somma(i)=sum(term(i,:));
end
p=K*(sum(somma));
%checked
return