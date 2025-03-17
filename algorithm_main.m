%% Matlab code for ``Dynamic opinion updating with endogenous networks''
%% Ugo Bolletta, Paolo Pin
%% European Economic Review (2025)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This program allows to perform numerical simulation for the model, after
% having selected a value (or a range of values) for the parameters V,f,n, 
% and with a uniform distribution of initial opinions. The latter can be
% changed to any type of discrete distribution. 
% The algorithm selects the payoff maximizing network in time t, given a
% distribution of opinions at t-1. Thus, opinions are updated. The outcome
% is a vector of opinions with size nx1. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ADDITIONAL FILES NEEDED:
% - net_sel.m (the algorithm to create the optimal set of links for an 
%   agent)
% - networkComponents.m (used to count the number of components at the end 
%   of the simulation)
% - diameter.m (calculates the diameter of the network - SOURCE: Matlab 
%   Tools for Network Analysis (2006-2011)Copyright (c) 2011,
%   Massachusetts Institute of Technology.)
% - simple_dijkstra.m (used in diameter.m- SOURCE: Matlab Tools for Network
%   Analysis (2006-2011)Copyright (c) 2011, Massachusetts Institute of 
%   Technology.)
% - polarizationrev.m is used to calculate measures of polarization
% - figure6.m to replicate the figure 6 in the paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTES: 
% In the ``Parameters'' section we propose several combination of
% parameters that allow to reproduce the figures in the paper. The
% initially uncommented section is left to the user to put in any
% combination of parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters %%
% The reader can manipulate the main parameters to play
% around with the model. For replication purposes leave this section
% commented and uncomment further sections. NOTE: the parameter ``n'' must
% be an odd number for the simulations. Select one set of parameters at a
% time. 
%--------------------------------------------------------------------------

% Figure 1a --- UNCOMMENT TO REPLICATE%
% rng=1;
% n=101; %number of nodes
% T=25; % iterations
% A=zeros(n,n);
% f=0.5;
% V=0.0350916;
% grid_fine=100;
%--------------------------------------------------------------------------

% Figure 1b --- UNCOMMENT TO REPLICATE%
rng=1;
n=101; %number of nodes
T=50; % iterations
A=zeros(n,n);
f=0.5;
V=0.0350917;
grid_fine=100;
%--------------------------------------------------------------------------

% Figure 4b --- UNCOMMENT TO REPLICATE%
% rng=1;
% n=101; %number of nodes
% T=25; % iterations
% A=zeros(n,n);
% f=0.5;
% V=0.01;
% grid_fine=100;
%--------------------------------------------------------------------------

% Figure 5 --- UNCOMMENT TO REPLICATE%
% rng=1;
% T=25; % iterations
% n=80;
% A=zeros(n,n);
% grid_fine=100;
% f=linspace(0.01,0.95,grid_fine); % flexibility
% V=linspace(0.005,0.1,grid_fine); % benefit
% Pol_1=zeros(length(V)*length(f),T);
% Pol_2=zeros(length(V)*length(f),T);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DISTRIBUTIONS
% Uncomment to produce a Uniform, Normal, or Bimodal distribution for given
% parameters. Select 1 of the following options.
%--------------------------------------------------------------------------

%% Uniform distribution
%--------------------------------------------------------------------------

S=zeros(n,T); %size of peer group
M=zeros(T,n);
C=zeros(n,T);
y=zeros(n,n);
B=zeros(n,T);
Z=zeros(length(f),length(V));
W=zeros(length(f),length(V));
uW=zeros(length(f),length(V));
Y=zeros(length(f),length(V));
grid=zeros(length(f), length(V)); 
for i=1:n
    M(1,i)=(i-1)/(n-1); % Creates a uniform distribution of initial opinions
end
%--------------------------------------------------------------------------

%% Normal Distribution
%--------------------------------------------------------------------------

% seg_a=linspace(0,0.2-0.2/11,11);
% seg_b=linspace(0.2,0.4,20);
% seg_c=linspace(0.4,0.6,42);
% seg_d=linspace(0.6,0.8,20);
% seg_e=linspace(0.8,1,11);
% 
% for i=1:length(seg_a)
% seg_1(i)=seg_a(1,i);
% end
% for i=1:length(seg_b)-1
% seg_2(i)=seg_b(1,i);
% end
% for i=1:length(seg_c)-1
% seg_3(i)=seg_c(1,i);
% end
% for i=1:length(seg_d)-1
% seg_4(i)=seg_d(1,i);
% end
% for i=1:length(seg_e)
% seg_5(i)=seg_e(1,i);
% end
% M=zeros(T,n);
% M(1,:)=[seg_1 seg_2 seg_3 seg_4 seg_5];
%--------------------------------------------------------------------------

%% Bimodal Distribution
%--------------------------------------------------------------------------

% seg_a=linspace(0,0.1-0.005,6);
% seg_b=linspace(0.1,0.2,10);
% seg_c=linspace(0.2,0.3,21);
% seg_d=linspace(0.3,0.4,10);
% seg_e=linspace(0.4,0.5-0.005,6);
% seg_f=linspace(0.5+0.005,0.6,6);
% seg_g=linspace(0.6,0.7,10);
% seg_h=linspace(0.7,0.8,21);
% seg_i=linspace(0.8,0.9,10);
% seg_l=linspace(0.9+0.005,1,6);
% 
% for i=1:length(seg_a)
% seg_1(i)=seg_a(1,i);
% end
% for i=1:length(seg_b)-1
% seg_2(i)=seg_b(1,i);
% end
% for i=1:length(seg_c)-1
% seg_3(i)=seg_c(1,i);
% end
% for i=1:length(seg_d)-1
% seg_4(i)=seg_d(1,i);
% end
% for i=1:length(seg_e)
% seg_5(i)=seg_e(1,i);
% end
% for i=1:length(seg_f)
% seg_6(i)=seg_f(1,i);
% end
% for i=1:length(seg_g)-1
% seg_7(i)=seg_g(1,i);
% end
% for i=1:length(seg_h)-1
% seg_8(i)=seg_h(1,i);
% end
% for i=1:length(seg_i)
% seg_9(i)=seg_i(1,i);
% end
% for i=1:length(seg_l)
% seg_10(i)=seg_l(1,i);
% end
% M=zeros(T,n);
% M(1,:)=[seg_1 seg_2 seg_3 seg_4 seg_5 seg_6 seg_7 seg_8 seg_9 seg_10];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Network formation
% The following cods calculates the endogenous network and the distribution
% of opinions over time for a given choice of parameters. There is the
% option to run the code for multiple values of ``f'' and ``V''. 
%--------------------------------------------------------------------------

for w=1:length(f)
for z=1:length(V)
for t=2:T
    a=M(t-1,:); % (vector of opinions)
    for i=1:n
       
        F=net_sel(a,i,n,f(w),V(z));
        h=sum(F>0);
        S(i,t)=h;
        if h>0
            m=mean(a(logical(F)));
        end
        M(t,i)=(1-f(w))*a(i)+f(w)*m;
        A(i,:)=F;
    end
%% Uncomment the following lines to generate networks for a specified time 

%                 if t==2 || t==3 || t==4 ||t==5 || t==6 || t==7
%                     [Z(w,z),sizes,member] = networkComponents(A);
%                      W(w,z)=diameter(A);
%                      uW(w,z)=diameter(((A'+A)>0)+0);
%                      G=digraph(A);
%                      deg_out=outdegree(G);
%                      deg_in=indegree(G);
%                      figure
%                      plot(G,'Layout','force')
%                 end
     if t==T
         Y(w,z)=diameter(A);
     end
    B(:,t)=M(t,:)-M(t-1,:);
    if t==2
        G_1=digraph(A);
        figure
        %plot(G,'Layout','force')
        deg_out=outdegree(G_1);
        deg_in=indegree(G_1);
        axisx=linspace(0,n,n);
        figure
        plot(axisx,deg_in)
        diam=diameter(A);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Polarization
% Introduces a measure of polarization based on Esteban and Ray (1994). We
% refer to this measure in Appendix D in the paper.
%--------------------------------------------------------------------------

alpha(1,1)=0.8;
alpha(1,2)=1;
alpha(1,3)=1.6;
K=4;
s=50;
for t=1:T
    for o=1:3
        p(t,o)=polarizationrev(M,alpha(o),t,K,n,s);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plotting results
% The code here allows to generate figures that reproduce the evolution of
% opinions over time. NOTE: this part is still part of the loop started in
% ``Network formation''.
%--------------------------------------------------------------------------

 x=sign(B);
 for i=1:n
     for j=1:T-1
        y(i,j)=x(i,j)*x(i,j+1);
     end
 end
 grid(w,z)=sum(sum((y(:)<0)));
 if  any(y(:)<0)==1
    grid(w,z)=1;
 else
    grid(w,z)=0;
 end
    figure
    for ii=1:n
    plot(M(:,ii));
    hold on
    end
%--------------------------------------------------------------------------
%% DO NOT UNCOMMENT --> THIS CLOSES THE LOOP
end 

datetime('now')
end
save('results')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plotting the evolution of polarization
%Plot the measured polarization over time for the current selection of parameters with 3 different values of
%alpha: 0.8, 1.0 and 1.1, respectively.
%--------------------------------------------------------------------------

axis=linspace(1,T,T);
figure
plot(axis, p(:,1))
hold on
plot(axis, p(:,2))
plot(axis,p(:,3))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%