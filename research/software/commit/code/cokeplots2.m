%This code computes the Ramsey outcome and the function lhat
%to verify implementability of the sovereign debt example

% Copyright by Marco Bassetto
% This code can be freely distributed and modified for research purposes only, 
% provided this copyright notice is included in the modified code. 
% Proper credit should be given in all publications arising from
% modifications of this code; this should include a citation of 
% "Equilibrium and Government Commitment" by Marco Bassetto

warning off optim:fmincon:SwitchingToMediumScale

%Parameter values
gam=0.5;
alph1=1;
alph2=2;
alph3=1;
B=.5;

%Distribution of the shock: taken to be uniform over
%(-qdl,+qdh)
qdl=0;
qdh=0.01;

%Construct a grid for L of size gridl
gridl=400;
l=linspace(0.01,5,gridl)';

%Compute the function lhat for the best strategy. 
%The Ramsey outcome is either the largest fixed point,
% or the default point
lhat=zeros(gridl,1);

options=optimset('TolX',1e-6,'TolFun',1e-8,'Display','off');

guess=1;
for i=1:gridl,
   L=l(i);
   lhat(i)=fsolve(@(x) cokefn2(x,alph1,alph2,alph3,gam,qdl,qdh,L,B),guess,options);   
   guess=lhat(i);
end;

ibest=max(find(lhat>l));
ibest=ibest+1;
lbest=l(ibest);

l1curve=(sin(1.2*(l-1.9)+.12)+1.3)/2;

% Plot the results
plot(l,l,'-',l,lhat,'-',l,(sin(1.2*(l-1.9)+.12)+1.3)/2,'-',lbest,lbest,'s');
title('Household best reply');
text(lbest-.1,lbest+.1,'R');
text(l(gridl/2)-.1,l1curve(gridl/2)+.1,'l_1');
text(l(gridl/2)-.1,lhat(gridl/2)+.1,'lhat');
xlabel('Labor supply of other households');
ylabel('Best reply');

% Save the results for PSTricks

%This creates the area to be shaded. Notice that it is augmented in a way
%that creates a polygon to be filled.
sdarea=[[0; l; l(gridl)] [0; lhat; 0]];
save sdarea1.dat sdarea -ascii;

l1curve=[l l1curve];

save l1curve.dat l1curve -ascii;

lbest=[lbest lbest; lbest lbest];

save lbest1.dat lbest -ascii;

% Now, consider the case that cannot be implemented

B=1;

%Compute the function lhat for the best strategy. 
%The Ramsey outcome is either the largest fixed point,
% or the default point
lhat=zeros(gridl,1);

guess=1;
for i=1:gridl,
   L=l(i);
   lhat(i)=fsolve(@(x) cokefn2(x,alph1,alph2,alph3,gam,qdl,qdh,L,B),guess,options);
   guess=lhat(i);
end;

ibest=max(find(lhat>l));
ibest=ibest+1;
lbest=l(ibest);

% Plot the results
plot(l,l,'-',l,lhat,'-',lbest,lbest,'s');
title('Household best reply');
text(lbest-.1,lbest+.1,'R');
text(l(gridl/2)-.1,lhat(gridl/2)+.1,'lhat');
xlabel('Labor supply of other households');
ylabel('Best reply');

sdarea=[[0; l; l(gridl)] [0; lhat; 0]];
save sdarea2.dat sdarea -ascii;

lbest=[lbest lbest; lbest lbest];

save lbest2.dat lbest -ascii;