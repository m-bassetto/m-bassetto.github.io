% This code calculates the optimal debt repayment and the efficiency wedges
% for linear preferences.

% Copyright by Marco Bassetto, Vadym Lepetyuk
% This code can be freely distributed and modified for research purposes only, 
% provided this copyright notice is included in the modified code. 
% Proper credit should be given in all publications arising from
% modifications of this code; this should include a citation of 
% "Politics and Efficiency of Separating Capital and Ordinary Government Budgets"
% by Marco Bassetto with Thomas J. Sargent

% Parameters
alph=.0452;
bet=.96;

n=0.02;
thet=[0.4; 0];
age=[1.02/1.42; 0.4/1.42];
delt=.06;

tcap=size(thet,1);

% Compute the expected present value of taxes
% per unit of public capital (formula (22) in QJE)
% >>> Qs = 1-bet*(1-delt) + QsC - QsB*x <<<
QsC=bet*(1-delt)*(1-thet/(1+n));
QsB=zeros(tcap,1);
for i=1:tcap,
   tempsum=0;
   for m=2:tcap-i,
       tempsum=tempsum+prod(thet((i):(i+m-1)))*...
           (1-alph)^(m-2)*bet^m/(1+n)^m;
   end;
   QsB(i)=1-bet*(alph+1/bet-delt)*thet(i)/(1+n)-...
       (-1+alph+1/bet)*(delt-alph)*tempsum;
end;

% Compute the amount of debt that
% matches incentives for each generation
optx=QsC./QsB;

[xsort,isort]=sort(optx);

votpower=cumsum(age(isort));
medvot=min(find(votpower>0.5));

optdebt=xsort(medvot);
fprintf('Optimal debt financing: %.4f\n',optdebt);

% Compute efficiency costs of having x=1 or x=0
dev0=QsC;
dev1=QsC-QsB;

pareto=1-bet*(1-delt);
perdev0=(pareto+dev0)/pareto;
perdev1=(pareto+dev1)/pareto;

[devsort0,isort0]=sort(perdev0);

votpower0=cumsum(age(isort0));
medvot0=min(find(votpower0>0.5));

bb=devsort0(medvot0);
fprintf('Efficiency wedge, balanced budget: %.4f\n',bb);

[devsort1,isort1]=sort(perdev1);

votpower1=cumsum(age(isort1));
medvot1=min(find(votpower1>0.5));

fd=devsort1(medvot1);
fprintf('Efficiency wedge, "golden rule":   %.4f\n',fd);

save linear;