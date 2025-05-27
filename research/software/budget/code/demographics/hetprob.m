% This code calculates the optimal debt repayment as a share of gross investment

% Copyright by Marco Bassetto, Vadym Lepetyuk
% This code can be freely distributed and modified for research purposes only, 
% provided this copyright notice is included in the modified code. 
% Proper credit should be given in all publications arising from
% modifications of this code; this should include a citation of 
% "Politics and Efficiency of Separating Capital and Ordinary Government Budgets"
% by Marco Bassetto with Thomas J. Sargent

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

% Plot the optimal debt financing for each age
if exist('doplot','var'),
    figure(1);
    plot(1:tcap,optx,1:tcap,optdebt*ones(tcap,1));
    title('Optimal Debt Profile');
    legend('optx','optdebt',2);
    xlabel('Age');
    ylabel('x');
end;

% Now, check that the incentives work indeed for the
% median voter.

dev=QsC-QsB*optdebt;

[devsort,idevsort]=sort(dev);

% check: medvot2 is the same as medvot
votpower2=cumsum(age(idevsort));
medvot2=min(find(votpower2>0.5));

% check: optdev is zero
optdev=devsort(medvot2);

if dev(1)>0,
   small=sort(isort(1:medvot-1));
   big=sort(isort(medvot+1:tcap));
else,
   big=sort(isort(1:medvot-1));
   small=sort(isort(medvot+1:tcap));
end;

% Now, compute efficiency costs of having x=1 or x=0

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

% Plot the efficiency costs against x

if exist('doplot','var'),
    xgrid=(0.7:0.005:1.3)*optdebt;
    devgrid=QsC*ones(size(xgrid))-QsB*xgrid;

    [devsortgrid,isortgrid]=sort(devgrid,1);

    votpowergrid=cumsum(age(isortgrid),1);
    
    medvotgrid=zeros(size(xgrid));
    meddevgrid=zeros(size(xgrid));
    for i=1:length(xgrid)
        medvotgrid(i)=min(find(votpowergrid(:,i)>0.5));
        meddevgrid(i)=devsortgrid(medvotgrid(i),i);
    end;

    perdevgrid=(pareto+meddevgrid)/pareto;

    figure(2);
    plot(xgrid,perdevgrid);
    title('Efficiency wedge');
    xlabel('x');
    ylabel('(Pareto choice = 1)');
    text(1.05*optdebt,1,sprintf('optdebt=%.4f',optdebt));
end;

% Compute the maturity structure that makes the golden rule efficient
% Root Finding Implementation: False Position Method

if fd<1,
    alphlow=alph;
    fdlow=fd;
    alphhigh=2*delt;  %initial guess for alphhigh
    alphcurr=alphhigh;
    for i=1:tcap,
        tempsum=0;
        for m=2:tcap-i,
            tempsum=tempsum+prod(thet((i):(i+m-1)))*...
                (1-alphcurr)^(m-2)*bet^m/(1+n)^m;
        end;
        dev1(i)=bet*(1-delt)*(1-thet(i)/(1+n))...
            -1+(bet*thet(i)/(1+n))*(alphcurr+1/bet-delt)+...
            (-1+alphcurr+1/bet)*(delt-alphcurr)*tempsum;
    end;
    perdev1=(pareto+dev1)/pareto;
    [devsort1,isort1]=sort(perdev1);
    votpower1=cumsum(age(isort1));
    medvot1=min(find(votpower1>0.5));
    fdhigh=devsort1(medvot1);
    if fdhigh<1,
        error('alpha is not high enough');
    end;
else,
    alphlow=delt;
    alphhigh=alph;
    fdhigh=fd;
    alphcurr=delt;
    for i=1:tcap,
        tempsum=0;
        for m=2:tcap-i,
            tempsum=tempsum+prod(thet((i):(i+m-1)))*...
                (1-alphcurr)^(m-2)*bet^m/(1+n)^m;
        end;
        dev1(i)=bet*(1-delt)*(1-thet(i)/(1+n))...
            -1+(bet*thet(i)/(1+n))*(alphcurr+1/bet-delt)+...
            (-1+alphcurr+1/bet)*(delt-alphcurr)*tempsum;
    end;
    perdev1=(pareto+dev1)/pareto;
    [devsort1,isort1]=sort(perdev1);
    votpower1=cumsum(age(isort1));
    medvot1=min(find(votpower1>0.5));
    fdlow=devsort1(medvot1);
    if fdlow>1,
        error('alpha is not low enough');
    end;
end;

epsi=1;
while epsi>1e-08,
    alphcurr=(alphhigh*(1-fdlow)+alphlow*(fdhigh-1))/(fdhigh-fdlow);
    %fprintf('Processing alph=%.4f...\n',alphcurr);
    for i=1:tcap,
        tempsum=0;
        for m=2:tcap-i,
            tempsum=tempsum+prod(thet((i):(i+m-1)))*...
                (1-alphcurr)^(m-2)*bet^m/(1+n)^m;
        end;
        dev1(i)=bet*(1-delt)*(1-thet(i)/(1+n))...
            -1+(bet*thet(i)/(1+n))*(alphcurr+1/bet-delt)+...
            (-1+alphcurr+1/bet)*(delt-alphcurr)*tempsum;
    end;
    perdev1=(pareto+dev1)/pareto;
    [devsort1,isort1]=sort(perdev1);
    votpower1=cumsum(age(isort1));
    medvot1=min(find(votpower1>0.5));
    fdev=devsort1(medvot1);
    epsi=abs(fdev-1);
    if fdev>1,
        alphhigh=alphcurr;
        fdhigh=fdev;
    else,
        alphlow=alphcurr;
        fdlow=fdev;
    end;
end;
fprintf('Optimal alph for "golden rule": %.4f\n\n',alphcurr);
    
% Plot the efficiency costs against alph

if exist('doplot','var'),
    alphgrid=[0.5*delt:.005:2*delt];
    nalphgrid=length(alphgrid);
    fd=zeros(1,nalphgrid);
    for j=1:nalphgrid,
        %fprintf('Processing alph=%.4f...\n',alphgrid(j));
        QsB=zeros(tcap,1);
        for i=1:tcap,
           tempsum=0;
           for m=2:tcap-i,
               tempsum=tempsum+prod(thet((i):(i+m-1)))*...
                   (1-alphgrid(j))^(m-2)*bet^m/(1+n)^m;
           end;
           dev1(i)=QsC(i)-1+(bet*thet(i)/(1+n))*(alphgrid(j)+1/bet-delt)+...
               (-1+alphgrid(j)+1/bet)*(delt-alphgrid(j))*tempsum;
        end;
        perdev1=(pareto+dev1)/pareto;
        [devsort1,isort1]=sort(perdev1);
        votpower1=cumsum(age(isort1));
        medvot1=min(find(votpower1>0.5));
        fdgrid(j)=devsort1(medvot1);
    end;

    figure(3);
    plot(alphgrid,fdgrid);
    title('Efficiency wedge'); 
    xlabel('alph');
    ylabel('(Pareto choice = 1)');
end;