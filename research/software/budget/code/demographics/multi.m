% This code computes demographic variables for all 50 US states,
% launches hetprob.m with the different parameterizations,
% and saves the output

% Copyright by Marco Bassetto, Vadym Lepetyuk
% This code can be freely distributed and modified for research purposes only, 
% provided this copyright notice is included in the modified code. 
% Proper credit should be given in all publications arising from
% modifications of this code; this should include a citation of 
% "Politics and Efficiency of Separating Capital and Ordinary Government Budgets"
% by Marco Bassetto with Thomas J. Sargent

% We stack all results; there are 53 sets of results (all states, US now
% and back in time, with complete population or only males voting), 
% times 2 values for the depreciation.

optdebtall=zeros(53,2); % Matrix of optimal values of debt
bball=zeros(53,2); % Wedge under balanced budget
fdall=zeros(53,2); % Wedge under golden rule
alphall=zeros(53,2); % Maturity of debt that makes the golden rule optimal

alph=1-.5^(1/15);
bet=0.96;

n=.0124;
load usdeathsnow.txt;
thet=[1-usdeathsnow; 0];
load usagenow.txt;
age=usagenow/sum(usagenow);
delt=.06;

hetprob;
optdebtall(51,1)=optdebt;
bball(51,1)=bb;
fdall(51,1)=fd;
alphall(51,1)=alphcurr;
save usnow;

delt=.03;
hetprob;
save usnowld;
optdebtall(51,2)=optdebt;
bball(51,2)=bb;
fdall(51,2)=fd;
alphall(51,2)=alphcurr;

n=.023;
load usdeaths1880.txt;
thet=1-usdeaths1880; % Note inconsistency
load usage1880.txt;
age=usage1880/sum(usage1880);
delt=.06;

hetprob;
optdebtall(52,1)=optdebt;
bball(52,1)=bb;
fdall(52,1)=fd;
alphall(52,1)=alphcurr;
save us1880;

delt=.03;
hetprob;
save us1880ld;
optdebtall(52,2)=optdebt;
bball(52,2)=bb;
fdall(52,2)=fd;
alphall(52,2)=alphcurr;

lbl=['AL' 'AK' 'AZ' 'AR' 'CA' 'CO' 'CT' 'DE' 'FL' 'GA' 'HI' 'ID' 'IL'...
    'IN' 'IA' 'KS' 'KY' 'LA' 'ME' 'MD' 'MA' 'MI' 'MN' 'MS' 'MO' 'MT'...
    'NE' 'NV' 'NH' 'NJ' 'NM' 'NY' 'NC' 'ND' 'OH' 'OK' 'OR' 'PA' 'RI'...
    'SC' 'SD' 'TN' 'TX' 'UT' 'VT' 'VA' 'WA' 'WV' 'WI' 'WY'];

load popgrowthallstates.txt;
load migrationallstates.txt;
load pop1990allstates.txt;
load pop2000allstates.txt;
mig=reshape(migrationallstates',8,18,50);
% Compute outmigration rate
outrate=squeeze(mig(6,4:18,:)./(mig(1,4:18,:)-mig(7,4:18,:)-mig(8,4:18,:)));
% First, compute the migration rate for the relevant groups.
% We will use interpolation to compute yearly numbers
% We have data aggregated for 5 years, and movers moved some time
% during the last 5 years. Hence, a person that is now 15-19 moved
% some time between when s/he was 10 and 18. Centering the interval,
% we attribute the migration rate of the 15-19 group to age 14, and so
% on. For the last group (85 and above), we center at 84, and keep the
% rate constant from then on. This is a bit arbitrary, but it is also
% irrelevant, since for this group mortality is a lot more important
% than mobility
outyearly=ones(73,50);
for ii=1:50,
    outyearly(1:67,ii)=1-(1-interp1(14:5:84,outrate(:,ii),18:84)).^(1/5);
    outyearly(68:72,ii)=1-(1-outrate(15,ii))^(1/5);
end;

% For the age structure, we use the average of the age structure in 2000 
% and 1990. Since we have data for 1990 that are somewhat aggregated by
% age, we use the age structure of 2000 within each cell

select=[1:8 10:51]; % This gets rid of DC and Puerto Rico

% This sums males and females aged 18 and above. Though we truncate at 90,
% we need to carry over the old for imputing numbers in 1990, when our data
% codes people 85 and above
agestate2000=(pop2000allstates(select,21:105)+pop2000allstates(select,125:209))';

agefrac2000=zeros(73,50);
for ii=1:50,
    agefrac2000(:,ii)=agestate2000(1:73,ii)/sum(agestate2000(1:73,ii));
end;

% For 1990, DC was cleaned out before saving
agestate1990=zeros(73,50);
agestate1990(1:4,:)=pop1990allstates(:,13:16)';

for ii=1:50,
    agestate1990(5:7,ii)=pop1990allstates(ii,17)*agefrac2000(5:7,ii)/...
        sum(agefrac2000(5:7,ii));
    for jj=1:7,
        agestate1990(3+5*jj:7+5*jj,ii)=pop1990allstates(ii,17+jj)*...
            agefrac2000(3+5*jj:7+5*jj,ii)/sum(agefrac2000(3+5*jj:7+5*jj,ii));
    end;
    agestate1990(43:44,ii)=pop1990allstates(ii,25)*agefrac2000(43:44,ii)/...
        sum(agefrac2000(43:44,ii));
    agestate1990(45:47,ii)=pop1990allstates(ii,26)*agefrac2000(45:47,ii)/...
        sum(agefrac2000(45:47,ii));
    for jj=1:4,
        agestate1990(43+5*jj:47+5*jj,ii)=pop1990allstates(ii,26+jj)*...
            agefrac2000(43+5*jj:47+5*jj,ii)/sum(agefrac2000(43+5*jj:47+5*jj,ii));
    end;
    agestate1990(68:73,ii)=pop1990allstates(ii,31)*agestate2000(68:73,ii)/...
        sum(agestate2000(68:85,ii));
    agefrac1990(:,ii)=agestate1990(:,ii)/sum(agestate1990(:,ii));
end;

agefrac=(agefrac2000+agefrac1990)/2;

for ii=1:50,
    n=popgrowthallstates(ii);
    % The migration rate is conditional on survival (because of the way
    % we computed the starting number of people). We now need to integrate
    % mortality
    thet=[1-usdeathsnow; 0].*(1-outyearly(:,ii));
    age=agefrac(:,ii);
    delt=.06;
    hetprob;
    optdebtall(ii,1)=optdebt;
    bball(ii,1)=bb;
    fdall(ii,1)=fd;
    alphall(ii,1)=alphcurr;

    delt=.03;
    hetprob;
    optdebtall(ii,2)=optdebt;
    bball(ii,2)=bb;
    fdall(ii,2)=fd;
    alphall(ii,2)=alphcurr;
    plot(18:90,perdev0,18:90,bb*ones(73,1));
    axis([18 90 1.25 3.5]);
    title(['Voting outcome under balanced budget in ' lbl(2*ii-1:2*ii)]);
    if ii>1,
        print -dps bbvot.ps -append
    else,
        print -dps bbvot.ps
    end;
end;
    
n=.023;
load usdm1880.txt;
thet=[1-usdm1880'; 0];
load usam1880.txt;
age=usam1880/sum(usam1880);
delt=.06;

hetprob;
optdebtall(53,1)=optdebt;
bball(53,1)=bb;
fdall(53,1)=fd;
alphall(53,1)=alphcurr;

delt=.03;
hetprob;
optdebtall(53,2)=optdebt;
bball(53,2)=bb;
fdall(53,2)=fd;
alphall(53,2)=alphcurr;

halflifeall=-log(2)./log(1-alphall); % Halflife of debt implied by alphall

save multi;