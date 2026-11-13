% This code computes demographic variables for 12 Eurozone coutries,
% launches hetprob.m and herprobnet.m with different parameterizations,
% and saves the output. Based on multi.m for U.S. states.

% Copyright by Marco Bassetto, Vadym Lepetyuk
% This code can be freely distributed and modified for research purposes only, 
% provided this copyright notice is included in the modified code. 
% Proper credit should be given in all publications arising from
% modifications of this code; this should include a citation of 
% "Government Investment and the European Stability and Growth Pact"
% by Marco Bassetto and Vadym Lepetyuk

optdebtall=zeros(12,4); % Matrix of optimal values of debt
bball=zeros(12,4); % Wedge under balanced budget
fdall=zeros(12,4); % Wedge under golden rule

alph=0;
bet=0.96;

% list of Eurosone countries (ordered as appears in Eurostat statistics)
%lbl={'be' 'de' 'ie' 'gr' 'es' 'fr' 'it' 'lu' 'nl' 'at' 'pt' 'fi'};
lbl={'Belgium' 'Germany' 'Ireland' 'Greece' 'Spain' 'France'...
    'Italy' 'Luxembourg' 'Netherlands' 'Austria' 'Portugal' 'Finland'};

load eu12population.txt;
load eu12mortality.txt;
load eu8migration.txt;
load eu12popgrowth.txt;

if (size(eu12population,1)~=73 | size(eu12population,2)~=12)
    error('Incorrect size: eu12population');
end;
if (size(eu12mortality,1)~=72 | size(eu12mortality,2)~=12)
    error('Incorrect size: eu12mortality');
end;
if (size(eu8migration,1)~=15 | size(eu8migration,2)~=12)
    error('Incorrect size: eu12migration');
end;
if (size(eu12popgrowth,1)~=1 | size(eu12popgrowth,2)~=12)
    error('Incorrect size: eu12popgrowth'); 
end;

% migration: linear interpolation of 5-year cells;
% note: migration data are only awailable for 8 countries
% Belgium, Germany, Spain, Italy, Luxembourg, Netherlands, Austria, Finland
migration=interp1(17:5:92,[eu8migration; eu8migration(15,:)],(18:90)');

% death probability
death=eu12mortality./eu12population(1:72,:);
death=[death; ones(1,12)];

% total survival probability
thetall=(1-death).*(1-migration);

% age profile; sum normalized to 1
ageall=eu12population./kron(sum(eu12population),ones(73,1));

% run the loop only for 8 countries w/ migration data
countryset = [1 2 5 7 8 9 10 12];
for ii=countryset,
    fprintf('\nCountry %i of 12: %s\n',ii,lbl{ii})
    
    n=eu12popgrowth(ii);
    thet=thetall(:,ii);
    age=ageall(:,ii);
    delt=.06;
    hetprob;
    optdebtall(ii,1)=optdebt;
    bball(ii,1)=bb;
    fdall(ii,1)=fd;
 
    delt=.03;
    hetprob;
    optdebtall(ii,2)=optdebt;
    bball(ii,2)=bb;
    fdall(ii,2)=fd;
    
    delt=.06;
    hetprobnet;
    optdebtall(ii,3)=optdebt;
    bball(ii,3)=bb;
    fdall(ii,3)=fd;
    
        if ii==2,
        pause;
    end;


    delt=.03;
    hetprobnet;
    optdebtall(ii,4)=optdebt;
    bball(ii,4)=bb;
    fdall(ii,4)=fd;

    if exist('doplot','var'),
        plot(18:90,perdev0,18:90,bb*ones(73,1));
        axis([18 90 1 3.5]);
        title(['Voting outcome under balanced budget in ' lbl{ii}]);
        if ii>1,
            print -dpsc bbvot.ps -append
        else,
            print -dpsc bbvot.ps
        end;
    end;
end;
    
%dosave=1;
if exist('dosave','var'),
    save eu8;
end;

%dotex=1;
if exist('dotex','var'),
    fid = fopen('eu8table.tex','w');
    fprintf(fid,'\\documentclass{article}\n');
    fprintf(fid,'\\begin{document}\n');
    fprintf(fid,'\\begin{table}\n');
    fprintf(fid,'\\begin{center}\n');
    fprintf(fid,'\\begin{tabular}{cccc}\n');
    fprintf(fid,'Country & SGP, no exclusions & Excluding gross investment & Excluding net investment\\\\\n');
    fprintf(fid,'\\multicolumn{4}{c}{Generic capital}\\\\\n');
    [lblabc,iiabc]=sort(lbl(countryset)); %list of countries in alphabetical order
    countrysetabc=countryset(iiabc);
    for ii=countrysetabc,
        fprintf(fid,'%-12s & %6.0f\\%% & %6.0f\\%% & %6.1f\\%%\\\\\n',...
            lbl{ii},(bball(ii,1)-1)*100,(fdall(ii,1)-1)*100,(fdall(ii,3)-1)*100);
    end;
    fprintf(fid,'\\multicolumn{4}{c}{Major infrastructure}\\\\\n');
    for ii=countrysetabc,
        fprintf(fid,'%-12s & %6.0f\\%% & %6.0f\\%% & %6.1f\\%%\\\\\n',...
            lbl{ii},(bball(ii,2)-1)*100,(fdall(ii,2)-1)*100,(fdall(ii,4)-1)*100);
    end;
    fprintf(fid,'\\end{tabular}\n');
    fprintf(fid,'\\end{center}\n');
    fprintf(fid,'\\caption{Efficiency wedge $\\tau$ in the baseline calibration}\n');
    fprintf(fid,'\\end{table}\n');
    fprintf(fid,'\\end{document}\n');
    fclose(fid);
end;