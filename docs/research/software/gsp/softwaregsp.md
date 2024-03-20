## Government Investment and the European Stability and Growth Pact (joint with Vadym Lepetyuk)
*Economic Perspectives*, Federal Reserve Bank of Chicago, 2007, vol. 31, n.3, pp.33-43.

The code for this project is written in Matlab, except for a Stata code that pulls data from the Eurobarometer survey.

Main code for the baseline exercise: it loads the data for all countries, and calls the subroutines for each country and each choice of the depreciation rate: [eu8.m](code/eu8.m)

Main code for the Eurobarometer exercise: [eu12.m](code/eu12.m)

Subroutine that computes wedges for the case in which gross investment is excluded from the deficit count: [hetprob.m](code/hetprob.m)

Same as above for net investment: [hetprobnet.m](code/hetprobnet.m)

### Files that contain demographic data needed as inputs
Mortality by age: [eu12mortality.txt](code/eu12mortality.txt)

Population structure by age: [eu12population.txt](code/eu12population.txt)

Population growth: [eu12popgrowth.txt](code/eu12popgrowth.txt)

Baseline emigration rates: [eu8migration.txt](code/eu8migration.txt)

Excel spreadsheet that contains emigration data sent to us by Anna Lööf  at Eurostat (note: all other data were downloaded from the Eurostat web site or their publications): [emigration.xls](code/emigration.xls)

Stata code that handles the Eurobarometer survey data: [eu12migration.do](code/eu12migration.do)

Note: the Eurobarometer data are not freely available, so we can not include them here. They can be downloaded from the archives mentioned here: [http://www.gesis.org/en/data_service/eurobarometer/staff/archives.htm](http://www.gesis.org/en/eurobarometer/home/)
