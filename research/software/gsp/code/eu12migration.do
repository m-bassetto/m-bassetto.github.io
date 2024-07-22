* This Stata file calculates out-of-country migration rates for Euro 12 countries
* based on Eurobarometer 64.1 survey (2005)

* Copyright by Marco Bassetto, Vadym Lepetyuk
* This code can be freely distributed and modified for research purposes only, 
* provided this copyright notice is included in the modified code. 
* Proper credit should be given in all publications arising from
* modifications of this code; this should include a citation of 
* "Government Investment and the European Stability and Growth Pact"
* by Marco Bassetto and Vadym Lepetyuk

* The relevant question of the survey is
* QA18. Do you think that in the next five years you are likely to move…?
* (1) In the same city/town/village
* (2) To another city/town/village but in the same region
* (3) To another region but in the same country
* (4) To another country in the European Union
* (5) To another country outside the European Union
* (6) You don't think you will move
* (7) Don't know

* Original dataset is restricted to Euro 12 countries.
* Original dataset is restricted to respondents residing in their home countries.
* Population weights are used in the calcucations.
* Respondents are aggregated into 10-year age groups.
* Annual probability of moving is calculated (by diving by 5).

* The survey data are available from ICPSR (http://www.icpsr.umich.edu)
* Univ of Minnesota is among the subscribing institutions.

* Dataset description:
* Papacostas, Antonis. EUROBAROMETER 64.1: MOBILITY, FOOD RISK, SMOKING, AIDS
* PREVENTION, AND MEDICAL ERRORS, SEPTEMBER-OCTOBER 2005 [Computer file].
* ICPSR04641-v1. Brussels, Belgium: TNS Opinion & Social [producer], 2006.
* Cologne, Germany: Zentralarchiv fur Empirische Sozialforschung/Ann Arbor,
* MI: Inter-university Consortium for Political and Social Research
* [distributors], 2007-02-28.

**version 9.0
clear
set mem 50m
set more off
use 04641-0001-Data.dta, clear

* Restrict the dataset to Euro 12 countries
keep if VEUZONE==1   /* Euro countries; note: this is not the same as VEU12 */

* Merge East and West Germany
replace COUNTRY=3 if COUNTRY==4
label define COUNTRY 3 "GERMANY", modify
* Merge Northern Ireland and Great Britain
**replace COUNTRY=16 if COUNTRY==17
**label define COUNTRY 16 "UNITED KINGDOM", modify

* Keep citizens residing in their home countries
drop if COUNTRY==1 & Q1_1==0     /*  BELGIUM  */
drop if COUNTRY==3 & Q1_3==0     /*  GERMANY  */
drop if COUNTRY==5 & Q1_4==0     /*  GREECE  */
drop if COUNTRY==6 & Q1_5==0     /*  SPAIN  */
drop if COUNTRY==7 & Q1_15==0     /*  FINLAND  */
drop if COUNTRY==8 & Q1_6==0     /*  FRANCE  */
drop if COUNTRY==9 & Q1_7==0     /*  IRELAND  */
drop if COUNTRY==10 & Q1_8==0     /*  ITALY  */
drop if COUNTRY==11 & Q1_9==0     /*  LUXEMBOURG  */
drop if COUNTRY==12 & Q1_10==0     /*  NETHERLANDS  */
drop if COUNTRY==13 & Q1_13==0     /*  AUSTRIA  */
drop if COUNTRY==14 & Q1_11==0     /*  PORTUGAL  */

* Define respondents that expect to migrate out of country
generate int migrate=(QA18_4|QA18_5)&(~QA18_7)
* Define respondents that expect to stay within the country
generate int stay=(QA18_1|QA18_2|QA18_3|QA18_6)&(~QA18_7)
* Define respondents that know whether they plan to migrate or stay
* (in other words, any answer except "don't know")
generate int respond=(~QA18_7)

* Define age, country
rename VD11 age
rename COUNTRY country

* Tabulate age by 5-year groups
**recode age 15/19=1 20/24=2 25/29=3 30/34=4 35/39=5 40/44=6 45/49=7 ///
**           50/54=8 55/59=9 60/64=10 65/69=11 70/74=12 75/79=13 80/max=14

* Tabulate age by 10-year groups
recode age 15/24=1 25/34=2 35/44=3 45/54=4 55/64=5 65/74=6 75/max=7
**recode age min/17=. 18/27=1 28/37=2 38/47=3 48/57=4 58/67=5 68/77=6 78/max=7
**recode age min/17=. 18/24=1 25/34=2 35/44=3 45/54=4 55/64=5 65/74=6 75/max=7
**recode age min/17=. 18/27=1 28/37=2 38/47=3 48/57=4 58/67=5 68/max=6
**drop if age==.

* Display number of respondents by age and country
tabulate country age

* Display number of respondents by age and country that expressed the desire to migrate out of country
tabulate country age if migrate

* Create national population weights
generate nw=W1
replace nw=W3 if country==3
replace nw=W4 if country==16

* Create Euro 12 population weights
**generate ew=W12

* EU12 unionwide summation by age
** collapse (mean) migrate stay respond [aweight=ew], by (age) fast

* EU12 country-by-country summation by age
collapse (mean) migrate stay respond [aweight=nw], by (age country) fast

* Out-of-country annual migration rate
generate rate=migrate/respond/5

* Note: (migrate+stay)/respond>=1 because some respondents expressed both the desire to migrate and to stay
* such respondents are counted as half toward the out-of-country annual migration rate
**generate rate=(migrate-0.5*(migrate+stay-respond))/respond/5

* Display out-of-country annual migration rate
table country age, c(mean rate)

* Graphs
**plot rate age

* Reshape data
drop migrate stay respond
reshape wide rate, i(age) j(country)
rename rate1 mr_BE  /*  BELGIUM  */
rename rate3 mr_DE  /*  GERMANY  */
rename rate5 mr_EL  /*  GREECE  */
rename rate6 mr_ES  /*  SPAIN  */
rename rate7 mr_FI  /*  FINLAND  */
rename rate8 mr_FR  /*  FRANCE  */
rename rate9 mr_IE  /*  IRELAND  */
rename rate10 mr_IT  /*  ITALY  */
rename rate11 mr_LU  /*  LUXEMBOURG  */
rename rate12 mr_NL  /*  NETHERLANDS  */
rename rate13 mr_AT  /*  AUSTRIA  */
rename rate14 mr_PT  /*  PORTUGAL  */

* Write .csv file
outsheet using migrate.csv, nolabel noquote comma replace
