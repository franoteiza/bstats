/*=====================================================================================*
bststatsFO:     Create bootstrapped t-values for non-standard MHT correction of p-values
Purpose:        Romano & Wolf R code for calculating adjusted p-values requires the input
                of bootstrapped t-values from the models being estimated. This code 
                creates these files ready for upload to R to run the MHT correction code

Notes:          This code is 99% based on the rwolf command by Damian Clarke (@damiancclarke)
                but stops short of actually calculating MHT adjusted p-values. This is 
                because it is meant for use in non-standard examples, where the outcome
                variable is the same one but what changes in each hypotheses tested are
                the independent variables. 

Date created:   06/06/2018
Last edit:      02/09/2020
Created by:     Francisco Oteiza
*=====================================================================================*/


cap program drop bststatsFO
program bststatsFO, eclass
vers 14.2
#delimit ;
syntax varlist(min=1 fv ts numeric) [if] [in] [pweight fweight aweight iweight],
indepvar(varlist max=1)
file(name)
[
 method(name)
 controls(varlist fv ts)
 seed(numlist integer >0 max=1)
 reps(integer 100)
 Verbose
 strata(varlist)
 otherendog(varlist)
 cluster(varlist)
 iv(varlist)
 indepexog
 bl(name)
 *
 ]
;
#delimit cr
cap set seed `seed'
if `"`method'"'=="" local method regress
if `"`method'"'=="ivreg2"|`"`method'"'=="ivreg" {
    dis as error "To estimate IV regression models, specify method(ivregress)"
    exit 200
}
if `"`method'"'!="ivregress"&length(`"`indepexog'"')>0 {
    dis as error "indepexog argument can only be specified with method(ivregress)"
    exit 200
}
if `"`method'"'=="ivregress" {
    local ivr1 "("
    local ivr2 "=`iv')"
    local method ivregress 2sls
    if length(`"`iv'"')==0 {
        dis as error "Instrumental variable(s) must be included when specifying ivregress"
        dis as error "Specify the IVs using iv(varlist)"
        exit 200
    }    
}
else {
    local ivr1
    local ivr2
    local otherendog
}

local bopts
if length(`"`strata'"')!=0  local bopts `bopts' strata(`strata')
if length(`"`cluster'"')!=0 local bopts `bopts' cluster(`cluster')

if length(`"`verbose'"')==0 local q qui

*------------------------------------------------------------------------------------
*--- Run original regressions and collect their point estimates and (unadjusted) s.e.
*------------------------------------------------------------------------------------    
local j=0
local cand
local wt [`weight' `exp']

tempname nullvals
tempfile nullfile
file open `nullvals' using "`nullfile'", write all

foreach var of varlist `varlist' {
    local ++j
    file write `nullvals' "b`j'; se`j';"
}

local j=0

`q' dis "Displaying original (uncorrected) models:"
foreach var of varlist `varlist' {
    local ++j
    local Xv `controls'
    if length(`"`bl'"')!=0 local Xv `controls' `var'`bl' 
    if length(`"`indepexog'"')==0 {
        `q' `method' `var' `ivr1'`indepvar' `otherendog'`ivr2' `Xv' `if' `in' `wt', `options'
    }
    else {
        `q' `method' `var' `ivr1'`otherendog'`ivr2' `indepvar' `Xv' `if' `in' `wt', `options'
    }
    if _rc!=0 {
        dis as error "Your original `method' does not work."
        dis as error "Please test the `method' and try again."
        exit _rc
    }
    local t`j' = abs(_b[`indepvar']/_se[`indepvar'])
    local n`j' = e(N)-e(rank)
    if `"`method'"'=="areg" local n`j' = e(df_r)
    local cand `cand' `j'
    
    if `j'==1 file write `nullvals' _n "`= _b[`indepvar']';`= _se[`indepvar']'"
    else file write `nullvals' ";`= _b[`indepvar']';`= _se[`indepvar']'"
}

*-------------------------------------------------------------------------------
*--- Run bootstrap reps to create null Studentized distribution
*-------------------------------------------------------------------------------  
dis "Running `reps' bootstrap replications for each variable.  This may take some time"
forvalues i=1/`reps' {
    if length(`"`verbose'"')!=0 dis "Bootstrap sample `i'."
    local j=0
    preserve
    bsample `if' `in', `bopts'
    
    foreach var of varlist `varlist' {
        local ++j
        local Xv `controls'
        if length(`"`bl'"')!=0 local Xv `controls' `var'`bl' 
        if length(`"`indepexog'"')==0 {
            qui `method' `var' `ivr1'`indepvar'  `otherendog'`ivr2' `Xv' `if' `in' `wt', `options'
        }
        else {
            qui `method' `var' `ivr1'`otherendog'`ivr2' `indepvar' `Xv' `if' `in' `wt', `options'
        }
        if `j'==1 file write `nullvals' _n "`= _b[`indepvar']';`= _se[`indepvar']'"
        else file write `nullvals' ";`= _b[`indepvar']';`= _se[`indepvar']'"
    }
    restore
}

preserve
file close `nullvals'
qui insheet using `nullfile', delim(";") names clear // removed qui

*-------------------------------------------------------------------------------
*--- Create null t-distribution
*-------------------------------------------------------------------------------
forval x=1/`j' {
    qui sum b`x'
    gen t_`x'=abs((b`x'-r(mean))/r(sd))
    replace t_`x' = abs((b`x')/se`x') if _n==1
}

keep t_*
export delimited using  "`file'_tnull.csv",  replace
restore

end
