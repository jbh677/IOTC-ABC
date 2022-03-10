//////////////////////////////////////////////////////
// projc.cpp /////////////////////////////////////////
//////////////////////////////////////////////////////
// R. Hillary CSIRO 2021 /////////////////////////////
//////////////////////////////////////////////////////

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
#include <iostream>
#include <R.h>
#include "assert.h"

using namespace Rcpp;

/*
 * Project forward in time for fixed exploitation rate
 * inputs:
 * 1.  Ninit: initial population
 * 2.  SSBinit: duh
 * 3.  epsinit: autocorrelation bollox
 * 4.  B0
 * 5.  steepness
 * 6.  sigmar
 * 7.  rho
 * 8.  M
 * 9.  weight-at-age
 * 10.  maturity-at-age
 * 11.  selectivity-at-age
 * 11. fleet catch share
 * 12. fixed catch biomass
 * 13. dimensions
 */
//[[Rcpp::export]]
RcppExport SEXP projc(SEXP Ninit_,SEXP SSBinit_,SEXP epsinit_,SEXP  B0_,SEXP hh_,SEXP SPR_,SEXP sigmar_,SEXP rho_,SEXP M_,SEXP wt_,SEXP mat_,SEXP sel_,SEXP cshare_,SEXP cfix_,SEXP dm_)
{
  IntegerVector dm = as<Rcpp::IntegerVector>(dm_);
  int ny = dm[0];
  int na = dm[1];
  int nf = dm[2];
  int a,f,i,y;

  // read ins

  NumericVector Ninit = as<Rcpp::NumericVector>(Ninit_);
  NumericVector M = as<Rcpp::NumericVector>(M_);
  NumericVector wt = as<Rcpp::NumericVector>(wt_); 
  NumericVector mat = as<Rcpp::NumericVector>(mat_);
  NumericMatrix sel = as<Rcpp::NumericMatrix>(sel_); 
  NumericVector cshare = as<Rcpp::NumericVector>(cshare_); 
  double SSBinit = as<double>(SSBinit_);
  double epsinit = as<double>(epsinit_);
  double B0 = as<double>(B0_);
  double hh = as<double>(hh_);
  double SPR = as<double>(SPR_);
  double sigmar = as<double>(sigmar_);
  double rho = as<double>(rho_); 
  double cfix = as<double>(cfix_);
  double alp = 4.*hh/(SPR*(1.-hh));
  double bet = (5.*hh-1.)/(B0*(1.-hh));
  double eps1,eps2,xisum,xsum;

  // memory allocation

  NumericMatrix N(na,ny);
  NumericVector SSB(ny);
  NumericMatrix X(nf,ny);

  // init popn

  for(a=0;a<na;a++) N(a,0) = Ninit(a);
  SSB(0) = SSBinit;

  // set up the fixed catch-derived exploitation rats

  for(f=0;f<nf;f++) {
      
    xsum = 0.; 
    for(a=0;a<na;a++) xsum += N(a,0)*sel(f,a)*wt(a);
    X(f,0) = cshare(f)*cfix/xsum;

  }

  // loop through the years

  eps1 = epsinit;
  for(y=1;y<ny;y++) {

    // recruitment

    if(y == 1) N(0,y) = alp*SSB(y-1)/(1.+bet*SSB(y-1))*exp(eps1);
    else {

      eps2 = rho*eps1+R::rnorm(-sigmar*sigmar/2.,sqrt(1.-rho*rho)*sigmar);
      N(0,y) = alp*SSB(y-1)/(1.+bet*SSB(y-1))*exp(eps2); 

    }

    eps1 = eps2; 

    // rest of the ages

    for(a=1;a<na;a++) {

      for(xisum=0.,f=0;f<nf;f++) xisum += X(f,y-1)*sel(f,a-1);
      xisum = xisum < 0.99 ? xisum : 0.99;
      N(a,y) = N(a-1,y-1)*exp(-M(a-1))*(1.-xisum);

    }

    // plus group

    for(xisum=0.,f=0;f<nf;f++) xisum += X(f,y-1)*sel(f,na-1); 
    xisum = xisum < 0.99 ? xisum : 0.99; 
    N(na-1,y) += N(na-1,y-1)*exp(-M(na-1))*(1.-xisum);

    // SSB

    SSB(y) = 0.;
    for(a=0;a<na;a++) SSB(y) += N(a,y)*wt(a)*mat(a); 

    // set up the fixed catch-derived exploitation rates

    for(f=0;f<nf;f++) {
      
      xsum = 0.; 
      for(a=0;a<na;a++) xsum += N(a,y)*sel(f,a)*wt(a);
      X(f,y) = cshare(f)*cfix/xsum;

    }

  }

  List res = Rcpp::List::create(Named("N")=N,Named("eps")=eps2);
 
  return Rcpp::wrap(res); 
}
