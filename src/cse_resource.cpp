#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace Rcpp;
static double addlog(double a,double b) {
 if(!std::isfinite(a))return b;if(!std::isfinite(b))return a;
 double m=std::max(a,b);return m+std::log(std::exp(a-m)+std::exp(b-m));
}
static double one(double p,double q,double r,bool reverse) {
 if(reverse) {
  if(p==0||p==1)return R_PosInf;
  return addlog(r*std::log(q)+(1-r)*std::log(p),
                r*std::log1p(-q)+(1-r)*std::log1p(-p));
 }
 return addlog(p==0?R_NegInf:r*std::log(p)+(1-r)*std::log(q),
               p==1?R_NegInf:r*std::log1p(-p)+(1-r)*std::log1p(-q));
}
static double resource_moment(NumericMatrix low,NumericMatrix up,NumericMatrix center,
 NumericMatrix enroll,int cell,double r,bool reverse) {
 int J=low.ncol();std::vector<double> psi(J);
 for(int j=0;j<J;++j)psi[j]=std::max(0.,std::max(one(low(cell,j),center(cell,j),r,reverse),one(up(cell,j),center(cell,j),r,reverse)));
 double out=0;
 for(int k=0;k<enroll.nrow();++k) {
  double value=0;
  for(int j=0;j<J;++j)if(enroll(k,j)>0)value+=enroll(k,j)*psi[j];
  out=std::max(out,value);
 };return out;
}
// [[Rcpp::export]]
List cse_resource_bounds(NumericMatrix low,NumericMatrix up,NumericMatrix center,
 NumericMatrix enroll,NumericVector orders,NumericMatrix plus,NumericMatrix minus) {
 int R=low.nrow(),P=plus.ncol();NumericMatrix out(R,P);
 for(int cell=0;cell<R;++cell) {
  std::vector<double> ub(P,1.),lb(P,0.),aff(P,1.);
  for(double r:orders) {
   if(r<=1)continue;
   double mf=resource_moment(low,up,center,enroll,cell,r,false);
   double mr=resource_moment(low,up,center,enroll,cell,r,true);
   for(int m=0;m<P;++m) {
    double a=plus(cell,m),b=minus(cell,m);
    if(a==0)ub[m]=0;else ub[m]=std::min(ub[m],std::exp((1-1/r)*std::log(a)+mf/r));
    double complement=1-std::exp((1-1/r)*std::log1p(-b)+mf/r);
    double inverse=b==0||!std::isfinite(mr)?0:std::exp(r/(r-1)*std::log(b)-mr/(r-1));
    lb[m]=std::max(lb[m],std::max(complement,inverse));
    double f=(a-b+1)/2;
    aff[m]=std::min(aff[m],-1+2*std::exp((1-1/r)*std::log(f)+mf/r));
   }
  }
  for(int m=0;m<P;++m)out(cell,m)=std::min(aff[m],ub[m]-lb[m]);
 };return List::create(_["upper"]=out);
}
// [[Rcpp::export]]
NumericVector cse_resource_moments(NumericMatrix low,NumericMatrix up,NumericMatrix center,
 NumericMatrix enroll,double order) {
 return NumericVector::create(resource_moment(low,up,center,enroll,0,order,false),
                               resource_moment(low,up,center,enroll,0,order,true));
}
