#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace Rcpp;
static double la(double a,double b) {
 if(!std::isfinite(a))return b;if(!std::isfinite(b))return a;
 double m=std::max(a,b);return m+std::log(std::exp(a-m)+std::exp(b-m));
}
static double bm(double p,double p0,double q,bool reverse) {
 if(reverse) {
  if(p==0||p==1)return R_PosInf;
  return la(q*std::log(p0)+(1-q)*std::log(p),
            q*std::log1p(-p0)+(1-q)*std::log1p(-p));
 }
 double a=p==0?R_NegInf:q*std::log(p)+(1-q)*std::log(p0);
 double b=p==1?R_NegInf:q*std::log1p(-p)+(1-q)*std::log1p(-p0);
 return la(a,b);
}



static double moment(NumericMatrix lo,NumericMatrix hi,NumericMatrix center,
 NumericVector nmax,int r,double q,bool reverse) {
 int J=lo.ncol();std::vector<double> v;
 for(int j=0;j<J;++j){v.push_back(lo(r,j));v.push_back(hi(r,j));}
 std::sort(v.begin(),v.end());v.erase(std::unique(v.begin(),v.end()),v.end());
 int V=v.size();std::vector<double> prev(V,R_NegInf),cur(V,R_NegInf);
 for(int j=0;j<J;++j){double best=R_NegInf;
  for(int t=0;t<V;++t){if(j>0)best=std::max(best,prev[t]);
   if(v[t]>=lo(r,j)&&v[t]<=hi(r,j)&&(j==0||best>R_NegInf))
    cur[t]=(j==0?0:best)+nmax[j]*bm(v[t],center(r,j),q,reverse);
   else cur[t]=R_NegInf;
  }prev.swap(cur);
 }return std::max(0.,*std::max_element(prev.begin(),prev.end()));
}
// [[Rcpp::export]]
List cse_bounds(NumericMatrix low,NumericMatrix up,NumericMatrix center,
 NumericVector nmax,NumericVector orders,NumericMatrix plus,NumericMatrix minus) {
 int R=low.nrow(),P=plus.ncol();NumericMatrix out(R,P);
 for(int r=0;r<R;++r){
  std::vector<double> ub(P,1.),lb(P,0.),aff(P,1.);
  for(double q:orders){if(q<=1)continue;
   double mf=moment(low,up,center,nmax,r,q,false);
   double mr=moment(low,up,center,nmax,r,q,true);
   for(int m=0;m<P;++m){double a=plus(r,m),b=minus(r,m);
    if(a==0)ub[m]=0;else ub[m]=std::min(ub[m],std::exp((1-1/q)*std::log(a)+mf/q));
    double complement=1-std::exp((1-1/q)*std::log1p(-b)+mf/q);
    double inverse=b==0||!std::isfinite(mr)?0:std::exp(q/(q-1)*std::log(b)-mr/(q-1));
    lb[m]=std::max(lb[m],std::max(complement,inverse));
    double f=(a-b+1)/2;
    aff[m]=std::min(aff[m],-1+2*std::exp((1-1/q)*std::log(f)+mf/q));
   }
  }
  for(int m=0;m<P;++m)out(r,m)=std::min(aff[m],ub[m]-lb[m]);
 }
 return List::create(_["upper"]=out);
}
