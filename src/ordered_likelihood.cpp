#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
using namespace Rcpp;

static inline double score(int x,int n,double p) {
  if(n==0) return 0.0;
  if(p<=0.0) return x==0 ? 0.0 : -std::numeric_limits<double>::infinity();
  if(p>=1.0) return x==n ? 0.0 : -std::numeric_limits<double>::infinity();
  return x*std::log(p)+(n-x)*std::log1p(-p);
}

// [[Rcpp::export]]
List box_dp_extrema(IntegerMatrix N,IntegerMatrix X,NumericVector logA,
                    NumericMatrix low,NumericMatrix up,IntegerVector states) {
  const int J=N.ncol(),K=N.nrow(),R=low.nrow(),S=states.size();
  if(X.nrow()!=K || X.ncol()!=J || logA.size()!=K ||
     up.nrow()!=R || low.ncol()!=J || up.ncol()!=J)
    stop("dimension mismatch");
  NumericMatrix U(S,R),L(S,R);
  const double inf=std::numeric_limits<double>::infinity();
  for(int r=0;r<R;++r) {
    for(int h=0;h<S;++h) {
      const int s=states[h]-1;
      if(s<0 || s>=K) stop("state out of range");
      std::vector<double> v;
      v.reserve(2*J+J*(J+1)/2);
      for(int k=0;k<J;++k) {
        v.push_back(low(r,k));v.push_back(up(r,k));
      }
      for(int i=0;i<J;++i) {
        int nx=0,nn=0;
        for(int k=i;k<J;++k) {
          nx+=X(s,k);nn+=N(s,k);
          if(nn>0) v.push_back((double)nx/nn);
        }
      }
      std::sort(v.begin(),v.end());
      v.erase(std::unique(v.begin(),v.end()),v.end());
      const int V=v.size();
      std::vector<double> maxprev(V,-inf),minprev(V,inf),
                          maxcur(V,-inf),mincur(V,inf);
      std::vector<int> prevValid(V,0),curValid(V,0);
      for(int k=0;k<J;++k) {
        double bestmax=-inf,bestmin=inf;bool seen=false;
        for(int t=0;t<V;++t) {
          if(k>0 && prevValid[t]) {
            seen=true;bestmax=std::max(bestmax,maxprev[t]);
            bestmin=std::min(bestmin,minprev[t]);
          }
          if(v[t]>=low(r,k)-1e-14 && v[t]<=up(r,k)+1e-14 && (k==0 || seen)) {
            curValid[t]=1;
            const double f=score(X(s,k),N(s,k),v[t]);
            maxcur[t]=(k==0 ? f : bestmax+f);
            mincur[t]=(k==0 ? f : bestmin+f);
          } else {
            curValid[t]=0;maxcur[t]=-inf;mincur[t]=inf;
          }
        }
        maxprev.swap(maxcur);minprev.swap(mincur);prevValid.swap(curValid);
      }
      double maxval=-inf,minval=inf;bool any=false;
      for(int t=0;t<V;++t)if(prevValid[t]) {any=true;maxval=std::max(maxval,maxprev[t]);minval=std::min(minval,minprev[t]);}
      if(!any)stop("empty monotone box");
      if(maxval==-inf || minval==inf) stop("empty monotone box");
      U(h,r)=std::exp(logA[s]+maxval);
      L(h,r)=std::exp(logA[s]+minval);
    }
  }
  return List::create(_["U"]=U,_["L"]=L);
}
