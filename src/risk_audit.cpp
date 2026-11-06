#include <Rcpp.h>
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <limits>
#include <numeric>
#include <vector>
using namespace Rcpp;





using u128 = unsigned __int128;
using i128 = __int128;
static constexpr uint64_t S = uint64_t(1) << 63;

struct Rat { uint64_t n, d; };
static bool less_rat(const Rat &a, const Rat &b) {
  return u128(a.n) * b.d < u128(b.n) * a.d;
}
static bool eq_rat(const Rat &a, const Rat &b) {
  return u128(a.n) * b.d == u128(b.n) * a.d;
}



static bool sum_le_target(const Rat &a,const Rat &b,int tn,int td) {
  uint64_t D=std::max(a.d,b.d);
  return (u128(a.n)*(D/a.d)+u128(b.n)*(D/b.d))*td<=u128(tn)*D;
}
static bool sum_ge_target(const Rat &a,const Rat &b,int tn,int td) {
  uint64_t D=std::max(a.d,b.d);
  return (u128(a.n)*(D/a.d)+u128(b.n)*(D/b.d))*td>=u128(tn)*D;
}
static Rat rat_double(double x) {
  if (!std::isfinite(x) || x < 0 || x > 1) stop("bad box endpoint");
  if (x == 0) return {0,1};
  if (x == 1) return {1,1};
  uint64_t bits;
  std::memcpy(&bits, &x, sizeof(bits));
  const int raw = int((bits >> 52) & 2047);
  if (raw == 0) stop("subnormal box endpoint not supported");
  uint64_t n = (uint64_t(1) << 52) | (bits & ((uint64_t(1) << 52)-1));
  int shift = 1075 - raw;
  if (shift < 0 || shift > 63) stop("box endpoint dyadic denominator too large");
  while (shift && (n & 1u) == 0u) { n >>= 1; --shift; }
  if(n>(uint64_t(1)<<53)) stop("endpoint numerator exceeds checked range");
  return {n, uint64_t(1) << shift};
}
static uint64_t qfloor(const Rat &p) {
  return uint64_t((u128(p.n)*S)/p.d);
}
static uint64_t qceil(const Rat &p) {
  const u128 x=u128(p.n)*S;
  return uint64_t(x/p.d + (x%p.d != 0));
}
static uint64_t mlo(uint64_t a,uint64_t b) {
  return uint64_t((u128(a)*b)>>63);
}
static uint64_t mhi(uint64_t a,uint64_t b) {
  const u128 x=u128(a)*b;
  return uint64_t((x+S-1)>>63);
}
struct Iv {uint64_t lo,hi;};
static Iv likelihood_factor(int n,int x,const Rat &p) {
  uint64_t pl=qfloor(p),pu=qceil(p),ql=S-pu,qu=S-pl;
  Iv out{S,S};
  for(int k=0;k<x;++k) {out.lo=mlo(out.lo,pl);out.hi=mhi(out.hi,pu);}
  for(int k=x;k<n;++k) {out.lo=mlo(out.lo,ql);out.hi=mhi(out.hi,qu);}
  return out;
}

// [[Rcpp::export]]
DataFrame fixedpoint_interval_audit(IntegerMatrix N,IntegerMatrix X,
    NumericVector coeff,NumericMatrix low,NumericMatrix up,
    LogicalMatrix feasible,IntegerVector states,IntegerMatrix h,
    int target_num,int target_den,int cap_num=20,int cap_den=1000,int max_cells=0,int endpoint_num=0) {
  const int K=N.nrow(),J=N.ncol(),R=low.nrow(),P=h.ncol();
  const int Sstates=states.size(),P0=5*J-2;
  if(target_num<1 || target_den<1 || target_den>100 || target_num>target_den || cap_num<0 || cap_den<1) stop("bad rational target or cap");
  if(J<2 || X.nrow()!=K || X.ncol()!=J || coeff.size()!=K ||
     up.nrow()!=R || low.ncol()!=J || up.ncol()!=J ||
     feasible.nrow()!=R || feasible.ncol()!=J ||
     h.nrow()!=Sstates || P%P0!=0) stop("dimension mismatch");
  long double coefficient_total=0;
  for(int ss=0;ss<Sstates;++ss) {
    double a=coeff[states[ss]-1];
    if(a<1 || a!=std::floor(a) || a>=9007199254740992.0) stop("coefficient range");
    coefficient_total+=a;
  }
  if(coefficient_total>=9007199254740992.0L || cap_den>10000) stop("audit accumulator range");
  const int Ruse=(max_cells>0 ? std::min(R,max_cells):R);
  std::vector<Rat> elo(Ruse*J),eup(Ruse*J),endpoints;
  endpoints.reserve(2*Ruse*J);
  for(int r=0;r<Ruse;++r) for(int j=0;j<J;++j) {
    elo[r*J+j]=rat_double(low(r,j));eup[r*J+j]=rat_double(up(r,j));
    if(less_rat(eup[r*J+j],elo[r*J+j])) stop("inverted box");
    endpoints.push_back(elo[r*J+j]); endpoints.push_back(eup[r*J+j]);
  }



  for(int r=0;r<Ruse;++r) for(int j=0;j<J;++j) {
    bool possible=true;
    if(j>0) possible &= sum_le_target(elo[r*J+j-1],elo[r*J+j],target_num,target_den);
    if(j<J-1) possible &= sum_ge_target(eup[r*J+j],eup[r*J+j+1],target_num,target_den);
    if(possible && !feasible(r,j)) stop("stored label mask omits necessary cell");
  }
  std::sort(endpoints.begin(),endpoints.end(),less_rat);
  endpoints.erase(std::unique(endpoints.begin(),endpoints.end(),eq_rat),
                  endpoints.end());
  std::vector<int> elid(Ruse*J),euid(Ruse*J);
  auto index_of=[](const std::vector<Rat>&v,const Rat&q) {
    auto it=std::lower_bound(v.begin(),v.end(),q,less_rat);
    if(it==v.end() || !eq_rat(*it,q)) stop("candidate index not found");
    return int(it-v.begin());
  };
  for(int r=0;r<Ruse;++r) for(int j=0;j<J;++j) {
    elid[r*J+j]=index_of(endpoints,elo[r*J+j]);
    euid[r*J+j]=index_of(endpoints,eup[r*J+j]);
  }
  std::vector<i128> sum(size_t(P)*Ruse,0);
  for(int ss=0;ss<Sstates;++ss) {
    const int s=states[ss]-1;
    if(s<0 || s>=K || coeff[s]<1) stop("bad state or coefficient");
    std::vector<std::pair<int,int>> active;
    for(int m=0;m<P;++m) {
      if(h(ss,m)<-1 || h(ss,m)>1) stop("contrast coefficient not signed unit");
      if(h(ss,m)) active.emplace_back(m,h(ss,m));
    }
    if(active.empty()) continue;
    std::vector<Rat> v=endpoints;
    for(int i=0;i<J;++i) {
      int nx=0,nn=0;
      for(int k=i;k<J;++k) {
        nx+=X(s,k);nn+=N(s,k);
        if(nn>0) v.push_back({uint64_t(nx),uint64_t(nn)});
      }
    }
    std::sort(v.begin(),v.end(),less_rat);
    v.erase(std::unique(v.begin(),v.end(),eq_rat),v.end());
    const int V=v.size();
    std::vector<Iv> f(size_t(J)*V);
    for(int i=0;i<J;++i) for(int t=0;t<V;++t) {
      if(X(s,i)<0 || N(s,i)<X(s,i)) stop("invalid counts");
      f[i*V+t]=likelihood_factor(N(s,i),X(s,i),v[t]);
    }
    std::vector<int> endpoint_map(endpoints.size());
    for(size_t t=0;t<endpoints.size();++t)
      endpoint_map[t]=index_of(v,endpoints[t]);
    std::vector<int> pool;
    for(int i=0;i<J;++i) {
      int nx=0,nn=0;
      for(int k=i;k<J;++k) {
        nx+=X(s,k);nn+=N(s,k);
        if(nn>0) pool.push_back(index_of(v,{uint64_t(nx),uint64_t(nn)}));
      }
    }
    std::sort(pool.begin(),pool.end());
    pool.erase(std::unique(pool.begin(),pool.end()),pool.end());
    for(int r=0;r<Ruse;++r) {
      std::vector<int> ids=pool;
      ids.reserve(pool.size()+2*J);
      for(int j=0;j<J;++j) {
        ids.push_back(endpoint_map[elid[r*J+j]]);
        ids.push_back(endpoint_map[euid[r*J+j]]);
      }
      std::sort(ids.begin(),ids.end());
      ids.erase(std::unique(ids.begin(),ids.end()),ids.end());
      const int T=ids.size();
      std::vector<uint64_t> prevU(T),prevL(T),curU(T),curL(T);
      std::vector<uint8_t> prevValid(T,0),curValid(T,0);
      for(int i=0;i<J;++i) {
        uint64_t bestU=0,bestL=S+1;
        bool seen=false;
        for(int t=0;t<T;++t) {
          if(i>0 && prevValid[t]) {
            seen=true;
            bestU=std::max(bestU,prevU[t]);
            bestL=std::min(bestL,prevL[t]);
          }
          const int id=ids[t];
          const bool within=!less_rat(v[id],elo[r*J+i]) &&
                            !less_rat(eup[r*J+i],v[id]);
          if(within && (i==0 || seen)) {
            const Iv &g=f[i*V+id];
            curU[t]=(i==0 ? g.hi:mhi(bestU,g.hi));
            curL[t]=(i==0 ? g.lo:mlo(bestL,g.lo));
            curValid[t]=1;
          } else curValid[t]=0;
        }
        prevU.swap(curU);prevL.swap(curL);prevValid.swap(curValid);
      }
      uint64_t upper=0,lower=S+1;
      bool seen=false;
      for(int t=0;t<T;++t) if(prevValid[t]) {
        seen=true;upper=std::max(upper,prevU[t]);
        lower=std::min(lower,prevL[t]);
      }
      if(!seen) stop("empty ordered cell");
      const i128 U=i128(coeff[s])*upper,L=i128(coeff[s])*lower;
      for(auto e:active) sum[size_t(e.first)*Ruse+r]+=(e.second>0?U:-L);
    }
  }
  CharacterVector policy(P),metric(P);
  IntegerVector argmax(P);
  NumericVector bound_pp(P),slack_probability(P);
  LogicalVector passes(P);
  for(int m=0;m<P;++m) {
    const int local=m%P0;
    policy[m]=std::to_string(m/P0+1);
    if(local<J) metric[m]="C"+std::to_string(local+1);
    else if(local<2*J) metric[m]="O"+std::to_string(local-J+1);
    else if(local<3*J) metric[m]="B"+std::to_string(local-2*J+1);
    else if(local<4*J-1) metric[m]="Tp"+std::to_string(local-3*J+1);
    else metric[m]="Tm"+std::to_string(local-(4*J-1)+1);
    i128 mx=std::numeric_limits<int64_t>::min();int who=-1;
    for(int r=0;r<Ruse;++r) {
      if(local<3*J && !feasible(r,local%J)) continue;
      if(who<0 || sum[size_t(m)*Ruse+r]>mx) {
        mx=sum[size_t(m)*Ruse+r];who=r;
      }
    }
    if(who<0) stop("no admissible label cell");
    const int numerator=(local==0 || local==J-1 ? endpoint_num:cap_num);
    const i128 cap_scaled=i128(numerator)*S;
    passes[m]=(mx*cap_den<=cap_scaled);
    argmax[m]=who+1;
    bound_pp[m]=100.0*double((long double)mx/(long double)S);
    slack_probability[m]=double((long double)(cap_scaled-mx*cap_den)/
                                 ((long double)cap_den*(long double)S));
  }
  return DataFrame::create(_["policy"]=policy,_["metric"]=metric,
    _["argmax_cell"]=argmax,_["bound_pp"]=bound_pp,
    _["slack_probability"]=slack_probability,_["pass"]=passes);
}
