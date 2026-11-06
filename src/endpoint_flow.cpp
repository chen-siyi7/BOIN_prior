#include <Rcpp.h>
#include <algorithm>
#include <cstdint>
#include <functional>
#include <queue>
#include <vector>
using namespace Rcpp;
using Integer = int64_t;
using Wide = unsigned __int128;
struct Edge { int to, reverse; Integer capacity; };

// [[Rcpp::export]]
DataFrame endpoint_flow_audit(IntegerMatrix n, IntegerMatrix x,
    NumericVector coefficient, IntegerVector reference, IntegerVector action,
    IntegerVector capacity, IntegerMatrix pairs) {
  int k=n.nrow(), j=n.ncol();
  if(x.nrow()!=k || x.ncol()!=j || coefficient.size()!=k ||
     reference.size()!=k || action.size()!=k || capacity.size()!=k ||
     pairs.ncol()!=3) stop("Incompatible endpoint inputs.");

  for(int e=0; e<pairs.nrow(); ++e) {
    int endpoint=pairs(e,0), a=pairs(e,1)-1, b=pairs(e,2)-1, prefix=0;
    if((endpoint!=1 && endpoint!=j) || a<0 || a>=k || b<0 || b>=k ||
       reference[a]!=endpoint || reference[b]==endpoint ||
       capacity[a]<=0 || capacity[b]<=0) stop("Invalid endpoint pair.");
    for(int d=0; d<j; ++d) {
      if(n(a,d)!=n(b,d)) stop("An endpoint pair has different enrollment.");
      prefix += x(a,d)-x(b,d);
      if(d<j-1 && prefix<0) stop("The endpoint prefix inequality failed.");
    }
    if(prefix!=0) stop("An endpoint pair has different total DLTs.");
    double ca=coefficient[a], cb=coefficient[b];
    if(ca<1 || cb<1 || ca!=std::floor(ca) || cb!=std::floor(cb) ||
       ca>=9007199254740992.0 || cb>=9007199254740992.0)
      stop("The path coefficient is outside the exact integer range.");
    if(Wide(uint64_t(ca))*uint64_t(capacity[b]) !=
       Wide(uint64_t(cb))*uint64_t(capacity[a]))
      stop("The exact endpoint coefficient ratio failed.");
  }
  IntegerVector endpoint_out(2), losses(2), gains(2);
  NumericVector required_out(2), verified_out(2);
  LogicalVector passed(2,true);
  for(int pass=0; pass<2; ++pass) {
    int endpoint=pass==0 ? 1 : j, source=k, sink=k+1;
    std::vector<std::vector<Edge>> graph(k+2);
    std::vector<bool> bad(k,false),good(k,false);
    auto add=[&](int a,int b,Integer cap) {
      int ai=graph[a].size(), bi=graph[b].size();
      graph[a].push_back({b,bi,cap});
      graph[b].push_back({a,ai,0});
    };
    Integer demand=0;
    for(int s=0; s<k; ++s) {
      bad[s]=reference[s]==endpoint && action[s]!=endpoint;
      good[s]=reference[s]!=endpoint && action[s]==endpoint;
      if(bad[s]) {
        if(capacity[s]<=0) stop("An endpoint departure has no capacity.");
        demand+=capacity[s]; add(source,s,capacity[s]); ++losses[pass];
      }
      if(good[s]) {
        ++gains[pass];
        if(capacity[s]>0) add(s,sink,capacity[s]);
      }
    }
    for(int e=0; e<pairs.nrow(); ++e) {
      int a=pairs(e,1)-1,b=pairs(e,2)-1;
      if(pairs(e,0)==endpoint && bad[a] && good[b])
        add(a,b,std::min(capacity[a],capacity[b]));
    }
    Integer flow=0;
    while(true) {
      std::vector<int> level(k+2,-1),used(k+2,0);
      std::queue<int> queue;
      level[source]=0;queue.push(source);
      while(!queue.empty()) {
        int a=queue.front();queue.pop();
        for(const auto &e:graph[a]) if(e.capacity>0 && level[e.to]<0) {
          level[e.to]=level[a]+1;queue.push(e.to);
        }
      }
      if(level[sink]<0) break;
      std::function<Integer(int,Integer)> send=[&](int a,Integer amount) {
        if(a==sink) return amount;
        for(int &i=used[a]; i<int(graph[a].size()); ++i) {
          Edge &e=graph[a][i];
          if(e.capacity>0 && level[e.to]==level[a]+1) {
            Integer sent=send(e.to,std::min(amount,e.capacity));
            if(sent>0) {
              e.capacity-=sent;
              graph[e.to][e.reverse].capacity+=sent;
              return sent;
            }
          }
        }
        return Integer(0);
      };
      while(Integer added=send(source,demand-flow)) {
        flow+=added;
        if(flow==demand) break;
      }
      if(flow==demand) break;
    }
    if(flow!=demand) stop("The endpoint flow does not cover every departure.");
    endpoint_out[pass]=endpoint;
    required_out[pass]=double(demand);
    verified_out[pass]=double(flow);
  }
  return DataFrame::create(_["endpoint"]=endpoint_out,
      _["changed_reference_states"]=losses,_["new_endpoint_states"]=gains,
      _["required_integer_flow"]=required_out,
      _["verified_integer_flow"]=verified_out,_["pass"]=passed);
}
