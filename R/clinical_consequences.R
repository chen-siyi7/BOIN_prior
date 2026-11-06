out <- output
o <- clinical$model
J <- 3L
b <- o$reference
K <- length(b)
fits <- clinical
actions <- c(list(BOIN = b), fits$policies)
P <- o$planning
M <- exact_mass(o$catalog, P)
n <- as.matrix(o$catalog[paste0("n", 1:J)])
x <- as.matrix(o$catalog[paste0("x", 1:J)])
priors<-fits$priors
clinical_metrics <- lapply(actions, function(action) {
 selected <- vapply(1:J, function(j)colSums(M[action==j,,drop=FALSE]), numeric(nrow(P)))
 cbind(PCS=rowSums(selected*(col(selected)==o$labels)),
   A=rowSums(selected*(abs(P-.35)<=.05)),R=rowSums(selected*(P>.40)),
   B=rowSums(selected*(P<.30)),N=colSums(M[action<0,,drop=FALSE]),DLT=rowSums(selected*P))
})
historical_rows <- list()
for(prior in names(priors))for(rule in names(actions)) {
 w<-priors[[prior]];value<-clinical_metrics[[rule]];delta<-value-clinical_metrics$BOIN
 bracket<-o$labels>0;selected<-sum(w*(1-value[,'N']))
 historical_rows[[length(historical_rows)+1L]]<-data.frame(prior=prior,rule=rule,
   PCS=100*sum(w*value[,'PCS'])/sum(w[bracket]),
   PCS_gain=100*sum(w*delta[,'PCS'])/sum(w[bracket]),
   within=100*sum(w*value[,'A']),toxic=100*sum(w*value[,'R']),
   under=100*sum(w*value[,'B']),none=100*sum(w*value[,'N']),
   utility_gain=100*sum(w*(delta[,'A']-delta[,'R']-.25*delta[,'B']-.5*delta[,'N'])),
   mean_DLT=100*sum(w*value[,'DLT'])/selected,cohort20_DLT=20*sum(w*value[,'DLT'])/selected)
}
historical<-do.call(rbind,historical_rows)
compare_numeric_table(historical,read_reference_csv('historical_summary.csv'),c('prior','rule'),1e-8)
write.csv(historical,file.path(out,'historical_summary.csv'),row.names=FALSE)



expansion<-function(p) {
 ep<-ed<-list()
 key<-function(n,x)paste(n,x,sep='/')
 for(nn in seq(18,0,by=-3))for(xx in 0:nn) {
  if(nn==18||(nn>=3&&pbeta(.35,xx+1,nn-xx+1,lower.tail=FALSE)>.95)) {
   ep[[key(nn,xx)]]<-rep(0,length(p));ed[[key(nn,xx)]]<-rep(0,length(p))
  }else {
   vp<-rep(3,length(p));vd<-3*p
   for(y in 0:3) {
    q<-dbinom(y,3,p)
    vp<-vp+q*ep[[key(nn+3,xx+y)]];vd<-vd+q*ed[[key(nn+3,xx+y)]]
   }
   ep[[key(nn,xx)]]<-vp;ed[[key(nn,xx)]]<-vd
  }
 }
 list(patients=ep,DLT=ed)
}
cont<-lapply(1:J,function(j)expansion(P[,j]))
result<-observed<-list();baseline<-list()
for(rule in names(actions)) {
 a<-actions[[rule]];expn<-expd<-dist<-over<-rep(0,nrow(P));change<-rep(0,nrow(P))
 for(s in which(a>0)) {
  j<-a[s];key<-paste(n[s,j],x[s,j],sep='/');ms<-M[s,]
  expn<-expn+ms*cont[[j]]$patients[[key]];expd<-expd+ms*cont[[j]]$DLT[[key]]
  dist<-dist+ms*abs(P[,j]-.35);over<-over+ms*pmax(P[,j]-.35,0)
 }
 if(rule!='BOIN')change<-colSums(M[a!=b,,drop=FALSE])
 for(prior in names(priors)) {
  w<-priors[[prior]];qrec<-sum(w*colSums(M[a>0,,drop=FALSE]));rr<-data.frame(rule=rule,prior=prior,change_percent=100*sum(w*change),
   expected_expansion_patients=sum(w*expn),expected_expansion_DLT=sum(w*expd),
   target_distance_pp=100*sum(w*dist)/qrec,above_target_excess_pp=100*sum(w*over)/qrec)
  result[[length(result)+1]]<-rr
 }
 for(arm in c('MN','DD')) {
  at<-which(n[,1]==0&n[,2]==6&n[,3]==12&x[,1]==0&x[,2]==1&x[,3]==if(arm=='MN')5 else 3)
  stopifnot(length(at)==1)
  observed[[length(observed)+1]]<-data.frame(source=arm,rule=rule,BOIN=b[at],recommendation=a[at],dose_mg=c(150,200,250)[a[at]])
 }
}
result<-do.call(rbind,result)
for(prior in names(priors)) {
 ii<-result$prior==prior;bb<-result[ii&result$rule=='BOIN',]
 result$patient_change[ii]<-result$expected_expansion_patients[ii]-bb$expected_expansion_patients
 result$DLT_change[ii]<-result$expected_expansion_DLT[ii]-bb$expected_expansion_DLT
 result$distance_reduction_pp[ii]<-bb$target_distance_pp-result$target_distance_pp[ii]
 result$above_excess_reduction_pp[ii]<-bb$above_target_excess_pp-result$above_target_excess_pp[ii]
}
stopifnot(all(result$expected_expansion_patients>=0),all(result$expected_expansion_DLT>=0),
 all(result$expected_expansion_DLT<=result$expected_expansion_patients+1e-10))
write.csv(result,file.path(out,'protocol_consequences.csv'),row.names=FALSE)
write.csv(do.call(rbind,observed),file.path(out,'published_counts_lookup.csv'),row.names=FALSE)

g<-if(file.exists(file.path(out,'clinical_exact_gains.json')))
 jsonlite::read_json(file.path(out,'clinical_exact_gains.json'),simplifyVector=TRUE)else
 read_reference_json('cost_audited_gains.json')
ceil<-read_reference_json('clinical_gain_bounds.json')$results
if(file.exists(file.path(out,'gain_ceiling_checks.csv'))) {
 verified<-read.csv(file.path(out,'gain_ceiling_checks.csv'))
 ceil$directed_upper_pp<-verified$upper_pp[match(ceil$rule,verified$rule)]
}
quad<-if(file.exists(file.path(out,'continuous_prior_gains.csv')))
 read.csv(file.path(out,'continuous_prior_gains.csv'))else read_reference_csv('continuous_prior_gains.csv')
quad<-quad[quad$order==256,]
settings<-list(
 MN_fixed=list(rule='primary_MN',priors='MN_0.5',nominal=TRUE,bound='primary_MN'),
 MN_flexible=list(rule='expanded_MN',priors='MN_0.5',nominal=TRUE,bound='expanded_MN'),
 DD_fixed=list(rule='primary_DD',priors='DD_0.5',nominal=TRUE,bound='primary_DD'),
 DD_flexible=list(rule='expanded_DD',priors='DD_0.5',nominal=TRUE,bound='expanded_DD'),
 MN_preferences=list(rule='uncertain_MN',priors=grep('^MN_',names(priors),value=TRUE),nominal=FALSE,bound='expanded_MN'),
 sources_preferences=list(rule='uncertain_sources',priors=names(priors),nominal=FALSE,bound='expanded_DD'))
thresholds<-frontier<-review<-list()
for(key in names(settings)) {
 s<-settings[[key]];z<-g$results[g$results$rule==s$rule&g$results$prior%in%s$priors,]
 qq<-quad[quad$rule==s$rule&quad$prior%in%s$priors,]
 if(s$nominal) {
  z<-z[as.numeric(z$lambda)==1&as.numeric(z$underdose)==.25&as.numeric(z$none)==.5,]
  qq<-qq[qq$lambda==1&qq$underdose==.25&qq$none==.5,]
 }
 lower<-min(as.numeric(z$lower_pp));up<-ceil$directed_upper_pp[ceil$rule==s$bound]
 absolute<-historical[historical$rule==s$rule&historical$prior%in%s$priors,]
 thresholds[[key]]<-data.frame(specification=key,rule=s$rule,gain_lower_pp=lower,gain_upper_pp=up,
  minimum_within_percent=min(absolute$within),maximum_toxic_percent=max(absolute$toxic),maximum_none_percent=max(absolute$none),
  continuous_gain=min(qq$gain_pp),min_change_percent=min(result$change_percent[result$rule==s$rule&result$prior%in%s$priors]),
  max_change_percent=max(result$change_percent[result$rule==s$rule&result$prior%in%s$priors]))
 for(gamma in c(0,.001,.0025,.005,.01,.025,.05,.1,.2,.5))frontier[[length(frontier)+1]]<-data.frame(specification=key,
  required_gain_pp=gamma,status=if(lower>=gamma)'gain requirement met'else if(as.numeric(up)<gamma)'requirement unattainable'else'not resolved')
 for(cost in c(0,1,5,10,25,50,100)) {
  pp<-result$change_percent[match(paste(s$rule,z$prior),paste(result$rule,result$prior))]/100
  
  review[[length(review)+1]]<-data.frame(specification=key,review_cost_pp=cost,
   available_fixed_burden_pp=min(as.numeric(z$lower_pp)-cost*pp))
 }
}
write.csv(do.call(rbind,thresholds),file.path(out,'benefit_thresholds.csv'),row.names=FALSE)
write.csv(do.call(rbind,frontier),file.path(out,'adoption_frontier.csv'),row.names=FALSE)
write.csv(do.call(rbind,review),file.path(out,'review_burden_frontier.csv'),row.names=FALSE)
saveRDS(list(settings=settings,results=result),file.path(out,'clinical_analysis.rds'))
cat('Protocol-linked consequences and clinical benefit thresholds complete.\n')

cases<-list()
for(pp in c(.05,.35,.8)) {
 cc<-expansion(pp)
 for(nn in seq(0,18,by=3))for(xx in 0:nn) {
  key<-paste(nn,xx,sep='/')
  cases[[length(cases)+1]]<-data.frame(p=pp,n=nn,x=xx,
    patients=cc$patients[[key]],DLT=cc$DLT[[key]])
 }
}
write.csv(do.call(rbind,cases),file.path(out,'continuation_cases.csv'),row.names=FALSE)
