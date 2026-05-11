# Pointwise Overlap Indicator for Propensity Scores based on Nethery et al. (2019)
# -------------------------------------------------
# For each unit i, decides whether unit i lies in the region of overlap
# (RO) on the propensity-score scale. A unit at PS value k is "in RO" if,
# in *both* treatment arms, there exists a window of b + 1 consecutive
# PS values containing k whose range is smaller than a. In words: each
# arm has at least b nearby observations packed within distance a of k.
#
# Args:
#   ps : numeric vector of estimated propensity scores, length n
#   E  : 0/1 treatment indicator, length n
#   a  : tolerance on the PS distance defining "nearby"
#   b  : minimum number of neighbours required on each side within
#        that tolerance
#
# Returns:
#   Integer vector of length n with 1 = in the overlap region, 0 = not.



pw_overlap<-function(ps,E,a,b){
	
	# Sorted PS within each treatment arm
  ps1<-ps[which(E==1)]
  ps0<-ps[which(E==0)]
  ps1<-ps1[order(ps1)]
  ps0<-ps0[order(ps0)]
  
  RO<-rep(0,length(ps))
  for (k in ps){
    cte<-0   # how many treatment arms satisfy the local-density condition
    for (e in 0:1){
      cth<-0
      temp<-get(paste('ps',e,sep=''))
      fooless<-temp[which(temp<k)]
      foomore<-temp[which(temp>k)]
      
      if (k %in% temp){
      	# k itself is in this arm: take b neighbours on each side,
      	# padding with +/- Inf when fewer than b exist
        if (length(fooless)>=b & length(foomore)>=b){
          allcheck<-c(fooless[(length(fooless)-(b-1)):length(fooless)],k,foomore[1:b])
        } else if (length(fooless)>=b & length(foomore)<b){
          allcheck<-c(fooless[(length(fooless)-(b-1)):length(fooless)],k,foomore,rep(Inf,b-length(foomore)))
        } else if (length(fooless)<b & length(foomore)>=b){
          allcheck<-c(rep(-Inf,b-length(fooless)),fooless,k,foomore[1:b])
        } else{
          allcheck<-c(rep(-Inf,b-length(fooless)),fooless,k,foomore,rep(Inf,b-length(foomore)))
        }
        
      	# Slide a window of length b+1 across allcheck; count windows
      	# whose span is below a
        for (h in 1:(b+1)){
          if (abs(allcheck[h]-allcheck[h+b])<a) cth<-cth+1
        }
      } else{
      	# k is not in this arm: need b+1 neighbours on each side
        if (length(fooless)>=(b+1) & length(foomore)>=(b+1)){
          allcheck<-c(fooless[(length(fooless)-b):length(fooless)],k,foomore[1:(b+1)])
        } else if (length(fooless)>=(b+1) & length(foomore)<(b+1)){
          allcheck<-c(fooless[(length(fooless)-b):length(fooless)],k,foomore,rep(Inf,(b+1)-length(foomore)))
        } else if (length(fooless)<(b+1) & length(foomore)>=(b+1)){
          allcheck<-c(rep(-Inf,(b+1)-length(fooless)),fooless,k,foomore[1:(b+1)])
        } else{
          allcheck<-c(rep(-Inf,(b+1)-length(fooless)),fooless,k,foomore,rep(Inf,(b+1)-length(foomore)))
        }
        
        for (h in 1:(b+2)){
          if (abs(allcheck[h]-allcheck[h+b+1])<a) cth<-cth+1
        }
      }
      if (cth>0) cte<-cte+1
    }
    
    # Mark k as in RO only if both arms satisfied the condition
    if (cte==2) RO[which(ps==k)]<-1
  }
  
  return(RO)
}