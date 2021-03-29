#############################
#code associated with Table 1
#############################

#BB, BG, GB, GG
none<-c(0,0,0,0)
match<-c(1,1,-1,-1) 
silver<-c(-1,1,1,-1)
year2<-c(0,0,1,1)
year1<-c(0,1,0,1)
year<-c("00","01","10","11")
mat<-rbind(none,match,silver,year2,year1,year)
mat<-as.data.frame(t(mat));mat

mat$seed.neg<-c(0.8,1.25,0.8,1.25)
mat$seed.pos<-c(1.25,0.8,1.25,0.8)
mat$seed.none<-c(1,1,1,1)

##

nyear<-500; pfav<-0.5

pairs<-c("nonevsnone","nonevsmatch","nonevssilver","matchvsnone","matchvsmatch","matchvssilver","silvervsnone","silvervsmatch","silvervssilver")
pairs<-factor(pairs,levels=c("nonevsnone","nonevsmatch","nonevssilver","matchvsnone","matchvsmatch","matchvssilver","silvervsnone","silvervsmatch","silvervssilver"))


#Vmicr & Vocto
 
g1<-1; g2<-1
d1<-0; d2<-0
c<-1
 
#wet year
lambdai.g<-236; lambdaj.g<-924
aii.g<-0.099; aij.g<-0.083
ajj.g<-0.456; aji.g<-0.302

#dry year
lambdai.b<-153; lambdaj.b<-1127
aii.b<-0.082; aij.b<-0.002
ajj.b<-0.515; aji.b<-0.574

ri=NA;rj=NA


summ.sp1<-data.frame(actual=c(lambdai.g,lambdai.b,aii.g,aii.b,aij.g,aij.b,ri))
summ.sp2<-data.frame(actual=c(lambdaj.g,lambdaj.b,ajj.g,ajj.b,aji.g,aji.b,rj))


for (i in 1:length(pairs)){

  S.comp<-data.frame()
  
for (r in 1:100){
  
  me.typei<-sapply(strsplit(as.character(pairs[i]),"vs"), `[`, 1)
  me.typej<-sapply(strsplit(as.character(pairs[i]),"vs"), `[`, 2)
  
  #if loop is replicating each level 100 times, k loop is the amount of autocor, k=1 is always repeating, k=11 #is never repeating
    
      S<-data.frame(sp1=rep(NA,nyear),sp2=rep(NA,nyear),year=rep(NA,nyear),year.comb=rep(NA,nyear))  #matrix of N; each row is a year
      
      #species 1 is the invader, species 2 is the resident as set up here, 10 is a bit higher than carrying #capacity most of the time, 0.1 is 1/100th of that, to "win" for the results matrix, one species has to be #at least 100 times more abundant than the other after 500 time steps
      S[1:2,1]<-2
      S[1:2,2]<-2
      
      S[1,3]<-rbinom(1,1,prob=pfav) #random 1st year condition
      
      pswitch<-0.5
      
      bin<-rbinom(nyear,1,prob=pswitch) 
      
          for(j in 2:nyear){
          S[j,3]<-abs(S[j-1,3]-bin[j])
          S[j,4]<-paste(S[j-1,3],S[j,3],sep="")
          }  
      
      for(t in 3:nyear){
        
        tmp1<-S[t-1,3]
        tmp2<-S[t,3]
        tmp3<-S$year.comb[t]
        
        aii<-ifelse(tmp2==0,aii.b,aii.g); ajj<-ifelse(tmp2==0,ajj.b,ajj.g)
        aij<-ifelse(tmp2==0,aij.b,aij.g); aji<-ifelse(tmp2==0,aji.b,aji.g)  
        
        lambdai<-ifelse(tmp2==0,lambdai.b,lambdai.g);lambdaj<-ifelse(tmp2==0,lambdaj.b,lambdaj.g)
        
        tradeoff.typei<-"seed.pos"
        tradeoff.typej<-"seed.pos"
        
        #effect of current environment on # of seeds produced
        #trade-off removed in revision
        #si.o<-1/mat[mat$year1==tmp2,tradeoff.typei][1]
        #sj.o<-1/mat[mat$year1==tmp2,tradeoff.typej][1]
        si.o<-1
        sj.o<-1
        
        
        #effect of past environment on offspring success
        si.m<-mat[mat$year==tmp3,tradeoff.typei][1]
        sj.m<-mat[mat$year==tmp3,tradeoff.typej][1]
        
        formi<-as.numeric(as.character(mat[mat$year==S$year.comb[t],colnames(mat)==me.typei]))
        formj<-as.numeric(as.character(mat[mat$year==S$year.comb[t],colnames(mat)==me.typej]))
        
        Mi<-si.m^formi
        Mj<-sj.m^formj
        
        #lambda is modified by: si.o = negative inverse of seed size in current conditions + Mi = maternal effects due to seed size under different forms
        lambdai<-Mi*lambdai*si.o; lambdaj<-Mj*lambdaj*sj.o 
        
        S[t,1] <- as.integer((1-g1)*(1-d1)*S[t-1,1] + (g1*lambdai*S[t-1,1])/(c+aii*g1*S[t-1,1]+ aij*g2*S[t-1,2]))
        S[t,2] <- as.integer((1-g2)*(1-d2)*S[t-1,2] + (g2*lambdaj*S[t-1,2])/(c+ajj*g2*S[t-1,2]+ aji*g1*S[t-1,1]))

      }
      
      S$density.1<-c(2,S$sp1[1:nrow(S)-1])
      S$density.2<-c(2,S$sp2[1:nrow(S)-1])

      S$rep<-r
      S$time<-1:500
      
      S<-subset(S,time>2)   
      
      S$good<-S$year
      S$bad<-S$good+1
      S$bad<-ifelse(S$bad==2,0,1)
      
      S$gg<-ifelse(S$year.comb=="11",1,0)
      S$gb<-ifelse(S$year.comb=="10",1,0)
      S$bg<-ifelse(S$year.comb=="01",1,0)
      S$bb<-ifelse(S$year.comb=="00",1,0)      
      
      S.comp<-rbind(S.comp,S)
      
}

  
S<-S.comp
  
lm.sp1<-nls(sp1 ~ ((x1*density.1*good) + (x2*density.1*bad))/ (1 + y1*density.1*good + y2*density.1*bad + z1*density.2*good + z2*density.2*bad),data=S,start=list(x1=236,x2=153,y1=0.099,y2=0.088,z1=0.083,z2=0.002),nls.control(maxiter = 10000000, printEval = FALSE, warnOnly = TRUE) )
lm.sp2<-nls(sp2 ~ ((x1*density.2*good) + (x2*density.2*bad))/ (1 + y1*density.2*good + y2*density.2*bad + z1*density.1*good + z2*density.1*bad),data=S,start=list(x1=924,x2=1127,y1=0.456,y2=0.515,z1=0.302,z2=0.574),nls.control(maxiter = 10000000, printEval = FALSE, warnOnly = TRUE) )


summ.sp1[,1+i]<-round(c(summary(lm.sp1)$coefficients[,1],cor(S$sp1,fitted.values(lm.sp1))^2*100),3)
summ.sp2[,1+i]<-round(c(summary(lm.sp2)$coefficients[,1],cor(S$sp2,fitted.values(lm.sp2))^2*100),3)

}

summ.sp<-rbind(summ.sp1,summ.sp2)
colnames(summ.sp)<-as.character(c("actual",as.character(pairs)))
summ.sp$sp<-rep(c("vm","vo"),each=7)

write.csv(summ.sp,"summary of simulation.csv")


