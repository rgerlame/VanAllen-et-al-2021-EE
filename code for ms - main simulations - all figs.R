
require(ggplot2); require(cowplot)

##############################
#setting up parameter values
##############################

g1<-1; g2<-1
d1<-0; d2<-0
c<-1

lambdai.g<-60; lambdaj.g<-60
lambdai.b<-10; lambdaj.b<-10


#list of values to simulate different ecological scenarios; species j is dominant in scenario 

aii.g<-c(1,1,1.25,1.25,1.25)
ajj.g<-c(1,1,0.8,1.25,1)

aii.b<-c(1,1,1.25,1.25,1.25)
ajj.b<-c(1,1,0.8,1.25,1)


aij.g<-c(0,1,0.8,0.8,1)
aji.g<-c(0,1,1.25,0.8,0.8)

aij.b<-c(0,1,0.8,0.8,1)
aji.b<-c(0,1,1.25,0.8,0.8)


scenario<-c("no comp","neutral","fd","nd","fd+nd")

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


#simulation

nyear<-100; pfav<-0.5

pairs<-c("nonevsnone","nonevsmatch","nonevssilver","matchvsnone","matchvsmatch","matchvssilver","silvervsnone","silvervsmatch","silvervssilver")
pairs<-factor(pairs,levels=c("nonevsnone","nonevsmatch","nonevssilver","matchvsnone","matchvsmatch","matchvssilver","silvervsnone","silvervsmatch","silvervssilver"))



####################
#Part 1 - simulations
####################

#if a = 2, neutrality, a = 1 is no competition

for(a in 1:5){

win<-matrix(NA,nrow=500,ncol=11)

for (i in 1:length(pairs)) {
  
me.typei<-sapply(strsplit(as.character(pairs[i]),"vs"), `[`, 1)
me.typej<-sapply(strsplit(as.character(pairs[i]),"vs"), `[`, 2)

#if loop is replicating each level 100 times, k loop is the amount of autocor, k=1 is always repeating, k=11 #is never repeating
  
for(k in 1:11){
      
    for(f in 1:500){
    
    S<-data.frame(sp1=rep(NA,nyear),sp2=rep(NA,nyear),year=rep(NA,nyear),year.comb=rep(NA,nyear))  #matrix of N; each row is a year
    
    S[1:2,1]<-2
    S[1:2,2]<-2
    
    S[1,3]<-rbinom(1,1,prob=pfav) #random 1st year condition
    
    #as k increases, the probability of switching also increases
    pswitch<-((k-1)/10)
    
    bin<-rbinom(nyear,1,prob=pswitch) #calcs when things switch, but when pswitch = 1, no switching - always good or bad
    
    for(j in 2:nyear){
    S[j,3]<-abs(S[j-1,3]-bin[j])
    S[j,4]<-paste(S[j-1,3],S[j,3],sep="")
    }  
    
    plot(1:100,S$sp1,col="white",ylim=c(0,100),xlim=c(0,100),xlab="Year",ylab="Population size (N)")
    
      for(t in 3:nyear){
      
      tmp1<-S[t-1,3]
      tmp2<-S[t,3]
      tmp3<-S$year.comb[t]
      
      aii<-ifelse(tmp2==0,aii.b[a],aii.g[a]); ajj<-ifelse(tmp2==0,ajj.b[a],ajj.g[a])
      aij<-ifelse(tmp2==0,aij.b[a],aij.g[a]); aji<-ifelse(tmp2==0,aji.b[a],aji.g[a])  
      
      lambdai<-ifelse(tmp2==0,lambdai.b,lambdai.g);lambdaj<-ifelse(tmp2==0,lambdaj.b,lambdaj.g)
      
      tradeoff.typei<-"seed.pos"
      tradeoff.typej<-"seed.pos"
      
      #here is where a tradeoff could be included if wanted - if the code below were not commented out, there would be a tradeoff
      #effect of current environment on # of seeds produced
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
    
    lines(1:100,S$sp1,lwd=2,col="#545454")
    lines(1:100,S$sp2+1.5,lwd=2, col="#6BAA75")
    text(1:100,rep(100,100),labels="+",col=ifelse(S$year=="1","black","white"))
    text(1:100,rep(95,100),labels="-",col=ifelse(S$year=="0","black","white"))
    
    #for replicate 50 (random replicate), print the plot so we get one plot per autocorrelation x ME type
    if(f==50) {dev.copy(pdf,height=3.5,width=10,paste("./dynamics plots/",scenario[a],"/",pairs[i]," k=",k," 500 sm neutral.pdf",sep=""));dev.off()}

    
    if(S[t,1]!=0 & S[t,2]!=0) win[f,k]<-3
    if(S[t,1]!=0 & S[t,2]==0) win[f,k]<-1
    if(S[t,1]==0 & S[t,2]!=0) win[f,k]<-2
    
    print(t)
    } 
  }

#1 = 1 wins, 2 = 2 wins, 3 = coexist
prop.1<- sapply(1:ncol(win),function(k) sum(win[,k]==1))/500
prop.2<- sapply(1:ncol(win),function(k) sum(win[,k]==2))/500
prop.3<- sapply(1:ncol(win),function(k) sum(win[,k]==3))/500

plot(x=1:11,y=prop.1+0.01,ylim=c(0,1),type="l",col="#545454",xlab="Autocorrelation",ylab="Proportion outcomes",lwd=2,main=paste("sp1=",me.typei,", sp2=",me.typej,sep=""))
points(x=1:11,y=prop.2,ylim=c(0,1),type="l",col="#6BAA75",lwd=2)
points(x=1:11,y=prop.3,ylim=c(0,1),type="l",col="#84DD63",lwd=2)

dev.copy(pdf,paste("./summary plots/",scenario[a],"/",pairs[i]," 500.pdf",sep=""));dev.off()

assign(paste("win",pairs[i],sep=""),win)

write.csv(win,paste("./summary plots/",scenario[a],"/",pairs[i]," 500.csv",sep=""))


}
}


#################
#plots for ms
#################

#figure 3

  for(i in c(1,7,4,3,9,6,2,8,5)){
    
    ms<-read.csv(paste("./summary plots/neutral/",pairs[i]," 500 sm neutral.csv",sep=""))[,2:12]


    for (k in 1:11) {
    win=ms

    win<-win[,seq(from = 1, to = 11, by=1)]

    #1 = 1 wins, 2 = 2 wins, 3 = coexist
    prop.1<- sapply(1:ncol(win),function(k) sum(win[,k]==1))/500
    prop.2<- sapply(1:ncol(win),function(k) sum(win[,k]==2))/500
    prop.3<- sapply(1:ncol(win),function(k) sum(win[,k]==3))/500

    prop.1<-as.data.frame(cbind(seq(from=1,to=0,by=-0.1),prop.1,"spi"));colnames(prop.1)<-c("k","prop","sp")
    prop.2<-as.data.frame(cbind(seq(from=1,to=0,by=-0.1),prop.2+0.02,"spj"));colnames(prop.2)<-c("k","prop","sp")
    prop.3<-as.data.frame(cbind(seq(from=1,to=0,by=-0.1),prop.3,"both"));colnames(prop.3)<-c("k","prop","sp")

    prop<-rbind(prop.1,prop.2,prop.3)
    prop$prop<-as.numeric(as.character(prop$prop))
    prop$sp<-factor(prop$sp,levels=c("spi","spj","both"))
    
      }

    plot<-ggplot(data=prop,aes(x=k,y=prop, fill=sp,group=sp, col=sp)) + theme_classic() + theme(legend.text.align = 0, legend.key.size = unit(0.5, "cm"),legend.position="none",axis.text=element_text(size=7),panel.border = element_rect(colour = "black", fill=NA)) +
      geom_line(size=1.5) + 
      scale_color_manual(values=c("#545454","#6BAA75","#84DD63"),name="Outcome", labels=c(expression(italic("i")~"wins"), expression(italic("j")~"wins"), "neutrality")) +
      labs(x=ifelse(i==1|i==9|i==5,expression("Autocorrelation ("~italic("k")~")"),""),y=ifelse(i==1|i==9|i==5,"Proportion outcomes","")) +
      scale_y_continuous(breaks = seq(0, 1, by = 0.2)) + 
      ylim(0,1.02) 
    
    assign(paste("plot.",pairs[i],sep=""),value=plot)

  }

  grid<-plot_grid(plot.nonevsnone,plot.nonevssilver,plot.nonevsmatch,
          NULL,plot.silvervssilver+theme(legend.position=c(-1, 0.5)),plot.silvervsmatch,
          NULL,NULL,plot.matchvsmatch, align="hv",
          nrow=3, ncol=3,labels=c("A","B","C","","D","E","","","F"), hjust=c(-5,-5,-5,-5,-5,-5,-5,-5,-5.5), vjust=2.5)
  grid
  ggsave("figure 3.pdf",w=8, h=8)

