Jfunc = function(alpha,q,Jmax) {
  (alpha*q/(sqrt(1+(alpha*alpha*q*q)/(Jmax*Jmax))))
}

q = 1:2000

plot(q,Jfunc(.55,q,200),ylim=c(0,250))
abline(h=200)
lines(q,Jfunc(.85,q,200),ylim=c(0,250),col="pink")
lines(q,Jfunc(.25,q,200),ylim=c(0,250),col="blue")


cctr = function(Kc,Ko,J,V,G,O) {
  ((Kc*J*(Ko+O)) - (8*Ko*G*V)) / ((4*V-J)*Ko)
}

citr = function(Kc,Ko,J,V,G,O) {
  k2 = Kc*(1+(O/Ko))
  ((V/J)*8*G-k2) / (1-(4*(V/J)))
}

O is 21
J

cctr(Kc=27,Ko=16,J=300,V=200,G=3.74,O=21)
citr(Kc=27,Ko=16,J=300,V=200,G=3.74,O=21)



##examine Jmax/Vcmax ratio in light of Aci curves

#kinetic constants at 25 C
R=0.008314 #(kJ mol^-1 K^-1)
Kc=exp(35.9774-80.99/(R*(25+273.15))) #Michaelis-Menten constant for Rubisco for O2 (Pa)
Ko=exp(12.3772-23.72/(R*(25+273.15))) #Michaelis-Menten constant for Rubisco for CO2 (kPa)
Gs=exp(11.187-24.46/(R*(25+273.15))) #Photorespiration compensation point (Pa)
O=21 #oxygen (O2) partial pressure (kPa)  
R; Kc; Ko; Gs; O

ACI = function(Vcmax,Jmax,Ci,Gs=Gs,Kc=Kc,Ko=Ko,O=O) {
  k2 = Kc*(1+(O/Ko))
  Ac = (Vcmax*(Ci-Gs)) / (Ci+k2)
  Aj = (Jmax*(Ci-Gs)) / (4*Ci + 8*Gs)
  out = rep(0,length(Ac))
  for(i in 1:length(Ac)) {
     out[i] = min(Ac[i],Aj[i]) 
  }
  return(out)
}

Ci = 0:60
curve = ACI(Vcmax=50,Jmax=200,Ci=Ci,Gs=Gs,Kc=Kc,Ko=Ko,O=O)
plot(Ci,curve,type="l",lwd=2)

#Vcmax is typically between 50 and 200 for our species
#Jmax/Vcmax ration is between 1 and 2.5

curve = ACI(Vcmax=50,Jmax=50,Ci=Ci,Gs=Gs,Kc=Kc,Ko=Ko,O=O)
plot(Ci,curve,type="l",lwd=2,col="yellow",ylim=c(-5,60))
curve = ACI(Vcmax=50,Jmax=1.7*50,Ci=Ci,Gs=Gs,Kc=Kc,Ko=Ko,O=O)
lines(Ci,curve,type="l",lwd=2,col="orange")
curve = ACI(Vcmax=50,Jmax=2.5*50,Ci=Ci,Gs=Gs,Kc=Kc,Ko=Ko,O=O)
lines(Ci,curve,type="l",lwd=2,col="red")

curve = ACI(Vcmax=150,Jmax=150,Ci=Ci,Gs=Gs,Kc=Kc,Ko=Ko,O=O)
lines(Ci,curve,type="l",lwd=2,col="cyan")
curve = ACI(Vcmax=150,Jmax=1.7*150,Ci=Ci,Gs=Gs,Kc=Kc,Ko=Ko,O=O)
lines(Ci,curve,type="l",lwd=2,col="blue")
curve = ACI(Vcmax=150,Jmax=2.5*150,Ci=Ci,Gs=Gs,Kc=Kc,Ko=Ko,O=O)
lines(Ci,curve,type="l",lwd=2,col="purple")


##Simulate limitations on photosynthesis
#Ac or Aj, each limited by either substrate or enzymes: 4 possible limitations
  #ie for Ac, limited by Rubisco (Ac limitation if CO2=infinity) or CO2?
    #or rather than Ci = infinity, Ci=40?
  #for Aj, limited by chl/thylakoid proteins (Aj limitation if q=2000) or q?

#***add alpha limitation too? ie alpha -> 1? or something like 0.8?

#kinetic constants at 25 C
R=0.008314 #(kJ mol^-1 K^-1)
Kc=exp(35.9774-80.99/(R*(25+273.15))) #Michaelis-Menten constant for Rubisco for O2 (Pa)
Ko=exp(12.3772-23.72/(R*(25+273.15))) #Michaelis-Menten constant for Rubisco for CO2 (kPa)
Gs=exp(11.187-24.46/(R*(25+273.15))) #Photorespiration compensation point (Pa)
O=21 #oxygen (O2) partial pressure (kPa)  
R; Kc; Ko; Gs; O

Ac = function(Vcmax,Ci,Gs,Kc,Ko,O) {
  k2 = Kc*(1+(O/Ko))
  (Vcmax*(Ci-Gs)) / (Ci+k2)
}  
  
Aj = function(alpha,Jmax,q,Gs,Ci) {
  J = (alpha*q/(sqrt(1+(alpha*alpha*q*q)/(Jmax*Jmax))))
  return( (J*(Ci-Gs)) / (4*Ci + 8*Gs) )
}  

N = 100
Ci.vec = seq(0,40,length=N)
q.vec = seq(0,2000,length=N)
out = matrix(0,nrow=N*N,ncol=7)
Vcmax = 150
Jmax = 250
alpha = .5
nx = 1
for(i in 1:N) {
  for(j in 1:N) {
    Ci = Ci.vec[i]; q = q.vec[j]
    Ac1 = Ac(Vcmax=300,Ci=Ci,Gs,Kc,Ko,O)            #Ac at ambient CO2 but high Vcmax
    Ac2 = Ac(Vcmax=Vcmax,Ci=40,Gs,Kc,Ko,O)          #Ac at Vcmax but max CO2 
    Aj1 = Aj(alpha=alpha,Jmax=Jmax,q=q,Gs,Ci=Ci)    #Aj at ambient conditions
    Aj2 = Aj(alpha=alpha,Jmax=Jmax,q=2000,Gs,Ci=Ci) #Aj at ambient CO2 but max light (alpha irrelevant)
    Aj3 = Aj(alpha=.8,Jmax=Jmax,q=q,Gs,Ci=Ci)       #Aj at ambient CO2 but max alpha (=.8)
    out[nx,] = c(Ci,q,Ac1,Ac2,Aj1,Aj2,Aj3)
    nx = nx + 1
  }    
}    
minout = unlist(apply(out[,3:7],1,function(x)which(x==min(x))))
hist(minout)

#plot
cols = c("red2","red4","palegreen1","springgreen4","darkgreen")

plot(out[,2],out[,1],type="n",)
points(out[,2],out[,1],col=cols[minout],pch=19)
    
  #interpretation on which of the above are the minimum:
  #Ac1: (red): Anet limited by CO2 (would increase at higher CO2)
  #Ac2: (darkred): Anet limited by Rubisco content (increase in CO2 would not change Anet)
  #Aj1 (palegreen): Anet limited by chl content (if alpha higher, Anet increases at same q)
  #Aj2 (medium green): Anet limited by Jmax (photons cannot get higher, alpha irrelevant at q=2000)
  #Aj3 (darkgreen): Anet limited by photons (alpha at max)











