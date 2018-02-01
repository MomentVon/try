##############################p3 Kurve######
hwscheiteln=54
hwstichproben=24
tausend8=16.23*10^8  #m3
tausend24=35.8*10^8  #m3
dltat=3*3600  #s
constMaxp3=33333  #fur p3 unendlich
constMaxnp3=99999  #die Numer von integrate Funktion
constMax=9999
constFlnq=0.618
HWSCHEITEL=array(0,dim = c(hwscheiteln,5))  #Hochwasserscheitel Stichprobe
HWSCHEITEL[,1] <- sort(scan("input/hws.prn"),TRUE)
HWSCHEITEL[,2] <- c(1:hwscheiteln)
HWSCHEITEL[,3] <- HWSCHEITEL[,2]/(hwscheiteln+1)
HWSTICHPROBE1997=array(scan("input/hwstich.prn"),c(hwstichproben))
p3durchschnitthw=mean(HWSCHEITEL[,1])
HWSCHEITEL[,4]=(HWSCHEITEL[,1]/p3durchschnitthw-1)
p3abweichungskoffizient=(sum((HWSCHEITEL[,4])^2)/(hwscheiteln-1))^0.5
bestp3skewkoffizient=p3skewkoffizient0=(sum((HWSCHEITEL[,4])^3))/((hwscheiteln-3)*p3abweichungskoffizient^3 )

##############################Funktion##########
probemaxx <- function(A,x){
  zusammern=length(A)
  BEWEGWN=array(A[1:x])
  mark=1
  addition=sum(BEWEGWN)
  for(i in 2:(zusammern-x)){
    BEWEGWN=A[i:(i+x-1)]
    if((sum(BEWEGWN)) > addition ){
      mark=i
    }
  }
  return(A[mark:(mark+x-1)])
}

absolutabweichungfunk <- function(A,B){  ##A,B sind moglichkeit
  return(sum(abs(A-B)))
}

beweichen618funk <- function(paramekl,paramegr,DATENORI,abweichfk,processfk,umkriesn){
  paramelk=paramekl
  paramerc=paramegr
  abweischunglk=abweichfk(DATENORI,processfk(paramelk)) 
  abweischungrc=abweichfk(DATENORI,processfk(paramerc)) 
  for(i in 1:umkriesn){
    parametemrc=paramelk+0.618*(paramerc - paramelk)
    parametemlk=paramerc+0.618*(paramelk - paramerc)
    abweischungtemrc=abweichfk(DATENORI,processfk(parametemrc)) 
    abweischungtemlk=abweichfk(DATENORI,processfk(parametemlk)) 
    if(abweischungtemrc<abweischungtemlk){
      paramelk=parametemlk
    }
    if(abweischungtemlk<abweischungtemrc){
      paramerc=parametemrc
    }
  }
  return(mean(paramelk,paramerc))
}

p3moglichkeitdichtfunk <- function(x){  
  p3a0=p3durchschnitthw*(1-2*p3abweichungskoffizient/bestp3skewkoffizient)
  p3af=4/bestp3skewkoffizient^2
  p3bt=2/(p3durchschnitthw*p3abweichungskoffizient*bestp3skewkoffizient)
  return(p3bt^p3af/(gamma(p3af))*(x-p3a0)^(p3af-1)*exp(-p3bt*(x-p3a0))) 
}

p3moglichkeitfunk <- function(x){
  return(integrate(p3moglichkeitdichtfunk,x,constMaxp3)$value)
}
p3funk <- function(p3skewkoffizient,A=HWSCHEITEL[,1]){
  p3moglichkeitdichtfunk <- function(x){  
    p3a0=p3durchschnitthw*(1-2*p3abweichungskoffizient/p3skewkoffizient)
    p3af=4/p3skewkoffizient^2
    p3bt=2/(p3durchschnitthw*p3abweichungskoffizient*p3skewkoffizient)
    return(p3bt^p3af/(gamma(p3af))*(x-p3a0)^(p3af-1)*exp(-p3bt*(x-p3a0))) 
  }
  p3moglichkeitfunk <- function(x){
    return(integrate(p3moglichkeitdichtfunk,max(x,2200),constMaxp3)$value)
  }
  p3moglichkeitlistfunk <- function(A){  ##A ist Hwscheitel, aber B ist die Moglichkeit von diese Hwscheitel
    kreisn=length(A)
    B=A
    for(i in 1:kreisn){
      B[i]=p3moglichkeitfunk(A[i])
    }
    return(B)
  }
  
  return(p3moglichkeitlistfunk(A))
  
}


##############################Funktion machen###########
bestp3skewkoffizient=beweichen618funk(p3skewkoffizient0*0.6,p3skewkoffizient0*1.1,HWSCHEITEL[,3],absolutabweichungfunk,p3funk,10)

entwerfenhw=beweichen618funk(HWSCHEITEL[1,1],45000, 0.001,absolutabweichungfunk,p3moglichkeitfunk,10)

###
hw8=sum(probemaxx(HWSTICHPROBE1997,8))  #m3/s
tausendselbst8=hw8*dltat
hw24=sum(HWSTICHPROBE1997)
tausendselbst24=hw24*dltat
tausendselbst8/10^8
tausendselbst24/10^8

p3moglichkeitfunk(entwerfenhw)
