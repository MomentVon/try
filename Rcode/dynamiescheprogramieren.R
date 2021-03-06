#Hochwasser Regulieren durch dynamischen Programieren
##########################Parameter dynamiesche Progrmierung############
constMin=0.0
constMax=999.9
constBrech=-1.0
variblen=7  #Nutzeffekt, Entscheidungvariable, Zustandvariable Anfang, Zustandvariable End
##########################Parameter Talsprren##############
totraum=11.0  #10^9m3
vorhwbegrenzraum=26.0  #10^9m3
nutzraum=29.8  #10^9m3
hwschutzraum=30.0  #10^8m3
ganzraum=32.0 #10^8m3
hwschutzraumabfmax=11.0  #10^m3/s
ganzraumabfmax=13.0  #10^m3/s
fristhwanf=6  #MOnat, Juli
fristhwend=9  #September
##########################Parameter Hochwasser regieren#################
hwdltat=3*3600  #s
hwdltatvarible=3*3.6/100  #Masseinheit wechseln
hwabfmax=22.0  #10^3m3/s
hwabfmin=0.26  #10^3m3/s
hwschrittn=24
hwzustandvariblekl=vorhwbegrenzraum
hwzustandvariblegr=35.0
hwzustandvaribleprazision=0.1
hwzustanddiskretn=floor(abs((hwzustandvariblekl - hwzustandvariblegr) / hwzustandvaribleprazision))
HWZUFLUSS=array(scan("input/ihw.prn"),dim = c(1,hwschrittn))/1000.0  #10^3m3/s
HWWHDIS=array(seq(hwzustandvariblekl,constMax,by=hwzustandvaribleprazision),dim=c(hwzustanddiskretn))
HWZUSTANDANF=array(26.0,dim=c(hwzustanddiskretn)) #10^8m3
HWNUTZEFFEKANF=array(0,dim=c(hwzustanddiskretn))

##########################Parameter Wasserkraft##################
wkdltats=10*24*3600  #s
wkdltath=10*24  #h
wkwechsel1=3.6  #cl=q/wr(hsl) 
wkwechsel2=wkdltats/10^8 #Masseinheit wechseln
wkabfmax=11000  #m3/s
wkabfmin=200  #m3/s
wkmincl=176  #MW 保证出力
wkschrittn=36
wkzustandvariblekl=totraum
wkzustandvariblegr=nutzraum
wkzustandvaribleprazision=0.3
wkzustanddiskretn=ceiling(abs((wkzustandvariblekl - wkzustandvariblegr) / wkzustandvaribleprazision))
wkfristhwanf=3*(fristhwanf-1)+1
wkfristhwend=3*fristhwend
WKCLXZ=read.table("input/xzcl.prn",TRUE)  #st:m,ll:m3/s,cl:10^6W(MW),hsl:m3/kWh
WKWSW=read.table("input/wsw.prn",TRUE)  #ll:m3/s,wsw:m
WKKRQX=read.table("input/krqx.prn",TRUE)  #kr:m3,sw:m
WKZUFLUSS=array(scan("input/iwk.prn"),dim = c(1,wkschrittn))  #m3/s
WKWHDIS=array(seq(wkzustandvariblekl,constMax,by=wkzustandvaribleprazision),dim=c(wkzustanddiskretn))
WKZUSTANDANFANG=WKWHDIS  #10^8m3
WKNUTZEFFEKANF=array(0,dim=c(wkzustanddiskretn))
##########################Funktion andere kleine Tool#############
interpolation <- function(bktprm,BKTPRMLST,ZLPRMLST){
  
  n=length(BKTPRMLST)
  if(bktprm < BKTPRMLST[1]){
    zp=(bktprm - BKTPRMLST[1])*(ZLPRMLST[1] - ZLPRMLST[2])/(BKTPRMLST[1] - BKTPRMLST[2]+0.0000000001) + ZLPRMLST[1]
    
  }
  if(bktprm > BKTPRMLST[n]){
    zp=(bktprm - BKTPRMLST[n-1])*(ZLPRMLST[n-1] - ZLPRMLST[n])/(BKTPRMLST[n-1] - BKTPRMLST[n]+0.0000000001) + ZLPRMLST[n-1]
    
  }
  if(bktprm >= BKTPRMLST[1]) {
    for(i in 1:n){
      if(bktprm < BKTPRMLST[i]){
        zp=(bktprm - BKTPRMLST[i-1])*(ZLPRMLST[i-1] - ZLPRMLST[i])/(BKTPRMLST[i-1] - BKTPRMLST[i]+0.000000001) + ZLPRMLST[i-1]
        break()
      }
    }
  }
  return(zp)
}

##########################Funktion Hochwasser#####################################################
hwwasserhaltausgleichzustandwendfunk <- function(whanf,whend,ii){
  return((whanf - whend + HWZUFLUSS[ii]*hwdltatvarible)/hwdltatvarible)
}

tausendjahrhwregierennutzeffkfunk <- function(VQVV){
  VQVV[1]=max(VQVV[3],VQVV[4])
  return(VQVV)
}

hwregierenurteilenfunk <- function(VQVV,ii=0){
  vmax=VQVV[1]
  v=VQVV[3]
  q=VQVV[2]
  if(q<hwabfmin){
    return(constBrech)
  }
  if(v==constBrech){
    return(constBrech)
  }
  if(v<hwschutzraum & q>hwschutzraumabfmax){
    return(constBrech)
  }
  if((v>=hwschutzraum & v<ganzraum) & q>ganzraumabfmax){
    return(constBrech)
  }
  if(v>=ganzraum  & q>hwabfmax){
    return(constBrech)
  }
  else return(VQVV)
}

tausendjahrhwregierenzielfunk <- function(VQVVDIS,VORCONTROL){
  for(i in 1:hwzustanddiskretn){
    VQVVDIS[i,1]=ifelse(VQVVDIS[i,3]==constBrech,constMax,max(VQVVDIS[i,4],VORCONTROL[i]))
  }
  return(VQVVDIS[(which.min(VQVVDIS[,1])),])
}
##########################Funktion Wasserkraft############
wkwasserhaltausgleichzustandwendfunk <- function(whanf,whend,ii){
  return((whanf - whend + WKZUFLUSS[ii]*wkwechsel2)/wkwechsel2) ####
}

wknutzeffkfunk <- function(EQVV){
  q=EQVV[2]
  vanf=EQVV[3]
  hwsw=interpolation(q,WKWSW$ll,WKWSW$wsw)
  hanf=interpolation(vanf,WKKRQX$kr,WKKRQX$sw)
  hst=hanf-hwsw
  EQVV[5]=hst
  hsl=interpolation(hst,WKCLXZ$st,WKCLXZ$hsl)
  EQVV[6]=hsl
  cl=q/hsl*wkwechsel1  #MW
  clmax=interpolation(hst,WKCLXZ$st,WKCLXZ$cl)
  cl=min(cl,clmax)
  EQVV[7]=cl
  EQVV[1]=max((cl*wkdltath),constMin) #MWh
  return(EQVV)
}

wkregierenurteilenfunk <- function(EQVV,ii){
  e=EQVV[2]
  q=EQVV[2]
  vanf=EQVV[3]
  hwsw=interpolation(q,WKWSW$ll,WKWSW$wsw)
  hanf=interpolation(vanf,WKKRQX$kr,WKKRQX$sw)
  hst=hanf-hwsw
  hsl=interpolation(hst,WKCLXZ$st,WKCLXZ$hsl)
  cl=q/hsl*wkwechsel1  #MW
  clmax=interpolation(hst,WKCLXZ$st,WKCLXZ$cl)
  cl=min(cl,clmax)
  
  
  if(vanf==constBrech){
    return(constBrech)
  }
  if(ii >= wkfristhwanf & ii <= wkfristhwend & vanf > vorhwbegrenzraum){
    return(constBrech)
  }
  if(cl<wkmincl){
    return(constBrech)
  }
  if((vanf>=hwschutzraum & vanf<ganzraum) & q>ganzraumabfmax){
    return(constBrech)
  }
  if(vanf>=ganzraum  & q>hwabfmax){
    return(constBrech)
  }
  else return(EQVV)
}

wkregierenzielfunk <- function(EQVVDIS,VORCONTROL){
  for(i in 1:wkzustanddiskretn){
    EQVVDIS[i,1]=ifelse((EQVVDIS[i,3]==constBrech),constMin,(ifelse((EQVVDIS[,1]==VORCONTROL[]),EQVVDIS[i,1],(EQVVDIS[i,1]+VORCONTROL[i]))))
  }
  mark=max(EQVVDIS[,1])
  if(sum(EQVVDIS[,1]==mark)==1){
    markfolge=which.max(EQVVDIS[,1])
  }
  if(sum(EQVVDIS[,1]==mark)!=1){
    markfolge=which.max(EQVVDIS[,4]*(EQVVDIS[,1]==mark))
  }
  return(EQVVDIS[markfolge,])
}

##########################HauptFunktion##########################
dynamisheprogramierenfunk <- function(schrittn,zustanddiskretn,ZUSTANDANF,NUTZEFFEKANF,ZUSTANDDISKRET,zielfk,nutzeffekfk,zustandwendfk,urteilenfk){
  CONTROL<- array(constBrech, dim=c(zustanddiskretn,variblen,schrittn))
  for(i in 1:schrittn){
    CONTROL[,4,]=ZUSTANDDISKRET
  }
  j=1
  for(i in 1:zustanddiskretn){
    TEM=array(constBrech,dim=c(zustanddiskretn,variblen)) 
    TEM[,3]=ZUSTANDANF #CONTROL[,4,j-1]
    TEM[,4]=CONTROL[i,4,j] #
    for(k in 1:zustanddiskretn){ #
      TEM[k,2]=zustandwendfk(TEM[k,3],TEM[k,4],j) #
      TEM[k,]=nutzeffekfk(TEM[k,])
      TEM[k,]=urteilenfk(TEM[k,],j)
    }
    CONTROL[i,,j]=zielfk(TEM,NUTZEFFEKANF)#####
    CONTROL[i,,j]=urteilenfk(CONTROL[i,,j],j)
  }
  
  for(j in 2:schrittn ){
    for(i in 1:zustanddiskretn){
      TEM=array(constBrech,dim=c(zustanddiskretn,variblen)) 
      TEM[,3]=CONTROL[,4,j-1] #
      TEM[,4]=CONTROL[i,4,j] #
      for(k in 1:zustanddiskretn){ #
        TEM[k,2]=zustandwendfk(TEM[k,3],TEM[k,4],j) #
        TEM[k,]=nutzeffekfk(TEM[k,])
        TEM[k,]=urteilenfk(TEM[k,],j)
      }
      CONTROL[i,,j]=zielfk(TEM,CONTROL[,1,j-1])#####
      CONTROL[i,,j]=urteilenfk(CONTROL[i,,j],j)
    }
  }
  return(CONTROL)
}

wegsuchen <- function(CONTROL,schrittn,zielfk){
  BESTEND=zielfk(CONTROL[,,schrittn],CONTROL[,1,schrittn])
  bestendzustand=BESTEND[4]
  bestendnutzeffekt=BESTEND[1]
  temmaker=bestendzustand
  BESTWEG <- array(constBrech,dim=c(schrittn+1,variblen))
  BESTWEG[schrittn+1,1]=bestendzustand
  BESTWEG[schrittn+1,2]=bestendnutzeffekt
  for(i in schrittn:1){
    temmakerfolge=which(CONTROL[,4,i]==temmaker)
    temmaker=CONTROL[temmakerfolge,3,i]
    BESTWEG[i,]=CONTROL[temmakerfolge,,i]
  }
  return(BESTWEG)
}

keismainfunk <- function(umkreisn,dynamisheprogramierenfunk,schrittn,zustanddiskretn,ZUSTANDANF,NUTZEFFEKANF,ZUSTANDDISKRET,zielfk,nutzeffekfk,zustandwendfk,urteilenfk){
  CONTROLCONTROL=array(constBrech,dim=c(zustanddiskretn,variblen,schrittn,umkreisn))
  BESTWEG <- array(constBrech,dim=c(schrittn+1,2,umkreisn))
  
  for(i in 1:umkreisn){
    # minus=min(NUTZEFFEKANF[which(NUTZEFFEKANF!=constBrech)])
    # NUTZEFFEKANF=NUTZEFFEKANF-minus
    CONTROLCONTROL[,,,i]=dynamisheprogramierenfunk(schrittn,zustanddiskretn,ZUSTANDANF,NUTZEFFEKANF,ZUSTANDDISKRET,zielfk,nutzeffekfk,zustandwendfk,urteilenfk)
    NUTZEFFEKANF=CONTROLCONTROL[,1,schrittn,i]
    ZUSTANDANF=CONTROLCONTROL[,,schrittn,i]
    if(length(which(ZUSTANDANF!=constBrech))==0){
      print("tut mir leid")
      break()
    }
    BESTWEG[,,i]=wegsuchen(CONTROLCONTROL[,,,i],schrittn,zielfk)
    if(i>1 & (sum((BESTWEG[,1,i]-BESTWEG[,1,i-1]))==0) ) {
      
      return(BESTWEG)
    }   
  }
  return(BESTWEG)
}

##########################Funktion machen################
# HWCONTROL=dynamisheprogramierenfunk(hwschrittn,hwzustanddiskretn,HWZUSTANDANF,HWNUTZEFFEKANF,HWWHDIS,tausendjahrhwregierenzielfunk,tausendjahrhwregierennutzeffkfunk,hwwasserhaltausgleichzustandwendfunk,hwregierenurteilenfunk)
# HWWEG=wegsuchen(HWCONTROL,hwschrittn,tausendjahrhwregierenzielfunk)
# AA1=dynamisheprogramierenfunk(wkschrittn,wkzustanddiskretn,WKZUSTANDANFANG,WKNUTZEFFEKANF,WKWHDIS,wkregierenzielfunk,wknutzeffkfunk,wkwasserhaltausgleichzustandwendfunk,wkregierenurteilenfunk)
# AA2=dynamisheprogramierenfunk(wkschrittn,wkzustanddiskretn,AA1[,4,36],AA1[,1,36],WKWHDIS,wkregierenzielfunk,wknutzeffkfunk,wkwasserhaltausgleichzustandwendfunk,wkregierenurteilenfunk)
# WKWEG2=wegsuchen(AA2,wkschrittn,wkregierenzielfunk)
