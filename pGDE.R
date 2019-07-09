############
###Descriptoion:
###INPUT: expr-a matrix list of transcriptomics data; lab-a list of sample label ("treatment" or "control") vectors 
###nper-the number of label permutations in permutation test. When set to NULL, alternatively a asymptotical significance analysis will be used.  
pGDE<-function(expr,lab,nper=NULL){
  studys=NULL
  for(i in 1:length(expr)){
    studys[[i]]=cbind(expr[[i]][,c(which(lab[[i]]=="treatment"))],expr[[i]][,c(which(lab[[i]]=="control"))])
  }
  samples=array(0,c(length(expr),2))
  for(i in 1:length(expr)){
    samples[i,1]=length(which(lab[[i]]=="treatment"))
    samples[i,2]=length(which(lab[[i]]=="control"))
  }
  expr=studys
  res=regulation(expr,samples)#observed pd
  pud=res$pud
  pval=res$pval
  pm=NULL;
  pn=list();
  if(!is.null(nper)){
  for(j in 1:nper){
    a=regulation(expr,samples,1)#permutated pd
    pm=c(pm,a$pud)
    pn[[j]]=a$pud
    row.names(pn[[j]])=rownames(a);
    print(paste("permutation test:",j))
  }
  pm=abs(pm);
  pval=sapply(abs(pud),function(x){length(pm[pm>=x])/length(pm)})
  }
  pval_BH=p.adjust(pval,method = "BH",n=length(pval))
  return(list(pud=pud,pval=pval,pval_BH=pval_BH,pn=pn))
}


regulation<-function(expr,samples,miss=0){##miss=1: permuatation test
  n=nrow(expr[[1]])  ## number of genes (features)
  s=nrow(samples)  ## number of studies

up1=NULL;
down1=NULL;  
up2=NULL;
down2=NULL;  

  if(miss==0){
    for(i in 1:s){#ÑÐ¾¿
        s1=samples[i,1]
        s2=samples[i,1]+samples[i,2]
       upt=compUp(expr[[i]][,1:s1],expr[[i]][,(s1+1):s2])
       downt=compDown(expr[[i]][,1:s1],expr[[i]][,(s1+1):s2])
        up1=cbind(up1,upt$up1); 
		up2=cbind(up2,upt$up2); 
        down1=cbind(down1,downt$down1);
        down2=cbind(down2,downt$down2);
    }

  }else{
    for(i in 1:s){#studies
        s1=samples[i,1]
        s2=samples[i,1]+samples[i,2]
        label=c(1:s2)
        label=sample(label)###label permuatation
            upt=compUp(expr[[i]][,label[c(1:s1)]],expr[[i]][,label[c((s1+1):s2)]])
       downt=compDown(expr[[i]][,label[c(1:s1)]],expr[[i]][,label[c((s1+1):s2)]])
        up1=cbind(up1,upt$up1); 
		up2=cbind(up2,upt$up2); 
        down1=cbind(down1,downt$down1);
        down2=cbind(down2,downt$down2);
    }
  }
	pu=(ncol(up1)*rowMeans(up1)+ncol(up2)*rowMeans(up2))/(ncol(up1)+ncol(up2))
	pd=(ncol(down1)*rowMeans(down1)+ncol(down2)*rowMeans(down2))/(ncol(down1)+ncol(down2))
	pud=pu-pd;
	n=ncol(up1);
	m=ncol(up2);
	p=0.5;
	Z=(m+n)*(pud+1)/2-(m+n)*p;
	p.value <- 2 * (1 - pnorm(abs(Z), 0,((m+n)*p*(1-p))^0.5))
    return(list(pu=pu,pd=pd,pud=pud,pval=p.value))
}




compUp<-function(d1,d2){#d1-matrix of treatment data£¬d2-matrix of control data£»
if(class(d1)=="numeric") d1=matrix(d1,nrow=1);
if(class(d2)=="numeric") d2=matrix(d2,nrow=1);
up1=apply(d1,2,function(x){rowSums(x>=d2)/ncol(d2)});
up2=apply(d2,2,function(x){rowSums(x<=d1)/ncol(d1)});

  return(list(up1=up1,up2=up2))
}
compDown<-function(d1,d2){#d1-matrix of treatment data£¬d2-matrix of control data£»
down1=apply(d1,2,function(x){rowSums(x<=d2)/ncol(d2)});
down2=apply(d2,2,function(x){rowSums(x>d1)/ncol(d1)});

  return(list(down1=down1,down2=down2))
}
