pdf(file='transcript_pooled_screen_1_counts_mageck.pdf',width=4.5,height=4.5);
gstable=read.table('transcript_pooled_screen_1_counts_mageck.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("ENST00000493536","ENST00000451311","ENST00000293831","ENST00000245458","ENST00000597224","ENST00000330752","ENST00000392367","ENST00000359611","ENST00000225792","ENST00000579618")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='sym-01-01-TP6_vs_Sym-pool-1-1,Sym-pool-1-2,Sym-pool-1-3 neg.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1579.0103584784952,1704.511318020994,1623.3613216124863,148.21882484565762),c(1872.403651129218,1918.857891792196,1857.9970834330072,132.06677342016928),c(2518.0999133014775,2608.6434084495213,2635.8742870278193,274.5848742333016),c(3251.5831449282837,3564.081859621155,3277.742286313109,303.0884943959281),c(2702.914585837366,2908.5005834592343,2809.2661381358657,100.71279124128017),c(2315.9588652153498,2591.541288201819,2442.598049324814,567.2220412362666),c(2937.398201617274,2951.8259547534135,2771.8834913373416,1231.3563910254632),c(2762.979354411529,2767.123056078229,2702.685825986883,599.5261440872433),c(2698.2942190239683,2594.961712251359,2643.0326662020047,177.67256568037163),c(2263.979738564631,2264.320720795782,2100.586599891512,527.3169730085896))
targetgene="ENST00000493536"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2905.0556339234936,2825.2702649204166,2752.7944802061807,289.7868049867024),c(2694.8289439139207,2848.073091917353,2647.009543520997,106.41351527380547),c(3018.254620851725,3151.350690976607,3102.7596842774665,23.75301680218872),c(2958.1898522775614,3080.661927286104,3031.971267999411,24.703137474276268),c(3064.458288985697,3041.8971213913123,2958.001349866162,419.95333706269656),c(3167.261450583785,3173.0133766236963,3124.2348218000225,251.78197810320043),c(3154.5554418469424,3285.8873702585315,3207.749245498852,569.1222825804417))
targetgene="ENST00000451311"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(3087.560123052683,3230.020444116037,2983.45336470771,31.35398217888911),c(2853.076507272775,2780.804752276391,2736.091595466415,123.51568737138135),c(2343.681066095733,2451.3039021706604,2425.0997891212496,591.9251787105429),c(3296.631721358907,3362.276840698268,3311.9434312564395,358.1954933770059),c(2930.467651397178,3048.737969490393,3028.7897661442175,314.4899424609786),c(3549.5968043924036,3635.910764661505,3527.4901819457996,0.0),c(3435.242725760823,3712.3002351012415,3466.2462712333245,0.0))
targetgene="ENST00000293831"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2437.2434940670264,2390.876410628779,2596.1055138379006,40.8551888997646),c(2555.062847808655,2613.2039738489084,2798.130881642688,131.11665274808175),c(2563.1484897321,2478.667294566984,2324.882480682654,198.5752204662977),c(3122.212874153162,3039.616838691619,3016.0637587234432,45.60579226020234),c(1590.5612755119882,1800.2831914081266,1702.8988679923239,228.97908197309926),c(3615.437031483314,3667.834722457216,3617.367609355016,158.67015223862066),c(2745.6529788612897,2822.989982220723,2803.698509889277,148.21882484565762))
targetgene="ENST00000245458"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2769.909904631625,2878.8569083632174,2768.7019894821483,14.251810081313232),c(3483.7565773014935,3746.5044755966464,3674.634642748499,63.65808502986577),c(1277.5314239043275,1268.9773223795091,1309.1880134121277,61.75784368569067),c(1457.7257296268185,1534.630256893818,1506.441128434125,79.8101364553541),c(2571.2341316555453,2728.358250183437,2657.3494245503757,314.4899424609786),c(3125.67814926321,3272.20567406037,3242.7457659059805,1.9002413441750976),c(3153.400350143593,3180.994366072624,3190.250985295288,3.800482688350195))
targetgene="ENST00000597224"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(3237.7220444880923,3364.557123397962,3316.7156840392295,115.91472199468096),c(2948.9491186507666,2832.1111130194977,2805.2892608168736,2878.865636425273),c(2988.222236564643,2751.1610771803735,2650.1910453761902,51.306516292727636),c(3428.312175540727,3509.355074828508,3306.3758030098506,442.75623319279777),c(2488.0675290143954,2487.7884253657585,2479.9806961233376,387.6492342117199),c(3153.400350143593,3139.949277478139,3129.8024500466113,190.9742550895973),c(3176.5021842105793,3340.614155051179,3220.475252919626,17.102172097575878))
targetgene="ENST00000330752"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2804.5626557321043,2667.9307586415557,2584.9702573447234,109.26387729006811),c(2944.3287518373695,3046.4576867906994,2847.4441603981877,0.0),c(659.5573626124512,697.766506106253,726.1777984479177,11.401448065050586),c(3332.4395641627348,3373.6782541967364,3276.1515353855125,89.31134317622958),c(898.6613452057566,970.2602887196427,938.5430472820842,113.06435997841831),c(2976.67131953115,2986.030195248818,2873.691550703534,85.5108604878794),c(3503.3931362584317,3486.5522478315715,3588.734092658274,437.05550916027244))
targetgene="ENST00000392367"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1162.0222535693974,1167.5047422431423,1081.7106307657918,324.94126985394166),c(3030.9606295885674,3117.1464504812025,3082.875297682507,459.85840529037364),c(3309.337730095749,3513.9156402278954,3416.9329924778253,1043.2324979521286),c(2959.3449439809106,2769.4033387779227,2830.7412756584217,97.86242922501752),c(3015.9444374450263,3065.8400897380957,3016.8591341872416,1158.197099274722),c(2953.5694854641642,3058.9992416390146,3041.5157735649914,550.1198691386908),c(2630.14380852636,2803.6075792733272,2643.828041665803,398.1005616046829))
targetgene="ENST00000359611"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2348.30143290913,2454.7243262202005,2267.6154472891712,1104.0402209657318),c(3319.7335554258925,3323.5120348034766,2989.8163684180968,95.01206720875489),c(3045.9768217321084,3225.45987871665,3130.5978255104096,101.66291191336772),c(2133.4543760861598,2038.5727335261117,2090.246718862133,279.33547759373937),c(2314.803773512,2462.7053156691286,2281.1368301737434,935.8688620062355),c(1768.4453978277807,1849.3092694515399,1710.8526226303077,432.3049057998347),c(2284.7713892249185,2305.3658093902673,2172.965767097164,101.66291191336772))
targetgene="ENST00000225792"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2624.368350009613,2540.2349274587123,2438.6211720058222,2494.066764229816),c(3290.85626284216,3440.946593837699,3238.768888586989,158.67015223862066),c(3255.0484200383316,3248.2627057135865,3298.422048371867,54.15687830899028),c(2951.2593020574654,2976.909064450044,3082.079922218709,925.4175346132726),c(4298.096228162751,3948.309494519533,4259.235608640306,447.5068365532355),c(1853.922183875629,1795.7226260087396,1706.0803698475174,3219.9589577047027),c(2884.263983263206,3027.0752838433036,2886.417558124308,1.9002413441750976))
targetgene="ENST00000579618"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=9
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("ENST00000374080","ENST00000575354","ENST00000016171","ENST00000248846","ENST00000283027","ENST00000564367","ENST00000308388","ENST00000261745","ENST00000472484","ENST00000314133")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='sym-01-01-TP6_vs_Sym-pool-1-1,Sym-pool-1-2,Sym-pool-1-3 pos.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2618.592891492867,2713.5364126354284,2732.9100936112213,2618.5325722732846),c(3380.953415703406,3637.050906011352,3481.358405045494,5911.650821728728),c(2530.8059220383197,2531.1137966599376,2497.478956326902,1306.4159241203797),c(2970.8958610144036,3076.101361886717,3080.4891712911117,4017.110201586156),c(2973.2060444211024,3121.7070158805896,3164.0035949899416,5829.940443929199),c(933.3140963062357,956.5785925214809,988.6517015013819,789.550278504753),c(2275.530655598124,2349.8313220342934,2207.962287504293,4220.4360254128915))
targetgene="ENST00000374080"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2905.0556339234936,3028.2154251931506,3199.7954908608685,3923.048255049489),c(2490.377712421094,2384.035562529698,2573.039625387748,2857.9629816393467),c(2854.231598976124,2989.4506192983586,2768.7019894821483,4636.5888797872385),c(2909.6760007368907,3122.8471572304366,2967.5458554317424,3351.0756104527845),c(2611.662341272771,2615.484256548602,2663.712428260763,3590.506019818847),c(2529.6508303349706,2651.9687797437005,2704.27657691448,3529.698296805244),c(3015.9444374450263,3073.8210791870233,3038.334271709798,4743.952515733131))
targetgene="ENST00000575354"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2357.5421665359245,2501.47012156392,2324.882480682654,2957.7256522085395),c(2972.050952717753,2947.2653893540264,2907.0973201830657,3603.807709228073),c(2221.241345540707,2385.175703879545,2377.377261293347,3131.5977352005607),c(3106.0415903062717,3249.4028470634335,3187.0694834400942,1444.1834215730742),c(2523.875371818224,2675.9117480904833,2743.2499746406,3713.0715865181405),c(3093.3355815694295,3168.4528112243092,3049.469528202975,3906.8962036240005),c(3869.5572062201604,3521.896629676823,3692.1329029520634,4397.1584704211755))
targetgene="ENST00000016171"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1974.0517210239564,2016.9100478790224,2027.4120572220615,3263.66450862073),c(1805.4083323349582,1852.7296935010804,1716.4202508768963,2248.935630831228),c(2765.289537818228,2867.455494864749,2705.071952378278,2074.113427167119),c(3303.5622715790023,3347.4550031502595,3279.333037240706,4057.965390485921),c(2296.3223062584116,2266.6010034954757,2286.7044584203322,3227.559923081403),c(2023.7206642679764,2067.076267272282,1829.3635667362657,3262.7143879486425),c(1776.5310397512258,1836.767714603225,1729.1462582976703,2565.325814636382))
targetgene="ENST00000248846"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(3215.7753021244557,3394.200798493979,3253.0856469353594,3543.950106886557),c(3265.4442453684755,3154.7711150261475,3122.644070872426,3811.884136415246),c(3273.5298872919207,3245.982423013893,3325.4648141410116,3544.9002275586445),c(3177.6572759139285,2946.1252480041794,2841.0811566878006,4204.283973987403),c(2215.4658870239605,2407.9785308764813,2222.279045852664,3207.607388967565),c(2065.3039655885514,2142.3255963621723,1966.9635219733848,2260.3370788962784),c(1883.9545681627108,1942.800860138979,1887.4259755935473,2188.127907817625))
targetgene="ENST00000283027"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2530.8059220383197,2521.992665861163,2562.6997443583687,3561.052278984133),c(883.6451530622157,890.4503942303654,928.9985417165036,1165.7980646514225),c(2739.8775203445434,2875.436484313677,2814.0383909186558,3035.6355473197186),c(2750.273345674687,2815.0089927717954,2678.8245620729317,3409.983092122213),c(2162.3316686698927,2270.021427545016,2336.0177371758314,2505.468212294866),c(2612.81743297612,2577.859592003657,2647.009543520997,4616.6363456734),c(2417.606935110088,2463.845457018975,2284.318332028937,3020.4336165663176))
targetgene="ENST00000564367"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2817.2686644689466,2959.806944202341,2857.7840414275665,2719.245363514565),c(5098.574778583817,4747.5485807621535,4964.733645029466,8162.486693904132),c(3079.474481129238,3079.5217859362574,3094.0105541756843,4560.579226020234),c(3557.6824463158487,3601.7065241661003,3672.248516357104,4515.923554432119),c(2966.2754942010065,2934.7238345057112,2849.8302867895827,4132.07480290875))
targetgene="ENST00000308388"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(3673.191616650779,3809.2122498382214,3620.5491112102095,5591.460155235225),c(3950.4136254546115,4135.292675894411,4215.489958131395,7058.4464729384),c(2879.6436164498086,3008.8330222457544,2996.9747475922823,2491.216402213553),c(2734.1020618277967,2665.650475941862,2689.1644431023105,3619.0096399814734),c(3063.303197282348,2980.329488499584,3064.5816620151445,2658.4376405009616),c(2303.2528564785075,2186.7911090061984,2292.272086666921,3475.5414184962533),c(3032.1157212919165,3162.752104475075,3155.2544648881594,3610.4585539326854))
targetgene="ENST00000261745"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(3504.548227961781,3602.8466655159473,3522.717929163009,3694.0691730763897),c(2649.7803674832976,2821.849840870876,2768.7019894821483,1917.3435162726735),c(4533.734935646009,4585.648509083905,4423.87832964657,4746.802877749394),c(3112.9721405263676,3198.096486320326,3043.9018999563864,2867.464188360222),c(2671.7271098469346,2796.766731174246,2640.6465398106097,5167.7063354841775),c(2361.0074416459724,2536.8145034091717,2428.2812909764434,3559.1520376399576),c(3184.5878261340245,3234.5810095154247,3027.994390680419,1133.4939618004457))
targetgene="ENST00000472484"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2559.683214622052,2609.783549799368,2471.2315660215554,3096.4432703333214),c(2494.9980792344913,2388.596127929085,2341.58536542242,4059.865631830096),c(2647.4701840765993,2771.6836214776163,2558.722867039377,4862.7175997440745),c(1747.6537471674933,1840.1881386527655,1791.980919937742,1800.478673605905),c(2148.470568229701,2564.177895805495,2364.6512538725733,2132.0707881644594),c(1970.5864459139084,1813.9648876062886,1847.6572024036284,2240.38454478244),c(3022.8749876651223,2979.1893471497374,3152.0729630329656,2632.7843823545977))
targetgene="ENST00000314133"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=15
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("ENST00000493536","ENST00000597224","ENST00000245458","ENST00000451311","ENST00000293831","ENST00000392367","ENST00000325495","ENST00000251453","ENST00000579618","ENST00000225792")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='sym-01-02-TP6_vs_Sym-pool-1-1,Sym-pool-1-2,Sym-pool-1-3 neg.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1519.1592782133414,1618.378621513376,1581.5128898800026,356.85767116981754),c(1801.4317410269396,1821.8937926468307,1810.1000052717716,98.1089065753879),c(2422.6534209985984,2476.822933794384,2567.924408163806,478.6852144996948),c(3128.334578032594,3383.980983846698,3193.2457712863747,521.8100085987664),c(2600.462846392991,2761.5276678800146,2736.8464120804356,267.37372341424395),c(2228.1743619734816,2460.585021203949,2379.6306148071963,294.3267197261637),c(2826.0585548621266,2802.6637131091175,2700.4274479332726,1571.89874491116),c(2658.2509096461686,2627.294257132417,2633.0136206821403,502.4038512541842),c(2596.0176107581315,2463.832603722036,2574.898252362199,309.4203976608388),c(2178.165461081309,2149.8996269736217,2046.435836439533,678.1373872079009))
targetgene="ENST00000493536"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2664.9187630984584,2733.3819527232604,2697.3279616228756,101.34326613281827),c(3351.7076686842997,3557.18538481134,3579.906688508384,118.5931837724469),c(1229.1076530387386,1204.8531142102927,1275.4386167283117,23.71863675448938),c(1402.4718427982714,1457.0820231150528,1467.6067679729176,65.7653110010842),c(2473.7736307994865,2590.4883219274307,2588.845940758985,265.2174837092904),c(3007.2019069826642,3106.8539423032703,3159.151421872009,19.406157344582223),c(3033.873320791823,3020.251741820949,3108.009897750461,0.0))
targetgene="ENST00000597224"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2344.8617973885516,2270.0601801428425,2529.1808292838455,20.48427719705901),c(2458.215306077477,2481.1530438185,2725.9982099940466,42.04667424659481),c(2465.994468438482,2353.4147981070764,2264.949621322512,136.92122126455234),c(3003.8679802565193,2886.0183310733514,2938.313022256232,116.43694406749333),c(1530.272367300491,1709.3109320198132,1659.0000476399243,214.54585064288122),c(3478.3968842778045,3482.490986895338,3524.1159349212403,435.5604204006232),c(2641.5812760154445,2680.3381049278387,2731.422311037241,300.79543884102446))
targetgene="ENST00000245458"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2794.941905418108,2682.5031599398967,2681.830530070891,459.2790571551126),c(2592.6836840319866,2704.153710060477,2578.772610250195,518.5756490413361),c(2903.8501784721734,2992.106026664195,3022.774024214547,1.07811985247679),c(2846.062115218996,2924.9893212903958,2953.8104538082166,107.811985247679),c(2948.3025348207716,2888.1833860854094,2881.7473970914893,1019.9013804430434),c(3047.2090276964022,3012.674049278746,3043.6955568097255,523.96624830372),c(3034.984629700538,3119.8442723756184,3125.0570724576432,1064.1042943945918))
targetgene="ENST00000451311"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2970.5287129950707,3066.8004245801967,2906.543287574664,43.12479409907161),c(2744.933004525935,2640.284587204765,2665.5582269413076,51.749752918885925),c(2254.8457757826404,2327.43413796238,2362.583440100014,808.5898893575926),c(3171.6756254724773,3192.3736152795623,3226.565249123141,627.4657541414919),c(2819.3907014098368,2894.6785511215835,2950.7109674978196,443.10725936796075),c(3415.052276481052,3452.180216726526,3436.5554466525286,282.467401348919),c(3305.0326945182715,3524.70955963047,3376.8903351773893,2.15623970495358))
targetgene="ENST00000293831"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2698.258030359907,2533.114364107893,2518.332627197456,66.84343085356099),c(2832.7264083144164,2892.5134961095255,2774.040247805198,47.43727350897876),c(634.5573868762384,662.5068336897566,707.4577503480854,26.95299631191975),c(3206.126201642641,3203.1988903398524,3191.6960281311763,285.70176090634936),c(864.5983309802338,921.2309076306909,914.3484615670765,119.6713036249237),c(2863.843057758435,2835.1395382899877,2799.6110098659724,241.49884695480097),c(3370.599920132454,3310.369113436725,3496.2205581276685,48.51539336145555))
targetgene="ENST00000392367"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(3398.3826428503276,3211.8591103880844,3330.398040521436,2173.489622593209),c(2344.8617973885516,2311.196225371945,2396.677789514379,346.07647264504965),c(2207.059492707898,2225.676552395653,2312.2167875560644,100.26514628034148),c(3297.253532157267,3212.9416378941137,3177.7483397343904,233.95200798746345),c(2544.8974009572435,2651.109862265055,2633.7884922597395,491.6226527294163),c(2801.609758870398,2845.964813350278,2661.6838690533114,154.171138904181),c(2819.3907014098368,2789.6733830367693,2913.517131773057,340.68587338266565))
targetgene="ENST00000325495"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2998.3114357129443,2731.2168977112024,2798.061266710774,7.54683896733753),c(2826.0585548621266,3057.057677025936,3032.8473547233366,112.12446465758617),c(2379.3123735587155,2589.4057944214014,2595.044913379779,294.3267197261637),c(2137.0470314588556,2901.1737161577576,2787.2130646243845,19.406157344582223),c(3268.3595005306784,3345.0099936296533,3209.518074415958,3285.0311904967793),c(3091.6613840450004,2973.7030590617014,3041.370942076928,0.0),c(2273.7380272307946,2291.7107302634226,2195.9860509161817,631.778233551399))
targetgene="ENST00000251453"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2524.8938406003745,2411.8712834326434,2375.7562569192005,2532.50353346798),c(3166.1190809289023,3267.068013195564,3155.277063984013,192.98345359334542),c(3131.668504758739,3084.120864676661,3213.392432303954,40.96855439411802),c(2839.394261766706,2826.4793182417557,3002.627363196967,1063.026174542115),c(4135.180449328342,3748.7927533784755,4149.437298043808,488.38829317198594),c(1783.6507984875004,1704.9808219956972,1662.099533950321,1768.1165580619358),c(2774.938345061239,2874.110528507032,2812.00895510756,18.32803749210543))
targetgene="ENST00000579618"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2259.2910114175006,2330.6817204804674,2209.1588677353684,1319.618699431591),c(3193.9018036467764,3155.567680074576,2912.742260195458,307.2641579558852),c(2930.5215922813322,3062.4703145560807,3049.8945294305195,5.390599262383951),c(2052.587554396519,1935.5591807798771,2036.362505930743,480.84145420464836),c(2227.063053064767,2338.2594130226703,2222.331684554555,678.1373872079009),c(1701.4139392425939,1755.8596147790608,1666.7487634159165,728.8090202743101),c(2198.169021438178,2188.8706171906665,2116.9491500010618,186.51473447848468))
targetgene="ENST00000225792"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=21
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("ENST00000564367","ENST00000559835","ENST00000558494","ENST00000333646","ENST00000575354","ENST00000262306","ENST00000355710","ENST00000314133","ENST00000425967","ENST00000602385")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='sym-01-02-TP6_vs_Sym-pool-1-1,Sym-pool-1-2,Sym-pool-1-3 pos.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2434.8778189944633,2394.5508433361792,2496.636223024678,3805.763079243069),c(850.1513151669394,845.45398220866,905.0500026358858,1585.9143029933582),c(2636.0247314718695,2730.1343702051736,2741.495641546031,3918.965663753132),c(2646.026511650304,2672.760412385636,2609.767473354164,4163.698870265363),c(2080.370277114393,2155.312264503767,2275.7978234089014,3106.063294985632),c(2513.780751513225,2447.594691131601,2578.772610250195,5135.084857346951),c(2325.9695459403974,2339.3419405286995,2225.431170864952,3665.6074984210863))
targetgene="ENST00000564367"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(3233.9089243605144,3154.4851525685467,3268.4083143134985,4850.461216293079),c(851.2626240756543,882.2599174136465,845.3848911607461,1591.3049022557423),c(2104.819073106122,2067.627536515417,2139.420425751439,2365.3949563340775),c(3948.48055266423,4057.313092596745,4084.3480855254743,6908.592014671271),c(3058.322116783552,2910.916463712019,3133.580659811235,4322.182488579451),c(3982.9311288343933,4132.007490512747,4135.489609647023,5743.1444541438605),c(2992.75489116937,3272.4806507257094,3017.349923171352,5027.272872099272))
targetgene="ENST00000559835"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2044.8083920355143,2163.972484551999,1934.0794576876465,3165.359886871856),c(3399.493951759043,3588.5786824861816,3507.8436317916567,4822.430100128682),c(3866.2436934193233,3756.370445920679,3913.1014668760477,4386.869679728059),c(862.3757131628039,855.1967297629211,857.0079648247344,1603.1642206329868),c(2611.5759354801407,2651.109862265055,2503.6100672230714,5859.581398211354),c(4543.03081882673,4363.668376802955,4591.888968852962,6704.827362553157),c(2681.588396729183,2796.1685480729434,2755.443329942817,4704.915036208712))
targetgene="ENST00000558494"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1928.1209566204443,1845.7093977794689,1836.4456389101451,3885.5439483263517),c(2850.5073508538558,2795.086020566914,2777.914605693194,4739.414871487969),c(1828.1031548360984,1845.7093977794689,1971.273293412409,1453.305561138713),c(1897.0043071764255,2074.122701551591,1990.6450828523894,2493.6912187788157),c(2601.574155301706,2527.7017265777476,2762.4171741412097,4684.430759011653),c(925.7203209595563,988.3476130044899,929.8458931190607,805.3555298001622),c(3168.3416987463324,3153.402625062518,3132.8057882336357,4012.7620909186126))
targetgene="ENST00000333646"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2794.941905418108,2875.1930560130613,3117.3083566816513,3554.5611536159768),c(2395.9820071894396,2263.5650151066684,2506.709553533468,3295.8123890215475),c(2746.04431343465,2838.387120808075,2697.3279616228756,5121.069299264753),c(2799.3871410529678,2965.0428390134693,2891.04585602268,4221.91734229911),c(2512.6694426045096,2483.318098830558,2595.044913379779,3314.140426513653),c(2433.766510085748,2517.958979023487,2634.563363837339,3099.5945758707717),c(2901.6275606547433,2918.4941562542217,2960.00942642901,3939.449940950191))
targetgene="ENST00000575354"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1965.9054595167527,1908.4959931291517,2010.0168722923697,4348.057365038894),c(2074.813732570818,2062.2148989852717,2180.4886193641973,1447.914961876329),c(1105.7523641713788,1130.1587162942906,1166.9565958644214,1814.4757117184377),c(1790.31865193979,1677.9176343449717,1704.717470718278,2463.503862909465),c(2239.287451060631,2115.2587467806934,2096.0276174058827,2175.6458622981622),c(1534.7176029353507,1543.684223597374,1435.8370332913496,2624.1437209285073),c(3066.1012791445564,3227.014495472491,3088.63810831048,5909.174911425286))
targetgene="ENST00000262306"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(3545.0754188007018,3630.797255221313,3307.1518931934597,3418.7180522039016),c(2737.15384216493,2737.7120627473764,2699.6525763556733,4557.212616419391),c(2353.7522686582715,2456.254911179833,2262.6250065897148,2046.2714800009476),c(2664.9187630984584,2665.1827198434326,2723.673595261249,4403.04147751521),c(2459.326614986192,2434.6043610592524,2465.6413599207094,3709.810412372635),c(2210.3934194340422,2280.8854552031326,2211.483482468166,3903.871985818457),c(2258.1797025087853,2157.477319515825,2213.0332256233646,2810.658455406992))
targetgene="ENST00000355710"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2462.660541712337,2477.905461300413,2407.525991600768,2740.5806649960005),c(2400.4272428242994,2267.8951251307844,2281.2219244520957,3775.575723373719),c(2547.1200187746736,2631.624367156533,2492.761865136682,5233.193763922339),c(1681.4103788857246,1747.1993947308288,1745.7856643310365,2379.410514416276),c(2067.0345702098134,2434.6043610592524,2303.6932002024732,3066.172860443991),c(1895.8929982677107,1722.3012620921613,1800.0266747629819,3088.8133773460036),c(2908.295414107033,2828.6443732538137,3070.816062025698,1817.7100712758681))
targetgene="ENST00000314133"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1348.0177062712385,1356.4069650543545,1493.952401611291,1630.1172169449067),c(2817.168083592407,2902.256243663787,2840.679203478731,3881.2314689164446),c(2077.036350388248,2356.6623806251637,2319.1906317544576,3209.562800823404),c(2778.2722717873835,2656.5224997952005,2787.9879362019838,3287.187430201733),c(2364.865357745421,2338.2594130226703,2462.5418736103125,3651.591940338888),c(3442.8349991989257,3494.3987894616575,3373.015977289393,3896.3251468511194),c(2053.698863305234,2254.9047950584363,2326.9393475304496,2772.9242605703043))
targetgene="ENST00000425967"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1683.6329967031545,1761.272252309206,1792.2779589869897,2647.8623576829964),c(1790.31865193979,1842.4618152613818,1872.8646030573082,2783.705459095072),c(2078.147659296963,2222.428969877566,2050.310194327529,2554.065930517516),c(2890.514471567594,2853.542505892481,2847.6530476771236,4121.652196018768),c(2028.13875840479,1897.6707180688616,2013.1163586027667,2285.614087250795),c(2130.379178006566,2165.0550120580283,2110.750177380268,2684.5184326672074),c(2898.293633928599,2836.222065796017,2953.8104538082166,4564.759455386729))
targetgene="ENST00000602385"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=27
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("ENST00000493536","ENST00000293831","ENST00000451311","ENST00000597224","ENST00000245458","ENST00000392367","ENST00000325495","ENST00000330752","ENST00000359611","ENST00000225792")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='sym-01-03-TP6_vs_Sym-pool-1-1,Sym-pool-1-2,Sym-pool-1-3 neg.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1552.7739700632592,1678.1202308749926,1614.9022222385577,59.776098831812156),c(1841.2923229499218,1889.1480592392056,1848.3153312833272,52.68401930939377),c(2476.2598791060022,2568.2535707304232,2622.1391300825967,292.8015688541307),c(3197.555761322659,3508.898890779416,3260.662448723712,445.78785569487036),c(2658.0037234440574,2863.4680327505725,2794.6274615122907,59.776098831812156),c(2277.477549361254,2551.416243999236,2429.8700266999563,126.64427718604271),c(2888.591225947965,2906.122593802914,2757.4396102407513,1576.4679624118596),c(2717.0704728539254,2724.2794651060917,2688.6025238444972,491.3797954818457),c(2653.460127335606,2554.7837093454737,2629.2602079856574,457.94570630473044),c(2226.3620931411765,2229.2620592091876,2089.640749109275,673.7475546297472))
targetgene="ENST00000493536"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(3036.258099472635,3180.0097752968923,2967.907023820103,71.93394944167225),c(2805.6705969687273,2737.7493264910413,2721.8342207254473,121.578506098601),c(2304.7391260119625,2413.3501648035012,2412.4629473813634,449.8404725648237),c(3241.85582338006,3310.2184353514067,3294.68537648278,873.3389354749505),c(2881.775831785288,3001.534111946308,3013.007183872821,596.7478341006332),c(3490.617710317773,3579.6156630504024,3509.1089444527206,1.0131542174883417),c(3378.163706633601,3654.8223891163716,3448.1841668376455,0.0))
targetgene="ENST00000293831"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2856.7860531888055,2781.526375992128,2738.4500691659227,273.55163872185227),c(2650.0524302542676,2803.9761449670445,2633.2163623762467,281.656872461759),c(2968.1041578458644,3102.558072333431,3086.591655537782,6.07892530493005),c(2909.0374084359964,3032.9637885111906,3016.1721073852923,0.0),c(3013.540118930378,2994.799181253833,2942.587635720331,230.9991615873419),c(3114.635132343421,3123.8853528596014,3107.9548892469643,253.28855437208543),c(3102.14024304518,3235.011709285437,3191.03413144934,864.2205475175555))
targetgene="ENST00000451311"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2723.8858670166023,2834.2833330831813,2754.27468672828,0.0),c(3425.871465772341,3688.4970425787465,3655.4866569045253,92.1970337914391),c(1256.3043239868066,1249.3296434540914,1302.366025382002,39.51301448204533),c(1433.5045722164105,1510.8694520118663,1498.5912831552318,26.342009654696884),c(2528.5112343531932,2686.114857848734,2643.502363791779,354.6039761209196),c(3073.7427673673587,3221.5418479004875,3225.848290086526,4.052616869953367),c(3101.004344018067,3131.7427720008222,3173.627052130747,0.0))
targetgene="ENST00000597224"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2396.746947208103,2353.858277019973,2582.577586176704,140.8284362308795),c(2512.608647973613,2572.743524525407,2783.5502292186407,331.30142911868774),c(2520.5599411634034,2440.2898875734004,2312.7678567385124,446.80100991235867),c(3070.3350702860203,2992.5542043563414,3000.347489822935,174.26252540799476),c(1564.1329603343877,1772.409260569641,1694.025310050344,265.44640498194553),c(3555.363954863205,3611.0453396152852,3598.5180336800395,314.0778074213859),c(2700.031987447233,2779.2813990946365,2789.0888453654657,165.1441374505997))
targetgene="ENST00000245458"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2757.9628378299876,2626.6229700652057,2571.500353883054,92.1970337914391),c(2895.406620110642,2999.289135048816,2832.6065436619483,6.07892530493005),c(648.5983444814345,686.9629306324384,722.3937917216086,38.499860264556986),c(3277.068693220558,3321.443319838865,3259.079986967476,272.5384845043639),c(883.7294430937935,955.2376698826881,933.6524361790779,54.71032774437045),c(2927.2117928698017,2939.797247265288,2858.7171626398376,7.092079522418391),c(3445.1817492332593,3432.569676264701,3570.033722067796,49.64455665692874))
targetgene="ENST00000392367"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(3473.57922491108,3330.4232274288315,3400.7103141505736,2080.0056085035653),c(2396.746947208103,2396.512838072314,2447.277106018549,249.23593750213206),c(2255.8954678461105,2307.8362506213944,2361.032940303702,37.486706047068644),c(3370.2124134438113,3331.545715877577,3244.8378311613546,64.84186991925387),c(2601.208772088415,2748.9742109784997,2689.393754722615,413.3669207352434),c(2863.6014473514824,2951.022131752746,2717.878066334858,167.17044588557638),c(2881.775831785288,2892.6527324179638,2975.0281017231637,106.38119283627587))
targetgene="ENST00000325495"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(3183.924972997305,3312.4634122488983,3299.4327617514873,226.94654471738852),c(2899.9502162190934,2788.261306684603,2790.6713071217014,2253.254979694072),c(2938.5707831409304,2708.56462682365,2636.381285888718,63.828715701765525),c(3371.348312470924,3455.0194452396167,3289.146760335955,628.1556148427718),c(2446.7265044010683,2449.2697951633672,2467.0578779714956,392.09068216798823),c(3101.004344018067,3091.333187845973,3113.4935053937893,317.11727007385093),c(3123.722324560324,3288.8911548252363,3203.693825499226,3.039462652465025))
targetgene="ENST00000330752"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1142.7144212755222,1149.428171515714,1076.0739942402931,359.6697472083613),c(2980.5990471441055,3068.8834188710566,3066.8108835848357,792.2865980758831),c(3254.350712678301,3459.5093990346,3399.127852394338,1737.559482992506),c(2910.1733074631093,2726.5244420035833,2815.9906952214733,209.72292302008674),c(2965.8323597916383,3018.371438677495,3001.138720701053,1123.588027194571),c(2904.4938123275447,3011.63650798502,3025.6668779227066,1126.627489847036),c(2586.4420847359484,2760.1990954659577,2630.0514388637753,334.34089177115277))
targetgene="ENST00000359611"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2309.282722120414,2416.7176301497384,2255.7992335140266,1180.3246633739182),c(3264.5738039223165,3272.053828094049,2974.236870845046,153.99944105822794),c(2995.3657344965723,3175.519821501909,3114.284736271907,128.6705856210194),c(2098.0055030774247,2007.0093463575163,2079.354747693743,295.84103150659575),c(2276.3416503341414,2424.575049290959,2269.25015844203,840.9180005153236),c(1739.061410509766,1820.676263865711,1701.9376188315225,822.6812246005335),c(2246.8082756292074,2269.6716433640368,2161.6427590180006,89.15757113897406))
targetgene="ENST00000225792"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=33
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("ENST00000575354","ENST00000348750","ENST00000215909","ENST00000425967","ENST00000563280","ENST00000275036","ENST00000222005","ENST00000374080","ENST00000563694","ENST00000262306")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='sym-01-03-TP6_vs_Sym-pool-1-1,Sym-pool-1-2,Sym-pool-1-3 pos.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2856.7860531888055,2981.329319868883,3183.1218226681613,5368.704198470722),c(2448.998302455294,2347.123346327498,2559.6318907112855,3664.5788046553316),c(2806.80649599584,2943.1647126115254,2754.27468672828,4557.167670262561),c(2861.3296492972568,3074.495861114786,2952.0824062577453,4145.827057962294),c(2568.2677003021427,2574.9885014228985,2649.832210816722,4406.207691856798),c(2487.6188693771305,2610.9081317827645,2690.184985600733,4172.169067616991),c(2965.8323597916383,3026.228857818716,3022.501954410235,5392.006745472954))
targetgene="ENST00000575354"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(3722.3411118487934,3794.010956760853,3648.3655790014645,5872.241844562428),c(3328.1841494406362,3639.1075508339304,3489.328172499774,6405.160962961296),c(1840.156423922809,1907.1078744191386,1921.8998029482884,3276.540739357297),c(843.9729771448439,906.9706665866181,871.9364276858846,1254.284921250567),c(1548.2303739548079,1471.5823563057627,1485.1403582272283,1470.0867695755837),c(3254.350712678301,3263.0739205040827,3198.9464402305184,4187.366380879316),c(2378.5725627742977,2379.675511341127,2199.621841167658,3159.0148501286494))
targetgene="ENST00000348750"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2685.2653000947657,2887.0402901742345,2759.813302875105,3231.9619537878098),c(2952.2015714662843,3305.7284815564235,3138.0216626154433,1394.100203263958),c(4016.5389598710203,4030.8560194462198,3892.064689461766,2794.2793318328463),c(3106.683839153631,3229.3992670417083,3188.6604388149863,5990.780888008564),c(2584.1702866817227,2669.277531117547,2570.709123004936,3455.8690358527333),c(2070.7439264267164,1993.5394849725665,2018.4299700786676,3410.277096065758),c(2024.17206631509,2127.1156103733183,2115.751368087165,3947.2488313345793))
targetgene="ENST00000215909"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1377.845519887881,1406.4780262785055,1525.4931330112393,1843.9406758287819),c(2879.5040337310625,3009.3915310875286,2900.6523991800846,2817.581878835078),c(2122.9952816739074,2443.657352919638,2368.154018206763,3025.278493420188),c(2839.7475677821126,2754.5866532222285,2846.8486994680698,3959.4066819444392),c(2417.1931296961343,2424.575049290959,2514.5317306585675,3211.698869438043),c(3519.015185995594,3623.392712551489,3444.2280124470562,4737.5091209754855),c(2099.1414021045375,2338.1434387375316,2376.0663269879415,3075.9362042946054))
targetgene="ENST00000425967"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2577.3548925190453,2656.930158181343,2631.633900620011,3296.8038237070637),c(2730.7012611792798,2782.648864440874,2858.7171626398376,3447.763802112827),c(2781.8167173993575,2768.0565146071785,2750.3185323376906,4783.101060762461),c(2109.3644933485534,2218.037174721729,2164.016451652354,1636.2440612436717),c(3179.3813768888535,3149.7025871807555,3376.18215692892,3996.893387991508),c(1917.3975577664826,1999.1519272162955,1969.3736556353601,2214.7551194295147),c(3022.6273111472806,3277.6662703377783,3380.929542197627,4697.9961064934405),c(2309.282722120414,2488.5568908694704,2352.3294006444057,2833.792346314892),c(2861.3296492972568,2815.2010294545025,2833.397774540066,3784.1310023189562),c(2063.9285322640394,2148.4428908994887,2195.665686777069,2968.541857240841))
targetgene="ENST00000563280"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1107.5015514350239,1013.6070692174704,1051.5458370186395,1270.4953887303805),c(2327.4571065542195,2206.812290234271,2283.4923142481516,3147.8701537362776),c(1023.4450234286734,1084.3238414884568,1083.9863030214717,1571.402191324418),c(2098.0055030774247,2038.439022922399,2152.1479884805863,2571.3854039854114),c(2255.8954678461105,2301.1013199289196,2239.974615951669,2560.2407075930396),c(1508.4739080058582,1554.646501512953,1469.3157406648709,1993.8875000170565),c(2847.6988609719024,2976.8393660739,3069.1845762191892,3512.6056720320807))
targetgene="ENST00000275036"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2140.0337670806002,2250.589339735358,2153.730450236822,2697.0165269539657),c(2458.0854946721965,2422.3300723934676,2490.794804315032,2771.9899390481028),c(2743.196150477521,2692.849788541209,2713.9219119442687,3446.7506478953383),c(2216.1390018971606,2622.1330162702225,2277.1624672232088,3266.4091971824137),c(1776.5460784044897,1864.4533133667978,1876.0084120174522,2327.215237570721),c(2559.18050808524,2806.221121864536,2760.6045337532228,4526.7730437379105),c(2646.644733172929,2707.4421383749045,2783.5502292186407,3302.882749011994))
targetgene="ENST00000222005"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2575.0830944648196,2671.5225080150385,2718.669297212976,2977.660245198236),c(3324.7764523592973,3580.738151499148,3463.217553521885,7218.723799604434),c(2488.7547684042434,2491.924356215708,2484.4649572900885,1405.2448996563298),c(2921.5322977342375,3028.4738347162074,3064.437190950482,3041.4889609000015),c(2923.804095788463,3073.37337266604,3147.5164331528576,3075.9362042946054),c(917.8064139071788,941.7678084977383,983.4999815005033,741.6288872014661),c(2237.721083412305,2313.4486928651236,2196.4569176551868,3838.8413300633265))
targetgene="ENST00000374080"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(3145.304406075468,3133.987748898314,3352.445230585384,4516.641501563027),c(1569.8124554699518,1635.4656698226515,1593.5389885293753,2713.226994433779),c(2891.9989230293036,3091.333187845973,2872.959318445959,2279.5969893487686),c(1752.69219883512,1644.4455774126181,1592.7477576512574,883.4704776498339),c(1559.5893642259364,1480.5622638957293,1590.3740650169038,2442.7148183643917),c(2391.067452072539,2383.042976687364,2383.97863576912,3477.1452744199887),c(2397.8828462352158,2326.918554250073,2414.045409137599,2047.5846735439386))
targetgene="ENST00000563694"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2009.405378962623,1978.9471351388709,2052.4528978377357,2913.8315294964705),c(2120.7234836196817,2138.3404948607763,2226.5236910236654,1122.5748729770826),c(1130.219531977281,1171.8779404906302,1191.593702445501,1676.7702299432055),c(1829.9333326787935,1739.8570955560124,1740.7079318592978,1048.6146151004336),c(2288.836539632383,2193.3424288493215,2140.2795253088184,3004.002254852933),c(1568.676556442839,1600.6685279115313,1466.1508171523994,2642.306199209595),c(3133.9454158043395,3346.1380657112727,3153.8462801778005,4950.271506648038))
targetgene="ENST00000262306"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-01-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=39
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("ENST00000493536","ENST00000293831","ENST00000196551","ENST00000245458","ENST00000355968","ENST00000330752","ENST00000359611","ENST00000392367","ENST00000251453","ENST00000443387")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='sym-02-01-TP6_vs_Sym-pool-1-1,Sym-pool-1-2,Sym-pool-1-3 neg.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1340.3133910562644,1448.6552111878425,1388.4496708278873,14.084791169447252),c(1589.354796563427,1630.8272377452433,1589.1320093355928,112.67832935557801),c(2137.441984274072,2217.07232320922,2254.444982422155,115.80828294878852),c(2760.045498041978,3029.094441587422,2803.430227085607,104.85344537255176),c(2294.318460184096,2471.9193603613285,2402.7458291837816,26.604605542289253),c(1965.8583387474837,2202.537321090278,2089.1371578208928,101.72349177934126),c(2493.3554889949382,2508.7413657293137,2370.772710845266,921.7713332004923),c(2345.3033148548534,2351.763342844745,2311.5884279633324,125.19814372842),c(2290.396548286345,2205.4443215140664,2260.567494444424,31.299535932105),c(1921.7368298977894,1924.434280547863,1796.6171389791527,90.7686542031045))
targetgene="ENST00000493536"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2620.8176256718325,2745.17740019743,2551.7269550589936,25.039628745684002),c(2421.78059686099,2363.391344539898,2340.160150733921,39.124419915131256),c(1989.3898101339873,2083.3503037149576,2074.171017322013,436.6285262528648),c(2798.2841390450467,2857.5814165839115,2832.682228969781,159.62763325373552),c(2487.472621148312,2591.1063777366494,2590.5028645333637,143.97786528768302),c(3013.0088154468913,3090.1414504869763,3017.037868751436,0.0),c(2915.9414959775645,3155.064459951582,2964.6563770053567,0.0))
targetgene="ENST00000293831"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1727.6021909591352,1719.0062506001557,1724.5075529390958,37.559443118526005),c(2164.895367558326,2094.009305268848,2125.8722299545066,31.299535932105),c(1665.8320785695635,1762.6112569569802,1648.9965713311117,53.209211084578506),c(2022.7260612648672,1948.6592840794324,1970.0883129434403,151.80274927070926),c(2660.0367446493383,2783.9374058479407,2685.0616613217403,43.819350304947),c(1809.9623408118978,1852.728270094418,1803.419930115007,61.03409506760475),c(2617.8761917485194,2711.262395253233,2752.409293566699,15.6497679660525))
targetgene="ENST00000196551"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2068.8085260634366,2031.9932962280307,2220.431026742883,1.5649767966052501),c(2168.817279456077,2220.948323774271,2393.221921593585,7.82488398302625),c(2175.6806252771403,2106.606307105264,1988.4558490102474,154.93270286391976),c(2650.231964904962,2583.3543766065472,2579.6183987159966,39.124419915131256),c(1350.1181708006409,1530.0512230539152,1456.4775821864316,1.5649767966052501),c(3068.8960599898373,3117.273454442334,3093.909408586591,25.039628745684002),c(2330.5961452382885,2399.244349766621,2397.9838753886834,40.68939671173651))
targetgene="ENST00000245458"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2400.2100814233618,2336.2593405845405,2310.2278697361617,538.352018032206),c(2082.535217705564,2104.6683068227385,2165.3284185424623,10.95483757623675),c(2772.7917117096677,3000.9934374908016,2850.3694859230027,17.21474476265775),c(1844.2790699172153,1957.3802853507973,1883.6928655180893,50.079257491368004),c(2916.921973952002,2933.1634276024074,2719.0756170010122,640.0755098115472),c(2841.4251699203032,2992.2724362194367,2976.2211219363094,25.039628745684002),c(2047.2380106258086,2189.940319253862,2100.02162363826,15.6497679660525))
targetgene="ENST00000355968"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2748.2797623487268,2859.519416866437,2836.7639036512937,6.2599071864210005),c(2503.1602687393147,2406.996350896723,2399.3444336158545,1294.2358107925418),c(2536.496519870195,2338.197340867066,2266.6900064666934,18.779721559263002),c(2910.0586281309384,2982.582434806809,2827.920275174683,250.39628745684),c(2111.949556938693,2114.358308235366,2121.1102761594084,123.63316693181476),c(2676.704870214778,2668.626389037671,2676.8983119587147,28.169582338894504),c(2696.3144297035315,2839.170413899919,2754.4501309074553,0.0))
targetgene="ENST00000330752"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(986.3608422842736,992.2561446530774,925.1795944762013,89.20367740649925),c(2572.7742049243875,2649.2463862124155,2636.761844257174,247.2663338636295),c(2809.069396763861,2986.45843537186,2922.4790719630596,236.31149628739277),c(2511.9845705092534,2353.7013431272703,2421.1133652505887,126.76312052502526),c(2560.027991256698,2605.641379855591,2580.298677829582,741.7990015908886),c(2507.0821806370654,2599.8273790080143,2601.3873303507307,123.63316693181476),c(2232.548347794524,2382.771347365154,2261.2477735580096,97.0285613895255))
targetgene="ENST00000359611"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2380.600521934609,2267.460330554884,2210.907119152687,9.389860779631501),c(2499.238356841564,2589.168377454124,2435.3992266358828,1.5649767966052501),c(559.8529234038969,593.0280864528158,621.0948307035087,23.47465194907875),c(2828.6789562526137,2867.271417996539,2802.069668858436,90.7686542031045),c(762.81186411249,824.6191202146181,802.7293540308217,106.41842216915701),c(2526.6917401258183,2537.811369967197,2457.8484373842025,29.73455913549975),c(2973.7896964693855,2963.202431981553,3069.419360497515,0.0))
targetgene="ENST00000392367"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2645.3295750327734,2444.787356405971,2456.4878791570313,29.73455913549975),c(2493.3554889949382,2736.456398926065,2662.6124505734206,9.389860779631501),c(2099.2033432710036,2317.848337900548,2278.2547513976456,118.938236541999),c(1885.4591448435965,2596.920378584226,2446.9639715668354,350.554802439576),c(2883.585722821122,2994.2104365019622,2817.7160884709015,2962.5010759737384),c(2727.689724885536,2661.8433880488315,2670.09552082286,4.694930389815751),c(2006.0579356994272,2051.373299053286,1927.911007901143,7.82488398302625))
targetgene="ENST00000251453"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2107.047167066505,2240.3283265995265,2130.6341837496047,234.7465194907875),c(2537.4769978446325,2489.3613629040583,2566.0128164442876,982.8054282680971),c(3672.8704922434285,3615.3395270513984,3626.567954523992,64.16404866081525),c(2742.3968945021006,2819.790411074663,2770.09655051992,336.4700112701288),c(1908.9906162301002,2169.591316287344,2010.2247806449816,23.47465194907875),c(2131.559116427446,2074.6293024435927,2138.117253999045,679.1999297266785),c(2986.535910137075,2962.2334318402905,3139.4881091968155,78.24883983026251))
targetgene="ENST00000443387"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=45
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("ENST00000330436","ENST00000258704","ENST00000325656","ENST00000445274","ENST00000432701","ENST00000458333","ENST00000400385","ENST00000314133","ENST00000419003","ENST00000376708")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='sym-02-01-TP6_vs_Sym-pool-1-1,Sym-pool-1-2,Sym-pool-1-3 pos.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1791.3332592975823,1737.4172532841483,1750.3581592553426,1156.5178526912798),c(2263.9236429765288,2161.839315157242,2385.7388513441456,3607.2715161751016),c(2360.9909624458555,2247.111327588366,2293.2208918965252,2882.6872593468706),c(2387.4638677556723,2190.909319395125,2252.404145081399,2837.3029322453185),c(2686.5096499591546,2817.8524107921376,2629.278774007734,3358.4402055148666),c(2077.6328278333754,2040.7142974993956,2060.565435050304,2519.612642534453),c(1564.8428472024857,1544.5862251728568,1501.3760036830708,1924.9214598244578))
targetgene="ENST00000330436"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2385.502911806797,2270.3673309786723,2432.678110181541,3663.6106808528907),c(1910.9515721789753,2229.669325045636,2227.914096992323,2084.549093078193),c(2928.6877096452536,2945.7604294388234,2805.4710644263632,3147.168337973158),c(2447.2730241963686,2539.7493702497227,2445.6034133396643,3641.701005700417),c(2085.476651628877,2140.521312049461,2172.1312096783167,2727.754556482951),c(2984.5749541881996,2979.6754343830207,3012.2759149563376,3895.2272467504677),c(2518.847916330317,2719.014396383335,2687.1024986624966,3957.826318614678))
targetgene="ENST00000258704"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2229.606913871211,2200.5993208077525,2043.5584572106682,2723.059626093135),c(2919.863407875315,2698.665393416817,2843.566694787148,3444.5139293281554),c(1763.8798760133282,1863.3872716483083,1693.2147137141656,2222.2670511794554),c(1833.4938121984012,1835.286267551688,1780.2904402531021,2050.119603552878),c(2677.685348189216,2593.044378019175,2712.9531049787433,2912.4218184823703),c(2120.773858708632,2199.6303206664898,2088.456878707307,2486.7481298057423),c(1639.359173259747,1687.9982460797469,1699.3372257364344,1982.825601298852))
targetgene="ENST00000325656"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2100.1838212454413,2217.07232320922,2225.1929805379814,2939.0264240246597),c(2906.136716233188,2932.1944274611446,2866.696184649053,3951.5664114282567),c(1986.4483762106743,2190.909319395125,2066.007667958988,3372.524996684314),c(1910.9515721789753,2043.6212979231839,1883.0125864045037,2472.663338636295),c(3291.4645601871835,3436.074500917786,3309.557887593176,3754.379335055995),c(2545.3208216401335,2684.1303912978756,2642.204077165857,2834.172978652108),c(2130.5786384530084,2223.8553241980594,1997.9797566004436,2438.23384911098))
targetgene="ENST00000445274"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2063.9061361912486,2096.9163056926363,1951.040497763048,2303.6458446029283),c(3070.8570159387127,3129.87045627875,2878.2609295800057,4208.2226060715175),c(2826.7180003037383,2851.767415736335,2607.5098423729996,3975.0410633773354),c(2460.999715838496,2656.029387201255,2475.535694337424,2652.635670245899),c(2190.387794893705,2263.584329989833,2284.3772634199145,2632.2909718900305),c(2301.181806005159,2435.0973549933433,2484.379322814035,2856.0826538045817),c(3190.475328820106,3252.9334742191218,3168.7401110809897,2856.0826538045817),c(4132.714662254685,3887.6285667462366,3938.1357885461243,4169.098186156386),c(2179.602537174891,2076.567302726118,2208.8662818119305,2325.5555197554017),c(2071.7499599867497,2216.1033230679573,2039.4767825291556,2416.324173958506))
targetgene="ENST00000432701"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1422.6735409090268,1415.7092063849084,1423.143905620745,1701.1297779099068),c(2050.1794445491214,2113.3893080941034,2111.5863685692125,2230.0919351624816),c(2440.4096783753052,2453.508357677336,2578.938119602411,4011.035529699256),c(2838.48373599699,3099.831451899604,2968.7380516868693,3378.784903870735),c(1908.9906162301002,1964.1632863396367,1957.8432888989025,2449.1886866872164),c(2307.0646738517853,2330.445339736964,2334.037638711652,2735.5794404659773),c(1929.5806536932907,2100.7923062576874,2092.53855338882,2765.313999601477))
targetgene="ENST00000458333"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2314.9084976472864,2518.4313671419413,2395.943038047927,3134.648523600316),c(3264.9916548773667,3495.1835095348147,3308.8776084795904,3799.7636621575475),c(2096.261909347691,2072.691302161067,2114.3074850235544,2748.0992548388194),c(1955.0730810286695,1974.8222878935271,1930.6321243554848,1973.4357405192204),c(3256.167353107428,2933.1634276024074,3136.0867136288884,3780.9839405982843),c(2963.0044387505714,2834.325413193605,2868.0567428762242,3449.2088597179713),c(2044.2965767024955,2196.7233202427014,2125.1919508409214,2082.984116281588),c(1926.6392197699777,2074.6293024435927,1928.5912870147286,2995.365588702449),c(2618.856669722957,2809.1314095207726,2659.8913341190787,2496.137990585374),c(1817.8061646073988,1829.4722667041115,1878.2506326094058,2630.7259950934254))
targetgene="ENST00000400385"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2172.739191353827,2218.041323350483,2113.627205909969,1873.2772255364844),c(2117.8324247853193,2030.0552959455051,2002.7417103955415,2413.1942203652957),c(2247.2555174110885,2355.639343409796,2188.457908404367,4594.771874833014),c(1483.463175324161,1563.9662279981121,1532.6688429080011,2358.4200324841117),c(1823.6890324540248,2179.2813176999716,2022.4698046895194,1754.3389889944854),c(1672.695424390627,1541.6792247490685,1580.288380858982,748.0589087773095),c(2565.910859103324,2531.9973691196205,2695.9461271391074,1344.3150682839098))
targetgene="ENST00000314133"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2626.700493518458,2845.9534148887583,2685.7419404353254,3067.3545213462903),c(2719.8459010900347,2802.3484085319333,2831.32167074261,4785.699044018855),c(2297.2598941074084,2355.639343409796,2304.785636827478,2154.9730489254293),c(1970.760728619672,1917.6512795590236,1804.780488342178,2516.482688941242),c(2005.0774577249897,2081.412303432432,1976.2108249657094,3226.9821546000258),c(2551.2036894867592,2760.6814024576342,2658.530775891908,3366.265089497893),c(2278.630812593093,2209.3203220791174,2202.7437697896617,2985.9757279228174),c(2014.8822374693661,2070.7533018785416,1943.5574275136082,1696.4348475200911),c(2532.574607972444,2657.9673874837804,2602.7478885779014,2854.517677007976),c(2702.197297550157,2960.295431557765,2989.1464250944327,3403.824532616419))
targetgene="ENST00000419003"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2931.6291435685666,3142.467458115166,2876.9003713528346,2632.2909718900305),c(1749.1727063967635,1771.3322582283452,1817.7057915003013,2707.4098581270828),c(1348.1572148517657,1341.096195507675,1329.265387945954,1588.4514485543289),c(2148.2272419928863,2311.0653369117085,2139.4778122262155,2680.8052525847934),c(2257.060297155465,2343.04234157338,2262.6083317851803,3192.55266507471),c(3470.892029509273,3570.7655205533106,3524.5260874861756,3931.2217130723884),c(2451.1949360941194,2603.7033795730654,2551.046675945408,3241.066945769473))
targetgene="ENST00000376708"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=51
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("ENST00000493536","ENST00000196551","ENST00000245458","ENST00000293831","ENST00000330752","ENST00000355968","ENST00000316292","ENST00000359611","ENST00000251453","ENST00000225792")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='sym-02-02-TP6_vs_Sym-pool-1-1,Sym-pool-1-2,Sym-pool-1-3 neg.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1291.462540996714,1388.7345919997804,1336.4554426954833,27.98490735458799),c(1531.4270511745965,1563.371450391726,1529.6226918846883,43.726417741543735),c(2059.5379219991487,2125.36772340836,2170.0212332644933,148.6698203212487),c(2659.449197443855,2903.8022304958617,2698.448250538014,57.71887141883773),c(2210.6966685678935,2369.673541265177,2312.768556394143,19.23962380627924),c(1894.2080429395842,2111.433931515385,2010.9038042713516,59.46792812849948),c(2402.479328276989,2404.9724807273788,2281.9927573707782,1026.6962885714468),c(2259.823261202736,2254.4875282832554,2225.0247889658262,136.42642335361646),c(2206.917699903675,2114.22068989398,2175.914471375351,78.70755193477872),c(1851.6946454671247,1844.8340466298087,1729.3379834192901,62.966041547822975))
targetgene="ENST00000493536"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1664.6356965883028,1647.903121209104,1659.9287345580842,36.73019090289674),c(2085.9907026486794,2007.3949520478432,2046.2632329364942,3.4981134193234986),c(1605.1169401268596,1689.704496888027,1587.2454645241799,15.741510386955744),c(1949.0030885707542,1868.0570331180993,1896.313063226908,38.47924761255848),c(2563.0854965062804,2668.7856072343607,2584.51231372811,0.0),c(1743.9940385368939,1776.0940066244684,1735.8860257646868,0.0),c(2522.46158336593,2599.1166477694887,2649.337932947538,0.0))
targetgene="ENST00000196551"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1993.405970375323,1947.944106637819,2137.28102153751,19.23962380627924),c(2089.769671312898,2129.083401246486,2303.601297110588,29.73396406424974),c(2096.3828664752805,2019.4709050217543,1913.9927775594792,75.20943851545522),c(2553.6380748457336,2476.499279111314,2483.0176573744598,8.745283548308747),c(1300.9099626572606,1466.763826600437,1401.935866149451,8.745283548308747),c(2957.0429797510715,2988.33390131324,2978.0496586864565,45.475474451205486),c(2245.652128711916,2300.004581800305,2308.1849267523658,55.96981470917598))
targetgene="ENST00000245458"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2525.295809864094,2631.6288288530955,2456.170683758333,24.48679393526449),c(2333.513150154999,2265.634561797635,2252.526566816493,50.72264458019073),c(1916.8818549248958,1997.1768379929952,1996.4981111114787,426.76983715746684),c(2696.2941419199865,2739.3834861587643,2726.60483262322,122.43396967632245),c(2396.8108752806606,2483.930634787567,2493.4945251270947,206.38869174008641),c(2903.192676285956,2962.324156446354,2904.0567801834727,26.235850644926238),c(2809.6632018465452,3024.5617602349735,2853.6368541239176,0.0))
targetgene="ENST00000293831"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2648.1122914511993,2741.2413250778272,2730.533658030458,80.45660864444046),c(2411.9267499375355,2307.435937476558,2309.494535221445,1166.6208253443867),c(2444.0479835833935,2241.482655849813,2181.807709486208,29.73396406424974),c(2803.9947488502175,2859.214096438344,2722.021202981442,227.3773722560274),c(2034.9746256817277,2026.9022606980072,2041.6796032947166,83.95472206376397),c(2579.146113329209,2558.244191550097,2576.6546629136337,68.21321167680823),c(2598.0409566503026,2721.7340164276634,2651.302345651157,0.0))
targetgene="ENST00000330752"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2312.7288225017965,2239.6248169307496,2223.715180496747,594.6792812849948),c(2006.6323607000882,2017.613066102691,2084.241878539796,15.741510386955744),c(2671.7308456025658,2876.8635661694448,2743.6297427212517,27.98490735458799),c(1777.0600143488068,1876.417308253884,1813.1529254403688,78.70755193477872),c(2810.6079440126,2811.839204002231,2617.2525254550937,746.8472150255669),c(2737.8627972263916,2868.50329103366,2864.768526111092,52.471701289852476),c(1972.6216427221207,2099.357978541474,2021.3806720239866,29.73396406424974))
targetgene="ENST00000355968"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2444.9927257494483,2509.011460194921,2390.035456069825,3.4981134193234986),c(2294.778721346758,2357.597588291266,2215.857529682271,69.96226838646997),c(2149.288427774341,2118.8652871916383,2148.412693524684,38.47924761255848),c(2090.7144134789523,2094.713381243816,2108.469635217764,6.996226838646997),c(1943.3346355744263,2099.357978541474,2020.7258677894467,99.6962324507197),c(2386.4187114540596,2452.3473731634917,2511.8290436942057,1990.4265355950706),c(2733.1390863961183,2699.4399493989044,2595.643985715284,208.13774844974816))
targetgene="ENST00000316292"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(950.4106190509834,951.2135265603847,890.5337589739623,82.20566535410222),c(2479.0034437274157,2539.6658023594646,2538.0212130757927,190.64718135313066),c(2706.686305746588,2862.92977427647,2813.0389915824576,236.12265580433615),c(2420.429429432027,2256.345367202319,2330.448270726715,17.490567096617493),c(2466.721795568705,2497.8644266805413,2483.672461608999,552.7019202531128),c(2415.705718601754,2492.2909099233516,2503.971392879729,160.91321728888093),c(2151.1779121064505,2284.2129509882675,2176.56927560989,111.93962941835196))
targetgene="ENST00000359611"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2548.9143640154603,2343.6637963982917,2364.498090922778,27.98490735458799),c(2402.479328276989,2623.268553717311,2562.9037739883006,24.48679393526449),c(2022.6929775230174,2221.9753471996487,2192.939381473382,173.1566142565132),c(1816.7391853231024,2489.504151544757,2355.330831639222,307.83398090046785),c(2778.4867103667416,2870.3611299527233,2712.199139463347,2630.581291331271),c(2628.2727059640515,2551.7417553333757,2570.106620568237,0.0),c(1932.942471747825,1966.5224958284516,1855.715200685448,0.0))
targetgene="ENST00000251453"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1920.6608235891144,1999.96359637159,1866.8468726726226,620.915131929921),c(2715.1889852410795,2707.800224534689,2461.40911763465,69.96226838646997),c(2491.2850918861263,2627.913151014969,2577.309467148173,31.483020773911488),c(1744.9387807029486,1660.9079936425467,1720.8255283702742,131.1792532246312),c(1893.2633007735294,2006.4660325883115,1877.978544659797,472.2453116086723),c(1446.4002562296776,1506.7073633602968,1408.4839084948478,269.35473328790937),c(1868.7000044561084,1878.275147172947,1788.9251687624007,68.21321167680823))
targetgene="ENST00000225792"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=57
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("ENST00000355710","ENST00000286307","ENST00000368748","ENST00000273860","ENST00000414069","ENST00000376708","ENST00000374597","ENST00000524377","ENST00000445201","ENST00000471622")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='sym-02-02-TP6_vs_Sym-pool-1-1,Sym-pool-1-2,Sym-pool-1-3 pos.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(3013.727509714351,3115.5958672690726,2794.7044730153466,3319.7096349380004),c(2326.899954992616,2349.2373131554814,2281.337953136239,2572.8624199124333),c(2000.9639077037602,2107.718253677259,1912.0283648558602,2366.473728172347),c(2265.4917141990636,2286.9997093668626,2301.636884406969,2697.0454462984176),c(2090.7144134789523,2089.1398644866263,2083.587074305256,2399.70580565592),c(1879.0921682827097,1957.2333012331353,1868.8112853762416,2348.983161075729),c(1919.71608142306,1851.33648284653,1870.1208938453208,2174.077490109554))
targetgene="ENST00000355710"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2346.739540479764,2309.2937763956215,2457.480292227412,2245.788815205686),c(1773.2810456845882,1719.4299195930391,1721.4803326048138,2877.1982873935776),c(1891.37381644142,1874.5594693348207,1912.6831690904,3055.602071779076),c(2028.3614305193453,2243.340494768876,2156.9251485737,1591.641605792192),c(1843.1919659726327,1941.4416704210978,1954.5906401009393,1873.2397360477335),c(2091.659155645007,2007.3949520478432,1963.7578993844947,1656.3567040496766),c(2307.0603695054688,2213.615072063864,2267.5870642109057,2251.0359853346713))
targetgene="ENST00000286307"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2650.946517949363,2640.918023448412,2681.423340439982,2805.486962297446),c(2515.8483882035475,2601.9034061480834,2562.248969753761,2784.498281781505),c(2496.9535448824545,2303.720259638432,2594.9891814807447,2850.9624367486513),c(3111.03595281798,3113.7380283500092,3309.3806013635335,3534.8436102263954),c(1831.8550599799769,2055.698763943488,1929.7080791884316,2114.609561981055),c(2667.0071347722924,2796.976492649725,2689.935795488998,3078.339809004679),c(2825.7238186694744,3010.627968341999,2800.5977111262036,3986.100241319127))
targetgene="ENST00000368748"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2548.9143640154603,2587.9696142551093,2397.8931068843017,2254.534098753995),c(852.1574337812992,906.6253925028667,918.6903410591685,1091.4113868289317),c(1299.0204783251513,1327.425907670693,1327.9429876464674,1703.5812352105438),c(2137.0067796156304,2281.426192609673,2154.305931635541,2705.790729846726),c(2746.365476720883,2636.2734261507535,2605.4660492333796,3099.3284895206198),c(2195.580793911019,2300.9335012598367,2316.697381801381,2964.651122876665),c(3142.212444297784,3062.64745807577,3105.081680187154,3468.379455259249))
targetgene="ENST00000273860"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2205.9729577376206,2383.6073331581515,2313.423360628683,2859.70772029696),c(1243.2806905279265,1372.0140417282112,1331.8718130537054,1518.1812239863984),c(2124.7251314569203,2127.2255623274227,2182.462513720747,2639.3265748795798),c(2575.3671446649905,2659.496412639044,2613.978504282395,3090.583205972311),c(2883.3530907988084,2949.319284012912,2933.5229707377584,3515.603986420116),c(3140.3229599656747,3065.4342164543646,3109.6653098289316,3741.232301966482),c(1845.081450304742,1908.0005698779592,1995.1885026423995,1962.4416282404827))
targetgene="ENST00000414069"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2824.7790765034197,3012.4858072610623,2769.167107868299,2357.7284446240383),c(1685.4200242415054,1698.0647720238117,1749.6369146900201,2167.0812632709076),c(1299.0204783251513,1285.62453199177,1279.487474290531,1264.5680010854448),c(2069.93008582575,2215.4729109829273,2059.359317627288,2861.456777006622),c(2174.796466257817,2246.127253147471,2177.87888407897,2413.698259333214),c(3344.38726783348,3423.0682083740408,3392.5407391500726,4757.434250279958),c(2361.8554151366384,2496.0065877614784,2455.515879523793,3926.632313190627))
targetgene="ENST00000376708"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2346.739540479764,2460.7076482992766,2497.4233505343327,3622.2964457094827),c(2954.2087532529076,3091.4439613212503,2944.6546427249327,3905.643632674686),c(2615.046315639286,2541.523641278528,2509.8646309905866,3398.417186872779),c(2136.062037449576,2092.8555423247526,2035.7863651838595,2039.4001234655998),c(1808.2365058286105,1985.100885019084,2044.2988202328752,2200.3133407544806),c(839.8757856225886,891.7626811503607,819.8149016436771,963.7302470236239),c(1715.6517735552543,1790.0277985174428,1696.5977716923062,1820.768034757881))
targetgene="ENST00000374597"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2386.4187114540596,2309.2937763956215,2368.426916330016,2210.807681012451),c(1758.1651710277138,1790.9567179769742,1714.9322902594172,2098.868051594099),c(3029.78812653728,3055.216102399517,3055.971362596678,3795.453059965996),c(1858.3078406295072,2097.5001396224106,1851.1315710436702,2327.9944805597884),c(3222.51552841243,3276.2989337680438,3242.5905694404864,3921.385143061642),c(2279.6628466898837,2286.9997093668626,2402.476736526079,2467.9190173327283),c(2254.154808206408,2411.4749169441,2279.37354043262,3160.545474358781))
targetgene="ENST00000524377"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2715.1889852410795,2832.275432111927,2697.793446303474,2793.2435653298135),c(2039.698336512001,1970.2381736665782,2036.441169418399,2282.5190061085827),c(2264.5469720330093,2483.930634787567,2392.654673007984,2891.190741070872),c(2111.498741132155,2070.5614752959937,2135.97141306843,2565.8661930737862),c(1880.0369104487643,1834.6159325749607,1897.6226716959873,2375.2190117206555),c(2192.7465674128553,2403.114641808316,2386.761434897127,2126.852958948687),c(2540.4116845209683,2635.3445066912223,2398.547911118841,2709.2888432660498),c(2281.552331021993,2302.7913401789,2183.772122189827,2474.9152441713754),c(2318.3972754981246,2457.9208899206815,2465.337943041888,2744.2699774592847),c(2304.2261430073045,2511.798218573516,2424.740080500428,3288.226614164089))
targetgene="ENST00000445201"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(3079.8594613381765,2951.1771229319747,2984.597701031853,4124.275721382405),c(2575.3671446649905,2752.3883585922067,2619.2169381587128,2630.581291331271),c(1058.1112259812141,1151.860129819216,1098.76150555758,1596.888775921177),c(945.68690822071,996.7305800774344,945.5373146752953,982.9698708299031),c(2178.5754349220356,2342.73487693876,2281.9927573707782,3167.541701197428),c(1733.6018747102928,1785.3832012197845,1766.0070205535121,1694.8359516622352),c(1016.5425706748092,1021.8114054847882,996.6120449693902,1213.8453565052541))
targetgene="ENST00000471622"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=63
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("ENST00000493536","ENST00000196551","ENST00000245458","ENST00000293831","ENST00000355968","ENST00000330752","ENST00000316292","ENST00000392367","ENST00000419869","ENST00000371646")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='sym-02-03-TP6_vs_Sym-pool-1-1,Sym-pool-1-2,Sym-pool-1-3 neg.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1249.1074898277332,1328.453549459583,1286.049671811552,13.434903667897398),c(1481.2020782814598,1495.5099155454702,1471.9314225143485,40.3047110036922),c(1991.9929245241099,2033.1115191729268,2088.1766841663316,211.11991478124483),c(2572.2293956584263,2777.7563850238503,2596.6735411736427,205.3620989235745),c(2138.1942400855123,2266.81271215478,2225.5401473975508,5.757815857670313),c(1832.0852356288258,2019.7825537937338,1935.0605301975872,72.9323341971573),c(2323.687159204042,2300.579424448736,2195.9250888110037,713.9691663511189),c(2185.7096676429683,2156.62659835345,2141.105725044416,101.72141348550888),c(2134.5392071964775,2022.4483468695723,2093.8476528318406,92.12505372272501),c(1790.9661156271814,1764.7550162051716,1664.114249512155,97.88286958039534))
targetgene="ENST00000493536"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1610.0419876199458,1576.3723055125754,1597.3228407850484,105.55995739062242),c(2017.5781547473553,1920.2596122957582,1969.0863421906415,9.596359762783855),c(1552.4752196176435,1616.3592016501548,1527.3808939104367,17.27344757301094),c(1885.0832125198342,1786.969958503827,1824.7916950349113,44.14325490880574),c(2479.026056988032,2552.9411689614594,2487.0348136404677,0.0),c(1686.7976782896822,1698.998787001152,1670.4153258071651,40.3047110036922),c(2439.7344534309054,2486.2963420654937,2549.4154689610673,9.596359762783855))
targetgene="ENST00000196551"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1928.0298489659963,1863.389360011201,2056.671302691281,76.77087810227084),c(2021.2331876363905,2036.6659099407118,2216.718640584537,17.27344757301094),c(2027.6294951922018,1931.811382291059,1841.8046010314386,92.12505372272501),c(2469.8884747654447,2369.001446728594,2389.368131067812,38.38543905113542),c(1258.245072050321,1403.0957555830644,1349.0604347616522,30.70835124090834),c(2860.0632356699375,2858.6187749909554,2865.729498970572,90.20578177016824),c(2172.003294309087,2200.1678852588143,2221.129393991044,57.57815857670314))
targetgene="ENST00000245458"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2442.4757280976814,2517.397261283611,2363.5337182582707,9.596359762783855),c(2256.982808979152,2167.2897706568046,2167.570245483458,32.62762319346511),c(1854.0154329630363,1910.48503768435,1921.1981623485651,564.2659540516908),c(2607.865966326518,2620.474593549371,2623.768169242186,191.9271952556771),c(2318.2046098704895,2376.1102282641637,2399.449853139828,218.7970025914719),c(2807.979017001188,2833.7380396164613,2794.5273368369585,24.950535383238027),c(2717.51695299757,2893.2740849768575,2746.009049365381,24.950535383238027))
targetgene="ENST00000293831"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2236.880128089459,2142.4090352823105,2139.845509785414,500.92997961731726),c(1940.822464077619,1930.0341869071667,2005.6325847017,21.111991478124484),c(2584.1082525477905,2751.98705195741,2640.1509676092123,40.3047110036922),c(1718.7792160687388,1794.9673377313427,1744.7680260882837,63.33597443437345),c(2718.430711219829,2689.7852135211756,2518.540195115518,568.1044979568043),c(2648.071328105904,2743.989672729894,2756.7208790668983,24.950535383238027),c(1907.9271680763034,2008.230783798433,1945.1422522696034,30.70835124090834))
targetgene="ENST00000355968"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2561.2642969913213,2622.2517889332635,2627.5488150191923,67.174518339487),c(2332.8247414266298,2207.276666794384,2222.389609250046,907.8156335593528),c(2363.892520983428,2144.186230666203,2099.5186214973496,1.9192719525567712),c(2712.0344036640176,2735.1036958104323,2619.357415835679,213.0391867338016),c(1968.235210745382,1938.9201638266286,1964.6755887841346,107.47922934317918),c(2494.559946766431,2447.1980436198605,2479.473522086456,74.85160614971407),c(2512.8351112116065,2603.5912374023933,2551.3057918495706,0.0))
targetgene="ENST00000330752"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2364.8062792056867,2400.1023659467114,2299.892847678669,101.72141348550888),c(2219.518721866543,2255.2609421594793,2132.284218231402,72.9323341971573),c(2078.7999556386926,2026.8913353293035,2067.3831323927984,1.9192719525567712),c(2022.1469458586491,2003.787795338702,2028.946566993237,30.70835124090834),c(1879.6006631862817,2008.230783798433,1944.5121446401024,92.12505372272501),c(2308.153269425643,2345.8979067379923,2417.0928667658563,1984.5271989437015),c(2643.50253699461,2582.2648927956843,2497.746643341985,437.5940051829438))
targetgene="ENST00000316292"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2218.604963644284,2079.31859915413,2047.8497958782673,9.596359762783855),c(2329.1697085375945,2374.3330328802713,2255.785313613599,24.950535383238027),c(521.7559449097554,543.8217874710801,575.2882657344178,9.596359762783855),c(2636.19247121654,2629.3605704688334,2595.4133259146406,11.515631715340627),c(710.9038969173199,756.1966358462241,743.5270028111862,38.38543905113542),c(2354.75493876084,2327.237355207122,2276.5788653871323,61.41670248181668),c(2771.4286881108374,2717.331741971508,2843.045624308536,23.031263430681253))
targetgene="ENST00000392367"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2007.526814302509,1917.5938192199196,1979.7981718921587,128.59122082130366),c(2755.8947983324383,2743.1010750379482,2633.2197836847013,38.38543905113542),c(1940.822464077619,2068.6554268507753,2026.426136475233,477.898716186636),c(1705.9866009571163,1649.2373162521644,1675.456186843173,546.9925064786798),c(1225.3497760490052,1138.2936433830941,1260.2152590020105,174.65374768266616),c(1575.3191751741126,1698.110189309206,1619.3766078175836,0.0),c(1599.9906471750992,1574.5951101286828,1637.6497290731127,99.8021415329521))
targetgene="ENST00000419869"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1431.8591342794864,1661.6776839394113,1572.1185356050082,15.35417562045417),c(2416.8904978744363,2451.6410320795912,2405.7509294348383,178.49229158777973),c(1559.7852853957136,1611.0276154984774,1500.9163734713945,487.49507594941986),c(2711.120645441759,2632.0263635446718,2729.626250998355,667.9066394897563),c(2600.555900548448,2560.9385481889753,2517.279979856516,401.1278380843652),c(2448.8720356534927,2478.298962837978,2416.462759136355,635.2790162962913),c(2548.4716818796983,2531.6148243547505,2614.946662429172,1053.6803019536674))
targetgene="ENST00000371646"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=69
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("ENST00000598996","ENST00000368748","ENST00000293525","ENST00000308388","ENST00000284142","ENST00000414790","ENST00000258704","ENST00000524320","ENST00000440444","ENST00000477539")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='sym-02-03-TP6_vs_Sym-pool-1-1,Sym-pool-1-2,Sym-pool-1-3 pos.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2266.12039120174,2324.5715621312834,2295.482094272162,2543.035337137722),c(2415.9767396521775,2453.418227463484,2434.1057727623834,3312.663390112987),c(2140.935514752289,2185.950322187675,2125.983141936392,2535.358249327495),c(1814.7238294059093,1830.5112454091911,1742.8777031997806,2201.4049295826167),c(1979.2003094124873,1902.4876584568342,1889.062673244014,2366.462317502499),c(2492.732430321914,2670.236064298359,2616.836985317675,2201.4049295826167),c(1708.7278756238925,1738.0970854467853,1726.4949048327544,2040.1860855678478),c(2494.559946766431,2538.72360589032,2567.0584825870956,2928.808999601633),c(1821.1201369617206,1821.625268489729,1875.8304130244928,2084.3293404766537),c(1917.064750298891,2076.652806078291,2078.0949620943156,2577.582232283744))
targetgene="ENST00000598996"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2564.0055716580973,2526.283238203073,2580.290742806617,2976.790798415552),c(2433.3381458750937,2488.962135141332,2465.6111542374338,2886.585016645384),c(2415.0629814299186,2203.722276026599,2497.116535712484,2067.0558929036424),c(3009.0058258981167,2978.5794634036934,3184.5639594980807,3792.48137825218),c(1771.7771929597473,1966.466692276961,1856.9271841394625,2326.1576064988067),c(2579.5394614364964,2675.567650450036,2588.4821419901295,3435.4967950766204),c(2733.050842775969,2879.9451195976644,2694.9703313757996,4529.48180803398))
targetgene="ENST00000368748"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2656.2951521062328,2639.1351450802417,2778.1445384699323,3491.1556817007668),c(2778.7387538889075,2642.6895358480265,2769.3230316569184,3312.663390112987),c(1573.491658729595,1643.017132408541,1566.4475669394992,1880.8865135056358),c(2468.060958320927,2466.747192842677,2542.4842850365562,3406.707715788269),c(2130.884174307442,2278.3644821500807,2280.9896187936392,2118.8762356226753),c(1732.4855894026205,1658.1232931716265,1812.8196500743923,2055.540261188302),c(1866.808048074659,1796.7445331152353,1913.636870794553,2276.2565357323306))
targetgene="ENST00000293525"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2228.65630408913,2306.799608292359,2263.976712797112,2704.2541811524907),c(4033.3287930501933,3700.1207892640155,3933.1318233452753,3913.3955112632566),c(2436.07942054187,2400.1023659467114,2451.1186787589104,2999.8220618462333),c(2814.375324556999,2807.080108858075,2909.2069254061416,3347.2102852590087),c(2346.531114760511,2287.2504590695426,2257.675636502102,2734.9625323933988))
targetgene="ENST00000308388"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2084.282504972245,2112.1967137561396,2149.297124227929,1971.092295275804),c(2678.225349440443,2710.222960435938,2689.929470339792,2237.871096681195),c(1772.690951182006,1677.6724423944431,1660.96371136465,2573.74368837863),c(1881.4281796307991,1957.580715357499,1867.6390138409797,3324.1790218283277),c(2332.8247414266298,2388.5505959514107,2439.1466337983916,1932.7068562246686),c(2643.50253699461,2546.720985117836,2604.864940357156,2448.99101146244),c(2321.8596427595244,2439.2006643923446,2474.4326610504477,2498.892082228916))
targetgene="ENST00000284142"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1804.672488961063,1756.7576369776557,1803.368035631877,1834.8239866442732),c(2570.4018792139086,2656.907098919166,2770.5832469159204,3210.941976627478),c(2059.6110329712587,2040.2203007084966,2029.576674622738,2230.194008870968),c(2701.983063219171,2790.1967527110974,2702.531622929812,2951.840263032314),c(1890.565761853387,1829.622647717245,1896.623964798026,2113.118419765005),c(1260.9863467170971,1311.570193312605,1255.1743979660025,1368.4409021729778),c(1912.4959591875972,1736.319890062893,1786.35512963535,1817.5505390712624),c(2429.683112986059,2447.1980436198605,2618.727308206178,2577.582232283744),c(1890.565761853387,1850.948992323954,1892.2132113915188,2228.2747369184112),c(1856.7567076298126,1747.8716600581936,1722.0841514262474,1974.9308391809175))
targetgene="ENST00000414790"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2223.1737547555776,2081.9843922299683,2253.264883095595,3055.4809484703796),c(1780.914775182335,2044.6632891682275,2063.6024866157923,2308.884158925796),c(2729.395809886934,2701.3369835164763,2598.563864062146,1850.1781622647275),c(2280.74052275788,2329.0145505910145,2265.236928056114,2675.465101864139),c(1943.5637387443953,1962.9123015091764,2011.93366099671,2182.2122100570487),c(2781.480028555684,2732.4379027345935,2790.1165834304516,3001.74133379879),c(2347.44487298277,2493.405123601063,2488.925136528971,2804.056322685443))
targetgene="ENST00000258704"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1468.409463169837,1340.8939171468298,1482.0131445863644,1656.3316950564936),c(1715.1241831797038,1884.71570461791,1918.0476242010602,2477.7800907507917),c(1165.9554916021855,1148.0682179945024,1136.0840559903124,1604.5113523374607),c(1891.4795200756457,1979.7956576561544,1923.7185928665692,2372.2201333601693),c(2911.2336961164287,3043.4470949157667,3068.6241556698956,3283.8743108246354),c(1766.2946436261948,1754.0918439018171,1739.0970574227747,1740.7796609689915),c(1801.017456072028,1827.8454523333526,1844.3250315494424,2337.673238214147),c(1845.7916089627074,1713.2163500722916,1764.3013626028148,1777.24582806757),c(2737.619633887263,2715.5545465876157,2698.1208695233045,3138.009642430321),c(2300.843203647573,2216.162643713846,2145.516478450923,1729.2640292536507))
targetgene="ENST00000524320"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2595.0733512148954,2507.6226866722027,2644.5617210157193,3118.816922904753),c(2691.9317227743245,2815.966085777537,2841.785409049534,2299.287799163012),c(2447.0445192089755,2393.882182103088,2376.135870848291,3019.0147813718013),c(2016.6643965250967,1959.3579107413916,1981.0583871511606,2258.98308815932),c(1737.054380513914,1810.9620961863745,1953.3336514531163,2166.8580344365946),c(1959.0976285227944,1963.8008992011225,1901.034718204533,1767.6494683047863),c(1978.2865511902285,1956.6921176655528,1900.404610575032,2262.8216320644333),c(1887.8244871866107,1882.9385092340176,1904.815363981539,2264.74090401699),c(2215.8636889775075,2448.975239003753,2282.2498340526413,2051.7017172831884),c(2585.9357689923077,2501.402502828579,2489.555244158472,2535.358249327495))
targetgene="ENST00000440444"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1767.2084018484534,1695.4443962333673,1713.8927522427343,1401.068525366443),c(2479.939815210291,2427.648894397044,2488.29502889947,2541.116065185165),c(2798.8414347786,2833.7380396164613,2776.8843232109302,3314.582662065544),c(1736.1406222916555,1603.0302362709615,1621.2669307060867,1963.415207465577),c(2477.198540543515,2378.7760213400024,2530.5122400760374,2024.8319099473936),c(2566.746846324874,2528.0604335869652,2662.2047346417476,3289.6321266823056),c(2023.060704080908,2193.0591037232443,2183.322936220983,2105.441331954778),c(3065.65883567816,3166.073576404344,3089.4177074434288,3602.4734549490595),c(2411.4079485408834,2571.6017204923296,2508.458473043502,1905.8370488888738),c(2333.7384996488886,2349.4522975057776,2378.656301366295,3013.2569655141306))
targetgene="ENST00000477539"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-02-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=75
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("ENST00000493536","ENST00000293831","ENST00000451311","ENST00000245458","ENST00000392367","ENST00000293308","ENST00000330752","ENST00000355968","ENST00000225792","ENST00000196551")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='sym-03-01-TP6_vs_Sym-pool-1-1,Sym-pool-1-2,Sym-pool-1-3 neg.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1285.5216893161644,1375.5807806235184,1329.8406806599514,173.61097422504548),c(1524.3823397084875,1548.563514240389,1522.0518520439227,207.27024473806452),c(2050.063849823876,2105.236672954254,2159.28075242875,747.5901134996857),c(2647.2154758046836,2876.29800684222,2685.0923297401564,662.5561669404798),c(2200.5272516458117,2347.2284758331743,2301.321550265041,46.06005438623656),c(1885.4945040811335,2091.4348591018443,2000.9508722717837,510.20367935523575),c(2391.427692707393,2382.193070925946,2270.698075502171,2194.938745559504),c(2249.4278572379408,2233.133481319919,2214.0120690262197,471.2297871822663),c(2196.765666600263,2094.1952218723263,2165.1448220641932,550.9491120815219),c(1843.176672318714,1827.3601540590685,1720.7786563561644,614.7245720009264))
targetgene="ENST00000493536"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2513.6792066877156,2606.70257625848,2444.0139113941586,37.2023516196526),c(2322.7787656261344,2244.174932401847,2241.377727324955,86.80548711252274),c(1908.0640143544238,1978.2599855120832,1986.6164798295893,685.586194133598),c(2683.8909299987804,2713.436603383783,2713.1095513317186,437.5705166692473),c(2385.78531513907,2460.403349422935,2481.1530190852986,648.3838425139454),c(2889.837711242555,2934.2656250223413,2889.6832036878413,58.46083825945409),c(2796.7384813652325,2995.9137268964387,2839.512830140161,69.09008157935484))
targetgene="ENST00000293831"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2365.096597388554,2280.059648418113,2255.0605564743223,306.4765157238048),c(2193.944477816102,2298.4620668879925,2168.402638528328,604.0953286810256),c(2457.2554310044898,2543.2142325373943,2541.748405318212,102.74935209237387),c(2408.3548254123607,2486.166735280767,2483.759272256607,56.6892977061373),c(2494.8712814599735,2454.882623881971,2423.163886023694,343.6788673434574),c(2578.5665487234255,2560.6965300837805,2559.3406142245412,648.3838425139454),c(2568.222189848167,2651.7885015096854,2627.754759971379,1156.8159813158643))
targetgene="ENST00000451311"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1984.2361115267788,1929.4935765669018,2126.702587787399,63.77545991940447),c(2080.156530188263,2108.91715664823,2292.1996641654623,170.06789311841192),c(2086.7393040179727,2000.342887675939,1904.5195049333845,393.28200283632754),c(2541.8910945293287,2453.0423820349833,2470.7280064000665,173.61097422504548),c(1294.9256519300352,1452.8709381970136,1394.9970099426534,410.9974083694955),c(2943.4402981416197,2960.029010880173,2963.309855777295,513.7467604618693),c(2235.321913317134,2278.2194065711246,2296.7606072152516,224.98565027123243))
targetgene="ENST00000245458"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2283.2821226478763,2153.0829609759417,2117.5807016878207,166.52481201177832),c(2397.0700702757154,2458.563107575947,2332.5965883207377,90.34856821915632),c(536.9662652520335,563.1140051783233,594.8772863510709,164.75327145846154),c(2713.0432141017804,2722.637812618723,2683.7892031545025,357.8511917699917),c(731.628291359163,783.0229058933874,768.8446855358857,579.2937609345905),c(2423.401165594554,2409.7966986307656,2354.0981769840296,239.15797469776675),c(2852.2218607870714,2813.729784044628,2939.853577235522,304.704975170488))
targetgene="ENST00000392367"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1737.852291043359,1797.9162845072608,1801.572504666715,168.2963525650951),c(2155.388231099231,2122.71897050064,2083.0478471679885,138.1801631587097),c(2082.977718972424,2119.038486806664,2035.483726791616,614.7245720009264),c(1942.8586760257465,2018.7453061458189,1957.2961316523733,157.66710924519438),c(2993.281299995136,3117.3696887976457,2938.5504506498683,783.0209245660214),c(1957.90501620794,2057.3903849325666,1901.9132517620765,2249.856502712324),c(2665.0830047710383,2698.714668607879,2531.3233926329794,379.1096784097932))
targetgene="ENST00000293308"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2635.9307206680382,2715.276845230771,2717.018931088681,246.2441369110339),c(2400.8316553212635,2285.5803739590765,2298.0637338009055,1854.80295932268),c(2432.8051282084252,2220.251788391003,2171.0088916996365,38.97389217296939),c(2791.0961037969096,2832.132202514508,2708.5486082819293,969.0326826642845),c(2025.6135470278111,2007.703855063891,2031.5743470346538,914.114925511464),c(2567.28179358678,2534.0130233024547,2563.9015572743306,180.69713643831264),c(2586.0897188145223,2695.954305837397,2638.179772656611,21.25848663980149))
targetgene="ENST00000330752"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2302.0900478756184,2218.4115465440154,2212.708942440566,1482.7794431261539),c(1997.4016591861982,1998.502645828951,2073.9259610684103,115.1501359655914),c(2659.440627202716,2849.614500060894,2730.050196945221,76.176243792622),c(1768.885367669133,1858.6442654578643,1804.178757838023,267.5026235508354),c(2797.6788776266194,2785.2060354163145,2604.298481429606,1922.121500348718),c(2725.2683654998127,2841.333411749448,2850.5894061182203,251.55875857098428),c(1963.5473937762627,2079.4732870964226,2011.375884957016,389.73892172969397))
targetgene="ENST00000355968"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1911.8255993999721,1981.0203482825652,1857.606947849839,1309.1684689011083),c(2702.6988552265225,2682.152491984987,2449.2264177367747,184.24021754494623),c(2479.82494127778,2603.022092564504,2564.5531205671577,125.77937928549214),c(1736.9118947819718,1645.1762112072581,1712.3083335494132,393.28200283632754),c(1884.5541078197464,1987.4611947470232,1868.6835238278984,1023.950439817105),c(1439.7466761836486,1492.4361379072554,1401.5126428709239,903.4856821915632),c(1860.1038050236818,1860.4845073048523,1780.0709160034235,380.88121896311003))
targetgene="ENST00000225792"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1656.9782125640684,1632.2945182783421,1651.7129473165,763.5339784795368),c(2076.3949451427147,1988.381315670517,2036.1352900844429,230.3002719311828),c(1597.733248096681,1673.699959835572,1579.3894218127007,295.84727240390407),c(1940.0374872415853,1850.3631771464184,1886.927296027055,1341.0561988608106),c(2551.2950571431993,2643.5074131982396,2571.720316788255,784.7924651193383),c(1735.9714985205846,1759.271205720513,1727.2942892844346,1181.6175490622993),c(2510.8580179035544,2574.4983439361904,2636.22508277813,86.80548711252274))
targetgene="ENST00000196551"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=81
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("ENST00000456354","ENST00000321253","ENST00000292807","ENST00000284142","ENST00000414069","ENST00000317450","ENST00000563280","ENST00000314133","ENST00000334269","ENST00000553704")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='sym-03-01-TP6_vs_Sym-pool-1-1,Sym-pool-1-2,Sym-pool-1-3 pos.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2554.116245927361,2729.078659083181,2640.7860258279193,3000.9896973186433),c(1834.7131059662302,2037.1477246156987,1930.5820366464654,2301.2311787585113),c(3148.4466831240074,3170.7367023602974,3222.6320463224497,3932.8200283632755),c(2786.394122489974,2966.469857344631,2952.2332797992353,3573.197296039967),c(1759.4814050552623,1752.830359256055,1705.792700621143,2494.329099070041),c(1879.8521265128109,1966.2984135066613,2008.118068492881,2469.527531323606),c(2768.5265935236193,2912.1827228584852,2733.308013409356,3309.237753595765))
targetgene="ENST00000456354"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2539.069905745167,2740.1201101651086,2653.1657283916325,3362.383970195269),c(2415.8779955034574,2463.1637121934173,2354.0981769840296,3218.8891853766086),c(1630.6471172452295,1718.7858850867774,1763.130270389921,2301.2311787585113),c(1853.521031193972,1961.6978088891915,1985.9649165367623,1879.604527069115),c(2324.6595581489087,2278.2194065711246,2340.415347834662,2830.9218042002317),c(2312.4344067508764,2339.867508445222,2264.1824425739005,2763.6032631741937),c(2483.5865263233286,2423.5985124831755,2555.431234467579,3008.0758595319107))
targetgene="ENST00000321253"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1804.620425601843,1984.7008319765412,1966.4180177519515,2786.633290367312),c(2370.7389749568765,2569.89773931872,2488.971778599223,3068.3082383446817),c(1774.5277452374557,1871.52595838678,1767.0396501468829,2253.3995838189576),c(3667.545419409686,3901.312715614527,3827.9343453587526,5911.630826418131),c(1923.1103545366175,1900.969827938588,1851.0913149215687,1884.9191487290652),c(2345.348275899425,2483.406372510285,2475.9405127426826,3227.746888143193),c(2035.0175096416822,2177.0061049867854,2033.5290369131349,132.8655414987593))
targetgene="ENST00000292807"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2145.043872223973,2187.1274351452194,2222.482391832971,2129.3917450867825),c(2756.301442125587,2806.3688166566762,2781.523697078556,2797.2625336872125),c(1824.368747090972,1737.1883035566573,1717.5208398920292,3050.5928328115137),c(1936.275902196037,2027.0263944572648,1931.2335999392924,3287.979266955964),c(2400.8316553212635,2473.285042351851,2522.2015065334012,1934.5222842219355),c(2720.566384192877,2637.0665667337817,2693.562652546908,2790.1763714739454),c(2389.5469001846186,2525.7319349910085,2558.6890509317145,2749.430938747659))
targetgene="ENST00000284142"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2195.825270338876,2361.030289685584,2301.9731135578677,2600.6215322690487),c(1237.5614799854222,1359.0186040006265,1325.2797376101623,1496.9517675526881),c(2114.9511918595854,2107.0769148012423,2171.6604549924637,2333.1189087182133),c(2563.520208541232,2634.3062039632996,2601.040664965471,2972.6450484655747),c(2870.089389753426,2921.3839320934253,2919.0035518650575,3247.2338342296775),c(3125.877172850717,3036.3990475301744,3094.2740776355263,3675.9466481323407),c(1836.5938984890042,1889.92837685666,1985.3133532439354,2180.7664211329693))
targetgene="ENST00000414069"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2977.2945635515553,3063.0825543115,2953.5364063848897,3362.383970195269),c(2490.169300153038,2588.3001577886002,2458.99986712918,3227.746888143193),c(2541.8910945293287,2593.820883329564,2425.770139195002,3038.192048938296),c(2357.5734272974573,2294.7815831940165,2401.010734067575,2586.4492078425146),c(2898.301277595039,3008.7954198253547,3047.361520551981,3438.5602139878906),c(1684.2497041442944,1571.5665373277386,1612.6191497468787,2058.530122954111),c(2765.705404739458,2874.4577649952316,2796.5096528135773,3192.316077076857))
targetgene="ENST00000317450"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2133.7591170873275,2177.9262259102798,2167.0995119426743,2186.0810427929196),c(2260.712612374586,2280.9797693416067,2354.0981769840296,2868.124155819884),c(2303.0304441370054,2269.0181973361846,2264.8340058667277,2613.0223161422664),c(1746.3158573958428,1818.1589448241286,1782.0256058819045,2271.1149893521256),c(2632.16913562249,2581.8593113241423,2780.220570492902,4055.056326542134),c(1587.3888892214231,1638.7353647428001,1621.7410358464572,2134.706366746733),c(2502.3944515510702,2686.7530966024574,2784.129950249864,2859.2664530533),c(1911.8255993999721,2039.9080873861808,1937.0976695747356,1806.9713643831265),c(2368.858182434102,2307.6632761229325,2333.248151613565,2947.8434807191397),c(1708.700006940359,1761.111447567501,1808.0881375949853,1844.1737160027792))
targetgene="ENST00000563280"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2083.9181152338115,2106.156793877748,2024.4071508135567,2147.1071506199505),c(2031.2559245961338,1927.6533347199138,1918.202334082752,2919.498831866071),c(2155.388231099231,2236.813965013895,2096.079113024529,3810.5837301844167),c(1422.8195434786808,1485.0751705193034,1467.97209873928,1806.9713643831265),c(1749.137046180004,2069.3519569379887,1937.0976695747356,1984.1254197148055),c(1604.316021926391,1463.9123892789416,1513.5815292371715,1640.4465523713482),c(2461.017016050038,2404.2759730898015,2582.145329473487,1877.8329865157982))
targetgene="ENST00000314133"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2067.9313787902306,2118.11836588317,2163.190132185712,2313.6319626317286),c(1924.9911470593916,2075.792803402446,2083.0478471679885,2593.5353700557816),c(1579.8657191303262,1588.1287139506305,1555.933143270928,1856.5744998759967),c(1826.2495396137463,2152.162840052448,1929.9304733536385,2168.365637259752),c(2382.023730093522,2563.456892854262,2464.212373471796,2760.06018206756),c(2554.116245927361,2722.637812618723,2718.3220576743347,3535.9949444203144),c(2104.6068329843274,2123.639091424134,2024.4071508135567,2611.2507755889496))
targetgene="ENST00000334269"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1563.8789826867455,1570.6464164042447,1520.0971621654417,1714.8512556106534),c(1911.8255993999721,1967.2185344301554,1935.1429796962545,2269.343448798809),c(2980.1157523357165,3072.28376354644,2960.703602605987,3282.664645296013),c(2057.5870199149726,1982.8605901295532,2041.3477964270592,2217.968772752622),c(2217.4543843507795,2184.3670723747377,2197.071423412717,2154.1933128332175),c(1934.3951096732626,1995.7422830584692,1962.5086379949894,2152.421772279901),c(2579.5069449848124,2605.782455334986,2572.371880081082,2788.4048309206287),c(2425.281958117328,2478.805767892815,2453.135797493737,2776.004047047411),c(1108.7271921753897,1185.115749460262,1181.2842498953905,1328.655414987593),c(1395.5480518984548,1297.3705021265291,1327.2344274886434,1541.240281385608))
targetgene="ENST00000553704"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-01-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=87
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("ENST00000493536","ENST00000293831","ENST00000451311","ENST00000245458","ENST00000392367","ENST00000196551","ENST00000330752","ENST00000443387","ENST00000225792","ENST00000534624")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='sym-03-02-TP6_vs_Sym-pool-1-1,Sym-pool-1-2,Sym-pool-1-3 neg.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1303.4047158157284,1395.297829565009,1352.5080208967324,183.53878416885522),c(1545.5881816659075,1570.7600315437523,1547.9954614477053,180.1399177953579),c(2078.5825021787036,2135.412330464709,2196.0860270709313,657.6806432717312),c(2684.0411668041515,2917.525762689109,2730.8601441035935,844.6182938140837),c(2231.1390161000763,2380.8727513179515,2340.5479322916503,33.98866373497319),c(1911.723815077202,2121.4126866898096,2035.05738960013,407.86396481967824),c(2424.695093137818,2416.3385155476976,2309.4024756614954,1903.3651691584985),c(2280.7198831245223,2265.1423627787804,2251.7502474312087,560.8129516270576),c(2227.325103252042,2124.212615444789,2202.0500506809612,382.3724670184484),c(1868.8172955368161,1853.5528357967275,1750.1095948987115,662.7789428319771))
targetgene="ENST00000493536"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2548.647260698933,2644.066054286067,2485.672506802373,110.46315713866287),c(2355.0911836611917,2276.3420777987,2279.582357611347,241.31951251830964),c(1934.6072921654081,2006.6156077356316,2020.4786652200573,701.8659061271964),c(2721.2268170724865,2752.3299661452916,2759.354923573735,365.3781351509618),c(2418.9742238657664,2495.6698302721297,2523.444656332561,814.0284964526079),c(2930.038545502365,2976.324266543688,2938.938301164629,32.28923054822453),c(2835.6442025135157,3038.85600873824,2887.912765834375,78.17392659043833))
targetgene="ENST00000293831"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2397.9977032015777,2312.7411516134393,2293.4984127014163,278.70704262678015),c(2224.4646686160163,2331.407343313306,2205.363397130978,696.7676065669503),c(2491.4385679784186,2579.667692921528,2585.0729003028678,47.58412922896246),c(2441.8577009539727,2521.802498651942,2526.095333492574,44.18526285546515),c(2529.5776964587617,2490.0699727621695,2464.467089522267,572.7089839342982),c(2614.4372573275255,2597.4005750364013,2602.9649711329566,630.4897122837526),c(2603.948996995431,2689.7982239507396,2672.545246583303,1157.314000175837))
targetgene="ENST00000451311"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2011.839027338103,1957.1501997309856,2162.9525625707665,112.16259032541151),c(2109.093804962978,2139.1455688046826,2331.2705622316043,164.84501911461996),c(2115.768152447038,2029.015037775471,1936.9823346796416,401.0662320726836),c(2577.2516070591905,2488.203353592183,2512.8419476925083,122.35918944590348),c(1312.9394979358142,1473.6958347044474,1418.7749498970622,582.9055830547902),c(2984.386803586854,3002.4569349235007,3013.8199309350016,559.1135184403089),c(2266.4177099443937,2310.8745324434526,2335.9092472616276,226.0246138375717))
targetgene="ENST00000245458"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2315.0450987568315,2183.9444288843615,2153.67519251072,112.16259032541151),c(2430.4159624098693,2493.803211102143,2372.356058211809,79.87335977718699),c(544.436059056899,571.1854660159099,605.0170617730116,125.7580558194008),c(2750.7846416447524,2761.663061995225,2729.534805523587,370.47643471120773),c(741.8060489426749,794.2464568293127,781.9497622038923,633.8885786572499),c(2457.11335234611,2444.3378030974973,2394.224144781918,343.2855037232292),c(2891.899417022022,2854.060710909563,2989.9638364948833,390.86963295219164))
targetgene="ENST00000392367"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1680.0286095591173,1655.6912037781442,1679.8666501583618,812.3290632658592),c(2105.279892114944,2016.882013170558,2070.841531260308,248.11724526530426),c(1619.9594822025767,1697.6901351028434,1606.3103589679956,321.1928722954966),c(1967.0255513736997,1876.8855754215604,1919.0902638495527,1198.1003966578048),c(2586.7863891792767,2681.3984376857998,2615.555687643019,1108.030437760126),c(1760.120779367838,1784.487926507222,1756.7362877987443,1028.157077982939),c(2545.7868260629075,2611.400218811301,2681.159947353346,30.58979736147587))
targetgene="ENST00000196551"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2672.5994282600486,2754.1965853152783,2763.330939313755,234.521779771315),c(2434.2298752579036,2318.341009123399,2337.234585841634,2042.7186904718885),c(2466.6481344661956,2252.0760285888737,2208.0140742909907,57.78072834945442),c(2829.923333241464,2872.7269026094295,2754.7162385437123,990.7695478744685),c(2053.7920686664806,2036.4815144554175,2066.202846230285,800.4330309586186),c(2602.9955187834225,2570.334597071595,2607.60365616298,122.35918944590348),c(2622.065083023594,2734.5970840304185,2683.147955223356,6.797732746994638))
targetgene="ENST00000330752"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2049.0246776064378,2157.8117605045486,2075.480216290331,977.1740823804791),c(2467.601612678204,2397.672323847831,2499.588561892442,1721.525818176392),c(3571.7293821841395,3482.1780616100655,3532.6899850075847,180.1399177953579),c(2666.878558987997,2715.9308923305525,2698.389348893432,1194.7015302843076),c(1856.4220787807046,2089.6801608000364,1958.1877519597472,243.0189457050583),c(2072.861632906652,1998.2158214706917,2082.7695784803673,1373.142014892917),c(2904.2946337781336,2853.1274013245697,3058.2187733652227,480.9395918498706))
targetgene="ENST00000443387"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1938.4212050134424,2009.4155364906114,1889.2701457994042,1320.4595861037083),c(2740.296381312658,2720.597440255519,2490.9738611223993,244.71837889180696),c(2514.3220450666245,2640.332815946094,2608.266325452983,134.2552217531441),c(1761.0742575798465,1668.7575379680507,1741.4948941286686,310.99627317500466),c(1910.7703368651935,2015.9487035855645,1900.5355237294602,949.9831513925006),c(1459.7751425851354,1513.8281468591601,1425.4016427970953,785.1381322778807),c(1885.9799033529705,1887.151980856487,1810.4125002890116,343.2855037232292))
targetgene="ENST00000225792"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1767.7486050639066,1637.025012078278,1693.782705248431,147.85068724713338),c(2370.346835053329,2508.736164462036,2304.101121341469,234.521779771315),c(2851.8533321176615,2862.460497174503,2905.8048366644643,304.19854042801),c(2884.2715913259535,3036.05607998326,2983.9998128848533,56.08129516270576),c(3310.4763520937886,3285.2497391764755,3129.787056685579,428.25716306066215),c(1268.1260219714109,1259.967939740978,1423.4136349270852,1436.021042802617),c(2505.7407411585473,2604.8670517163478,2848.8152777241803,3543.3181943709546))
targetgene="ENST00000534624"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=93
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("ENST00000314133","ENST00000292807","ENST00000563027","ENST00000334269","ENST00000524320","ENST00000284142","ENST00000308388","ENST00000306704","ENST00000286307","ENST00000221476")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='sym-03-02-TP6_vs_Sym-pool-1-1,Sym-pool-1-2,Sym-pool-1-3 pos.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2112.9077178110124,2136.3456400497025,2058.9134840402485,2289.136502550444),c(2059.512937938532,1955.283580560999,1950.898389769711,3628.289853708388),c(2185.3720619236647,2268.8756011187534,2131.8071059406116,3823.724670184484),c(1442.612534768981,1506.3616701792137,1492.9939103774316,1425.8244436821253),c(1773.469474335958,2099.0132566499697,1970.1157991798066,1906.764035531996),c(1626.6338296866368,1484.8955497243674,1539.3807606776625,1594.0683291702426),c(2495.252480826453,2438.737945587537,2626.1583962830723,1923.7583673994825))
targetgene="ENST00000314133"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1829.7246888444643,2013.1487748305847,1999.935917229955,3227.223621635704),c(2403.718572473629,2606.733670886334,2531.3966878126002,3352.981677455105),c(1799.2133860601898,1898.3516958764067,1797.1591144889458,2487.9701854000373),c(3718.5650268334607,3957.232640371664,3893.182078769379,5699.898908355003),c(1949.8629435575453,1928.2176025961928,1882.6434528993711,2185.471078158776),c(2377.9746607493976,2519.0025698969625,2518.1433020125346,3155.8474277922605),c(2063.3268507865664,2208.2104780941877,2068.190854100295,227.72404702432036))
targetgene="ENST00000292807"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2296.929012728668,2418.2051347176844,2331.2705622316043,2917.926781647448),c(2790.8307265491126,2914.725833934129,2852.7912934642004,3409.062972617811),c(1950.816421769554,1935.6840792761393,2021.8040038000638,2280.639336616701),c(1692.4238263152288,1625.825297058358,1662.637248618276,1988.3368284959315),c(2277.8594484884966,2107.4130429149095,2259.039609621245,2532.1554482555025),c(1557.0299202100105,1556.7603877688528,1737.5188783886488,1889.7697036645093),c(2065.2338072105836,2019.6819419255378,2123.1924051705687,2406.3973924361017))
targetgene="ENST00000563027"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2096.6985882068666,2148.478664654616,2200.062042810951,2151.4824144238028),c(1951.7698999815625,2105.5464237449232,2118.5537201405455,2530.456015068754),c(1601.8433961744138,1610.892343698465,1582.4542645278768,2120.8926170623267),c(1851.6546877206617,2183.011119299368,1962.8264369897704,2725.89083154485),c(2415.160311017732,2600.200503791381,2506.215254792475,3001.1990077981327),c(2589.6468238153025,2761.663061995225,2764.6562778937614,3201.7321238344744),c(2133.884238475201,2154.0785221645756,2058.9134840402485,2516.8605495747647))
targetgene="ENST00000334269"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1532.2394866977875,1408.3641637549154,1558.5981700877583,1386.737480386906),c(1789.678603940104,1979.5496297708253,2017.1653187700408,2418.2934247433423),c(1216.6381985229475,1205.8359838113654,1194.7927298759473,1475.1080060978363),c(1973.6998988577598,2079.4137553651103,2023.1293423800705,2448.8832221048183),c(3037.781583459335,3196.5853286021106,3227.199442316064,3689.4694484313395),c(1843.0733838125843,1842.3531207768076,1828.967240409104,2311.229133978177),c(1879.3055558689105,1919.817816331253,1939.633011839655,2414.894558369845),c(1926.025988257331,1799.420879867115,1855.474012009236,2455.6809548518127),c(2856.6207231777043,2852.194091739577,2837.5498997941245,2797.2670253882934),c(2400.8581378376034,2327.6741049733323,2256.3889324612314,2600.132775725449))
targetgene="ENST00000524320"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2174.8838015915703,2218.4768835291143,2260.3649482012515,2029.1232249778993),c(2794.644639397147,2846.5942342296166,2828.9351990240816,2596.7339093519518),c(1849.7477312966446,1762.0884964673824,1746.796248448695,2928.12338076794),c(1963.2116385256654,2056.081015740277,1964.1517755697769,3045.3842706535975),c(2434.2298752579036,2508.736164462036,2565.1928216027686,2125.990916622573),c(2758.412467340821,2674.865270590846,2739.4748448736364,2486.2707522132887),c(2422.7881367138007,2561.934810806655,2602.3023018429535,2722.491965171352))
targetgene="ENST00000284142"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2325.533359088926,2422.8716826426507,2380.9707589818518,2962.1120445029132),c(4208.65282780587,3886.301111912172,4136.38170820059,4420.225718733263),c(2541.972913214873,2520.8691890669493,2577.783538112831,3089.5695335090627),c(2936.7128929864252,2948.324978993888,3059.5441119452294,3572.208558545682),c(2448.5320484380327,2402.338871772798,2374.3440660818187,3344.4845115213616))
targetgene="ENST00000308388"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2289.3011870325995,2207.277168509195,2270.9676568413042,2489.669618586786),c(1792.5390385761298,1790.087784017182,1671.9146186783223,2012.1288931104127),c(1630.447742534671,1633.2917737383048,1579.8035873678637,2063.1118887128723),c(2179.651192651613,2147.5453550696225,2340.5479322916503,2554.248079683235),c(2310.2777076967886,2282.8752448936534,2325.969207911578,2552.5486464964865),c(2034.722504426309,2120.479377104816,1977.405161369843,2386.0041941951176),c(3054.944191275489,3027.6562937183203,3183.463269175846,3725.1575453530613))
targetgene="ENST00000306704"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2368.439878629312,2320.207628293386,2486.9978453823796,1840.486141248798),c(1789.678603940104,1727.5560418226296,1742.157563418672,3024.991072412614),c(1908.8633804411763,1883.4187425165137,1935.656996099635,2904.331316153459),c(2047.1177211824204,2253.9426477588604,2182.832641270865,2369.0098623276313),c(1860.235991628739,1950.6170326360325,1978.067830659846,2193.9682440925194),c(2111.0007613869952,2016.882013170558,1987.3452007198923,1949.2498652007123),c(2328.3937937249516,2224.0767410390745,2294.823751281423,2438.6866229843263))
targetgene="ENST00000286307"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2013.7459837621202,2277.275387383693,2186.8086570108853,2603.531642098946),c(2981.5263689508283,3175.119208147264,3044.3027182751534,3867.909933039949),c(1154.66211474239,1177.836696261566,1152.381895315736,1303.4652542362217),c(1592.308614054328,1625.825297058358,1610.2863747080155,1779.3065465258464),c(1882.1659905049362,1791.9544031871685,1809.087161709005,2195.667677279268),c(2319.8124898168744,2531.1355945018754,2500.2512311824453,2764.977794840069),c(1611.3781782944995,1599.6926286785454,1649.38386281821,1816.6940766343168))
targetgene="ENST00000221476"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-02-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=99
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("ENST00000493536","ENST00000293831","ENST00000451311","ENST00000245458","ENST00000392367","ENST00000443387","ENST00000330752","ENST00000196551","ENST00000355968","ENST00000293308")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='sym-03-03-TP6_vs_Sym-pool-1-1,Sym-pool-1-2,Sym-pool-1-3 neg.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1268.4737710107154,1367.7297876602763,1318.440063840436,186.36518607162532),c(1504.166776011975,1539.725239218893,1509.003424366124,126.07056704845242),c(2022.8769720580538,2093.2212402452924,2140.769412820777,551.7871195454004),c(2612.1094845612024,2859.8818168735947,2662.0732499198616,597.4648612296223),c(2171.3450067045164,2333.831898542719,2281.5925063617933,14.616877338951005),c(1860.4900591634853,2079.4981988975305,1983.7968819470748,471.3942941811699),c(2359.713825662216,2368.596936623716,2251.231564176345,1929.4278087415325),c(2219.5971179646167,2220.388090067887,2195.0315222585996,475.04851351590764),c(2167.6333058383548,2082.242807167083,2146.5832102605436,528.034693869605),c(1818.7334244191675,1816.930674443685,1706.0265598248857,548.1329002106627))
targetgene="ENST00000493536"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2480.344103812467,2591.825075880644,2423.0615773961176,23.752425675795383),c(2291.9752848547673,2231.3665231460964,2222.1625769775114,140.68744438740342),c(1882.7602643604546,1966.9692598458823,1969.5853770943115,615.735957903311),c(2648.2985680062775,2697.9499289700034,2689.8502821320803,369.07615280851286),c(2354.1462743629736,2446.3608375943672,2459.8822945146403,820.3722406486252),c(2851.5141904286234,2917.518590534195,2864.9101828183902,76.73860602949277),c(2759.6495939911247,2978.8148418875317,2815.169915833719,49.33196101895964))
targetgene="ENST00000293831"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2333.731919599085,2267.0464306502777,2235.728104336967,367.249043141144),c(2164.8495301887338,2285.34381911396,2149.8130977270807,654.1052609180574),c(2424.6685908200434,2528.699085680939,2519.95820139223,40.196412682115266),c(2376.416479559943,2471.977181443523,2462.46620448787,56.640399688435146),c(2461.785599481659,2440.8716210552625,2402.3902976102804,507.9364875285474),c(2544.3709437537536,2546.0816047214375,2537.39959371153,756.4234022907145),c(2534.1637663718093,2636.6536776166663,2605.227230508809,1300.9020831666394))
targetgene="ENST00000451311"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1957.9222069002262,1918.4811804171234,2108.470538155406,113.28079937687029),c(2052.5705789873464,2096.880717938029,2272.54882145549,116.93501871160804),c(2059.066055503129,1988.9261260023015,1888.192212937577,491.49250052222754),c(2508.1818603086786,2439.041882208894,2449.5466546217217,84.04704469896828),c(1277.7530231761193,1444.5788192077434,1383.0378131711777,526.2075842022361),c(2904.405927771426,2943.1349343833504,2937.9056395621283,540.8244615411871),c(2205.6782397165107,2265.2166918039093,2277.0706639086416,241.17847609269157))
targetgene="ENST00000245458"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2253.0024257600708,2140.794450250867,2099.4268532491024,147.9958830568789),c(2365.281376961458,2444.531098747999,2312.5994260405496,115.10790904423916),c(529.8452986445636,559.9000869886884,589.7774513896709,151.65010239161668),c(2677.0642497190297,2707.0986232018445,2660.7812949332465,268.5851211032247),c(721.9258184684247,778.5538791296958,762.2534421027509,610.2546289012045),c(2391.2632830245893,2396.04301931924,2333.916683319694,292.3375467790201),c(2814.3971817670076,2797.6706960970737,2914.6504498030613,286.85621777691347))
targetgene="ENST00000392367"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1994.1112903453015,2115.178106401712,2023.201509038827,1264.359889819262),c(2401.4704604065337,2350.2995481600333,2436.6271047555733,1436.1081985519363),c(3476.007861160307,3413.377817899994,3443.706016821835,224.7344890863717),c(2595.4068306634754,2662.270021465822,2630.420352747798,988.4663300465617),c(1806.6703966041425,2048.39263850927,1908.8634927234145,266.75801143585585),c(2017.3094207588115,1958.7354350372252,2030.3072614652087,1330.1358378445414),c(2826.460209582033,2796.7558266738897,2981.186131613725,480.5298425180143))
targetgene="ENST00000443387"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2600.974381962718,2699.779667816372,2693.7261470919248,160.78565072846106),c(2368.99307782762,2272.5356471893824,2278.362618895256,1819.8012286994),c(2400.542535189993,2207.579918143309,2152.3970077003105,67.60305769264839),c(2754.082042691882,2815.9680845607563,2685.3284396789286,1048.7609490697346),c(1998.7509164280036,1996.2450813877745,2014.1578241325233,824.0264599833629),c(2533.235841155269,2519.550391449098,2541.9214361646823,195.5007344084697),c(2551.794345486077,2680.5674099295047,2615.5628704017277,7.308438669475502))
targetgene="ENST00000330752"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1635.0042315441701,1622.978356728649,1637.5529455342996,778.348718299141),c(2048.8588781211847,1977.0328235009079,2018.6796665856752,191.84651507373195),c(1576.5449429021253,1664.1474807719349,1565.8494437771765,336.1881787958731),c(1914.3097217228278,1839.802410023288,1870.750820618277,1251.5701221476797),c(2517.4611124740823,2628.419852808009,2549.673166084371,873.3584210023225),c(1712.949949733563,1749.230337128059,1712.4863347579599,1098.0929100886942),c(2477.5603281628455,2559.8046460691994,2613.624937921805,38.369303014746386))
targetgene="ENST00000196551"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2271.5609300908786,2205.7501792969406,2193.739567271985,1635.2631522951438),c(1970.9131599317918,1987.0963871559334,2056.1463611975055,109.62658004213253),c(2624.1725123762276,2833.3506036012545,2706.645696958073,71.25727702738615),c(1745.4273323124767,1848.0362348319452,1788.711678968235,308.78153378534),c(2760.577519207665,2769.3097439783655,2581.972040749742,1741.2355130025385),c(2689.127277534055,2825.1167787925974,2826.151533219945,153.47721205898554),c(1937.5078521363378,2067.6048963961366,1994.1325218399934,319.7441917895532))
targetgene="ENST00000355968"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1714.8058001666436,1787.6548529017925,1786.1277689950055,215.59894074952732),c(2126.8045963105774,2110.603759285791,2065.190046103809,133.37900571792792),c(2055.3543546369674,2106.9442815930543,2018.0336890923677,814.8909116465185),c(1917.093497372449,2007.2235144659842,1940.516389895478,175.40252806741205),c(2953.5859642480664,3099.577605747837,2913.358494816446,915.3819433518066),c(1931.9403008370953,2045.6480302397176,1885.6083029643476,2006.1664147710253),c(2629.7400636754696,2683.312018199057,2509.6225614993114,422.0623331622103))
targetgene="ENST00000293308"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=105
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("ENST00000284142","ENST00000376708","ENST00000215375","ENST00000400385","ENST00000334269","ENST00000306704","ENST00000292807","ENST00000562866","ENST00000203166","ENST00000229379")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='sym-03-03-TP6_vs_Sym-pool-1-1,Sym-pool-1-2,Sym-pool-1-3 pos.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2116.5974189286335,2174.64461890868,2203.429229671596,1965.9700020889102),c(2719.748809679888,2790.3517407116005,2757.6779189293593,2634.6921403459187),c(1800.1749200883598,1727.2734709716399,1702.7966723583488,3020.2122801607516),c(1910.5980208566664,2015.4573392746413,1914.6772901631812,3124.3575312007774),c(2368.99307782762,2459.169009518945,2500.5788765930074,2590.8415083290656),c(2684.487651451353,2622.0157668457205,2670.470957332858,2609.112605002754),c(2357.8579752291353,2511.3165666404407,2536.7536162182228,2737.0102817185757))
targetgene="ENST00000284142"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2774.496397455771,2966.9215393861377,2731.8388191970625,2123.1014334826336),c(1655.4185863080586,1672.381305580592,1726.0518621174158,1967.797111756279),c(1275.8971727430385,1266.1792816868378,1262.2400219226909,1627.954713625668),c(2033.084149439998,2181.963574294153,2031.5992164518234,2400.8221029227025),c(2136.0838484759815,2212.1542652592298,2148.5211427404656,2656.617456354345),c(3284.8552665529864,3371.293824433524,3346.8093928257226,4100.034093575757),c(2319.813041350979,2458.2541400957607,2422.4155999028103,3614.0229220556357))
targetgene="ENST00000376708"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(3091.846821512585,3177.3415067184883,3005.733276359407,4522.096426737967),c(1623.8691289456854,1665.9772196183033,1754.474871822942,2194.3587105100196),c(1145.0597172108432,1063.9931391631446,1091.0559861962256,1101.747129423432),c(2335.5877700321657,2286.2586885371443,2342.3143907326908,2698.640978703829),c(2533.235841155269,2546.0816047214375,2439.211014728803,3080.5068991839244),c(1686.968043670432,1846.206495985577,1903.695672776955,2265.6159875374055),c(1926.372749537853,1851.6957125246818,1873.3347305915065,2205.3213685142327))
targetgene="ENST00000215375"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2190.8314362518645,2377.7456308555575,2275.132731428719,2541.509547310106),c(3089.990971079504,3299.9340094251615,3142.0345274472716,3961.173758855722),c(1983.9041129633572,1956.9056961908568,2007.698049199449,2541.509547310106),c(1850.282881781541,1864.5038844492597,1833.2841260064467,1960.4886730868036),c(3081.6396441306406,2769.3097439783655,2977.956244147188,3365.5360072934686),c(2804.1900043850633,2675.993062813584,2723.441111784066,2671.2343336932963),c(1934.7240764867165,2074.008982358426,2018.0336890923677,2424.574528598498),c(1823.3730505018696,1958.7354350372252,1831.3461935265245,2578.0517406574836),c(2478.488253379386,2652.2064578107966,2525.771998831997,3117.0490925313015),c(1720.3733514658861,1727.2734709716399,1783.5438590217757,2554.299314981688))
targetgene="ENST00000400385"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2040.5075511723212,2106.0294121698703,2144.645277780621,2280.232864876357),c(1899.4629182581816,2063.9454187034003,2065.190046103809,2764.4169267291086),c(1558.9143637878578,1579.0646244158106,1542.5942540181095,1922.119370072057),c(1802.0307705214404,2139.879580827683,1913.3853351765663,2249.172000531086),c(2350.434573496812,2548.82621299099,2443.0868796886475,3018.3851704933822),c(2520.2448881237037,2707.0986232018445,2695.01810207854,3062.2358025102353),c(2076.6966346173963,2111.518628708975,2007.0520717061418,2762.5898170617397))
targetgene="ENST00000334269"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2227.94844491348,2163.6661858304706,2213.764869564515,2616.4210436722296),c(1744.4994070959362,1754.7195536671638,1629.8012156146108,1969.624221423648),c(1586.7521202840696,1601.0214905722298,1540.01034404488,1726.6186356635874),c(2121.2370450113353,2105.1145427466863,2281.5925063617933,2709.6036367080424),c(2248.362799677369,2237.770609108385,2267.3810015090303,2682.1969916975095),c(1980.1924120971958,2078.5833294743466,1927.5968400293295,2470.25227028272),c(2973.072393795415,2967.836408809322,3103.275877848827,3526.32165802193))
targetgene="ENST00000306704"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1780.6884905410116,1973.3733458081713,1949.5600748017816,2861.253739099659),c(2339.2994708983274,2555.2302989532786,2467.634024434329,3212.058795234483),c(1750.994883611719,1860.844406756523,1751.8909618497123,2093.8676788047314),c(3618.908344507527,3879.0463543007168,3795.1177731810694,5678.656846182465),c(1897.6070678251008,1890.1202282984154,1835.222058486369,2031.7459501141896),c(2314.245490051737,2469.2325731739706,2454.714474568181,2810.0946684133305),c(2008.0301685934076,2164.5810552536545,2016.0957566124455,168.09408939793656))
targetgene="ENST00000292807"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1647.9951845757355,1689.7638246210906,1540.01034404488,1686.4222229814723),c(1652.6348106584376,1622.978356728649,1772.5622416355498,2406.303431924809),c(2390.3353578080487,2499.4232641390468,2364.92360299845,3219.367233903959),c(2284.551883122444,2369.5118060469003,2325.518975906698,2475.7335992848266),c(1141.3480163446816,1180.1815559075294,1191.1824976588753,1596.8938492803973),c(2924.820282535314,3085.854564400075,2988.2918840401067,3162.7268342155235),c(2009.8860190264882,1929.459613495333,1960.5416921880078,2199.840039512126),c(1699.9589967019974,1803.2076330959228,1730.5737045705675,2150.5080784931665),c(1518.0856542600807,1567.171321914417,1463.1390223412973,1684.5951133141034),c(2368.99307782762,2320.108857194957,2340.3764582527683,2848.463971428077))
targetgene="ENST00000562866"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2763.3612948572863,2578.102034532882,2667.8870473596285,2861.253739099659),c(627.2774463813047,687.9818062344667,696.3637377853945,765.5589506275588),c(1773.2650888086885,1901.098661376625,1772.5622416355498,2168.779175166855),c(1638.7159324103316,1773.9318115540307,1741.5553219567937,1960.4886730868036),c(2242.7952483781264,2293.5776439226174,2216.3487795377446,2810.0946684133305),c(1912.4538712897472,1880.05666464339,1835.8680359796765,2099.349007806838),c(2296.614910937469,2349.3846787368493,2258.3373166027263,2579.8788503248525))
targetgene="ENST00000203166"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2515.605262041002,2443.616229324815,2499.9328990997,2906.931480783881),c(1923.5889738882317,2070.3495046656894,1985.734814426997,1814.3198996972935),c(2160.209904106032,2042.9034219701653,2059.3762486640426,2705.9494173733046),c(2359.713825662216,2235.940870262017,2343.6063457193054,2747.972939722789),c(2552.722270702617,2465.573095481234,2643.3399026139464,3091.4695571881375),c(2603.758157612339,2486.615092214469,2543.213391151297,3053.100254173391),c(2546.2267941868345,2376.830761432373,2470.8639119008662,2780.8609137354288))
targetgene="ENST00000229379"
collabel=c("Sym-pool-1-1","Sym-pool-1-2","Sym-pool-1-3","sym-03-03-TP6")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



dev.off()
Sweave("transcript_pooled_screen_1_counts_mageck_summary.Rnw");
library(tools);

texi2dvi("transcript_pooled_screen_1_counts_mageck_summary.tex",pdf=TRUE);

