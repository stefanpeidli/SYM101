pdf(file='gene_SYM101_counts_mageck.pdf',width=4.5,height=4.5);
gstable=read.table('gene_SYM101_counts_mageck.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("RAN","SNHG16","EIF4A1","RPS5","RPS29","ECT2","RPLP1","PPIA","HNRNPA1","FBL")
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
targetmat=list(c(2380.600521934609,2267.460330554884,2210.907119152687,9.389860779631501),c(2499.238356841564,2589.168377454124,2435.3992266358828,1.5649767966052501),c(559.8529234038969,593.0280864528158,621.0948307035087,23.47465194907875),c(2828.6789562526137,2867.271417996539,2802.069668858436,90.7686542031045),c(762.81186411249,824.6191202146181,802.7293540308217,106.41842216915701),c(2526.6917401258183,2537.811369967197,2457.8484373842025,29.73455913549975),c(2973.7896964693855,2963.202431981553,3069.419360497515,0.0),c(2879.6638109233713,2945.7604294388234,3100.7121997224453,50.079257491368004),c(2421.78059686099,2422.5003531569273,2390.5008051392438,7.82488398302625),c(3649.339020856925,3715.1465416014635,3746.977357628615,26.604605542289253),c(2155.0905878139497,2321.724338465599,2361.2488032550696,2173.7527704846925),c(782.4214236012429,813.9601186607275,793.2054464406256,114.24330615218327),c(961.8488929233324,992.2561446530774,920.4176406811032,28.169582338894504),c(2144.3053300951356,2098.854305975162,2036.7556660748137,67.29400225402576),c(2535.516041895757,2616.3003814094814,2497.304625972158,1.5649767966052501),c(1434.4392766022786,1431.2132086451127,1487.770421411362,1746.5141050114592),c(3368.922320167758,3625.029528464026,3619.084884274552,15.6497679660525),c(2714.9435112178467,2923.4734261897797,2835.403345424123,139.28293489786725),c(963.8098488722077,962.2171402739315,948.3090843381063,25.039628745684002),c(587.306306688151,600.780087582918,591.1625497057493,20.344698355868253),c(608.8768221257792,581.4000847576625,600.6864572959454,15.6497679660525))
targetgene="RAN"
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
targetmat=list(c(1340.3133910562644,1448.6552111878425,1388.4496708278873,14.084791169447252),c(1589.354796563427,1630.8272377452433,1589.1320093355928,112.67832935557801),c(2137.441984274072,2217.07232320922,2254.444982422155,115.80828294878852),c(2760.045498041978,3029.094441587422,2803.430227085607,104.85344537255176),c(2294.318460184096,2471.9193603613285,2402.7458291837816,26.604605542289253),c(1965.8583387474837,2202.537321090278,2089.1371578208928,101.72349177934126),c(2493.3554889949382,2508.7413657293137,2370.772710845266,921.7713332004923),c(2345.3033148548534,2351.763342844745,2311.5884279633324,125.19814372842),c(2290.396548286345,2205.4443215140664,2260.567494444424,31.299535932105),c(1921.7368298977894,1924.434280547863,1796.6171389791527,90.7686542031045))
targetgene="SNHG16"
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
targetgene="EIF4A1"
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
targetgene="RPS5"
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
targetgene="RPS29"
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
targetmat=list(c(1236.3827257658738,1291.6771883032736,1281.6458499949729,953.0708691325973),c(2405.1124712955498,2556.2223726511897,2544.2438848095535,302.04052174481325),c(1860.9471954826554,1957.3802853507973,1960.5644053532442,1430.3887920971986),c(2629.6419274417713,2594.9823783017005,2566.693095557873,636.9455562183368),c(1254.0313293057513,1195.7461743182594,1190.4884487745237,46.9493038981575),c(2771.81123373523,3009.7144387621665,2844.927253014319,244.136380270419),c(3093.4080093507787,3415.7254979512672,3172.1415066489167,2608.816319940952),c(1260.894675126815,1278.1111863255949,1244.910777861359,982.8054282680971),c(2487.472621148312,2459.3223585249125,2478.937089905351,145.54284208428825),c(2248.235995385526,2251.9563282946797,2182.335396382098,1010.9750106069915),c(1870.7519752270318,1989.3572900124686,1921.788495878874,336.4700112701288),c(2758.084542093103,2778.123405000364,2693.9052897983506,50.079257491368004),c(1204.0269526094314,1183.1491724818434,1240.8291031798465,43.819350304947),c(1798.196605118646,1902.1472772988193,1889.1350984267729,666.6801153538365))
targetgene="ECT2"
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
targetmat=list(c(2228.6264358967733,2265.5223302723584,2188.457908404367,162.75758684694603),c(2805.14748486611,2968.047432687867,2939.4860498026956,40.68939671173651),c(2523.7503062025053,2650.2153863536782,2609.550679713756,1054.7943609119386),c(1149.120186040923,1125.009164006077,1153.0730975273243,262.916101829682),c(2287.455114363032,2241.2973267407892,2249.002749513472,2327.120496552007),c(1832.5133342239635,1775.2082587933962,1897.9787269033836,812.2229574381248),c(2457.077803940745,2512.6173662943647,2555.1283506269206,1305.1906483687785),c(2092.3399974499403,1982.5742890236293,2027.912037598203,18.779721559263002),c(2255.0993412065895,2144.397312614512,2072.130179981257,7.82488398302625),c(2862.015207383494,2755.8364017513204,2746.9670606580153,56.33916467778901),c(2866.917597255682,2950.6054301451372,2775.538783428604,107.98339896576226),c(2825.7375223293006,2880.837419974218,3036.7659630454136,120.50321333860425),c(1098.1353313701654,1137.606165842493,1125.8619329839066,223.79168191455076),c(2121.7543366830696,2057.1872999008624,2003.4219895091271,1513.3325623172768))
targetgene="RPLP1"
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
targetmat=list(c(2400.2100814233618,2336.2593405845405,2310.2278697361617,538.352018032206),c(2082.535217705564,2104.6683068227385,2165.3284185424623,10.95483757623675),c(2772.7917117096677,3000.9934374908016,2850.3694859230027,17.21474476265775),c(1844.2790699172153,1957.3802853507973,1883.6928655180893,50.079257491368004),c(2916.921973952002,2933.1634276024074,2719.0756170010122,640.0755098115472),c(2841.4251699203032,2992.2724362194367,2976.2211219363094,25.039628745684002),c(2047.2380106258086,2189.940319253862,2100.02162363826,15.6497679660525),c(2269.8065108231544,2084.3193038562204,2102.7427400926017,528.9621572525746),c(2065.8670921401235,2073.66030230233,2045.5992945514245,472.62299257478554),c(2820.8351324571126,3000.9934374908016,2900.7101403283255,2945.286331211081),c(2158.0320217372628,2108.5443073877896,2195.2606995402216,1518.0274927070927),c(3092.427531376341,3251.964474077859,3092.5488503594197,2763.749022804872),c(2994.379733932576,2940.9154287325096,3036.7659630454136,1309.8855787585944),c(1436.400232551154,1485.4772165558277,1460.5592568679442,45.38432710155225))
targetgene="PPIA"
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
targetgene="HNRNPA1"
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
targetmat=list(c(2760.045498041978,2749.053400762481,2834.042787196952,76.68386303365726),c(2595.3251983364535,2602.7343794318026,2433.3583892951265,181.537308406209),c(1053.0333445460337,1076.5591569429384,1106.8141178035144,1070.444128877991),c(1646.2225190808106,1628.8892374627178,1721.1061573711686,1975.0007173158256),c(1817.8061646073988,1978.6982884585782,1868.7267250192094,2128.3684433831404),c(1299.1333161298833,1243.227181240135,1277.5641753134603,56.33916467778901),c(2716.904467166722,2788.7824065542545,2897.9890238739836,18.779721559263002),c(2351.186182701479,2446.7253566884965,2368.051594390924,89.20367740649925),c(2957.1215709039457,3184.1344641894652,3142.8895047647425,1082.9639432508332),c(1084.4086397280385,1078.497157225464,1119.7394209616377,1128.3482703523853),c(1237.3632037403115,1304.2741901396896,1288.4486411308274,54.77418788118376),c(2182.5439710982037,2318.8173380418107,2272.8125184889623,389.6792223547073),c(2653.173398828275,2781.0304054241524,2773.4979460878476,26.604605542289253),c(2676.704870214778,2703.510394123131,2728.5995245912086,64.16404866081525))
targetgene="FBL"
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
startindex=9
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("CYP2A13","SPDYE1","CABP4","FNDC7","TUSC8","SLC22A6","AATBC","BIN1","COX8A","LINC01433")
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
targetgene="CYP2A13"
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
targetgene="SPDYE1"
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
targetgene="CABP4"
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
targetgene="FNDC7"
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
targetgene="TUSC8"
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
targetgene="SLC22A6"
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
targetgene="AATBC"
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
targetmat=list(c(3213.0263222321714,3275.2204774681654,3136.7669927424736,3006.3204262786853),c(2872.800465102308,3141.4984579739034,3000.0308909118,3054.834706973448),c(1588.3743185889894,1646.3312400054479,1652.397966899039,1402.219209758304),c(728.4951350071723,782.9521141403189,749.6675831711573,870.127098912519),c(1336.391479158514,1270.3591851954927,1276.883896199875,1447.6035368598564),c(2809.069396763861,2816.883410650875,2750.3684562259427,2795.0485587369767),c(2053.1208784724345,2054.280299477074,1891.1759357675292,2133.063373772956),c(2612.973801876331,2675.4093900265107,2701.3883600477907,3089.264196498764),c(1974.6826405174224,2046.5282983469722,1905.4617971528235,1605.6661933169867),c(919.6883400225136,911.8291329282674,974.8399697679386,1067.3141752847805),c(808.8943289110593,878.8831281253332,854.4305666633153,768.4036071331778),c(1189.3197829928667,1231.5991795449818,1176.2025873892294,1472.6431656055404),c(2145.2858080695732,2161.839315157242,2235.397167241763,2957.806145583923),c(1668.7735124928765,1683.153245373433,1867.3661667920387,2688.6301365678196),c(1098.1353313701654,1186.0561729056317,1187.767332320182,1380.3095346058305),c(2326.674233340538,2243.2353270233148,2349.684058324117,1345.8800450805152),c(895.1763906615723,866.2861262889172,855.7911248904862,1122.0883631659644),c(731.4365689304852,820.743119649567,776.878747714575,837.2625861838088),c(772.6166438568664,909.8911326457419,816.3349363025305,741.7990015908886),c(2284.5136804397193,2616.3003814094814,2513.6313246982086,2877.992328957055),c(1156.9640098364243,1170.5521706454272,1194.5701234560363,1516.4625159104874),c(2213.9192662802084,2297.49933493403,2230.6352134466647,2494.5730137887685),c(1088.330551625789,1093.0321593444055,1152.392818413739,1342.7500914873046),c(2470.8044955828723,2383.7403475064166,2308.8673115089905,962.4607299122289),c(916.7469060992006,897.2941308093259,913.6148495452488,1111.1335255897277),c(743.202304623737,718.0291046757133,751.7084205119136,902.9916116412293),c(1236.3827257658738,1251.9481825115,1186.406774093011,1428.8238153005934),c(818.6991086554359,912.7981330695302,834.0221932557521,744.9289551840991))
targetgene="BIN1"
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
targetgene="COX8A"
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
targetgene="LINC01433"
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
startindex=15
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("RAN","SNHG16","RPS5","RPS29","EIF4A1","RPLP1","HNRNPA1","EEF1A1","PPIA","FBL")
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
targetmat=list(c(2293.8339791807034,2173.671535304004,2128.113762253954,40.22830432222023),c(2408.147781273317,2482.072795868504,2344.199159652048,36.73019090289674),c(539.4477768172083,568.4987092333549,597.8362661347262,1.7490567096617493),c(2725.581149067681,2748.6726807540804,2697.1386420689346,20.98868051594099),c(735.009405190522,790.5104600614135,772.6689967568202,139.92453677293994),c(2434.600561922847,2432.840064513328,2365.807699391857,66.46415496714647),c(2865.40298964377,2840.6357072477113,2954.476706243028,34.981134193234986),c(2774.707741702523,2823.915156976142,2984.597701031853,71.71132509613172),c(2333.513150154999,2322.2986488290644,2300.982080172429,59.46792812849948),c(3516.3303420554275,3561.477207844253,3606.6617238445474,66.46415496714647),c(2076.5432809881327,2225.691025037775,2272.825498087223,2551.873739396492),c(753.904248511615,780.2923460065656,763.5017374732647,157.41510386955744),c(926.792064899617,951.2135265603847,885.9501293321846,22.737737225602743),c(2066.1511171615316,2012.0395493455012,1960.4838782117965,20.98868051594099),c(2443.103241417339,2508.0825407353896,2403.7863449951587,27.98490735458799),c(1382.157788937961,1372.0140417282112,1432.0568609382763,1766.5472767583667),c(3246.1340825637963,3475.087698107812,3483.5585277510877,71.71132509613172),c(2615.991057805341,2802.550009406915,2729.2240495613787,83.95472206376397),c(928.6815492317263,922.4170233149043,912.7971029483114,8.745283548308747),c(565.9005574667386,575.930064909608,569.0248798149803,0.0),c(586.684885119941,557.3516757189755,578.1921390985358,0.0))
targetgene="RAN"
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
targetmat=list(c(1291.462540996714,1388.7345919997804,1336.4554426954833,27.98490735458799),c(1531.4270511745965,1563.371450391726,1529.6226918846883,43.726417741543735),c(2059.5379219991487,2125.36772340836,2170.0212332644933,148.6698203212487),c(2659.449197443855,2903.8022304958617,2698.448250538014,57.71887141883773),c(2210.6966685678935,2369.673541265177,2312.768556394143,19.23962380627924),c(1894.2080429395842,2111.433931515385,2010.9038042713516,59.46792812849948),c(2402.479328276989,2404.9724807273788,2281.9927573707782,1026.6962885714468),c(2259.823261202736,2254.4875282832554,2225.0247889658262,136.42642335361646),c(2206.917699903675,2114.22068989398,2175.914471375351,78.70755193477872),c(1851.6946454671247,1844.8340466298087,1729.3379834192901,62.966041547822975))
targetgene="SNHG16"
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
targetgene="RPS5"
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
targetgene="RPS29"
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
targetgene="EIF4A1"
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
targetmat=list(c(2147.398943442232,2171.813696384941,2106.505222514145,85.70377877342571),c(2702.9073370823694,2845.2803045453693,2829.4090974459496,5.247170128985248),c(2431.766335424683,2540.5947218189963,2511.8290436942057,895.5170353468156),c(1107.237818616056,1078.4754925162174,1109.8931775447545,239.62076922365966),c(2204.083473405511,2148.5907098966504,2164.782799388176,2578.1095900414184),c(1765.723108356151,1701.7804498619382,1826.9038143657021,874.5283548308746),c(2367.5238681329665,2408.6881585655055,2459.444704931031,1290.803851730371),c(2016.0797823606347,1900.5692142017062,1951.9714231627806,50.72264458019073),c(2172.9069819257074,2055.698763943488,1994.5336984078597,15.741510386955744),c(2757.702382713539,2641.8469429079437,2644.0994990712206,150.41887703091044),c(2762.4260935438124,2828.5597542738,2671.601276921887,143.42265019226343),c(2722.7469225695168,2761.6775531875232,2923.0461029851235,125.93208309564595),c(1058.1112259812141,1090.5514454901286,1083.7010081631674,220.38114541738042),c(2044.4220473422745,1972.0960125856413,1928.3984707193522,1362.5151768265027))
targetgene="RPLP1"
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
targetgene="HNRNPA1"
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
targetgene="EEF1A1"
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
targetmat=list(c(2312.7288225017965,2239.6248169307496,2223.715180496747,594.6792812849948),c(2006.6323607000882,2017.613066102691,2084.241878539796,15.741510386955744),c(2671.7308456025658,2876.8635661694448,2743.6297427212517,27.98490735458799),c(1777.0600143488068,1876.417308253884,1813.1529254403688,78.70755193477872),c(2810.6079440126,2811.839204002231,2617.2525254550937,746.8472150255669),c(2737.8627972263916,2868.50329103366,2864.768526111092,52.471701289852476),c(1972.6216427221207,2099.357978541474,2021.3806720239866,29.73396406424974),c(2187.0781144165276,1998.1057574525269,2023.9998889621452,407.53021335118757),c(1990.571743877159,1987.887643397679,1968.9963332608122,411.0283267705111),c(2718.023211739244,2876.8635661694448,2792.0852560771878,2787.9963952008284),c(2079.3775074862965,2021.3287439408175,2113.0532648595413,1154.3774283767546),c(2979.716791736383,3117.453706188136,2976.740050217377,2639.3265748795798),c(2885.242575130918,2819.270559678484,2923.0461029851235,1770.0453901776903),c(1384.0472732700703,1424.0335314619822,1405.864691556689,8.745283548308747))
targetgene="PPIA"
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
targetmat=list(c(2659.449197443855,2635.3445066912223,2727.914441092299,125.93208309564595),c(2500.732513546673,2495.0776683019467,2342.234746948429,69.96226838646997),c(1014.6530863426999,1032.029519539636,1065.3664895960565,1024.947231861785),c(1586.2220968057666,1561.513611472663,1656.6547133853858,1803.2774676612635),c(1751.551975865331,1896.8535363635797,1798.7472322804958,1766.5472767583667),c(1251.7833700224185,1191.8036665790758,1229.7223524655155,89.20189219274921),c(2617.8805421374504,2673.4302045320187,2789.466039139029,17.490567096617493),c(2265.4917141990636,2345.521635317355,2279.37354043262,131.1792532246312),c(2849.342372820841,3052.429344020922,3025.1955635733134,893.7679786371539),c(1044.884835656449,1033.8873584586993,1077.8077700523104,1056.4302526356967),c(1192.2646135609752,1250.3255925295682,1240.1992202181505,103.19434587004321),c(2102.996061637663,2222.9042666591804,2187.700947597065,300.8377540618209),c(2556.4723013438975,2665.9988488557656,2669.636864218268,40.22830432222023),c(2579.146113329209,2591.6852920932356,2626.419784738649,0.0))
targetgene="FBL"
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
startindex=21
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("OCIAD2","RET","STARD8","LSM11","SPRR2G","SPRN","FOXB2","PCBP4","PCAT1","NTRK1")
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
targetmat=list(c(2548.9143640154603,2587.9696142551093,2397.8931068843017,2254.534098753995),c(852.1574337812992,906.6253925028667,918.6903410591685,1091.4113868289317),c(1299.0204783251513,1327.425907670693,1327.9429876464674,1703.5812352105438),c(2137.0067796156304,2281.426192609673,2154.305931635541,2705.790729846726),c(2746.365476720883,2636.2734261507535,2605.4660492333796,3099.3284895206198),c(2195.580793911019,2300.9335012598367,2316.697381801381,2964.651122876665),c(3142.212444297784,3062.64745807577,3105.081680187154,3468.379455259249),c(1752.4967180313859,1846.6918855488718,1807.9144915640513,1866.2435092090866),c(878.6102144308296,934.4929762888155,923.9287749354859,1140.3849746994606),c(2728.415375565845,2793.2608148115983,2911.914430997949,3167.541701197428),c(2270.215425029337,2317.654051531406,2277.409127729001,2562.3680796544627),c(1296.1862518269872,1312.563196318187,1341.039072337261,1752.5548230810728),c(2096.3828664752805,2250.771850445129,2206.6902703987153,2266.777495721627),c(3021.2854470427883,2868.50329103366,2881.138631974584,3410.660583840411))
targetgene="OCIAD2"
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
targetmat=list(c(3013.727509714351,3115.5958672690726,2794.7044730153466,3319.7096349380004),c(2326.899954992616,2349.2373131554814,2281.337953136239,2572.8624199124333),c(2000.9639077037602,2107.718253677259,1912.0283648558602,2366.473728172347),c(2265.4917141990636,2286.9997093668626,2301.636884406969,2697.0454462984176),c(2090.7144134789523,2089.1398644866263,2083.587074305256,2399.70580565592),c(1879.0921682827097,1957.2333012331353,1868.8112853762416,2348.983161075729),c(1919.71608142306,1851.33648284653,1870.1208938453208,2174.077490109554))
targetgene="RET"
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
targetmat=list(c(2346.739540479764,2460.7076482992766,2497.4233505343327,3622.2964457094827),c(2954.2087532529076,3091.4439613212503,2944.6546427249327,3905.643632674686),c(2615.046315639286,2541.523641278528,2509.8646309905866,3398.417186872779),c(2136.062037449576,2092.8555423247526,2035.7863651838595,2039.4001234655998),c(1808.2365058286105,1985.100885019084,2044.2988202328752,2200.3133407544806),c(839.8757856225886,891.7626811503607,819.8149016436771,963.7302470236239),c(1715.6517735552543,1790.0277985174428,1696.5977716923062,1820.768034757881),c(1675.9726025809587,1672.9839466164578,1762.7329993808137,1423.7321616646639),c(1934.8319560799343,1985.100885019084,1999.7721322841771,2121.605788819702),c(2215.420379398167,2260.061045040445,2233.537244014842,3164.0435877781047),c(2463.8875690705413,2577.7515002002615,2642.789890602141,2383.964295268964),c(1846.9709346368513,1960.9489790712619,1945.4233808173838,1943.2020044342034),c(2207.8624420697297,2325.085407207659,2242.704503298398,3024.1190510051647),c(844.599496452862,892.6916006098922,871.5444361723116,1017.9510050231381))
targetgene="STARD8"
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
targetgene="LSM11"
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
targetgene="SPRR2G"
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
targetgene="SPRN"
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
targetgene="FOXB2"
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
targetmat=list(c(2085.0459604826247,2038.9782136719184,2108.469635217764,2611.341667524992),c(3329.2713931766057,3425.8549667526354,3404.3272153717867,3031.1152778438113),c(1772.3363035185337,1934.9392342043764,1951.9714231627806,2191.568057206172),c(2872.9609269722073,2930.740894822279,2933.5229707377584,2614.8397809443154),c(2286.276041852266,2322.2986488290644,2297.053254765191,2238.792588367039),c(958.9132985454753,998.5884189964976,1073.8789446450721,1182.3623357313425),c(1791.2311468396267,1906.142730958896,1870.7756980798606,1579.3982088245596),c(1978.2900957184486,2067.7747169173986,2059.359317627288,2322.7473104308033),c(2974.0483387400554,3016.2014850991886,2964.298769761123,3505.1096461621455),c(2461.998084738432,2495.0776683019467,2452.241858351095,2709.2888432660498),c(2131.3383266193027,2268.42132017623,2124.184936846716,2231.796361528392),c(1064.7244211435966,1160.2204049550005,1101.3807224957386,1560.1585850182803),c(945.68690822071,946.5689292627266,922.6191664664066,1002.2094946361824),c(2545.1353953512416,2520.1584937093003,2509.209826756047,2801.9888488781226),c(3079.8594613381765,2951.1771229319747,2984.597701031853,4124.275721382405),c(2575.3671446649905,2752.3883585922067,2619.2169381587128,2630.581291331271),c(1058.1112259812141,1151.860129819216,1098.76150555758,1596.888775921177),c(945.68690822071,996.7305800774344,945.5373146752953,982.9698708299031),c(2178.5754349220356,2342.73487693876,2281.9927573707782,3167.541701197428),c(1733.6018747102928,1785.3832012197845,1766.0070205535121,1694.8359516622352),c(1016.5425706748092,1021.8114054847882,996.6120449693902,1213.8453565052541))
targetgene="PCBP4"
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
targetmat=list(c(2833.2817559979117,2901.944391576799,2839.2311609640446,2882.4454575225627),c(2369.413352465076,2482.072795868504,2539.9856257794117,2523.888832041904),c(3060.9646180170835,3157.397242947996,3023.8859551042337,3714.9964513215555),c(1813.9049588249384,1930.294636906718,1793.5087984041786,1551.4133014699717),c(2684.9572359273307,2582.3960974979195,2539.9856257794117,3442.1436046143226),c(2334.457892321054,2250.771850445129,2304.910905579667,2912.1794215868126),c(2796.4368115217803,3034.7798742898212,3028.4695847460116,3393.1700167437934),c(3264.0841837188345,3130.4585786215785,3127.345024161503,4071.8040200925525),c(1964.1189632276287,1750.084261757583,1965.722312088114,1806.775581080587),c(2552.693332679679,2412.403836403632,2481.0532446708407,2553.622796106154),c(1518.2006608498314,1401.7394644332232,1540.0995596373232,1542.6680179216628),c(1773.2810456845882,1970.2381736665782,1993.2240899387805,1993.9246490143942),c(1205.4910038857404,1200.1639417148604,1180.6120348750399,1304.796305407665),c(1955.6162837331367,2069.632555836462,1999.1173280496375,2641.0756315892413),c(3009.9485410501325,3181.549148895818,3188.8966222082327,4463.592723056784),c(1826.186606983649,1833.6870131154292,1807.2596873295117,1992.1755923047324),c(1862.0868092937258,1910.787328256554,1916.611994497638,2327.9944805597884),c(1908.379175430404,1790.9567179769742,1833.4518567110988,2242.2907017863627),c(2830.447529499748,2838.7778683286483,2803.871732298902,2611.341667524992),c(2378.8607741256224,2316.7251320718747,2229.608418607604,2385.713351978626))
targetgene="PCAT1"
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
targetgene="NTRK1"
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
startindex=27
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("RAN","SNHG16","RPS5","RPS29","EIF4A1","PPIA","HNRNPA1","RHEB","FBL","EEF1A1")
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
targetmat=list(c(2218.604963644284,2079.31859915413,2047.8497958782673,9.596359762783855),c(2329.1697085375945,2374.3330328802713,2255.785313613599,24.950535383238027),c(521.7559449097554,543.8217874710801,575.2882657344178,9.596359762783855),c(2636.19247121654,2629.3605704688334,2595.4133259146406,11.515631715340627),c(710.9038969173199,756.1966358462241,743.5270028111862,38.38543905113542),c(2354.75493876084,2327.237355207122,2276.5788653871323,61.41670248181668),c(2771.4286881108374,2717.331741971508,2843.045624308536,23.031263430681253),c(2683.7078987739956,2701.3369835164763,2872.0305752655822,49.90107076647605),c(2256.982808979152,2221.4942298655233,2214.1982100665327,17.27344757301094),c(3401.008103247127,3406.8835509217665,3470.632823291537,88.28650981761147),c(2008.4405725247677,2129.0800699031174,2187.1035819979893,2107.3606039073347),c(729.1790613624953,746.4220612348158,734.7054959981722,47.981798813919276),c(896.3968160358495,909.9240365529183,852.5356227148602,26.869807335794796),c(1998.3892320799214,1924.7026007554894,1886.5422427260098,78.69015005482761),c(2362.978762761169,2399.213768254765,2313.1251078981904,1.9192719525567712),c(1336.8282791645747,1312.4587910045511,1378.0453857186985,1046.0032141434403),c(3139.67325168112,3324.243965570769,3352.1725889453483,23.031263430681253),c(2530.1965174345232,2680.8992366017133,2626.28859976019,255.26316969005057),c(898.224332480367,882.3775081025858,878.3700355244014,42.22398295624897),c(547.3411751330009,550.9305690066498,547.5635300363735,9.596359762783855),c(567.4438560226937,533.1586151677255,556.3850368493877,24.950535383238027))
targetgene="RAN"
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
targetmat=list(c(1249.1074898277332,1328.453549459583,1286.049671811552,13.434903667897398),c(1481.2020782814598,1495.5099155454702,1471.9314225143485,40.3047110036922),c(1991.9929245241099,2033.1115191729268,2088.1766841663316,211.11991478124483),c(2572.2293956584263,2777.7563850238503,2596.6735411736427,205.3620989235745),c(2138.1942400855123,2266.81271215478,2225.5401473975508,5.757815857670313),c(1832.0852356288258,2019.7825537937338,1935.0605301975872,72.9323341971573),c(2323.687159204042,2300.579424448736,2195.9250888110037,713.9691663511189),c(2185.7096676429683,2156.62659835345,2141.105725044416,101.72141348550888),c(2134.5392071964775,2022.4483468695723,2093.8476528318406,92.12505372272501),c(1790.9661156271814,1764.7550162051716,1664.114249512155,97.88286958039534))
targetgene="SNHG16"
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
targetgene="RPS5"
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
targetgene="RPS29"
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
targetgene="EIF4A1"
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
targetmat=list(c(2236.880128089459,2142.4090352823105,2139.845509785414,500.92997961731726),c(1940.822464077619,1930.0341869071667,2005.6325847017,21.111991478124484),c(2584.1082525477905,2751.98705195741,2640.1509676092123,40.3047110036922),c(1718.7792160687388,1794.9673377313427,1744.7680260882837,63.33597443437345),c(2718.430711219829,2689.7852135211756,2518.540195115518,568.1044979568043),c(2648.071328105904,2743.989672729894,2756.7208790668983,24.950535383238027),c(1907.9271680763034,2008.230783798433,1945.1422522696034,30.70835124090834),c(2115.350284529043,1911.3736353762963,1947.6626827876073,447.1903649457277),c(1925.28857429922,1901.599060764888,1894.733641909523,568.1044979568043),c(2628.88240543847,2751.98705195741,2686.7789321922864,2802.137050732886),c(2011.181847191544,1933.5885776749515,2033.3573203997441,1310.8627435962746),c(2881.993433004148,2982.1338541714786,2864.46928371157,3316.5019340181007),c(2790.6176107782712,2696.893995056745,2812.8004580924876,1189.9486105851981),c(1338.6557956090921,1362.2202617535388,1352.8410805386584,0.0))
targetgene="PPIA"
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
targetgene="HNRNPA1"
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
targetmat=list(c(2408.6666738741073,2319.2399759796062,2354.082103815756,1900.0792330312036),c(1568.9228676183013,1642.1285347165947,1655.292742699141,1583.3993608593362),c(2976.110529896801,2968.804888792285,2965.286504431731,32.62762319346511),c(1859.4979822965888,1910.48503768435,1827.3121255529154,23.031263430681253),c(2490.904913877396,2496.070916676902,2582.8111733246205,2740.7203482510695),c(2564.0055716580973,2607.145628170178,2580.290742806617,136.26830863153074),c(1884.1694542975756,1979.7956576561544,1952.7035438236155,24.950535383238027),c(1982.8553423015223,1949.5833361299833,1882.761596949004,2109.2798758598915),c(2092.506328972574,2145.9634260500957,2216.0885329550356,376.17730270112713),c(2726.654535220158,2704.891374284261,2723.955282332846,347.3882234127756),c(2387.6502347621554,2440.977859776237,2298.6326324196675,84.44796591249793),c(1892.3932782979043,1796.7445331152353,1830.4626637004203,28.78907928835157),c(2681.8803823294784,2461.4156066909995,2626.28859976019,383.8543905113542),c(2909.406179671911,2852.398591147332,2793.8972292074573,2199.48565763006))
targetgene="RHEB"
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
targetmat=list(c(2572.2293956584263,2520.9516520513957,2625.028384501188,42.22398295624897),c(2418.7180143189535,2386.7734005675184,2253.894990725096,94.04432567528178),c(981.3763307059147,987.2320357522385,1025.1851131981357,957.7167043258288),c(1534.2000551724682,1493.7327201615778,1594.1723026375435,1729.2640292536507),c(1694.1077440677523,1814.5164869541593,1730.9056582392616,1892.4021452209763),c(1210.729644492865,1140.0708387669865,1183.342128202888,69.09379029204376),c(2532.0240338790404,2557.3841574211906,2684.2585016742823,42.22398295624897),c(2191.1922169765207,2243.7091721641787,2193.4046582929996,46.06252686136251),c(2755.8947983324383,2919.9320157352436,2911.0972482946445,907.8156335593528),c(1010.6165938181953,989.009231136131,1037.1571581586547,990.344327519294),c(1153.1628764905627,1196.0524933595977,1193.423850274904,86.3672378650547),c(2034.0258027480131,2126.414276827279,2105.1895901628586,293.64860874118597),c(2472.6297494322207,2550.2753758856206,2568.9488054755984,61.41670248181668),c(2494.559946766431,2479.1875605299238,2527.361701928532,46.06252686136251))
targetgene="FBL"
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
targetgene="EEF1A1"
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
startindex=33
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("FENDRR","SPRR2G","KRT86","PCBP4","GMPPB","ASB17","HOXA11-AS","H19","SPDYE1","PCAT1")
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
targetgene="FENDRR"
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
targetgene="SPRR2G"
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
targetgene="KRT86"
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
targetmat=list(c(2016.6643965250967,1950.4719338219295,2028.946566993237,2556.470240805619),c(3220.0839752398915,3277.14828789762,3275.9295657757266,3299.22848644509),c(1714.210424957445,1850.948992323954,1878.3508435424967,1681.2822304397316),c(2778.7387538889075,2803.5257180902904,2822.8821801645036,2951.840263032314),c(2211.294897866214,2221.4942298655233,2210.4175642895266,2571.824416426073),c(927.4645955926475,955.242518842175,1033.3765123816486,859.8338347454335),c(1732.4855894026205,1823.4024638736214,1800.2174974843722,1971.092295275804),c(1913.409717409856,1978.018462272262,1981.6884947806616,1735.0218451113212),c(2876.5108836705954,2885.2767057493415,2852.497238751051,2491.214994418689),c(2381.253927206344,2386.7734005675184,2359.753072481265,2629.4025750027763),c(2061.438549415776,2169.955563732643,2044.0691501012611,2621.7254871925493),c(1029.8055164856294,1109.8585172408154,1059.841032820691,1282.073664307923),c(914.6719804810248,905.4810480931873,887.8216499669164,940.4432567528179),c(2461.6646507651158,2410.7655382500657,2414.5724362478522,2836.6839458789077),c(2978.8518045635774,2823.074867313107,2872.0305752655822,3120.73619485731),c(2490.904913877396,2632.914961236618,2520.430518004021,3022.853325276915),c(1023.409208929818,1101.8611380132995,1057.320602302687,1358.844542410194),c(914.6719804810248,953.4653234582826,909.8754169994517,974.9901518988397),c(2107.1264605287142,2241.04337908834,2195.9250888110037,2541.116065185165),c(1676.7463378448356,1707.8847639206142,1699.4002767642112,1842.5010744545004),c(983.2038471504322,977.4574611408302,959.0238121005301,794.5785883585032))
targetgene="PCBP4"
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
targetgene="GMPPB"
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
targetgene="ASB17"
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
targetmat=list(c(2579.5394614364964,2329.0145505910145,2332.6584444127216,1592.9957206221202),c(1815.6375876281681,1932.6999799830053,1904.815363981539,1558.4488254760981),c(1586.2842738412178,1528.38803014748,1682.3873707676842,1946.141759892566),c(2027.6294951922018,1982.461450731993,2004.3723694426978,1886.644329363306),c(1176.0068320470318,1201.3840795112749,1245.0926758939866,1303.1856557860476),c(764.8156320305872,789.9633481401801,781.3334605812465,909.7349055119096),c(678.0086009160044,754.4194404623317,814.0990573152989,957.7167043258288),c(1925.28857429922,2017.1167607178952,2005.0024770721989,1857.8552500749545),c(1292.9678844961538,1324.899158691798,1282.8991336640468,1710.0713097280832),c(2054.128483637706,2038.4431053246042,2114.641204605374,2143.8267710059135),c(2187.537184087486,2170.8441614245894,2222.389609250046,1953.8188477027932),c(2024.8882205254256,2004.6763930306481,2011.93366099671,2318.4805186885797),c(2554.8679894355096,2570.7131228003836,2594.153110655639,2723.4469006780582),c(1876.8593885195055,1782.5269700440958,1865.1185833229756,2076.6522526664266),c(1183.316897825102,1265.363113331402,1190.9034197569001,1320.4591033590586),c(804.1072355877141,813.0668881307815,844.9743311608481,884.7843701286715),c(762.0743573638108,797.0721296757497,811.5786267972949,921.2505372272502),c(1645.6785582880375,1708.7733616125604,1700.0303843937122,1303.1856557860476),c(1582.6292409521827,1595.921454735392,1560.7765982739902,1740.7796609689915),c(2144.590547641324,2376.1102282641637,2333.9186596717236,2151.5038588161406))
targetgene="HOXA11-AS"
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
targetgene="H19"
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
targetgene="SPDYE1"
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
targetmat=list(c(2740.3609085540393,2775.979189639958,2732.146681516359,2884.665744692827),c(2291.705621424985,2374.3330328802713,2444.1874948343993,2571.824416426073),c(2960.576640118402,3020.3435549251653,2909.8370330356424,3471.962962175199),c(1754.4157867368308,1846.506003864223,1725.8647972032536,1750.3760207317753),c(2596.900867659413,2470.301583610462,2444.1874948343993,3199.4263449121377),c(2257.896567201411,2153.072207585665,2217.9788558435384,3498.832769510994),c(2704.7243378859475,2903.0486595882658,2914.2477864421494,2291.6107113527846),c(3157.0346579040365,2994.574221858725,3009.3940384968014,2919.212639838849),c(1899.7033440759747,1674.1180516266584,1891.5831037620178,1942.3032159874524),c(2468.974716543186,2307.6882059843056,2387.4778081793092,1996.042830659042),c(1468.409463169837,1340.8939171468298,1482.0131445863644,1656.3316950564936),c(1715.1241831797038,1884.71570461791,1918.0476242010602,2477.7800907507917),c(1165.9554916021855,1148.0682179945024,1136.0840559903124,1604.5113523374607),c(1891.4795200756457,1979.7956576561544,1923.7185928665692,2372.2201333601693),c(2911.2336961164287,3043.4470949157667,3068.6241556698956,3283.8743108246354),c(1766.2946436261948,1754.0918439018171,1739.0970574227747,1740.7796609689915),c(1801.017456072028,1827.8454523333526,1844.3250315494424,2337.673238214147),c(1845.7916089627074,1713.2163500722916,1764.3013626028148,1777.24582806757),c(2737.619633887263,2715.5545465876157,2698.1208695233045,3138.009642430321),c(2300.843203647573,2216.162643713846,2145.516478450923,1729.2640292536507))
targetgene="PCAT1"
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



dev.off()
Sweave("gene_SYM101_counts_mageck_summary.Rnw");
library(tools);

texi2dvi("gene_SYM101_counts_mageck_summary.tex",pdf=TRUE);

