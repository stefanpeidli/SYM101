pdf(file='gene_pooled_screen_1_counts_filtered_mageck.pdf',width=4.5,height=4.5);
gstable=read.table('gene_pooled_screen_1_counts_filtered_mageck.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("RAN","FBL","TMSB4X","SNHG16","EIF4A1","RAD51","RPS12","RPS29","HNRNPA1","RPS16")
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
targetmat=list(c(4299.2513198661,4371.301935312703,4380.928054601457,52.256636964815186),c(3968.8950927082,4265.268789776948,4231.3974674073625,120.6653253551187),c(3392.5043327368985,3466.0297035343287,3625.3213639929995,121.61544602720625),c(3332.4395641627348,3373.6782541967364,3276.1515353855125,89.31134317622958),c(2987.067144861294,3078.3816445864104,2919.82332760384,7.60096537670039),c(2944.3287518373695,3046.4576867906994,2847.4441603981877,0.0),c(2853.076507272775,2850.353374617047,2794.9493797874948,28.503620162626465),c(2976.67131953115,2986.030195248818,2873.691550703534,85.5108604878794),c(3503.3931362584317,3486.5522478315715,3588.734092658274,437.05550916027244),c(2804.5626557321043,2667.9307586415557,2584.9702573447234,109.26387729006811),c(3198.448926574216,3439.806452487852,3315.124933111633,478.8608187321246),c(2526.1855552249226,2469.546163768209,2381.354138612339,377.19790681875685),c(1135.4551443923635,1132.160360397891,1108.7533965349367,108.31375661798056),c(1133.1449609856647,1167.5047422431423,1076.1430025192033,111.16411863424321),c(898.6613452057566,970.2602887196427,938.5430472820842,113.06435997841831),c(659.5573626124512,697.766506106253,726.1777984479177,11.401448065050586),c(691.8999303062317,706.8876369050276,691.1812780407892,9.501206720875487),c(921.7631792727426,957.7187338713277,927.4077907889069,120.6653253551187),c(717.3119477799163,684.0848099080912,702.3165345339664,19.002413441750974),c(1689.8991620000281,1683.9887737237511,1739.4861393270492,1559.1480228956675))
targetgene="RAN"
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
targetmat=list(c(3483.7565773014935,3746.5044755966464,3674.634642748499,63.65808502986577),c(3200.7591099809147,3281.3268048591444,3388.2994757810834,11.401448065050586),c(3125.67814926321,3272.20567406037,3242.7457659059805,1.9002413441750976),c(3251.5831449282837,3234.5810095154247,3313.534182184036,9.501206720875487),c(3153.400350143593,3180.994366072624,3190.250985295288,3.800482688350195),c(3057.527738765601,3062.419665688555,2845.0580340067927,0.0),c(2769.909904631625,2878.8569083632174,2768.7019894821483,14.251810081313232),c(2571.2341316555453,2728.358250183437,2657.3494245503757,314.4899424609786),c(1457.7257296268185,1534.630256893818,1506.441128434125,79.8101364553541),c(1530.4965069378245,1462.8013518534683,1493.7151210133509,76.0096537670039),c(1240.56848939715,1266.6970396798156,1294.0758795999584,56.05711965316538),c(1277.5314239043275,1268.9773223795091,1309.1880134121277,61.75784368569067))
targetgene="FBL"
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
targetmat=list(c(3018.254620851725,3151.350690976607,3102.7596842774665,23.75301680218872),c(2958.1898522775614,3080.661927286104,3031.971267999411,24.703137474276268),c(3167.261450583785,3173.0133766236963,3124.2348218000225,251.78197810320043),c(2694.8289439139207,2848.073091917353,2647.009543520997,106.41351527380547),c(2905.0556339234936,2825.2702649204166,2752.7944802061807,289.7868049867024),c(3064.458288985697,3041.8971213913123,2958.001349866162,419.95333706269656),c(3154.5554418469424,3285.8873702585315,3207.749245498852,569.1222825804417))
targetgene="TMSB4X"
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
targetmat=list(c(3251.5831449282837,3564.081859621155,3277.742286313109,303.0884943959281),c(2702.914585837366,2908.5005834592343,2809.2661381358657,100.71279124128017),c(2698.2942190239683,2594.961712251359,2643.0326662020047,177.67256568037163),c(2518.0999133014775,2608.6434084495213,2635.8742870278193,274.5848742333016),c(1872.403651129218,1918.857891792196,1857.9970834330072,132.06677342016928),c(2762.979354411529,2767.123056078229,2702.685825986883,599.5261440872433),c(2315.9588652153498,2591.541288201819,2442.598049324814,567.2220412362666),c(1579.0103584784952,1704.511318020994,1623.3613216124863,148.21882484565762),c(2263.979738564631,2264.320720795782,2100.586599891512,527.3169730085896),c(2937.398201617274,2951.8259547534135,2771.8834913373416,1231.3563910254632))
targetgene="SNHG16"
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
targetmat=list(c(3549.5968043924036,3635.910764661505,3527.4901819457996,0.0),c(3435.242725760823,3712.3002351012415,3466.2462712333245,0.0),c(3087.560123052683,3230.020444116037,2983.45336470771,31.35398217888911),c(3296.631721358907,3362.276840698268,3311.9434312564395,358.1954933770059),c(2853.076507272775,2780.804752276391,2736.091595466415,123.51568737138135),c(2930.467651397178,3048.737969490393,3028.7897661442175,314.4899424609786),c(2343.681066095733,2451.3039021706604,2425.0997891212496,591.9251787105429))
targetgene="EIF4A1"
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
targetmat=list(c(4387.038289320647,4335.957553467451,4511.3696306643915,210.92678920343585),c(3621.2124900000604,3711.160093751395,3556.9190741063394,110.21399796215566),c(3118.7475990431144,3232.300726815731,3296.035921980472,319.2405458214164),c(2745.6529788612897,2700.9948577871137,2689.1644431023105,158.67015223862066),c(3621.2124900000604,3777.2882920425104,3735.878553460974,829.4553467324301),c(2659.0211011100923,2850.353374617047,2778.246495047729,404.7514063092958),c(3302.407179875653,3538.9987499245253,3348.5307025911648,948.2204307433738),c(3107.1966820096213,3250.54298841328,3243.541141369779,936.8189826783231),c(2806.8728391388026,2852.6336573167405,2849.0349113257844,761.9967790142141),c(3322.0437388325913,3300.7092078065402,3333.4185687789954,1125.8929964237452),c(2840.3704985359327,2833.2512543693447,2724.160963509439,1073.6363594589302),c(3904.209957320639,3989.3545831140186,4003.920084761027,2280.289613010117),c(3205.379476794312,3535.578325874985,3310.352680328843,2012.3555834814283),c(3103.7314068995734,2948.405530703873,3097.987431494676,2345.847939384158))
targetgene="RAD51"
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
targetmat=list(c(3973.515459521597,3777.2882920425104,3848.026493856545,54.15687830899028),c(3497.617677741685,3577.763555819317,3529.0809328733963,290.7369256587899),c(2789.5464635885633,2869.7357775644427,2955.615223474767,10.451327392963037),c(2895.814900296699,2970.0682163509628,3000.9516249112744,115.91472199468096),c(2857.696874086172,2925.6027037069366,2957.2059744023636,342.0434419515176))
targetgene="RPS12"
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
targetmat=list(c(3615.437031483314,3667.834722457216,3617.367609355016,158.67015223862066),c(3122.212874153162,3039.616838691619,3016.0637587234432,45.60579226020234),c(2745.6529788612897,2822.989982220723,2803.698509889277,148.21882484565762),c(2437.2434940670264,2390.876410628779,2596.1055138379006,40.8551888997646),c(2555.062847808655,2613.2039738489084,2798.130881642688,131.11665274808175),c(2563.1484897321,2478.667294566984,2324.882480682654,198.5752204662977),c(1590.5612755119882,1800.2831914081266,1702.8988679923239,228.97908197309926))
targetgene="RPS29"
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
targetmat=list(c(3176.5021842105793,3340.614155051179,3220.475252919626,17.102172097575878),c(3237.7220444880923,3364.557123397962,3316.7156840392295,115.91472199468096),c(3153.400350143593,3139.949277478139,3129.8024500466113,190.9742550895973),c(2988.222236564643,2751.1610771803735,2650.1910453761902,51.306516292727636),c(3428.312175540727,3509.355074828508,3306.3758030098506,442.75623319279777),c(2488.0675290143954,2487.7884253657585,2479.9806961233376,387.6492342117199),c(2948.9491186507666,2832.1111130194977,2805.2892608168736,2878.865636425273))
targetgene="HNRNPA1"
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
targetmat=list(c(3213.465118717757,3131.968288029211,3121.8486954086275,0.0),c(2937.398201617274,3219.759171967416,3113.0995653068453,76.95977443909145),c(2221.241345540707,3055.578817589474,2860.96554328276,18.052292769663428),c(3116.4374156364156,2876.5766256635234,2872.100799775937,24.703137474276268),c(2473.051336870855,2727.21810883359,2663.712428260763,412.35237168599616),c(2363.317625052671,2413.6792376257154,2254.0940644045986,413.3024923580837))
targetgene="RPS16"
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
targetgenelist=c("PCBP4","MED12","ZFAND6","PCAT1","CIC","LINC01089","BIN1","MBD1","NAA25","LSM11")
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
targetmat=list(c(2605.8868827560245,2784.2251763259314,2580.198004561933,4319.248575309997),c(1156.246795052651,1161.8040354939083,1120.6840284919124,1863.186637963683),c(1156.246795052651,1223.3716683856364,1148.5221697248555,1863.186637963683),c(2119.5932756459683,2191.3516744055855,2145.127625864221,2766.751397118942),c(2418.7620268134374,2537.9546447590183,2501.455833645894,3035.6355473197186),c(3148.779983330196,3378.238819596124,3181.501855193506,3773.8793095317437),c(3512.633869885226,3597.145958766713,3563.2820778167265,4162.478664415551),c(2166.95203548329,2374.9144317309233,2371.01425758296,2837.0603268534205),c(2795.3219221053096,2850.353374617047,2790.1771270047043,3211.4078716559147),c(1293.7027077512178,1413.7752738100553,1334.6400282536756,1404.2783533453971),c(3765.598952918723,3622.229068463343,3625.3213639929995,3731.1238792878044),c(1301.7883496746629,1424.0365459586765,1337.8215301088692,1391.926784608259))
targetgene="PCBP4"
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
targetmat=list(c(2973.2060444211024,3121.7070158805896,3164.0035949899416,5829.940443929199),c(2275.530655598124,2349.8313220342934,2207.962287504293,4220.4360254128915),c(3380.953415703406,3637.050906011352,3481.358405045494,5911.650821728728),c(2004.0841053110382,1943.9410014888258,1885.039849202152,2998.580841108304),c(2530.8059220383197,2531.1137966599376,2497.478956326902,1306.4159241203797),c(3293.166446248859,3321.231752103783,3215.7030001368357,4413.310521846664),c(1900.125852009601,1943.9410014888258,2023.4351799030696,1158.197099274722),c(2704.069677540715,2662.2300518923216,2835.5135284412117,3636.1118120790493),c(2970.8958610144036,3076.101361886717,3080.4891712911117,4017.110201586156),c(2962.8102190909585,2943.844965304486,2851.4210377171794,3777.679792220094))
targetgene="MED12"
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
targetmat=list(c(2612.81743297612,2577.859592003657,2647.009543520997,4616.6363456734),c(2125.3687341627146,2279.1425583437904,1985.2571576407477,3221.859199048878),c(3110.6619571196693,3446.647300586933,3097.1920560308777,4368.654850258549),c(2530.8059220383197,2521.992665861163,2562.6997443583687,3561.052278984133),c(4722.014883291945,4595.909781232526,4713.394998469179,6179.5848512574175),c(3178.8123676172777,3065.8400897380957,3216.498375600634,4247.039404231343),c(2787.2362801818645,2944.985106654333,2828.3551492670267,3672.216397618376),c(2750.273345674687,2815.0089927717954,2678.8245620729317,3409.983092122213),c(2417.606935110088,2463.845457018975,2284.318332028937,3020.4336165663176),c(4139.848664803897,4351.919532365307,4244.918850291935,5044.190648112796),c(3361.3168567464677,3322.3718934536296,3354.893706301552,4010.4593568815435),c(4104.040822000068,4273.249779225876,4192.424069681242,3806.1834123827207),c(3533.4255205455133,3779.568574742204,3600.66472461525,3500.2445559705297))
targetgene="ZFAND6"
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
targetmat=list(c(3990.841835071837,3842.276348983779,3798.713215101046,5739.678980080882),c(3680.122166870875,3904.984123225354,3873.4785086980933,5605.711965316538),c(2391.039825933054,2540.2349274587123,2428.2812909764434,3641.8125361115744),c(2854.231598976124,2762.5624906788416,2799.7216325702852,4059.865631830096),c(2232.7922625742,2250.63902459762,2195.236280083519,3246.562336523154),c(3419.0714419139326,3724.8417899495566,3678.611520067491,2324.945284598232),c(3742.497118851737,3875.340448129337,3673.0438918209024,5005.235700557207),c(3282.770620918715,3169.592952574156,3085.261424073902,2553.924366571331),c(2333.285240765589,2198.1925225046666,2227.051298635454,2555.8246079155065),c(2276.6857473014734,2345.270756634906,2328.063982537848,2376.2518008909597))
targetgene="PCAT1"
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
targetmat=list(c(2854.231598976124,2989.4506192983586,2768.7019894821483,4636.5888797872385),c(3015.9444374450263,3073.8210791870233,3038.334271709798,4743.952515733131),c(2611.662341272771,2615.484256548602,2663.712428260763,3590.506019818847),c(2529.6508303349706,2651.9687797437005,2704.27657691448,3529.698296805244),c(2905.0556339234936,3028.2154251931506,3199.7954908608685,3923.048255049489))
targetgene="CIC"
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
targetmat=list(c(2963.9653107943077,3095.4837648341127,3189.4556098314893,5111.649215831013),c(3044.821730028759,3155.911256375994,3061.4001601599507,5057.492337522022),c(2384.1092757129586,2442.182771371886,2504.6373355010874,1577.2003156653311),c(2066.4590572919005,2165.1284233591086,2147.5137522556165,2751.5494663655413))
targetgene="LINC01089"
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
targetmat=list(c(1871.2485594258685,1937.100153389745,1931.9670015662564,3600.0072265397225),c(1965.9660791005113,1980.425524683924,2183.3056481265435,3044.186633368506),c(3384.4186908134534,3696.338256203386,3507.6057953508403,4923.525322757678),c(3309.337730095749,3314.390904004702,3215.7030001368357,4242.288800870905),c(2527.340646928272,2543.6553515082524,2613.603774041465,3248.4625778673294),c(2608.197066162723,2703.2751404868072,2608.036145794876,3172.4529241003256),c(2326.354690545493,2407.9785308764813,2227.846674099252,2809.5068273628817),c(1282.1517907177247,1286.0794426272114,1347.3660356744497,1610.4545391883953),c(910.2122622392496,1070.5927275061626,954.4505565580517,1185.750598765261),c(1363.0082099521758,1377.290750614957,1396.679314429949,1180.0498747327356),c(2741.0326120478926,2639.4272248953853,2747.2268519595923,3020.4336165663176),c(2910.83109244024,2804.747720623174,2699.5043241316894,2793.3547759373937))
targetgene="BIN1"
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
targetmat=list(c(2239.722812794296,2193.631957105279,2207.1669120404945,3212.3579923280026),c(2252.428821531138,2360.0925941829146,2265.2293208977762,3239.9114918185414),c(1189.7444544497805,1264.416756980122,1250.3302290910478,1696.9155203483622),c(1173.5731706028903,1289.499866676752,1240.7857235254671,1676.9629862345237),c(2732.9469701244475,2929.023127756477,2875.282301631131,3401.432006073425),c(2325.199598842144,2438.7623473223452,2308.179595942888,1862.2365172915956),c(2982.4467780478967,3176.433800673237,3206.9538700350536,3726.373275927366),c(3288.5460794354617,3364.557123397962,3342.1676988807776,3877.4424627892868),c(2805.7177474354535,2882.2773324127575,2936.526212343606,3318.771507601808),c(5114.746062430707,4792.014093406179,4454.102597270908,5017.587269294345))
targetgene="MBD1"
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
targetmat=list(c(3950.4136254546115,4135.292675894411,4215.489958131395,7058.4464729384),c(3673.191616650779,3809.2122498382214,3620.5491112102095,5591.460155235225),c(2303.2528564785075,2186.7911090061984,2292.272086666921,3475.5414184962533),c(2734.1020618277967,2665.650475941862,2689.1644431023105,3619.0096399814734),c(2879.6436164498086,3008.8330222457544,2996.9747475922823,2491.216402213553))
targetgene="NAA25"
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
targetmat=list(c(2168.107127186639,2110.4016385664613,2091.0420943259314,3196.205940902514),c(2312.493590105302,2300.80524399088,2323.2917297550575,3447.037798333627),c(2479.9818870909503,2753.441359880067,2619.966777751852,3753.9267754179054),c(2557.3730312153534,2463.845457018975,2385.331015931331,2337.29685333537))
targetgene="LSM11"
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
targetgenelist=c("RAN","FBL","RAD51","RPS29","SNHG16","RHEB","MYC","EIF4A1","RPS16","TMSB4X")
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
targetmat=list(c(4136.2917582370565,4150.41045811524,4267.992649416489,70.07779041099135),c(3818.4574103445802,4049.7354000545415,4122.316792827836,22.640516902012592),c(3263.914264895818,3290.883618328203,3531.8646506972327,223.17080946269556),c(3206.126201642641,3203.1988903398524,3191.6960281311763,285.70176090634936),c(2873.8448379368697,2922.8242662783377,2844.5535613667266,71.15591026346814),c(2832.7264083144164,2892.5134961095255,2774.040247805198,47.43727350897876),c(2744.933004525935,2706.318765072535,2722.8987236836497,0.0),c(2863.843057758435,2835.1395382899877,2799.6110098659724,241.49884695480097),c(3370.599920132454,3310.369113436725,3496.2205581276685,48.51539336145555),c(2698.258030359907,2533.114364107893,2518.332627197456,66.84343085356099),c(3077.214368231706,3265.9854856895354,3229.6647354335378,474.3727350897876),c(2430.432583359603,2344.7545780588443,2319.965503332057,286.77988075882615),c(1092.4166572667993,1074.9498134868109,1080.170979173309,134.76498155959877),c(1090.1940394493693,1108.5081661737104,1048.401244491741,125.06190288730765),c(864.5983309802338,921.2309076306909,914.3484615670765,119.6713036249237),c(634.5573868762384,662.5068336897566,707.4577503480854,26.95299631191975),c(665.6740363202571,671.1670537379887,673.3634009337198,31.265475721826913),c(886.8245091545328,909.3231050643718,903.5002594806874,112.12446465758617),c(690.1228323119861,649.5165036174084,684.2116030201089,24.796756606966174),c(1625.844933449977,1598.8931264048538,1694.6441402094883,2052.7401991158085))
targetgene="RAN"
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
targetmat=list(c(3351.7076686842997,3557.18538481134,3579.906688508384,118.5931837724469),c(3079.436986049136,3115.5141623515024,3300.952920572666,9.703078672291111),c(3007.2019069826642,3106.8539423032703,3159.151421872009,19.406157344582223),c(3128.334578032594,3071.1305346043127,3228.1149922783393,128.29626244473803),c(3033.873320791823,3020.251741820949,3108.009897750461,0.0),c(2941.634681368482,2907.6688811939316,2771.7156330724006,2.15623970495358),c(2664.9187630984584,2733.3819527232604,2697.3279616228756,101.34326613281827),c(2473.7736307994865,2590.4883219274307,2588.845940758985,265.2174837092904),c(1402.4718427982714,1457.0820231150528,1467.6067679729176,65.7653110010842),c(1472.4843040473133,1388.882790235225,1455.2088227313302,63.609071296130615),c(1193.54576795986,1202.6880591982347,1260.7160567539265,28.031116164396543),c(1229.1076530387386,1204.8531142102927,1275.4386167283117,23.71863675448938))
targetgene="FBL"
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
targetmat=list(c(4220.751235299394,4116.85210542834,4395.07158814276,156.32737860913457),c(3483.953428821379,3523.6270321244406,3465.2256950236997,62.53095144365383),c(3000.5340535303744,3068.9654795922547,3211.0678175711564,340.68587338266565),c(2641.5812760154445,2564.5076617827344,2619.8408038629536,80.85898893575926),c(3483.953428821379,3586.4136274741236,3639.5717999835238,673.8249077979938),c(2558.233107861823,2706.318765072535,2706.626420554066,476.5289747947412),c(3177.232170016052,3360.1653787140594,3262.209341692705,832.308526112082),c(2989.420964443225,3086.285919688719,3159.9262934496082,999.4171032459844),c(2700.480648177337,2708.483820084593,2775.5899909603963,699.6997842574368),c(3196.124421464206,3133.9171299539958,3247.48678171832,1204.2598752165745),c(2732.7086065300705,2690.0808524821,2653.9351532773194,760.074495996137),c(3756.2241114565427,3787.7637435955203,3900.70352163446,4405.197717220164),c(3083.882221683996,3356.9177961959726,3225.0155059679423,1210.7285943314353),c(2986.08703771708,2799.41613059103,3018.1247947489514,3211.7190405283577))
targetgene="RAD51"
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
targetmat=list(c(3478.3968842778045,3482.490986895338,3524.1159349212403,435.5604204006232),c(3003.8679802565193,2886.0183310733514,2938.313022256232,116.43694406749333),c(2641.5812760154445,2680.3381049278387,2731.422311037241,300.79543884102446),c(2344.8617973885516,2270.0601801428425,2529.1808292838455,20.48427719705901),c(2458.215306077477,2481.1530438185,2725.9982099940466,42.04667424659481),c(2465.994468438482,2353.4147981070764,2264.949621322512,136.92122126455234),c(1530.272367300491,1709.3109320198132,1659.0000476399243,214.54585064288122))
targetgene="RPS29"
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
targetmat=list(c(3128.334578032594,3383.980983846698,3193.2457712863747,521.8100085987664),c(2600.462846392991,2761.5276678800146,2736.8464120804356,267.37372341424395),c(2596.0176107581315,2463.832603722036,2574.898252362199,309.4203976608388),c(2422.6534209985984,2476.822933794384,2567.924408163806,478.6852144996948),c(1801.4317410269396,1821.8937926468307,1810.1000052717716,98.1089065753879),c(2658.2509096461686,2627.294257132417,2633.0136206821403,502.4038512541842),c(2228.1743619734816,2460.585021203949,2379.6306148071963,294.3267197261637),c(1519.1592782133414,1618.378621513376,1581.5128898800026,356.85767116981754),c(2178.165461081309,2149.8996269736217,2046.435836439533,678.1373872079009),c(2826.0585548621266,2802.6637131091175,2700.4274479332726,1571.89874491116))
targetgene="SNHG16"
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
targetmat=list(c(3619.5331156846037,3616.724397642936,3646.5456441819165,475.4508549422644),c(3261.6916470783885,2998.601191700369,3229.6647354335378,337.4515138252353),c(3316.145783605421,3295.213728352319,3349.7698299614167,734.1996195366941),c(2903.8501784721734,2973.7030590617014,2826.731515081945,225.32704916764914),c(3118.3327978541593,3176.1357026891274,3173.099110268795,791.339971717964),c(2291.518969770234,2411.8712834326434,2401.3270189799746,354.70143146486396),c(2261.5136292349302,2327.43413796238,2247.1275750377304,152.01489919922741),c(2301.5207499486687,2188.8706171906665,2251.001932925726,128.29626244473803),c(2544.8974009572435,2614.303927060069,2725.2233384164474,605.9033570919561),c(2929.4102833726174,2825.3967907357264,2894.920213910676,2168.099023330825),c(3538.407565348412,3474.913294353135,3435.7805750749294,3907.1063453758875),c(1908.117396263575,2000.510831141618,2035.5876343531438,2399.8947916133347))
targetgene="RHEB"
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
targetmat=list(c(3525.0718584438323,3429.4471390999165,3382.3144362205835,95.95266687043431),c(3269.470809439393,3469.50065682299,3474.5241539548906,141.2337006744595),c(2953.859079364346,2954.2175639531793,2974.7319864033952,264.1393638568136),c(2774.938345061239,2869.780418482916,2786.4381930467853,88.40582790309679),c(2639.3586581980144,2598.0660144696335,2589.620812336584,251.2019256270921),c(1945.9018991598834,2031.9041288164594,1953.4512471276269,1063.026174542115),c(1832.5483904709583,1886.8454430085715,1942.603045041238,142.3118205269363),c(2619.3550978411454,2801.581185603088,2721.348980528451,1016.667020885613),c(1962.5715327906078,2030.8216013104304,1993.744569162786,1312.0718604642536),c(984.6196931214488,974.2747554261126,1015.8566382325739,235.03012783994023),c(1024.626813835187,938.5513477271552,988.7361330166013,253.35816533204567))
targetgene="MYC"
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
targetmat=list(c(3415.052276481052,3452.180216726526,3436.5554466525286,282.467401348919),c(3305.0326945182715,3524.70955963047,3376.8903351773893,2.15623970495358),c(2970.5287129950707,3066.8004245801967,2906.543287574664,43.12479409907161),c(3171.6756254724773,3192.3736152795623,3226.565249123141,627.4657541414919),c(2744.933004525935,2640.284587204765,2665.5582269413076,51.749752918885925),c(2819.3907014098368,2894.6785511215835,2950.7109674978196,443.10725936796075),c(2254.8457757826404,2327.43413796238,2362.583440100014,808.5898893575926))
targetgene="EIF4A1"
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
targetmat=list(c(3091.6613840450004,2973.7030590617014,3041.370942076928,0.0),c(2826.0585548621266,3057.057677025936,3032.8473547233366,112.12446465758617),c(2137.0470314588556,2901.1737161577576,2787.2130646243845,19.406157344582223),c(2998.3114357129443,2731.2168977112024,2798.061266710774,7.54683896733753),c(2379.3123735587155,2589.4057944214014,2595.044913379779,294.3267197261637),c(2273.7380272307946,2291.7107302634226,2195.9860509161817,631.778233551399))
targetgene="RPS16"
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
targetmat=list(c(2903.8501784721734,2992.106026664195,3022.774024214547,1.07811985247679),c(2846.062115218996,2924.9893212903958,2953.8104538082166,107.811985247679),c(3047.2090276964022,3012.674049278746,3043.6955568097255,523.96624830372),c(2592.6836840319866,2704.153710060477,2578.772610250195,518.5756490413361),c(2794.941905418108,2682.5031599398967,2681.830530070891,459.2790571551126),c(2948.3025348207716,2888.1833860854094,2881.7473970914893,1019.9013804430434),c(3034.984629700538,3119.8442723756184,3125.0570724576432,1064.1042943945918))
targetgene="TMSB4X"
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
targetgenelist=c("ZFAND6","MED12","MBD1","ELOB","BIN1","PCBP4","COX8A","SHTN1","TERC","PCAT1")
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
targetmat=list(c(2513.780751513225,2447.594691131601,2578.772610250195,5135.084857346951),c(2044.8083920355143,2163.972484551999,1934.0794576876465,3165.359886871856),c(2992.75489116937,3272.4806507257094,3017.349923171352,5027.272872099272),c(2434.8778189944633,2394.5508433361792,2496.636223024678,3805.763079243069),c(4543.03081882673,4363.668376802955,4591.888968852962,6704.827362553157),c(3058.322116783552,2910.916463712019,3133.580659811235,4322.182488579451),c(2681.588396729183,2796.1685480729434,2755.443329942817,4704.915036208712),c(2646.026511650304,2672.760412385636,2609.767473354164,4163.698870265363),c(2325.9695459403974,2339.3419405286995,2225.431170864952,3665.6074984210863),c(3982.9311288343933,4132.007490512747,4135.489609647023,5743.1444541438605),c(3233.9089243605144,3154.4851525685467,3268.4083143134985,4850.461216293079),c(3948.48055266423,4057.313092596745,4084.3480855254743,6908.592014671271),c(3399.493951759043,3588.5786824861816,3507.8436317916567,4822.430100128682))
targetgene="ZFAND6"
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
targetmat=list(c(2860.50913103229,2963.9603115074406,3082.4391356896863,4614.352968600661),c(2189.2785501684584,2231.089189925798,2151.0434994154275,3577.2016705179894),c(3252.8011758086686,3453.2627442325547,3391.6128951517744,4951.804482425897),c(1928.1209566204443,1845.7093977794689,1836.4456389101451,3885.5439483263517),c(2434.8778189944633,2403.2110633844113,2433.0967536615426,1163.2913208224566),c(3168.3416987463324,3153.402625062518,3132.8057882336357,4012.7620909186126),c(1828.1031548360984,1845.7093977794689,1971.273293412409,1453.305561138713),c(2601.574155301706,2527.7017265777476,2762.4171741412097,4684.430759011653),c(2858.2865132148604,2920.6592112662797,3001.0776200417686,3504.9676404020447),c(2850.5073508538558,2795.086020566914,2777.914605693194,4739.414871487969))
targetgene="MED12"
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
targetmat=list(c(2154.827973998295,2082.782921599823,2150.2686278378283,3021.9699464924424),c(2167.0523719941593,2240.831937480059,2206.834253002571,2458.1132636470816),c(1144.6481759764022,1200.5230041861766,1218.0981199859696,2182.114581413023),c(1129.0898512543927,1224.3386093188149,1208.799661054779,2307.1764843003307),c(2629.35687801958,2781.0131629885373,2801.1607530211704,3211.7190405283577),c(2237.0648332432015,2315.526335396061,2248.677318192929,3688.248015323099),c(2869.39960230201,3015.921631796833,3124.2822008800445,3761.5601652915207),c(3163.8964631114727,3194.5386702916203,3256.010369071911,3531.9206367139645),c(2699.369339268622,2736.6295352413476,2860.8258644963103,3459.6866065980194),c(4920.875847789814,4549.863107839946,4339.280834555617,5432.645936630545))
targetgene="MBD1"
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
targetmat=list(c(2074.813732570818,2062.2148989852717,2180.4886193641973,1447.914961876329),c(1438.0337278771499,1500.3831233562134,1444.360620644941,2617.6750018136463),c(1534.7176029353507,1543.684223597374,1435.8370332913496,2624.1437209285073),c(2494.8885000650707,2455.172383673804,2652.385410122121,1664.617052224164),c(1965.9054595167527,1908.4959931291517,2010.0168722923697,4348.057365038894),c(2132.601795823996,2231.089189925798,2202.959895114575,1775.6633970292733),c(3066.1012791445564,3227.014495472491,3088.63810831048,5909.174911425286))
targetgene="ELOB"
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
targetmat=list(c(1800.3204321182245,1839.2142327432948,1882.1630619884988,2817.1271745218523),c(1891.4477626328508,1880.3502779723974,2127.0224805098514,3706.5760528152045),c(3256.1351025348135,3509.5541745460637,3417.1836572125485,3551.3267940585465),c(3183.9000234683417,3146.907460026344,3132.8057882336357,4758.821028832552),c(2431.5438922683184,2415.11886595073,2546.228003991028,3032.7511450172105),c(2509.3355158783647,2566.6727167947924,2540.8039029478336,2517.409855533305),c(2238.1761421519163,2286.298092733278,2170.4152888554077,2912.00172153981),c(1233.5528886735983,1221.0910268007278,1312.6324524530742,1830.6475095055896),c(875.7114200673833,1016.4933281612442,929.8458931190607,1470.5554787783417),c(1311.3445122836451,1307.693227283049,1360.6744902642256,2087.240034395066),c(2637.1360403805847,2506.0511764571675,2676.4064290276965,2858.0957289159705),c(2800.4984499616826,2663.0176648313745,2629.914134371744,3721.6697307498794))
targetgene="BIN1"
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
targetmat=list(c(2507.112898060935,2643.5321697228524,2513.683397731861,3884.465828473875),c(1112.4202176236683,1103.0955286435653,1091.7940528372972,1664.617052224164),c(1112.4202176236683,1161.552013969132,1118.9145580532697,1633.351576502337),c(2039.2518474919395,2080.617866587765,2089.828644785089,2055.974558673239),c(2327.0808548491127,2409.7062284205854,2436.9711115495384,2255.426731381445),c(3029.4280851569633,3207.5290003639684,3099.486310396869,3240.828276545231),c(3379.4903914021734,3415.374281521539,3471.4246676444936,4962.585680950665),c(2084.815512749253,2254.9047950584363,2309.8921728232667,2666.190395175102),c(2689.3675590901876,2706.318765072535,2718.249494218054,4163.698870265363),c(1244.6659777607479,1342.3341074759774,1300.2345072114867,2037.6465211811333),c(3622.867042410748,3439.189886654178,3531.8646506972327,3209.562800823404),c(1252.4451401217525,1352.0768550302384,1303.3339935218835,2007.4591653117832))
targetgene="PCBP4"
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
targetmat=list(c(2547.1200187746736,2631.624367156533,2492.761865136682,5233.193763922339),c(2400.4272428242994,2267.8951251307844,2281.2219244520957,3775.575723373719),c(2462.660541712337,2477.905461300413,2407.525991600768,2740.5806649960005),c(1895.8929982677107,1722.3012620921613,1800.0266747629819,3088.8133773460036),c(2908.295414107033,2828.6443732538137,3070.816062025698,1817.7100712758681),c(1681.4103788857246,1747.1993947308288,1745.7856643310365,2379.410514416276))
targetgene="COX8A"
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
targetmat=list(c(3325.036254875141,3450.015161714468,3414.8590424797508,5335.615149907634),c(2930.5215922813322,2819.9841532055816,2946.8366096098234,3858.590952014432),c(1829.2144637448134,1861.9473103699042,1788.4036010989935,3044.610463394455),c(2104.819073106122,1921.4863232015,1949.5768892396309,2859.1738487684474))
targetgene="SHTN1"
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
targetmat=list(c(2890.514471567594,2853.542505892481,2847.6530476771236,4121.652196018768),c(2898.293633928599,2836.222065796017,2953.8104538082166,4564.759455386729),c(1790.31865193979,1842.4618152613818,1872.8646030573082,2783.705459095072))
targetgene="TERC"
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
targetmat=list(c(3839.572279610164,3648.1176953177774,3700.786654613862,4400.885237810257),c(3540.6301831658416,3707.656708149373,3773.624582908188,3835.950435112419),c(2300.409441039954,2411.8712834326434,2365.6829264104103,3634.3420226992594),c(2746.04431343465,2622.964147108301,2727.547953149245,4776.07094647218),c(2148.160120546005,2136.9092969012736,2138.6455541738396,2696.377751044452),c(3289.474369796262,3536.6173621967887,3583.78104639638,2676.97159369987),c(3600.640864236449,3679.510992992619,3578.3569453531854,4197.120585692144),c(3158.3399185678977,3009.426466760659,3005.726849507364,3626.795183731922),c(2244.843995604206,2087.113031623939,2169.6404172778084,2672.659114289963),c(2190.3898590771732,2226.759079901682,2268.049107632909,3121.1569729203075))
targetgene="PCAT1"
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
targetgenelist=c("RAN","FBL","EIF4A1","TMSB4X","RAD51","SNHG16","MYC","RPS29","RPS12","HNRNPM")
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
targetmat=list(c(4227.8161789140095,4303.620712491453,4358.099676673188,76.99972052911397),c(3902.9490571597357,4199.229286758092,4209.34827158703,119.55219766362431),c(3336.135442630426,3412.364884187276,3606.430342461218,221.88077362994682),c(3277.068693220558,3321.443319838865,3259.079986967476,272.5384845043639),c(2937.434884113817,3030.718811613699,2904.608553570674,94.22334222641578),c(2895.406620110642,2999.289135048816,2832.6065436619483,6.07892530493005),c(2805.6705969687273,2806.221121864536,2780.385305706169,47.61824822195206),c(2927.2117928698017,2939.797247265288,2858.7171626398376,7.092079522418391),c(3445.1817492332593,3432.569676264701,3570.033722067796,49.64455665692874),c(2757.9628378299876,2626.6229700652057,2571.500353883054,92.1970337914391),c(3145.304406075468,3386.5476498661224,3297.8502999952516,388.03806529803484),c(2484.211172295792,2431.309979983434,2368.9452490848807,238.0912411097603),c(1116.5887436519267,1114.6310296045938,1102.9758440963005,236.06493267478362),c(1114.316945597701,1149.428171515714,1070.5353780934681,253.28855437208543),c(883.7294430937935,955.2376698826881,933.6524361790779,54.71032774437045),c(648.5983444814345,686.9629306324384,722.3937917216086,38.499860264556986),c(680.4035172405942,695.9428382224049,687.5796330844227,49.64455665692874),c(906.4474236360504,942.8902969464841,922.5752038854279,51.670865091905426),c(705.3932958370768,673.4930692474886,698.6568653780727,40.52616869953367),c(1661.8202766660922,1657.915438797568,1730.4219304437656,2129.650165160494))
targetgene="RAN"
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
targetmat=list(c(3425.871465772341,3688.4970425787465,3655.4866569045253,92.1970337914391),c(3147.576204129694,3230.521755490454,3370.643540782095,0.0),c(3073.7427673673587,3221.5418479004875,3225.848290086526,4.052616869953367),c(3197.555761322659,3184.4997290918755,3296.267838239016,0.0),c(3101.004344018067,3131.7427720008222,3173.627052130747,0.0),c(3006.724724767701,3015.0039733312574,2830.2328510275947,5.065771087441709),c(2723.8858670166023,2834.2833330831813,2754.27468672828,0.0),c(2528.5112343531932,2686.114857848734,2643.502363791779,354.6039761209196),c(1433.5045722164105,1510.8694520118663,1498.5912831552318,26.342009654696884),c(1505.0662109245197,1440.1526797408799,1485.9315891053461,17.22362169730181),c(1219.9555551191957,1247.0846665565998,1287.3326386977626,23.302547002231858),c(1256.3043239868066,1249.3296434540914,1302.366025382002,39.51301448204533))
targetgene="FBL"
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
targetmat=list(c(3490.617710317773,3579.6156630504024,3509.1089444527206,1.0131542174883417),c(3378.163706633601,3654.8223891163716,3448.1841668376455,0.0),c(3036.258099472635,3180.0097752968923,2967.907023820103,71.93394944167225),c(3241.85582338006,3310.2184353514067,3294.68537648278,873.3389354749505),c(2805.6705969687273,2737.7493264910413,2721.8342207254473,121.578506098601),c(2881.775831785288,3001.534111946308,3013.007183872821,596.7478341006332),c(2304.7391260119625,2413.3501648035012,2412.4629473813634,449.8404725648237))
targetgene="EIF4A1"
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
targetmat=list(c(2968.1041578458644,3102.558072333431,3086.591655537782,6.07892530493005),c(2909.0374084359964,3032.9637885111906,3016.1721073852923,0.0),c(3114.635132343421,3123.8853528596014,3107.9548892469643,253.28855437208543),c(2650.0524302542676,2803.9761449670445,2633.2163623762467,281.656872461759),c(2856.7860531888055,2781.526375992128,2738.4500691659227,273.55163872185227),c(3013.540118930378,2994.799181253833,2942.587635720331,230.9991615873419),c(3102.14024304518,3235.011709285437,3191.03413144934,864.2205475175555))
targetgene="TMSB4X"
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
targetmat=list(c(4314.144504974585,4268.823570580333,4487.8615406845165,183.38091336538983),c(3561.043449998769,3653.699900667626,3538.384486943082,81.05233739906734),c(3066.9273732046818,3182.254752194384,3278.8607589204225,365.74867251329135),c(2700.031987447233,2659.1751350788345,2675.1515989164936,275.57794715682894),c(3561.043449998769,3718.8042306948832,3716.411434519601,613.971455797935),c(2614.8395604137695,2806.221121864536,2763.769457265694,452.87993521728873),c(3247.535318515624,3484.204144907008,3331.0819968762016,695.0237931970024),c(3055.568382933553,3200.214567374317,3226.639520964644,766.9577426386746),c(2760.2346358842137,2808.4660987620277,2834.189005418184,1058.746157275317),c(3266.8456019765426,3249.6040591191327,3316.0486101919623,960.470198178948),c(2793.175707670486,2789.383795133349,2709.9657575536794,1112.4433308021992),c(3839.3387116414165,3927.587082161605,3983.0562404453203,4387.970915942008),c(3152.119800238145,3480.836679560771,3293.1029147265444,1640.2966781136251),c(3052.1606858522146,2902.755128456676,3081.844270269075,2623.0562690773168))
targetgene="RAD51"
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
targetmat=list(c(3197.555761322659,3508.898890779416,3260.662448723712,445.78785569487036),c(2658.0037234440574,2863.4680327505725,2794.6274615122907,59.776098831812156),c(2653.460127335606,2554.7837093454737,2629.2602079856574,457.94570630473044),c(2476.2598791060022,2568.2535707304232,2622.1391300825967,292.8015688541307),c(1841.2923229499218,1889.1480592392056,1848.3153312833272,52.68401930939377),c(2717.0704728539254,2724.2794651060917,2688.6025238444972,491.3797954818457),c(2277.477549361254,2551.416243999236,2429.8700266999563,126.64427718604271),c(1552.7739700632592,1678.1202308749926,1614.9022222385577,59.776098831812156),c(2226.3620931411765,2229.2620592091876,2089.640749109275,673.7475546297472),c(2888.591225947965,2906.122593802914,2757.4396102407513,1576.4679624118596))
targetgene="SNHG16"
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
targetmat=list(c(3603.0717140019447,3556.04340562674,3453.7227829844705,95.23649644390412),c(3341.81493776599,3597.575478230335,3547.879257480496,14.184159044836782),c(3019.219614065942,3063.2709766273274,3037.5353410944745,260.38063389450383),c(2836.339870700774,2975.716877625154,2845.266237711834,223.90708206492351),c(2697.760189393007,2693.9722769899545,2644.2935946698967,99.28911331385748),c(1988.9591964745916,2106.9108182958935,1994.6930437351316,195.53876397524994),c(1873.0974957090816,1956.4973661639547,1983.6158114414816,125.63112296855437),c(2677.314006904976,2905.0001053541678,2778.8028439499335,1410.3106707437717),c(2005.9976818812845,2105.788329847148,2035.8370493972607,1284.6795477752173),c(1006.4065380219807,1010.239603871233,1037.3036812125179,374.8670604706864),c(1047.298902998043,973.1974850626211,1009.6106004783927,352.5776676859429))
targetgene="MYC"
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
targetmat=list(c(3555.363954863205,3611.0453396152852,3598.5180336800395,314.0778074213859),c(3070.3350702860203,2992.5542043563414,3000.347489822935,174.26252540799476),c(2700.031987447233,2779.2813990946365,2789.0888453654657,165.1441374505997),c(2396.746947208103,2353.858277019973,2582.577586176704,140.8284362308795),c(2512.608647973613,2572.743524525407,2783.5502292186407,331.30142911868774),c(2520.5599411634034,2440.2898875734004,2312.7678567385124,446.80100991235867),c(1564.1329603343877,1772.409260569641,1694.025310050344,265.44640498194553))
targetgene="RPS29"
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
targetmat=list(c(3907.492653268187,3718.8042306948832,3827.9749883342192,338.3935086411061),c(3439.5022540976947,3522.3687521643656,3510.6914062089563,431.60369665003356),c(2743.196150477521,2825.303425493215,2940.2139430859775,13.171004827348442),c(2847.6988609719024,2924.0824089828466,2985.314103138696,92.1970337914391),c(2810.2141930771786,2880.3053594817598,2941.796404842213,408.3011496478017))
targetgene="RPS12"
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
targetmat=list(c(3370.2124134438113,3331.545715877577,3244.8378311613546,64.84186991925387),c(2863.6014473514824,2951.022131752746,2717.878066334858,167.17044588557638),c(2881.775831785288,2892.6527324179638,2975.0281017231637,106.38119283627587),c(2396.746947208103,2396.512838072314,2447.277106018549,249.23593750213206),c(2601.208772088415,2748.9742109784997,2689.393754722615,413.3669207352434),c(2255.8954678461105,2307.8362506213944,2361.032940303702,37.486706047068644),c(3473.57922491108,3330.4232274288315,3400.7103141505736,2080.0056085035653))
targetgene="HNRNPM"
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
targetgenelist=c("CIC","MED12","ATRAID","BIN1","RET","SMIM31","ZFAND6","ELOB","FARS2","TCEAL4")
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
targetmat=list(c(2806.80649599584,2943.1647126115254,2754.27468672828,4557.167670262561),c(2965.8323597916383,3026.228857818716,3022.501954410235,5392.006745472954),c(2568.2677003021427,2574.9885014228985,2649.832210816722,4406.207691856798),c(2487.6188693771305,2610.9081317827645,2690.184985600733,4172.169067616991),c(2856.7860531888055,2981.329319868883,3183.1218226681613,5368.704198470722))
targetgene="CIC"
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
targetmat=list(c(2923.804095788463,3073.37337266604,3147.5164331528576,3075.9362042946054),c(2237.721083412305,2313.4486928651236,2196.4569176551868,3838.8413300633265),c(3324.7764523592973,3580.738151499148,3463.217553521885,7218.723799604434),c(1970.7848120407862,1913.8428051116136,1875.2171811393343,3389.0008574985027),c(2488.7547684042434,2491.924356215708,2484.4649572900885,1405.2448996563298),c(3238.4481262987215,3269.8088511965575,3198.9464402305184,3712.197052877284),c(1868.55389960063,1913.8428051116136,2012.8913539318426,857.128467995137),c(2659.1396224711702,2621.010527821477,2820.7380804901804,2706.1349149113607),c(2921.5322977342375,3028.4738347162074,3064.437190950482,3041.4889609000015),c(2913.5810045444478,2898.265174661693,2836.5626980525376,5958.359953048937))
targetgene="MED12"
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
targetmat=list(c(2477.395778133115,2564.886105384186,2547.7634275395176,3572.3817708638926),c(2627.3344497120106,2755.7091416709745,2625.304053595068,2740.582158305964),c(3069.1991712589074,3112.6604683721434,3243.255369405119,5264.349314069424),c(2762.5064339384394,2777.036422197145,2874.541780202195,4113.406123002667),c(2496.7060615940336,2845.5082175706398,2784.3414600967585,3462.9611153751516),c(2332.000702662671,2582.845920564119,2366.571556450527,3088.0940549044653),c(3366.804716362473,3439.3046069571756,3199.7376711086363,4866.179706596505),c(2592.1215798715125,2704.074673028667,2618.9742065701253,3275.5275851398087))
targetgene="ATRAID"
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
targetmat=list(c(1840.156423922809,1907.1078744191386,1921.8998029482884,3276.540739357297),c(1933.3001441460624,1949.7624354714796,2171.9287604335327,2395.0965701424398),c(3328.1841494406362,3639.1075508339304,3489.328172499774,6405.160962961296),c(3254.350712678301,3263.0739205040827,3198.9464402305184,4187.366380879316),c(2485.347071322905,2504.271729151912,2599.9846654952967,2538.964469025784),c(2564.860003220804,2661.420111976326,2594.4460493484717,2645.3456618620603),c(2287.70064060527,2370.69560375116,2216.237689608133,3672.6840383952385),c(1260.8479200952581,1266.1669701852786,1340.3451075316593,883.4704776498339),c(895.0884333649219,1054.0166533723198,949.4770537414352,1210.7192898985684),c(1340.3608519931572,1355.9660460849439,1389.4014219749668,1705.138548032879),c(2695.488391338781,2598.5607588465605,2732.9114530190977,5057.665853701802),c(2862.4655483243696,2761.3215839147033,2685.437600332026,3032.3705729426065))
targetgene="BIN1"
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
targetmat=list(c(2513.744547000726,2524.476521229337,2517.696654171039,3422.434946675618),c(2308.146823093301,2237.1194783504084,2259.755387904616,2995.8970211130263),c(2723.8858670166023,2763.5665608121953,2781.176536584287,3843.907101150768),c(2797.7193037789375,2838.7732868781645,2756.6483793626335,4010.0643928188565))
targetgene="RET"
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
targetmat=list(c(2621.6549545764465,2629.9904354114433,2862.673317030427,3685.855043222587),c(2258.167265900336,2094.56344535969,2117.3338298434005,2941.186693368656),c(2221.8184970327247,2251.711828184104,2265.294004051441,3759.815301099236),c(2159.3440505415183,2234.8745014529163,2040.5844346659678,2431.57012197202),c(2531.9189314345317,2660.29762352758,2595.2372802265895,2903.6999873215873),c(2720.478169935264,2674.8899733612757,2753.483455850162,2948.2787728910744),c(3378.163706633601,3715.4367653486456,3478.250940206124,4957.363586170456))
targetgene="SMIM31"
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
targetmat=list(c(2569.4035993292555,2537.9463826142864,2633.2163623762467,2911.8052210614937),c(2090.054209887635,2243.854409042883,1974.9122717821851,2129.650165160494),c(3058.976080014892,3393.282580558597,3081.053039390957,5895.54439156466),c(2488.7547684042434,2482.9444486257416,2549.3458892957533,2683.845522126617),c(4643.555222837311,4524.750936894378,4688.834183726454,7422.367797319591),c(3125.9941226145497,3018.371438677495,3199.7376711086363,4276.52395201829),c(2740.924352423295,2899.3876631104386,2813.6170025871197,2594.687950987643),c(2704.5755835556843,2771.4239799534157,2664.8655975009615,2911.8052210614937),c(2377.436663747185,2425.697537739705,2272.4150819545016,4097.1956555228535),c(4071.0621131724365,4284.538408862773,4222.799196515033,4448.760168991308),c(3305.4661688983792,3270.931339645303,3337.4118439011445,3358.6062309738527),c(4035.8492433319384,4207.086705899313,4170.577958559254,3842.89394693328),c(3474.7151239381933,3721.049207592375,3581.902185239564,3374.816698453666))
targetgene="ZFAND6"
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
targetmat=list(c(2120.7234836196817,2138.3404948607763,2226.5236910236654,1122.5748729770826),c(1469.8533410840216,1555.7689899616987,1474.8543568116959,2671.687671516757),c(1568.676556442839,1600.6685279115313,1466.1508171523994,2642.306199209595),c(2550.0933158683374,2545.8038017555073,2708.3832957974437,2349.5046303554645),c(2009.405378962623,1978.9471351388709,2052.4528978377357,2913.8315294964705),c(2179.7902330295497,2313.4486928651236,2249.4693864890833,2299.8600736985354),c(3133.9454158043395,3346.1380657112727,3153.8462801778005,4950.271506648038))
targetgene="ELOB"
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
targetmat=list(c(2543.27792170566,2679.379927156259,2651.414672572958,3952.314602422021),c(2350.1750870964765,2408.8602110085176,2213.8639969737796,2966.5155488058645),c(2847.6988609719024,3026.228857818716,3001.138720701053,3993.8539253390427),c(2327.4571065542195,2421.2075839447216,2292.195853907448,3342.3957634940393))
targetgene="FARS2"
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
targetmat=list(c(2627.3344497120106,2753.464164773483,2626.886515351304,4891.508562033713),c(3131.673617750114,3184.4997290918755,3012.215952994703,625.1161521903068),c(2321.7776114186554,2497.5367984594373,2415.627870893835,2336.333625528116),c(2605.7523681968664,2778.158910645891,2754.27468672828,4129.61659048248),c(3446.317648260372,3547.0634980367736,3504.3615591840135,5486.23008769937))
targetgene="TCEAL4"
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



dev.off()
Sweave("gene_pooled_screen_1_counts_filtered_mageck_summary.Rnw");
library(tools);

texi2dvi("gene_pooled_screen_1_counts_filtered_mageck_summary.tex",pdf=TRUE);

