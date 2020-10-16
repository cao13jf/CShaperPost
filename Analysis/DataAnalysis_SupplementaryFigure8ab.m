  % Mission 1 : Data Analysis of Supplementary Figure 8ab
  
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200403 - Outlier\WorkSpace_Ratio') ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_CellName') ; Lineage={}         ;
  % AB ¡ú MS ¡ú E ¡ú C ¡ú D ¡ú P3 ¡ú P4
for k0=1:length(CellName)-2
    lineage=[] ;
if  size(CellName{k0,1},2)==1
    for k1=1:size(CellName{k0,1},1)
    for k2=size(CellName{k0,1},2)
    if  isnumeric(CellName{k0,1}{k1,k2})==0
        for i=1:size(Ratio,1)
        if  strcmp(CellName{k0,1}{k1,k2},Ratio{i,1})==1
            lineage=[lineage;str2num(Ratio{i,2})]     ;
        end
        end
    end
    end
    end
end
if  size(CellName{k0,1},2)~=1
    for k1=1:size(CellName{k0,1},1)
    for k2=size(CellName{k0,1},2)-1
    if  isnumeric(CellName{k0,1}{k1,k2})==0
        for i=1:size(Ratio,1)
        if  strcmp(CellName{k0,1}{k1,k2},Ratio{i,1})==1
            lineage=[lineage;str2num(Ratio{i,2})]     ;
        end
        end
    end
    end
    end
end
        Lineage{end+1,1}=lineage                      ;
end
        save('WorkSpace_Lineage','Lineage','-v7.3')   ;

    load('WorkSpace_Lineage') ; LINE=1.25 ; WIDTH=0.2 ; figure ; set(gcf,'position',[50,50,500,500])  ;
for k0=1:size(Lineage,1)-3
    data=sort(Lineage{k0,1})' ;
    Q1=data(1,round(length(data)*0.25)) ; Q2=data(1,round(length(data)*0.50)) ; Q3=data(1,round(length(data)*0.75)) ; IQR=Q3-Q1 ;
if  k0==1
    patch([k0-WIDTH,k0+WIDTH,k0+WIDTH,k0-WIDTH],[Q1,Q1,Q3,Q3],'b','facealpha',0.6,'edgecolor','none') ; hold on     ;
end
if  k0~=1
    patch([k0-WIDTH,k0+WIDTH,k0+WIDTH,k0-WIDTH],[Q1,Q1,Q3,Q3],'r','facealpha',0.6,'edgecolor','none') ; hold on     ;
end
    plot([k0-WIDTH,k0+WIDTH],[Q2,Q2],'k-','linewidth',LINE)                   ; hold on  ;
    plot([k0-WIDTH,k0+WIDTH],[Q1,Q1],'k-','linewidth',LINE)                   ; hold on  ;
    plot([k0-WIDTH,k0+WIDTH],[Q3,Q3],'k-','linewidth',LINE)                   ; hold on  ;
    plot([k0-WIDTH,k0-WIDTH],[Q1,Q3],'k-','linewidth',LINE)                   ; hold on  ;
    plot([k0+WIDTH,k0+WIDTH],[Q1,Q3],'k-','linewidth',LINE)                   ; hold on  ;
    plot([k0-WIDTH,k0+WIDTH],[Q3+1.5*IQR,Q3+1.5*IQR],'k-','linewidth',LINE)   ; hold on  ;
    plot([k0-WIDTH,k0+WIDTH],[Q1-1.5*IQR,Q1-1.5*IQR],'k-','linewidth',LINE)   ; hold on  ;
    plot([k0,k0],[Q3,Q3+1.5*IQR],'k-','linewidth',LINE)                       ; hold on  ;
    plot([k0,k0],[Q1,Q1-1.5*IQR],'k-','linewidth',LINE)                       ; hold on  ;
    outlier=data(1,find(data>Q3+1.5*IQR)) ;
    plot(k0*ones(1,length(outlier)),outlier,'.','markersize',30,'color',[0.5,0.5,0.5])   ; hold on ;
    outlier=data(1,find(data<Q1-1.5*IQR)) ;
    plot(k0*ones(1,length(outlier)),outlier,'.','markersize',30,'color',[0.5,0.5,0.5])   ; hold on ;
if  k0==1
    AB=data ;
end
if  k0~=1
    P1=data ; [p,h,stats]=ranksum(AB,P1,'alpha',0.05,'tail','left') ; Lineage{k0,2}=p    ;
end
end
    Xaxis={'AB';'MS';'E';'C'}   ;
    y=ylabel('Cell shape irregularity \it\eta') ; set(y,'Fontname','Arial','Fontsize',18)          ;
    set(gca,'xtick',[1,2,3,4],'FontSize',18,'Fontname','Arial')               ; clear x            ;
    set(gca,'ytick',[2.35,2.45,2.55,2.65]  ,'FontSize',18,'Fontname','Arial') ; clear y            ;
    set(gca,'XTickLabel',Xaxis) ;   set(gca,'FontSize',18,'Fontname','Arial') ; axis([0.5,4.5,2.35,2.75]) ;