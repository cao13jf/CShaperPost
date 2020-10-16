  % Mission 1 : Data Analysis of Figure 2e

    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200315 - QualityControl_Final\WorkSpace_WTinfo') ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_CellName') ; CellCycle=CellName               ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_AveCycle') ; load('WorkSpace_EMBRYO')         ;
for k0=1:length(CellName)
for k1=1:size(CellName{k0,1},1)
for k2=1:size(CellName{k0,1},2)
if  isnumeric(CellName{k0,1}{k1,k2})==0
    CellCycle{k0,1}{k1,k2}=[] ;
end
end
end
end
for I=1:size(EMBRYO,1)
for i=1:size(WTinfo,1)
if  strcmp(['CD',EMBRYO{I,1},'.csv'],WTinfo{i,1})==1
    for k0=1:length(AveCycle)
    for k1=1:size(AveCycle{k0,1},1)
    for k2=1:size(AveCycle{k0,1},2)
    if  isnumeric(CellName{k0,1}{k1,k2})==0 && AveCycle{k0,1}(k1,k2)~=0
        CellCycle{k0,1}{k1,k2}=[CellCycle{k0,1}{k1,k2},EMBRYO{I,3}{k0,1}(k1,k2)*WTinfo{i,2}]     ;
    end
    end
    end
    end
end
end
end
        save('WorkSpace_CellCycle','CellCycle','-v7.3')             ;
for I =1:size(EMBRYO,1)
    cycle=[] ; cycle_ave=[] ;
for k0=1:length(AveCycle)
for k1=1:size(AveCycle{k0,1},1)
for k2=1:size(AveCycle{k0,1},2)
if  AveCycle{k0,1}(k1,k2)~=0
    cycle=[cycle,CellCycle{k0,1}{k1,k2}(I)] ; cycle_ave=[cycle_ave,mean(CellCycle{k0,1}{k1,k2})] ;
end
end
end
end 
    b =sum(cycle.*cycle_ave)/sum(cycle_ave.^2)     ;
    yy=sum(cycle.*cycle_ave)/sum(cycle_ave.^2)*cycle_ave ; y=cycle  ;
    SSE=sum((yy-y).^2)   ; SST=sum((y-mean(y)).^2) ; EMBRYO{I,20}=b ; EMBRYO{I,21}=1-SSE/SST     ; clear cycle_i ;
    clear k0 ; clear k1  ; clear k2  ; clear nn    ; clear b ; clear yy ; clear y  ; clear SSE   ; clear cycle_j ;
end
    save('WorkSpace_EMBRYO','EMBRYO','-v7.3')      ;
  
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200324 - LostRatio\WorkSpace_CellCycle') ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200324 - LostRatio\WorkSpace_EMBRYO')    ;
    Xmax=[] ; Ymax=[] ; Threshold=[] ; c=hsv(size(EMBRYO,1)) ;
    figure ; set(gcf,'position',[0,0,550,550]) ; load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200315 - QualityControl_Final\WorkSpace_WTinfo')                           ;
for I=1:size(EMBRYO,1)
for i=1:size(WTinfo,1)
if  strcmp(['CD',EMBRYO{I,1},'.csv'],WTinfo{i,1})==1
    Data=cell2mat(table2cell(readtable(['E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200318 - DataUpdate4\RevisionData\LostNucleusRatio\nucelus_lost_',EMBRYO{I,1},'.csv']))) ;
    N=intersect(find(Data(:,1)*WTinfo{i,2}>=WTinfo{i,3}),find(Data(:,1)*WTinfo{i,2}<=WTinfo{i,4})) ; Data=Data(N,:)       ;
    plot((Data(:,1)*WTinfo{i,2}-WTinfo{i,3})/EMBRYO{I,20},Data(:,2),'-' ,'linewidth',1.5,'color',c(I,:)) ; hold on        ;
    plot((Data(:,1)*WTinfo{i,2}-WTinfo{i,3})/EMBRYO{I,20},Data(:,3),'.','markersize',7.5,'color',c(I,:)) ; hold on        ;
    Xmax=[Xmax,max((Data(:,1)*WTinfo{i,2}-WTinfo{i,3})/EMBRYO{I,20})] ; Ymax=[Ymax,max(Data(:,2))]       ;
end
end
end
    x=xlabel('Developmental time (min)')                       ; set(x,'Fontname','Arial','Fontsize',18) ;
    y=ylabel('Cell number')                                    ; set(y,'Fontname','Arial','Fontsize',18) ;
    set(gca,'xtick',[0,50,100,150,200         ],'FontSize',18,'Fontname','Arial') ; clear x ;
    set(gca,'ytick',[0,100,200,300,400],'FontSize',18,'Fontname','Arial')         ; clear y ;
    axis([0,210,0,400]) ; plot(0,0,'k.','markersize',40)       ; hold on          ;

    figure ; set(gcf,'position',[100,100,235,250]) ; LostNum200=[] ; LostNum350=[]          ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200324 - LostRatio\WorkSpace_EMBRYO')            ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200315 - QualityControl_Final\WorkSpace_WTinfo') ;
for I=1:size(EMBRYO,1)
for i=1:size(WTinfo,1)
if  strcmp(['CD',EMBRYO{I,1},'.csv'],WTinfo{i,1})==1
    Data=cell2mat(table2cell(readtable(['E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200318 - DataUpdate4\RevisionData\LostNucleusRatio\nucelus_lost_',EMBRYO{I,1},'.csv']))) ;
    N=intersect(find(Data(:,1)*WTinfo{i,2}>=WTinfo{i,3}),find(Data(:,1)*WTinfo{i,2}<=WTinfo{i,4})) ; Data=Data(N,:)               ;
    plot(Data(:,2),Data(:,3)./Data(:,2)*100,'.','markersize',7.5,'color',c(I,:)) ; hold on ;
    Xmax=[Xmax,max((Data(:,1)*WTinfo{i,2}-WTinfo{i,3})/EMBRYO{I,20})] ; Ymax=[Ymax,max(Data(:,2))] ;
    LostNum200=[LostNum200;Data(find(Data(:,2)<=200),2:3)]            ;
    LostNum350=[LostNum350;Data(find(Data(:,2)<=350),2:3)]            ;
end
end
end
    axis([-50,400,0,30]) ; plot(0,0,'k.','markersize',30) ; hold on              ;
    x=xlabel('Cell number')         ; set(x,'Fontname','Arial','Fontsize',18)    ;
    y=ylabel('Cell loss ratio (%)') ; set(y,'Fontname','Arial','Fontsize',18)    ;
    set(gca,'xtick',[0,200,400],'FontSize',18,'Fontname','Arial') ; clear x      ;
    set(gca,'ytick',[0,15,30]  ,'FontSize',18,'Fontname','Arial') ; clear y      ;