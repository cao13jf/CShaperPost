  % Mission 1 : Data Output of Figure 2e
  
    SourceData={}     ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200324 - LostRatio\WorkSpace_CellCycle') ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200324 - LostRatio\WorkSpace_EMBRYO')    ;
    Xmax=[] ; Ymax=[] ; Threshold=[] ; c=hsv(size(EMBRYO,1))    ;
    figure ; set(gcf,'position',[0,0,550,550]) ; load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200315 - QualityControl_Final\WorkSpace_WTinfo') ;
for I=1:size(EMBRYO,1)
    SourceData{2,5*I-4}='Developmental Time (Normalized / min)' ;
    SourceData{2,5*I-3}='Total Cell Number'    ; SourceData{2,5*I-2}='Lost Cell Number' ;
    SourceData{2,5*I-1}='Cell Loss Ratio (%)'  ;
if  I+3<10
    III=['0',num2str(I+3)] ;
end
if  I+3>=10
    III=[    num2str(I+3)] ;
end
    SourceData{1,5*I-4}=['Sample ',III]        ; I
for i=1:size(WTinfo,1)
if  strcmp(['CD',EMBRYO{I,1},'.csv'],WTinfo{i,1})==1
    Data=cell2mat(table2cell(readtable(['E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200318 - DataUpdate4\RevisionData\LostNucleusRatio\nucelus_lost_',EMBRYO{I,1},'.csv']))) ;
    N=intersect(find(Data(:,1)*WTinfo{i,2}>=WTinfo{i,3}),find(Data(:,1)*WTinfo{i,2}<=WTinfo{i,4})) ; Data=Data(N,:) ;
    plot((Data(:,1)*WTinfo{i,2}-WTinfo{i,3})/EMBRYO{I,20},Data(:,2),'-' ,'linewidth',1.5,'color',c(I,:)) ; hold on  ;
    plot((Data(:,1)*WTinfo{i,2}-WTinfo{i,3})/EMBRYO{I,20},Data(:,3),'.','markersize',7.5,'color',c(I,:)) ; hold on  ;
    Xmax=[Xmax,max((Data(:,1)*WTinfo{i,2}-WTinfo{i,3})/EMBRYO{I,20})] ; Ymax=[Ymax,max(Data(:,2))]       ;
    for t=1:size(Data,1)
        SourceData{2+t,5*I-4}=sprintf('%2.4f',(Data(t,1)*WTinfo{i,2}-WTinfo{i,3})/EMBRYO{I,20})          ;
        SourceData{2+t,5*I-3}=Data(t,2) ; SourceData{2+t,5*I-2}=Data(t,3)                                ;
        SourceData{2+t,5*I-1}=Data(t,3)./Data(t,2)*100         ;
    end
end
end
end
    x=xlabel('Developmental time (min)')                       ; set(x,'Fontname','Arial','Fontsize',18) ;
    y=ylabel('Cell number')                                    ; set(y,'Fontname','Arial','Fontsize',18) ;
    set(gca,'xtick',[0,50,100,150,200] ,'FontSize',18,'Fontname','Arial') ; clear x ;
    set(gca,'ytick',[0,100,200,300,400],'FontSize',18,'Fontname','Arial') ; clear y ;
    axis([0,210,0,400]) ; plot(0,0,'k.','markersize',40)       ; hold on  ;
    xlswrite('Figure 2e.xls',SourceData,'Sheet1')              ;