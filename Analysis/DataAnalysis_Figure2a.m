  % Mission 1 : Data Analysis of Figure 2a
  
    clear all ; clc ; figure ; set(gcf,'position',[50,50,1400,400]) ; Data=zeros(7,3,7) ; CELLNUM=[]  ; % Method ~ Embryo ~ Time Point
    File=dir('D:\STORAGE\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200530 - CAO FigureUpdate_DATA') ; File([1,2,7,11],:)=[]          ;
    File=[File([1,4,5,6,7],:);File([2,3],:)] ; X=0 ; c=[[238,0,0];[255,192,0];[0,139,69];[112,197,184];[59,73,146];[156,0,232];[256,256,256]]/256              ;
for I=1:size(Data,1)
    file=dir(['D:\STORAGE\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200530 - CAO FigureUpdate_DATA\',File(I).name,'\*.mat'])         ;
for i=1:size(Data,2)
    load(['D:\STORAGE\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200530 - CAO FigureUpdate_DATA\',File(I).name,'\',file(3*i-2).name]) ;
    CellNum=cell_num ; Dice=DICES ; CELLNUM=[CELLNUM,CellNum]    ;
for j=1:size(Data,3)
    Data(I,i,j)=Dice(j)*100       ;
end
end
end
for j=1:size(Data,3)
    X=X+2 ; XXX=X ;
for I=1:size(Data,1)
    X=X+1 ; value=[] ; VALUE=[]   ;
for i=1:size(Data,2)
    value=[value,Data(I,i,j)]     ;
end
if  I~=length(File)
    patch([X+0.05,X+0.95,X+0.95,X+0.05],[0,0,mean(value),mean(value)],c(I,:),'facealpha',1,'edgecolor','none') ; hold on ;
    plot((X+0.5)*ones(1,length(value)),value,'.','markersize',35,'color',[0.5,0.5,0.5])   ; hold on            ;
    plot((X+0.5)*ones(1,2) ,mean(value)+std(value)*[-1,1],'k-','linewidth',1.5) ; hold on ;
    plot((X+0.5)+[-0.1,0.1],mean(value)+std(value)*[1,1] ,'k-','linewidth',1.5) ; hold on ;
    plot((X+0.5)+[-0.1,0.1],mean(value)-std(value)*[1,1] ,'k-','linewidth',1.5) ; hold on ;
end
if  I==length(File)
    plot([X,X+1,X+1,X,X]+0.1,[0,0,mean(value),mean(value),0],'k-','linewidth',1.5)            ; hold on ;
    plot((X+0.5)*ones(1,length(value))+0.1,value,'.','markersize',35,'color',[0.5,0.5,0.5])   ; hold on ;
    plot((X+0.5)*ones(1,2)+0.1,mean(value)+std(value)*[-1,1],'k-','linewidth',1.5)  ; hold on ;
    plot((X+0.5)+[-0.1,0.1]+0.1,mean(value)+std(value)*[1,1],'k-','linewidth',1.5)  ; hold on ;
    plot((X+0.5)+[-0.1,0.1]+0.1,mean(value)-std(value)*[1,1],'k-','linewidth',1.5)  ; hold on ;
end
if  I~=length(File)
for i=1:size(Data,2)
    VALUE=[VALUE,Data(length(File),i,j)] ;
end
    [p,h,stats]=ranksum(value,VALUE,'alpha',0.05,'tail','left')             ;
    plot([X-1,XXX+6]+1.5,(155-8.25*I)*[1,1]   ,'k-','linewidth',1.25)       ; hold on ;
    plot((XXX+6+0.1+1.5)*[1,1],(155-8.25*I)*[1,0.99],'k-','linewidth',1.25) ; hold on ;
    plot((X-1+1.5)*[1,1],(155-8.25*I)*[1,0.99],'k-','linewidth',1.25)       ; hold on ;
    
    [p,h,stats]=ranksum(value,VALUE,'alpha',0.01,'tail','left')             ;
if  h==1
    plot([0.8,0,-0.4]+(X-1+1.5+XXX+6+1.5)/2,(158.5-8.25*I)*ones(1,3),'k*','linewidth',0.8,'markersize',6) ; hold on ;
else if h==0
    [p,h,stats]=ranksum(value,VALUE,'alpha',0.05,'tail','left')             ;
if  h==1
    plot([0,0.6]-0.3+(X-1+1.5+XXX+6+1.5)/2,(158.5-8.25*I)*ones(1,2) ,'k*','linewidth',0.8,'markersize',6) ; hold on ;
else if h==0
    [p,h,stats]=ranksum(value,VALUE,'alpha',0.10,'tail','left')             ;
if  h==1
    plot([0.8]-0.8+(X-1+1.5+XXX+6+1.5)/2,(158.5-8.25*I)*ones(1,1)   ,'k*','linewidth',0.8,'markersize',6) ; hold on ;
else if h==0
    text([0.8]-1.4+(X-1+1.5+XXX+6+1.5)/2,(158.5-8.25*I)*ones(1,1)+0.6,'n.s.','Fontname','Arial','Fontsize',11,'FontWeight','Bold') ;
    end
end
end
end
end
end

end
end
end
    CELLNUM=CELLNUM(:,1:3) ; Xaxis=cell(1,size(Data,3)) ; timepoint=[24,34,44,54,64,74,84]       ;
for t=1:length(timepoint)
    Xaxis{1,t}=['~ ',num2str(round(mean(CELLNUM(t,:)))),' / ',num2str(timepoint(t))]             ;
end
    x=xlabel('Cell number / Time point') ; set(x,'Fontname','Arial','Fontsize',18)               ;
    y=ylabel('Dice ratio (%)          ') ; set(y,'Fontname','Arial','Fontsize',18)               ;
    set(gca,'xtick',9*[1:1:7]-2.5,'FontSize',18,'Fontname','Arial')         ; clear x            ;
    set(gca,'ytick',[0,20,40,60,80,100],'FontSize',18,'Fontname','Arial')   ; clear y            ;
    set(gca,'XTickLabel',Xaxis) ; set(gca,'FontSize',18,'Fontname','Arial') ; axis([0,65,0,160]) ;