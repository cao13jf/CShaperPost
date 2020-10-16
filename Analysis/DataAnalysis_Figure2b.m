  % Mission 1 : Data Analysis of Figure 2b
  
    clear all ; clc ; figure ; set(gcf,'position',[50,50,800,400]) ; X=0     ; % Method ~ Embryo ~ Time Point
    File=dir('D:\STORAGE\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200530 - CAO FigureUpdate_DATA') ; File([1,2,7,11],:)=[]            ;
    File=[File([1,4,5,6,7],:);File([2,3],:)] ; X=0 ; c=[[238,0,0];[255,192,0];[0,139,69];[112,197,184];[59,73,146];[156,0,232];[0,0,0]]/256   ;
    Embryo=dir('E:\Project 12 - Wnt & Membrane\CJF\Paper Preparation\NC Figure\Figure - Comparison Between Algorithms - B\HDistance\CShaper') ; Embryo(1:2,:)=[] ;
for i=[1,3,2]
    X=X+5     ;
for I=length(File)
    load(['D:\STORAGE\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200530 - CAO FigureUpdate_DATA\',File(I).name,'\',Embryo(i).name])     ;
    CellNum=cell_nums' ; Distance=HDistances ; clear cell_nums     ; clear HDistances ; RRR=sort(Distance) ; XXX=X+4.5 ;
end
for I=1:length(File)    
    load(['D:\STORAGE\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200530 - CAO FigureUpdate_DATA\',File(I).name,'\',Embryo(i).name])     ;
    CellNum=cell_nums' ; Distance=HDistances ; clear cell_nums     ; clear HDistances ; R=sort(Distance)   ; X=X+4.5   ;
    plot(X,mean(R),'.','markersize',25,'color',c(I,:)) ; hold on   ; length(Distance)
    plot([X,X],mean(R)+std(R)*[-1,1],'-','linewidth',2,'color',c(I,:)) ; hold on      ;
    plot([X-0.5,X+0.5],(mean(R)+std(R)*[1])*[1,1],'-','linewidth',2,'color',c(I,:))   ; hold on ;
    plot([X-0.5,X+0.5],(mean(R)-std(R)*[1])*[1,1],'-','linewidth',2,'color',c(I,:))   ; hold on ;
if  I~=length(File)
    plot([XXX+27,X],[57.5,57.5]-I*6+45,'k-','linewidth',1.5)         ; hold on ;
    plot((XXX+27)*[1,1],(57.5-I*6)*[1,0.98]+45,'k-','linewidth',1.5) ; hold on ;
    plot(X*[1,1],(57.5-I*6)*[1,0.98]+45  ,'k-','linewidth',1.5)      ; hold on ;
    [p,h,stats]=ranksum(R,RRR,'alpha',0.01,'tail','right')           ;
if  h==1
    plot([0,1.6,3.2]-1.6+(XXX+27+X)/2,([57.5]-I*6+47.5)*ones(1,3),'k*','linewidth',1,'markersize',6.5) ; hold on ;
else if h==0
    [p,h,stats]=ranksum(R,RRR,'alpha',0.05,'tail','right')           ;
if  h==1
    plot([0,1.6]-0.8+(XXX+27+X)/2,([57.5]-I*6+47.5)*ones(1,2)   ,'k*','linewidth',1,'markersize',6.5) ; hold on ;
else if h==0
    [p,h,stats]=ranksum(R,RRR,'alpha',0.10,'tail','right')           ;
if  h==1
    plot(0+(XXX+27+X)/2,([57.5]-I*6+47.5)*ones(1,1)             ,'k*','linewidth',1,'markersize',6.5) ; hold on ;
end
end
end
end
end
end

end
end
    Xaxis={'Sample 02';'Sample 03';'Sample 04'} ;
    y=ylabel('Hausdorff distance (\mum)')       ; set(y,'Fontname','Arial','Fontsize',18) ;
    set(gca,'xtick',36.5*[1,2,3]-14,'FontSize',18,'Fontname','Arial')       ; clear x     ;
    set(gca,'ytick',[0,20,40,60,80,100],'FontSize',18,'Fontname','Arial')   ; clear y     ;
    set(gca,'XTickLabel',Xaxis) ; set(gca,'FontSize',18,'Fontname','Arial') ; axis([0,120,-15,105]) ;