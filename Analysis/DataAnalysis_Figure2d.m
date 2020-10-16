  % Mission 1 : Data Analysis of Figure 2d

    data=importdata('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200325 - CAO Data\Volume_Consistency.csv') ; A=data.data  ; B=data.textdata ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200331 - Volume&Surface_PLUS\WorkSpace_EMBRYO' )      ; DATA=cell(size(data,1)-1,3)    ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200315 - QualityControl_Final\WorkSpace_WTinfo')      ; B=B(2:end,:) ; ValidCell={}    ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_CellName')   ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_AveCycle')   ;
    figure ; set(gcf,'position',[50,50,400,400]) ; c=hsv(size(EMBRYO,1)) ; Xrange=[] ; Yrange=[] ;
for k0=1:length(AveCycle)
for k1=1:size(AveCycle{k0,1},1)
for k2=1:size(AveCycle{k0,1},2)
if  AveCycle{k0,1}(k1,k2)~=0
    ValidCell{end+1,1}=CellName{k0,1}{k1,k2}     ;
end
end
end
end
for I=1:size(A,1)
    Data{I,1}=A(I,1) ; Data{I,2}=A(I,2) ; Data{I,3}=B{I,1} ; Data{I,4}=B{I,2}        ;
end
for I=1:size(EMBRYO,1)
    I
for i=1:size(WTinfo,1)
if  strcmp(['CD',EMBRYO{I,1},'.csv'],WTinfo{i,1})==1
    for j=1:size(Data,1)
    if  strcmp(EMBRYO{I,1},Data{j,3})==1
        for J=1:length(ValidCell)
        if  strcmp(Data{j,4},ValidCell{J,1})==1
            Data{j,2}=(Data{j,2}*WTinfo{i,2}-WTinfo{i,3})/EMBRYO{I,20}    ; Xrange=[Xrange,Data{j,2}] ;
            x=200/100*cos([-pi:pi/180:pi]) ; y=1/100*sin([-pi:pi/180:pi]) ; Yrange=[Yrange,Data{j,1}] ;
            patch(Data{j,2}+x,Data{j,1}+y,'b','facealpha',0.2,'edgecolor','none')    ; hold on        ;
        end
        end
    end
    end
end
end
end
    axis([0,200,0,1.0]) ; axis square      ;
    x=xlabel('Developmental time (\itt_{c}\rm / min)')                ; set(x,'Fontname','Arial','Fontsize',15) ;
    y=ylabel({'Volume inconsistency';'coefficient (\it\rho_{c}\rm)'}) ; set(y,'Fontname','Arial','Fontsize',15) ;
    set(gca,'xtick',[0,50,100,150,200]     ,'FontSize',15,'Fontname','Arial') ; clear x ;
    set(gca,'ytick',[0,0.2,0.4,0.6,0.8,1.0],'FontSize',15,'Fontname','Arial') ; clear y ;

    data=importdata('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200325 - CAO Data\Volume_Consistency.csv') ; A=data.data  ; B=data.textdata ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200331 - Volume&Surface_PLUS\WorkSpace_EMBRYO' )      ; DATA=cell(size(data,1)-1,3)    ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200315 - QualityControl_Final\WorkSpace_WTinfo')      ; B=B(2:end,:) ; ValidCell={}    ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_CellName') ; Point=[]       ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_AveCycle') ;
for k0=1:length(AveCycle)
for k1=1:size(AveCycle{k0,1},1)
for k2=1:size(AveCycle{k0,1},2)
if  AveCycle{k0,1}(k1,k2)~=0
    ValidCell{end+1,1}=CellName{k0,1}{k1,k2}               ;
end
end
end
end
for I=1:size(A,1)
    Data{I,1}=A(I,1) ; Data{I,2}=A(I,2) ; Data{I,3}=B{I,1} ; Data{I,4}=B{I,2} ;
end
for I=1:size(EMBRYO,1)
    I
for i=1:size(WTinfo,1)
if  strcmp(['CD',EMBRYO{I,1},'.csv'],WTinfo{i,1})==1
    for j=1:size(Data,1)
    if  strcmp(EMBRYO{I,1},Data{j,3})==1
        for J=1:length(ValidCell)
        if  strcmp(Data{j,4},ValidCell{J,1})==1
            Data{j,2}=(Data{j,2}*WTinfo{i,2}-WTinfo{i,3})/EMBRYO{I,20} ; Point=[Point;[Data{j,1},Data{j,2}]] ;
        end
        end
    end
    end 
end
end
end

    step1=15 ; figure    ; set(gcf,'position',[200,200,400,100])       ;
for I=0:step1-1
    Ratio=length(intersect(find(Point(:,1)>1/step1*I),find(Point(:,1)<1/step1*(I+1))))/size(Point,1)         ;
    patch(1-[I,I+1,I+1,I]*1/step1,[0,0,Ratio,Ratio],'k','facealpha',0.40,'edgecolor','k') ; hold on          ;
end
    axis([0,1,0,0.8])    ;
    y=ylabel('Fraction') ; set(y,'Fontname','Arial','Fontsize',15)     ;
    set(gca,'xtick',[]         ,'FontSize',15,'Fontname','Arial')      ; clear x ;
    set(gca,'ytick',[0,0.4,0.8],'FontSize',15,'Fontname','Arial')      ; clear y ;
    set(gca,'YTickLabel',{'0.0';'0.4';'0.8'}) ;

    step2=15 ; figure    ; set(gcf,'position',[200,200,400,100])       ;
for I=0:step2-1
    Ratio=length(intersect(find(Point(:,2)>200/step2*I),find(Point(:,2)<200/step2*(I+1))))/size(Point,1)     ;
    patch(    [I,I+1,I+1,I]*200/step2,[0,0,Ratio,Ratio],'k','facealpha',0.40,'edgecolor','k') ; hold on      ;
end
    axis([0,200,0,0.4])  ;
    y=ylabel('Fraction') ; set(y,'Fontname','Arial','Fontsize',15)     ;
    set(gca,'xtick',[]             ,'FontSize',15,'Fontname','Arial')  ; clear x ;
    set(gca,'ytick',[0,0.2,0.4]    ,'FontSize',15,'Fontname','Arial')  ; clear y ;
    set(gca,'YTickLabel',{'0.0';'0.2';'0.4'})                          ;