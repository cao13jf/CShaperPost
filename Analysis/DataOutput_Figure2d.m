  % Mission 1 : Data Output of Figure 2d

    data=importdata('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200325 - CAO Data\Volume_Consistency.csv') ; A=data.data  ; B=data.textdata ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200331 - Volume&Surface_PLUS\WorkSpace_EMBRYO' )      ; DATA=cell(size(data,1)-1,3)    ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200315 - QualityControl_Final\WorkSpace_WTinfo')      ; B=B(2:end,:) ; ValidCell={}    ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_CellName')   ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_AveCycle')   ;
    figure ; set(gcf,'position',[50,50,400,400]) ; c=hsv(size(EMBRYO,1)) ; Xrange=[] ; Yrange=[] ;
    SourceData=cell(1,2)                         ;
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
    Data{I,1}=A(I,1) ; Data{I,2}=A(I,2) ; Data{I,3}=B{I,1} ; Data{I,4}=B{I,2} ;
end

for I=1:size(EMBRYO,1)
    SourceData{2,4*I-2}='Cell Lifespan (Averaged & Normalized / min)'         ; SourceData{2,4*I-3}='Cell Name' ;
    SourceData{2,4*I-1}='Volume Inconsistency' ; NUM=3     ;
if  I+3<10
    III=['0',num2str(I+3)] ;
end
if  I+3>=10
    III=[    num2str(I+3)] ;
end
    SourceData{1,4*I-2}=['Sample ',III]        ; I
for i=1:size(WTinfo,1)
if  strcmp(['CD',EMBRYO{I,1},'.csv'],WTinfo{i,1})==1
    for j=1:size(Data,1)
    if  strcmp(EMBRYO{I,1},Data{j,3})==1
        for J=1:length(ValidCell)
        if  strcmp(Data{j,4},ValidCell{J,1})==1
            Data{j,2}=(Data{j,2}*WTinfo{i,2}-WTinfo{i,3})/EMBRYO{I,20}     ; Xrange=[Xrange,Data{j,2}] ; NUM    ;
            x=200/100*cos([-pi:pi/180:pi])  ; y=1/100*sin([-pi:pi/180:pi]) ; Yrange=[Yrange,Data{j,1}] ;
            patch(Data{j,2}+x,Data{j,1}+y,'b','facealpha',0.2,'edgecolor','none') ; hold on            ;
            SourceData{NUM,4*I-3}=Data{j,4} ;
            SourceData{NUM,4*I-2}=sprintf('%2.4f',Data{j,2}) ;
            SourceData{NUM,4*I-1}=sprintf('%2.4f',Data{j,1}) ; NUM=NUM+1   ;
        end
        end
    end
    end
end
end
end
    axis([0,200,0,1.0]) ; axis square ;
    x=xlabel('Developmental time (\itt_{c}\rm / min)')                        ; set(x,'Fontname','Arial','Fontsize',15) ;
    y=ylabel({'Volume inconsistency';'coefficient (\it\rho_{c}\rm)'})         ; set(y,'Fontname','Arial','Fontsize',15) ;
    set(gca,'xtick',[0,50,100,150,200]     ,'FontSize',15,'Fontname','Arial') ; clear x ;
    set(gca,'ytick',[0,0.2,0.4,0.6,0.8,1.0],'FontSize',15,'Fontname','Arial') ; clear y ;
    xlswrite('Figure 2d.xls',SourceData,'Sheet1')                             ;