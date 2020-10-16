  % Figure 1ABCDEF : Structure
    
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200331 - Volume&Surface_PLUS\WorkSpace_CellVolume')             ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200319 - DataImport\WorkSpace_EMBRYO') ; c=hsv(12)              ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200319 - DataImport\WorkSpace_HLH1'  ) ; C=zeros(size(c))       ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200423 - StandardDataset_CONTACT\WorkSpace_Contact12Cell_TEST') ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200424 - QueryUpdate\WorkSpace_Query') ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_CellName') ; SSS=Sequence0          ;
for k1=1:size(C,1)
for k2=1:size(C,2)
    C(k1,k2)=c(size(C,1)+1-k1,k2) ;
end
end
for n=8
    Parameter=[] ; CELLNAME={}    ;
    load(['E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200319 - DataImport\Normalization5\WorkSpace_Criteria'  ,num2str(n)]) ;
    load(['E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200319 - DataImport\Normalization5\WorkSpace_Eggshell0_',num2str(n)]) ;
for I=1:size(EMBRYO,1)
    Parameter=[Parameter;Sequence0{20,1}{1,I}] ;
end
for k0=1:length(CellName)
for k1=1:size(CellName{k0,1},1)
for k2=1:size(CellName{k0,1},2)
if  Criteria{1,n}{k0,1}{k1,k2}>1-1e-5
    CELLNAME{end+1,1}=CellName{k0,1}{k1,k2}    ;
end
end
end
end
    Sequence0{8,1}=CELLNAME                    ;
    
for I=1:size(Sequence0{8,1},1)
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  isnumeric(CellName{k0,1}{k1,k2})==0 && strcmp(CellName{k0,1}{k1,k2},Sequence0{8,1}{I,1})==1
        Sequence0{8,1}{I,2}=(3/4/pi*mean(CellVolume{k0,1}{k1,k2}))^(1/3) ; Sequence0{8,1}{I,3}=[k0,k1,k2]          ;
    end
    end
    end
    end
end


    figure('Position',[0 0 800 450])  ; Sequence0{8,1}=SSS{8,n} ; Sequence0{10,1}=SSS{10,n} ;
for j1=1:length(Sequence0{8,1})
for j2=1:length(Sequence0{8,1})
if  Sequence0{10,1}(j1,j2)<17
    Sequence0{10,1}(j1,j2)=0 ;
end
if  Sequence0{10,1}(j1,j2)==17
    Sequence0{10,1}(j1,j2)=1 ;
end
end
end
for j1=1:length(Sequence0{8,1})
for j2=1:length(Sequence0{8,1})
if  j1~=j2
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  strcmp(CellName{k0,1}{k1,k2},Sequence0{8,1}{j1,1})==1
        g1=Criteria{4,n}{k0,1}{k1,k2} ;
    end
    if  strcmp(CellName{k0,1}{k1,k2},Sequence0{8,1}{j2,1})==1
        g2=Criteria{4,n}{k0,1}{k1,k2} ;
    end
    end
    end
    end
    if  Sequence0{10,1}(j1,j2)==1 || isnan(Sequence0{10,1}(j1,j2))==1
        plot3([g1(2),g2(2)],[g1(3),g2(3)],[g1(1),g2(1)],'-','linewidth',2,'color',[0.5,0.5,0.5]) ; hold on ;
    end
end
end
end
for j =1:length(CELLNAME)
for k0=1:length(CellName)
for k1=1:size(CellName{k0,1},1)
for k2=1:size(CellName{k0,1},2)
if  isnumeric(CellName{k0,1}{k1,k2})==0 && strcmp(CellName{k0,1}{k1,k2},CELLNAME{j,1})==1
    g=Criteria{4,n}{k0,1}{k1,k2}  ;
    [rx,ry,rz]=ellipsoid(Criteria{4,n}{k0,1}{k1,k2}(2),Criteria{4,n}{k0,1}{k1,k2}(3),Criteria{4,n}{k0,1}{k1,k2}(1),Criteria{3,n}{k0,1}{k1,k2},Criteria{3,n}{k0,1}{k1,k2},Criteria{3,n}{k0,1}{k1,k2},20) ;
    surf(rx,ry,rz,'FaceAlpha',1,'FaceColor',C(j,:),'EdgeColor','none') ; hold on ;
    view(-40,20) ; axis equal     ; grid off ;
end
end
end
end
end
    axis([-25,25,-8,8,-15,15]) ;
    set(gca,'FontSize',20,'Fontname','arial') ; h=findobj(gca) ; view(-42.5,15)  ; camlight ;    
    set(gcf,'position',[0,0,800,450]) ; set(gcf, 'PaperPositionMode', 'manual')  ;
    set(gcf,'PaperUnits','points')    ; set(gcf, 'PaperPosition', [0,0,800,450]) ;
  % print(gcf,'-r600','-dpng',['Connection.png'])              ; clf             ;
end