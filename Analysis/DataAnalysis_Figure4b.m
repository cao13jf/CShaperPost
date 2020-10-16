  % Mission 1 : Data Analysis of Figure 4b
  
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200331 - Volume&Surface_PLUS\WorkSpace_EMBRYO' )           ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200315 - QualityControl_Final\WorkSpace_WTinfo')           ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_CellName') ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_AveCycle') ; AVECYCLE=AveCycle ; Resolution=WTinfo{1,2} ;
for I=1:size(EMBRYO,1)
    EMBRYO{I,14}=[] ; EMBRYO{I,15}=[]                 ; EMBRYO{I,16}=[]      ;
end
for k0=1:length(AVECYCLE)
for k1=1:size(AVECYCLE{k0,1},1)
for k2=1:size(AVECYCLE{k0,1},2)
if  AVECYCLE{k0,1}(k1,k2)~=0
    avecycle=[]     ; AveCycle{k0,1}(k1,k2)=0         ;
    for I=1:size(EMBRYO,1)
        avecycle=[avecycle,EMBRYO{I,3}{k0,1}(k1,k2)*Resolution/EMBRYO{I,20}] ;
    end
        AveCycle{k0,1}(k1,k2)=mean(avecycle)          ;
end
end
end
end
        save('WorkSpace_AveCycle','AveCycle','-v7.3') ;
  
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200403 - Outlier\ServerComputation\Tool\WorkSpace_EMBRYO'  ) ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200403 - Outlier\ServerComputation\Tool\WorkSpace_CellName') ;
    N=cell(7,7) ; MORPHOLOGY=cell(length(CellName),size(EMBRYO,1)) ; load('WorkSpace_AveCycle') ; N(:,2)=AveCycle(:,1) ; Tree=cell(4,5)       ;
for I=1:size(EMBRYO,1)
    load(['E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200403 - Outlier\ServerComputation\Result\WorkSpace_Morphology_',num2str(I)]) ;
    MORPHOLOGY(:,I)=Morphology              ;
end
for k0=1:length(AveCycle)
    Morphology=cell(size(CellName{k0,1}))   ;
for k1=1:size(AveCycle{k0,1},1)
for k2=1:size(AveCycle{k0,1},2)
if  AveCycle{k0,1}(k1,k2)~=0
    t=linspace(0,AveCycle{k0,1}(k1,k2),501) ; morphology=cell(1,500) ; morphology_plus=[] ; k2
    for i=1:size(MORPHOLOGY,2)
        L=size(MORPHOLOGY{k0,i}{k1,k2},2)   ; seg=linspace(0,AveCycle{k0,1}(k1,k2),L+1)   ;
    for j=1:L
    for k=1:500
    if  t(1,k)>=seg(1,j) && t(1,k)<seg(1,j+1)
    if  isnan(MORPHOLOGY{k0,i}{k1,k2}(1,j))==0 && isnan(MORPHOLOGY{k0,i}{k1,k2}(2,j))==0
        morphology{1,k}=[morphology{1,k};(MORPHOLOGY{k0,i}{k1,k2}(2,j).^(1/2))/(MORPHOLOGY{k0,i}{k1,k2}(1,j).^(1/3))] ;
    end
    end
    end
    end
    end
    for k=1:500
    if  isempty(morphology{1,k})==0
        morphology_plus=[morphology_plus,mean(morphology{1,k})]  ;
    end
    end 
        Morphology{k1,k2}=morphology_plus   ;
end
end
end
        N{k0,4}=Morphology ; clear Express  ; clear L ; clear k1 ; clear k2 ; clear express ;
end
        save('WorkSpace_N','N','-v7.3')     ;

   load('WorkSpace_N')    ; Tree=cell(4,5)         ;
   load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200403 - Outlier\WorkSpace_EMBRYO')              ;
   load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200315 - QualityControl_Final\WorkSpace_WTinfo') ;
   load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200331 - Volume&Surface_PLUS\WorkSpace_EMBRYO' ) ;
  % ABcell: 64¡Á6
    CycleABa=zeros(128,7) ; CycleABp=zeros(128,7)  ;
    ExpreABa=cell( 128,7) ; ExpreABp=cell( 128,7)  ;
for k2=1:6
    CycleABa(1:2^(k2),k2) =N{1,2}(1:2^(k2),k2+1)                    ;
    CycleABp(1:2^(k2),k2) =N{1,2}(2^(k2)+1:2^(k2+1),k2+1)           ;
    ExpreABa(1:2^(k2),k2) =N{1,4}(1:2^(k2),k2+1)                    ;
    ExpreABp(1:2^(k2),k2) =N{1,4}(2^(k2)+1:2^(k2+1),k2+1)           ;
end
    Tree{1,1}=CycleABa    ; Tree{2,1}=CycleABp     ; clear CycleABa ; clear CycleABp ;
    Tree{1,3}=ExpreABa    ; Tree{2,3}=ExpreABp     ; clear ExpreABa ; clear ExpreABp ;
  % EMScell: 32¡Á5
    CycleEMS =zeros(64,6) ; ExpreEMS=cell(64,6)    ;
for k2=1:5
    CycleEMS( 1:2^(k2-1),k2)=N{2,2}(1:2^(k2-1),k2) ;
    ExpreEMS( 1:2^(k2-1),k2)=N{2,4}(1:2^(k2-1),k2) ;
end
for k2=1:4
    CycleEMS( 2^(k2-1)+1:2^(k2),k2)=N{3,2}(1:2^(k2-1),k2)     ;
    ExpreEMS( 2^(k2-1)+1:2^(k2),k2)=N{3,4}(1:2^(k2-1),k2)     ;
end
    Tree{3,1}=CycleEMS  ; Tree{3,3}=ExpreEMS ; clear CycleEMS ; clear ExpreEMS ;
  % P2cell: 16¡Á4
    CycleP2=zeros(16,4) ; Expre=cell(16,4)   ;
for k2=1:4
    CycleP2( 1:2^(k2-1),k2)=N{4,2}(1:2^(k2-1),k2)  ;
    ExpreP2( 1:2^(k2-1),k2)=N{4,4}(1:2^(k2-1),k2)  ;
end
for k2=2:4
    CycleP2( 2^(k2-1)+1:2^(k2-1)+2^(k2-2),k2)=N{5,2}(1:2^(k2-2),k2-1) ;
    ExpreP2( 2^(k2-1)+1:2^(k2-1)+2^(k2-2),k2)=N{5,4}(1:2^(k2-2),k2-1) ;
end
    CycleP2(2,1)=N{6,2} ; CycleP2(4,2)=N{7,2}(1,1) ;
    ExpreP2(2,1)=N{6,4} ; ExpreP2{4,2}=N{7,4}{1,1} ;
    Tree{4,3}=ExpreP2   ; clear ExpreP2            ; Tree{4,1}=CycleP2        ; clear CycleP2     ;
  % Initial DivTime
    T_ABa=[] ; T_ABp=[] ; T_EMS=[] ; T_P2 =[]      ;
for I=1:size(EMBRYO)
for i=1:size(WTinfo,1)
if  strcmp(['CD',EMBRYO{I,1},'.csv'],WTinfo{i,1})==1
    T_ABa=[T_ABa,(EMBRYO{I,2}{1,1}(1,1)-WTinfo{i,3}/WTinfo{i,2})*WTinfo{i,2}/EMBRYO{I,20}]        ;
    T_ABp=[T_ABp,(EMBRYO{I,2}{1,1}(2,1)-WTinfo{i,3}/WTinfo{i,2})*WTinfo{i,2}/EMBRYO{I,20}]        ;
    T_EMS=[T_EMS,(EMBRYO{I,2}{8,1}(1,1)-WTinfo{i,3}/WTinfo{i,2})*WTinfo{i,2}/EMBRYO{I,20}]        ;
    T_P2 =[T_P2 ,(EMBRYO{I,2}{9,1}(1,1)-WTinfo{i,3}/WTinfo{i,2})*WTinfo{i,2}/EMBRYO{I,20}]        ;
end
end
end
    Tree{1,5}=mean(T_ABa) ; Tree{2,5}=mean(T_ABp) ; Tree{3,5}=mean(T_EMS) ; Tree{4,5}=mean(T_P2 ) ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200403 - Outlier\WorkSpace_EMBRYO')                  ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200331 - Volume&Surface_PLUS\WorkSpace_CellVolume' ) ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200331 - Volume&Surface_PLUS\WorkSpace_CellSurface') ;
    Tree{1,7}=(mean(CellSurface{1,1}{1,1}./cell2mat(EMBRYO(:,10))')^(1/2))/(mean(CellVolume{1,1}{1,1}./cell2mat(EMBRYO(:,7))')^(1/3)) ;
    Tree{2,7}=(mean(CellSurface{1,1}{2,1}./cell2mat(EMBRYO(:,10))')^(1/2))/(mean(CellVolume{1,1}{2,1}./cell2mat(EMBRYO(:,7))')^(1/3)) ;
    Tree{3,7}=(mean(CellSurface{8,1}{1,1}./cell2mat(EMBRYO(:,10))')^(1/2))/(mean(CellVolume{8,1}{1,1}./cell2mat(EMBRYO(:,7))')^(1/3)) ;
    Tree{4,7}=(mean(CellSurface{9,1}{1,1}./cell2mat(EMBRYO(:,10))')^(1/2))/(mean(CellVolume{9,1}{1,1}./cell2mat(EMBRYO(:,7))')^(1/3)) ;    
    save('WorkSpace_Tree','Tree','-v7.3')            ; load('WorkSpace_Tree') ; Branch=cell(4,1)  ;
    
for n=1:4
    [R1,R2]=size(Tree{n,1}) ; branch=cell(R2,1)      ; clear R1 ; 
for i=1:R2
    i0=R2+1-i ; branch{i,1}=[(2^(i0-1)+1)/2:2^(i0-1) : (2^(i0-1)+1)/2+(2^i-1)*2^(i0-1)]           ;
end
if  n>1
    branch0=Branch{n-1,1}   ; Delta=branch0{end,1}(1,end)       ; clear branch0                   ;
for i=1:R2
    branch0=branch{i,1}     ; branch{i,1}=branch0+Delta         ; clear branch0                   ;
end
end
    Branch{n,1}=branch      ; clear Delta ;   clear branch      ;
end
    clear R2 ; clear n      ; clear i     ;   clear i0          ;
    save('WorkSpace_Branch','Branch','-v7.3')                   ;
  
    figure   ; MIN=[0,0,0]/256            ; MAX=[0,256,0]/256   ;
    c=[linspace(MIN(1),MAX(1),1000)',linspace(MIN(2),MAX(2),1000)',linspace(MIN(3),MAX(3),1000)'] ;
for i=1:size(c,1)
    patch([0,1,1,0],[i,i,i+1,i+1],c(i,:),'edgecolor','none') ; hold on  ;
end
    set(gca,'xtick',[],'xcolor','w') ; set(gca,'ytick',[],'ycolor','w') ;
    
    figure ; MIN=[0,256,0]/256       ; MAX=[0,0,0]/256       ;
    c=[linspace(MIN(1),MAX(1),1000)',linspace(MIN(2),MAX(2),1000)',linspace(MIN(3),MAX(3),1000)'] ;
for i=1:size(c,1)
    patch([0,1,1,0],[i,i,i+1,i+1],c(i,:),'edgecolor','none') ; hold on  ;
end
    set(gca,'xtick',[],'xcolor','w') ; set(gca,'ytick',[],'ycolor','w') ;

    load('WorkSpace_Tree') ; load('WorkSpace_Branch')        ; Range=[] ;
    figure ; set(gcf,'position',[0,0,3000,800]) ; TREE=Tree  ; RANGE=[] ;
for n =1:size(Tree,1)
    Tree{n,7}=exp(Tree{n,7}).^2 ;
for k1=1:size(Tree{n,3},1)
for k2=1:size(Tree{n,3},2)
if  isempty(Tree{n,3}{k1,k2})==0
    RANGE=[RANGE,min(Tree{n,3}{k1,k2}),max(Tree{n,3}{k1,k2})]    ;
    Tree{n,3}{k1,k2}=exp(Tree{n,3}{k1,k2}).^2                    ;
    Range=[Range,min(Tree{n,3}{k1,k2}),max(Tree{n,3}{k1,k2})]    ;
end
end
end
end
    SHAPE{1,1}='Tetrahedron'  ; SHAPE{1,2}=sqrt(3)               ; SHAPE{1,3}=sqrt(2)/12           ;
    SHAPE{2,1}='Cube'         ; SHAPE{2,2}=6                     ; SHAPE{2,3}=1                    ;
    SHAPE{3,1}='Octahedron'   ; SHAPE{3,2}=2*sqrt(3)             ; SHAPE{3,3}=sqrt(2)/3            ;
    SHAPE{4,1}='Dodecahedron' ; SHAPE{4,2}=3*sqrt(25+10*sqrt(5)) ; SHAPE{4,3}=(15+7*sqrt(5))/4     ;
    SHAPE{5,1}='Icosahedron'  ; SHAPE{5,2}=5*sqrt(3)             ; SHAPE{5,3}=5/12*(3+sqrt(5))     ;
    SHAPE{6,1}='Sphere'       ; SHAPE{6,2}=4*pi                  ; SHAPE{6,3}=4*pi/3               ;
    ETA=(cell2mat(SHAPE(:,2)).^(1/2))./(cell2mat(SHAPE(:,3)).^(1/3))                               ;
    MAX=max(Range)*1.001      ; MIN=min(Range)*0.999             ; c=hsv(10000) ;

    ColorMIN=[0,256,0]/256    ; ColorMAX=[256,0,256]/256         ; Nc=10000     ; Nn=1000 ; k=0.50 ; ccc=[]               ;
    c=[linspace(ColorMIN(1),ColorMAX(1),Nc)',linspace(ColorMIN(2),ColorMAX(2),Nc)',linspace(ColorMIN(3),ColorMAX(3),Nc)'] ;
for n=1:Nn
    ccc=[ccc;c(round((((n-1)/(Nn-1))^k)*(Nc-1)+1),:)]            ;
end
    figure ; c=ccc ;
for i=1:size(c,1)
    patch([0,1,1,0],[i,i,i+1,i+1],c(i,:),'edgecolor','none')     ; hold on      ;
end
    set(gca,'xtick',[],'xcolor','w') ; set(gca,'ytick',[],'ycolor','w')         ; NNN=size(c,1)    ;

    figure ; set(gcf,'position',[0,0,3000,800])                  ;
for n=1:length(Branch)
    branch=Branch{n,1} ; R2=length(branch)   ; T0=150            ;
for k2=1
    C=round((Tree{n,7}-MIN)/(MAX-MIN)*NNN)   ; 
    plot((branch{1,1}(1,1)+branch{1,1}(1,2))/2*ones(1,2),[T0-Tree{n,5},T0+7.5],'-','linewidth',2,'color',c(C,:)) ; hold on                 ;
end
for k2=1:size(Tree{n,1},2)
for k1=1:size(Tree{n,1},1)
if  Tree{n,1}(k1,k2)~=0
    nk1=k1 ; nk2=k2 ; t=[Tree{n,1}(nk1,nk2)] ; k2
    for k=2:k2
        nk1=ceil(nk1/2) ; nk2=nk2-1          ; t=[t,Tree{n,1}(nk1,nk2)]               ;
    end
    for k=1:length(Tree{n,3}{k1,k2})
        C=round((Tree{n,3}{k1,k2}(1,length(Tree{n,3}{k1,k2})+1-k)-MIN)/(MAX-MIN)*NNN) ;
        plot(branch{k2,1}(1,k1)*ones(1,2),T0-Tree{n,5}-sum(t)+t(1,1)/length(Tree{n,3}{k1,k2})*[k-1,k]   ,'-','linewidth',2,'color',c(C,:)) ;
        hold on ; clear nk1 ; clear nk2 ; clear k ;
    end
    if  (-1)^k1==1
        C=round((Tree{n,3}{k1,k2}(1,2)-MIN)/(MAX-MIN)*NNN)                            ;
        plot([branch{k2,1}(1,k1),branch{k2,1}(1,k1)-2^(R2-k2-1)],[-sum(t)+t(1,1)+T0-Tree{n,5}]*ones(1,2),'-','linewidth',2,'color',c(C,:)) ; hold on ;
    else if (-1)^k1==-1
        C=round((Tree{n,3}{k1,k2}(1,2)-MIN)/(MAX-MIN)*NNN)                            ;
        plot([branch{k2,1}(1,k1),branch{k2,1}(1,k1)+2^(R2-k2-1)],[-sum(t)+t(1,1)+T0-Tree{n,5}]*ones(1,2),'-','linewidth',2,'color',c(C,:)) ; hold on ;
        end
    end
end
end
end
        clear nk1 ; clear nk2 ; clear k  ; clear k1 ; clear k2 ; clear branch         ; clear R2 ; clear t    ;
        set(gca,'xtick',[],'xcolor','w') ; set(gca,'ytick',[],'ycolor','w')           ;
end
    Axis=T0-20*[0:1:8]        ;
for i=1:length(Axis)
    plot([-5,-4]-10,[Axis(i),Axis(i)]           ,'k-','linewidth',2) ; hold on ; clear T0        ;
end
    plot([-5,-5]-10,[min(Axis)-20/1.5,max(Axis)],'k-','linewidth',2) ; hold on ; clear i         ;
    fill([-5,-3.5,-6.5,-5]-10,min(Axis)-20+[-5,7.5,7.5,-5],'k') ; axis([-15-10,360,-70,160])     ; clear Axis ;

    figure ;
for i=log(sqrt(MIN))+0.001:0.001:log(sqrt(MAX))
    C=round((exp(i).^2-MIN)/(MAX-MIN)*NNN)     ;
    patch([0,1,1,0],[i,i,i+0.001,i+0.001],c(C,:),'edgecolor','none')    ; hold on ;
end
for i=[2.4,2.5,2.6,2.7,2.8]
    plot([0.75,0.99],[i,i],'k-','linewidth',1) ; hold on ;
end
    set(gca,'xtick',[],'xcolor','w') ; set(gca,'ytick',[],'ycolor','w') ; clear i ; clear c      ;