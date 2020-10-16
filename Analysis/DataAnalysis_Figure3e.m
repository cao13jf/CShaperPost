  % Mission 1 : Data Analysis of Figure 3e

    load('WorkSpace_EMBRYO') ; load('WorkSpace_CellSurface')    ; load('WorkSpace_Query')      ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200323 - MissingCell\WorkSpace_MissingCell') ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_CellName') ; c=hsv(size(EMBRYO,1))        ;
    Max=[]     ; figure('Position',[0 0 600 600]) ; CELLNAME={} ;
for I =1:size(EMBRYO,1)
    surface=[] ; I
for k0=1
for k1=1:size(CellName{k0,1},1)
for k2=2:size(CellName{k0,1},2)-1
if  isnumeric(CellName{k0,1}{k1,k2})==0 && length(CellSurface{k0,1}{k1,k2})==size(EMBRYO,1)
    surface=[surface,[mean(CellSurface{k0,1}{k1,k2});CellSurface{k0,1}{k1,k2}(I)]] ;
end
end
end
end
for k0=[2,3,4,5,7]
for k1=1:size(CellName{k0,1},1)
for k2=1:size(CellName{k0,1},2)-1
if  isnumeric(CellName{k0,1}{k1,k2})==0 && length(CellSurface{k0,1}{k1,k2})==size(EMBRYO,1)
    surface=[surface,[mean(CellSurface{k0,1}{k1,k2});CellSurface{k0,1}{k1,k2}(I)]] ;
end
end
end
end
    surface=[surface,[mean(CellSurface{ 6,1}{ 1, 1});CellSurface{ 6,1}{ 1, 1}(I)]] ; EMBRYO{I,14}=surface             ;
end
for I=1:size(EMBRYO,1)
    b =sum(EMBRYO{I,14}(2,:).*EMBRYO{I,14}(1,:))/sum(EMBRYO{I,14}(1,:).^2)         ;
    yy=sum(EMBRYO{I,14}(2,:).*EMBRYO{I,14}(1,:))/sum(EMBRYO{I,14}(1,:).^2)*EMBRYO{I,14}(1,:)   ; y=EMBRYO{I,14}(2,:)  ;
    SSE=sum((yy-y).^2)   ; SST=sum((y-mean(y)).^2) ; EMBRYO{I,15}=b ; EMBRYO{I,16}=1-SSE/SST   ; clear cycle_i        ;
    clear k0 ; clear k1  ; clear k2  ; clear nn    ; clear b ; clear yy ; clear y  ; clear SSE ; clear cycle_j        ;
end
    plot([0:1:1700],[0:1:1700],'--','linewidth',1.5,'color','k') ; hold on  ; clear SST        ;
for I=1:size(EMBRYO,1)
    plot(EMBRYO{I,14}(1,:),EMBRYO{I,14}(2,:)./EMBRYO{I,15},'.','markersize',20,'color',c(I,:)) ; hold on              ;
end
    axis([0,1700,0,1700])     ;
    x=xlabel({'\rm Cell surface area in average (\itS\rm_{A} / \mum^{2})'}) ; set(x,'Fontname','arial','Fontsize',18) ;
    y=ylabel({'\rm Cell surface area in samples (\itS\rm_{S} / \mum^{2})'}) ; set(y,'Fontname','arial','Fontsize',18) ;
    set(gca,'FontSize',18,'Fontname','arial') ; axis square                 ;
    set(gca,'xtick',[0,800,1600],'FontSize',18,'Fontname','arial')          ;
    set(gca,'ytick',[0,800,1600],'FontSize',18,'Fontname','arial')          ;

    load('WorkSpace_EMBRYO') ; load('WorkSpace_CellSurface')    ; load('WorkSpace_Query')      ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200323 - MissingCell\WorkSpace_MissingCell') ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_CellName') ; c=hsv(size(EMBRYO,1))        ;
    Max=[]     ; figure('Position',[0 0 200 200]) ; CELLNAME={} ;
for I =1:size(EMBRYO,1)
    surface=[] ; I
for k0=1
for k1=1:size(CellName{k0,1},1)
for k2=2:size(CellName{k0,1},2)-1
if  isnumeric(CellName{k0,1}{k1,k2})==0 && length(CellSurface{k0,1}{k1,k2})==size(EMBRYO,1)
    surface=[surface,[mean(CellSurface{k0,1}{k1,k2});CellSurface{k0,1}{k1,k2}(I)]] ;
end
end
end
end
for k0=[2,3,4,5,7]
for k1=1:size(CellName{k0,1},1)
for k2=1:size(CellName{k0,1},2)-1
if  isnumeric(CellName{k0,1}{k1,k2})==0 && length(CellSurface{k0,1}{k1,k2})==size(EMBRYO,1)
    surface=[surface,[mean(CellSurface{k0,1}{k1,k2});CellSurface{k0,1}{k1,k2}(I)]] ;
end
end
end
end
    surface=[surface,[mean(CellSurface{ 6,1}{ 1, 1});CellSurface{ 6,1}{ 1, 1}(I)]] ; EMBRYO{I,14}=surface            ;
end
for I=1:size(EMBRYO,1)
    b =sum(EMBRYO{I,14}(2,:).*EMBRYO{I,14}(1,:))/sum(EMBRYO{I,14}(1,:).^2)         ;
    yy=sum(EMBRYO{I,14}(2,:).*EMBRYO{I,14}(1,:))/sum(EMBRYO{I,14}(1,:).^2)*EMBRYO{I,14}(1,:)   ; y=EMBRYO{I,14}(2,:) ;
    SSE=sum((yy-y).^2)   ; SST=sum((y-mean(y)).^2) ; EMBRYO{I,15}=b ; EMBRYO{I,16}=1-SSE/SST   ; clear cycle_i       ;
    clear k0 ; clear k1  ; clear k2  ; clear nn    ; clear b ; clear yy  ; clear y ; clear SSE ; clear cycle_j       ;
    plot([0:1:1700],EMBRYO{I,15}*[0:1:1700],'-','linewidth',1.0,'color',c(I,:))    ; hold on   ; clear SST           ;
end
for I=1:size(EMBRYO,1)
    plot(EMBRYO{I,14}(1,:),EMBRYO{I,14}(2,:),'.','markersize',20/3*2,'color',c(I,:))           ; hold on             ;
end
    axis([0,1700,0,1700])            ;
    set(gca,'FontSize',14,'Fontname','arial') ; axis square              ;
    set(gca,'xtick',[0,800,1600],'FontSize',14,'Fontname','arial')       ;
    set(gca,'ytick',[0,800,1600],'FontSize',14,'Fontname','arial')       ;

    figure('Position',[0 0 250 150]) ; SURFACE=[]              ;
for I=1:size(EMBRYO,1)
    SURFACE=[SURFACE;EMBRYO{I,14}(2,:)/EMBRYO{I,15}]           ;
end
    histogram(std(SURFACE)./mean(SURFACE),9,'Normalization','probability','Facecolor','b','Facealpha',0.75) ;
    x=xlabel({'Variation coefficient';'of cell surface area'}) ; set(x,'Fontname','arial','Fontsize' ,14)   ;
    y=ylabel('Fraction') ; axis([0,0.3,0,0.45])                ; set(y,'Fontname','arial','Fontsize' ,14)   ;
    set(gca,'xtick',[0,0.10,0.20,0.30],'FontSize',14,'Fontname','arial') ;
    set(gca,'ytick',[0,0.20,0.40]     ,'FontSize',14,'Fontname','arial') ;
    set(gca,'FontSize',14,'Fontname','arial')                  ;