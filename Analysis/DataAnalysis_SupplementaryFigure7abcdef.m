  % Mission 1 : Data Analysis of Supplementary Figure 7abcdef

  % Volume Asymmetry
    load('WorkSpace_EMBRYO') ; load('WorkSpace_CellVolume') ; c=hsv(size(EMBRYO,1)) ; figure('Position',[0 0 550 575])        ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200323 - MissingCell\WorkSpace_MissingCell') ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200324 - LostRatio\WorkSpace_CellCycle')     ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_AveCycle') ; 
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_WTinfo'  ) ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_CellName') ;
for I =1:size(EMBRYO,1)
    volume=[] ; I
for k0=1
for k1=1:size(CellName{k0,1},1)
for k2=1:size(CellName{k0,1},2)-2
if  isnumeric(CellName{k0,1}{k1,k2})==0                   && length(CellVolume{k0,1}{k1,k2})==size(EMBRYO,1)
if  length(CellVolume{k0,1}{2*k1-1,k2+1})==size(EMBRYO,1) && length(CellVolume{k0,1}{2*k1,k2+1})==size(EMBRYO,1)
    if  mean(CellVolume{k0,1}{2*k1,k2+1})/mean(CellVolume{k0,1}{2*k1-1,k2+1})>1
        volume=[volume,[mean(CellVolume{k0,1}{2*k1,k2+1})/mean(CellVolume{k0,1}{2*k1-1,k2+1});CellVolume{k0,1}{2*k1,k2+1}(I)/CellVolume{k0,1}{2*k1-1,k2+1}(I)]]       ;
    end
    if  mean(CellVolume{k0,1}{2*k1,k2+1})/mean(CellVolume{k0,1}{2*k1-1,k2+1})<1
        volume=[volume,[mean(CellVolume{k0,1}{2*k1,k2+1})/mean(CellVolume{k0,1}{2*k1-1,k2+1});CellVolume{k0,1}{2*k1,k2+1}(I)/CellVolume{k0,1}{2*k1-1,k2+1}(I)].^(-1)] ;
    end
end
end
end
end
end
for k0=[2,3,4,5]
for k1=1:size(CellName{k0,1},1)
for k2=1:size(CellName{k0,1},2)-2
if  isnumeric(CellName{k0,1}{k1,k2})==0                   && length(CellVolume{k0,1}{k1,k2})==size(EMBRYO,1)
if  length(CellVolume{k0,1}{2*k1-1,k2+1})==size(EMBRYO,1) && length(CellVolume{k0,1}{2*k1,k2+1})==size(EMBRYO,1)
    if  mean(CellVolume{k0,1}{2*k1,k2+1})/mean(CellVolume{k0,1}{2*k1-1,k2+1})>1
        volume=[volume,[mean(CellVolume{k0,1}{2*k1,k2+1})/mean(CellVolume{k0,1}{2*k1-1,k2+1});CellVolume{k0,1}{2*k1,k2+1}(I)/CellVolume{k0,1}{2*k1-1,k2+1}(I)]]       ;
    end
    if  mean(CellVolume{k0,1}{2*k1,k2+1})/mean(CellVolume{k0,1}{2*k1-1,k2+1})<1
        volume=[volume,[mean(CellVolume{k0,1}{2*k1,k2+1})/mean(CellVolume{k0,1}{2*k1-1,k2+1});CellVolume{k0,1}{2*k1,k2+1}(I)/CellVolume{k0,1}{2*k1-1,k2+1}(I)].^(-1)] ;
    end
end
end
end
end
end
    if  mean(CellVolume{2,1}{1,1})/mean(CellVolume{3,1}{1,1})>1
        volume=[volume,[mean(CellVolume{2,1}{1,1})/mean(CellVolume{3,1}{1,1});CellVolume{2,1}{1,1}(I)/CellVolume{3,1}{1,1}(I)]]       ;
    end
    if  mean(CellVolume{2,1}{1,1})/mean(CellVolume{3,1}{1,1})<1
        volume=[volume,[mean(CellVolume{2,1}{1,1})/mean(CellVolume{3,1}{1,1});CellVolume{2,1}{1,1}(I)/CellVolume{3,1}{1,1}(I)].^(-1)] ;
    end
    if  mean(CellVolume{4,1}{1,1})/mean(CellVolume{6,1}{1,1})>1
        volume=[volume,[mean(CellVolume{4,1}{1,1})/mean(CellVolume{6,1}{1,1});CellVolume{4,1}{1,1}(I)/CellVolume{6,1}{1,1}(I)]]       ;
    end
    if  mean(CellVolume{4,1}{1,1})/mean(CellVolume{6,1}{1,1})<1
        volume=[volume,[mean(CellVolume{4,1}{1,1})/mean(CellVolume{6,1}{1,1});CellVolume{4,1}{1,1}(I)/CellVolume{6,1}{1,1}(I)].^(-1)] ;
    end
    if  mean(CellVolume{5,1}{1,1})/mean(CellVolume{7,1}{1,1})>1
        volume=[volume,[mean(CellVolume{5,1}{1,1})/mean(CellVolume{7,1}{1,1});CellVolume{5,1}{1,1}(I)/CellVolume{7,1}{1,1}(I)]]       ;
    end
    if  mean(CellVolume{5,1}{1,1})/mean(CellVolume{7,1}{1,1})<1
        volume=[volume,[mean(CellVolume{5,1}{1,1})/mean(CellVolume{7,1}{1,1});CellVolume{5,1}{1,1}(I)/CellVolume{7,1}{1,1}(I)].^(-1)] ;
    end
        EMBRYO{I,14}=volume           ;
end
for I=1:size(EMBRYO,1)
    b=1      ; yy=b*EMBRYO{I,14}(1,:) ; y=EMBRYO{I,14}(2,:)  ;
    SSE=sum((yy-y).^2)   ; SST=sum((y-mean(y)).^2) ; EMBRYO{I,15}=b ; EMBRYO{I,16}=1-SSE/SST    ; clear cycle_i ;
    clear k0 ; clear k1  ; clear k2   ; clear nn   ; clear b ; clear yy ; clear y   ; clear SSE ; clear cycle_j ;
end
    plot([0:1:5],[0:1:5],'k--','linewidth',1.5,'color','k')  ; hold on  ; clear SST ;
for I=1:size(EMBRYO,1)
    plot(EMBRYO{I,14}(1,:),EMBRYO{I,14}(2,:),'.','markersize',20,'color',c(I,:))    ; hold on   ;
end
    axis([0,4,0,4]) ;
    x=xlabel({'\rm Volume ratio between sister cells';'in average (\itV\rm_{A,1} / \itV\rm_{A,2})'}) ; set(x,'Fontname','arial','Fontsize',18) ;
    y=ylabel({'\rm Volume ratio between sister cells';'in samples (\itV\rm_{S,1} / \itV\rm_{S,2})'}) ; set(y,'Fontname','arial','Fontsize',18) ;
    set(gca,'FontSize',18,'Fontname','arial') ; axis square             ;
    set(gca,'xtick',[0,1,2,3,4],'FontSize',18,'Fontname','arial')       ;
    set(gca,'ytick',[0,1,2,3,4],'FontSize',18,'Fontname','arial')       ;
    
    figure('Position',[0 0 525 475]) ; VOLUME=[]         ;
for I=1:size(EMBRYO,1)
    VOLUME=[VOLUME;EMBRYO{I,14}(2,:)]                    ;
end
    histogram(std(VOLUME)./mean(VOLUME),11,'Normalization','probability','Facecolor','b','Facealpha',0.75) ;
    x=xlabel({'\rm Variation coefficient of';'volume ratio between sister cells'}) ; set(x,'Fontname','arial','Fontsize' ,18) ;
    y=ylabel('Fraction') ; axis([0,0.25,0,0.3])                                    ; set(y,'Fontname','arial','Fontsize' ,18) ;
    set(gca,'xtick',[0,0.1,0.2,0.3],'FontSize',18,'Fontname','arial')   ;
    set(gca,'ytick',[0,0.1,0.2,0.3],'FontSize',18,'Fontname','arial')   ;
    set(gca,'FontSize',18,'Fontname','arial')            ;

  % Surface Area Asymmetry
    load('WorkSpace_EMBRYO') ; load('WorkSpace_CellSurface') ; c=hsv(size(EMBRYO,1)) ; figure('Position',[0 0 550 575])       ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200323 - MissingCell\WorkSpace_MissingCell') ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200324 - LostRatio\WorkSpace_CellCycle')     ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_AveCycle') ; CellVolume=CellSurface       ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_WTinfo'  ) ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_CellName') ;
for I =1:size(EMBRYO,1)
    volume=[] ; I
for k0=1
for k1=1:size(CellName{k0,1},1)
for k2=1:size(CellName{k0,1},2)-2
if  isnumeric(CellName{k0,1}{k1,k2})==0                   && length(CellVolume{k0,1}{k1,k2})==size(EMBRYO,1)
if  length(CellVolume{k0,1}{2*k1-1,k2+1})==size(EMBRYO,1) && length(CellVolume{k0,1}{2*k1,k2+1})==size(EMBRYO,1)
    if  mean(CellVolume{k0,1}{2*k1,k2+1})/mean(CellVolume{k0,1}{2*k1-1,k2+1})>1
        volume=[volume,[mean(CellVolume{k0,1}{2*k1,k2+1})/mean(CellVolume{k0,1}{2*k1-1,k2+1});CellVolume{k0,1}{2*k1,k2+1}(I)/CellVolume{k0,1}{2*k1-1,k2+1}(I)]]       ;
    end
    if  mean(CellVolume{k0,1}{2*k1,k2+1})/mean(CellVolume{k0,1}{2*k1-1,k2+1})<1
        volume=[volume,[mean(CellVolume{k0,1}{2*k1,k2+1})/mean(CellVolume{k0,1}{2*k1-1,k2+1});CellVolume{k0,1}{2*k1,k2+1}(I)/CellVolume{k0,1}{2*k1-1,k2+1}(I)].^(-1)] ;
    end
end
end
end
end
end
for k0=[2,3,4,5]
for k1=1:size(CellName{k0,1},1)
for k2=1:size(CellName{k0,1},2)-2
if  isnumeric(CellName{k0,1}{k1,k2})==0                   && length(CellVolume{k0,1}{k1,k2})==size(EMBRYO,1)
if  length(CellVolume{k0,1}{2*k1-1,k2+1})==size(EMBRYO,1) && length(CellVolume{k0,1}{2*k1,k2+1})==size(EMBRYO,1)
    if  mean(CellVolume{k0,1}{2*k1,k2+1})/mean(CellVolume{k0,1}{2*k1-1,k2+1})>1
        volume=[volume,[mean(CellVolume{k0,1}{2*k1,k2+1})/mean(CellVolume{k0,1}{2*k1-1,k2+1});CellVolume{k0,1}{2*k1,k2+1}(I)/CellVolume{k0,1}{2*k1-1,k2+1}(I)]]       ;
    end
    if  mean(CellVolume{k0,1}{2*k1,k2+1})/mean(CellVolume{k0,1}{2*k1-1,k2+1})<1
        volume=[volume,[mean(CellVolume{k0,1}{2*k1,k2+1})/mean(CellVolume{k0,1}{2*k1-1,k2+1});CellVolume{k0,1}{2*k1,k2+1}(I)/CellVolume{k0,1}{2*k1-1,k2+1}(I)].^(-1)] ;
    end
end
end
end
end
end
    if  mean(CellVolume{2,1}{1,1})/mean(CellVolume{3,1}{1,1})>1
        volume=[volume,[mean(CellVolume{2,1}{1,1})/mean(CellVolume{3,1}{1,1});CellVolume{2,1}{1,1}(I)/CellVolume{3,1}{1,1}(I)]]       ;
    end
    if  mean(CellVolume{2,1}{1,1})/mean(CellVolume{3,1}{1,1})<1
        volume=[volume,[mean(CellVolume{2,1}{1,1})/mean(CellVolume{3,1}{1,1});CellVolume{2,1}{1,1}(I)/CellVolume{3,1}{1,1}(I)].^(-1)] ;
    end
    if  mean(CellVolume{4,1}{1,1})/mean(CellVolume{6,1}{1,1})>1
        volume=[volume,[mean(CellVolume{4,1}{1,1})/mean(CellVolume{6,1}{1,1});CellVolume{4,1}{1,1}(I)/CellVolume{6,1}{1,1}(I)]]       ;
    end
    if  mean(CellVolume{4,1}{1,1})/mean(CellVolume{6,1}{1,1})<1
        volume=[volume,[mean(CellVolume{4,1}{1,1})/mean(CellVolume{6,1}{1,1});CellVolume{4,1}{1,1}(I)/CellVolume{6,1}{1,1}(I)].^(-1)] ;
    end
    if  mean(CellVolume{5,1}{1,1})/mean(CellVolume{7,1}{1,1})>1
        volume=[volume,[mean(CellVolume{5,1}{1,1})/mean(CellVolume{7,1}{1,1});CellVolume{5,1}{1,1}(I)/CellVolume{7,1}{1,1}(I)]]       ;
    end
    if  mean(CellVolume{5,1}{1,1})/mean(CellVolume{7,1}{1,1})<1
        volume=[volume,[mean(CellVolume{5,1}{1,1})/mean(CellVolume{7,1}{1,1});CellVolume{5,1}{1,1}(I)/CellVolume{7,1}{1,1}(I)].^(-1)] ;
    end
        EMBRYO{I,14}=volume           ;
end
for I=1:size(EMBRYO,1)
    b=1      ; yy=b*EMBRYO{I,14}(1,:) ; y=EMBRYO{I,14}(2,:)  ;
    SSE=sum((yy-y).^2)   ; SST=sum((y-mean(y)).^2) ; EMBRYO{I,15}=b ; EMBRYO{I,16}=1-SSE/SST    ; clear cycle_i ;
    clear k0 ; clear k1  ; clear k2   ; clear nn   ; clear b ; clear yy ; clear y   ; clear SSE ; clear cycle_j ;
end
    plot([0:1:5],[0:1:5],'k--','linewidth',1.5,'color','k')  ; hold on  ; clear SST ;
for I=1:size(EMBRYO,1)
    plot(EMBRYO{I,14}(1,:),EMBRYO{I,14}(2,:),'.','markersize',20,'color',c(I,:))    ; hold on   ;
end
    axis([0,3,0,3]) ;
    x=xlabel({'\rm Surface area ratio between sister cells';'in average (\itS\rm_{A,1} / \itS\rm_{A,2})'}) ; set(x,'Fontname','arial','Fontsize',18) ;
    y=ylabel({'\rm Surface area ratio between sister cells';'in samples (\itS\rm_{S,1} / \itS\rm_{S,2})'}) ; set(y,'Fontname','arial','Fontsize',18) ;
    set(gca,'FontSize',18,'Fontname','arial') ; axis square             ;
    set(gca,'xtick',[0,1,2,3],'FontSize',18,'Fontname','arial')         ;
    set(gca,'ytick',[0,1,2,3],'FontSize',18,'Fontname','arial')         ;
    
    figure('Position',[0 0 525 475]) ; VOLUME=[]             ;
for I=1:size(EMBRYO,1)
    VOLUME=[VOLUME;EMBRYO{I,14}(2,:)]                        ;
end
    histogram(std(VOLUME)./mean(VOLUME),10,'Normalization','probability','Facecolor','b','Facealpha',0.75) ;
    x=xlabel({'\rm Variation coefficient of';'surface area ratio between sister cells'}) ; set(x,'Fontname','arial','Fontsize' ,18) ;
    y=ylabel('Fraction') ; axis([0,0.25,0,0.3])                                          ; set(y,'Fontname','arial','Fontsize' ,18) ;
    set(gca,'xtick',[0,0.1,0.2,0.3],'FontSize',18,'Fontname','arial')   ;
    set(gca,'ytick',[0,0.1,0.2,0.3],'FontSize',18,'Fontname','arial')   ;
    set(gca,'FontSize',18,'Fontname','arial')                ;

  % Symmetry VS Asymmetry
    load('WorkSpace_EMBRYO') ; load('WorkSpace_CellVolume')  ; c=hsv(size(EMBRYO,1)) ; figure('Position',[0 0 500 525])       ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200323 - MissingCell\WorkSpace_MissingCell') ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200324 - LostRatio\WorkSpace_CellCycle')     ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_AveCycle') ; 
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_WTinfo'  ) ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_CellName') ;
for I =1:size(EMBRYO,1)
    volume=[] ; I
for k0=1
for k1=1:size(CellName{k0,1},1)
for k2=1:size(CellName{k0,1},2)-2
if  isnumeric(CellName{k0,1}{k1,k2})==0                   && length(CellVolume{k0,1}{k1,k2})==size(EMBRYO,1)
if  length(CellVolume{k0,1}{2*k1-1,k2+1})==size(EMBRYO,1) && length(CellVolume{k0,1}{2*k1,k2+1})==size(EMBRYO,1)
    if  mean(CellVolume{k0,1}{2*k1,k2+1})/mean(CellVolume{k0,1}{2*k1-1,k2+1})>1
        volume=[volume,[mean(CellVolume{k0,1}{2*k1,k2+1})/mean(CellVolume{k0,1}{2*k1-1,k2+1});CellVolume{k0,1}{2*k1,k2+1}(I)/CellVolume{k0,1}{2*k1-1,k2+1}(I)]]       ;
    end
    if  mean(CellVolume{k0,1}{2*k1,k2+1})/mean(CellVolume{k0,1}{2*k1-1,k2+1})<1
        volume=[volume,[mean(CellVolume{k0,1}{2*k1,k2+1})/mean(CellVolume{k0,1}{2*k1-1,k2+1});CellVolume{k0,1}{2*k1,k2+1}(I)/CellVolume{k0,1}{2*k1-1,k2+1}(I)].^(-1)] ;
    end
end
end
end
end
end
for k0=[2,3,4,5]
for k1=1:size(CellName{k0,1},1)
for k2=1:size(CellName{k0,1},2)-2
if  isnumeric(CellName{k0,1}{k1,k2})==0                   && length(CellVolume{k0,1}{k1,k2})==size(EMBRYO,1)
if  length(CellVolume{k0,1}{2*k1-1,k2+1})==size(EMBRYO,1) && length(CellVolume{k0,1}{2*k1,k2+1})==size(EMBRYO,1)
    if  mean(CellVolume{k0,1}{2*k1,k2+1})/mean(CellVolume{k0,1}{2*k1-1,k2+1})>1
        volume=[volume,[mean(CellVolume{k0,1}{2*k1,k2+1})/mean(CellVolume{k0,1}{2*k1-1,k2+1});CellVolume{k0,1}{2*k1,k2+1}(I)/CellVolume{k0,1}{2*k1-1,k2+1}(I)]]       ;
    end
    if  mean(CellVolume{k0,1}{2*k1,k2+1})/mean(CellVolume{k0,1}{2*k1-1,k2+1})<1
        volume=[volume,[mean(CellVolume{k0,1}{2*k1,k2+1})/mean(CellVolume{k0,1}{2*k1-1,k2+1});CellVolume{k0,1}{2*k1,k2+1}(I)/CellVolume{k0,1}{2*k1-1,k2+1}(I)].^(-1)] ;
    end
end
end
end
end
end
    if  mean(CellVolume{2,1}{1,1})/mean(CellVolume{3,1}{1,1})>1
        volume=[volume,[mean(CellVolume{2,1}{1,1})/mean(CellVolume{3,1}{1,1});CellVolume{2,1}{1,1}(I)/CellVolume{3,1}{1,1}(I)]]       ;
    end
    if  mean(CellVolume{2,1}{1,1})/mean(CellVolume{3,1}{1,1})<1
        volume=[volume,[mean(CellVolume{2,1}{1,1})/mean(CellVolume{3,1}{1,1});CellVolume{2,1}{1,1}(I)/CellVolume{3,1}{1,1}(I)].^(-1)] ;
    end
    if  mean(CellVolume{4,1}{1,1})/mean(CellVolume{6,1}{1,1})>1
        volume=[volume,[mean(CellVolume{4,1}{1,1})/mean(CellVolume{6,1}{1,1});CellVolume{4,1}{1,1}(I)/CellVolume{6,1}{1,1}(I)]]       ;
    end
    if  mean(CellVolume{4,1}{1,1})/mean(CellVolume{6,1}{1,1})<1
        volume=[volume,[mean(CellVolume{4,1}{1,1})/mean(CellVolume{6,1}{1,1});CellVolume{4,1}{1,1}(I)/CellVolume{6,1}{1,1}(I)].^(-1)] ;
    end
    if  mean(CellVolume{5,1}{1,1})/mean(CellVolume{7,1}{1,1})>1
        volume=[volume,[mean(CellVolume{5,1}{1,1})/mean(CellVolume{7,1}{1,1});CellVolume{5,1}{1,1}(I)/CellVolume{7,1}{1,1}(I)]]       ;
    end
    if  mean(CellVolume{5,1}{1,1})/mean(CellVolume{7,1}{1,1})<1
        volume=[volume,[mean(CellVolume{5,1}{1,1})/mean(CellVolume{7,1}{1,1});CellVolume{5,1}{1,1}(I)/CellVolume{7,1}{1,1}(I)].^(-1)] ;
    end
        EMBRYO{I,14}=volume           ;
end
for I=1:size(EMBRYO,1)
    b=1      ; yy=b*EMBRYO{I,14}(1,:) ; y=EMBRYO{I,14}(2,:)  ;
    SSE=sum((yy-y).^2)   ; SST=sum((y-mean(y)).^2) ; EMBRYO{I,15}=b  ; EMBRYO{I,16}=1-SSE/SST   ; clear cycle_i ;
    clear k0 ; clear k1  ; clear k2   ; clear nn   ; clear b ; clear yy ; clear y   ; clear SSE ; clear cycle_j ;
end
    figure('Position',[0 0 575 550])  ; VOLUME=[]       ;
for I=1:size(EMBRYO,1)
    VOLUME=[VOLUME;EMBRYO{I,14}(2,:)]                   ;
end
    plot(mean(VOLUME),std(VOLUME)./mean(VOLUME),'.','markersize',30,'color',[0.6,0.6,0.6]) ; hold on ;
    p=polyfit(mean(VOLUME),std(VOLUME)./mean(VOLUME),1) ;
    plot([0.25:0.01:3.75],p(1)*[0.25:0.01:3.75]+p(2),'k--','linewidth',1.5,'color','k')    ; hold on ; clear SST              ;
    x=xlabel({'\rm Volume ratio between sister cells';'in average (\itV\rm_{A,1} / \itV\rm_{A,2})'}) ; set(x,'Fontname','arial','Fontsize',18) ;
    y=ylabel({'\rm Variation coefficient of';'volume ratio between sister cells'}) ; set(y,'Fontname','arial','Fontsize' ,18) ;
    axis([0,4,0,0.25])                                  ;
    set(gca,'xtick',[0,1,2,3,4]    ,'FontSize',18,'Fontname','arial') ;
    set(gca,'ytick',[0,0.1,0.2,0.3],'FontSize',18,'Fontname','arial') ;
    set(gca,'FontSize',18,'Fontname','arial')           ;
    
    load('WorkSpace_EMBRYO') ; load('WorkSpace_CellSurface') ; c=hsv(size(EMBRYO,1))           ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200323 - MissingCell\WorkSpace_MissingCell') ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200324 - LostRatio\WorkSpace_CellCycle')     ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_AveCycle') ; CellVolume=CellSurface       ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_WTinfo'  ) ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_CellName') ;
for I =1:size(EMBRYO,1)
    volume=[] ; I
for k0=1
for k1=1:size(CellName{k0,1},1)
for k2=1:size(CellName{k0,1},2)-2
if  isnumeric(CellName{k0,1}{k1,k2})==0                   && length(CellVolume{k0,1}{k1,k2})==size(EMBRYO,1)
if  length(CellVolume{k0,1}{2*k1-1,k2+1})==size(EMBRYO,1) && length(CellVolume{k0,1}{2*k1,k2+1})==size(EMBRYO,1)
    if  mean(CellVolume{k0,1}{2*k1,k2+1})/mean(CellVolume{k0,1}{2*k1-1,k2+1})>1
        volume=[volume,[mean(CellVolume{k0,1}{2*k1,k2+1})/mean(CellVolume{k0,1}{2*k1-1,k2+1});CellVolume{k0,1}{2*k1,k2+1}(I)/CellVolume{k0,1}{2*k1-1,k2+1}(I)]]       ;
    end
    if  mean(CellVolume{k0,1}{2*k1,k2+1})/mean(CellVolume{k0,1}{2*k1-1,k2+1})<1
        volume=[volume,[mean(CellVolume{k0,1}{2*k1,k2+1})/mean(CellVolume{k0,1}{2*k1-1,k2+1});CellVolume{k0,1}{2*k1,k2+1}(I)/CellVolume{k0,1}{2*k1-1,k2+1}(I)].^(-1)] ;
    end
end
end
end
end
end
for k0=[2,3,4,5]
for k1=1:size(CellName{k0,1},1)
for k2=1:size(CellName{k0,1},2)-2
if  isnumeric(CellName{k0,1}{k1,k2})==0                   && length(CellVolume{k0,1}{k1,k2})==size(EMBRYO,1)
if  length(CellVolume{k0,1}{2*k1-1,k2+1})==size(EMBRYO,1) && length(CellVolume{k0,1}{2*k1,k2+1})==size(EMBRYO,1)
    if  mean(CellVolume{k0,1}{2*k1,k2+1})/mean(CellVolume{k0,1}{2*k1-1,k2+1})>1
        volume=[volume,[mean(CellVolume{k0,1}{2*k1,k2+1})/mean(CellVolume{k0,1}{2*k1-1,k2+1});CellVolume{k0,1}{2*k1,k2+1}(I)/CellVolume{k0,1}{2*k1-1,k2+1}(I)]]       ;
    end
    if  mean(CellVolume{k0,1}{2*k1,k2+1})/mean(CellVolume{k0,1}{2*k1-1,k2+1})<1
        volume=[volume,[mean(CellVolume{k0,1}{2*k1,k2+1})/mean(CellVolume{k0,1}{2*k1-1,k2+1});CellVolume{k0,1}{2*k1,k2+1}(I)/CellVolume{k0,1}{2*k1-1,k2+1}(I)].^(-1)] ;
    end
end
end
end
end
end
    if  mean(CellVolume{2,1}{1,1})/mean(CellVolume{3,1}{1,1})>1
        volume=[volume,[mean(CellVolume{2,1}{1,1})/mean(CellVolume{3,1}{1,1});CellVolume{2,1}{1,1}(I)/CellVolume{3,1}{1,1}(I)]]       ;
    end
    if  mean(CellVolume{2,1}{1,1})/mean(CellVolume{3,1}{1,1})<1
        volume=[volume,[mean(CellVolume{2,1}{1,1})/mean(CellVolume{3,1}{1,1});CellVolume{2,1}{1,1}(I)/CellVolume{3,1}{1,1}(I)].^(-1)] ;
    end
    if  mean(CellVolume{4,1}{1,1})/mean(CellVolume{6,1}{1,1})>1
        volume=[volume,[mean(CellVolume{4,1}{1,1})/mean(CellVolume{6,1}{1,1});CellVolume{4,1}{1,1}(I)/CellVolume{6,1}{1,1}(I)]]       ;
    end
    if  mean(CellVolume{4,1}{1,1})/mean(CellVolume{6,1}{1,1})<1
        volume=[volume,[mean(CellVolume{4,1}{1,1})/mean(CellVolume{6,1}{1,1});CellVolume{4,1}{1,1}(I)/CellVolume{6,1}{1,1}(I)].^(-1)] ;
    end
    if  mean(CellVolume{5,1}{1,1})/mean(CellVolume{7,1}{1,1})>1
        volume=[volume,[mean(CellVolume{5,1}{1,1})/mean(CellVolume{7,1}{1,1});CellVolume{5,1}{1,1}(I)/CellVolume{7,1}{1,1}(I)]]       ;
    end
    if  mean(CellVolume{5,1}{1,1})/mean(CellVolume{7,1}{1,1})<1
        volume=[volume,[mean(CellVolume{5,1}{1,1})/mean(CellVolume{7,1}{1,1});CellVolume{5,1}{1,1}(I)/CellVolume{7,1}{1,1}(I)].^(-1)] ;
    end
        EMBRYO{I,14}=volume           ;
end
for I=1:size(EMBRYO,1)
    b=1      ; yy=b*EMBRYO{I,14}(1,:) ; y=EMBRYO{I,14}(2,:)         ;
    SSE=sum((yy-y).^2)   ; SST=sum((y-mean(y)).^2) ; EMBRYO{I,15}=b ; EMBRYO{I,16}=1-SSE/SST    ; clear cycle_i  ;
    clear k0 ; clear k1  ; clear k2   ; clear nn   ; clear b ; clear yy ; clear y   ; clear SSE ; clear cycle_j  ;
end
    figure('Position',[0 0 575 550])  ; VOLUME=[]       ;
for I=1:size(EMBRYO,1)
    VOLUME=[VOLUME;EMBRYO{I,14}(2,:)]                   ;
end
    plot(mean(VOLUME),std(VOLUME)./mean(VOLUME),'.','markersize',30,'color',[0.6,0.6,0.6]) ; hold on ;
    p=polyfit(mean(VOLUME),std(VOLUME)./mean(VOLUME),1) ;
    plot([0.25:0.01:2.75],p(1)*[0.25:0.01:2.75]+p(2),'k--','linewidth',1.5,'color','k')    ; hold on ; clear SST            ;
    x=xlabel({'\rm Surface area ratio between sister cells';'in average (\itS\rm_{A,1} / \itS\rm_{A,2})'}) ; set(x,'Fontname','arial','Fontsize',18) ;
    y=ylabel({'\rm Variation coefficient of';'surface area ratio between sister cells'})                   ; set(y,'Fontname','arial','Fontsize',18) ;
    axis([0,3,0,0.15])                                  ;
    set(gca,'xtick',[0,1,2,3]         ,'FontSize',18,'Fontname','arial') ;
    set(gca,'ytick',[0,0.05,0.10,0.15],'FontSize',18,'Fontname','arial') ; set(gca,'YTickLabel',{'0';'0.05';'0.10';'0.15'}) ;
    set(gca,'FontSize',18,'Fontname','arial')           ;