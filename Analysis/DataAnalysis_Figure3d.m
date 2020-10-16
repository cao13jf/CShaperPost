  % Mission 1 : Data Analysis of Figure 3d
  
    load('WorkSpace_EMBRYO') ; load('WorkSpace_CellVolume') ; load('WorkSpace_Query') ; Max=[] ; CELLNAME={} ; c=hsv(size(EMBRYO,1)) ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_CellName') ; figure('Position',[0 0 600 600])    ;
for I =1:size(EMBRYO,1)
    volume=[] ; I
for k0=1
for k1=1:size(CellName{k0,1},1)
for k2=2:size(CellName{k0,1},2)-1
if  isnumeric(CellName{k0,1}{k1,k2})==0 && length(CellVolume{k0,1}{k1,k2})==size(EMBRYO,1)
    volume=[volume,[mean(CellVolume{k0,1}{k1,k2});CellVolume{k0,1}{k1,k2}(I)]] ;
end
end
end
end
for k0=[2,3,4,5,7]
for k1=1:size(CellName{k0,1},1)
for k2=1:size(CellName{k0,1},2)-1
if  isnumeric(CellName{k0,1}{k1,k2})==0 && length(CellVolume{k0,1}{k1,k2})==size(EMBRYO,1)
    volume=[volume,[mean(CellVolume{k0,1}{k1,k2});CellVolume{k0,1}{k1,k2}(I)]] ;
end
end
end
end
    volume=[volume,[mean(CellVolume{ 6,1}{ 1, 1});CellVolume{ 6,1}{ 1, 1}(I)]] ; EMBRYO{I,14}=volume                 ;
end
for I=1:size(EMBRYO,1)
    b =sum(EMBRYO{I,14}(2,:).*EMBRYO{I,14}(1,:))/sum(EMBRYO{I,14}(1,:).^2)         ;
    yy=sum(EMBRYO{I,14}(2,:).*EMBRYO{I,14}(1,:))/sum(EMBRYO{I,14}(1,:).^2)*EMBRYO{I,14}(1,:)   ; y=EMBRYO{I,14}(2,:) ;
    SSE=sum((yy-y).^2)   ; SST=sum((y-mean(y)).^2) ; EMBRYO{I,15}=b ; EMBRYO{I,16}=1-SSE/SST   ; clear cycle_i       ;
    clear k0 ; clear k1  ; clear k2  ; clear nn    ; clear b ; clear yy ; clear y  ; clear SSE ; clear cycle_j       ;
end
    plot([0:1:7500],[0:1:7500],'k--','linewidth',1.5,'color','k')   ; hold on      ; clear SST ;
for I=1:size(EMBRYO,1)
    plot(EMBRYO{I,14}(1,:),EMBRYO{I,14}(2,:)/EMBRYO{I,15},'.','markersize',20,'color',c(I,:))  ; hold on             ;
end
    axis([0,4000,0,4000]) ;
    x=xlabel({'\rm Cell volume in average (\itV\rm_{A} / \mum^{3})'}) ; set(x,'Fontname','arial','Fontsize',18)      ;
    y=ylabel({'\rm Cell volume in samples (\itV\rm_{S} / \mum^{3})'}) ; set(y,'Fontname','arial','Fontsize',18)      ;
    set(gca,'FontSize',18,'Fontname','arial') ; axis square           ;
    set(gca,'xtick',[0,2000,4000],'FontSize',18,'Fontname','arial')   ;
    set(gca,'ytick',[0,2000,4000],'FontSize',18,'Fontname','arial')   ;

    load('WorkSpace_EMBRYO') ; load('WorkSpace_CellVolume') ; load('WorkSpace_Query') ; Max=[] ; CELLNAME={} ; c=hsv(size(EMBRYO,1)) ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_CellName') ; figure('Position',[0 0 200 200])    ;
for I =1:size(EMBRYO,1)
    volume=[] ; I
for k0=1
for k1=1:size(CellName{k0,1},1)
for k2=2:size(CellName{k0,1},2)-1
if  isnumeric(CellName{k0,1}{k1,k2})==0 && length(CellVolume{k0,1}{k1,k2})==size(EMBRYO,1)
    volume=[volume,[mean(CellVolume{k0,1}{k1,k2});CellVolume{k0,1}{k1,k2}(I)]] ;
end
end
end
end
for k0=[2,3,4,5,7]
for k1=1:size(CellName{k0,1},1)
for k2=1:size(CellName{k0,1},2)-1
if  isnumeric(CellName{k0,1}{k1,k2})==0 && length(CellVolume{k0,1}{k1,k2})==size(EMBRYO,1)
    volume=[volume,[mean(CellVolume{k0,1}{k1,k2});CellVolume{k0,1}{k1,k2}(I)]] ;
end
end
end
end
    volume=[volume,[mean(CellVolume{ 6,1}{ 1, 1});CellVolume{ 6,1}{ 1, 1}(I)]] ; EMBRYO{I,14}=volume                 ;
end
for I=1:size(EMBRYO,1)
    b =sum(EMBRYO{I,14}(2,:).*EMBRYO{I,14}(1,:))/sum(EMBRYO{I,14}(1,:).^2)         ;
    yy=sum(EMBRYO{I,14}(2,:).*EMBRYO{I,14}(1,:))/sum(EMBRYO{I,14}(1,:).^2)*EMBRYO{I,14}(1,:)   ; y=EMBRYO{I,14}(2,:) ;
    SSE=sum((yy-y).^2)   ; SST=sum((y-mean(y)).^2) ; EMBRYO{I,15}=b ; EMBRYO{I,16}=1-SSE/SST   ; clear cycle_i       ;
    clear k0 ; clear k1  ; clear k2  ; clear nn    ; clear b ; clear yy ; clear y  ; clear SSE ; clear cycle_j       ;
    plot([0:1:4000],EMBRYO{I,15}*[0:1:4000],'-','linewidth',1.0,'color',c(I,:))    ; hold on   ; clear SST           ;
end
for I=1:size(EMBRYO,1)
    plot(EMBRYO{I,14}(1,:),EMBRYO{I,14}(2,:)  ,'.','markersize',20/3*2,'color',c(I,:))         ; hold on             ;
end
    axis([0,4000,0,4000]) ; h=findobj(gca)         ; axis square        ;
    set(gca,'FontSize',14,'Fontname','arial')      ; axis square        ;
    set(gca,'xtick',[0,2000,4000],'FontSize',14,'Fontname','arial')     ;
    set(gca,'ytick',[0,2000,4000],'FontSize',14,'Fontname','arial')     ;

    figure('Position',[0 0 250 150]) ; VOLUME=[]         ;
for I=1:size(EMBRYO,1)
    VOLUME=[VOLUME;EMBRYO{I,14}(2,:)/EMBRYO{I,15}]       ;
end
    histogram(std(VOLUME)./mean(VOLUME),10,'Normalization','probability','Facecolor','b','Facealpha',0.75)           ; % axis square ;
    x=xlabel({'Variation coefficient';'of cell volume'}) ; set(x,'Fontname','arial','Fontsize' ,14) ;
    y=ylabel('Fraction') ; axis([0,0.3,0,0.45])          ; set(y,'Fontname','arial','Fontsize' ,14) ;
    set(gca,'xtick',[0,0.1,0.2,0.3],'FontSize',14,'Fontname','arial')   ;
    set(gca,'ytick',[0,0.2,0.4]    ,'FontSize',14,'Fontname','arial')   ;
    set(gca,'FontSize',14,'Fontname','arial')            ;