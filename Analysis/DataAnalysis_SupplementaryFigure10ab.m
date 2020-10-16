  % Mission 1 : Data Analysis of Supplementary Figure 10ab

    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200331 - Volume&Surface_PLUS\WorkSpace_EMBRYO') ; embryo=EMBRYO ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200324 - LostRatio\WorkSpace_EMBRYO')           ;
for k1=1:size(embryo,1)
for k2=1:size(embryo,2)
if  isempty(EMBRYO{k1,k2})==1 && isempty(embryo{k1,k2})==0
    EMBRYO{k1,k2}=embryo{k1,k2}               ;
end
end
end
    save('WorkSpace_EMBRYO','EMBRYO','-v7.3') ; clear k1 ; clear k2 ; clear embryo             ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200331 - Volume&Surface_PLUS\WorkSpace_CellVolume')             ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200323 - MissingCell\WorkSpace_MissingCell')    ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200324 - LostRatio\WorkSpace_CellCycle')        ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_AveCycle') ; 
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_WTinfo'  ) ;
    
  % AB + MS Cells
    c=hsv(length(AveCycle)) ; figure('Position',[100 100 200 200])        ; CYCLE=[]           ;
for k0=1:2
for k1=1:size(AveCycle{k0,1},1)
for k2=1:size(AveCycle{k0,1},2)
if  AveCycle{k0,1}(k1,k2)~=0 && length(CellVolume{k0,1}{k1,k2})==size(EMBRYO,1)
    cycle=CellCycle{k0,1}{k1,k2}./cell2mat(EMBRYO(:,20))' ; volume=CellVolume{k0,1}{k1,k2}./cell2mat(EMBRYO(:,15))'            ;
    loglog(mean(volume)+[-1,1]*std(volume),mean(cycle)*[ 1,1] ,'-','linewidth',1.5,'color',c(k0,:)) ; hold on                  ;
    loglog(mean(volume)*[ 1,1],mean(cycle)+[-1,1]*std(cycle)  ,'-','linewidth',1.5,'color',c(k0,:)) ; hold on                  ;
    loglog(mean(volume),mean(cycle),'.','markersize',15,'color',c(k0,:))  ; hold on ; CYCLE=[CYCLE,[mean(cycle);mean(volume)]] ;
end
end
end
end
    p=polyfit(log(CYCLE(2,:)),log(CYCLE(1,:)),1)                          ;
    plot([0:1:3500],exp(p(2)).*([0:1:3500].^(p(1))),'k--','linewidth',1.5)          ; hold on  ;
for k0=1:2
for k1=1:size(AveCycle{k0,1},1)
for k2=1:size(AveCycle{k0,1},2)
if  AveCycle{k0,1}(k1,k2)~=0 && length(CellVolume{k0,1}{k1,k2})==size(EMBRYO,1) 
    cycle=CellCycle{k0,1}{k1,k2}./cell2mat(EMBRYO(:,20))' ; volume=CellVolume{k0,1}{k1,k2}./cell2mat(EMBRYO(:,15))'            ;
    loglog(mean(volume)+[-1,1]*std(volume),mean(cycle)*[ 1,1] ,'-','linewidth',1.5,'color',c(k0,:)) ; hold on                  ;
    loglog(mean(volume)*[ 1,1],mean(cycle)+[-1,1]*std(cycle)  ,'-','linewidth',1.5,'color',c(k0,:)) ; hold on                  ;
    loglog(mean(volume),mean(cycle),'.','markersize',15,'color',c(k0,:)) ; hold on  ; CYCLE=[CYCLE,[mean(cycle);mean(volume)]] ;
end
end
end
end
    axis([-500,3500,10,100]) ; h=findobj(gca)               ;
    set(gca,'FontSize',14,'Fontname','arial') ; axis square ;
    set(gca,'xtick',[1,10,100,1000],'FontSize',14,'Fontname','arial')    ;
    set(gca,'ytick',[10,100       ],'FontSize',14,'Fontname','arial')    ;
        
    c=hsv(length(AveCycle)) ; figure('Position',[0 0 550 550]) ; CYCLE=[]           ;
    plot([0:1:3500],exp(p(2)).*([0:1:3500].^(p(1))),'k--','linewidth',1.5)          ; hold on     ;
for k0=1:2
for k1=1:size(AveCycle{k0,1},1)
for k2=1:size(AveCycle{k0,1},2)
if  AveCycle{k0,1}(k1,k2)~=0 && length(CellVolume{k0,1}{k1,k2})==size(EMBRYO,1)
    cycle=CellCycle{k0,1}{k1,k2}./cell2mat(EMBRYO(:,20))' ; volume=CellVolume{k0,1}{k1,k2}./cell2mat(EMBRYO(:,15))'            ;
    plot(mean(volume)+[-1,1]*std(volume),mean(cycle)*[ 1,1] ,'-','linewidth',1.5,'color',c(k0,:)) ; hold on                    ;
    plot(mean(volume)*[ 1,1],mean(cycle)+[-1,1]*std(cycle)  ,'-','linewidth',1.5,'color',c(k0,:)) ; hold on                    ;
    plot(mean(volume),mean(cycle),'.','markersize',15,'color',c(k0,:))   ; hold on  ; CYCLE=[CYCLE,[mean(cycle);mean(volume)]] ;
end
end
end
end
    axis([-50,3500,10,70]) ; h=findobj(gca)                 ;
    set(gca,'xtick',[0,1000,2000,3000],'FontSize',18,'Fontname','arial') ;
    set(gca,'ytick',[0,20,40,60      ],'FontSize',18,'Fontname','arial') ;
    x=xlabel({'\rm Cell volume (\itV\rm_ / \mum^{3})'})     ; set(x,'Fontname','arial','Fontsize',18) ;
    y=ylabel({'\rm Cell cycle duration (\itC\rm_ / min)'})  ; set(y,'Fontname','arial','Fontsize',18) ;
    set(gca,'FontSize',18,'Fontname','arial') ; axis square ;

  % C + P3 + P4 Cells
    c=hsv(length(AveCycle)) ; figure('Position',[100 100 200 200])       ; CYCLE=[]                   ;
for k0=[4,6,7]
for k1=1:size(AveCycle{k0,1},1)
for k2=1:size(AveCycle{k0,1},2)
if  AveCycle{k0,1}(k1,k2)~=0 && length(CellVolume{k0,1}{k1,k2})==size(EMBRYO,1)
    cycle=CellCycle{k0,1}{k1,k2}./cell2mat(EMBRYO(:,20))' ; volume=CellVolume{k0,1}{k1,k2}./cell2mat(EMBRYO(:,15))' ; valid=[]  ;
    loglog(mean(volume)+[-1,1]*std(volume),mean(cycle)*[ 1,1] ,'-','linewidth',1.5,'color',c(k0,:))   ; hold on     ;
    loglog(mean(volume)*[ 1,1],mean(cycle)+[-1,1]*std(cycle)  ,'-','linewidth',1.5,'color',c(k0,:))   ; hold on     ;
    loglog(mean(volume),mean(cycle),'.','markersize',15,'color',c(k0,:))   ; hold on ; CYCLE=[CYCLE,[mean(cycle);mean(volume)]] ;
end
end
end
end
    p=polyfit(log(CYCLE(2,:)),log(CYCLE(1,:)),1)                           ;
    plot([0:1:3500],exp(p(2)).*([0:1:3500].^(p(1))),'k--','linewidth',1.5) ; hold on ;
for k0=[4,6,7]
for k1=1:size(AveCycle{k0,1},1)
for k2=1:size(AveCycle{k0,1},2)
if  AveCycle{k0,1}(k1,k2)~=0 && length(CellVolume{k0,1}{k1,k2})==size(EMBRYO,1) 
    cycle=CellCycle{k0,1}{k1,k2}./cell2mat(EMBRYO(:,20))' ; volume=CellVolume{k0,1}{k1,k2}./cell2mat(EMBRYO(:,15))' ; valid=[]  ;
    loglog(mean(volume)+[-1,1]*std(volume),mean(cycle)*[ 1,1] ,'-','linewidth',1.5,'color',c(k0,:)) ; hold on       ;
    loglog(mean(volume)*[ 1,1],mean(cycle)+[-1,1]*std(cycle)  ,'-','linewidth',1.5,'color',c(k0,:)) ; hold on       ;
    loglog(mean(volume),mean(cycle),'.','markersize',15,'color',c(k0,:)) ; hold on ; CYCLE=[CYCLE,[mean(cycle);mean(volume)]]   ;
end
end
end
end
    axis([-500,3500,10,100]) ; h=findobj(gca)                  ;
    set(gca,'FontSize',14,'Fontname','arial') ; axis square    ;
    set(gca,'xtick',[1,10,100,1000],'FontSize',14,'Fontname','arial')    ;
    set(gca,'ytick',[10,100       ],'FontSize',14,'Fontname','arial')    ;
    
    c=hsv(length(AveCycle)) ; figure('Position',[0 0 550 550]) ; CYCLE=[]         ;
    plot([0:1:3500],exp(p(2)).*([0:1:3500].^(p(1))),'k--','linewidth',1.5)        ; hold on       ;
for k0=[4,6,7]
for k1=1:size(AveCycle{k0,1},1)
for k2=1:size(AveCycle{k0,1},2)
if  AveCycle{k0,1}(k1,k2)~=0 && length(CellVolume{k0,1}{k1,k2})==size(EMBRYO,1)
    cycle=CellCycle{k0,1}{k1,k2}./cell2mat(EMBRYO(:,20))' ; volume=CellVolume{k0,1}{k1,k2}./cell2mat(EMBRYO(:,15))' ; valid=[]  ;
    plot(mean(volume)+[-1,1]*std(volume),mean(cycle)*[ 1,1] ,'-','linewidth',1.5,'color',c(k0,:)) ; hold on         ;
    plot(mean(volume)*[ 1,1],mean(cycle)+[-1,1]*std(cycle)  ,'-','linewidth',1.5,'color',c(k0,:)) ; hold on         ;
    plot(mean(volume),mean(cycle),'.','markersize',15,'color',c(k0,:)) ; hold on ; CYCLE=[CYCLE,[mean(cycle);mean(volume)]]     ;
end
end
end
end
    axis([-50,3500,10,70]) ; h=findobj(gca)                ;
    set(gca,'xtick',[0,1000,2000,3000],'FontSize',18,'Fontname','arial') ;
    set(gca,'ytick',[0,20,40,60      ],'FontSize',18,'Fontname','arial') ;
    x=xlabel({'\rm Cell volume (\itV\rm_ / \mum^{3})'})    ; set(x,'Fontname','arial','Fontsize',18) ;
    y=ylabel({'\rm Cell cycle duration (\itC\rm_ / min)'}) ; set(y,'Fontname','arial','Fontsize',18) ;
    set(gca,'FontSize',18,'Fontname','arial') ; axis square              ;