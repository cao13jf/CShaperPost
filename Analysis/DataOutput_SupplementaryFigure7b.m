  % Mission 1 : Data Output of Supplementary Figure 7b

    load('D:\STORAGE\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUSPLUS\Revision61 - SizeAsymmetry_PLUSPLUSPLUS\WorkSpace_EMBRYO'    ) ;
    load('D:\STORAGE\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUSPLUS\Revision61 - SizeAsymmetry_PLUSPLUSPLUS\WorkSpace_CellVolume') ;
    SourceData={}     ; SourceData{1,1}='Volume Ratio Between Sister Cells' ; SourceData{2,1}='Cell 1'   ; SourceData{2,2}='Cell 2'     ;
    SourceDataPLUS={} ; SourceDataPLUS{1,3}='Variation Coefficient of Volume Ratio Between Sister Cells' ; SourceDataPLUS{1,1}='Cell 1' ; SourceDataPLUS{1,2}='Cell 2' ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200323 - MissingCell\WorkSpace_MissingCell') ; 
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200324 - LostRatio\WorkSpace_CellCycle'    ) ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_AveCycle') ; VOLUME=[]                    ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_WTinfo'  ) ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_CellName') ;
for I =1:size(EMBRYO,1)
    volume=[] ; NUM=2                 ; I
if  I+3<10
    III=['0',num2str(I+3)]            ;
end
if  I+3>=10
    III=[    num2str(I+3)]            ;
end
    SourceData{2,2+I}=['Sample ',III] ;
for k0=1
for k1=1:size(CellName{k0,1},1)
for k2=1:size(CellName{k0,1},2)-2
if  isnumeric(CellName{k0,1}{k1,k2})==0                   && length(CellVolume{k0,1}{k1,k2})==size(EMBRYO,1)
if  length(CellVolume{k0,1}{2*k1-1,k2+1})==size(EMBRYO,1) && length(CellVolume{k0,1}{2*k1,k2+1})==size(EMBRYO,1)
    if  mean(CellVolume{k0,1}{2*k1,k2+1})/mean(CellVolume{k0,1}{2*k1-1,k2+1})>1
        volume=[volume,[mean(CellVolume{k0,1}{2*k1,k2+1})/mean(CellVolume{k0,1}{2*k1-1,k2+1});CellVolume{k0,1}{2*k1,k2+1}(I)/CellVolume{k0,1}{2*k1-1,k2+1}(I)]]       ;
        SourceDataPLUS{NUM,1}=CellName{k0,1}{2*k1,k2+1}           ; SourceDataPLUS{NUM,2}=CellName{k0,1}{2*k1-1,k2+1} ;
        NUM=NUM+1 ; SourceData{NUM,1}=CellName{k0,1}{2*k1,k2+1}   ; SourceData{NUM,2}=CellName{k0,1}{2*k1-1,k2+1}     ;
    end
    if  mean(CellVolume{k0,1}{2*k1,k2+1})/mean(CellVolume{k0,1}{2*k1-1,k2+1})<1
        volume=[volume,[mean(CellVolume{k0,1}{2*k1,k2+1})/mean(CellVolume{k0,1}{2*k1-1,k2+1});CellVolume{k0,1}{2*k1,k2+1}(I)/CellVolume{k0,1}{2*k1-1,k2+1}(I)].^(-1)] ;
        SourceDataPLUS{NUM,1}=CellName{k0,1}{2*k1-1,k2+1}         ; SourceDataPLUS{NUM,2}=CellName{k0,1}{2*k1,k2+1}   ;
        NUM=NUM+1 ; SourceData{NUM,1}=CellName{k0,1}{2*k1-1,k2+1} ; SourceData{NUM,2}=CellName{k0,1}{2*k1,k2+1}       ;
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
        SourceDataPLUS{NUM,1}=CellName{k0,1}{2*k1,k2+1}           ; SourceDataPLUS{NUM,2}=CellName{k0,1}{2*k1-1,k2+1} ;
        NUM=NUM+1 ; SourceData{NUM,1}=CellName{k0,1}{2*k1,k2+1}   ; SourceData{NUM,2}=CellName{k0,1}{2*k1-1,k2+1}     ;
    end
    if  mean(CellVolume{k0,1}{2*k1,k2+1})/mean(CellVolume{k0,1}{2*k1-1,k2+1})<1
        volume=[volume,[mean(CellVolume{k0,1}{2*k1,k2+1})/mean(CellVolume{k0,1}{2*k1-1,k2+1});CellVolume{k0,1}{2*k1,k2+1}(I)/CellVolume{k0,1}{2*k1-1,k2+1}(I)].^(-1)] ;
        SourceDataPLUS{NUM,1}=CellName{k0,1}{2*k1-1,k2+1}         ; SourceDataPLUS{NUM,2}=CellName{k0,1}{2*k1,k2+1}   ;
        NUM=NUM+1 ; SourceData{NUM,1}=CellName{k0,1}{2*k1-1,k2+1} ; SourceData{NUM,2}=CellName{k0,1}{2*k1,k2+1}       ;
    end
end
end
end
end
end
    if  mean(CellVolume{2,1}{1,1})/mean(CellVolume{3,1}{1,1})>1
        volume=[volume,[mean(CellVolume{2,1}{1,1})/mean(CellVolume{3,1}{1,1});CellVolume{2,1}{1,1}(I)/CellVolume{3,1}{1,1}(I)]]       ;
        SourceDataPLUS{NUM,1}=CellName{2,1}{1,1}         ; SourceDataPLUS{NUM,2}=CellName{3,1}{1,1} ;
        NUM=NUM+1 ; SourceData{NUM,1}=CellName{2,1}{1,1} ; SourceData{NUM,2}=CellName{3,1}{1,1}     ;
    end
    if  mean(CellVolume{2,1}{1,1})/mean(CellVolume{3,1}{1,1})<1
        volume=[volume,[mean(CellVolume{2,1}{1,1})/mean(CellVolume{3,1}{1,1});CellVolume{2,1}{1,1}(I)/CellVolume{3,1}{1,1}(I)].^(-1)] ;
        SourceDataPLUS{NUM,1}=CellName{3,1}{1,1}         ; SourceDataPLUS{NUM,2}=CellName{2,1}{1,1} ;
        NUM=NUM+1 ; SourceData{NUM,1}=CellName{3,1}{1,1} ; SourceData{NUM,2}=CellName{2,1}{1,1}     ;
    end
        
    if  mean(CellVolume{4,1}{1,1})/mean(CellVolume{6,1}{1,1})>1
        volume=[volume,[mean(CellVolume{4,1}{1,1})/mean(CellVolume{6,1}{1,1});CellVolume{4,1}{1,1}(I)/CellVolume{6,1}{1,1}(I)]]       ;
        SourceDataPLUS{NUM,1}=CellName{4,1}{1,1}         ; SourceDataPLUS{NUM,2}=CellName{6,1}{1,1} ;
        NUM=NUM+1 ; SourceData{NUM,1}=CellName{4,1}{1,1} ; SourceData{NUM,2}=CellName{6,1}{1,1}     ;
    end
    if  mean(CellVolume{4,1}{1,1})/mean(CellVolume{6,1}{1,1})<1
        volume=[volume,[mean(CellVolume{4,1}{1,1})/mean(CellVolume{6,1}{1,1});CellVolume{4,1}{1,1}(I)/CellVolume{6,1}{1,1}(I)].^(-1)] ;
        SourceDataPLUS{NUM,1}=CellName{6,1}{1,1}         ; SourceDataPLUS{NUM,2}=CellName{4,1}{1,1} ;
        NUM=NUM+1 ; SourceData{NUM,1}=CellName{6,1}{1,1} ; SourceData{NUM,2}=CellName{4,1}{1,1}     ;
    end

    if  mean(CellVolume{5,1}{1,1})/mean(CellVolume{7,1}{1,1})>1
        volume=[volume,[mean(CellVolume{5,1}{1,1})/mean(CellVolume{7,1}{1,1});CellVolume{5,1}{1,1}(I)/CellVolume{7,1}{1,1}(I)]]       ;
        SourceDataPLUS{NUM,1}=CellName{5,1}{1,1}         ; SourceDataPLUS{NUM,2}=CellName{7,1}{1,1} ;
        NUM=NUM+1 ; SourceData{NUM,1}=CellName{5,1}{1,1} ; SourceData{NUM,2}=CellName{7,1}{1,1}     ;
    end
    if  mean(CellVolume{5,1}{1,1})/mean(CellVolume{7,1}{1,1})<1
        volume=[volume,[mean(CellVolume{5,1}{1,1})/mean(CellVolume{7,1}{1,1});CellVolume{5,1}{1,1}(I)/CellVolume{7,1}{1,1}(I)].^(-1)] ;
        SourceDataPLUS{NUM,1}=CellName{7,1}{1,1}         ; SourceDataPLUS{NUM,2}=CellName{5,1}{1,1} ;
        NUM=NUM+1 ; SourceData{NUM,1}=CellName{7,1}{1,1} ; SourceData{NUM,2}=CellName{5,1}{1,1}     ;
    end
        EMBRYO{I,14}=volume                              ;
end
for I=1:size(EMBRYO,1)
    for iii=1:length(EMBRYO{I,14}(2,:))
        SourceData{2+iii,2+I}=sprintf('%2.4f',EMBRYO{I,14}(2,iii))      ;
    end
end
for I=1:size(EMBRYO,1)
    VOLUME=[VOLUME;EMBRYO{I,14}(2,:)]                    ;
end
    VAR=std(VOLUME)./mean(VOLUME)                        ;
for iii=1:length(VAR)
    SourceDataPLUS{1+iii,3}=sprintf('%2.4f',VAR(1,iii))  ;
end
    xlswrite('Supplementary Figure 7b.xls',SourceDataPLUS,'Figure S7b') ;