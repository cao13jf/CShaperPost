  % Mission 1 : Data Output of Supplementary Figure 10a

    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200331 - Volume&Surface_PLUS\WorkSpace_EMBRYO') ; embryo=EMBRYO  ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200324 - LostRatio\WorkSpace_EMBRYO')           ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_CellName')   ;
for k1=1:size(embryo,1)
for k2=1:size(embryo,2)
if  isempty(EMBRYO{k1,k2})==1 && isempty(embryo{k1,k2})==0
    EMBRYO{k1,k2}=embryo{k1,k2} ;
end
end
end
    save('WorkSpace_EMBRYO','EMBRYO','-v7.3') ; clear k1 ; clear k2 ; clear embryo               ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200331 - Volume&Surface_PLUS\WorkSpace_CellVolume')     ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200323 - MissingCell\WorkSpace_MissingCell')            ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200324 - LostRatio\WorkSpace_CellCycle')                ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_AveCycle')   ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_WTinfo'  )   ;
    SourceData1={} ; SourceData1{1,1}='Cell Name' ; SourceData1{1,2}='Cell Cycle Duration (min)' ; SourceData1{1,3}='Cell Volume (¦Ìm3)' ; NUM1=1 ;
    SourceData2={} ; SourceData2{1,1}='Cell Name' ; SourceData2{1,2}='Cell Cycle Duration (min)' ; SourceData2{1,3}='Cell Volume (¦Ìm3)' ; NUM2=1 ;
    
  % AB + MS Cells
    c=hsv(length(AveCycle)) ; CYCLE=[]            ;
for k0=1:2
for k1=1:size(AveCycle{k0,1},1)
for k2=1:size(AveCycle{k0,1},2)
if  AveCycle{k0,1}(k1,k2)~=0 && length(CellVolume{k0,1}{k1,k2})==size(EMBRYO,1)
    cycle=CellCycle{k0,1}{k1,k2}./cell2mat(EMBRYO(:,20))'   ; volume=CellVolume{k0,1}{k1,k2}./cell2mat(EMBRYO(:,15))' ;
    NUM1=NUM1+1 ; SourceData1{NUM1,1}=CellName{k0,1}{k1,k2} ;
    SourceData1{NUM1,2}=sprintf('%2.2f',mean(cycle))        ; SourceData1{NUM1,3}=sprintf('%2.2f',mean(volume))       ;
end
end
end
end

  % C + P3 + P4 Cells
    c=hsv(length(AveCycle)) ; CYCLE=[]                      ;
for k0=[4,6,7]
for k1=1:size(AveCycle{k0,1},1)
for k2=1:size(AveCycle{k0,1},2)
if  AveCycle{k0,1}(k1,k2)~=0 && length(CellVolume{k0,1}{k1,k2})==size(EMBRYO,1)
    cycle=CellCycle{k0,1}{k1,k2}./cell2mat(EMBRYO(:,20))'   ; volume=CellVolume{k0,1}{k1,k2}./cell2mat(EMBRYO(:,15))' ; valid=[] ;
    NUM2=NUM2+1 ; SourceData2{NUM2,1}=CellName{k0,1}{k1,k2} ;
    SourceData2{NUM2,2}=sprintf('%2.2f',mean(cycle))        ; SourceData2{NUM2,3}=sprintf('%2.2f',mean(volume))                  ;
end
end
end
end
    xlswrite('Supplementary Figure 10a.xls',SourceData1,'Figure S10a') ;