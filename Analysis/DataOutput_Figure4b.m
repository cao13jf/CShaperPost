  % Mission 1 : Data Output of Figure 4e

    load('D:\STORAGE\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUSPLUS\Revision45 - Fig3de_Update (Label)\WorkSpace_EMBRYO'    ) ;
    load('D:\STORAGE\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUSPLUS\Revision45 - Fig3de_Update (Label)\WorkSpace_CellVolume') ;
    load('WorkSpace_N') ; SourceData={} ; SourceData{1,1}='Cell Name'  ; SourceData{1,2}='Cell Cycle Duration (min)'             ;
    SourceData{1,3}='Dynamics of Cell Shape Irregularity ¦Ç (Averaged & Interpolated)' ; NUM=1 ; load('WorkSpace_Tree')          ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_CellName') ;
    SourceData{2,1}='ABa' ; SourceData{2,2}=sprintf('%2.2f',Tree{1,5}) ; SourceData{2,3}=sprintf('%2.4f',Tree{1,7}) ; NUM=NUM+1  ;
    SourceData{3,1}='ABp' ; SourceData{3,2}=sprintf('%2.2f',Tree{2,5}) ; SourceData{3,3}=sprintf('%2.4f',Tree{2,7}) ; NUM=NUM+1  ;
    SourceData{4,1}='EMS' ; SourceData{4,2}=sprintf('%2.2f',Tree{3,5}) ; SourceData{4,3}=sprintf('%2.4f',Tree{3,7}) ; NUM=NUM+1  ;
    SourceData{5,1}='P2'  ; SourceData{5,2}=sprintf('%2.2f',Tree{4,5}) ; SourceData{5,3}=sprintf('%2.4f',Tree{4,7}) ; NUM=NUM+1  ;
for k0=1
for k1=1:size(CellName{k0,1},1)
for k2=2:size(CellName{k0,1},2)-1
if  isnumeric(CellName{k0,1}{k1,k2})==0 && length(CellVolume{k0,1}{k1,k2})==size(EMBRYO,1)
    NUM=NUM+1 ; SourceData{NUM,1}=CellName{k0,1}{k1,k2}  ; SourceData{NUM,2}=sprintf('%2.2f',N{k0,2}(k1,k2)) ; ETA=''            ;
    for i=1:length(N{k0,4}{k1,k2})
        ETA=[ETA,sprintf('%2.4f',N{k0,4}{k1,k2}(i)),','] ;
    end
        SourceData{NUM,3}=ETA(1:end-1)                   ;
end
end
end
end
for k0=[2,3,4,5,7]
for k1=1:size(CellName{k0,1},1)
for k2=1:size(CellName{k0,1},2)-1
if  isnumeric(CellName{k0,1}{k1,k2})==0 && length(CellVolume{k0,1}{k1,k2})==size(EMBRYO,1)
    NUM=NUM+1 ; SourceData{NUM,1}=CellName{k0,1}{k1,k2}  ; SourceData{NUM,2}=sprintf('%2.2f',N{k0,2}(k1,k2)) ; ETA=''            ;
    for i=1:length(N{k0,4}{k1,k2})
        ETA=[ETA,sprintf('%2.4f',N{k0,4}{k1,k2}(i)),','] ;
    end
        SourceData{NUM,3}=ETA(1:end-1)                   ;
end
end
end
end
    NUM=NUM+1 ; SourceData{NUM,1}=CellName{6,1}{1,1}     ; SourceData{NUM,2}=sprintf('%2.2f',N{6,2}(1,1))    ; ETA=''            ;
for i=1:length(N{6,4}{1,1})
    ETA=[ETA,sprintf('%2.4f',N{6,4}{1,1}(i)),',']        ;
end
    SourceData{NUM,3}=ETA(1:end-1)                       ;
    xlswrite('Figure 4b.xls',SourceData,'Figure 4b')     ;