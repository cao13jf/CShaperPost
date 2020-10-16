  % Mission 1 : Data Output of Figure 3d
  
    load('D:\STORAGE\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUSPLUS\Revision45 - Fig3de_Update (Label)\WorkSpace_EMBRYO'    ) ;
    load('D:\STORAGE\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUSPLUS\Revision45 - Fig3de_Update (Label)\WorkSpace_CellVolume') ;
    load('D:\STORAGE\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUSPLUS\Revision45 - Fig3de_Update (Label)\WorkSpace_Query'     ) ;
    Max=[] ; CELLNAME={}   ; c=hsv(size(EMBRYO,1)) ; SourceData1={} ; SourceData2={} ; SourceData3={}   ;
    SourceData1{1,1}='Cell Volume Before Normalization (¦Ìm3)'      ; SourceData1{2,1}='Cell Name'      ; CELLNAME={} ;
    SourceData2{1,1}='Cell Volume After Normalization (¦Ìm3)'       ; SourceData2{2,1}='Cell Name'      ; VOLUME=[]   ;
    SourceData3{1,1}='Variation Coefficient of Cell Volume (Normalized)' ; SourceData3{2,1}='Cell Name' ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_CellName')          ;
for I =1:size(EMBRYO,1)
    volume=[] ; NUM=2      ; I
if  I+3<10
    III=['0',num2str(I+3)] ;
end
if  I+3>=10
    III=[    num2str(I+3)] ;
end
    SourceData1{2,1+I}=['Sample ',III] ; SourceData2{2,1+I}=['Sample ',III]    ;   
    
for k0=1
for k1=1:size(CellName{k0,1},1)
for k2=2:size(CellName{k0,1},2)-1
if  isnumeric(CellName{k0,1}{k1,k2})==0 && length(CellVolume{k0,1}{k1,k2})==size(EMBRYO,1)
    volume=[volume,[mean(CellVolume{k0,1}{k1,k2});CellVolume{k0,1}{k1,k2}(I)]] ;
    NUM=NUM+1 ; SourceData1{NUM,1}=CellName{k0,1}{k1,k2} ; SourceData2{NUM,1}=CellName{k0,1}{k1,k2} ; SourceData3{NUM,1}=CellName{k0,1}{k1,k2} ;
end
end
end
end
for k0=[2,3,4,5,7]
for k1=1:size(CellName{k0,1},1)
for k2=1:size(CellName{k0,1},2)-1
if  isnumeric(CellName{k0,1}{k1,k2})==0 && length(CellVolume{k0,1}{k1,k2})==size(EMBRYO,1)
    volume=[volume,[mean(CellVolume{k0,1}{k1,k2});CellVolume{k0,1}{k1,k2}(I)]] ;
    NUM=NUM+1 ; SourceData1{NUM,1}=CellName{k0,1}{k1,k2} ; SourceData2{NUM,1}=CellName{k0,1}{k1,k2} ; SourceData3{NUM,1}=CellName{k0,1}{k1,k2} ;
end
end
end
end
    volume=[volume,[mean(CellVolume{ 6,1}{ 1, 1});CellVolume{ 6,1}{ 1, 1}(I)]] ;
    NUM=NUM+1 ; SourceData1{NUM,1}=CellName{6,1}{1,1}    ; SourceData2{NUM,1}=CellName{6,1}{1,1}    ; SourceData3{NUM,1}=CellName{6,1}{1,1}    ;
    EMBRYO{I,14}=volume ;
end
for I=1:size(EMBRYO,1)
    b =sum(EMBRYO{I,14}(2,:).*EMBRYO{I,14}(1,:))/sum(EMBRYO{I,14}(1,:).^2)        ;
    yy=sum(EMBRYO{I,14}(2,:).*EMBRYO{I,14}(1,:))/sum(EMBRYO{I,14}(1,:).^2)*EMBRYO{I,14}(1,:)  ; y=EMBRYO{I,14}(2,:) ;
    SSE=sum((yy-y).^2)   ; SST=sum((y-mean(y)).^2) ; EMBRYO{I,15}=b ; EMBRYO{I,16}=1-SSE/SST  ; clear cycle_i       ;
    clear k0 ; clear k1  ; clear k2  ; clear nn    ; clear b ; clear yy ; clear y ; clear SSE ; clear cycle_j       ;
end
for I=1:size(EMBRYO,1)
    for iii=1:length(EMBRYO{I,14}(2,:))
        SourceData1{2+iii,1+I}=sprintf('%2.2f',EMBRYO{I,14}(2,iii)) ;
        SourceData2{2+iii,1+I}=sprintf('%2.2f',EMBRYO{I,14}(2,iii)/EMBRYO{I,15})  ;
    end
end
for k1=3:size(SourceData1,1)
for k2=2:size(SourceData1,2)
if  isempty(SourceData1{k1,k2})==1
    k1
    k2
end
end
end
for I=1:size(EMBRYO,1)
    VOLUME=[VOLUME;EMBRYO{I,14}(2,:)/EMBRYO{I,15}]             ;
end
    VAR=std(VOLUME)./mean(VOLUME) ;
for iii=1:length(VAR)
    SourceData3{2+iii,2}=sprintf('%2.4f',VAR(1,iii))           ;
end
    xlswrite('Figure 3d.xls',SourceData1,'Main Graph'        ) ;
    xlswrite('Figure 3d.xls',SourceData2,'Top-Left Inset'    ) ;
    xlswrite('Figure 3d.xls',SourceData3,'Bottom-Right Inset') ; 