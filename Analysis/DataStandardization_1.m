  % Data Standardization
  
  
  
  % Mission 1 : Timing Setting
    
    File=dir('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200318 - DataUpdate4\RevisionData\LostCells') ; File([1,2,3,12,14],:)=[] ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 2 - WT Position\WorkSpace_Sequence0') ; Sequence0(:,58:end)=[] ; EMBRYO=cell(length(File),1)     ;
for I=1:length(File)
    EMBRYO{I,1}=File(I).name ;
end
    save('WorkSpace_EMBRYO','EMBRYO','-v7.3')  ;
for n=1:27
    tmax=[]     ; n
for I=1:size(EMBRYO,1)
    Data=importdata(['E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200318 - DataUpdate4\RevisionData\RawData\',EMBRYO{I,1},'\CD',EMBRYO{I,1},'.csv']) ;
    N=Data.data ; N=N(:,1)   ; C=Data.textdata ; C=C(2:size(C,1),1) ; K=[0] ; Tmax=zeros(Sequence0{3,n},1)                  ; I
    for k=1:Sequence0{3,n}
    for j=1:size(C,1)
        if strcmp(C{j}(1,1:strfind(C{j},':')-1),Sequence0{5,n}{k,1})==1
           K=[K,N(j,:)]            ;
        end
    end
           Tmax(k,1)=max(K); K=[0] ;
    end
           tmax=[tmax,min(Tmax)]   ;
end
           Sequence0{7,n}=tmax     ;
end
for n=28:57
    tmin=[]     ; n
for I=1:size(EMBRYO,1)
    Data=importdata(['E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200318 - DataUpdate4\RevisionData\RawData\',EMBRYO{I,1},'\CD',EMBRYO{I,1},'.csv']) ;
    N=Data.data ; N=N(:,1) ; C=Data.textdata ; C=C(2:size(C,1),1) ; K=[] ; Tmin=zeros(Sequence0{3,n},1) ; clear Data ; I
    for k=1:Sequence0{3,n}
    for j=1:size(C,1)
        if strcmp(C{j}(1,1:strfind(C{j},':')-1),Sequence0{5,n}{k,1})==1
           K=[K,N(j,:)]            ;
        end
    end
           Tmin(k,1)=min(K) ; K=[] ;
    end
           tmin=[tmin,max(Tmin)]   ;
end
           Sequence0{7,n}=tmin     ;
end
           save('WorkSpace_Sequence0','Sequence0','-v7.3')        ;



  % Mission 2 : Data Merging

    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_CellName') ; Position=CellName ;
    load('WorkSpace_Sequence0') ; load('WorkSpace_Query') ; load('WorkSpace_EMBRYO') ; Sequence0=Sequence0(1:9,:)  ; Sequence0{10,1}=[] ;
for n=1:size(Sequence0,2)
    n
  % Cell Position
    Position=cell(1,length(EMBRYO)) ; Segmentation1=cell(1,length(EMBRYO))           ; Segmentation2=cell(1,length(EMBRYO))             ;
for I=1:size(EMBRYO,1)
    t=Sequence0{7,n}(I)             ; I
    if  t<10
        T=['00',num2str(t)]         ;
    end
    if  t>=10 && t<=99
        T=['0' ,num2str(t)]         ;
    end
    if  t>=100
        T=[     num2str(t)]         ;
    end
    Data=readtable(['E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200318 - DataUpdate4\RevisionData\LostCells\',EMBRYO{I,1},'\',EMBRYO{I,1},'_',T,'_nucLoc.csv']) ;
    A=[] ; B={} ; position=CellName ;
    for i =1:size(Data,1)
    if  strcmp(table2cell(Data(i,8)),'mother')==0
        A=[A;cell2mat(table2cell(Data(i,3:5)))] ; B{end+1,1}=table2cell(Data(i,2)) ;
    end
    end
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
        position{k0,1}{k1,k2}=0                 ;
    if  isnumeric(CellName{k0,1}{k1,k2})==0
        for i=1:length(B)
        if  strcmp(CellName{k0,1}{k1,k2},B{i})==1
            position{k0,1}{k1,k2}=[A(i,3),A(i,1),A(i,2)] ;
        end
        end
    end
    end
    end
    end
        Position{1,I}=position ; clear position          ;
        
    Seg1=load_nii(['E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200318 - DataUpdate4\RevisionData\SegmentedCell\',EMBRYO{I,1},'LabelUnified\',EMBRYO{I,1},'LabelUnified\',EMBRYO{I,1},'_',T,'_segCell.nii']) ;
    Seg1=Seg1.img  ; Cell=unique(Seg1) ; segmentation1=cell(length(Cell)-1,2)           ;
    for i=2:length(Cell)
        [x,y,z]=ind2sub(size(Seg1),find(Seg1==Cell(i)))  ; segmentation1{i-1,1}=[z,x,y] ;
        for k0=1:length(CellName)
        for k1=1:size(CellName{k0,1},1)
        for k2=1:size(CellName{k0,1},2)
        if  isnumeric(CellName{k0,1}{k1,k2})==0 && length(Position{1,I}{k0,1}{k1,k2})==3
        if  Seg1(Position{1,I}{k0,1}{k1,k2}(2),Position{1,I}{k0,1}{k1,k2}(3),Position{1,I}{k0,1}{k1,k2}(1))==Cell(i)
            segmentation1{i-1,2}=CellName{k0,1}{k1,k2}   ;
        end
        end
        end
        end
        end
    end
        Segmentation1{1,I}=segmentation1 ; clear segmentation1 ;
        
end
        Sequence0{8,n}=Position ; Sequence0{9,n}=Segmentation1 ;
end
for n=1:size(Sequence0,2)
for I=1:size(EMBRYO,1)
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence0{8,n}{1,I}{k0,1}{k1,k2})==3
        Sequence0{8,n}{1,I}{k0,1}{k1,k2}     = Sequence0{8,n}{1,I}{k0,1}{k1,k2}.*[0.250,0.250,0.250] ;
        Sequence0{8,n}{1,I}{k0,1}{k1,k2}(:,3)=-Sequence0{8,n}{1,I}{k0,1}{k1,k2}(:,3)                 ;
    end
    end
    end
    end
    for i=1:size(Sequence0{9,n}{1,I},1)
        Sequence0{ 9,n}{1,I}{i,1}     = Sequence0{ 9,n}{1,I}{i,1}.*(ones(size(Sequence0{ 9,n}{1,I}{i,1},1),1)*[0.250,0.250,0.250]) ;
        Sequence0{ 9,n}{1,I}{i,1}(:,3)=-Sequence0{ 9,n}{1,I}{i,1}(:,3) ;
    end
end
end
        save('WorkSpace_Sequence00','Sequence0','-v7.3')               ;
    


  % Mission 3 : Data Input ~ hlh-1 & PHA-4

    File=dir('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200304 - EmbryoAddition\DataStorage\*.csv')                   ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 2 - WT Position\WorkSpace_Sequence0') ; File=[File(1:5,:);File(22:29,:);File(6:21,:)] ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_CellName')     ; SSS=Sequence0 ; SSS(:,58:end)=[]              ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_WTinfo')       ; HLH1=[]       ;
for i=1:length(File)
for I=1:size(WTinfo,1)
if  strcmp(File(i).name,WTinfo{I,1})==1
    HLH1=[HLH1;I] ;
end
end
end
for n=1:size(SSS,2)
    Position1=cell(1,length(HLH1))         ; load(['Merging\WorkSpace_Sequence0_',num2str(n)])     ; n
for i=HLH1'
    Data=importdata(['E:\Project 1 - C.elegans Resource\1-Time Normalization\Wildtype Data-20190131\',WTinfo{i,1}])      ; i
    A=Data.data     ; T=A(:,1)*WTinfo{i,2} ; A=[A(:,7)*WTinfo{i,9},A(:,8:9)*0.09] ; B=Data.textdata ; B=B(2:size(B,1),1) ;
    position1=CellName ; clear Data        ;
    for j=1:length(T)
    if  T(j,1)-WTinfo{i,3}==SSS{7,n}(1,i)
        for k0=1:length(CellName)
        for k1=1:size(CellName{k0,1},1)
        for k2=1:size(CellName{k0,1},2)
        if  strcmp(B{j}(1,1:strfind(B{j},':')-1),CellName{k0,1}{k1,k2})==1
            position1{k0,1}{k1,k2}=A(j,:)  ;
        end
        end
        end
        end
    end
    end
        for k0=1:length(CellName)
        for k1=1:size(CellName{k0,1},1)
        for k2=1:size(CellName{k0,1},2)
        if  isnumeric(cell2mat(position1{k0,1}(k1,k2)))==0 || length(cell2mat(position1{k0,1}(k1,k2)))==1
            position1{k0,1}{k1,k2}=0       ;
        end
        end
        end
        end
            Position1{1,find(HLH1==i)}=position1 ; clear position1 ; clear k0       ; clear k1 ; clear k2 ;       
end
            Sequence0{8,n}=[Sequence0{8,n},Position1]              ;
            clear Position1 ; clear Position2 ; clear Position3 ; clear A ; clear B ; clear T  ; clear i  ;
            save(['Addition\WorkSpace_Sequence00_',num2str(n)],'Sequence0','-v7.3') ;
end