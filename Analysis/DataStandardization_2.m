  % Data Standardization

  
  
  % Mission 1 : Existence Probability Of Cell

    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_CellName') ; load('WorkSpace_EMBRYO') ;
    load('WorkSpace_Sequence0') ; S=Sequence0   ; load('WorkSpace_HLH1') ; Criteria=cell(5,size(S,2))                     ;
for n=1:size(Sequence0,2)
    load(['Addition\WorkSpace_Sequence00_',num2str(n)]) ; Positive=CellName                    ; n
    S{12,n}=cell(1,size(EMBRYO,1)+length(HLH1)) ; S{13,n}=cell(1,size(EMBRYO,1)+length(HLH1))  ;
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
        Positive{k0,1}{k1,k2}=0 ;
        for I=1:length(Sequence0{8,n})
        if  length(Sequence0{8,n}{1,I}{k0,1}{k1,k2})==3
            Positive{k0,1}{k1,k2}=Positive{k0,1}{k1,k2}+1/length(Sequence0{8,n})               ;
        end
        end
    end
    end
    end
            Criteria{1,n}=Positive ; clear Positive                      ;
end
            save('Normalization1\WorkSpace_Criteria','Criteria','-v7.3') ;



  % Mission 2 : Cavity Importing

    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200322 - Cavity\WorkSpace_EMBRYO'       ) ; embryo=EMBRYO         ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200319 - DataImport\WorkSpace_EMBRYO'   ) ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200319 - DataImport\WorkSpace_Sequence0') ; SSS=Sequence0         ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200322 - Cavity\WorkSpace_Mission'      ) ;
    File=dir('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200405 - CavityResult\Result*'      ) ;
for n=1:57
    load(['E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200319 - DataImport\Addition\WorkSpace_Sequence00_',num2str(n)]) ; n
for I=1:size(EMBRYO,1)
for i=1:size(embryo,1)
if  strcmp(EMBRYO{I,1},embryo{i,1})==1
    t=SSS{7,n}(1,I)         ; I
    if  t<10
        T=['00',num2str(t)] ;
    end
    if  t>=10 && t<=99
        T=['0' ,num2str(t)] ;
    end
    if  t>=100
        T=[     num2str(t)] ;
    end
    for j=1:length(File)
        file=dir(['E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200405 - CavityResult\',File(j).name,'\Embryo',num2str(i),'\*.mat'])     ;
    for J=1:length(file)
    if  strcmp(['Seg_',num2str(i),'_',T,'.mat'],file(J).name)==1
        load(['E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200405 - CavityResult\',File(j).name,'\Embryo',num2str(i),'\',file(J).name]) ;
        [x,y,z]=ind2sub(size(Seg),find(Seg==0)) ;
    if  isempty([z,x,y])==0
        Sequence0{9,n}{1,I}{end+1,1}=[z,x,y]    ; Sequence0{9,n}{1,I}{end,2}='Cavity' ;
        Sequence0{9,n}{1,I}{end,1}=Sequence0{9,n}{1,I}{end,1}.*(ones(size(Sequence0{9,n}{1,I}{end,1},1),1)*[0.250,0.250,0.250]) ;
    end
    end
    end
    end
end
end
end
        save(['Cavity\WorkSpace_Sequence000_',num2str(n)],'Sequence0','-v7.3')        ;
end



  % Mission 3 : Cavity Data Merging

    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_CellName')      ; load('WorkSpace_EMBRYO')         ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 2 - WT Position\WorkSpace_Sequence0')  ; SSS=Sequence0 ; SSS(:,58:end)=[] ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_WTinfo')        ; load('WorkSpace_HLH1')           ;
    load('Normalization1\WorkSpace_Criteria')          ;
for i=1:(size(EMBRYO,1)+length(HLH1))
    S=SSS ; S(11,:)=[]     ; i
    
  % Merging Data of Embryo i
for n=1:size(S,2)
    load(['Cavity\WorkSpace_Sequence000_',num2str(n)]) ; Sequence0{9,n}=[Sequence0{9,n},cell(1,length(HLH1))]  ; n
for j=1:length(EMBRYO)+length(HLH1)
if  i~=j
    Sequence0{8,n}{1,j}=[] ; Sequence0{9,n}{1,j}=[]    ;
end
end
    S(:,n)=Sequence0(:,n)  ; clear Sequence0 ; clear n ; clear j ;
end
    Sequence0=S ; save(['Normalization1\WorkSpace_Sequence0_Embryo',num2str(i)],'Sequence0','-v7.3') ; clear S ; clear Sequence0       ;
end