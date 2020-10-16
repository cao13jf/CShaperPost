  % Mission 1 : Data Standardization of Cell-Cell Contact

  % Missing Cell
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200323 - MissingCell\WorkSpace_MissingCell') ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200319 - DataImport\WorkSpace_Query' )       ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200319 - DataImport\WorkSpace_EMBRYO')       ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_CellName') ; missingcell={}               ;
for I =1:size(EMBRYO,1)
for k0=1:length(CellName)
for k1=1:size(CellName{k0,1},1)
    k1
for k2=1:size(CellName{k0,1},2)
if  MissingCell{I,1}{k0,1}{k1,k2}==1
    for i=1:size(Query,1)
    if  strcmp(CellName{k0,1}{k1,k2},Query{i,1})==1
    if  isempty(missingcell)==1
        missingcell{end+1,1}=CellName{k0,1}{k1,k2} ; missingcell{end,2}=Query{i,2}             ;
    end
    if  isempty(missingcell)==0
        NNN=0 ;
    for j=1:size(missingcell,1)
    if  strcmp(CellName{k0,1}{k1,k2},missingcell{j,1})==1
        NNN=1 ;
    end
    end
    if  NNN==0
        missingcell{end+1,1}=CellName{k0,1}{k1,k2} ; missingcell{end,2}=Query{i,2}             ;
    end
    end
    end
    end
end
end
end
end
end
        MissingCell=missingcell   ; save('WorkSpace_MissingCell','MissingCell','-v7.3')        ;
  
    load('WorkSpace_MissingCell') ; load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_CellName') ; Num=0 ;
for k0=1:length(CellName)
for k1=1:size(CellName{k0,1},1)
for k2=1:size(CellName{k0,1},2)
if  isnumeric(CellName{k0,1}{k1,k2})==0
    for i=1:size(MissingCell)
    if  strcmp(CellName{k0,1}{k1,k2},MissingCell{i,1})==1
        CellName{k0,1}{k1,k2}=0   ;
    end
    end
if  isnumeric(CellName{k0,1}{k1,k2})==0
    Num=Num+1 ;
end
end
end
end
end
        
  % Checking for Signaling Cell Pair 
    load('WorkSpace_MissingCell')              ;
    ContactPair={{'P2'},{'EMS'};{'MS'},{'ABal'};{'C'},{'ABar'};{'P2'},{'ABp'};{'MS'},{'ABalp'};{'MS'},{'ABara'};{'ABalapa'},{'ABplaaa'};{'ABalapp'},{'ABplaaa'};{'MSapp'},{'ABplpapp'};{'MSappp'},{'ABplpppp'}} ;
for k1=1:size(ContactPair,1)
for k2=1:size(ContactPair,2)
for i =1:size(MissingCell,1)
if  strcmp(ContactPair{k1,k2},MissingCell{i,1})==1
    MissingCell{i,1}
end
end
end
end

  % Data Import
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200319 - DataImport\WorkSpace_Sequence0') ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200319 - DataImport\WorkSpace_Query'    ) ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200319 - DataImport\WorkSpace_EMBRYO'   ) ;
    load('WorkSpace_MissingCell') ; embryo=EMBRYO      ;
for I=1:size(embryo,1)
    I
  % Generation of WorkSpace_EMBRYO_I
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200319 - DataImport\WorkSpace_EMBRYO'   ) ;
for I=1:size(EMBRYO,1)
    EMBRYO{I,11}=[]               ;
for n=1:size(Sequence0,2)
    EMBRYO{I,11}=[EMBRYO{I,11},Sequence0{7,n}(1,I)]    ;
end
    EMBRYO{I,11}=[min(EMBRYO{I,11}):max(EMBRYO{I,11})] ;
end

  % Importing Contact Information
    contact=importdata(['E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200318 - DataUpdate4\RevisionData\ContactSurface\',EMBRYO{I,1},'_Stat.csv']) ;
    A=contact.data ; B=contact.textdata ; Dot=strfind(B{1,1},',') ;
for j=2:size(B,2)-1
    B{1,j}=B{1,1}(Dot(j-1)+1:Dot(j)-1)  ; B{5,j}=[]               ;
end
    B{1,end}=B{1,1}(Dot(end)+1:end)     ; B{1,1}='cell1'          ;
for j=2:size(B,2)
    I
    j
for i=1:size(Query,1)
if  strcmp(B{1,j},Query{i,1})==1
    B{3,j}=Query{i,2}          ;
end
if  strcmp(B{2,j},Query{i,1})==1
    B{4,j}=Query{i,2}          ;
end
end
if  isempty([find(cell2mat(MissingCell(:,2))==B{3,j}),find(cell2mat(MissingCell(:,2))==B{4,j})])==1
    B{5,j}=NaN                 ;
    for t=EMBRYO{I,11}
    if  t<10
        ttt=['00',num2str(t)]  ;
    end
    if  t>=10 && t<=99
        ttt=['0' ,num2str(t)]  ;
    end
    if  t>=100
        ttt=[     num2str(t)]  ;
    end
    if         (isnan(A(t,j))==0 && A(t,j)~=0)     && isnan(B{5,j})==1
        Data=table2cell(readtable(['E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200318 - DataUpdate4\RevisionData\LostCells\',EMBRYO{I,1},'\',EMBRYO{I,1},'_',ttt,'_nucLoc.csv'])) ;
        surface1=Data{find(cell2mat(Data(:,1))==B{3,j}),7} ; surface2=Data{find(cell2mat(Data(:,1))==B{4,j}),7} ;
        B{5,j}=0               ;
    if  A( t,j)/surface1>=1/48    || A( t,j)/surface2>=1/48
        B{5,j}=1               ;
    end
    end
    if  t>1 && (isnan(A(t-1,j))==0 && A(t-1,j)~=0) && (isnan(A(t,j))==0 && A(t,j)~=0)
        T=0                    ;
    for tt=t-1:t
    if  tt<10
        ttt=['00',num2str(tt)] ;
    end
    if  tt>=10 && tt<=99
        ttt=['0' ,num2str(tt)] ;
    end
    if  tt>=100
        ttt=[     num2str(tt)] ;
    end
        Data=table2cell(readtable(['E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200318 - DataUpdate4\RevisionData\LostCells\',EMBRYO{I,1},'\',EMBRYO{I,1},'_',ttt,'_nucLoc.csv'])) ;
        surface1=Data{find(cell2mat(Data(:,1))==B{3,j}),7} ; surface2=Data{find(cell2mat(Data(:,1))==B{4,j}),7} ;
    if  A(tt,j)/surface1>=1/48    || A(tt,j)/surface2>=1/48
        T=T+1                  ;
    end
    end
    if  T==1 && (isnan(B{5,j})==1 || B{5,j}==0 || B{5,j}==1)
        B{5,j}=1               ;
    end
    if  T==2 && (isnan(B{5,j})==1 || B{5,j}==0 || B{5,j}==1)
        B{5,j}=2               ;
    end
    end
    end
end
end
        EMBRYO{I,13}=B ; save(['WorkSpace_embryo_',num2str(I)],'EMBRYO','-v7.3') ;
end

  % Data Merging
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200319 - DataImport\WorkSpace_EMBRYO') ; embryo=EMBRYO ;
for I=1:size(EMBRYO,1)
    load(['WorkSpace_embryo_',num2str(I)]) ; embryo{I,11}=EMBRYO{I,11} ; embryo{I,13}=EMBRYO{I,13} ; I
end
    EMBRYO=embryo ; save('WorkSpace_EMBRYO','EMBRYO','-v7.3') ; clear  ; clc     ;
     
  % Connection Criteria
  % Num1 ~ Connection ; % Num2 ~ Surface 1/48 ; % Num3 ~ 2 TP % Num4 ~ Reproducible
    load('WorkSpace_EMBRYO') ; Num=[]      ; Connection=[]    ; CONNECTION={}    ;
for I=1:size(EMBRYO,1)
for i=1:size(EMBRYO{I,13},2)
if  isempty(EMBRYO{I,13}{5,i})==0 && isnan(EMBRYO{I,13}{5,i})==0
    X=cell2mat(EMBRYO{I,13}(3:4,i))'       ; Y=[X(2),X(1)]    ;
    if  isempty(Connection)==0
        x=Connection(:,1:2)-ones(size(Connection,1),1)*X      ; y=Connection(:,1:2)-ones(size(Connection,1),1)*Y ;
        z=[intersect(find(x(:,1)==0),find(x(:,2)==0)),intersect(find(y(:,1)==0),find(y(:,2)==0))]                ;
    if  isempty(z)==1
        Connection=[Connection;X]          ;
    end
    end
    if  isempty(Connection)==1
        Connection=[Connection;X]          ;
    end
end
end
end
        CONNECTION{end+1,1}=Connection     ; save('WorkSpace_CONNECTION','CONNECTION','-v7.3')                  ;
        
     Connection=[] ; 
for I=1:size(EMBRYO,1)
for i=1:size(EMBRYO{I,13},2)
if  isempty(EMBRYO{I,13}{5,i})==0 && isnan(EMBRYO{I,13}{5,i})==0 && EMBRYO{I,13}{5,i}>0
    X=cell2mat(EMBRYO{I,13}(3:4,i))'       ; Y=[X(2),X(1)] ;
    if  isempty(Connection)==0
        x=Connection(:,1:2)-ones(size(Connection,1),1)*X   ; y=Connection(:,1:2)-ones(size(Connection,1),1)*Y   ;
        z=[intersect(find(x(:,1)==0),find(x(:,2)==0)),intersect(find(y(:,1)==0),find(y(:,2)==0))]               ;
    if  isempty(z)==1
        Connection=[Connection;X]          ;
    end
    end
    if  isempty(Connection)==1
        Connection=[Connection;X]          ;
    end
end
end
end
    CONNECTION{end+1,1}=Connection         ; Connection=[] ;
for I=1:size(EMBRYO,1)
for i=1:size(EMBRYO{I,13},2)
if  isempty(EMBRYO{I,13}{5,i})==0 && isnan(EMBRYO{I,13}{5,i})==0 && EMBRYO{I,13}{5,i}>1
    X=cell2mat(EMBRYO{I,13}(3:4,i))'       ; Y=[X(2),X(1)] ;
    if  isempty(Connection)==0
        x=Connection(:,1:2)-ones(size(Connection,1),1)*X   ; y=Connection(:,1:2)-ones(size(Connection,1),1)*Y ;
        z=[intersect(find(x(:,1)==0),find(x(:,2)==0)),intersect(find(y(:,1)==0),find(y(:,2)==0))]             ;
    if  isempty(z)==1
        Connection=[Connection;X]          ;
    end
    end
    if  isempty(Connection)==1
        Connection=[Connection;X]          ;
    end
end
end
end
    CONNECTION{end+1,1}=Connection         ; Connection=[] ;
for I=1:size(EMBRYO,1)
for i=1:size(EMBRYO{I,13},2)
if  isempty(EMBRYO{I,13}{5,i})==0 && isnan(EMBRYO{I,13}{5,i})==0 && EMBRYO{I,13}{5,i}>1
    X=cell2mat(EMBRYO{I,13}(3:4,i))'       ; Y=[X(2),X(1)] ;
    if  isempty(Connection)==0
        x=Connection(:,1:2)-ones(size(Connection,1),1)*X   ; y=Connection(:,1:2)-ones(size(Connection,1),1)*Y ;
        z=[intersect(find(x(:,1)==0),find(x(:,2)==0)),intersect(find(y(:,1)==0),find(y(:,2)==0))]             ;
    if  isempty(z)==1
        Connection=[Connection;[X,1]]      ;
    end
    if  isempty(z)==0
        Connection(z,3)=Connection(z,3)+1  ;
    end
    end
    if  isempty(Connection)==1
        Connection=[Connection;[X,1]]      ;
    end
end
end
end
    CONNECTION{end+1,1}=Connection         ;
    save('WorkSpace_CONNECTION','CONNECTION','-v7.3')          ;
    
  % Contact Duration & Contact Area Surface
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200319 - DataImport\WorkSpace_Query' ) ;
    load('WorkSpace_EMBRYO')   ; load('WorkSpace_MissingCell') ; load('WorkSpace_CONNECTION')                           ;
    CONNECTION{1,2}=cell(size(CONNECTION{1,1},1),1)            ;
for i=1:size(CONNECTION{1,1},1)
    CONNECTION{1,2}{i,1}={}    ; i
for I=1:size(EMBRYO,1)
    contact=importdata(['E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200318 - DataUpdate4\RevisionData\ContactSurface\',EMBRYO{I,1},'_Stat.csv']) ; A=contact.data ;
for j=1:size(EMBRYO{I,13},2)
if  isempty(EMBRYO{I,13}{5,j})==0 && isnan(EMBRYO{I,13}{5,j})==0
if  (CONNECTION{1,1}(i,1)==EMBRYO{I,13}{3,j} && CONNECTION{1,1}(i,2)==EMBRYO{I,13}{4,j}) || (CONNECTION{1,1}(i,2)==EMBRYO{I,13}{3,j} && CONNECTION{1,1}(i,1)==EMBRYO{I,13}{4,j})
    Time=[]                    ;
    for t=EMBRYO{I,11}
    if  isnan(A(t,j))==0 && A(t,j)~=0
        Time=[Time,t]          ;
    end
    end
        t=Time(1)              ;
  while t<=max(EMBRYO{I,11}) && isempty(Time)==0
        t=Time(1)              ;
    if  ismember(t+1,Time)==0
    if  t<10
        ttt=['00',num2str(t)]  ;
    end
    if  t>=10 && t<=99
        ttt=['0' ,num2str(t)]  ;
    end
    if  t>=100
        ttt=[     num2str(t)]  ;
    end
        Data=table2cell(readtable(['E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200318 - DataUpdate4\RevisionData\LostCells\',EMBRYO{I,1},'\',EMBRYO{I,1},'_',ttt,'_nucLoc.csv'])) ;
        surface1=Data{find(cell2mat(Data(:,1))==EMBRYO{I,13}{3,j}),7} ; surface2=Data{find(cell2mat(Data(:,1))==EMBRYO{I,13}{4,j}),7} ;
        CONNECTION{1,2}{i,1}{end+1,1}=[t;A(t,j);surface1;surface2;A(t,j)/surface1;A(t,j)/surface2] ; Time(Time==t)=[]                 ;
    end
    if  isempty(Time)==0 && ismember(t+1,Time)==1
        t0=t                   ;
  while ismember(t0 ,Time)==1
        t0=t0+1                ;
    end
        connection=[]          ;
    for tt=t:1:t0-1
    if  tt<10
        ttt=['00',num2str(tt)] ;
    end
    if  tt>=10 && tt<=99
        ttt=['0' ,num2str(tt)] ;
    end
    if  tt>=100
        ttt=[     num2str(tt)] ;
    end
        Data=table2cell(readtable(['E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200318 - DataUpdate4\RevisionData\LostCells\',EMBRYO{I,1},'\',EMBRYO{I,1},'_',ttt,'_nucLoc.csv'])) ;
        surface1=Data{find(cell2mat(Data(:,1))==EMBRYO{I,13}{3,j}),7} ; surface2=Data{find(cell2mat(Data(:,1))==EMBRYO{I,13}{4,j}),7} ;
        connection=[connection,[t;A(t,j);surface1;surface2;A(t,j)/surface1;A(t,j)/surface2]] ;
    end
        CONNECTION{1,2}{i,1}{end+1,1}=connection ; Time(Time<=t0)=[]  ;
    end
    end
end
end
end
end
end
        save('WorkSpace_CONNECTION','CONNECTION','-v7.3')             ;

  % Contact Detail - Merging
    load('WorkSpace_CONNECTION') ; File=dir('WorkSpace_CONNECTION_*')         ;
    connection=CONNECTION ; connection{1,2}=cell(size(connection{1,1},1),1)   ;
for n=1:length(File)
    load(File(n).name)    ; n
for i=1:size(CONNECTION{1,2},1)
if  isempty(CONNECTION{1,2}{i,1})==0
    connection{1,2}{i,1}=CONNECTION{1,2}{i,1}                         ;
end
end
end
    CONNECTION=connection ; save('WorkSpace_CONNECTION','CONNECTION','-v7.3') ;
    
  % Signaling Cell
    ContactPair={{'P2'},{'EMS'};{'MS'},{'ABal'};{'C'},{'ABar'} ;{'P2'},{'ABp'};{'MS'},{'ABalp'};{'MS'},{'ABara'};{'ABalapa'},{'ABplaaa'};{'ABalapp'},{'ABplaaa'};{'MSapp'},{'ABplpapp'};{'MSappp'},{'ABplpppp'}} ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200319 - DataImport\WorkSpace_Query' ) ; load('WorkSpace_CONNECTION') ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200319 - DataImport\WorkSpace_EMBRYO')                                ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_CellName') ; load('WorkSpace_MissingCell')                         ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200315 - QualityControl_Final\WorkSpace_WTinfo')                      ;
for k1=1:size(ContactPair,1)
for k2=1:2
    ContactPair{k1,k2+4}=0              ;
    for i=1:size(Query,1)
    if  strcmp(Query{i,1}      ,ContactPair{k1,k2})==1
        ContactPair{k1,k2+2}=Query{i,2} ;
    end
    end
    for i=1:size(MissingCell,1)
    if  strcmp(MissingCell{i,1},ContactPair{k1,k2})==1
        ContactPair{k1,k2+4}=1          ;
    end
    end
end
end
for k1=1:size(ContactPair,1)
    ContactPair{k1,8}=[] ; ContactPair{k1,9}=[] ; ContactPair{k1,10}=[] ; ContactPair{k1,11}=[] ; ContactPair{k1,12}=[] ;
if  sum(cell2mat(ContactPair(k1,5:6)))==0
    X1=intersect(find(CONNECTION{4,1}(:,1)==ContactPair{k1,3}),find(CONNECTION{4,1}(:,2)==ContactPair{k1,4}))           ;
    X2=intersect(find(CONNECTION{4,1}(:,2)==ContactPair{k1,3}),find(CONNECTION{4,1}(:,1)==ContactPair{k1,4}))           ;
    X=[X1,X2] ;
if  isempty(X)==0
    ContactPair{k1,7}=CONNECTION{4,1}(X,3)      ;
    Y1=intersect(find(CONNECTION{1,1}(:,1)==CONNECTION{4,1}(X,1)),find(CONNECTION{1,1}(:,2)==CONNECTION{4,1}(X,2)))     ;
    Y2=intersect(find(CONNECTION{1,1}(:,2)==CONNECTION{4,1}(X,1)),find(CONNECTION{1,1}(:,1)==CONNECTION{4,1}(X,2)))     ;
    Y=[Y1,Y2] ; T=[] ; S1=[] ; S2=[]            ; s1=[] ; s2=[] ; C=[] ; Cmax=[]                      ;
for j =1:length(CONNECTION{1,2}{Y,1})
    T=[T,size(CONNECTION{1,2}{Y,1}{j,1},2)*WTinfo{end,2}]       ; C=[C,CONNECTION{1,2}{Y,1}{j,1}(2,:)*(0.25^2)]         ;
    Cmax=[Cmax,max(CONNECTION{1,2}{Y,1}{j,1}(2,:)*(0.25^2))]    ;
if  CONNECTION{1,1}(Y,1)==ContactPair{k1,3} && CONNECTION{1,1}(Y,2)==ContactPair{k1,4}
    S1=[S1,CONNECTION{1,2}{Y,1}{j,1}(3,:)*(0.25^2)] ; S2=[S2,CONNECTION{1,2}{Y,1}{j,1}(4,:)*(0.25^2)] ;
    s1=[s1,CONNECTION{1,2}{Y,1}{j,1}(5,:)*100]      ; s2=[s2,CONNECTION{1,2}{Y,1}{j,1}(6,:)*100]      ;
end
if  CONNECTION{1,1}(Y,2)==ContactPair{k1,3} && CONNECTION{1,1}(Y,1)==ContactPair{k1,4}
    S1=[S1,CONNECTION{1,2}{Y,1}{j,1}(4,:)*(0.25^2)] ; S2=[S2,CONNECTION{1,2}{Y,1}{j,1}(3,:)*(0.25^2)] ;
    s1=[s1,CONNECTION{1,2}{Y,1}{j,1}(6,:)*100]      ; s2=[s2,CONNECTION{1,2}{Y,1}{j,1}(5,:)*100]      ;
end
end
    ContactPair{k1,8} =[sprintf('%2.2f',sum(T)/size(EMBRYO,1))                      ] ;
    ContactPair{k1,9} =[sprintf('%2.2f',mean(S1)) ,' ¡À ',sprintf('%2.2f',std(S1))  ] ;
    ContactPair{k1,10}=[sprintf('%2.2f',mean(S2)) ,' ¡À ',sprintf('%2.2f',std(S2))  ] ;
    ContactPair{k1,11}=[sprintf('%2.2f',mean(s1)) ,' ¡À ',sprintf('%2.2f',std(s1))  ] ;
    ContactPair{k1,12}=[sprintf('%2.2f',mean(s2)) ,' ¡À ',sprintf('%2.2f',std(s2))  ] ;
    ContactPair{k1,13}=[sprintf('%2.2f',mean(C))  ,' ¡À ',sprintf('%2.2f',std(C))   ] ;
    ContactPair{k1,14}=[sprintf('%2.2f',max(Cmax)),' ¡À ',sprintf('%2.2f',std(Cmax))] ;
end
end
end

  % Cell Pair of MSappp ¡ú ABplpppp (10/11)
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200319 - DataImport\WorkSpace_Query') ;
    CellA='MSappp'     ; CellB='ABplpppp' ; ValidEmbryo=[] ; load('WorkSpace_EMBRYO') ;
for i=1:size(Query,1)
if  strcmp(Query{i,1},CellA)==1
    NumberA=Query{i,2} ;
end
if  strcmp(Query{i,1},CellB)==1
    NumberB=Query{i,2} ;
end
end 
for I=1:size(EMBRYO,1)
for i=2:size(EMBRYO{I,13},2)
if  (EMBRYO{I,13}{3,i}==NumberA && EMBRYO{I,13}{4,i}==NumberB) || (EMBRYO{I,13}{3,i}==NumberB && EMBRYO{I,13}{4,i}==NumberA)
if  isempty(EMBRYO{I,13}{5,i})==0 && isnan(EMBRYO{I,13}{5,i})==0 && EMBRYO{I,13}{5,i}>1
    ValidEmbryo=[ValidEmbryo,I] ;
end
end
end
end

  % Cell Pair of C ¡ú ABar (3/14)
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200319 - DataImport\WorkSpace_Query') ;
    CellA='C'          ; CellB='ABar' ; ValidEmbryo=[] ; load('WorkSpace_EMBRYO') ; load('WorkSpace_CONNECTION')       ;
for i=1:size(Query,1)
if  strcmp(Query{i,1},CellA)==1
    NumberA=Query{i,2} ;
end
if  strcmp(Query{i,1},CellB)==1
    NumberB=Query{i,2} ;
end
end
for I=1:size(EMBRYO,1)
for i=2:size(EMBRYO{I,13},2)
if  (EMBRYO{I,13}{3,i}==NumberA && EMBRYO{I,13}{4,i}==NumberB) || (EMBRYO{I,13}{3,i}==NumberB && EMBRYO{I,13}{4,i}==NumberA)
if  isempty(EMBRYO{I,13}{5,i})==0 && isnan(EMBRYO{I,13}{5,i})==0 && EMBRYO{I,13}{5,i}>1
    ValidEmbryo=[ValidEmbryo,I] ;
end
end
end
end

  % Cell Pair of MSapp ¡ú ABplpapp (12/7)
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200323 - MissingCell\WorkSpace_MissingCell')         ;
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200319 - DataImport\WorkSpace_EMBRYO') ; Invalid1=[] ; Invalid2=[] ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_CellName') ; CellA='MSapp' ; CellB='ABplpapp'     ;
for k0=1:length(CellName)
for k1=1:size(CellName{k0,1},1)
for k2=1:size(CellName{k0,1},2)
if  isnumeric(CellName{k0,1}{k1,k2})==0 && strcmp(CellName{k0,1}{k1,k2},CellA)==1
    for I=1:size(EMBRYO,1)
    if  MissingCell{I,1}{k0,1}{k1,k2}==1
        Invalid1=[Invalid1,I] ;
    end
    end
end
if  isnumeric(CellName{k0,1}{k1,k2})==0 && strcmp(CellName{k0,1}{k1,k2},CellB)==1
    for I=1:size(EMBRYO,1)
    if  MissingCell{I,1}{k0,1}{k1,k2}==1
        Invalid2=[Invalid2,I] ;
    end
    end
end
end
end
end

  % Data Generation
    load('WorkSpace_CONNECTION') ; load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200424 - QueryUpdate\WorkSpace_QUERY') ;
    TABLE1=cell(1,2) ; TABLE2=cell(1,2) ; TABLE3=cell(1,2) ; TABLE4=cell(1,2) ; 
for i=1:size(CONNECTION{1,1},1)
    i
for j=1:size(Query,1)
if  CONNECTION{1,1}(i,1)==Query{j,2}
    TABLE1{end,1}=Query{j,1}     ;
end
if  CONNECTION{1,1}(i,2)==Query{j,2}
    TABLE1{end,2}=Query{j,1}     ;
end
end
    TABLE1{end+1,1}={}           ;
end
for i=1:size(CONNECTION{2,1},1)
    i
for j=1:size(Query,1)
if  CONNECTION{1,1}(i,1)==Query{j,2}
    TABLE2{end,1}=Query{j,1}     ;
end
if  CONNECTION{1,1}(i,2)==Query{j,2}
    TABLE2{end,2}=Query{j,1}     ;
end
end
    TABLE2{end+1,1}={}           ;
end
for i=1:size(CONNECTION{3,1},1)
    i
for j=1:size(Query,1)
if  CONNECTION{1,1}(i,1)==Query{j,2}
    TABLE3{end,1}=Query{j,1}     ;
end
if  CONNECTION{1,1}(i,2)==Query{j,2}
    TABLE3{end,2}=Query{j,1}     ;
end
end
    TABLE3{end+1,1}={}           ;
end
for i=1:size(CONNECTION{4,1},1)
if  CONNECTION{4,1}(i,3)==17
    i
for j=1:size(Query,1)
if  CONNECTION{1,1}(i,1)==Query{j,2}
    TABLE4{end,1}=Query{j,1}     ;
end
if  CONNECTION{1,1}(i,2)==Query{j,2}
    TABLE4{end,2}=Query{j,1}     ;
end
end
    TABLE4{end+1,1}={}           ;
end    
end