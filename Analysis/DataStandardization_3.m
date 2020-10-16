  % Data Standardization

  
  
  % Mission 1 : Normalization

    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_CellName')     ; load('WorkSpace_EMBRYO')         ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 2 - WT Position\WorkSpace_Sequence0') ; SSS=Sequence0 ; SSS(:,58:end)=[] ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_WTinfo')       ; load('WorkSpace_HLH1')           ;
    load('Normalization1\WorkSpace_Criteria') ;
for n =1:size(Sequence0,2)
    load(['Normalization1\WorkSpace_Sequence000_',num2str(n)])            ; Sequence0{14,n}=[]     ;
    Sequence=Sequence0          ; Sequence(18:20,n)=Sequence(12:14,n)     ; n
for nn=1:size(Sequence,2)
    Sequence{12,nn}=[]          ; Sequence{13,nn}=[] ; Sequence{14,nn}=[] ; Sequence{20,nn}=[]     ;
end
for N=1:30
    N
  % Calculation for Average Position
    PositionAve=CellName            ;
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
        PositionAve{k0,1}{k1,k2}=[] ;
    end
    end
    end
for i=1:length(Sequence{18,n})
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,n}{1,i}{k0,1}{k1,k2})==3
        PositionAve{k0,1}{k1,k2}=[PositionAve{k0,1}{k1,k2};Sequence{18,n}{1,i}{k0,1}{k1,k2}] ;
    end
    end
    end
    end
end
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(PositionAve{k0,1}{k1,k2}) <3
        PositionAve{k0,1}{k1,k2}=NaN ;
    end
    if  length(PositionAve{k0,1}{k1,k2})>=3
        PositionAve{k0,1}{k1,k2}=mean(PositionAve{k0,1}{k1,k2}) ;
    end
    end
    end
    end
    
  % Optimization
for i=1:length(Sequence{18,n})
    Angle=pi/180*[-15:1:15] ; Distance=[-10:0.25:10]            ;
    
  % X Rotation
    Var=[] ;
for k=Angle
    var=[] ;
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,n}{1,i}{k0,1}{k1,k2})==3
        y= Sequence{18,n}{1,i}{k0,1}{k1,k2}(1,3)*cos(k)+Sequence{18,n}{1,i}{k0,1}{k1,k2}(1,1)*sin(k) ;   
        z=-Sequence{18,n}{1,i}{k0,1}{k1,k2}(1,3)*sin(k)+Sequence{18,n}{1,i}{k0,1}{k1,k2}(1,1)*cos(k) ;
        var=[var;Criteria{1,n}{k0,1}{k1,k2}*sum(([z,Sequence{18,n}{1,i}{k0,1}{k1,k2}(1,2),y]-PositionAve{k0,1}{k1,k2}).^2)] ;
    end
    end
    end
    end
        Var=[Var,sum(var)] ;
end
        k=Angle(1,find(Var(1,:)==min(Var))) ; k=mean(k) ;
       
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,n}{1,i}{k0,1}{k1,k2})==3
        y= Sequence{18,n}{1,i}{k0,1}{k1,k2}(1,3)*cos(k)+Sequence{18,n}{1,i}{k0,1}{k1,k2}(1,1)*sin(k) ;   
        z=-Sequence{18,n}{1,i}{k0,1}{k1,k2}(1,3)*sin(k)+Sequence{18,n}{1,i}{k0,1}{k1,k2}(1,1)*cos(k) ;
        Sequence{18,n}{1,i}{k0,1}{k1,k2}(1,3)=y ; Sequence{18,n}{1,i}{k0,1}{k1,k2}(1,1)=z            ;
    end
    end
    end
    end
    if  ismember(i,[1:1:size(EMBRYO,1)])==1
    for j=1:size(Sequence{19,n}{1,i},1)
        y= Sequence{19,n}{1,i}{j,1}(:,3)*cos(k)+Sequence{19,n}{1,i}{j,1}(:,1)*sin(k) ;
        z=-Sequence{19,n}{1,i}{j,1}(:,3)*sin(k)+Sequence{19,n}{1,i}{j,1}(:,1)*cos(k) ;
        Sequence{19,n}{1,i}{j,1}(:,3)=y ; Sequence{19,n}{1,i}{j,1}(:,1)=z            ;
    end
    end

  % Y Rotation
    Var=[] ;
for k=Angle
    var=[] ;
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,n}{1,i}{k0,1}{k1,k2})==3
        x= Sequence{18,n}{1,i}{k0,1}{k1,k2}(1,2)*cos(k)+Sequence{18,n}{1,i}{k0,1}{k1,k2}(1,1)*sin(k) ;   
        z=-Sequence{18,n}{1,i}{k0,1}{k1,k2}(1,2)*sin(k)+Sequence{18,n}{1,i}{k0,1}{k1,k2}(1,1)*cos(k) ;
        var=[var;Criteria{1,n}{k0,1}{k1,k2}*sum(([z,x,Sequence{18,n}{1,i}{k0,1}{k1,k2}(1,3)]-PositionAve{k0,1}{k1,k2}).^2)] ;
    end
    end
    end
    end
        Var=[Var,sum(var)] ;
end
        k=Angle(1,find(Var(1,:)==min(Var))) ; k=mean(k) ;
                
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,n}{1,i}{k0,1}{k1,k2})==3
        x= Sequence{18,n}{1,i}{k0,1}{k1,k2}(1,2)*cos(k)+Sequence{18,n}{1,i}{k0,1}{k1,k2}(1,1)*sin(k) ;   
        z=-Sequence{18,n}{1,i}{k0,1}{k1,k2}(1,2)*sin(k)+Sequence{18,n}{1,i}{k0,1}{k1,k2}(1,1)*cos(k) ;
        Sequence{18,n}{1,i}{k0,1}{k1,k2}(1,2)=x ; Sequence{18,n}{1,i}{k0,1}{k1,k2}(1,1)=z            ;
    end
    end
    end
    end
    if  ismember(i,[1:1:size(EMBRYO,1)])==1
    for j=1:size(Sequence{19,n}{1,i},1)
        x= Sequence{19,n}{1,i}{j,1}(:,2)*cos(k)+Sequence{19,n}{1,i}{j,1}(:,1)*sin(k) ;
        z=-Sequence{19,n}{1,i}{j,1}(:,2)*sin(k)+Sequence{19,n}{1,i}{j,1}(:,1)*cos(k) ;
        Sequence{19,n}{1,i}{j,1}(:,2)=x ; Sequence{19,n}{1,i}{j,1}(:,1)=z            ;
    end
    end
        
  % Z Rotation
    Var=[] ;
for k=Angle
    var=[] ;
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,n}{1,i}{k0,1}{k1,k2})==3
        x= Sequence{18,n}{1,i}{k0,1}{k1,k2}(1,2)*cos(k)+Sequence{18,n}{1,i}{k0,1}{k1,k2}(1,3)*sin(k) ;   
        y=-Sequence{18,n}{1,i}{k0,1}{k1,k2}(1,2)*sin(k)+Sequence{18,n}{1,i}{k0,1}{k1,k2}(1,3)*cos(k) ;
        var=[var;Criteria{1,n}{k0,1}{k1,k2}*sum(([Sequence{18,n}{1,i}{k0,1}{k1,k2}(1,1),x,y]-PositionAve{k0,1}{k1,k2}).^2)] ;
    end
    end
    end
    end
        Var=[Var,sum(var)] ;
end
        k=Angle(1,find(Var(1,:)==min(Var))) ; k=mean(k) ;
                
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,n}{1,i}{k0,1}{k1,k2})==3
        x= Sequence{18,n}{1,i}{k0,1}{k1,k2}(1,2)*cos(k)+Sequence{18,n}{1,i}{k0,1}{k1,k2}(1,3)*sin(k) ;   
        y=-Sequence{18,n}{1,i}{k0,1}{k1,k2}(1,2)*sin(k)+Sequence{18,n}{1,i}{k0,1}{k1,k2}(1,3)*cos(k) ;
        Sequence{18,n}{1,i}{k0,1}{k1,k2}(1,2)=x ; Sequence{18,n}{1,i}{k0,1}{k1,k2}(1,3)=y            ;
    end
    end
    end
    end
    if  ismember(i,[1:1:size(EMBRYO,1)])==1
    for j=1:size(Sequence{19,n}{1,i},1)
        x= Sequence{19,n}{1,i}{j,1}(:,2)*cos(k)+Sequence{19,n}{1,i}{j,1}(:,3)*sin(k) ;
        y=-Sequence{19,n}{1,i}{j,1}(:,2)*sin(k)+Sequence{19,n}{1,i}{j,1}(:,3)*cos(k) ;
        Sequence{19,n}{1,i}{j,1}(:,2)=x ; Sequence{19,n}{1,i}{j,1}(:,3)=y            ;
    end
    end
        
  % X Translation
    Var=[] ;
for k=Distance
    var=[] ;
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,n}{1,i}{k0,1}{k1,k2})==3
        var=[var;Criteria{1,n}{k0,1}{k1,k2}*sum((Sequence{18,n}{1,i}{k0,1}{k1,k2}+[0,k,0]-PositionAve{k0,1}{k1,k2}).^2)] ;
    end
    end
    end
    end
        Var=[Var,sum(var)] ;
end
        k=Distance(1,find(Var(1,:)==min(Var))) ; k=mean(k) ;
            
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,n}{1,i}{k0,1}{k1,k2})==3
        Sequence{18,n}{1,i}{k0,1}{k1,k2}=Sequence{18,n}{1,i}{k0,1}{k1,k2}+[0,k,0] ;
    end
    end
    end
    end
    if  ismember(i,[1:1:size(EMBRYO,1)])==1
    for j=1:size(Sequence{19,n}{1,i},1)
        Sequence{19,n}{1,i}{j,1}(:,2)=Sequence{19,n}{1,i}{j,1}(:,2)+k             ;
    end 
    end
    
  % Y Translation
    Var=[] ;
for k=Distance
    var=[] ;
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,n}{1,i}{k0,1}{k1,k2})==3
        var=[var;Criteria{1,n}{k0,1}{k1,k2}*sum((Sequence{18,n}{1,i}{k0,1}{k1,k2}+[0,0,k]-PositionAve{k0,1}{k1,k2}).^2)] ;
    end
    end
    end
    end
        Var=[Var,sum(var)] ;
end
        k=Distance(1,find(Var(1,:)==min(Var))) ; k=mean(k) ;
                
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,n}{1,i}{k0,1}{k1,k2})==3
        Sequence{18,n}{1,i}{k0,1}{k1,k2}=Sequence{18,n}{1,i}{k0,1}{k1,k2}+[0,0,k] ;
    end
    end
    end
    end
    if  ismember(i,[1:1:size(EMBRYO,1)])==1
    for j=1:size(Sequence{19,n}{1,i},1)
        Sequence{19,n}{1,i}{j,1}(:,3)=Sequence{19,n}{1,i}{j,1}(:,3)+k             ;
    end
    end
    
  % Z Translation
    Var=[] ;
for k=Distance
    var=[] ;
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,n}{1,i}{k0,1}{k1,k2})==3
        var=[var;Criteria{1,n}{k0,1}{k1,k2}*sum((Sequence{18,n}{1,i}{k0,1}{k1,k2}+[k,0,0]-PositionAve{k0,1}{k1,k2}).^2)] ;
    end
    end
    end
    end
        Var=[Var,sum(var)] ;
end
        k=Distance(1,find(Var(1,:)==min(Var))) ; k=mean(k) ;
                
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,n}{1,i}{k0,1}{k1,k2})==3
        Sequence{18,n}{1,i}{k0,1}{k1,k2}=Sequence{18,n}{1,i}{k0,1}{k1,k2}+[k,0,0] ;
    end
    end
    end
    end
    if  ismember(i,[1:1:size(EMBRYO,1)])==1
    for j=1:size(Sequence{19,n}{1,i},1)
        Sequence{19,n}{1,i}{j,1}(:,1)=Sequence{19,n}{1,i}{j,1}(:,1)+k             ;
    end
    end
end
end
        save(['Normalization2\WorkSpace_00',num2str(n)],'Sequence','-v7.3')       ;
end