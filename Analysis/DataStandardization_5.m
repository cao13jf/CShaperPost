  % Data Standardization

  
  
  % Mission 1 : Normalization
   
for n=1:57
    load(['Normalization4\WorkSpace_Eggshell_',num2str(n)]) ; load(['Normalization4\WorkSpace_Criteria',num2str(n)])    ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_CellName') ; load('WorkSpace_HLH1') ;
    load('WorkSpace_EMBRYO') ; Angle=pi/180*[-15:1:15] ; Distance=[-10:0.25:10] ; Scale=[0.75:0.025:1.25]               ; n

for I=1:length(HLH1)
    Sequence=cell(18,1) ; Sequence{18,1}=Sequence0{18,1}{1,I+size(EMBRYO,1)}    ; I

for N=1:30
    
  % X Rotation
    Var=[] ;
for k=Angle
    var=[] ;
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,1}{k0,1}{k1,k2})==3 && Criteria{1,n}{k0,1}{k1,k2}>0 && length(Criteria{2,n}{k0,1}{k1,k2})==3
        y= Sequence{18,1}{k0,1}{k1,k2}(1,3)*cos(k)+Sequence{18,1}{k0,1}{k1,k2}(1,1)*sin(k) ;
        z=-Sequence{18,1}{k0,1}{k1,k2}(1,3)*sin(k)+Sequence{18,1}{k0,1}{k1,k2}(1,1)*cos(k) ;
        var=[var;Criteria{1,n}{k0,1}{k1,k2}*sum(([z,Sequence{18,1}{k0,1}{k1,k2}(1,2),y]-Criteria{2,n}{k0,1}{k1,k2}).^2)] ;
    end
    end
    end
    end
        Var=[Var,sum(var)] ;
end
        k=Angle(1,find(Var(1,:)==min(Var))) ; k=k(1) ;
                
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,1}{k0,1}{k1,k2})==3
        y= Sequence{18,1}{k0,1}{k1,k2}(1,3)*cos(k)+Sequence{18,1}{k0,1}{k1,k2}(1,1)*sin(k) ;
        z=-Sequence{18,1}{k0,1}{k1,k2}(1,3)*sin(k)+Sequence{18,1}{k0,1}{k1,k2}(1,1)*cos(k) ;
        Sequence{18,1}{k0,1}{k1,k2}(1,3)=y  ; Sequence{18,1}{k0,1}{k1,k2}(1,1)=z           ;
    end
    end
    end
    end

  % Y Rotation
    Var=[] ;
for k=Angle
    var=[] ;
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,1}{k0,1}{k1,k2})==3 && Criteria{1,n}{k0,1}{k1,k2}>0 && length(Criteria{2,n}{k0,1}{k1,k2})==3
        x= Sequence{18,1}{k0,1}{k1,k2}(1,2)*cos(k)+Sequence{18,1}{k0,1}{k1,k2}(1,1)*sin(k) ;
        z=-Sequence{18,1}{k0,1}{k1,k2}(1,2)*sin(k)+Sequence{18,1}{k0,1}{k1,k2}(1,1)*cos(k) ;
        var=[var;Criteria{1,n}{k0,1}{k1,k2}*sum(([z,x,Sequence{18,1}{k0,1}{k1,k2}(1,3)]-Criteria{2,n}{k0,1}{k1,k2}).^2)] ;
    end
    end
    end
    end
        Var=[Var,sum(var)] ;
end
        k=Angle(1,find(Var(1,:)==min(Var))) ; k=k(1) ;
        
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,1}{k0,1}{k1,k2})==3
        x= Sequence{18,1}{k0,1}{k1,k2}(1,2)*cos(k)+Sequence{18,1}{k0,1}{k1,k2}(1,1)*sin(k) ;
        z=-Sequence{18,1}{k0,1}{k1,k2}(1,2)*sin(k)+Sequence{18,1}{k0,1}{k1,k2}(1,1)*cos(k) ;
        Sequence{18,1}{k0,1}{k1,k2}(1,2)=x  ; Sequence{18,1}{k0,1}{k1,k2}(1,1)=z           ;
    end
    end
    end
    end

  % Z Rotation
    Var=[] ;
for k=Angle
    var=[] ;
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,1}{k0,1}{k1,k2})==3 && Criteria{1,n}{k0,1}{k1,k2}>0 && length(Criteria{2,n}{k0,1}{k1,k2})==3
        x= Sequence{18,1}{k0,1}{k1,k2}(1,2)*cos(k)+Sequence{18,1}{k0,1}{k1,k2}(1,3)*sin(k) ;
        y=-Sequence{18,1}{k0,1}{k1,k2}(1,2)*sin(k)+Sequence{18,1}{k0,1}{k1,k2}(1,3)*cos(k) ;
        var=[var;Criteria{1,n}{k0,1}{k1,k2}*sum(([Sequence{18,1}{k0,1}{k1,k2}(1,1),x,y]-Criteria{2,n}{k0,1}{k1,k2}).^2)] ;
    end
    end
    end
    end
        Var=[Var,sum(var)] ;
end
        k=Angle(1,find(Var(1,:)==min(Var))) ; k=k(1) ;

    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,1}{k0,1}{k1,k2})==3
        x= Sequence{18,1}{k0,1}{k1,k2}(1,2)*cos(k)+Sequence{18,1}{k0,1}{k1,k2}(1,3)*sin(k) ;   
        y=-Sequence{18,1}{k0,1}{k1,k2}(1,2)*sin(k)+Sequence{18,1}{k0,1}{k1,k2}(1,3)*cos(k) ;
        Sequence{18,1}{k0,1}{k1,k2}(1,2)=x  ; Sequence{18,1}{k0,1}{k1,k2}(1,3)=y           ;
    end
    end
    end
    end
        
  % X Translation
    Var=[] ;
for k=Distance
    var=[] ;
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,1}{k0,1}{k1,k2})==3 && Criteria{1,n}{k0,1}{k1,k2}>0 && length(Criteria{2,n}{k0,1}{k1,k2})==3
        var=[var;Criteria{1,n}{k0,1}{k1,k2}*sum((Sequence{18,1}{k0,1}{k1,k2}+[0,k,0]-Criteria{2,n}{k0,1}{k1,k2}).^2)] ;
    end
    end
    end
    end
        Var=[Var,sum(var)] ;
end
        k=Distance(1,find(Var(1,:)==min(Var))) ; k=k(1) ;
            
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,1}{k0,1}{k1,k2})==3
        Sequence{18,1}{k0,1}{k1,k2}=Sequence{18,1}{k0,1}{k1,k2}+[0,k,0] ;
    end
    end
    end
    end
    
  % Y Translation
    Var=[] ;
for k=Distance
    var=[] ;
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,1}{k0,1}{k1,k2})==3 && Criteria{1,n}{k0,1}{k1,k2}>0 && length(Criteria{2,n}{k0,1}{k1,k2})==3
        var=[var;Criteria{1,n}{k0,1}{k1,k2}*sum((Sequence{18,1}{k0,1}{k1,k2}+[0,0,k]-Criteria{2,n}{k0,1}{k1,k2}).^2)] ;
    end
    end
    end
    end
        Var=[Var,sum(var)] ;
end
        k=Distance(1,find(Var(1,:)==min(Var))) ; k=k(1) ;
                
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,1}{k0,1}{k1,k2})==3
        Sequence{18,1}{k0,1}{k1,k2}=Sequence{18,1}{k0,1}{k1,k2}+[0,0,k] ;
    end
    end
    end
    end
    
  % Z Translation
    Var=[] ;
for k=Distance
    var=[] ;
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,1}{k0,1}{k1,k2})==3 && Criteria{1,n}{k0,1}{k1,k2}>0 && length(Criteria{2,n}{k0,1}{k1,k2})==3
        var=[var;Criteria{1,n}{k0,1}{k1,k2}*sum((Sequence{18,1}{k0,1}{k1,k2}+[k,0,0]-Criteria{2,n}{k0,1}{k1,k2}).^2)] ;
    end
    end
    end
    end
        Var=[Var,sum(var)] ;
end
        k=Distance(1,find(Var(1,:)==min(Var))) ; k=k(1) ;
                
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,1}{k0,1}{k1,k2})==3
        Sequence{18,1}{k0,1}{k1,k2}=Sequence{18,1}{k0,1}{k1,k2}+[k,0,0] ;
    end
    end
    end
    end
end

for N=1:30
    
  % X Scaling
    Var=[] ;
for k=Scale
    var=[] ;
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,1}{k0,1}{k1,k2})==3 && Criteria{1,n}{k0,1}{k1,k2}>0 && length(Criteria{2,n}{k0,1}{k1,k2})==3
        var=[var;Criteria{1,n}{k0,1}{k1,k2}*sum((Sequence{18,1}{k0,1}{k1,k2}.*[1,k,1]-Criteria{2,n}{k0,1}{k1,k2}).^2)] ;
    end
    end
    end
    end
        Var=[Var,sum(var)] ;
end
        k=Scale(1,find(Var(1,:)==min(Var))) ; k=k(1) ;
                
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,1}{k0,1}{k1,k2})==3
        Sequence{18,1}{k0,1}{k1,k2}=Sequence{18,1}{k0,1}{k1,k2}.*[1,k,1] ;
    end
    end
    end
    end

  % Y Scaling
    Var=[] ;
for k=Scale
    var=[] ;
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,1}{k0,1}{k1,k2})==3 && Criteria{1,n}{k0,1}{k1,k2}>0 && length(Criteria{2,n}{k0,1}{k1,k2})==3
        var=[var;Criteria{1,n}{k0,1}{k1,k2}*sum((Sequence{18,1}{k0,1}{k1,k2}.*[1,1,k]-Criteria{2,n}{k0,1}{k1,k2}).^2)] ;
    end
    end
    end
    end
        Var=[Var,sum(var)] ;
end
        k=Scale(1,find(Var(1,:)==min(Var))) ; k=k(1) ;
                
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,1}{k0,1}{k1,k2})==3
        Sequence{18,1}{k0,1}{k1,k2}=Sequence{18,1}{k0,1}{k1,k2}.*[1,1,k] ;
    end
    end
    end
    end
    
  % Z Scaling
    Var=[] ;
for k=Scale
    var=[] ;
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,1}{k0,1}{k1,k2})==3 && Criteria{1,n}{k0,1}{k1,k2}>0 && length(Criteria{2,n}{k0,1}{k1,k2})==3
        var=[var;Criteria{1,n}{k0,1}{k1,k2}*sum((Sequence{18,1}{k0,1}{k1,k2}.*[k,1,1]-Criteria{2,n}{k0,1}{k1,k2}).^2)] ;
    end
    end
    end
    end
        Var=[Var,sum(var)] ;
end
        k=Scale(1,find(Var(1,:)==min(Var))) ; k=k(1) ;
                
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,1}{k0,1}{k1,k2})==3
        Sequence{18,1}{k0,1}{k1,k2}=Sequence{18,1}{k0,1}{k1,k2}.*[k,1,1] ;
    end
    end
    end
    end
end

for N=1:30
    
  % X Rotation
    Var=[] ;
for k=Angle
    var=[] ;
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,1}{k0,1}{k1,k2})==3 && Criteria{1,n}{k0,1}{k1,k2}>0 && length(Criteria{2,n}{k0,1}{k1,k2})==3
        y= Sequence{18,1}{k0,1}{k1,k2}(1,3)*cos(k)+Sequence{18,1}{k0,1}{k1,k2}(1,1)*sin(k) ;   
        z=-Sequence{18,1}{k0,1}{k1,k2}(1,3)*sin(k)+Sequence{18,1}{k0,1}{k1,k2}(1,1)*cos(k) ;
        var=[var;Criteria{1,n}{k0,1}{k1,k2}*sum(([z,Sequence{18,1}{k0,1}{k1,k2}(1,2),y]-Criteria{2,n}{k0,1}{k1,k2}).^2)] ;
    end
    end
    end
    end
        Var=[Var,sum(var)] ;
end
        k=Angle(1,find(Var(1,:)==min(Var))) ; k=k(1) ;
                
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,1}{k0,1}{k1,k2})==3
        y= Sequence{18,1}{k0,1}{k1,k2}(1,3)*cos(k)+Sequence{18,1}{k0,1}{k1,k2}(1,1)*sin(k) ;   
        z=-Sequence{18,1}{k0,1}{k1,k2}(1,3)*sin(k)+Sequence{18,1}{k0,1}{k1,k2}(1,1)*cos(k) ;
        Sequence{18,1}{k0,1}{k1,k2}(1,3)=y  ; Sequence{18,1}{k0,1}{k1,k2}(1,1)=z           ;
    end
    end
    end
    end
    
  % Y Rotation
    Var=[] ;
for k=Angle
    var=[] ;
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,1}{k0,1}{k1,k2})==3 && Criteria{1,n}{k0,1}{k1,k2}>0 && length(Criteria{2,n}{k0,1}{k1,k2})==3
        x= Sequence{18,1}{k0,1}{k1,k2}(1,2)*cos(k)+Sequence{18,1}{k0,1}{k1,k2}(1,1)*sin(k) ;   
        z=-Sequence{18,1}{k0,1}{k1,k2}(1,2)*sin(k)+Sequence{18,1}{k0,1}{k1,k2}(1,1)*cos(k) ;
        var=[var;Criteria{1,n}{k0,1}{k1,k2}*sum(([z,x,Sequence{18,1}{k0,1}{k1,k2}(1,3)]-Criteria{2,n}{k0,1}{k1,k2}).^2)] ;
    end
    end
    end
    end
        Var=[Var,sum(var)] ;
end
        k=Angle(1,find(Var(1,:)==min(Var))) ; k=k(1) ;
                
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,1}{k0,1}{k1,k2})==3
        x= Sequence{18,1}{k0,1}{k1,k2}(1,2)*cos(k)+Sequence{18,1}{k0,1}{k1,k2}(1,1)*sin(k) ;   
        z=-Sequence{18,1}{k0,1}{k1,k2}(1,2)*sin(k)+Sequence{18,1}{k0,1}{k1,k2}(1,1)*cos(k) ;
        Sequence{18,1}{k0,1}{k1,k2}(1,2)=x  ; Sequence{18,1}{k0,1}{k1,k2}(1,1)=z           ;
    end
    end
    end
    end

  % Z Rotation
    Var=[] ;
for k=Angle
    var=[] ;
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,1}{k0,1}{k1,k2})==3 && Criteria{1,n}{k0,1}{k1,k2}>0 && length(Criteria{2,n}{k0,1}{k1,k2})==3
        x= Sequence{18,1}{k0,1}{k1,k2}(1,2)*cos(k)+Sequence{18,1}{k0,1}{k1,k2}(1,3)*sin(k) ;   
        y=-Sequence{18,1}{k0,1}{k1,k2}(1,2)*sin(k)+Sequence{18,1}{k0,1}{k1,k2}(1,3)*cos(k) ;
        var=[var;Criteria{1,n}{k0,1}{k1,k2}*sum(([Sequence{18,1}{k0,1}{k1,k2}(1,1),x,y]-Criteria{2,n}{k0,1}{k1,k2}).^2)] ;
    end
    end
    end
    end
        Var=[Var,sum(var)] ;
end
        k=Angle(1,find(Var(1,:)==min(Var))) ; k=k(1) ;
                
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,1}{k0,1}{k1,k2})==3
        x= Sequence{18,1}{k0,1}{k1,k2}(1,2)*cos(k)+Sequence{18,1}{k0,1}{k1,k2}(1,3)*sin(k) ;   
        y=-Sequence{18,1}{k0,1}{k1,k2}(1,2)*sin(k)+Sequence{18,1}{k0,1}{k1,k2}(1,3)*cos(k) ;
        Sequence{18,1}{k0,1}{k1,k2}(1,2)=x  ; Sequence{18,1}{k0,1}{k1,k2}(1,3)=y           ;
    end
    end
    end
    end
    
  % X Translation
    Var=[] ;
for k=Distance
    var=[] ;
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,1}{k0,1}{k1,k2})==3 && Criteria{1,n}{k0,1}{k1,k2}>0 && length(Criteria{2,n}{k0,1}{k1,k2})==3
        var=[var;Criteria{1,n}{k0,1}{k1,k2}*sum((Sequence{18,1}{k0,1}{k1,k2}+[0,k,0]-Criteria{2,n}{k0,1}{k1,k2}).^2)] ;
    end
    end
    end
    end
        Var=[Var,sum(var)] ;
end
        k=Distance(1,find(Var(1,:)==min(Var))) ; k=k(1) ;
            
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,1}{k0,1}{k1,k2})==3
        Sequence{18,1}{k0,1}{k1,k2}=Sequence{18,1}{k0,1}{k1,k2}+[0,k,0] ;
    end
    end
    end
    end
    
  % Y Translation
    Var=[] ;
for k=Distance
    var=[] ;
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,1}{k0,1}{k1,k2})==3 && Criteria{1,n}{k0,1}{k1,k2}>0 && length(Criteria{2,n}{k0,1}{k1,k2})==3
        var=[var;Criteria{1,n}{k0,1}{k1,k2}*sum((Sequence{18,1}{k0,1}{k1,k2}+[0,0,k]-Criteria{2,n}{k0,1}{k1,k2}).^2)] ;
    end
    end
    end
    end
        Var=[Var,sum(var)] ;
end
        k=Distance(1,find(Var(1,:)==min(Var))) ; k=k(1) ;
                
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,1}{k0,1}{k1,k2})==3
        Sequence{18,1}{k0,1}{k1,k2}=Sequence{18,1}{k0,1}{k1,k2}+[0,0,k] ;
    end
    end
    end
    end
    
  % Z Translation
    Var=[] ;
for k=Distance
    var=[] ;
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,1}{k0,1}{k1,k2})==3 && Criteria{1,n}{k0,1}{k1,k2}>0 && length(Criteria{2,n}{k0,1}{k1,k2})==3
        var=[var;Criteria{1,n}{k0,1}{k1,k2}*sum((Sequence{18,1}{k0,1}{k1,k2}+[k,0,0]-Criteria{2,n}{k0,1}{k1,k2}).^2)] ;
    end
    end
    end
    end
        Var=[Var,sum(var)] ;
end
        k=Distance(1,find(Var(1,:)==min(Var))) ; k=k(1) ;
                
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence{18,1}{k0,1}{k1,k2})==3
        Sequence{18,1}{k0,1}{k1,k2}=Sequence{18,1}{k0,1}{k1,k2}+[k,0,0] ;
    end
    end
    end
    end
end
        Sequence0{18,1}{1,I+size(EMBRYO,1)}=Sequence{18,1}                           ;
end
        save(['Normalization5\WorkSpace_Eggshell0_',num2str(n)],'Sequence0','-v7.3') ;
end



  % Mission 2 : Average Position & Positional Variation
  
for n=1:57
    load(['Normalization5\WorkSpace_Eggshell0_',num2str(n)])  ; load(['Normalization4\WorkSpace_Criteria',num2str(n)]) ;
    Positive=CellName ; Position=CellName     ; n
for k0=1:length(CellName)
for k1=1:size(CellName{k0,1},1)
for k2=1:size(CellName{k0,1},2)
    Positive{k0,1}{k1,k2}=0 ; Position{k0,1}{k1,k2}=NaN       ;
if  Criteria{1,n}{k0,1}{k1,k2}==1
    position=[] ;
    for I=1:size(EMBRYO,1)+length(HLH1)
        position=[position;Sequence0{18,1}{1,I}{k0,1}{k1,k2}] ;
    end
        Positive{k0,1}{k1,k2}=sqrt(sum(sum((position-ones(size(position,1),1)*mean(position))'.^2))/size(position,1))  ;
        Position{k0,1}{k1,k2}=mean(position) ;
end
end
end
end
        Criteria{3,n}=Positive ; Criteria{4,n}=Position       ;
        save(['Normalization5\WorkSpace_Criteria',num2str(n)],'Criteria','-v7.3') ;
end