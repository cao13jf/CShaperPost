  % Data Standardization

  
  
  % Mission 1 : Eggshell Reorientation

for n=1:57
    load(['Normalization2\WorkSpace_00',num2str(n)])   ; load('WorkSpace_EMBRYO')    ; Sequence0=Sequence(:,n)          ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_CellName') ; load('WorkSpace_HLH1') ;
    Sequence0{20,1}=cell(1,length(EMBRYO))             ;
for I=1:size(EMBRYO,1) 
    Phi=[0:pi/180:pi-pi/180] ; I
    
  % YZ Plane
    c=hsv(size(Sequence0{19,1}{1,I},1))     ; point=[] ; Range=zeros(1,length(Phi))  ; Step=[-20:0.01:20]               ;
for i=1:size(Sequence0{19,1}{1,I},1)
    point=[point;[Sequence0{19,1}{1,I}{i,1}(:,3),Sequence0{19,1}{1,I}{i,1}(:,1)]]    ;
end
for i=1:length(Phi)
    Point=point(:,1)*tan(Phi(i))-point(:,2) ;
    point1=point(Point>0,:) ; point2=point(Point<0,:)  ;
    Max1=max(abs((point1(:,1)*tan(Phi(i))-point1(:,2))/sqrt(tan(Phi(i)).^2+1)))      ;
    Max2=max(abs((point2(:,1)*tan(Phi(i))-point2(:,2))/sqrt(tan(Phi(i)).^2+1)))      ;
    Range(i)=Max1+Max2      ;
end
    phi=Phi(find(Range==min(Range)))        ; phi=phi(1) ;
if  phi~=0
if  abs(phi-pi/2)>1e-5
    Point=point(:,1)*tan(phi)-point(:,2)    ; k=phi-pi/2 ;
    point1=point(Point>0,:) ; point2=point(Point<0,:)    ;
    Point1=abs((point1(:,1)*tan(phi)-point1(:,2))/sqrt(tan(phi).^2+1))               ;
    Point2=abs((point2(:,1)*tan(phi)-point2(:,2))/sqrt(tan(phi).^2+1))               ;
    P1=point1(find(Point1==max(Point1)),:)  ; P2=point2(find(Point2==max(Point2)),:) ;
    L1=Step*tan(phi)+P1(2)-P1(1)*tan(phi) ; L2=Step*tan(phi)+P2(2)-P2(1)*tan(phi)    ;
    Point=point(:,1)*(-1/tan(phi))-point(:,2)            ;
    point1=point(Point>0,:) ; point2=point(Point<0,:)    ;
    Point1=abs((point1(:,1)*(-1/tan(phi))-point1(:,2))/sqrt((-1/tan(phi)).^2+1))     ;
    Point2=abs((point2(:,1)*(-1/tan(phi))-point2(:,2))/sqrt((-1/tan(phi)).^2+1))     ;
    P1=point1(find(Point1==max(Point1)),:)  ; P2=point2(find(Point2==max(Point2)),:) ;
    L3=Step*(-1/tan(phi))+P1(2)-P1(1)*(-1/tan(phi)) ; L4=Step*(-1/tan(phi))+P2(2)-P2(1)*(-1/tan(phi))         ;
    p1=find(abs(L1-L3)==min(abs(L1-L3))) ; p2=find(abs(L1-L4)==min(abs(L1-L4)))      ;
    p3=find(abs(L2-L3)==min(abs(L2-L3))) ; p4=find(abs(L2-L4)==min(abs(L2-L4)))      ;
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence0{18,1}{1,I}{k0,1}{k1,k2})==3
        Sequence0{18,1}{1,I}{k0,1}{k1,k2}=Sequence0{18,1}{1,I}{k0,1}{k1,k2}-[(L1(p1)+L1(p2)+L2(p3)+L2(p4))/4,0,(Step(p1)+Step(p2)+Step(p3)+Step(p4))/4] ;
    end
    end
    end
    end
    for j=1:size(Sequence0{19,1}{1,I},1)
        Sequence0{19,1}{1,I}{j,1}(:,3)=Sequence0{19,1}{1,I}{j,1}(:,3)-(Step(p1)+Step(p2)+Step(p3)+Step(p4))/4 ;
        Sequence0{19,1}{1,I}{j,1}(:,1)=Sequence0{19,1}{1,I}{j,1}(:,1)-(L1(p1)+L1(p2)+L2(p3)+L2(p4))/4         ;
    end
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence0{18,1}{1,I}{k0,1}{k1,k2})==3
        y= Sequence0{18,1}{1,I}{k0,1}{k1,k2}(1,3)*cos(k)+Sequence0{18,1}{1,I}{k0,1}{k1,k2}(1,1)*sin(k)        ;   
        z=-Sequence0{18,1}{1,I}{k0,1}{k1,k2}(1,3)*sin(k)+Sequence0{18,1}{1,I}{k0,1}{k1,k2}(1,1)*cos(k)        ;
        Sequence0{18,1}{1,I}{k0,1}{k1,k2}(1,3)=y ; Sequence0{18,1}{1,I}{k0,1}{k1,k2}(1,1)=z                   ;
    end
    end
    end
    end
    for j=1:size(Sequence0{19,1}{1,I},1)
        y= Sequence0{19,1}{1,I}{j,1}(:,3)*cos(k)+Sequence0{19,1}{1,I}{j,1}(:,1)*sin(k) ;
        z=-Sequence0{19,1}{1,I}{j,1}(:,3)*sin(k)+Sequence0{19,1}{1,I}{j,1}(:,1)*cos(k) ;
        Sequence0{19,1}{1,I}{j,1}(:,3)=y         ; Sequence0{19,1}{1,I}{j,1}(:,1)=z    ;
    end
end
end
if  phi==0 || abs(phi-pi/2)<1e-5
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence0{18,1}{1,I}{k0,1}{k1,k2})==3
        Sequence0{18,1}{1,I}{k0,1}{k1,k2}=Sequence0{18,1}{1,I}{k0,1}{k1,k2}-[(max(point(:,2))+min(point(:,2)))/2,0,(max(point(:,1))+min(point(:,1)))/2] ;
    end
    end
    end
    end
    for j=1:size(Sequence0{19,1}{1,I},1)
        Sequence0{19,1}{1,I}{j,1}(:,3)=Sequence0{19,1}{1,I}{j,1}(:,3)-(max(point(:,1))+min(point(:,1)))/2 ;
        Sequence0{19,1}{1,I}{j,1}(:,1)=Sequence0{19,1}{1,I}{j,1}(:,1)-(max(point(:,2))+min(point(:,2)))/2 ;
    end
end
        point=[] ;
    for i=1:size(Sequence0{19,1}{1,I},1)
        point=[point;[Sequence0{19,1}{1,I}{i,1}(:,3),Sequence0{19,1}{1,I}{i,1}(:,1)]]   ;
    end
        Sequence0{20,1}{1,I}=[Sequence0{20,1}{1,I},(max(point(:,1))-min(point(:,1)))/2] ; % Compression

  % XZ Plane
    Phi=[-pi/4:pi/180:pi/4] ; point=[]      ; Range=zeros(1,length(Phi))   ; Step=[-30:0.01:30] ;
for i=1:size(Sequence0{19,1}{1,I},1)
    point=[point;[Sequence0{19,1}{1,I}{i,1}(:,2),Sequence0{19,1}{1,I}{i,1}(:,1)]]               ;
end
for i=1:length(Phi)
    Point=point(:,1)*tan(Phi(i))-point(:,2) ;
    point1=point(Point>0,:) ; point2=point(Point<0,:)    ;
    Max1=max(abs((point1(:,1)*tan(Phi(i))-point1(:,2))/sqrt(tan(Phi(i)).^2+1)))      ;
    Max2=max(abs((point2(:,1)*tan(Phi(i))-point2(:,2))/sqrt(tan(Phi(i)).^2+1)))      ;
    Range(i)=Max1+Max2      ;
end
    phi=Phi(find(Range==min(Range)))        ; phi=phi(1) ;
if  phi~=0
    Point=point(:,1)*tan(phi)-point(:,2)    ; k=phi      ;
    point1=point(Point>0,:) ; point2=point(Point<0,:)    ;
    Point1=abs((point1(:,1)*tan(phi)-point1(:,2))/sqrt(tan(phi).^2+1))               ;
    Point2=abs((point2(:,1)*tan(phi)-point2(:,2))/sqrt(tan(phi).^2+1))               ;
    P1=point1(find(Point1==max(Point1)),:)  ; P2=point2(find(Point2==max(Point2)),:) ;
    L1=Step*tan(phi)+P1(2)-P1(1)*tan(phi) ; L2=Step*tan(phi)+P2(2)-P2(1)*tan(phi)    ;
    Point=point(:,1)*(-1/tan(phi))-point(:,2)            ;
    point1=point(Point>0,:) ; point2=point(Point<0,:)    ;
    Point1=abs((point1(:,1)*(-1/tan(phi))-point1(:,2))/sqrt((-1/tan(phi)).^2+1))     ;
    Point2=abs((point2(:,1)*(-1/tan(phi))-point2(:,2))/sqrt((-1/tan(phi)).^2+1))     ;
    P1=point1(find(Point1==max(Point1)),:)  ; P2=point2(find(Point2==max(Point2)),:) ;
    L3=Step*(-1/tan(phi))+P1(2)-P1(1)*(-1/tan(phi)) ; L4=Step*(-1/tan(phi))+P2(2)-P2(1)*(-1/tan(phi)) ;
    p1=find(abs(L1-L3)==min(abs(L1-L3))) ; p2=find(abs(L1-L4)==min(abs(L1-L4)))      ;
    p3=find(abs(L2-L3)==min(abs(L2-L3))) ; p4=find(abs(L2-L4)==min(abs(L2-L4)))      ;
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence0{18,1}{1,I}{k0,1}{k1,k2})==3
        Sequence0{18,1}{1,I}{k0,1}{k1,k2}=Sequence0{18,1}{1,I}{k0,1}{k1,k2}-[(L1(p1)+L1(p2)+L2(p3)+L2(p4))/4,(Step(p1)+Step(p2)+Step(p3)+Step(p4))/4,0] ;
    end
    end
    end
    end
    for j=1:size(Sequence0{19,1}{1,I},1)
        Sequence0{19,1}{1,I}{j,1}(:,2)=Sequence0{19,1}{1,I}{j,1}(:,2)-(Step(p1)+Step(p2)+Step(p3)+Step(p4))/4 ;
        Sequence0{19,1}{1,I}{j,1}(:,1)=Sequence0{19,1}{1,I}{j,1}(:,1)-(L1(p1)+L1(p2)+L2(p3)+L2(p4))/4         ;
    end    
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence0{18,1}{1,I}{k0,1}{k1,k2})==3
        x= Sequence0{18,1}{1,I}{k0,1}{k1,k2}(1,2)*cos(k)+Sequence0{18,1}{1,I}{k0,1}{k1,k2}(1,1)*sin(k) ;   
        z=-Sequence0{18,1}{1,I}{k0,1}{k1,k2}(1,2)*sin(k)+Sequence0{18,1}{1,I}{k0,1}{k1,k2}(1,1)*cos(k) ;
        Sequence0{18,1}{1,I}{k0,1}{k1,k2}(1,2)=x ; Sequence0{18,1}{1,I}{k0,1}{k1,k2}(1,1)=z            ;
    end
    end
    end
    end
    for j=1:size(Sequence0{19,1}{1,I},1)
        x= Sequence0{19,1}{1,I}{j,1}(:,2)*cos(k)+Sequence0{19,1}{1,I}{j,1}(:,1)*sin(k) ;
        z=-Sequence0{19,1}{1,I}{j,1}(:,2)*sin(k)+Sequence0{19,1}{1,I}{j,1}(:,1)*cos(k) ;
        Sequence0{19,1}{1,I}{j,1}(:,2)=x ; Sequence0{19,1}{1,I}{j,1}(:,1)=z            ;
    end
end
if  phi==0
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence0{18,1}{1,I}{k0,1}{k1,k2})==3
        Sequence0{18,1}{1,I}{k0,1}{k1,k2}=Sequence0{18,1}{1,I}{k0,1}{k1,k2}-[(max(point(:,2))+min(point(:,2)))/2,(max(point(:,1))+min(point(:,1)))/2,0] ;
    end
    end
    end
    end
    for j=1:size(Sequence0{19,1}{1,I},1)
        Sequence0{19,1}{1,I}{j,1}(:,2)=Sequence0{19,1}{1,I}{j,1}(:,2)-(max(point(:,1))+min(point(:,1)))/2  ;
        Sequence0{19,1}{1,I}{j,1}(:,1)=Sequence0{19,1}{1,I}{j,1}(:,1)-(max(point(:,2))+min(point(:,2)))/2  ;
    end
end

  % Searching for Axes
    Phi=[-pi/4:pi/180:pi/4] ; point=[] ; Range=zeros(1,length(Phi)) ;
for i=1:size(Sequence0{19,1}{1,I},1)
    point=[point;[Sequence0{19,1}{1,I}{i,1}(:,2),Sequence0{19,1}{1,I}{i,1}(:,1)]]                         ;
end
for i=1:length(Phi)
    point0=[point(:,1)*cos(Phi(i))+point(:,2)*sin(Phi(i)),-point(:,1)*sin(Phi(i))+point(:,2)*cos(Phi(i))] ;
    AXIS1=[1:1/50:1.5]*(max(point0(:,1))-min(point0(:,1)))/2        ; AXIS2=[] ;
for Axis=AXIS1
    AXIS2=[AXIS2,sqrt(max((point0(:,2).^2)./(1-(point0(:,1).^2)/(Axis.^2))))]  ;
end
    Range(i)=min(AXIS1.*AXIS2)       ;
end
    phi=Phi(find(Range==min(Range))) ; k=phi(1)  ; 
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence0{18,1}{1,I}{k0,1}{k1,k2})==3
        x= Sequence0{18,1}{1,I}{k0,1}{k1,k2}(1,2)*cos(k)+Sequence0{18,1}{1,I}{k0,1}{k1,k2}(1,1)*sin(k)    ;   
        z=-Sequence0{18,1}{1,I}{k0,1}{k1,k2}(1,2)*sin(k)+Sequence0{18,1}{1,I}{k0,1}{k1,k2}(1,1)*cos(k)    ;
        Sequence0{18,1}{1,I}{k0,1}{k1,k2}(1,2)=x ; Sequence0{18,1}{1,I}{k0,1}{k1,k2}(1,1)=z               ;
    end
    end
    end
    end
    for j=1:size(Sequence0{19,1}{1,I},1)
        x= Sequence0{19,1}{1,I}{j,1}(:,2)*cos(k)+Sequence0{19,1}{1,I}{j,1}(:,1)*sin(k) ;
        z=-Sequence0{19,1}{1,I}{j,1}(:,2)*sin(k)+Sequence0{19,1}{1,I}{j,1}(:,1)*cos(k) ;
        Sequence0{19,1}{1,I}{j,1}(:,2)=x ; Sequence0{19,1}{1,I}{j,1}(:,1)=z            ;
    end
    point=[] ;
for i=1:size(Sequence0{19,1}{1,I},1)
    point=[point;[Sequence0{19,1}{1,I}{i,1}(:,2),Sequence0{19,1}{1,I}{i,1}(:,3),Sequence0{19,1}{1,I}{i,1}(:,1)]] ;
end
    AXIS1=[1:1/50:1.5]*(max(point(:,1))-min(point(:,1)))/2 ; AXIS2=[]       ;
for Axis=AXIS1
    AXIS2=[AXIS2,sqrt(max((point(:,3).^2)./(1-(point(:,1).^2)/(Axis.^2))))] ;
end
    Sequence0{20,1}{1,I}=[Sequence0{20,1}{1,I},AXIS1(find(AXIS1.*AXIS2==min(AXIS1.*AXIS2)))] ;
    Sequence0{20,1}{1,I}=[Sequence0{20,1}{1,I},AXIS2(find(AXIS1.*AXIS2==min(AXIS1.*AXIS2)))] ;
    Sequence0{20,1}{1,I}=[Sequence0{20,1}{1,I},sqrt(max((point(:,2).^2)./(1-(point(:,1).^2)/(Sequence0{20,1}{1,I}(2).^2)-(point(:,3).^2)/(Sequence0{20,1}{1,I}(3).^2))))] ;
end
    save(['Normalization3\WorkSpace_Eggshell_',num2str(n)],'Sequence0','-v7.3')              ;
end



  % Mission 2 : Size Rescaling

for n=1:57
    load('WorkSpace_EMBRYO') ; load('WorkSpace_HLH1') ; load(['Normalization3\WorkSpace_Eggshell_',num2str(n)]) ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_CellName') ; Parameter=[]   ; 
for I=1:size(EMBRYO,1)
    Parameter=[Parameter;Sequence0{20,1}{1,I}]        ;
end
for I=1:size(EMBRYO,1)    
    K=[Sequence0{20,1}{1,I}(3),Sequence0{20,1}{1,I}(2),Sequence0{20,1}{1,I}(1)]./mean([Parameter(:,3),Parameter(:,2),Parameter(:,1)]) ; I
    for k0=1:length(CellName)
    for k1=1:size(CellName{k0,1},1)
    for k2=1:size(CellName{k0,1},2)
    if  length(Sequence0{18,1}{1,I}{k0,1}{k1,k2})==3
        Sequence0{18,1}{1,I}{k0,1}{k1,k2}=Sequence0{18,1}{1,I}{k0,1}{k1,k2}./K                 ;
    end
    end
    end
    end
    for j=1:size(Sequence0{19,1}{1,I},1)
        Sequence0{19,1}{1,I}{j,1}=Sequence0{19,1}{1,I}{j,1}./(ones(size(Sequence0{19,1}{1,I}{j,1},1),1)*K)         ;
    end
    
  % Figure Plotting
    figure  ; c=hsv(size(Sequence0{19,1}{1,I},1)) ;
    plot(mean(Parameter(:,2))*cos([0:pi/180:2*pi]),mean(Parameter(:,3))*sin([0:pi/180:2*pi]),'k-','linewidth',2.5) ; hold on ;
for i=1:size(Sequence0{19,1}{1,I},1)
    plot(Sequence0{19,1}{1,I}{i,1}(:,2),Sequence0{19,1}{1,I}{i,1}(:,1),'.','markersize',2.5,'color',c(i,:))        ;
    hold on ; axis equal ; axis([-1.05*mean(Parameter(:,2)),1.05*mean(Parameter(:,2)),-1.05*mean(Parameter(:,3)),1.05*mean(Parameter(:,3))]) ;
end
    figure  ; c=hsv(size(Sequence0{19,1}{1,I},1)) ;
    plot(mean(Parameter(:,4))*cos([0:pi/180:2*pi]),mean(Parameter(:,3))*sin([0:pi/180:2*pi]),'k-','linewidth',2.5) ; hold on ;
for i=1:size(Sequence0{19,1}{1,I},1)
    plot(Sequence0{19,1}{1,I}{i,1}(:,3),Sequence0{19,1}{1,I}{i,1}(:,1),'.','markersize',2.5,'color',c(i,:))        ;
    hold on ; axis equal ; axis([-1.05*mean(Parameter(:,1)),1.05*mean(Parameter(:,1)),-1.05*mean(Parameter(:,3)),1.05*mean(Parameter(:,3))]) ;
end
    plot(mean(Parameter(:,1))*[-1,1,1,-1,-1],mean(Parameter(:,3))*[-1,-1,1,1,-1],'k-','linewidth',2.5) ; hold on   ;
end
    Sequence0{21,1}=Parameter ; save(['Normalization3\WorkSpace_Eggshell_',num2str(n)],'Sequence0','-v7.3')        ;
end



  % Mission 3 : Average Position

for n=1:57
    load('Normalization1\WorkSpace_Criteria') ; load(['Normalization3\WorkSpace_Eggshell_',num2str(n)]) ; load('WorkSpace_EMBRYO') ;
    load('E:\Project 1 - C.elegans Resource\11-Patch\Problem 1 - WT Cycle\WorkSpace_CellName') ; PositionAve=CellName              ;
for k0=1:length(CellName)
for k1=1:size(CellName{k0,1},1)
for k2=1:size(CellName{k0,1},2)
    PositionAve{k0,1}{k1,k2}=NaN              ;
if  Criteria{1,n}{k0,1}{k1,k2}>0
    g=[] ;
    for I=1:size(EMBRYO,1)
    if  length(Sequence0{18,1}{1,I}{k0,1}{k1,k2})==3
        g=[g;Sequence0{18,1}{1,I}{k0,1}{k1,k2}] ;
    end
    end
        PositionAve{k0,1}{k1,k2}=mean(g)        ;
end
end
end
end
        Criteria{2,n}=PositionAve ; save(['Normalization3\WorkSpace_Criteria',num2str(n)],'Criteria','-v7.3')                      ;
end