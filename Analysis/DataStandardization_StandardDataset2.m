  % Mission 1 : Data Standardization of Standard Dataset 2

  % Data Matrixing
    load('Tool/WorkSpace_Query')                      ;
for n=1:57
    load(['RawData/WorkSpace_Eggshell0_',num2str(n)]) ;
for I=1:length(Sequence0{19,1})
if  isempty(Sequence0{19,1}{1,I})==0
    Seg=zeros(256,114,184) ; save(['Process/WorkSpace_Process_',num2str(n),'_',num2str(I)],'n','-v7.3')       ;
for i=[1:1:size(Sequence0{19,1}{1,I},1)]
    dot=round((Sequence0{19,1}{1,I}{i,1}/0.25+ones(size(Sequence0{19,1}{1,I}{i,1},1),1)*[184/2,256/2,114/2])) ;
    for J=1:size(Query,1)
    if  strcmp(Sequence0{19,1}{1,I}{i,2},Query{J,1})==1
        for j=1:size(Sequence0{19,1}{1,I}{i,1},1)
        for x=dot(j,2)-1:dot(j,2)+1
        for y=dot(j,3)-1:dot(j,3)+1
        for z=dot(j,1)-1:dot(j,1)+1
            Seg(x,y,z)=Query{J,2}                       ;
        end
        end
        end
        end
    end
    end
end
        Cell=unique(Seg) ; Cell(Cell==0)=[] ; Seg0=Seg  ;
    for i=1:length(Cell)
        [X,Y,Z]=ind2sub(size(Seg0),find(Seg0==Cell(i))) ;
    for j=1:size([X,Y,Z],1)
        for x=X(j)-1:X(j)+1
        for y=Y(j)-1:Y(j)+1
        for z=Z(j)-1:Z(j)+1
        if  Seg0(x,y,z)==0
            Seg(x,y,z)=0                    ;
        end
        end
        end
        end
    end
    end
        save(['Result/Seg_',num2str(n),'_',num2str(I)],'Seg','-v7.3') ;
end
end
end

  % Correlation of Temporal Information
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200413 - PositionalVariability\WorkSpace_S') ;
    Sequence0=S ; Sequence0(2:end,:)=[] ; load('WorkSpace_Frame') ;
for i=1:size(Frame,1)
  % First Moment
if  strcmp(Frame{i,2},'/')==0
    for j=28:57
    if  strcmp(Frame{i,1},Sequence0{1,j})==1
        Sequence0{2,j}=str2num(Frame{i,2}(strfind(Frame{i,2},'[')+1:strfind(Frame{i,2},']')-1)) ;
    end
    end
end
end
for i=1:size(Frame,1)
  % Last  Moment
if  strcmp(Frame{i,3},'/')==0
    for j=1:27
    if  strcmp(Frame{i,1},Sequence0{1,j})==1
        Sequence0{2,j}=str2num(Frame{i,3}(strfind(Frame{i,3},'[')+1:strfind(Frame{i,3},']')-1)) ;
    end
    end
end
end
    S=Sequence0 ;
    
  % Data Generation
    load('E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200319 - DataImport\WorkSpace_EMBRYO') ;
for n=1:size(Sequence0,2)
    n
for I=1:size(EMBRYO,1)
    if  I+4<10
        III=['0',num2str(I+4)] ;
    end
    if  I+4>=10
        III=[    num2str(I+4)] ;
    end
    Path1=['E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200410 - Standardization\ServerComputation_PLUS\Result\Seg_',num2str(n),'_',num2str(I),'.mat']  ;
    Path2=['E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200430 - DatasetReorganization\Standard Dataset 2\Seg_',num2str(Sequence0{2,n}),'_',III,'.mat'] ;
    copyfile(Path1,Path2)      ;
end
end