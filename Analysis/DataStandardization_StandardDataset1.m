  % Mission 1 : Data Standardization of Standard Dataset 1

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
for n=1:size(S,2)
    load(['E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200319 - DataImport\Normalization5\WorkSpace_Eggshell0_',num2str(n)]) ;
    Sequence0(1:17,:)=[] ; Sequence0(end-1:end,:)=[] ; Dataset=Sequence0 ; n
    save(['E:\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200430 - DatasetReorganization\Standard Dataset 1\WorkSpace_Dataset_',num2str(S{2,n})],'Dataset','-v7.3') ;
end