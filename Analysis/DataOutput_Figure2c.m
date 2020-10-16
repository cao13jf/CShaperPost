  % Mission 1 : Data Output of Figure 2c
  
    load('D:\STORAGE\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200530 - CAO FigureUpdate_DATA\F1_Score_Summary.mat') ; Data=F1_matrix ;
    Name={'3DUNet';'B-CShaper';'CellProfiler';'CShaper';'FusionNet';'RACE';'SingleCellDetector'} ; figure('Position',[0 0 400 400])            ;
    File=dir('D:\STORAGE\Project 12 - Wnt & Membrane\CJF\Paper Revision - PLUSPLUS\Revision20200530 - CAO FigureUpdate_DATA') ; File([1,2,7,11],:)=[]           ;
    File=[File([1,4,5,6,7],:);File([2,3],:)] ; c=[[238,0,0];[255,192,0];[0,139,69];[112,197,184];[59,73,146];[156,0,232];[0,0,0]]/256          ;
    SourceData=cell(2+length(File),length([0.5:0.05:0.95])+1)   ; R0=1 ; C0=-3    ; SourceData{1,1}='F1 Score' ; SourceData{3,1}='Method'      ;
    SourceData{1,2}='IoU Threshold' ; Threshold=[0.5:0.05:0.95] ;
for iii=1:length(Threshold)
    SourceData{2,1+iii}=sprintf('%2.2f',Threshold(iii))         ;
end
for I=1:length(File)
    SourceData{3+I,1}=File(I).name  ;
for i=1:length(Name)
if  strcmp(File(I).name,Name{i,1})==1
    plot([0.5:0.05:0.95],Data(:,i),'-','linewidth',1.5,'color',c(I,:)) ; hold on  ;
    for iii=1:length(Threshold)
        SourceData{3+I,1+iii}=sprintf('%2.4f',Data(iii,i))      ;
    end
end
end
end
for I=1:length(File)
for i=1:length(Name)
if  strcmp(File(I).name,Name{i,1})==1
    plot([0.5:0.05:0.95],Data(:,i),'.','markersize',20,'color',c(I,:)) ; hold on  ;
end
end
end
    axis([0.5,1.0,0,1])                           ;
    x=xlabel({'\rm IoU threshold'})               ; set(x,'Fontname','arial','Fontsize',18) ;
    y=ylabel({'\it F\rm_{1} score'})              ; set(y,'Fontname','arial','Fontsize',18) ;
    set(gca,'FontSize',18,'Fontname','arial')     ; axis square                   ;
    set(gca,'xtick',[0.5,0.6,0.7,0.8,0.9,1.0],'FontSize',18,'Fontname','arial')   ;
    set(gca,'ytick',[0.0,0.2,0.4,0.6,0.8,1.0],'FontSize',18,'Fontname','arial')   ;
    xlswrite('Figure 2c.xls',SourceData,'Sheet1') ;