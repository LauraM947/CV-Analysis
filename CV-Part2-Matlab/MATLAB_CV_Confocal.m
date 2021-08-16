
%% Coefficient of variation for confocal images
% Code by Chiara VICARIO 2019-02
% Further annotated & Readme file by Victoria NEGUEMBOR 2021-06


clear all
close all 

edges = 0:1:255;

%% Parameter to add

% add the correct categories number and names
categ = 2; 
Folder = 'C:\Users\';
categname = {'Condition1','Condition2'};

for m = 1:categ
   
    % load .tif files 
    I{m,1} = uipickfiles;
    
end

% % choose colors by checking/unchecking prefered settings

% colors = [[0.6314    0.0039    0.7569];  [0.5176    0.8902    0.0039]];
% % magenta Condition1    %% green Condition2
 colors = colormap(jet(2));

for m = 1:categ
    
    data = I {m,1};
   
    
for k = 1: length(data)
    
    IM = imread(data{1,k});
    IM = double(IM);
    IM(IM==0) = NaN;
    
% figure, imagesc(IM)
% STD = stdfilt(IM);
% figure, imagesc(STD)

    [N{m}(:,k),~] = histcounts(IM(:),edges,'normalization','probability');
    figure(1)
    plot(edges(:,1:end-1),N{m}(:,k), 'color',colors(m,:))
    hold all

    num = nanstd(IM(:));
    den = nanmean(IM(:));

    CV(m,k) = num./den;

end
end


xlabel('pixel intensity','fontsize',14)
ylabel('probability','fontsize',14)
set(gca,'FontSize',14)
set(gcf,'color','w');
grid on 
%Save figure of individual cells
saveas(gcf,strcat(Folder,'\CV_Histplot_individualcells'),'fig')

%Uncomment next line if you also want to save a .bmp of the figure. 
%saveas(gcf,strcat(Folder,'\CV_Histplot_individualcells'),'bmp')


% Save output analysis result
% ATT no headers!!! 
CVtable = array2table(CV');
CVtable.Properties.VariableNames  = categname;
writetable(CVtable,strcat(Folder,'\CVanalysis_result.xlsx'));

%[file,path] = uiputfile('*.xls');
%filename = fullfile(path,file);
%xlswrite(filename,CV');
%save(Folder,'\CV_Confocal_analysis.mat');
%saveas(gcf,'Density_Results_Histplot','fig')
%saveas(gcf,strcat(Folder,'\Density_Results_Histplot_individualcells','fig')



% Compute meand and SE of px histogram
N_mean = cellfun(@(x)mean(x,2), N,'UniformOutput',false);
N_std = cellfun(@(x)std(x,1,2)./sqrt(size(x,2)), N,'UniformOutput',false);
% Plot shade error bar
for m =1 :categ
    figure(2)
    hold all
%     plot(edges(:,1:end-1),N_mean{1,m},'color',colors(m,:),'LineWidth',2)
    
    lineProps.width = 2;
    lineProps.edgestyle = '-';
    lineProps.col = {colors(m,:)};
    mseb(edges(:,1:end-1),N_mean{1,m}',N_std{1,m}',lineProps,0.9);  

end
hold on


legend('Condition1','Condition2')


xlabel('pixel intensity','fontsize',20)
ylabel('probability','fontsize',20)
set(gca,'FontSize',18)
set(gcf,'color','w');
grid on 
box off
xlim([0 255])


%Save the other figure
saveas(gcf,strcat(Folder,'\CV_Histplot_smooth'),'fig')

%Uncomment next line if you also want to save a .bmp of the figure. 
%saveas(gcf,strcat(Folder,'\CV_Histplot_smooth'),'bmp')

%Save Workspace
save(strcat(Folder,'\CVanalysis_workspace.mat'));
%saveas(gcf,'Density_Results_Histplot_smooth','fig')


%Uncomment this part (CTRL+T or command+T for Mac) if you want to plot the boxplot in matlab 
boxplotpregunta = input('Do Boxplot (Y/N)>','s');

if boxplotpregunta=='Y';

    CV_2 = reshape(CV',[1,numel(CV)]);

    CV_3 = CV_2;

    CV_3(CV_2==0) = nan;

group =0 ;
for k = 1 : size(CV,1)
    group = horzcat(group,k*ones(1,size(CV,2)));
    
end
group = group(2:end);
group(isnan(CV_3)) = nan;

group = group(~isnan(group));
CV_4 = CV_3(~isnan(CV_3));

figure(3)
a = boxplot(CV_4,group,'Labels',categname,'Colors',colors);
set(a,'LineWidth',2);
% set(gca,'linewidth',3)
set(gca,'box','off')
set(gca,'FontSize',14)
set(gcf,'color','w');
ylabel('Coefficient of Variation')

for m = 1 :categ
    
    figure(3)
    hold all
    scatter(m*ones(1,size(CV,2))+0.4,CV(m,:),50,colors(m,:),'filled');
    
end

%Save figure boxplot
saveas(gcf,strcat(Folder,'\Coefficient_Variation_Boxplot'),'fig')

else
end


