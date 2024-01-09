% Hasreet Gill 03/23/2022
%% Quantification using Delaunay Triangulation and Voronoi Cells
% This code uses Delaunay triangle edge lengths and Voronoi cell areas to
% quantify spacing variability in 2D patterns.
clear
close all

%% Choose image stack
[filename, pathname] = uigetfile({'*.jpg;*.tif;*.png;*.gif;*.lsm','All Image Files';...
      '*.*','All Files' },'mytitle');
names={'UV1','UV2','UV3','rhodopsin1','rhodopsin2','rhodopsin3','blue1','blue2','blue3','green1','green2','green3','red1','red2','red3'};
% names={'hindgut1','hindgut2','hindgut3','hindgut4','hox1','hox2','hox3','hox4','mid1','mid2','mid3','mid4'}
% names={'ordered','clustered'}
lim=length(imread(filename,1));

%% SKIP IF LOADING STRUCT
% Identify and store center points
for i=2 %1:length(names)
    im{i}=imread(filename,i);

    figure, imshow(im{i})

    % Choose points either automatically or manually
    choice=questdlg('how will you choose points?', ...
	'segmentation type', ...
    'manual','auto','auto');
    
    switch choice
    case 'auto'
        close all

        % Binarize, extract center points & edge points
        im2=imgaussfilt(im{i},1);
        bw1=im2bw(im2);
        se=strel('disk',2);
        bw2=imclose(bw1,se);
        rp=regionprops(bw2,'centroid');
        centroids{i}=cat(1,rp.Centroid);

        figure, imshow(im{i})
        hold on
        scatter(centroids{i}(:,1),centroids{i}(:,2),'r.')

        choice=questdlg('edit points?', ...
        'edit', ...
        'no','yes','yes');
        switch choice
            case 'no'
                continue
            case 'yes'
            % Delete wrong points
            msgbox('choose incorrect points');
            pause(60)
            centroids{i}=setdiff(centroids{i},vertcat(cursor_info.Position),'rows');
            close all
    
            % Add correct points
            figure, imshow(im{i})
            hold on
            scatter(centroids{i}(:,1),centroids{i}(:,2),'r.')
            msgbox('choose missing points');
            [misspts(:,1),misspts(:,2)]=getpts;
            centroids{i}=vertcat(centroids{i},misspts);
        end

    case 'manual'
        % Manually choose points, double click final point
        [centroids{i}(:,1),centroids{i}(:,2)]=getpts
    end
    clearvars -except centroids filename names im lim
end

%% Calculate & store triangle edge lengths and cell areas
close all
imstruct=DTVoronoi(im,names,centroids,lim);
imstruct.centroids=centroids;
%%
for y=1:3:length(names)
a=[imstruct.cellareas{1,y};imstruct.cellareas{1,y+1};imstruct.cellareas{1,y+2}];
% a=datasample(a,350,'Replace',false);
e=[imstruct.edges{1,y};imstruct.edges{1,y+1};imstruct.edges{1,y+2}];
% e=datasample(e,350,'Replace',false);
Na=floor(0.75*length(a));
Ne=floor(0.75*length(e));
    for x=1:10000
    sampa=datasample(a,Na,'Replace',false);
    sampe=datasample(e,Ne,'Replace',false);
    imstruct.areasampvar{y}(x,:)=var(sampa);
    imstruct.edgesampvar{y}(x,:)=var(sampe);
    end
end
save(['imstruct_',filename,'scale_06.29.22.mat'],'imstruct')

% RESTART HERE IF LOADING STRUCT

% % Plot all points on overlaid image
% imblank=uint8(zeros(lim,lim));
% figure, imshow(imblank)
% hold on
% scatter(centroids{1}(:,1),centroids{1}(:,2),'m.') % UV
% scatter(centroids{2}(:,1),centroids{2}(:,2),'b.') % Blue
% scatter(centroids{3}(:,1),centroids{3}(:,2),'w.') % Rhodopsin
% scatter(centroids{4}(:,1),centroids{4}(:,2),'g.') % Green
% scatter(centroids{5}(:,1),centroids{5}(:,2),'r.') % Red

% % Calculate minimum distance from all points overlaid 
% % SOME OPSINS OVERLAP 100% - NOT THE BEST METHOD...
% 
allcents=vertcat(imstruct.centroids{:});
% mindistall=sort(pdist(allpts));
% mindistall=mindistall(1);

for i=1:3:length(names)
    dist1=sort(pdist(imstruct.centroids{i}));
    dist2=sort(pdist(imstruct.centroids{i+1}));
    dist3=sort(pdist(imstruct.centroids{i+2}));
    mindist{i}=(dist1(1)+dist2(1)+dist3(1))/3;
    mindist{i+1}=mindist{i};
    mindist{i+2}=mindist{i};
end

% Find minimum distance from manually chosen points
overlay=imread('BLUE_3crop.tif');
figure, imshow(overlay)
[pts(:,1),pts(:,2)]=getpts;
mandist=pdist(pts);
mandist=mandist(1);

% Random image, points, triangles and cells
randpts=[]
randstruct=[]
D=[]
rname={'randUV1','randUV2','randUV3','randrho1','randrho2','randrho3','randblue1','randblue2','randblue3','randgreen1','randgreen2','randgreen3','randred1','randred2','randred3'};
for i=1:length(rname);
    randim{i}=uint8(zeros(lim,lim));
    numpts{i}=length(imstruct.centroids{i});
    [randpts{i}(:,1),randpts{i}(:,2),D] = GetPointsRandom(numpts{i},lim,lim,mindist{i});
end
% totalpts=length(allcents);
meanpts=floor(length(allcents)/length(names));


% randim{6}=uint8(zeros(lim,lim)); % Mean number of points with no minimum distance between them
% rname{6}={'mean_nomin'};
% [randpts{6}(:,1),randpts{6}(:,2),D] = GetPointsRandom(meanpts,lim,lim,1); 
% 
% randim{7}=uint8(zeros(lim,lim)); % Total number of points with no minimum distance between them
% rname{7}={'total_nomin'};
% [randpts{7}(:,1),randpts{7}(:,2),D] = GetPointsRandom(totalpts,lim,lim,1); 

randim{16}=uint8(zeros(lim,lim)); % Mean number of points with no minimum distance between them
rname{16}={'globrand1'};
[randpts{16}(:,1),randpts{16}(:,2),D] = GetPointsRandom(meanpts,lim,lim,mandist); 
randim{17}=uint8(zeros(lim,lim)); % Mean number of points with no minimum distance between them
rname{17}={'globrand2'};
[randpts{17}(:,1),randpts{17}(:,2),D] = GetPointsRandom(meanpts,lim,lim,mandist); 
randim{18}=uint8(zeros(lim,lim)); % Mean number of points with no minimum distance between them
rname{18}={'globrand3'};
[randpts{18}(:,1),randpts{18}(:,2),D] = GetPointsRandom(meanpts,lim,lim,mandist); 

% randim{9}=uint8(zeros(lim,lim)); % Total number of points with no minimum distance between them
% rname{9}={'total_manualmin'};
% [randpts{9}(:,1),randpts{9}(:,2),D] = GetPointsRandom(totalpts,lim,lim,mandist); 
% 
% randim{10}=uint8(zeros(lim,lim)); % Manual minimum distance, 100 points
% rname{10}={'100pts'};
% pts100=100;
% [randpts{10}(:,1),randpts{10}(:,2),D] = GetPointsRandom(pts100,lim,lim,mandist);
% 
% randim{11}=uint8(zeros(lim,lim)); % Manual minimum distance, 500 points
% rname{11}={'500pts'};
% pts500=500;
% [randpts{11}(:,1),randpts{11}(:,2),D] = GetPointsRandom(pts500,lim,lim,mandist);

randstruct=DTVoronoi(randim,rname,randpts,lim);
randstruct.centroids=randpts;

% Na=floor(0.75*(length(randstruct.cellareas{1})));
% Ne=floor(0.75*(length(randstruct.edges{1})));
for i=1:3:length(rname)
%    
    ra=[randstruct.cellareas{i};randstruct.cellareas{i+1};randstruct.cellareas{i+2}];
    re=[randstruct.edges{i};randstruct.edges{i+1};randstruct.edges{i+2}];
%     ra=datasample(ra,350,'Replace',false);
%     re=datasample(re,350,'Replace',false);
    rNa=floor(0.75*length(ra));
    rNe=floor(0.75*length(re));
    for x=1:10000
        rsampa=datasample(ra,rNa,'Replace',false);
        rsampe=datasample(re,rNe,'Replace',false);
        randstruct.areasampvar{i}(x,:)=var(rsampa);
        randstruct.edgesampvar{i}(x,:)=var(rsampe);
    end
end
save(['randstruct_',filename,'scale_07.01.22.mat'],'randstruct')

% Boxplot
namesplot={'UV','UV','UV','Rhodopsin','Rhodopsin','Rhodopsin','Blue','Blue','Blue','Green','Green','Green','Red','Red','Red'}
namesplot2={'UV','UV','UV','Green','Green','Green','Blue','Blue','Blue','Rhodopsin','Rhodopsin','Rhodopsin','Red','Red','Red'}
for i=1:15
vcellnames2{i}=repmat({namesplot2{i}},length(imstruct.cellareas{i}),1);
cellareas2{i}=imstruct.cellareas{i};
end
figure,
vcelldata=vertcat(cellareas2{:});
vcellxlabels=vertcat(vcellnames2{:});
boxplot(vcelldata,vcellxlabels)
ylabel('Voronoi cell areas (pix^2)')
% saveas(gcf,'Voronoi_cell_areas_boxplot.tif')

%%
close all
v=violinplot(vcelldata,vcellxlabels,'ShowBox',false,'ShowWhiskers',false,'ShowData',false, ...
    'ViolinColor',[0.4940 0.1840 0.5560; 0 1 0; 0 0 1;0.8500 0.3250 0.0980;  1 0 0], ...
    'GroupOrder',{'UV','Green','Blue','Rhodopsin','Red'},'EdgeColor',[0 0 0], ...
    'MedianColor',[0.3 0.3 0.3])
ylabel('Voronoi cell areas')
x0=10;
y0=10;
width=1000;
height=600;
set(gcf,'position',[x0,y0,width,height])
ax=gca;
ax.FontSize=24
ax.LineWidth=1.8
ax.Box='off'
print -bestfit
saveas(gcf,'Voronoi_cell_areas_violinplot.pdf')
%%

% Just histogram
% Voronoi cells
close all
colors={[0.4940 0.1840 0.5560],[],[],[0.8500 0.3250 0.0980],[],[],'blue',[],[],'green',[],[],'red',[],[]};
figure
histogram(randstruct.areasampvar{16},'facecolor','black','facealpha',.5,'edgecolor','none')
hold on
for z=1:3:15
histogram(imstruct.areasampvar{z},'facecolor',colors{z},'facealpha',.3,'edgecolor','none')
pd=fitdist(imstruct.areasampvar{z},'Normal');
% 
% histogram(randstruct.areasampvar{z},'facecolor','black','facealpha',.3,'edgecolor','none')
% pd=fitdist(randstruct.areasampvar{z},'Normal');

% legend(namesplot{z},[namesplot{z}, ' Random'])

xlabel('Variance')
ylabel('Frequency')
title('Voronoi cell areas')

x0=10;
y0=10;
width=1000;
height=400;
set(gcf,'position',[x0,y0,width,height])

hold on
% saveas(gcf,[namesplot{z},'_vs_',namesplot{z},'random.tif'])
end
% legend('Random','UV','Rhodopsin','Blue','Green','Red')
%% 

% saveas(gcf,'all_opsins_nosubsamp_scale.tif')
%%
%%
close all
% s = settings;
% s.matlab.fonts.editor.code.Name.TemporaryValue = 'Arial'
for z=1:3:length(names);
    figure
%     histogram(randstruct.areasampvar{16},30,'facecolor','black','edgecolor','none')
%     pd=fitdist(randstruct.areasampvar{16},'Normal');
%     hold on
    histogram(randstruct.areasampvar{z},30,'facecolor',colors{z},'facealpha',.1,'edgecolor','none')
    pd=fitdist(randstruct.areasampvar{z},'Normal');
    hold on
    histogram(imstruct.areasampvar{z},30,'facecolor',colors{z},'facealpha',.3,'edgecolor','none')
    pd=fitdist(imstruct.areasampvar{z},'Normal');
    
  
    
    histogram(randstruct.areasampvar{z},30,'facecolor','black','facealpha',.3,'edgecolor','none')
    pd=fitdist(randstruct.areasampvar{z},'Normal');

    ax=gca;
    ax.FontSize=24
%     ax.FontName='Times'
    ax.LineWidth=1.8
    ax.Box='off'
%     if z==1
%     lgd=legend('Global','Cell-specific','FontSize',20,'Box','off','FontName','Times')
%     lgd.Title.String = 'Random Distribution';
%         lgd.Title.FontSize = 22;
%     lgd.Title.FontWeight='normal'
%     else
% %             legend('Random','FontSize',14,'Box','off','FontName','Times')
% 
%     end
    if z==13
        xlabel('Variance','FontSize',26,'FontName','Arial')
         x0=10;
    y0=10;
    width=750;
    height=295;
    set(gcf,'position',[x0,y0,width,height])
                    
    else
        set(gca,'xticklabel',{[]})
         x0=10;
    y0=10;
    width=1000;
    height=250;
    set(gcf,'position',[x0,y0,width,height])
    end
    ylabel('Frequency','FontSize',26,'FontName','Arial')
   
    xlim([0 300000]) 
    ylim([0 1200])
    hold off
    saveas(gcf,[namesplot{z}, '_spec_random.pdf'])
end
%%



figure
for z=6:9; %length(rname);
histogram(randstruct.areasampvar{z},'facealpha',.3,'edgecolor','none')
pd=fitdist(randstruct.areasampvar{z},'Normal');
hold on
end
legend('mean min-','total min-','mean min+','total min+')
xlabel('Variance')
ylabel('Frequency')
title('Voronoi cell areas')
hold off
saveas(gcf,'global_random_options.tif')
% % Triangles
% figure
% histogram(randstruct.edgesampvar{1},'facecolor','black','facealpha',.5,'edgecolor','none')
% hold on
% for z=1:length(colors)
% histogram(imstruct.edgesampvar{z},'facecolor',colors{z},'facealpha',.3,'edgecolor','none')
% pd=fitdist(imstruct.edgesampvar{z},'Normal');
% end
% legend(['random',names])
% xlabel('Variance')
% ylabel('Frequency')
% title('Delaunay triangulation edge lengths')
% saveas(gcf,'Voronoi cells variance distributions.tif')

for i=1:5;
figure
histogram(imstruct.areasampvar{i},'facecolor',colors{i},'facealpha',.3,'edgecolor','none')
pd=fitdist(imstruct.areasampvar{i},'Normal');
hold on
histogram(randstruct.areasampvar{i},'facecolor','black','facealpha',.3,'edgecolor','none')
pd=fitdist(randstruct.areasampvar{i},'Normal');
legend(names{i},['random ', names{i}])
xlabel('Variance')
ylabel('Frequency')
title('Voronoi cell areas')
x0=10;
y0=10;
width=1000;
height=400;
set(gcf,'position',[x0,y0,width,height])
hold off
saveas(gcf,[names{i},' vs random ',names{i},'.tif'])
end

%% Histogram + fit
% Voronoi cells
hr=histfit(randstruct.areasampvar{1},[],'normal')
hr(1).FaceColor='black';
hr(1).FaceAlpha=0.3;
hr(1).EdgeColor='none';
hr(2).Color='black';
hr(2).LineWidth=1;
hold on

% h1=histfit(imstruct.areasampvar{1},[],'normal')
% h1(1).FaceColor=colors{1};
% h1(1).FaceAlpha=0.3;
% h1(1).EdgeColor='none';
% h1(2).Color=colors{1};
% h1(2).LineWidth=1;

h2=histfit(imstruct.areasampvar{2},[],'normal')
h2(1).FaceColor=colors{2};
h2(1).FaceAlpha=0.3;
h2(1).EdgeColor='none';
h2(2).Color=colors{2};
h2(2).LineWidth=1;
% 
% h3=histfit(imstruct.areasampvar{3},[],'normal')
% h3(1).FaceColor=colors{3};
% h3(1).FaceAlpha=0.3;
% h3(1).EdgeColor='none';
% h3(2).Color=colors{3};
% h3(2).LineWidth=1;
% 
% h5=histfit(imstruct.areasampvar{5},[],'normal')
% h5(1).FaceColor=colors{5};
% h5(1).FaceAlpha=0.3;
% h5(1).EdgeColor='none';
% h5(2).Color=colors{5};
% h5(2).LineWidth=1;
% 
% h4=histfit(imstruct.areasampvar{4},[],'normal')
% h4(1).FaceColor=colors{4};
% h4(1).FaceAlpha=0.3;
% h4(1).EdgeColor='none';
% h4(2).Color=colors{4};
% h4(2).LineWidth=1;

legend('random',' ','blue',' ')
%'blue',' ','rhodopsin',' ','red',' ','green',' ')
xlabel('Variance')
ylabel('Frequency')
title('Voronoi cell areas')

x0=10;
y0=10;
width=1000;
height=400;
set(gcf,'position',[x0,y0,width,height])
hold off
saveas(gcf,'blue_vs_random_withfit.tif')

% % Triangles
% hr=histfit(randstruct.edgesampvar{1},[],'normal')
% hr(1).FaceColor='black';
% hr(1).FaceAlpha=0.3;
% hr(1).EdgeColor='none';
% hr(2).Color='black';
% hr(2).LineWidth=1;

% hold on
% 
% h1=histfit(imstruct.edgesampvar{1},[],'normal')
% h1(1).FaceColor=colors{1};
% h1(1).FaceAlpha=0.3;
% h1(1).EdgeColor='none';
% h1(2).Color=colors{1};
% h1(2).LineWidth=1;

% h2=histfit(imstruct.edgesampvar{2},[],'normal')
% h2(1).FaceColor=colors{2};
% h2(1).FaceAlpha=0.3;
% h2(1).EdgeColor='none';
% h2(2).Color=colors{2};
% h2(2).LineWidth=1;

% h3=histfit(imstruct.edgesampvar{3},[],'normal')
% h3(1).FaceColor=colors{3};
% h3(1).FaceAlpha=0.3;
% h3(1).EdgeColor='none';
% h3(2).Color=colors{3};
% h3(2).LineWidth=1;

% h5=histfit(imstruct.edgesampvar{5},[],'normal')
% h5(1).FaceColor=colors{5};
% h5(1).FaceAlpha=0.3;
% h5(1).EdgeColor='none';
% h5(2).Color=colors{5};
% h5(2).LineWidth=1;

% h4=histfit(imstruct.edgesampvar{4},[],'normal')
% h4(1).FaceColor=colors{4};
% h4(1).FaceAlpha=0.3;
% h4(1).EdgeColor='none';
% h4(2).Color=colors{4};
% h4(2).LineWidth=1;

% legend('random',' ','UV',' ','blue',' ','rhodopsin',' ','red',' ','green',' ')
% xlabel('Variance')
% ylabel('Frequency')
% title('Delaunay triangulation edge lengths')
% saveas(gcf,'Delaunay triangles variance distributions_with fit.tif')

%% t-tests
% Each vs. random
[t.UVr,p.UV]=ttest2(randstruct.areasampvar{1},imstruct.areasampvar{1},'Vartype','unequal');
[t.bluer,p.bluep]=ttest2(randstruct.areasampvar{1},imstruct.areasampvar{2},'Vartype','unequal');
[t.rhor,p.rhop]=ttest2(randstruct.areasampvar{1},imstruct.areasampvar{3},'Vartype','unequal');
[t.greenr,p.greenr]=ttest2(randstruct.areasampvar{1},imstruct.areasampvar{4},'Vartype','unequal');
[t.redr,p.redr]=ttest2(randstruct.areasampvar{1},imstruct.areasampvar{5},'Vartype','unequal');
[t.globrand,p.globrand]=ttest2(randstruct.areasampvar{2},randstruct.areasampvar{8});

% Pairwise
[t.rhoblue,p.rhoblue]=ttest2(imstruct.areasampvar{3},imstruct.areasampvar{2},'Vartype','unequal');
[t.rhogreen,p.rhogreen]=ttest2(imstruct.areasampvar{3},imstruct.areasampvar{4},'Vartype','unequal');
[t.bluegreen,p.bluegreen]=ttest2(imstruct.areasampvar{2},imstruct.areasampvar{4});
[t.blueUV,p.blueUV]=ttest2(imstruct.areasampvar{1},imstruct.areasampvar{2},'Vartype','unequal');
[t.greenred,p.greenred]=ttest2(imstruct.areasampvar{4},imstruct.areasampvar{5},'Vartype','unequal');


%% Quantification using Point Pattern Analysis
close all
    for n=1:length(imstruct.centroids)
%         figure
        for o=1:100
%             randmat=lim*rand(size(imstruct.centroids{n}));
              randmat=[];
%                 [randmat(:,1),randmat(:,2),D] = GetPointsRandom(length(imstruct.centroids{n}),lim,lim,mandist);
                                [randmat(:,1),randmat(:,2),D] = GetPointsRandom(meanpts,lim,lim,mindist{n});

            xK=[0:2:70];
%                 [ordX,ordY] = meshgrid(0:lim/(sqrt(meanpts)):lim)
%                 ordered=[ordX(:),ordY(:)];
            boxrand=[min(randmat(:,1)), max(randmat(:,1)), min(randmat(:,2)), max(randmat(:,2))];
            boxdata=[min(imstruct.centroids{n}(:,1)), max(imstruct.centroids{n}(:,1)), min(imstruct.centroids{n}(:,2)), max(imstruct.centroids{n}(:,2))];
%             boxdata=[min(ordered(:,1)), max(ordered(:,1)), min(ordered(:,2)), max(ordered(:,2))];
method=1;
            [K,CSR]=kfunction(imstruct.centroids{n}, randmat, xK, boxrand, boxdata, method);
            imstruct.CSR{n,:}(o,:)=CSR;
%             plot(xK',CSR,'b')
%             title(names{n});
%             hold on
% %             plot(xK',K,'red')
        end
          imstruct.ripK{:,n}=K;
%         final.Kmean{n}=mean(horzcat(imstruct.ripK{n,:}),2);
%         final.CSRmean{n}=mean(horzcat(imstruct.CSR{n,:}),2);
% % %         final.Kstd{n}=std(horzcat(imstruct.ripK{n,:}),0,2);
%         final.CSRstd{n}=std(horzcat(imstruct.CSR{n,:}),0,2);
        
%         figure,
%         curve1 = final.CSRmean{n} + final.CSRstd{n};
%         curve2 = final.CSRmean{n} - final.CSRstd{n};
%         x2 = [xK', fliplr(xK')];
%         inBetween = [curve1, fliplr(curve2)];
%         fill(x2, inBetween, 'g');
%         hold on
%         plot(xK', final.CSRmean{n}, 'b', 'LineWidth', 2)
%         hold on

%             saveas(gcf,[names{n},'_mindist_csrplot'],'tiff')
    end
    %%
%     close all
%     s = settings;
% s.matlab.fonts.editor.code.Name.TemporaryValue = 'Times New Roman'
% 
%     ppcolors={[0.4940 0.1840 0.5560]; [0.4940 0.1840 0.5560];[0.4940 0.1840 0.5560];[0.8500 0.3250 0.0980]; [0.8500 0.3250 0.0980];[0.8500 0.3250 0.0980];[0 0 1]; [0 0 1];[0 0 1];[0 1 0];[0 1 0];[0 1 0]; [1 0 0]; [1 0 0]; [1 0 0]};
% %    figure
%        CSRy=[imstruct.CSR{1}',imstruct.CSR{2}',imstruct.CSR{3}']';
%  shadedErrorBar(xK, mean(CSRy), [2.*std(CSRy);2.*std(CSRy)])
% %          s.patch.FaceColor = [0.7 0.7 0.7];
%      plot(xK,mean(CSRy,1),'LineWidth',1.2,'Color','k')

    for i=1:%:3:length(imstruct.centroids)
    figure

%     datay=[imstruct.ripK{i},imstruct.ripK{i+1},imstruct.ripK{i+2}]';
    datay1=[imstruct.ripK{1}]';
    datay2=[imstruct.ripK{2}]';
%     CSRy=[imstruct.CSR{i}',imstruct.CSR{i+1}',imstruct.CSR{i+2}']';
        CSRy=[imstruct.CSR{1}']';

%     SEM = std(CSRy)/sqrt(length(CSRy));               % Standard Error
%     ts = tinv([0.025  0.975],length(CSRy)-1);      % T-Score
%     CI = [(mean(CSRy)+(ts(1).*SEM));(mean(CSRy)+(ts(2).*SEM))];                      % Confidence Intervals
%      s=shadedErrorBar(xK, datay,{@mean,@std},'lineprops','w')
%      s.patch.FaceColor = ppcolors{i};
%      hold on
%      plot(xK,mean(datay,1),'Color',ppcolors{i},'LineWidth',1.8)
     plot(xK,mean(datay1,1),'Color','r','LineWidth',2.2)
     hold on
     plot(xK,mean(datay2,1),'Color','b','LineWidth',2.6)
%      plot(xK,mean(datay,1),'Color',ppcolors{i},'LineWidth',1.8)

%     shadedErrorBar(xK, mean(CSRy), [2.*std(CSRy);2.*std(CSRy)])
%          s.patch.FaceColor = [0.7 0.7 0.7];
     plot(xK,mean(CSRy,1),'LineWidth',2.2,'Color','k')
         xlim([0 15])
    ylim([0 455])

xlabel('Distance')
ylabel('K')
x0=10;
y0=10;
width=700;
height=500;
set(gcf,'position',[x0,y0,width,height])
ax=gca;
ax.FontSize=26
ax.LineWidth=2.6
ax.Box='off'
saveas(gcf,'PPschematic.tif')
% hold on
%     end
