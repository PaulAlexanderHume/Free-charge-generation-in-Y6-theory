%% Process Combined VOTCA Output
clear all
h_0=-6.70; e_0=-3.24; IE_0=-h_0; EA_0=e_0;% h_0=0; e_0=0;
E_n_1=importdata('E_neutral_1_768.txt');    E_n_2=importdata('E_neutral_769_1536.txt');  
E_e_1=importdata('E_electron_1_768.txt');   E_e_2=importdata('E_electron_769_1536.txt');
E_h_1=importdata('E_hole_1_768.txt');       E_h_2=importdata('E_hole_769_1536.txt'); 

dh_1=E_h_1-E_n_1;  IE1=IE_0+dh_1;          dh_2=E_h_2-E_n_2;  IE2=IE_0+dh_2; %IE2(202)=[];
de_1=E_e_1-E_n_1;  EA1=EA_0+de_1;          de_2=E_e_2-E_n_2;  EA2=EA_0+de_2; %EA2(202)=[];
IE=[IE1; IE2];
EA=[EA1; EA2];

%% Plotting Distributions (No distinction between positions)
binwidth=0.05;
ImageDPI=600; ImageSizeX=8.95; ImageSizeY=5;
ImageFontSize=10; linewidth=1;
h=figure('Units','centimeters','InnerPosition',[10 5 ImageSizeX ImageSizeY],...
           'PaperPosition',[10 5 ImageSizeX ImageSizeY])
h=histogram(EA,'Orientation','horizontal','BinWidth',binwidth,'Normalization','Probability','FaceColor',[150 0 0]/255)
hold on
h=histogram(-IE,'Orientation','horizontal','BinWidth',binwidth,'Normalization','Probability','FaceColor',[57 82 200]/255)
h=line([0 0.2],[-4.1 -4.1],'Color',[0 0 0],'LineWidth',1.0,'LineStyle','--','Color',[0 0 0])
% h=text(35,-5.4,'\downarrow','HorizontalAlignment','center')
% h=text(35,-5.1,'PESA','HorizontalAlignment','center')
% h=text(35,-4.35,'\uparrow','HorizontalAlignment','center')
% h=text(35,-4.65,'IPES','HorizontalAlignment','center')
% h=line([0 40],[-5.8 -5.8],'Color',[0 0 0],'LineWidth',1.0,'LineStyle','--','Color',[0.5 0.5 0.5])
h=line([0 0.2],[-5.65 -5.65],'Color',[0 0 0],'LineWidth',1.0,'LineStyle','--','Color',[0 0 0])
% xticks([0 10 20 30 pi 2*pi 3*pi])
% xticklabels({'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'})
xlabel('DOS','fontsize',ImageFontSize)
ylabel('Energy (eV)','fontsize',ImageFontSize)
legend({'EA','-IE'},'Location','best')
legend('boxoff')
h=text(-5,-3,'\bf b) \rm')
ylim([-6.5,-3])

%% Means and Onsets
'Mean IE'
mean(IE)
'Mean EA'
mean(EA)
'Mean IE - 2*sigma'
mean(IE)-2*std(IE)
'Mean EA - 2*sigma'
mean(EA)-2*std(EA)

%% Plotting Distributions (Distinguishing between positions)
binwidth=0.05;
ImageDPI=600; ImageSizeX=8.95; ImageSizeY=5;
ImageFontSize=10; linewidth=1;
map=colormap(lines(4))
h=figure('Units','centimeters','InnerPosition',[10 5 ImageSizeX ImageSizeY],...
           'PaperPosition',[10 5 ImageSizeX ImageSizeY])
h=histogram(-IE1,'Orientation','horizontal','BinWidth',binwidth,'Normalization','Probability','FaceColor',[150 0 0]/255)
hold on
h=histogram(EA1,'Orientation','horizontal','BinWidth',binwidth,'Normalization','Probability','FaceColor',[map(3,:)])
h=histogram(-IE2,'Orientation','horizontal','BinWidth',binwidth,'Normalization','Probability','FaceColor',[57 82 200]/255)
h=histogram(EA2,'Orientation','horizontal','BinWidth',binwidth,'Normalization','Probability','FaceColor',[map(1,:)])
xlabel('DOS','fontsize',ImageFontSize)
ylabel('Energy (eV)','fontsize',ImageFontSize)
legend({'Position 1 -IE','Position 1 EA','Position2 -IE','Position 2 EA'},'Location','best')
legend('boxoff')
h=text(-4,-5.5,'\bf c) \rm')
ylim([-6.5,-3])

%% Plotting Distributions (Distinguishing between positions, No Edge Molecules)
binwidth=0.05;
ImageDPI=600; ImageSizeX=8; ImageSizeY=5;
ImageFontSize=10; linewidth=0.75;
map=colormap(lines(4))
h=figure('Units','centimeters','InnerPosition',[10 5 ImageSizeX ImageSizeY],...
           'PaperPosition',[10 5 ImageSizeX ImageSizeY])
h=histogram(-IE1([1:337 433:768]),'Orientation','horizontal','BinWidth',binwidth,'Normalization','Probability','FaceColor',[150 0 0]/255)
hold on
% h=histogram(EA1([1:337 433:768]),'Orientation','horizontal','BinWidth',binwidth,'Normalization','Probability','FaceColor',[map(3,:)])
h=histogram(-IE2([1:337 433:768]),'Orientation','horizontal','BinWidth',binwidth,'Normalization','Probability','FaceColor',[57 82 200]/255)
% h=histogram(EA2([1:337 433:768]),'Orientation','horizontal','BinWidth',binwidth,'Normalization','Probability','FaceColor',[map(1,:)])
xlabel('DOS','fontsize',ImageFontSize)
ylabel('Energy (eV)','fontsize',ImageFontSize)
% legend({'Position 1 -IE','Position 1  EA','Position2 -IE','Position 2  EA'},'Location','best')
% h=text(-4,-5.5,'\bf c) \rm')
xlim([0,0.22])
ylim([-6.5,-5.25])
h=line([0 0.22],[-5.65 -5.65],'Color',[0 0 0],'LineWidth',0.75,'LineStyle','--','Color',[0 0 0])
legend({'Position 1','Position 2'},'Location','best')
legend('boxoff')
ax=gca;
ax.LineWidth=1;

%% Layer-by-Layer Box and Whisker Plot (Order is Z2 Z1 Y1 Y2)
meanIE=zeros(32,1);   stdIE=zeros(32,1);
meanEA=zeros(32,1);   stdEA=zeros(32,1);
IE3=IE1(385:768,1);   IE4=IE2(385:768,1);
EA3=EA1(385:768,1);   EA4=EA2(385:768,1);
layers=zeros(32,1);

color1=[0,0,0.5]; color2=[0.5,0,0]

for n=1:8
    indices=[(1+(n-1)*48):(n*48)];
    meanIE(3+(n-1)*4)=mean(IE1(indices)); stdIE(3+(n-1)*4)=std(IE1(indices));
    meanIE(2+(n-1)*4)=mean(IE3(indices)); stdIE(2+(n-1)*4)=std(IE3(indices));
    meanIE(4+(n-1)*4)=mean(IE2(indices)); stdIE(4+(n-1)*4)=std(IE2(indices));
    meanIE(1+(n-1)*4)=mean(IE4(indices)); stdIE(1+(n-1)*4)=std(IE4(indices));
    
    meanEA(3+(n-1)*4)=mean(EA1(indices)); stdEA(3+(n-1)*4)=std(EA1(indices));
    meanEA(2+(n-1)*4)=mean(EA3(indices)); stdEA(2+(n-1)*4)=std(EA3(indices));
    meanEA(4+(n-1)*4)=mean(EA2(indices)); stdEA(4+(n-1)*4)=std(EA2(indices));
    meanEA(1+(n-1)*4)=mean(EA4(indices)); stdEA(1+(n-1)*4)=std(EA4(indices));
    
    layers((1+(n-1)*4):(n*4))=[(1+(n-1)*4):(n*4)]';
end

ImageDPI=600; ImageSizeX=8; ImageSizeY=5;
ImageFontSize=10; linewidth=1;
map=colormap(lines(4))
h=figure('Units','centimeters','InnerPosition',[10 5 ImageSizeX ImageSizeY],...
           'PaperPosition',[10 5 ImageSizeX ImageSizeY])

plot(layers,-meanIE,'Marker','sq','MarkerSize',5,'MarkerEdgeColor',color1,'MarkerFaceColor',color1,...
     'LineWidth',1,'Color',color1)
hold on
plot(layers, meanEA,'Marker','sq','MarkerSize',5,'MarkerEdgeColor',color2,'MarkerFaceColor',color2,...
     'LineWidth',1,'Color',color2)
% errorbar(layers,-meanIE,stdIE,'Marker','sq','MarkerSize',5,'MarkerEdgeColor',color1,'MarkerFaceColor',color1,...
%      'LineWidth',1,'Color',color1)
% hold on
% errorbar(layers, meanEA,stdEA,'Marker','sq','MarkerSize',5,'MarkerEdgeColor',color2,'MarkerFaceColor',color2,...
%      'LineWidth',1,'Color',color2)
xlim([0,33])
ylim([-6.5,-3.25])
legend({'-IE','EA'},'Location','best')
legend('boxoff')
ax=gca;
ax.LineWidth=1;
xlabel('Layer','fontsize',ImageFontSize)
ylabel('Energy (eV)','fontsize',ImageFontSize)