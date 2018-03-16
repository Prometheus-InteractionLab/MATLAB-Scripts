%          Name:  Script_Sandia.m
%    Written by:  Telops
% Last modified:  17 April 2015 by Ethan Eagle
%
% Reads IR camera files .hcc and outputs movie and png frames.  adjust
% clims to get the best image quality.  hsvlinearcmp.mat is the best
% colormap and can be found in the 906-165-image\e$ folder.

close all 
clear all

    [filename2,pathname2] = uigetfile('*.mat','Select colormap mat file'); 
    addpath(pathname2)
    load(filename2)


[filename,pathname] = uigetfile('*.hcc','Select IR .hcc files');
filename = cellstr(filename);

%load all data
[D,H,DH] = readIRCam(fullfile(pathname,filename));
% or you can choose the frames that you want to load, below frames 1 to 51.
%[D,H,DH] = readIRCam(fullfile(pathname,filename),'frames',1:51);

%D is a matrix frame-by-pixel
%H is the header (information of each frame)
dayname = '160505'
runname = 'h2'
%time vector
qq=1;
time=zeros(1,size(D,1));
t0=double(H(1).POSIXTime)+double(H(1).SubSecondTime)*1e-7;
for ii=1:size(D,1);
time(qq)=(double(H(ii).POSIXTime)+double(H(ii).SubSecondTime)*1e-7) - t0;
qq=qq+1;
end


% %Remove Badpixel : i will send you the Filterbadpixel code if you want
% later next week.
% RawData=round((D(1,:)-H(1).DataOffset)./H(1).DataLSB);
% Badpix=find(RawData==65534);
% CleanData = filterBadPixels(Badpix, H(1).Width, H(1).Height, D(:,:), 'recurse',true,'method','medoid-sum');

%Make an image
figure,imagesc(formImage(H(1),D(1,:)))%with bad pixel
nrepeat = 3; % increasing this number slows the effective framerate of the finished movie
%adjust the plot limits to make scene visible
clims=[0 2000];

dribbleObj1 = VideoWriter(fullfile(pathname,[dayname runname num2str(clims(2))]),'MPEG-4');
open(dribbleObj1);
cmap = colormap(eval(filename2(1:end-4)));
%video
t0=double(H(1).POSIXTime)+double(H(1).SubSecondTime)*1e-7; 
qq=1;

for jj=1:45;
    img(:,:,jj) = formImage(H(jj),D(jj,:));
    img(img>62500)=0; % removes 'bad' pixels
    imagesc(img(:,:,jj),clims),axis('image'),colorbar %colormap(cmap)
    set(gcf,'color','w');
    drawnow
    if exist(fullfile(pathname,[ dayname runname]))
    else
        mkdir(fullfile(pathname,[ dayname runname]))
    end
    
    saveas(gcf,fullfile(pathname,[ dayname runname], ['IR' num2str(jj) '.png']))
    title(['Frame ',num2str(jj)]) 
    drawnow
    M1l(qq)=getframe(gcf);
    for k = 1:nrepeat
      writeVideo(dribbleObj1,M1l(qq))
    end
    if qq == 1
        saveas(gcf,fullfile(pathname,[dayname runname], ['IR' num2str(jj) '.fig']),'fig')
    end

    qq = qq + 1;
      

end
close(dribbleObj1)
N = size(img);
save(fullfile(pathname,[dayname runname],[dayname runname 'IR.mat']),'img','N','cmap');
