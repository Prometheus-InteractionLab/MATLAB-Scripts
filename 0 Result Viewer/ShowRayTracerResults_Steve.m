function varargout = ShowRayTracerResults(varargin)
% SHOWRAYTRACERRESULTS M-file for ShowRayTracerResults.fig
%      SHOWRAYTRACERRESULTS, by itself, creates a new SHOWRAYTRACERRESULTS or raises the existing
%      singleton*.
%
%      H = SHOWRAYTRACERRESULTS returns the handle to a new SHOWRAYTRACERRESULTS or the handle to
%      the existing singleton*.
%
%      SHOWRAYTRACERRESULTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SHOWRAYTRACERRESULTS.M with the given input arguments.
%
%      SHOWRAYTRACERRESULTS('Property','Value',...) creates a new SHOWRAYTRACERRESULTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ShowRayTracerResults_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ShowRayTracerResults_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ShowRayTracerResults

% Last Modified by GUIDE v2.5 22-Jul-2015 08:56:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ShowRayTracerResults_OpeningFcn, ...
                   'gui_OutputFcn',  @ShowRayTracerResults_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before ShowRayTracerResults is made visible.
function ShowRayTracerResults_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ShowRayTracerResults (see VARARGIN)

% Choose default command line output for ShowRayTracerResults
handles.output = hObject;

%Load and display the piston shape, but hide the handles to it so it isn't
%changed.
set(handles.ShowRayTracerResults_fig,'Name','Loading geometry...')
pause(0.01)

bowl='test';
switch bowl
    case 'GM'
        surfaces=getGMgeometry(handles);
    case 'Ford'
        surfaces=getFordgeometry(handles);
    case 'test'
        surfaces=getTestgeometry(handles);
end

%Now plot all the surfaces
set(handles.ShowRayTracerResults_fig,'CurrentAxes',handles.axes1);
setappdata(handles.ShowRayTracerResults_fig,'surfaces',surfaces);
for s=1:length(surfaces)
   h=plot(surfaces{s}.x,surfaces{s}.y,'k-','Linewidth',2); 
   hold all
   set(h,'HandleVisibility','off')
end
grid on
set(handles.axes1,'XLim',[-60 60]);
set(handles.axes1,'YLim',[-60 72]);
set(handles.ShowRayTracerResults_fig,'Name','Ray Tracer Viewer')

setappdata(handles.ShowRayTracerResults_fig,'fileNum',0);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ShowRayTracerResults wait for user response (see UIRESUME)
% uiwait(handles.ShowRayTracerResults_fig);

% --- Outputs from this function are returned to the command line.
function varargout = ShowRayTracerResults_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Loading, saving%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function load_but_Callback(hObject, eventdata, handles)
    dirIn = getappdata(handles.ShowRayTracerResults_fig,'dirName');
    dirName = uigetdir(dirIn,'Pick the data file directory');
    if isequal(dirName,0)
       return 
    end
    %Save the directory name
    setappdata(handles.ShowRayTracerResults_fig,'dirName',dirName);
    %Look at all .mat files in the directory
    D=dir(fullfile(dirName,'*.mat'));
    %Loop through each filename, save the name and crank angle in
    %cells/matrices.
    filenames=cell(length(D),1);
    % hVector=zeros(length(D),1);
    CAvector=zeros(length(D),1);
    for idx = 1:length(D)
        filename=D(idx).name;
        filenames{idx}=filename;
        filenameString=strrep(filename,'p','.');
        temp1=strfind(filenameString,'h');
    %     temp2=strfind(filenameString,'R');
    %     temp2=temp2(end);
        temp3=strfind(filenameString,'CA_');
    %     hVector(idx)=str2double(filenameString(temp1+1:temp2-2));
        CAvector(idx)=str2double(filenameString(temp3+3:temp1-2));
    end
    %Look for unique crank angles and display them in the listbox
    CAvector=unique(CAvector);
    set(handles.crankAngle_listbox,'String',{CAvector});
    CA=CAvector(get(handles.crankAngle_listbox,'Value'));
    %For the selected crank angle, display all the files
    [filenames4ThisCA nameCell] = getFilesFor_CA(dirName,CA);
    %Save the filenames
    setappdata(handles.ShowRayTracerResults_fig,'filenames',...
        filenames4ThisCA);
    %Set the filenames for the currently selected crank angle in the file 
    %listbox
    set(handles.files_listbox,'String',nameCell);

    % hVector=unique(hVector);
    % set(handles.crankAngle_listbox,'String',{hVector});
    % h=hVector(get(handles.crankAngle_listbox,'Value'));
    % [filenames radiiCell] = getRadiiFor_h(dirName,h);
    % set(handles.files_listbox,'String',radiiCell);
    setappdata(handles.ShowRayTracerResults_fig,'fileNum',0);
    redraw(handles)

function [filenames nameCell] = getFilesFor_CA(dirName,CA)
    %Look at .mat files in the directory for the selected h value, see 
    %which radii are available
    CAString=sprintf('%06.2f',CA);
    CAString=strrep(CAString,'.','p');
    D=dir(fullfile(dirName,['*',CAString,'*.mat']));
    filenames=cell(length(D),1);
    nameCell=cell(length(D),1);
%     radii=zeros(length(D),1);
    for idx=1:length(D)
        filename=D(idx).name;
        filenames{idx}=filename;
        filenameString=strrep(filename,'p','.');
        temp1=strfind(filenameString,'_h')+1;
        nameCell{idx}=filenameString(temp1:end-4);
    end
 
function surfaces=getFordgeometry(handles)
    Fordfilename=['C:\Users\sbusch\Documents\Data\Software\RayTracer\',...
        'Bowl form\FordBowl\RZcartesianCord.xlsx'];
    GMfilename=['C:\Users\sbusch\Documents\Data\Software\RayTracer\',...
            'Bowl form\FullBowlProfile_FP48_B2.xlsx'];
    dt=0.0005; %Spline resolution
        %Get bowl coordinates
        disp('Reading bowl coordinates...');
        pause(0.01)
        [numeric, txt, raw] = xlsread(Fordfilename,'RZcartesianCord',...
            '$H$3:$I$1846');
        x=numeric(:,1);
        x=[-flipud(x);x;x(end);x(end)]; %mirror to make the whole profile
        y=numeric(:,2);
        y=[flipud(y);y;y(end);y(end)]; %mirror to make the whole profile
     %Create bowl surface using coordinates from file, 3rd order B-splines,
    %step size of 0.0001, index of refraction on one side 1, other side 1.4
        surfaces{1}=createSurface(x,y,3,dt,'refracting',1,1.46);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Get coordinates for the side of the piston
    set(handles.ShowRayTracerResults_fig,'Name',['Reading side ',...
        'coordinates...'])
    pause(0.01)
    [numeric, txt, raw] = xlsread(GMfilename,'Full Piston','$C$19:$D$120');
    x=numeric(:,1);
    y=numeric(:,2);
    surfaces{end+1}=createSurface(x,y,3,2*dt,'absorbing',1.46,1);
    % surfaces{end+1}=createSurface(x+0.05,y,3,0.0001,'absorbing',1.46,1);
    surfaces{end+1}=createSurface(-x,y,3,2*dt,'absorbing',1.46,1);
    % surfaces{end+1}=createSurface(-x-0.05,y,3,0.0001,'absorbing',1.46,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Get bottom surface of piston
    set(handles.ShowRayTracerResults_fig,'Name',['Reading bottom ',...
        'coordinates...'])
    pause(0.01)
    [numeric, txt, raw] = xlsread(GMfilename,'Full Piston','$E$19:$F$280');
    x=numeric(:,1);
    x=[x;-flipud(x)]; %mirror to make the whole profile
    x=[x;x(end)];
    y=numeric(:,2);
    y=[y;flipud(y)]; %mirror to make the whole profile
    y=[y;y(end)];
    surfaces{end+1}=createSurface(x,y,3,dt,'refracting',1.46,1);
    
function surfaces=getGMgeometry(handles)
    GMfilename=['C:\Users\sbusch\Documents\Data\Software\RayTracer\',...
            'Bowl form\FullBowlProfile_FP48_B2.xlsx'];
    [numeric, txt, raw] = xlsread(GMfilename,'Full Piston','$A$19:$B$407');
    x=numeric(:,1);
    x=[x;-flipud(x)]; %mirror to make the whole profile
    y=numeric(:,2);
    y=[y;flipud(y)]; %mirror to make the whole profile
    %Create bowl surface using coordinates from file, 3rd order B-splines,
    %step size of 0.0001, index of refraction on one side 1, other side 1.4
    surfaces{1}=createSurface(x,y,3,0.0001,'refracting',1,1.46);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Get coordinates for the side of the piston
    set(handles.ShowRayTracerResults_fig,'Name',['Reading side ',...
        'coordinates...'])
    pause(0.01)
   [numeric, txt, raw] = xlsread(GMfilename,'Full Piston','$C$19:$D$120');
    x=numeric(:,1);
    y=numeric(:,2);
    surfaces{end+1}=createSurface(x,y,3,0.0001,'absorbing',1.46,1);
    surfaces{end+1}=createSurface(-x,y,3,0.0001,'absorbing',1.46,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Get bottom surface of piston
    set(handles.ShowRayTracerResults_fig,...
        'Name','Reading bottom coordinates...')
    pause(0.01)
    [numeric, txt, raw] = xlsread(GMfilename,'Full Piston','$E$19:$F$280');
    x=numeric(:,1);
    x=[x;-flipud(x)]; %mirror to make the whole profile
    y=numeric(:,2);
    y=[y;flipud(y)]; %mirror to make the whole profile
    surfaces{end+1}=createSurface(x,y,3,0.0001,'refracting',1.46,1);
    
    function surfaces=getTestgeometry(handles)
        testfilename=['\\906-165-press\e$\Data\matlabs\0 Ray Tracer\test1.xlsx'];
        dt=0.001; %Spline resolution
        %Get bowl coordinates
        disp('Reading bowl coordinates...');
        pause(0.01)
        [numeric, txt, raw] = xlsread(testfilename,'Sheet1',...
            '$A$14:$B$102'); %#ok<*NASGU>
        x=numeric(:,1);
        x=[x;-flipud(x)]; %mirror to make the whole profile
        y=numeric(:,2);
        y=[y;flipud(y)]; %mirror to make the whole profile
        %Create bowl surface using coordinates from file, 3rd order B-splines,
        %step size of 0.0001, index of refraction on one side 1, other side 1.4
        surfaces{1}=createSurface(x,y,3,dt,'refracting',1,1.46);
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Get coordinates for the side of the piston
    set(handles.ShowRayTracerResults_fig,'Name',['Reading side ',...
        'coordinates...'])
    pause(0.01)
   [numeric, txt, raw] = xlsread(testfilename,'Sheet1','$C$14:$D$86');
    x=numeric(:,1);
    y=numeric(:,2);
    surfaces{end+1}=createSurface(x,y,3,0.0001,'absorbing',1.46,1);
    surfaces{end+1}=createSurface(-x,y,3,0.0001,'absorbing',1.46,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Get bottom surface of piston
    set(handles.ShowRayTracerResults_fig,...
        'Name','Reading bottom coordinates...')
    pause(0.01)
    [numeric, txt, raw] = xlsread(testfilename,'Sheet1','$E$14:$F$83');
    x=numeric(:,1);
    x=[x;-flipud(x)]; %mirror to make the whole profile
    y=numeric(:,2);
    y=[y;flipud(y)]; %mirror to make the whole profile
    surfaces{end+1}=createSurface(x,y,3,0.0001,'refracting',1.46,1);
    
function saveaxes_but_Callback(hObject, eventdata, handles)
    figure
    newax=axes;
%     bowl='Ford';
%     switch bowl
%         case 'GM'
%             surfaces=getGMgeometry(handles);
%         case 'Ford'
%             surfaces=getFordgeometry(handles);
%     end

    %Now plot all the surfaces
%     set(handles.ShowRayTracerResults_fig,'CurrentAxes',handles.axes1);
%     setappdata(handles.ShowRayTracerResults_fig,'surfaces',surfaces);
%     for s=1:length(surfaces)
%        h=plot(surfaces{s}.x,surfaces{s}.y,'k-','Linewidth',2); 
%        hold all
%        set(h,'HandleVisibility','off')
%     end
    grid on
    
    a=get(handles.axes1,'Children');
%     colormap(getHSVlinearMap());
    colormap(flipud(gray(256)));
    newHandle=copyobj(a,newax);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Determine the mapping function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function getMappingfcn_but_Callback(hObject, eventdata, handles)
    %Loop through each file that the user has selected. Bin the image as is
    %done in the redraw function. Sum over the columns and 
    dirName=getappdata(handles.ShowRayTracerResults_fig,'dirName');
    %These are the filenames for the selected crank angle.
    filenames=getappdata(handles.ShowRayTracerResults_fig,'filenames');
    dx=str2double(get(handles.dx_ctl,'String')); %mm
    dy=str2double(get(handles.dy_ctl,'String'));
    %For each file, save the R and potential "R" coordinates
    R=zeros(length(filenames),1);
    Rmap=cell(length(filenames),1); %Cell because R may map to several "R"s
    RhighestVector=zeros(length(filenames),1);
    for idx=1:length(filenames)
        set(handles.files_listbox,'Value',idx)
        filename=filenames{idx};
        temp1=strfind(filename,'R');
        temp1=temp1(end);
        radiusString=filename(temp1+1:end-4);
        radiusString=strrep(radiusString,'p','.');
        R(idx)=str2double(radiusString);
        
        %Load data
        set(handles.ShowRayTracerResults_fig,'Name',['Loading ',...
            filename, '...'])
        drawnow
        S=load(fullfile(dirName,filename));
        X=S.X;
        Y=S.Y;
        clear S
        yMin=-70;
        yMax=20;
        X=X(Y < yMax & Y > yMin); %limit range of y data
        Y=Y(Y < yMax & Y > yMin);
        
        %Bin the data
        set(handles.ShowRayTracerResults_fig,'Name',...
            sprintf('Binning %0.0f points...',numel(X)));
        drawnow
        try
            [im,xVector,yVector]=binData(X,Y,dx,dy,'fast');%#ok<NASGU>
        catch
            continue
        end
        if isempty(im)
            %it didn't work
           continue 
        end
        %sum each column and display (hP = horizontal Profile)
        hP=sum(im,1);
        %Get local maxima of the horizontal profile
        limit_mm=0.5; %Peaks must be separated by at least this many mm.
        separationLimit=round(limit_mm/dx); %separation limit in pixels
        threshold=0.1; %Look for peaks at least this much of the global max
        [peaks, peakIndices]=findMaxima(hP,threshold,separationLimit);

        %Now look around each peak and get the integral over the fwhm, as
        %well as the centroid
%         peakAreas=zeros(length(peaks),1);
%         peakCentroids=zeros(length(peaks),1);
%         for k=1:length(peaks)
%             peakHeight=peaks(k);
%             iLow=find(hP(1:peakIndices(k))<peakHeight/2,1,'last');
%             iHigh=find(hP(peakIndices(k):end)<peakHeight/2,1,'first')+...
%                 peakIndices(k)-1;
%             %Integrate over the FWHM
%             centSum=iLow*hP(iLow);
%             centSumDenom=hP(iLow);
%             for j=iLow+1:iHigh
%                peakAreas(k) = peakAreas(k)+( (hP(j-1)+hP(j))/2 ) * dx;
%                centSum = centSum + j*hP(j);
%                centSumDenom = centSumDenom + hP(j);
%             end
%             peakCentroids(k) = round(centSum/centSumDenom);
%         end
%         %Which peak has the higher area?
%         [maxArea, iLargestArea] = max(peakAreas);
%         %Where are the peaks located?
%         
% %         [maxRow ImaxRow]=max(im,[],1);
% %         [maxVal iColMax]=max(maxRow);
% %         [maxCol ImaxCol]=max(im,[],2);
% %         [maxVal iRowMax]=max(maxCol);
%         
        
%         
%         iAllPeaks = peakCentroids;
%         iLargestArea = peakCentroids(iLargestArea);
%         try
%             xAllCentroids=xVector(iAllPeaks);
%         catch
%         end
        xAllPeaks=xVector(peakIndices);
%         [temp iHighestPeak]=max(hP);
%         xHighestPeak=xVector(iHighestPeak);
%         xLargestArea=xVector(iLargestArea);
        
        Rmap{idx}=xAllPeaks;
        
        [~,sortIndex] = sort(im(:),'descend');
        maxIndices = sortIndex(1:5); %five brightest pixels
        [I,J]=ind2sub(size(im),maxIndices);
        I=median(I);
        J=median(J);
        RhighestVector(idx)=xVector(J);
%         RhighestVector(idx)=xVector(iColMax);
    end
    assignin('base','R',R);
    assignin('base','Rmap',Rmap);
    assignin('base','RhighestVector',RhighestVector);
    f2=figure(2);
%     for r=1:numel(R)
%         plot(RhighestVector(r),R(r),'Color',[0 0 0.8],'Marker','o',...
%             'MarkerSize',12,'LineStyle','none')
%         for n=1:length(Rmap{r}) 
%            plot(Rmap{r}(n),R(r),'Color',[0.8 0 0],'Marker','+',...
%                'LineStyle','none')
%            hold all
%         end
%     end
    plot(RhighestVector,R,'Color',[0 0 0.8],'Marker','o',...
            'MarkerSize',12,'LineStyle','none')
    xlabel('R(distorted)')
    ylabel('R(undistorted)')
    %Allow the user to remove data points as desired. Make a button so that
    %when the user clicks it, the data are saved
    uicontrol('Style','pushbutton','String','Save','Position',...
        [84 340 100 40],'Callback',{@(hObject,...
        eventdata)ButtonCallback(hObject,eventdata,guidata(hObject))});
    
function smartMapping_but_Callback(hObject, eventdata, handles)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Bin data according to dx and dy specified by the user
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dx=str2double(get(handles.dx_ctl,'String')); %mm
    dy=str2double(get(handles.dy_ctl,'String'));
    yMin=-80;
    yMax=20;
    numObjects=2; %How many centroids to look for
    numBrightest=5; %How many of the brightest pixels are considered?
    
    dirName=getappdata(handles.ShowRayTracerResults_fig,'dirName');
    %These are the filenames for the selected crank angle:
    filenames=getappdata(handles.ShowRayTracerResults_fig,'filenames');
    R=zeros(length(filenames),1);   %initialize unwarped R vector
    Rwarped=zeros(length(filenames),numObjects+1); %multiple warped radii
    
    for idx=1:length(filenames)
        set(handles.files_listbox,'Value',idx)
        filename=filenames{idx};
        temp1=strfind(filename,'R');
        temp1=temp1(end);
        radiusString=filename(temp1+1:end-4);
        radiusString=strrep(radiusString,'p','.');
        R(idx)=str2double(radiusString);
        %Load data
        set(handles.ShowRayTracerResults_fig,'Name',['Loading ',...
            filename, '...'])
        pause(0.001)
        S=load(fullfile(dirName,filename));
        X=S.X;
        Y=S.Y;
        clear S
    
        set(handles.ShowRayTracerResults_fig,'Name',...
            sprintf('Binning %0.0f points...',numel(X)));
        pause(0.01)
        try
            [im,xVector,yVector] = ...
                binnedImageFromIntersecs(X,Y,dx,dy,yMin,yMax);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Find the largest objects based on morphological ops
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [bw Areas xVec yVec] = getLargestObjects(im,numObjects,...
                xVector,yVector);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Find the median of the numBrightest brightest points
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [xVal yVal] = getBrightestPix(im,numBrightest,xVector,yVector);
            Rwarped(idx,:)=[xVal,xVec];
        catch %#ok<CTCH>
        end
    end
    figure(1)
    plot(Rwarped,R,'Color',[0 0 0.8],'Marker','o','MarkerSize',12,...
        'LineStyle','none')
    xlabel('R(distorted)')
    ylabel('R(undistorted)')
    %Allow the user to remove data points as desired. Make a button so that
    %when the user clicks it, the data are saved
    uicontrol('Style','pushbutton','String','Save','Position',...
        [84 340 100 40],'Callback',{@(hObject,...
        eventdata)ButtonCallback(hObject,eventdata,guidata(hObject))});

function rayBundles_but_Callback(hObject, eventdata, handles)
    dirName=getappdata(handles.ShowRayTracerResults_fig,'dirName');
    %These are the filenames for the selected crank angle:
    filenames=getappdata(handles.ShowRayTracerResults_fig,'filenames');
    dx=str2double(get(handles.dx_ctl,'String')); %mm
    dy=str2double(get(handles.dy_ctl,'String'));
    
    nBundlesMax=3;
    RrealSpace=zeros(length(filenames),1);
    RvirtualImage=zeros(length(filenames),nBundlesMax);
    xVirtualImage=zeros(length(filenames),nBundlesMax);
    yVirtualImage=zeros(length(filenames),nBundlesMax);
    xRayOrigin=zeros(length(filenames),1);
    yRayOrigin=zeros(length(filenames),1);
    
    idx0=get(handles.files_listbox,'Value');
    
    fileNum=getappdata(handles.ShowRayTracerResults_fig,'fileNum');
    setappdata(handles.ShowRayTracerResults_fig,'fileNum',fileNum+1);
    for idx=1:length(filenames)
        
        set(handles.files_listbox,'Value',idx)
        filename=filenames{idx};
        temp1=strfind(filename,'R');
        temp1=temp1(end);
        radiusString=filename(temp1+1:end-4);
        radiusString=strrep(radiusString,'p','.');
        RrealSpace(idx)=str2double(radiusString);
        %Load data
        set(handles.ShowRayTracerResults_fig,'Name',['Loading ',...
            filename, '...'])
        drawnow
        S=load(fullfile(dirName,filename));
        
        finalRays=S.finalrays; %grab all of the final rays
        intersections=S.intersections;
        numUsefulRays=numel(S.usefulrays); %how many of them are useful?
        usefulRayIndices=zeros(numUsefulRays,1);%Initialize array of indices
        clear S
        %loop though final rays and find the ones with positive power.
        %Store their indices in usefulRayIndices
        cnt=1;
        for ray=1:length(finalRays)
           if finalRays{ray}.power > 0
              usefulRayIndices(cnt) = ray;
              cnt=cnt+1;
           end
        end
        %Compute the spacing between rays. If the value is 1, then the rays
        %are adjacent
        raySpacing=[1;diff(usefulRayIndices)];
        %Find the locations where the ray spacing is greater than 1
        breaks=find(raySpacing > 1);
        if isempty(breaks)
%             breaks=1;
        end
        %Get the ray origin from any ray's intersections
        temp = intersections{1};
        xRayOrigin(idx) = temp(1,1);
        yRayOrigin(idx) = temp(2,1);
        
        numBundles=numel(breaks)+1;
        numRays=numel(breaks);
        if numBundles == 1
            bundles(1).indices = usefulRayIndices;

        elseif numBundles == 2
            bundles(1).indices=...
                usefulRayIndices(1:breaks(1)-1);
                
            bundles(2).indices =...
                usefulRayIndices(breaks(1):end);
            
            numRays(1)=numel(bundles(1).indices); 
            numRays(2)=numel(bundles(2).indices); 
            
        elseif numBundles == 3
            bundles(1).indices=...
                usefulRayIndices(1:breaks(1)-1);
            bundles(2).indices =...
                usefulRayIndices(breaks(1):breaks(2)-1);
            bundles(3).indices =...
                usefulRayIndices(breaks(2):end);
            
            numRays(1)=numel(bundles(1).indices); 
            numRays(2)=numel(bundles(2).indices); 
            numRays(3)=numel(bundles(3).indices);
            
        elseif numBundles == 4
            bundles(1).indices=...
                usefulRayIndices(1:breaks(1)-1);
            bundles(2).indices =...
                usefulRayIndices(breaks(1):breaks(2)-1);
            bundles(3).indices =...
                usefulRayIndices(breaks(2):breaks(3)-1);
            bundles(4).indices =...
                usefulRayIndices(breaks(3):end);
            
            numRays(1)=numel(bundles(1).indices); 
            numRays(2)=numel(bundles(2).indices); 
            numRays(3)=numel(bundles(3).indices);
            numRays(4)=numel(bundles(4).indices);
        else
            beep
            errordlg('Too many bundles! Add some code...','','modal')
            keyboard
        end
        %Limit locations of intersections
        yMin=-60;
        yMax=10;
        %Go through each bundle and compute its intersections

        for b=1:numBundles
            [X,Y] = findCrossings(finalRays(bundles(b).indices));
            
            X=X(Y < yMax & Y > yMin); %limit range of y data
            Y=Y(Y < yMax & Y > yMin);
            if isempty(X) || isempty(Y)
               numRays(b)=0;
               continue; 
            end
            
            bundles(b).X=X; %#ok<AGROW>
            bundles(b).Y=Y; %#ok<AGROW>
            [im,xVector,yVector] = binData(X,Y,dx,dy,'fast');
            bundles(b).xVector=xVector; %#ok<AGROW>
            bundles(b).yVector=yVector; %#ok<AGROW>
            if isempty(im) || isempty(xVector)
                bundles(b).Xmean=[]; %#ok<AGROW>
                bundles(b).Ymean=[]; %#ok<AGROW>
            else
                [~,sortIndex] = sort(im(:),'descend');
                if numel(sortIndex) < 5
                    maxIndices = sortIndex(1:end);
                else
                    maxIndices = sortIndex(1:5); %five brightest pixels
                end
                [I,J]=ind2sub(size(im),maxIndices);
                I=median(I);
                J=median(J);
                bundles(b).Xmean=xVector(J); %#ok<AGROW>
                bundles(b).Ymean=yVector(round(I)); %#ok<AGROW>

            end
        end
        %What is the largest bundle (LB)?
        [~,iLB]=max(numRays);        
  
        try
            RvirtualImage(idx,1)=bundles(iLB).Xmean;
            xVirtualImage(idx,1)=bundles(iLB).Xmean;
            yVirtualImage(idx,1)=bundles(iLB).Ymean;
            altNumRays=numRays;
            altNumRays(altNumRays==max(numRays))=0;
            [~,i2ndLB]=max(altNumRays); 
            RvirtualImage(idx,2)=bundles(i2ndLB).Xmean;
            xVirtualImage(idx,2)=bundles(i2ndLB).Xmean;
            yVirtualImage(idx,2)=bundles(i2ndLB).Ymean;
            altNumRays(altNumRays==max(altNumRays))=0;
            [~,i3rdLB]=max(altNumRays); 
            RvirtualImage(idx,3)=bundles(i3rdLB).Xmean;
            xVirtualImage(idx,3)=bundles(i3rdLB).Xmean;
            yVirtualImage(idx,3)=bundles(i3rdLB).Ymean;
%             [~,iFourthLargestBundle]=max(altNumRays); 
%             RvirtualImage(idx,4)=bundles(iFourthLargestBundle).Xmean;
        catch
            %No point to calculate
        end
        
    end
    set(handles.files_listbox,'Value',idx0);
    f2=figure(3);
    plot(RvirtualImage(:,1),RrealSpace,'Color',[0.8 0 0],'Marker','o',...
            'MarkerSize',12,'LineStyle','none')
    hold all
    plot(RvirtualImage(:,2),RrealSpace,'Color',[0 0.6 0],'Marker','o',...
        'MarkerSize',10,'LineStyle','none')
    plot(RvirtualImage(:,3),RrealSpace,'Color',[0 0 0.8],'Marker','+',...
        'MarkerSize',6,'LineStyle','none')
%     plot(RvirtualImage(:,4),RrealSpace,'Color',[0 0 0],'Marker','+',...
%         'MarkerSize',8,'LineStyle','none')
    xlabel('RvirtualImage')
    ylabel('RrealSpace')
    %Allow the user to remove data points as desired. Make a button so that
    %when the user clicks it, the data are saved
    uicontrol('Style','pushbutton','String','Save','Position',...
        [84 390 100 20],'Callback',{@(hObject,...
        eventdata)RayBundleButtonCallback(hObject,eventdata,...
        handles)});
    %Create plot to show relationship between ray origins and virtual image
    %locations - put the head at zero and move the piston
    f4=figure(4);
    set(f4,'Units','Pixels','Position',[400 400 660 660]);
    ax3=axes('Parent',f4,'Units','Pixels','Position',[50 50 600 600]);
    han = get(ax3,'Children');
    delete(han);
    %Crank angle
    CAvector=get(handles.crankAngle_listbox,'String');
    CA=str2double(CAvector{get(handles.crankAngle_listbox,'Value')});
    %Piston position
    ImportPistonPosition('Piston Position Pin Offset.TXT'); %Kan's function
    PistonP=zeros(2*size(data,1),1);
    PistonP(1801:(1800+size(data,1)),1)=data(:,2);
    PistonP(1:1800,1)=data(1801:end,2);
    PistonP((1801+size(data,1)):end,1)=data(1:1800,2);
    SquishHeight=0.028*25.4;
    %PistonP is s, the distance between the top of the piston and the firedeck
    PistonP=PistonP+SquishHeight; %Function of crank angle

    CAD=(-360:0.1:359.9)'; %Crank angle vector with 0.1° resolution
    
    h=interp1(CAD,PistonP,CA,'linear'); 
    %Plot the piston surfaces
    surfaces = getappdata(handles.ShowRayTracerResults_fig,'surfaces');
    for s=1:length(surfaces)
       plot(surfaces{s}.x,surfaces{s}.y-h,'k-','Linewidth',2); 
       hold all
    end
    grid on
    set(ax3,'XLim',[-45 45]);
    set(ax3,'YLim',[-88 2]);
    %Plot virtual images
    color1 = 0.4*[1 1 1];
    color2 = [0.8 0.2 0.2];
    color3 = [0.15 0.2 0.8];
    plot(xVirtualImage(:,1),yVirtualImage(:,1)-h,'Color',color1,...
        'Marker','^','MarkerFaceColor',color1,'LineWidth',1,...
        'LineStyle','none')
    plot(xVirtualImage(:,2),yVirtualImage(:,2)-h,'Color',color2,...
        'Marker','+','MarkerFaceColor',color2,'LineWidth',1,...
        'LineStyle','none')
    plot(xVirtualImage(:,3),yVirtualImage(:,3)-h,'Color',color3,...
        'Marker','v','MarkerFaceColor',color3,'LineWidth',1,...
        'LineStyle','none')
    %Plot ray origins
    plot(xRayOrigin,yRayOrigin-h,'Color',[0 0.4 0],'LineWidth',1,...
        'Marker', 'o','MarkerFaceColor','none','MarkerSize',4);
    %Plot Firedeck
    plot([-41;41],[0;0],'Color',0.35*[1 1 1],'LineWidth',2)
    cnt=0;
    for x = -41:0.5:41
        plot([x, x+1],[0 1],'Color',0.35*[1 1 1],'LineWidth',2)
        cnt = cnt+1;
    end
    
    text(-40,-2,sprintf('%1.1f CAD ATDC',CA),'FontName','Arial',...
        'FontWeight','Bold','FontSize',14,'BackgroundColor',[1 1 1]);
    xlabel('x / mm','FontName','Arial','FontWeight','Bold','FontSize',14);
    ylabel('z / mm','FontName','Arial','FontWeight','Bold','FontSize',14);
    
return    

function RayBundleButtonCallback(hObject, eventdata, handles) 
    hFig=get(hObject,'Parent');
    hAx=get(hFig,'CurrentAxes');
    plot=get(hAx,'Children');
    Rdistorted=get(plot,'XData'); %This is a 3x1 cell, and each element
    %of the cell contains a vector of virtual image radii. 
    %If any of the data points have been deleted, the x-value becomes zero
    R=get(plot,'YData'); %This is a 3x1 cell, and each element of the cell
    %contains a vector of realspace radii. If a data pointa has been
    %deleted, the R (y) value has been changed to NaN.
    
    %loop through the radii. For each realspace radius (R), check to see 
    %how many unique virtual image radii there are (use the one associated
    %with the largest ray bundle at that point). If the realspace radius
    %has been changed to NaN, use another element. If all three are NaN,
    %don't save that point.
    rRealspaceVec=[];
    rVirtualImageVec=[];
    
    for b = 1:numel(R{1})
        validPt=false;
        cnt=0;
        while(~validPt && cnt < 3)
            cnt=cnt+1;
            Rtemp=R{cnt}(b);
            validPt = ~isnan(Rtemp);
            if validPt
                ptNum=cnt;
            else
                ptNum=[];
            end
        end
        if validPt
            %Add the R value to the realspace vector
            rRealspaceVec(end+1,1)=Rtemp;
            %Look at virtual image radii. Are they the same?
            rVirtualImageVec(end+1,1) = Rdistorted{ptNum}(b);
        end
    end
    
    R = rRealspaceVec;
    Rdistorted = rVirtualImageVec;
    
    
%     %only take values, not NaN
%     idx=~isnan(R);
%     Rdistorted=Rdistorted(idx);
%     R=R(idx);
    
    %Name the file and save it
    fileNum=getappdata(handles.ShowRayTracerResults_fig,'fileNum');
    handles=guidata(findobj('Tag','ShowRayTracerResults_fig'));
    CAcell=get(handles.crankAngle_listbox,'String');
    CA=CAcell{get(handles.crankAngle_listbox,'Value')};
    dx=str2double(get(handles.dx_ctl,'String')); %mm
    dy=str2double(get(handles.dy_ctl,'String'));
    filename2save=sprintf('%02.0f_CA %s_dx%3.2f_dy%3.2f',fileNum,CA,dx,dy);
    filename2save=strrep(filename2save,'.','p');
    filename2save=[filename2save,'.mat'];
    
    dirName=getappdata(handles.ShowRayTracerResults_fig,'dirName');
    Rmap=[];
    save(fullfile(dirName,filename2save),'Rdistorted','R','Rmap')

function ButtonCallback(hObject, eventdata, handles) 
    hFig=get(hObject,'Parent');
    hAx=get(hFig,'CurrentAxes');
    plot=get(hAx,'Children');
    Rdistorted=get(plot,'XData');
    R=get(plot,'YData');
    Rmap=[];
    %only take values, not NaN
    idx=~isnan(R);
    Rdistorted=Rdistorted(idx);
    R=R(idx);
    handles=guidata(findobj('Tag','ShowRayTracerResults_fig'));
    CAcell=get(handles.crankAngle_listbox,'String');
    CA=CAcell{get(handles.crankAngle_listbox,'Value')};
    dx=str2double(get(handles.dx_ctl,'String')); %mm
    dy=str2double(get(handles.dy_ctl,'String'));
    filename2save=sprintf('CA %s_dx%3.2f_dy%3.2f',CA,dx,dy);
    filename2save=strrep(filename2save,'.','p');
    filename2save=[filename2save,'.mat'];
    dirName=getappdata(handles.ShowRayTracerResults_fig,'dirName');
    save(fullfile(dirName,filename2save),'Rdistorted','R','Rmap')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Changing what to display, plotting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function crankAngle_listbox_Callback(hObject, eventdata, handles)
    dirName=getappdata(handles.ShowRayTracerResults_fig,'dirName');
    %See which crank angle is selected
    CAvector=get(handles.crankAngle_listbox,'String');
    CA=str2double(CAvector{get(handles.crankAngle_listbox,'Value')});
    %For the selected crank angle, display all the files
    [filenames4ThisCA nameCell] = getFilesFor_CA(dirName,CA);
    %Save the filenames
    setappdata(handles.ShowRayTracerResults_fig,'filenames',...
        filenames4ThisCA);
    %Set the filenames for the currently selected crank angle in the file 
    %listbox
    set(handles.files_listbox,'String',nameCell);
    
    redraw(handles)

function files_listbox_Callback(hObject, eventdata, handles)
    redraw(handles)

function redraw(handles)
    %Based on piston position and radial position, build up the filename to
    %load and display
    filenames=getappdata(handles.ShowRayTracerResults_fig,'filenames');
    filename=filenames{get(handles.files_listbox,'Value')};
    dirName=getappdata(handles.ShowRayTracerResults_fig,'dirName');
    %Load data
    set(handles.ShowRayTracerResults_fig,'Name','Loading file...')
    pause(0.01)
    S=load(fullfile(dirName,filename));
    X=S.X;
    Y=S.Y;
    usefulrays=S.usefulrays;
    startPoint=S.intersections{1,:};
    intersections=S.intersections;
    finalrays=S.finalrays;
    clear S
    set(handles.ShowRayTracerResults_fig,'Name','Plotting...')
    pause(0.01)
    butTag=get(get(handles.displayType_panel,'SelectedObject'),'Tag');
    switch butTag
        case 'intersec_but'
            set(handles.ShowRayTracerResults_fig,'Renderer','Painters')
            %Clear children from axes1
            h=get(handles.axes1,'Children');
            delete(h);
            plot(X,Y,'Color',[1 0 0],'LineStyle','none','Marker','o',...
                'MarkerFaceColor',[1 0 0],'MarkerEdgeColor','none',...
                'MarkerSize',1);
            hold on
            x0=zeros(length(usefulrays),1);
            y0=zeros(length(usefulrays),1);
            for idx=1:length(usefulrays)
               ray=usefulrays{idx}; 
               x0(idx)=ray.start(1);
               y0(idx)=ray.start(2);
            end
            plot(x0,y0,'Color',[0 0 1],'LineStyle','none','Marker','o',...
                'MarkerFaceColor','none','MarkerEdgeColor',[0 0 1],...
                'MarkerSize',3);
            plot(startPoint(1),startPoint(2),'Color',[0.7 0.7 0.7],...
                'Marker','+','LineWidth',2,'MarkerEdgeColor',...
                [0.6 0.6 0.6],'MarkerSize',20);
    case 'intensity_but'
        set(handles.ShowRayTracerResults_fig,'Renderer','Painters')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Bin data according to dx and dy specified by the user
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dx=str2double(get(handles.dx_ctl,'String')); %mm
        dy=str2double(get(handles.dy_ctl,'String'));
        yMin=-80;
        yMax=20;
        X=X(Y < yMax & Y > yMin); %limit range of y data
        Y=Y(Y < yMax & Y > yMin);
        
        if isempty(X) || isempty(Y)
            set(handles.ShowRayTracerResults_fig,'Name',...
            'No intersections to plot!');
            pause(0.5)
            set(handles.ShowRayTracerResults_fig,...
                'Name','Ray Tracer Viewer')
           return 
        end
        set(handles.ShowRayTracerResults_fig,'Name',...
            sprintf('Binning %0.0f points...',numel(X)));
        drawnow
        [im,xVector,yVector] = binData(X,Y,dx,dy, 'fast'); 
        if isempty(im)
            return
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Scale axes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        xlim=get(handles.axes1,'XLim');
        ylim=get(handles.axes1,'YLim');
        h=get(handles.axes1,'Children');
        delete(h);
        %Show intensity map
        imagesc(im,'XData',[min(xVector),max(xVector)],'YData',...
            [min(yVector),max(yVector)],'AlphaData',1)
%         colormap(GetLaVisionMap);
%         colormap(getHSVlinearMap());
        colormap(flipud(gray(256)));
        hold all
        %Plot the piston contour
        surfaces=getappdata(handles.ShowRayTracerResults_fig,'surfaces');
        for s=1:length(surfaces)
           h=plot(surfaces{s}.x,surfaces{s}.y,'Color',[0 0 0],...
               'Linewidth',2); 
        end
        %Plot the start point
        plot(startPoint(1),startPoint(2),'Color',[0.7 0.7 0.7],...
            'Marker','x','LineWidth',1,'MarkerEdgeColor',...
            [0 0 0],'MarkerSize',16);
        grid on
        set(handles.axes1,'XLim',xlim);
        set(handles.axes1,'YLim',ylim);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Find peaks and locate the mapped radii
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        set(handles.ShowRayTracerResults_fig,'Name',...
            'Calculating maxima...');
        drawnow

        usefulRayIndices=zeros(numel(usefulrays),1);%Initialize
        cnt=1;
        for ray=1:length(finalrays)
           if finalrays{ray}.power > 0
              usefulRayIndices(cnt) = ray;
              cnt=cnt+1;
           end
        end
        %Compute the spacing between rays. If the value is 1, then the rays
        %are adjacent
        raySpacing=[1;diff(usefulRayIndices)];
        %Find the locations where the ray spacing is greater than 1
        breaks=find(raySpacing > 1);
        if isempty(breaks)
            breaks=1;
        end
        %Get the ray origin from any ray's intersections
        temp = intersections{1};
        xRayOrigin = temp(1,1);
        yRayOrigin = temp(2,1);

        numBundles=numel(breaks);
        numRays=numel(breaks);
        if numBundles ==1
            bundles(1).indices = usefulRayIndices;

        else
            for b = 1:numBundles
                if b == 1
                    bundles(b).indices=...
                        usefulRayIndices(1:breaks(b)-1);%#ok<AGROW>
                else
                    bundles(b).indices=...
                        usefulRayIndices(breaks(b-1):breaks(b)-1);  %#ok<AGROW>
                end
                numRays(b)=numel(bundles(b).indices); 
            end
        end
        %Limit locations of intersections
        yMin=-60;
        yMax=10;
        %Go through each bundle and compute its intersections
        for b=1:numBundles
            [X,Y] = findCrossings(finalrays(bundles(b).indices));

            X=X(Y < yMax & Y > yMin); %limit range of y data
            Y=Y(Y < yMax & Y > yMin);
            if isempty(X) || isempty(Y)
               numRays(b)=0;
               continue; 
            end

            bundles(b).X=X; %#ok<AGROW>
            bundles(b).Y=Y; %#ok<AGROW> 
            bundles(b).xVector=xVector; %#ok<AGROW>
            bundles(b).yVector=yVector; %#ok<AGROW>
            if isempty(im) || isempty(xVector)
                bundles(b).Xmean=[]; %#ok<AGROW>
                bundles(b).Ymean=[]; %#ok<AGROW>
            else
                [~,sortIndex] = sort(im(:),'descend');
                maxIndices = sortIndex(1:5); %five brightest pixels
                [I,J]=ind2sub(size(im),maxIndices);
                I=median(I);
                J=median(J);
                bundles(b).Xmean=xVector(J); %#ok<AGROW>
                bundles(b).Ymean=yVector(I); %#ok<AGROW>
            end
        end
        %What is the largest bundle (LB)?
        [~,iLB]=max(numRays);   
        xVal = bundles(iLB).Xmean;
        yVal = bundles(iLB).Ymean;
  
        plot(xVal,yVal,'Marker','+','Color',[0.8 0.13 0.15],'LineStyle',...
            'none','MarkerSize',16)
        
        usefulIntersections=cell(numel(usefulrays),1);
        cnt=1;
        for r=1:numel(finalrays)
            if finalrays{r}.power > 0;
                usefulIntersections{cnt} = intersections{r};
                cnt=cnt+1;               
            end
        end
        usefulIntersections = {usefulIntersections{1};...
            usefulIntersections{end}};
        usefulrays = {usefulrays{1},usefulrays{end}};
        smartRayPlot(usefulIntersections, usefulrays,'r',1); hold all;

    case 'thresholded_but'
        set(handles.ShowRayTracerResults_fig,'Renderer','Painters')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Bin data according to dx and dy specified by the user
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dx=str2double(get(handles.dx_ctl,'String')); %mm
        dy=str2double(get(handles.dy_ctl,'String'));
        yMin=-80;
        yMax=20;
        set(handles.ShowRayTracerResults_fig,'Name',...
            sprintf('Binning %0.0f points...',numel(X)));
            pause(0.01)
        [im,xVector,yVector] = ...
            binnedImageFromIntersecs(X,Y,dx,dy,yMin,yMax);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Scale axes, remove previous plot
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        xlim=get(handles.axes1,'XLim');
        ylim=get(handles.axes1,'YLim');
        h=get(handles.axes1,'Children');
        delete(h);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Find the largest objects based on morphological ops
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        numObjects=5;
        %TUNA: Maybe not the largest, but the brightest total intensity?
        [bw Areas xVec yVec] = getLargestObjects(im,numObjects,...
            xVector,yVector);
        
        
        imagesc(bw.*im,'XData',[min(xVector),max(xVector)],'YData',...
            [min(yVector),max(yVector)],'AlphaData',1)
        colormap(GetLaVisionMap);
        hold all
        %Plot the centroids with stars
        plot(xVec, yVec,'r*','MarkerSize',12)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Find the median of the 5 brightest points
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        numBrightest=5;
        [xVal yVal] = getBrightestPix(im,numBrightest,xVector,yVector);
        %plot the median with an X
        plot(xVal,yVal,'Marker','x','Color',[0.8 0 0],'LineStyle',...
            'none','MarkerSize',16)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Plot the piston contour
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        surfaces=getappdata(handles.ShowRayTracerResults_fig,'surfaces');
        for s=1:length(surfaces)
           plot(surfaces{s}.x,surfaces{s}.y,'Color',[0.7 0.7 0.7],...
               'Linewidth',2); 
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Plot the start point
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        plot(startPoint(1),startPoint(2),'Color',[0.7 0.7 0.7],...
            'Marker','+','LineWidth',2,'MarkerEdgeColor',...
            [0.6 0.6 0.6],'MarkerSize',20);
        
        
        grid on
        set(handles.axes1,'XLim',xlim);
        set(handles.axes1,'YLim',ylim);
    case 'rays_but'
        %Plot the rays that end up being useful
        xlim=get(handles.axes1,'XLim');
        ylim=get(handles.axes1,'YLim');
        h=get(handles.axes1,'Children');
        delete(h);
        set(handles.ShowRayTracerResults_fig,'Renderer','painters')
        set(handles.axes1,'DrawMode','fast')
        usefulIntersections=cell(numel(usefulrays),1);
        cnt=1;
        for r=1:numel(finalrays)
            if finalrays{r}.power >0;
                usefulIntersections{cnt} = intersections{r};
                cnt=cnt+1;               
            end
        end
        smartRayPlot(usefulIntersections, usefulrays,'r'); hold all;
        set(handles.axes1,'XLim',[-50 50],'XTick',[-45, -30, -15,...
            0, 15, 30, 45]);
        set(handles.axes1,'YLim',[-80 20]);
        %Plot the initial rays on each end of the ray fan to help the user
        %know if the ray fan angle was too small or not.
        try
            endPoint1 = intersections{1,1}(:,1:2);
            endPoint2 = intersections{end,1}(:,1:2);
            hold all
            plot(endPoint1(1,:),endPoint1(2,:), 'Color', [0.6 0 0],...
                'LineWidth',2)
            plot(endPoint2(1,:),endPoint2(2,:), 'Color', [0.6 0 0],...
                'LineWidth',2)
        catch
        end
        drawnow
        %Find and plot intersections
%         [X,Y] = findCrossings(usefulrays);
%         plot(X,Y,'LineStyle','none','Color',[0 0 0],'Marker',...
%             'o','MarkerSize',4,'MarkerFaceColor',[0 0 0]);
        
    otherwise
            
    end
    
    
    set(handles.ShowRayTracerResults_fig,'Name','Ray Tracer Viewer')
    
function smartRayPlot(intersec,rays,color,varargin)
    %How many points total?
    a=0;
    for r = 1:numel(intersec)
       a=a+numel(intersec{r})/2; 
    end
    concatPts = zeros(a+numel(intersec),2);
    concatPts2 = zeros(a+2*numel(intersec),2);
%     concatPts3 = zeros(a+2*numel(intersec),2);
    cnt=1;
    cnt2=1;
%     cnt3=1;
    for r=1:numel(intersec)
        curPoint = intersec{r}';
        concatPts(cnt:cnt+size(curPoint,1),:) =...
            [curPoint; nan, nan];
        concatPts2(cnt2:cnt2+2,:) = [curPoint(end,:); ...
            (rays{r}.t_ray(30))'; nan,nan];
%         concatPts3(cnt3:cnt3+2,:) = [curPoint(end,:);...
%             (rays{r}.t_ray(-50))'; nan, nan];
        
        cnt=cnt+size(curPoint,1)+1;
        cnt2=cnt2+3;
%         cnt3=cnt3+3;
    end
    if isempty(concatPts)
        return
    end
    plot(concatPts(1,1),concatPts(1,2),'Color',color,'Marker','o',...
        'MarkerFaceColor',color,'MarkerSize',4,'LineWidth',0.25)
    plot(concatPts(:,1),concatPts(:,2),'Color',color,...
        'Marker','.','MarkerSize',1,'LineWidth',0.25)
    plot(concatPts2(:,1),concatPts2(:,2),'Color','b')
%     plot(concatPts3(:,1),concatPts3(:,2),'Color','k')
    if nargin == 4
        backproj = varargin{1};
    else
        backproj = false;
    end
    if backproj
        [X,Y] = findCrossings(rays);
        for r = 1:numel(intersec)
            ray = rays{r}; 
            start = ray.start;
            fin = [X(1),Y(1)];
            plot([start(1),fin(1)],[start(2),fin(2)],'Color',[0 0 0])
        end
    end
return
 
function plotTools_panel_SelectionChangeFcn(hObject, eventdata, handles)
    switch get(eventdata.NewValue,'Tag')
        case 'zoom_but'
            pan off
            zoom on
        case 'pan_but'
            zoom off
            pan on
        otherwise
            zoom off
            pan off
    end

function displayType_panel_SelectionChangeFcn(hObject, eventdata, handles)
redraw(handles)

function dx_ctl_Callback(hObject, eventdata, handles)
    minVal=0.005;
    maxVal=30;
    newVal=checkval(get(hObject,'String'),minVal,maxVal);
    set(hObject,'String',sprintf('%1.3f',newVal));
    redraw(handles)
    
function dy_ctl_Callback(hObject, eventdata, handles)
    minVal=0.005;
    maxVal=30;
    newVal=checkval(get(hObject,'String'),minVal,maxVal);
    set(hObject,'String',sprintf('%1.3f',newVal));
    redraw(handles)
    
function jetmap=GetLaVisionMap()
    %idx 0-16
    r1=linspace(0,0,17);
    g1=linspace(0,0,17);
    b1=linspace(1,224,17)/255;
    %idx 16-83
    r2=linspace(0,0,68);            r2=r2(2:end);
    g2=linspace(0,254,68)/255;      g2=g2(2:end);
    b2=linspace(224,255,68)/255;    b2=b2(2:end);
    %idx 83-166
    r3=linspace(0,255,84)/255;      r3=r3(2:end);
    g3=linspace(254,254,84)/255;    g3=g3(2:end);
    b3=linspace(255,0,84)/255;      b3=b3(2:end);
    %idx 166-251
    r4=linspace(255,255,86)/255;    r4=r4(2:end);
    g4=linspace(254,0,86)/255;      g4=g4(2:end);
    b4=linspace(0,0,86);            b4=b4(2:end);
    %
    jetmap=[r1 r2 r3 r4; g1 g2 g3 g4; b1 b2 b3 b4]';
    jetmap(end+1,:)=[252 252 251]/255;
    jetmap(end+1,:)=[255 255 255]/255;
    jetmap(end+1,:)=[255 255 255]/255;
    jetmap(end+1,:)=[255 255 255]/255;
    
function Map = getHSVlinearMap()
    Map = colormap([1,1,1;1,0.772549033164978,0.698039233684540;...
        0,0.749019622802734,0.749019622802734;...
        0,0.698039233684540,0.776470601558685;...
        0,0.643137276172638,0.803921580314636;...
        0,0.580392181873322,0.831372559070587;...
        0,0.513725519180298,0.858823537826538;...
        0,0.443137258291245,0.890196084976196;...
        0,0.364705890417099,0.917647063732147;...
        0,0.278431385755539,0.945098042488098;...
        0,0.188235297799110,0.972549021244049;...
        0,0.0941176488995552,1;...
        0,0.0196078438311815,1;...
        0.0549019612371922,0,1;...
        0.133333340287209,0,1;...
        0.207843139767647,0,1;...
        0.282352954149246,0,1;...
        0.352941185235977,0,1;...
        0.427450984716415,0,1;...
        0.501960813999176,0,1;...
        0.576470613479614,0,1;...
        0.647058844566345,0,1;...
        0.721568644046783,0,1;...
        0.796078443527222,0,1;...
        0.866666674613953,0,1;...
        0.941176474094391,0,1;...
        1,0,1;...
        0.992156863212585,0,0.925490200519562;...
        0.980392158031464,0,0.831372559070587;...
        0.968627452850342,0,0.737254917621613;...
        0.956862747669220,0,0.647058844566345;...
        0.945098042488098,0,0.556862771511078;...
        0.937254905700684,0,0.466666668653488;...
        0.925490200519562,0,0.384313732385635;...
        0.913725495338440,0,0.298039227724075;...
        0.901960790157318,0,0.215686276555061;...
        0.890196084976196,0,0.137254908680916;...
        0.886274516582489,0,0.109803922474384;...
        0.882352948188782,0,0.0862745121121407;...
        0.874509811401367,0,0.0156862754374743;...
        0.862745106220245,0.0509803928434849,0;...
        0.854901969432831,0.129411771893501,0;...
        0.850980401039124,0.211764708161354,0;...
        0.843137264251709,0.286274522542954,0;...
        0.839215695858002,0.364705890417099,0;...
        0.831372559070587,0.439215689897537,0;...
        0.823529422283173,0.513725519180298,0;...
        0.819607853889465,0.588235318660736,0;...
        0.811764717102051,0.658823549747467,0;...
        0.807843148708344,0.729411780834198,0;...
        0.800000011920929,0.800000011920929,0;...
        0.674509823322296,0.764705896377564,0.0274509806185961;...
        0.560784339904785,0.729411780834198,0.0509803928434849;...
        0.462745100259781,0.698039233684540,0.0745098069310188;...
        0.376470595598221,0.662745118141174,0.0941176488995552;...
        0.305882364511490,0.627451002597809,0.113725490868092;...
        0.243137255311012,0.592156887054443,0.125490203499794;...
        0.192156866192818,0.556862771511078,0.141176477074623;...
        0.149019613862038,0.521568655967712,0.149019613862038;...
        0.156862750649452,0.490196079015732,0.200000002980232;...
        0.160784319043159,0.454901963472366,0.235294118523598;...
        0.164705887436867,0.419607847929001,0.258823543787003;...
        0.164705887436867,0.384313732385635,0.274509817361832;...
        0,0,0]);
        
function val_out=checkval(str_in,min,max)
%First check to see if str_in is a number when converted with str2double.
%Then see if it is within the limits of min and max and coerce the value if
%need be.  Finally, output the value in val_out.

%Change commas to periods for easier data entry
str_in=strrep(str_in,',','.');

if isnan(str2double(str_in))
    val = min;
%     beep
%     errordlg('Wert muss numerisch sein!',...
%         'Fehler', 'modal');
    val_out=val;
else
    val = str2double(str_in);

if val > max
    val = max;
%     beep
%     errordlg(strcat(['Max Wert ist ' num2str(max) '!']),...
%         'Max Wert überschreiten', 'modal');
end

if val < min
    val = min;
%     beep
%     errordlg(strcat(['Min Wert ist ' num2str(min) '!']),...
%     'Min Wert unterschreiten', 'modal');
end

   val_out=val;     
end

function [vals, indices]=findMaxima(hP,threshold,separationLimit)
    %Find all local maxima
    diff_hP=diff(hP);
    if numel(hP)==1
        vals=hP;
        indices=1;
    elseif isempty(hP)
        vals=1;
        indices=1;
    else
        iAllMaxima = find([hP(1)>hP(2),...
            (diff_hP(2:end)<0 & diff_hP(1:(end-1))>0),...
            hP(end)>hP(end-1)]);
        %look at peaks, see which ones are greater than 10% of the maximum.
        hPmaxima=hP(iAllMaxima);
        iPeaksHighEnoughMaxima=iAllMaxima(hPmaxima>threshold*max(hP));
        %Look at indices - if several maxima are too close together, take
        %the highest one    
        iPeaksSeparated=[];
        isClean=false;
        current=1;
        while ~isClean
            %Look at each peak. If the next one in line is within the spacing,
            %remove the next one.
            try
                spacing = abs(iPeaksHighEnoughMaxima(current+1)-...
                    iPeaksHighEnoughMaxima(current));
                if spacing < separationLimit %one of them has to go
                    [val,iToKeep]=max(hP(iPeaksHighEnoughMaxima(current:...
                        current+1)));
                    if iToKeep==1 %remove the current+1 value
                        iPeaksHighEnoughMaxima=...
                            [iPeaksHighEnoughMaxima(1:current),...
                            iPeaksHighEnoughMaxima(current+2:end)];
                    else %remove the current value
                        if current == 1
                            iPeaksHighEnoughMaxima =...
                                iPeaksHighEnoughMaxima(2:end);
                        else
                            iPeaksHighEnoughMaxima=...
                                [iPeaksHighEnoughMaxima(1:current-1),...
                                iPeaksHighEnoughMaxima(current+1:end)];
                        end
                    end
                else
                    current = current+1;
                end
               if current == length(iPeaksHighEnoughMaxima)
                  isClean=true; 
               end
            catch
                indices=[];
                vals=[];
                return
            end
        end


    %     for idx=1:length(iPeaksHighEnoughMaxima)
    %        spacings=abs(iPeaksHighEnoughMaxima(idx)-...
    %            iPeaksHighEnoughMaxima);
    %        closeIndices=find(spacings<separationLimit);
    %        %find highest peak of the closeIndices values
    %        [temp,iHighestOfClosePeaks]=...
    %            max(hP(iPeaksHighEnoughMaxima(closeIndices)));
    %       iPeaksSeparated(end+1)=closeIndices(iHighestOfClosePeaks);%#ok<AGROW>
    %     end
        indices = iPeaksHighEnoughMaxima;
        vals = hP(indices);
    end
    
function [im, xVector, yVector] = binData(X,Y,dx,dy,type)
    ncols=ceil((max(X)-min(X))/dx);
    nrows=ceil((max(Y)-min(Y))/dy);
    try
        colBins=linspace(min(X),max(X),ncols);
        rowBins=linspace(min(Y),max(Y),nrows);
    catch err %#ok<NASGU>
        im=[];
        xVector=[];
        yVector=[];
        return
    end
switch type
    case 'slow'
        %don't use this. It's really slow.
        im=zeros(nrows,ncols);
        for idx=1:numel(X)
            xi = find((X(idx) >= colBins), 1, 'last');
            yi = find((Y(idx) >= rowBins), 1, 'last');
            if ~isempty(yi) && ~isempty(xi)
                im(yi,xi)=im(yi,xi)+1;
            end
        end
        xVector=[min(X),max(X)];
        yVector=[min(Y),max(Y)];
    case 'fast'
        %[h2D, binC_x, binC_y, Edges_MinMaxN_XY] = 
        %hist2(Data_XY, Edges_MinMaxN_XY, BIN_MODE, EPSILON, GRID_ADJ_MODE)
        if ncols == 1
            edgesMinMaxN_XY=[mean(X)-dx/2,mean(X)+dx/2,2;min(Y),max(Y),nrows];
        else
            edgesMinMaxN_XY=[min(X),max(X),ncols;min(Y),max(Y),nrows];
        end
%         edgesMinMaxN_XY=[ncols;nrows];
    if isempty(X) || isempty(edgesMinMaxN_XY)
        xVector = [];
        yVector = [];
        im = [];
    else
        [im, xVector, yVector]=hist2([X',Y'],edgesMinMaxN_XY);
        im=double(im);
    end
    otherwise
        
end

function [bw, Areas, xVec, yVec] =...
    getLargestObjects(im,numObjects,xVector,yVector)
    try
        %Threshold image
        bw=(im>0);
        se = strel('disk',2);
        bw = imclose(bw,se);
        bw=medfilt2(bw,[3,3]);
        cc = bwlabel(bw);
        stats=regionprops(cc,'Area','Centroid','PixelIdxList','PixelList');
        Areas=zeros(numel(stats),1);
        Centroids=zeros(numel(stats),2);
        Objectnums=zeros(numel(stats),1);
        xbar=zeros(numel(stats),1);
        ybar=zeros(numel(stats),1);
        for obj=1:numel(stats)
            Areas(obj)=stats(obj).Area;
            Centroids(obj,:)=stats(obj).Centroid;
            Objectnums(obj)=obj;
            idx=stats(obj).PixelIdxList;
            pixelValues = double(im(idx));
            sumPixelValues = sum(pixelValues);
            x = stats(obj).PixelList(:,1);
            y = stats(obj).PixelList(:,2);
            xbar(obj) = sum(x .* pixelValues) / sumPixelValues;
            ybar(obj) = sum(y .* pixelValues) / sumPixelValues;
        end
        %Sort
        [Areas I]=sort(Areas,1,'descend');
        Centroids=Centroids(I,:);
        Objectnums=Objectnums(I);
        xbar=xbar(I);
        ybar=ybar(I);
        %Take the 4 largest objects (or the largest 1, 2, or 3 if there
        %aren't 4)
        if numel(stats) > numObjects
            Areas=Areas(1:numObjects);
            Centroids=Centroids(1:numObjects,:);
            Objectnums=Objectnums(1:numObjects);
            xbar=xbar(1:numObjects);
            ybar=ybar(1:numObjects);
        else
            Areas=Areas(1:end);
            Centroids=Centroids(1:end,:);
            Objectnums=Objectnums(1:end);
            xbar=xbar(1:end);
            ybar=ybar(1:end);
        end
        bw=ismember(cc,Objectnums);
        xVec=linspace(min(xVector),max(xVector),size(bw,2));
        yVec=linspace(min(yVector),max(yVector),size(bw,1));
    %     xVec=xVec(round(Centroids(:,1)));
    %     yVec=yVec(round(Centroids(:,2)));
        xVec=xVec(round(xbar));
        yVec=yVec(round(ybar));
    catch
        bw = zeros(size(im));
        Areas = [];
        xVec = [];
        yVec = [];
    end
    
function [im,xVector,yVector] = ...
    binnedImageFromIntersecs(X,Y,dx,dy,yMin,yMax) 
    %Limit range based on y criteria
    X=X(Y < yMax & Y > yMin);
    Y=Y(Y < yMax & Y > yMin);
    [im,xVector,yVector] = binData(X,Y,dx,dy, 'fast');
    
function [xVal yVal] = getBrightestPix(im,numBrightest,xVector,yVector)
    [sortedValues,sortIndex] = sort(im(:),'descend');
    if numel(sortIndex)>numBrightest
        maxIndices = sortIndex(1:numBrightest);
    else
        maxIndices = sortIndex(1:end);
    end
    [I,J]=ind2sub(size(im),maxIndices);
    xVal=xVector(median(J));
    yVal=yVector(median(I));
    
function demoVid_but_Callback(hObject, eventdata, handles)
    %Loop through all files, show all useful rays, save to gif
    [saveFile, savePath] = uiputfile('*.gif', 'Save as...');
    if isequal(saveFile,0) || isequal(savePath,0)
       return
    end
    
    %Loop through each file for the given crank angle
    dirName=getappdata(handles.ShowRayTracerResults_fig,'dirName');
    %These are the filenames for the selected crank angle.
    filenames=getappdata(handles.ShowRayTracerResults_fig,'filenames');
    %Loop through each file, display the rays, and save each image to a
    %gif/video/something
    xVar=[];
    yVar=[];
    for idx=1:length(filenames)-1
        %Load the file
        set(handles.ShowRayTracerResults_fig,'Name',...
            sprintf('Loading %s...',filenames{idx}))
        pause(0.01)
        S=load(fullfile(dirName,filenames{idx}));
        X=S.X;
        Y=S.Y;
        usefulrays=S.usefulrays;
        startPoint=S.intersections{1,:};
        intersections=S.intersections;
        finalrays=S.finalrays;
        clear S
        set(handles.ShowRayTracerResults_fig,'Name','Plotting...')
        pause(0.01)
        
        %Plot the rays that end up being useful
        xlim=get(handles.axes1,'XLim');
        ylim=get(handles.axes1,'YLim');
        h=get(handles.axes1,'Children');
        delete(h);

        plot(intersections{1}(1,1),intersections{1}(2,1),...
            'Color',[1 0 0],'Marker','+','MarkerSize',20,'LineWidth',2);
        for r=1:numel(finalrays)
            if finalrays{r}.power > 0; %plot rays that escape unabsorbed
                plotpath(intersections{r},'r');
%                 plotpath([intersections{r}(:,end),...
%                     finalrays{r}.t_ray(-100)],'k');
                plotpath([intersections{r}(:,end),...
                    finalrays{r}.t_ray(100)],'b');
            end
        end
        %Find the crossings and plot them too.
        [X,Y] = findCrossings(usefulrays);
        plot(X,Y,'LineStyle','none','Color',[0 0 0],'Marker',...
            'o','MarkerSize',4,'MarkerFaceColor',[0 0 0]);
        %Draw a line to show the mapping function
        xVar=[xVar;intersections{1}(1,1)]; %#ok<AGROW>
        yVar=[yVar;mean(X)]; %#ok<AGROW>
        plot(xVar,yVar,'Color',[0.4 0.4 0.4],'LineWidth',2)
        
        set(handles.axes1,'XLim',[-50 50],'XTick',[-45, -30, -15,...
            0, 15, 30, 45]);
        set(handles.axes1,'YLim',[-80 20]);
        drawnow
        frame = getframe(handles.axes1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if idx == 1;
          imwrite(imind,cm,fullfile(savePath,saveFile),'gif',...
              'Loopcount',1,'DelayTime',0.1);
        else
          imwrite(imind,cm,fullfile(savePath,saveFile),'gif',...
              'WriteMode','append','DelayTime',0.1);
        end
        
    end
    
function demoVid2_but_Callback(hObject, eventdata, handles)
    %For current CA and file (radial position), go through the useful rays
    %and display them one by one, saving a frame to a gif each time
    %Based on piston position and radial position, build up the filename to
    %load and display
    
    %Where to save?
    [saveFile, savePath] = uiputfile('*.gif', 'Save as...');
    if isequal(saveFile,0) || isequal(savePath,0)
%        return
        Write=false;
    else
        Write=true;
    end
    
    filenames=getappdata(handles.ShowRayTracerResults_fig,'filenames');
    filename=filenames{get(handles.files_listbox,'Value')};
    dirName=getappdata(handles.ShowRayTracerResults_fig,'dirName');
    %Load data
    set(handles.ShowRayTracerResults_fig,'Name','Loading file...')
    pause(0.01)
    S=load(fullfile(dirName,filename));
    X=S.X;
    Y=S.Y;
    usefulrays=S.usefulrays;
    startPoint=S.intersections{1,:};
    intersections=S.intersections;
    finalrays=S.finalrays;
    clear S
    %Clear the axes
    xlim=get(handles.axes1,'XLim');
    ylim=get(handles.axes1,'YLim');
    h=get(handles.axes1,'Children');
    delete(h);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %First, plot a single ray piece by piece
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rayNum=1735; %GM, CA-40, R8, h10.426
    %start with starting point
    plot(intersections{rayNum}(1,1),intersections{rayNum}(2,1),'Color',...
        [1 0 0],'Marker','+','MarkerSize',20,'LineWidth',2);
    frame = getframe(handles.axes1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if Write
        imwrite(imind,cm,fullfile(savePath,saveFile),'gif',...
            'Loopcount',1,'DelayTime',1);
    else
        pause(1)
    end
    %First leg
    plot(intersections{rayNum}(1,1:2),intersections{rayNum}(2,1:2),...
        'Color',[1 0 0],'Marker','o','MarkerSize',8,'LineWidth',1.5);
    frame = getframe(handles.axes1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if Write
        imwrite(imind,cm,fullfile(savePath,saveFile),'gif',...
              'WriteMode','append','DelayTime',0.5);
    else
        pause(0.5)
    end
    %Second Leg
    plot(intersections{rayNum}(1,2:3),intersections{rayNum}(2,2:3),...
        'Color',[1 0 0],'Marker','o','MarkerSize',8,'LineWidth',1.5);
    frame = getframe(handles.axes1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if Write
        imwrite(imind,cm,fullfile(savePath,saveFile),'gif',...
              'WriteMode','append','DelayTime',0.5);
    else
        pause(0.5)
    end
    %Last Leg
    endPt = finalrays{rayNum}.t_ray(40);
    plot([intersections{rayNum}(1,end),endPt(1)],...
        [intersections{rayNum}(2,end),endPt(2)],...
        'Color',[0 0 1],'Marker','o','MarkerSize',8,'LineWidth',1.5);
    set(handles.axes1,'XLim',[-50 50],'XTick',[-45, -30, -15,...
        0, 15, 30, 45]);
    set(handles.axes1,'YLim',[-80 20]);
    frame = getframe(handles.axes1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if Write
        imwrite(imind,cm,fullfile(savePath,saveFile),'gif',...
              'WriteMode','append','DelayTime',2.5);
    else
        pause(2.5)
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %Plot a second ray faster
    %%%%%%%%%%%%%%%%%%%%%%%%%
    rayNum2=1805;
    %First leg
    plot(intersections{rayNum2}(1,1:2),intersections{rayNum2}(2,1:2),...
        'Color',[1 0 0],'Marker','o','MarkerSize',8,'LineWidth',1.5);
    frame = getframe(handles.axes1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if Write
        imwrite(imind,cm,fullfile(savePath,saveFile),'gif',...
              'WriteMode','append','DelayTime',0.5);
    else
        pause(0.5)
    end
    %Second Leg
    plot(intersections{rayNum2}(1,2:3),intersections{rayNum2}(2,2:3),...
        'Color',[1 0 0],'Marker','o','MarkerSize',8,'LineWidth',1.5);
    frame = getframe(handles.axes1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if Write
        imwrite(imind,cm,fullfile(savePath,saveFile),'gif',...
              'WriteMode','append','DelayTime',0.5);
    else
        pause(0.5)
    end
    %Last Leg
    endPt = finalrays{rayNum2}.t_ray(40);
    plot([intersections{rayNum2}(1,end),endPt(1)],...
        [intersections{rayNum2}(2,end),endPt(2)],...
        'Color',[0 0 1],'Marker','o','MarkerSize',8,'LineWidth',1.5);
    set(handles.axes1,'XLim',[-50 50],'XTick',[-45, -30, -15,...
        0, 15, 30, 45]);
    set(handles.axes1,'YLim',[-80 20]);
    frame = getframe(handles.axes1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if Write
        imwrite(imind,cm,fullfile(savePath,saveFile),'gif',...
              'WriteMode','append','DelayTime',5);
    else
        pause(5)
    end
    %Plot the intersection of the two rays
    [X,Y] = findCrossings({finalrays{rayNum},finalrays{rayNum2}});
    plot(X,Y,'LineStyle','none','Color',[0 0 0],'Marker',...
        'o','MarkerSize',8,'MarkerFaceColor',[0 0 0]);
    endPt2 = finalrays{rayNum2}.t_ray(-40);
    plot([intersections{rayNum2}(1,end),endPt2(1)],...
        [intersections{rayNum2}(2,end),endPt2(2)],...
        'Color',[0 0 0],'LineWidth',1.5);
    endPt = finalrays{rayNum}.t_ray(-40);
    plot([intersections{rayNum}(1,end),endPt(1)],...
        [intersections{rayNum}(2,end),endPt(2)],...
        'Color',[0 0 0],'LineWidth',1.5);
    frame = getframe(handles.axes1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if Write
        imwrite(imind,cm,fullfile(savePath,saveFile),'gif',...
              'WriteMode','append','DelayTime',5);
    else
        pause(5)
    end    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Clear plot, make ray fan (one total ray path at a time)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h=get(handles.axes1,'Children');
    delete(h);
    plot(intersections{rayNum}(1,1),intersections{rayNum}(2,1),'Color',...
        [1 0 0],'Marker','+','MarkerSize',20,'LineWidth',2);
    for r=2:60:numel(finalrays)
            plotpath(intersections{r},'r');
            
            set(handles.axes1,'XLim',[-50 50],'XTick',[-45, -30, -15,...
                0, 15, 30, 45]);
            set(handles.axes1,'YLim',[-80 20]);
            drawnow
            frame = getframe(handles.axes1);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            if Write
                imwrite(imind,cm,fullfile(savePath,saveFile),'gif',...
                      'WriteMode','append','DelayTime',0.1);
            else
                pause(0.1)
            end    
    end
    drawnow
    if Write
        imwrite(imind,cm,fullfile(savePath,saveFile),'gif',...
              'WriteMode','append','DelayTime',5);
    else
        pause(5)
    end    
    %Clear plot, plot only useful rays
    h=get(handles.axes1,'Children');
    delete(h);
    for r=rayNum:5:rayNum2
        if finalrays{r}.power > 0;
            plotpath(intersections{r},'r');
            plotpath([intersections{r}(:,end),...
                finalrays{r}.t_ray(40)],'b');
        end
        set(handles.axes1,'XLim',[-50 50],'XTick',[-45, -30, -15,...
            0, 15, 30, 45]);
        set(handles.axes1,'YLim',[-80 20]);
        drawnow
        frame = getframe(handles.axes1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if Write
            imwrite(imind,cm,fullfile(savePath,saveFile),'gif',...
                  'WriteMode','append','DelayTime',0.1);
        else
            pause(0.1)
        end    
    end
    drawnow
    if Write
        imwrite(imind,cm,fullfile(savePath,saveFile),'gif',...
              'WriteMode','append','DelayTime',2);
    else
        pause(2)
    end   
    %Plot the intersections of the useful rays    
    [X,Y] = findCrossings(finalrays(rayNum:5:rayNum2));
    plot(X,Y,'LineStyle','none','Color',[0 0 0],'Marker',...
        'o','MarkerSize',6,'MarkerFaceColor',[0 0 0]);
    drawnow
    frame = getframe(handles.axes1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if Write
        imwrite(imind,cm,fullfile(savePath,saveFile),'gif',...
              'WriteMode','append','DelayTime',5);
    else
        pause(5)
    end   
   
function demoImage_but_Callback(hObject, eventdata, handles)
    %Loop through all files, show all useful rays, save to gif
%     [saveFile, savePath] = uiputfile('*.gif', 'Save as...');
%     if isequal(saveFile,0) || isequal(savePath,0)
%        return
%     end
    
    %Loop through each file for the given crank angle
    dirName=getappdata(handles.ShowRayTracerResults_fig,'dirName');
    %These are the filenames for the selected crank angle.
    filenames=getappdata(handles.ShowRayTracerResults_fig,'filenames');
    %Loop through each file, display the rays, and save each image to a
    %gif/video/something
    xVar=[];
    yVar=[];
    dx=str2double(get(handles.dx_ctl,'String')); %mm
    dy=str2double(get(handles.dy_ctl,'String'));
    h=get(handles.axes1,'Children');
    delete(h);
    cnt=1;
    cnt3 = 1;
    for idx=[6,11,16,21,26,31,34,37,40]
        %Load the file
        set(handles.ShowRayTracerResults_fig,'Name',...
            sprintf('Loading %s...',filenames{idx}))
        pause(0.01)
        S=load(fullfile(dirName,filenames{idx}));
        X=S.X;
        Y=S.Y;
        usefulrays=S.usefulrays;
        startPoint=S.intersections{1,:};
        intersections=S.intersections;
        finalrays=S.finalrays;
        clear S
        set(handles.ShowRayTracerResults_fig,'Name','Plotting...')
        pause(0.01)
        
        %Plot the useful rays for each position at the given crank angle
        xlim=get(handles.axes1,'XLim');
        ylim=get(handles.axes1,'YLim');
        usefulIntersections=cell(numel(usefulrays),1);
        jj=1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        numUsefulRays=numel(usefulrays); %how many of them are useful?
        usefulRayIndices=zeros(numUsefulRays,1);%Initialize array of indices
        %loop though final rays and find the ones with positive power.
        %Store their indices in usefulRayIndices
        cnt2=1;
        for ray=1:length(finalrays)
           if finalrays{ray}.power > 0
              usefulRayIndices(cnt2) = ray;
              cnt2=cnt2+1;
           end
        end
        %Compute the spacing between rays. If the value is 1, then the rays
        %are adjacent
        raySpacing=[1;diff(usefulRayIndices)];
        %Find the locations where the ray spacing is greater than 1
        breaks=find(raySpacing > 1);
        if isempty(breaks)
%             breaks=1;
        end
        %Get the ray origin from any ray's intersections
        temp = intersections{1};
        xRayOrigin(idx) = temp(1,1);
        yRayOrigin(idx) = temp(2,1);
        
        numBundles=numel(breaks)+1;
        numRays=numel(breaks);
        if numBundles ==1
            bundles(1).indices = usefulRayIndices;

        else
            for b = 1:numBundles
                if b == 1
                    bundles(b).indices=...
                        usefulRayIndices(1:breaks(b)-1);%#ok<AGROW>
                elseif b == numBundles
                    bundles(b).indices=...
                        usefulRayIndices(breaks(b-1):end); %#ok<AGROW>
                else
                    bundles(b).indices=...
                        usefulRayIndices(breaks(b-1):breaks(b)-1);  %#ok<AGROW>
                end
                numRays(b)=numel(bundles(b).indices); 
            end
        end
        %Limit locations of intersections
        yMin=-60;
        yMax=10;
        %Go through each bundle and compute its intersections
        for b=1:numBundles
            [X,Y] = findCrossings(finalrays(bundles(b).indices));
            
            X=X(Y < yMax & Y > yMin); %limit range of y data
            Y=Y(Y < yMax & Y > yMin);
            if isempty(X) || isempty(Y)
               numRays(b)=0;
               continue; 
            end
            
            bundles(b).X=X; %#ok<AGROW>
            bundles(b).Y=Y; %#ok<AGROW>
            [im,xVector,yVector] = binData(X,Y,dx,dy,'fast');
            bundles(b).xVector=xVector; %#ok<AGROW>
            bundles(b).yVector=yVector; %#ok<AGROW>
            if isempty(im) || isempty(xVector)
                bundles(b).Xmean=[]; %#ok<AGROW>
                bundles(b).Ymean=[]; %#ok<AGROW>
            else
                [~,sortIndex] = sort(im(:),'descend');
                if numel(sortIndex) < 5
                    maxIndices = sortIndex(1:end);
                else
                    maxIndices = sortIndex(1:5); %five brightest pixels
                end
                [I,J]=ind2sub(size(im),maxIndices);
                I=median(I);
                J=median(J);
                bundles(b).Xmean=xVector(J); %#ok<AGROW>
                bundles(b).Ymean=yVector(round(I)); %#ok<AGROW>

            end
        end
        %What is the largest bundle (LB)?
        [~,iLB]=max(numRays);        
  
        try
            RvirtualImage(cnt,1)=bundles(iLB).Xmean;
            xVirtualImage(cnt,1)=bundles(iLB).Xmean;
            yVirtualImage(cnt,1)=bundles(iLB).Ymean;
            altNumRays=numRays;
            altNumRays(altNumRays==max(numRays))=0;
            [~,i2ndLB]=max(altNumRays); 
            RvirtualImage(cnt,2)=bundles(i2ndLB).Xmean;
            xVirtualImage(cnt,2)=bundles(i2ndLB).Xmean;
            yVirtualImage(cnt,2)=bundles(i2ndLB).Ymean;
            altNumRays(altNumRays==max(altNumRays))=0;
            [~,i3rdLB]=max(altNumRays); 
            RvirtualImage(cnt,3)=bundles(i3rdLB).Xmean;
            xVirtualImage(cnt,3)=bundles(i3rdLB).Xmean;
            yVirtualImage(cnt,3)=bundles(i3rdLB).Ymean;
            if ~isequal(xVirtualImage(cnt,1),xVirtualImage(cnt,2))
                XX(cnt3) = xVirtualImage(cnt,1);
                YY(cnt3) = yVirtualImage(cnt,1);
                cnt3 = cnt3+1;
                XX(cnt3) = xVirtualImage(cnt,2);
                YY(cnt3) = yVirtualImage(cnt,2);
            else
                XX(cnt3) = xVirtualImage(cnt,1);
                YY(cnt3) = yVirtualImage(cnt,1);
            end
            cnt = cnt+1;
            cnt3 = cnt3+1;
        catch
            %No point to calculate
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for r=1:numel(finalrays)
            if finalrays{r}.power >0;
                usefulIntersections{jj} = intersections{r};
                jj=jj+1;               
            end
        end
        smartRayPlot(usefulIntersections, usefulrays,'r'); hold all;
%         smartRayPlot(usefulIntersections, finalrays,'r'); hold all;
        drawnow
%         for r=1:10:numel(finalrays)
%             if finalrays{r}.power > 0; %plot rays that escape unabsorbed
%                 plotpath(intersections{r},'r');
% %                 plotpath([intersections{r}(:,end),...
% %                     finalrays{r}.t_ray(-100)],'k');
%                 plotpath([intersections{r}(:,end),...
%                     finalrays{r}.t_ray(100)],'b');
%             end
%         end
%         [X,Y] = findCrossings(usefulrays);
%         dx=str2double(get(handles.dx_ctl,'String')); %mm
%         dy=str2double(get(handles.dy_ctl,'String'));
%         yMin=-80;
%         yMax=20;
%         X=X(Y < yMax & Y > yMin); %limit range of y data
%         Y=Y(Y < yMax & Y > yMin);
%         
%         set(handles.ShowRayTracerResults_fig,'Name',...
%         sprintf('Binning %0.0f points...',numel(X)));
%         drawnow
%         [im,xVector,yVector] = binData(X,Y,dx,dy, 'fast'); 
%         if isempty(im)
%             continue
%         end
%         numBrightest=5;
%         try
%             [xVal,yVal] = getBrightestPix(im,numBrightest,xVector,yVector);
%             XX(cnt)=xVal;
%             YY(cnt)=yVal;
%             cnt=cnt+1;
%         catch err
%             disp(err.message);
%         end
        
        
        drawnow
        frame = getframe(handles.axes1);
        im = frame2im(frame);
        [imInd,map] = rgb2ind(im,256);
        if idx == 1;
%           imwrite(imInd,map,fullfile(savePath,saveFile),'gif',...
%               'Loopcount',1,'DelayTime',0.1);
        else
%           imwrite(imInd,map,fullfile(savePath,saveFile),'gif',...
%               'WriteMode','append','DelayTime',0.1);
        end
        
    end
    plot(XX,YY,'Marker','o','Color',[0 0 0],'LineStyle',...
        'none','MarkerSize',5,'MarkerFaceColor',[0 0 0])
    figure(1)
    ax=axes;
    
    bowl='GM';
    switch bowl
        case 'GM'
            surfaces=getGMgeometry(handles);
        case 'Ford'
            surfaces=getFordgeometry(handles);
    end

    %Now plot all the surfaces
    set(handles.ShowRayTracerResults_fig,'CurrentAxes',handles.axes1);
    setappdata(handles.ShowRayTracerResults_fig,'surfaces',surfaces);
    for s=1:length(surfaces)
       h=plot(surfaces{s}.x,surfaces{s}.y,'k-','Linewidth',2); 
       hold all
       set(h,'HandleVisibility','off')
    end
    grid on
    set(handles.axes1,'XLim',[-45 45]);
    set(handles.axes1,'YLim',[-60 30]);
    newHandle = copyobj(get(handles.axes1,'Children'),ax);
    set(ax,'Units','Pixels','Position',[100 100 700 1050])
    set(ax,'XLim',[-42 42],'XTick',[-45, -30, -15,...
            0, 15, 30, 45],'YLim',[-64 62]);
    xlabel('r (mm)','FontName','Arial','FontWeight','bold','FontSize',14)
    ylabel('z (mm)','FontName','Arial','FontWeight','bold','FontSize',14)

function showHist_but_Callback(hObject, eventdata, handles)
    dirName=getappdata(handles.ShowRayTracerResults_fig,'dirName');
    %These are the filenames for the selected crank angle:
    filenames=getappdata(handles.ShowRayTracerResults_fig,'filenames');
    dx=str2double(get(handles.dx_ctl,'String')); %mm
    dy=str2double(get(handles.dy_ctl,'String'));
    
    nBundlesMax=3;
    RvirtualImage=zeros(nBundlesMax,1);
    xVirtualImage=zeros(length(filenames),nBundlesMax);
    yVirtualImage=zeros(length(filenames),nBundlesMax);
        
    fileNum=getappdata(handles.ShowRayTracerResults_fig,'fileNum');
    setappdata(handles.ShowRayTracerResults_fig,'fileNum',fileNum+1);

        
    idx = get(handles.files_listbox,'Value');
    filename=filenames{idx};
    temp1=strfind(filename,'R');
    temp1=temp1(end);
    radiusString=filename(temp1+1:end-4);
    radiusString=strrep(radiusString,'p','.');
    RrealSpace=str2double(radiusString);
    %Load data
    set(handles.ShowRayTracerResults_fig,'Name',['Loading ',...
        filename, '...'])
    drawnow
    S=load(fullfile(dirName,filename));

    finalRays=S.finalrays; %grab all of the final rays
    intersections=S.intersections;
    numUsefulRays=numel(S.usefulrays); %how many of them are useful?
    usefulRayIndices=zeros(numUsefulRays,1);%Initialize array of indices
    clear S
    %loop though final rays and find the ones with positive power.
    %Store their indices in usefulRayIndices
    cnt=1;
    for ray=1:length(finalRays)
       if finalRays{ray}.power > 0
          usefulRayIndices(cnt) = ray;
          cnt=cnt+1;
       end
    end
    %Compute the spacing between rays. If the value is 1, then the rays
    %are adjacent
    raySpacing=[1;diff(usefulRayIndices)];
    %Find the locations where the ray spacing is greater than 1
    breaks=find(raySpacing > 1);
    if isempty(breaks)
        breaks=1;
    end
    %Get the ray origin from any ray's intersections
    temp = intersections{1};
    xRayOrigin = temp(1,1);
    yRayOrigin = temp(2,1);

    numBundles=numel(breaks);
    numRays=numel(breaks);
    if numBundles ==1
        bundles(1).indices = usefulRayIndices;

    else
        for b = 1:numBundles
            if b == 1
                bundles(b).indices=...
                    usefulRayIndices(1:breaks(b)-1);%#ok<AGROW>
            else
                bundles(b).indices=...
                    usefulRayIndices(breaks(b-1):breaks(b)-1);  %#ok<AGROW>
            end
            numRays(b)=numel(bundles(b).indices); 
        end
    end
    %Limit locations of intersections
    yMin=-60;
    yMax=10;
    %Go through each bundle and compute its intersections
    for b=1:numBundles
        [X,Y] = findCrossings(finalRays(bundles(b).indices));

        X=X(Y < yMax & Y > yMin); %limit range of y data
        Y=Y(Y < yMax & Y > yMin);
        if isempty(X) || isempty(Y)
           numRays(b)=0;
           continue; 
        end

        bundles(b).X=X; %#ok<AGROW>
        bundles(b).Y=Y; %#ok<AGROW>
        [im,xVector,yVector] = binData(X,Y,dx,dy,'fast');
        bundles(b).xVector=xVector; %#ok<AGROW>
        bundles(b).yVector=yVector; %#ok<AGROW>
        if isempty(im) || isempty(xVector)
            bundles(b).Xmean=[]; %#ok<AGROW>
            bundles(b).Ymean=[]; %#ok<AGROW>
        else
            [~,sortIndex] = sort(im(:),'descend');
            maxIndices = sortIndex(1:5); %five brightest pixels
            [I,J]=ind2sub(size(im),maxIndices);
            I=median(I);
            J=median(J);
            bundles(b).Xmean=xVector(J); %#ok<AGROW>
            bundles(b).Ymean=yVector(I); %#ok<AGROW>
        end
    end
    %What is the largest bundle (LB)?
    [~,iLB]=max(numRays);        

    try
        RvirtualImage(1)=bundles(iLB).Xmean;
        xVirtualImage(1)=bundles(iLB).Xmean;
        yVirtualImage(1)=bundles(iLB).Ymean;
        altNumRays=numRays;
        altNumRays(altNumRays==max(numRays))=0;
        [~,i2ndLB]=max(altNumRays); 
        RvirtualImage(2)=bundles(i2ndLB).Xmean;
        xVirtualImage(2)=bundles(i2ndLB).Xmean;
        yVirtualImage(2)=bundles(i2ndLB).Ymean;
        altNumRays(altNumRays==max(altNumRays))=0;
        [~,i3rdLB]=max(altNumRays); 
        RvirtualImage(3)=bundles(i3rdLB).Xmean;
        xVirtualImage(3)=bundles(i3rdLB).Xmean;
        yVirtualImage(3)=bundles(i3rdLB).Ymean;
    catch
        %No point to calculate
    end
    %We are interested only in the radial (x) distribution.
    radialHist = sum(im,1)/(sum(sum(im)))*100;
    figure(42)
    plot(xVector,radialHist,'Marker','o','Color',[0.22 0.31 0.64])
    xlabel('Radial position of virtual image (mm)','FontName','Arial',...
        'FontWeight','bold','FontSize',14)
    ylabel('Frequency (%)','FontName','Arial','FontWeight','bold',...
        'FontSize',14)
    set(gca,'FontName','Arial','FontWeight','bold','FontSize',14)
    hold all
    plot([xVirtualImage(1),xVirtualImage(1)],[0 max(radialHist)+1],...
        'LineWidth',2,'Color',[0 0 0])
    
    
        

function ImportPistonPosition(fileToRead1)
    %IMPORTFILE(FILETOREAD1)
    %  Imports data from the specified file
    %  FILETOREAD1:  file to read

    %  Auto-generated by MATLAB on 26-Jan-2014 14:58:10

    % Import the file
    newData1 = importdata(fileToRead1);

    % Create new variables in the base workspace from those fields.
    vars = fieldnames(newData1);
    for i = 1:length(vars)
        assignin('caller', vars{i}, newData1.(vars{i}));
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create Functions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function files_listbox_CreateFcn(hObject, eventdata, handles)

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function crankAngle_listbox_CreateFcn(hObject, eventdata, handles)

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function dx_ctl_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function dy_ctl_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Deprecated%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [maxel,IJ]= max2(M,userows,usecols)
% finds the location of the single overall maximum element in a 2-d array
% usage: [maxel,IJ] = max2(M)
% usage: [maxel,IJ] = max2(M,userows,usecols)
%
% The location in a 2-d array of the overall
% maximum element (or the first incidence of
% several, if the maximum is not unique), where
% you may restrict the search to a set of
% specified rows and/or columns.
%
% Note that max2 does NOT convert the matrix to
% linear indexing, so that really huge arrays
% can be worked with.
%
% arguments: (input)
%  M - an (nxm) 2-dimensional numeric array (or
%      vector) that max is able to operate on. M
%      may contain inf or -inf elements.
%
%  userows - (OPTIONAL) a list of the rows to be
%      searched for the maximum. The search will
%      be restricted to this set of rows. If empty.
%      there will be no row restriction.
%
%      userows must be a list of integers
%
%  usecols - (OPTIONAL) a list of the columns to be
%      searched for the maximum. The search will
%      be restricted to this set of columnss. If
%      empty. there will be no column restriction.
%
% arguments: (output)
%  maxel - overall maximum element found. If the
%      maximum was ot unique, then this is the
%      first element identified. Ties will be
%      resolved in a way consistent with find.
%
%  IJ - a (1x2) row vector, comtaining respectively
%      the row and column indices of the maximum as
%      found.
%
% Example:
%  M = magic(4)
% ans =
%    16     2     3    13
%     5    11    10     8
%     9     7     6    12
%     4    14    15     1
%
% % the overall maximum
%  [maxel,IJ] = max2(M)
% maxel =
%      16
% IJ =
%      1     1
%
%
% % a restricted maximum
%  [maxel,IJ] = max2(M,[1 2 3],[2 3])
% maxel =
%      11
% IJ =
%     2     2
%
%
% See also: max2, max, min, find
% 
%
% Author: John D'Errico
% e-mail: woodchips@rochester.rr.com
% Release: 1.0
% Release date: 2/16/09

% check the arguments
if (nargin<1) || (nargin>3)
  error('max2 may have 1, 2, or 3 arguments only')
end

if length(size(M)) > 2
  error('M must be a 2-d array or a vector')
end
[n,m] = size(M);

% default for userows?
if (nargin<2) || isempty(userows)
  userows = 1:n;
else
  userows = unique(userows);
  if ~isnumeric(userows) || any(diff(userows)==0) || ...
      any(userows<1) || any(userows>n) || any(userows~=round(userows))
    error('userows must be a valid set of indices into the rows of M')
  end
end

% default for usecols?
if (nargin<3) || isempty(usecols)
  usecols = 1:m;
else
  usecols = unique(usecols);
  if ~isnumeric(usecols) || any(diff(usecols)==0) || ...
      any(usecols<1) || any(usecols>m) || any(usecols~=round(usecols))
    error('usecols must be a valid set of indices into the columns of M')
  end
end

% restrict the search
Muse = M(userows,usecols);

% The maximum down the rows
[maxrows,rowind] = max(Muse,[],1);

% find the best of these maxima
% across the columns
[maxel,colind] = max(maxrows,[],2);
rowind = rowind(colind);

% package the row and column indices
% together, in terms of the original
% matrix in case there was a restiction.
IJ = [userows(rowind),usecols(colind)];
