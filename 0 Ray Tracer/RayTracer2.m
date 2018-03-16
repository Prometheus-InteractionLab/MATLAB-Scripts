%% Raytracerv2 Program
%by Steve Busch, Sandia National Laboratories
% Modified, Jan 2016 by Ethan Eagle
%
%This program is used to perform raytracing for a user-specified geometry
%and user specified ray fans. The first part of the program reads in the
%desired surface geometry from an excel file, and creates an array of
%corresponding surface objects. The ray fans originate from a user-
%specified plane. The piston position as a function of crank angle is
%loaded from an external file, so that the distance between the piston and
%the ray-origin plane is known for a given crank angle. 
%
%Then, parameters are defined for the simulation, such as which rays are to
%be traced, how many rays are to be traced, which crank angles are to be 
%simulated, and if the output is to be plotted as the ray tracing
%progresses. 
%
%The ray tracing is performed for each crank angle and radial position, and
%uses the array of surfaces, the ray fan for the current crank angle and
%radius, and a boolean specifying whether or not the output should be
%plotted. For each crank angle and radius, the results are saved in a .mat
%file that can be viewed with ShowRayTracerResults.m

clear planes rays surfaces
close all
pause(0.1)

core = 42; %This parameter is used with multiple instances of this program
          %to perform multiple simulations in parallel without using
          %parallel computing tools. Do a search for "core" in this file to
          %see how it is used (there's only one place).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Section 1: Define surfaces in environment%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bowl = 'Contour';
sheet = 'Contoured'; % 'Contoured' for contour HD
switch sheet
    case 'Flat'
        Top = '$A$14:$B$102';
        Side='$C$14:$D$86';
        Bottom = '$E$14:$F$83';
    case 'Contoured'
        Top = '$A$17:$B$136';
        Side = '$C$17:$D$145';
        Bottom = '$E$17:$F$130';
end

GMfilename=['C:\Users\sbusch\Documents\Data\Software\RayTracer\',...
            'Bowl form\FullBowlProfile_FP48_B2.xlsx'];
Fordfilename=['C:\Users\sbusch\Documents\Data\Software\RayTracer\',...
            'Bowl form\FordBowl\RZcartesianCord.xlsx'];
testfilename=['\\906-165-press\e$\Data\matlabs\0 Ray Tracer\test1.xlsx'];
contourfilename= ['\\906-165-press\e$\Data\matlabs\0 Ray Tracer\PistonContour.xlsx'];
switch bowl
    case 'Contour'
        dt=0.001; %Spline resolution
        %Get bowl coordinates
        disp('Reading bowl coordinates...');
        pause(0.01)
        [numeric, txt, raw] = xlsread(contourfilename,sheet,Top); %#ok<*NASGU>
        x=numeric(:,1);
        y=numeric(:,2);
        x1 = [-flipud(x);x]
        y1 = [flipud(y);y];
        %Create bowl surface using coordinates from file, 3rd order B-splines,
        %step size of 0.0001, index of refraction on one side 1, other side 1.4
        surfaces{1}=createSurface(x1,y1,3,dt,'refracting',1,1.46);
        %surfaces{end+1}=createSurface(-x,y,3,dt,'refracting',1,1.46);%mirror to make the whole profile
    case 'test'
        dt=0.001; %Spline resolution
        %Get bowl coordinates
        disp('Reading bowl coordinates...');
        pause(0.01)
        [numeric, txt, raw] = xlsread(testfilename,sheet,Top); %#ok<*NASGU>
        x=numeric(:,1);
        x1=[x;-flipud(x);x(end)]; %mirror to make the whole profile
        y=numeric(:,2);
        y1=[y;flipud(y);y(end)]; %mirror to make the whole profile
        %Create bowl surface using coordinates from file, 3rd order B-splines,
        %step size of 0.0001, index of refraction on one side 1, other side 1.4
        surfaces{1}=createSurface(x1,y1,3,dt,'refracting',1,1.46);
    case 'GM'
        dt=0.001; %Spline resolution
        %Get bowl coordinates
        disp('Reading bowl coordinates...');
        pause(0.01)
        [numeric, txt, raw] = xlsread(GMfilename,'Full Piston',...
            '$A$21:$B$410'); %#ok<*NASGU>
        x=numeric(:,1);
        x=[x;-flipud(x)]; %mirror to make the whole profile
        y=numeric(:,2);
        y=[y;flipud(y)]; %mirror to make the whole profile
        %Create bowl surface using coordinates from file, 3rd order B-splines,
        %step size of 0.0001, index of refraction on one side 1, other side 1.4
        surfaces{1}=createSurface(x,y,3,dt,'refracting',1,1.46);
    case 'Ford'
        dt=0.001; %Spline resolution
        %Get bowl coordinates
        disp('Reading bowl coordinates...');
        pause(0.01)
        [numeric, txt, raw] = xlsread(Fordfilename,'RZcartesianCord',...
            '$H$3:$I$1846'); %#ok<*ASGLU>
        x=numeric(:,1);
        x=[-flipud(x);x;x(end);x(end)]; %mirror to make the whole profile
        y=numeric(:,2);
        y=[flipud(y);y;y(end);y(end)]; %mirror to make the whole profile
        %Create bowl surface using coordinates from file, 3rd order B-splines,
        %step size of 0.0001, index of refraction on one side 1, other side 1.4
        surfaces{1}=createSurface(x,y,3,dt,'refracting',1,1.46);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get coordinates for the side of the piston
disp('Reading side coordinates...')
[numeric, txt, raw] = xlsread(contourfilename,sheet,Side);
x=numeric(:,1);
y=numeric(:,2);
surfaces{end+1}=createSurface(x,y,3,2*dt,'absorbing',1.46,1);
% surfaces{end+1}=createSurface(x+0.05,y,3,0.0001,'absorbing',1.46,1);
surfaces{end+1}=createSurface(-x,y,3,2*dt,'absorbing',1.46,1);
% surfaces{end+1}=createSurface(-x-0.05,y,3,0.0001,'absorbing',1.46,1);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get bottom surface of piston
disp('Reading bottom coordinates...')
[numeric, txt, raw] = xlsread(contourfilename,sheet,Bottom);
x=numeric(:,1);
x=[-flipud(x);x]; %mirror to make the whole profile
y=numeric(:,2);
y=[flipud(y);y]; %mirror to make the whole profile
surfaces{end+1}=createSurface(x,y,3,dt,'refracting',1.46,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make absorbing cylinder to simulate mirror
disp('Generating Bowl Rim...')
MirrorD=1.96 * 25.4;  %cylinder diameter in mm
BowlD = 138.58;
BowlRimD = 99.5;
x = (46.2:1:BowlD/2)';
y = -10.386*ones(length(x),1);
surfaces{end+1}=createSurface(x,y,3,0.02,'absorbing',1,1);
surfaces{end+1}=createSurface(-x,y,3,0.02,'absorbing',1,1);

disp('Simulating turning mirror and absorbing sidewalls...')
xMirror= (-MirrorD/2:-1:-BowlD/2)';
yMirror= (-310*ones(length(xMirror),1));
ySideWall=(-310:2:90)';
xSideWall=-BowlD/2*ones(length(ySideWall),1);
% surfaces{end+1}=createSurface(x,y,1,1,'absorbing',1,1);
% surfaces{end+1}=createSurface(-x,y,1,1,'absorbing',1,1);
%xnew=-(D/2:0.1:30)';
%x=[x;xnew];
%y=[y;(-60*ones(size(xnew)))];

%x=[x;(-30*ones(11,1))];
%y=[y;linspace(-60,-50,11)'];
surfaces{end+1}=createSurface([xMirror;xSideWall],[yMirror;ySideWall],3,0.02,'absorbing',1,1);
surfaces{end+1}=createSurface(-[xMirror;xSideWall],[yMirror;ySideWall],3,0.02,'absorbing',1,1);
%%
disp('Simulating absorbing cylinder walls...')
surfaces{end+1}=createSurface(x,y,1,1,'absorbing',1,1);
surfaces{end+1}=createSurface(-x,y,1,1,'absorbing',1,1);

disp('Surfaces splined.')



%Array of ray fans above piston CA=360
CAD=(0:0.1:719.9)'; %Crank angle vector with 0.1° resolution
CArad = CAD.*pi./180.;
SquishHeight=2.54; %mm
Stroke = 152.4; %mm
PistonP = 152.4*(1-cos(CArad));
%PistonP is s, the distance between the top of the piston and the firedeck
PistonP=PistonP+SquishHeight; %Function of crank angle


%% break%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Section 2: Define Rays and set up simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%Define the ray fans as having a linear or random distribution of angles
mode='linear';
%Define the path in which the simulation results are to be saved. This
%should be a folder in which you can store a few GB.
pathname2save=['\\906-165-press\e$\Data\matlabs\0 Ray Tracer\',...
        'Contour'];
    if ~exist(pathname2save)
        mkdir(pathname2save)
    end
%Define the type of ray tracing as nonsequential
type='nonsequential';
%Which crank angles are to be simulated? This determines the vertical
%position of the origin of each ray fan. Specify how many rays are in each
%fan, and if the output should be plotted to make sure the simulation is
%running properly (you probably want to do this with ~10-100 rays first).
%Define a startX vector that indicates which radial positions will be
%simulated for each crank angle.
%Finally, specify the start and end angle for the ray fan in degrees (0°
%is a horizontal line pointing right)
% radialDistance=[0:0.5:22.6,23:0.5:41];
radialDistance=(9.5:0.5:45);
% radialDistance=[0:0.5:22.64,22.64,22.75:0.5:41];
thetaInj=15.5*pi/180; %down angle of the injected spray
switch core
    case 1
        CA_imaged=[355 356 357 358 359];
        numRaysPerFan=3000;
        showplot=false;
        angle1=-60;
        angle2=-120;
    case 2
        CA_imaged=[360 361 362 363 364];
        numRaysPerFan=3000;
        showplot=false;
        angle1=-60;
        angle2=-120;
    case 3
        CA_imaged=[365 366 367 368 369];
        numRaysPerFan=3000;
        showplot=false;
        angle1=-60;
        angle2=-120;
    case 4
        CA_imaged=[370 371 372 373 374];
        numRaysPerFan=3000;
        showplot=false;
        angle1=-60;
        angle2=-120;
    case 5
       CA_imaged=[375 376 377 378 379];
        numRaysPerFan=3000;
        showplot=false;
        angle1=-60;
        angle2=-120;
    case 42
        CA_imaged = [360];
        numRaysPerFan = 30;
        showplot = false;
        angle1 = -60;
        angle2 = -120;
    otherwise
        return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Section 3: Perform Raytracing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Loop through for all crank angles (startY values). For each startY, make a
%ray fan at each x position (this is the radial position defined above).
%Run the raytracing algorithm for each ray fan.

tic
for idxCA=1:length(CA_imaged) %for each CA
    
    
    startY=zeros(numel(radialDistance),1);
    hRays=zeros(numel(radialDistance),1); 
     for rIdx=1:numel(radialDistance) 
         
    hRays(rIdx)=1.2+radialDistance(rIdx)*tan(thetaInj); %in mm;
    startY(rIdx)=interp1(CAD,PistonP-hRays(rIdx),...
                     CA_imaged(idxCA),'linear'); %Along jet axes in bowl
    startX=radialDistance; %the radial coordinate
    %startY=

%         if CA_imaged(idxCA) < 6
%             %Along jet axis in bowl, halfway between piston and head in
%             %squish region
%             if radialDistance(rIdx) <= 49
%                 %Define vertical distance between origin of rays and the firedeck
%                
%                 %For the crank angles to be simulated, calculate the distance 
%                 %between the piston and the origins of the ray fans
%                 startY(rIdx)=interp1(CAD,PistonP-hRays(rIdx),...
%                     CA_imaged(idxCA),'linear'); %Along jet axes in bowl
%             else
%                 startY(rIdx)=interp1(CAD,PistonP,...
%                     CA_imaged(idxCA),'linear')/2; %h/2 above squish region
%             end       
% 
%         else
%             if radialDistance(rIdx) <= 49 %constant distance above the cylinder
%                 %1 mm above piston bowl surface
%                 [numeric, txt, raw] = xlsread(testfilename,'Sheet1',...
%                     '$A$14:$B$63'); %#ok<*NASGU>
%                 x=flipud(numeric(:,1));
%                 y=flipud(numeric(:,2))+1; %Follow the bowl profile
%                 [x,iX,temp] = unique(x);
%                 y=y(iX);
%                 try
%                 startY(rIdx)=spline(x,y,radialDistance(rIdx)); 
%                 catch
%                     keyboard
%                 end
%             else
             %end 
        %end
    end
    for idx=1:length(startX) %For each position along the jet axis
        clc
        close all
        clear rays prevrays finalrays surfaceNums intersections
        %Create the ray fan 
        %Average the angles to get the center angle
        centerAngle=(angle2+angle1)/2; %in degrees
        %Take the difference of the angles to get the fanAngle
        fanAngle=abs(angle2-angle1); %in degrees
        %Where do the rays in this bundle originate?
        startPoint=[startX(idx);startY(idx)];
        %Make the ray fan. 100 units of power are divided equally among the
        %rays. At the moment, power is used to determine whether a not a
        %ray is absorbed and makes it to the mirror underneath the piston.
        %A ray with a power of 0 is said to be absorbed.
        rays=createrays(numRaysPerFan,startPoint,centerAngle(1),...
            fanAngle(1),mode,100,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Perform the ray tracing. Sequential tracing means that the order in 
    %which the surfaces are defined is important, because a plane can only 
    %be interacted with once. The advantage is better computational time.
    %Non-sequential tracing means that each ray is traced until it has no
    %more interactions. A ray can interact multiple times with a surface 
    %and internal reflection is possible. It's slower but necessary for our
    %piston geometry.
    
    switch type
        %This is the heart of the ray tracer. It requires the array of
        %surfaces, the array of rays that make up the ray bundle, and an
        %input as to whether or not the rays should be plotted as they are
        %computed. The result is 
        case 'sequential'
            %This function may not be compatible with the rest of this
            %code; it has not been thoroughly tested.
            [power, intersections, finalrays, prevrays]=...
                sequential_Surfacetrace(surfaces,rays);
        case 'nonsequential'
            [power, intersections, finalrays, prevrays, surfaceNums]=...
                nonSequential_Surfacetrace(surfaces,rays,showplot);
    end

    %Plot final rays traced backwards to visualize mapping of the virtual
    %image
    cnt=1;
    usefulrays={};
    for r=1:length(rays)
        %Rays are "useful" if they haven't been absorbed by either the
        %piston skirt or the imaginary cylinder that represents the mirror
        if finalrays{r}.power > 0
            if showplot
                plotpath([intersections{r}(:,end),...
                    finalrays{r}.t_ray(-100)],'k');
            end
            usefulrays{cnt}=finalrays{r};  %#ok<SAGROW>
            cnt=cnt+1;
        end
    end
    %Look to see where the useful rays cross each other.
    %The cloud of crossing points shows us where/how the virtual image of
    %the target point forms (this is how Steve understands it ;-) ).
    [X,Y] = findCrossings(usefulrays);
    
    %Save the data for the current run in a file.
    filename2save=sprintf('Raytrace_CA_%06.2f_h%06.3f_R%06.3f',...
        CA_imaged(idxCA),startY(idx),startX(idx));
    filename2save=strrep(filename2save,'.','p');%this is clever replacement for matlab numbers to a p.
    filename2save=strrep(filename2save,'-','m');%this is clever replacement for matlab numbers to a p.
    filename2save=[filename2save,'.mat']; %#ok<AGROW>
    
    save(fullfile(pathname2save,filename2save),'X','Y',...
        'intersections','finalrays','usefulrays')

        if showplot
            %Not sure how useful this stuff is. It's left over from an old
            %program...
            Xm=mean(X);
            Ym=mean(X);
            R=sqrt((X-Xm).^2+(Y-Ym).^2);
            Rlimit=60;
            Xm=mean(X((R<Rlimit)));
            Ym=mean(Y((R<Rlimit)));
            plot(X(R<Rlimit),Y(R<Rlimit),'LineStyle','none','Color',...
                [0 0.7 0],'Marker','o','MarkerSize',2,'MarkerEdgeColor',...
                [0 0.7 0],'MarkerFaceColor',[0 0.7 0])
            plot([startPoint(1),Xm],[startPoint(2),Ym],'c-+',...
                'LineWidth',2,'MarkerSize',8)
            %Make histogram
            [n,xout]=hist(X(R<Rlimit),82*10);
            % bar(xout,n)
            h = findobj(gca,'Type','patch');
            set(h,'FaceColor','r','EdgeColor','none')
            [maxN, newIdx]=max(n);

            endX(idx)=xout(newIdx);
        end
    end
end
toc
