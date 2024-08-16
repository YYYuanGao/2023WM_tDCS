% BY Yuan Gao: 20230418 11:20
% list all the parameters used in this experiment

Param = struct;

%% Screen Settings
Param.Settings.ViewDistance    = 670;              % 1100 mm 
Param.Settings.ScrnResolution  = [0 0 1920 1080]%[0 0 1920 1080];   % rect_computer = [0 0 40 30];
Param.Settings.SquareSize      = Param.Settings.ScrnResolution(3);             % 400 mm 
Param.Settings.SquareLength    = 540;              % 154 mm,532mm  
Param.Settings.PixelPerDegree  = 2*Param.Settings.ViewDistance*tan(1/2/180*pi)*Param.Settings.SquareSize/Param.Settings.SquareLength;       

%% Keys for response
Param.Keys.Space      = 32;  
Param.Keys.Esc        = 27;
Param.Keys.Left       = 37;  
Param.Keys.Right      = 39;
Param.Keys.Up         = 38;
Param.Keys.Down       = 40;
Param.Keys.Trigger1   = 83;  % 's'       

%% parameters for stimuli
% Locations
Param.Stimuli.Eccentricity     = 7.5;    % degree
Param.Stimuli.Theta            = (0:35)*2*pi/36-pi/36; % [10:10:360]
Param.Stimuli.Locations(:,1)   = cos(Param.Stimuli.Theta) * Param.Stimuli.Eccentricity * Param.Settings.PixelPerDegree + Param.Settings.ScrnResolution(3)/2;
Param.Stimuli.Locations(:,2)   = sin(Param.Stimuli.Theta) * Param.Stimuli.Eccentricity * Param.Settings.PixelPerDegree + Param.Settings.ScrnResolution(4)/2;
Param.Stimuli.Location_used    = [1 2 3 4 5 6 7];
Param.Stimuli.LocationsNum     = 7; 

%% parameters for fixation
Param.Fixation.CrossColor    = [1,1,1]*255;
Param.Fixation.CrossSize     = 0.3*Param.Settings.PixelPerDegree;
Param.Fixation.CrossWidth    = 0.1*Param.Settings.PixelPerDegree;
% +
Param.Fixation.CrossLoc      = [Param.Settings.ScrnResolution(3)/2-Param.Fixation.CrossSize, Param.Settings.ScrnResolution(3)/2+Param.Fixation.CrossSize, Param.Settings.ScrnResolution(3)/2, Param.Settings.ScrnResolution(3)/2;...
                                     Param.Settings.ScrnResolution(4)/2, Param.Settings.ScrnResolution(4)/2, Param.Settings.ScrnResolution(4)/2-Param.Fixation.CrossSize, Param.Settings.ScrnResolution(4)/2+Param.Fixation.CrossSize] + [offset',offset',offset',offset'];
% X
Param.Fixation.CrossLoc2     = [Param.Settings.ScrnResolution(3)/2-Param.Fixation.CrossSize*sqrt(2)/2, Param.Settings.ScrnResolution(3)/2+Param.Fixation.CrossSize*sqrt(2)/2, Param.Settings.ScrnResolution(3)/2-Param.Fixation.CrossSize*sqrt(2)/2, Param.Settings.ScrnResolution(3)/2+Param.Fixation.CrossSize*sqrt(2)/2;...
                                     Param.Settings.ScrnResolution(4)/2-Param.Fixation.CrossSize*sqrt(2)/2, Param.Settings.ScrnResolution(4)/2+Param.Fixation.CrossSize*sqrt(2)/2, Param.Settings.ScrnResolution(4)/2+Param.Fixation.CrossSize*sqrt(2)/2, Param.Settings.ScrnResolution(4)/2-Param.Fixation.CrossSize*sqrt(2)/2]+[offset',offset',offset',offset'];  

Param.Fixation.OvalSize      = 0.6*Param.Settings.PixelPerDegree;
Param.Fixation.OvalColor     = [0,0,0]; 
Param.Fixation.Offset        = 0.5*Param.Settings.PixelPerDegree; 
Param.Fixation.FontColor     = [1,1,1]*255;

%% parameters for trials
Param.Trial.ITI              = 0.5;     
Param.Trial.StiDura          = 0.3; 
Param.Trial.MaskDura         = 0.3;
Param.Trial.Delay            = [3 4 5];
Param.Trial.Feedback         = 0.3;

Param.Trial.TestNum          = 80;  
Param.Trial.Practice         = 12;
Param.Trial.minirun          = 40;
%% masks
Param.Stimuli.OuterRadius      = 0.3;
Param.Stimuli.OuterSize        = round(Param.Stimuli.OuterRadius*Param.Settings.PixelPerDegree);   % radius pixels

Param.Stimuli.GratingOri       = [15 75 135];
Param.Stimuli.GratingContrast  = 0.5;
Param.Stimuli.Spatial_freq     = 1.25/Param.Settings.PixelPerDegree;       % 0.02 % 6 bars in total
%% sets for stimuli 
Param.Stimuli.Color(1,:)     = [0.89,0.09,0.05]*255; % red
Param.Stimuli.Color(2,:)     = [0.98,0,0.81]*255; % pink
Param.Stimuli.Color(3,:)     = [1,1,0]*255; % yellow
Param.Stimuli.Color(4,:)     = [0.4940,0.1840,0.5560]*255; % purple
Param.Stimuli.Color(5,:)     = [0.1,0.06,1]*255; % blue
Param.Stimuli.Color(6,:)     = [0,1,0]*255; % green 
Param.Stimuli.Color(7,:)     = [0,1,1]*255; % cyan
% Param.Stimuli.Color(8,:)     = [1,0.38,0]*255; % orange
Param.Stimuli.Color(8,:)     = [0,0,0]*255; 

Param.Stimuli.Color_used     = [1 2 3 4 5 6 7 8];
Param.Stimuli.OvalSize       = 0.2 * 2 * Param.Settings.PixelPerDegree; % diameter
%% sets for wheel
Param.Stimuli.WheelOuter     = (Param.Stimuli.Eccentricity+0.2) * 2 * Param.Settings.PixelPerDegree; %diameter
Param.Stimuli.WheelInner     = (Param.Stimuli.Eccentricity-0.2) * 2 * Param.Settings.PixelPerDegree; %diameter
Param.Stimuli.Width          = 0.1 * Param.Settings.PixelPerDegree;
Param.Stimuli.Height         = 0.1 * Param.Settings.PixelPerDegree;

%% open window 
screens = Screen('Screens');
screenNumber = max(screens);	
Screen('Preference', 'SkipSyncTests', 1);

white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);
gray = (white+black)/2;
% gray = [70 70 70];

wnd = Screen('OpenWindow',screenNumber,gray,Param.Settings.ScrnResolution); 
Screen('BlendFunction',wnd, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
Screen('TextFont',wnd,'Arial');
Screen('TextSize',wnd,30);

RefreshDur = Screen('GetFlipInterval',wnd);
RefreshRate = 1./RefreshDur;
Slack = RefreshDur/2;

