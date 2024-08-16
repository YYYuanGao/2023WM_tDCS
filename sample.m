% ovals

color_temp = Shuffle(Param.Stimuli.Color_used);
location_temp = randperm(length(Param.Stimuli.Location_used), 7);
Screen('FillOval', wnd, Param.Stimuli.Color(color_temp(1),:), [Param.Stimuli.Locations(location_temp(1),1)-Param.Stimuli.OvalSize/2, Param.Stimuli.Locations(location_temp(1),2)-Param.Stimuli.OvalSize/2, Param.Stimuli.Locations(location_temp(1),1)+Param.Stimuli.OvalSize/2, Param.Stimuli.Locations(location_temp(1),2)+Param.Stimuli.OvalSize/2]);
Screen('FillOval', wnd, Param.Stimuli.Color(color_temp(2),:), [Param.Stimuli.Locations(location_temp(2),1)-Param.Stimuli.OvalSize/2, Param.Stimuli.Locations(location_temp(2),2)-Param.Stimuli.OvalSize/2, Param.Stimuli.Locations(location_temp(2),1)+Param.Stimuli.OvalSize/2, Param.Stimuli.Locations(location_temp(2),2)+Param.Stimuli.OvalSize/2]);
Screen('FillOval', wnd, Param.Stimuli.Color(color_temp(3),:), [Param.Stimuli.Locations(location_temp(3),1)-Param.Stimuli.OvalSize/2, Param.Stimuli.Locations(location_temp(3),2)-Param.Stimuli.OvalSize/2, Param.Stimuli.Locations(location_temp(3),1)+Param.Stimuli.OvalSize/2, Param.Stimuli.Locations(location_temp(3),2)+Param.Stimuli.OvalSize/2]);

% Param.Stimuli.Locations2(:,1)   = (cos(Param.Stimuli.Theta+Param.Stimuli.AngleJitter)*Param.Stimuli.Eccentricity*Param.Settings.PixelPerDegree+Param.Settings.ScrnResolution(3)/2)+ Param.Stimuli.RadiusJitter;
% Param.Stimuli.Locations2(:,2)   = (sin(Param.Stimuli.Theta+Param.Stimuli.AngleJitter)*Param.Stimuli.Eccentricity*Param.Settings.PixelPerDegree+Param.Settings.ScrnResolution(4)/2)+ Param.Stimuli.RadiusJitter;
% Param.Stimuli.Location2_used    = [1 2 3 4 5 6 7 8];
% location_temp = randperm(length(Param.Stimuli.Location2_used), 3);
% Screen('FillOval', wnd, Param.Stimuli.Color(color_temp(1),:), [Param.Stimuli.Locations2(location_temp(1),1)-Param.Stimuli.OvalSize/2, Param.Stimuli.Locations2(location_temp(1),2)-Param.Stimuli.OvalSize/2, Param.Stimuli.Locations2(location_temp(1),1)+Param.Stimuli.OvalSize/2, Param.Stimuli.Locations2(location_temp(1),2)+Param.Stimuli.OvalSize/2]);
% Screen('FillOval', wnd, Param.Stimuli.Color(color_temp(2),:), [Param.Stimuli.Locations2(location_temp(2),1)-Param.Stimuli.OvalSize/2, Param.Stimuli.Locations2(location_temp(2),2)-Param.Stimuli.OvalSize/2, Param.Stimuli.Locations2(location_temp(2),1)+Param.Stimuli.OvalSize/2, Param.Stimuli.Locations2(location_temp(2),2)+Param.Stimuli.OvalSize/2]);
% Screen('FillOval', wnd, Param.Stimuli.Color(color_temp(3),:), [Param.Stimuli.Locations2(location_temp(3),1)-Param.Stimuli.OvalSize/2, Param.Stimuli.Locations2(location_temp(3),2)-Param.Stimuli.OvalSize/2, Param.Stimuli.Locations2(location_temp(3),1)+Param.Stimuli.OvalSize/2, Param.Stimuli.Locations2(location_temp(3),2)+Param.Stimuli.OvalSize/2]);

% fixation
% Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);

for tempi = 1:36
    Screen('DrawDots', wnd, Param.Stimuli.Locations(tempi,:), 3, white);
end

vbl = Screen('Flip',wnd);

%% close all
is_true = 0;
while ~is_true
    [~,~,keyCode] = KbCheck;
    if keyCode(Param.Keys.Esc)
        is_true = 1;
    end
end

reset_test_gamma;
ShowCursor;
Screen('CloseAll');

delete *.asv

