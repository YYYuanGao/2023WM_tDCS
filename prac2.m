% practice.
% By @Yuan 20231030 08:39
%% check parameters used
if exist([CurrDir '\Results\prac\' SubjID '\' SubjID '_results_sess' num2str(sess_num) '_run' num2str(run_num) '.mat'],'file')
    disp(' ');
    disp([SubjID '_run' num2str(run_num) ' has been test, please enter a new run number!']);
    disp(' ');
    abort;
end

%%
trial_index = randperm(Param.Trial.Practice);
title = {'trial_num' 'loc1' 'loc2' 'loc3' 'loc4' 'loc5' 'loc6' 'loc7' 'set_size' ...
    'ovalcentre_x1' 'ovalcentre_x2' 'ovalcentre_x3' 'ovalcentre_x4' 'ovalcentre_x5' 'ovalcentre_x6' 'ovalcentre_x7'...
    'ovalcentre_y1' 'ovalcentre_y2' 'ovalcentre_y3' 'ovalcentre_y4' 'ovalcentre_y5' 'ovalcentre_y6' 'ovalcentre_y7'...
    'ori1' 'ori2' 'ori3' 'ori4' 'ori5' 'ori6' 'ori7' 'color1' 'color2' 'color3' 'color4' 'color5' 'color6' 'color7' ...
    'color_test' 'reported_x' 'reported_y' 'ori_report' 'error' 'trial_onset' 'fix_onset_delay' 'ITI' 'stim_dur' 'delay_dur'  ...
    'acc' 'trial_dur' 'RT' 'delta_x' 'delta_y'};

n_title = numel(title);
n_trial = Param.Trial.Practice;
results = zeros(n_trial,n_title);
results = array2table(results, 'VariableNames', title);

% start
DrawFormattedText(wnd,'Press space to start!','center','center', black);
Screen('Flip',wnd);
is_true = 0;
while (is_true == 0)
    [ifkey,RT_time,keyCode] = KbCheck;
    if keyCode(Param.Keys.Space)
        is_true = 1;
    elseif keyCode(Param.Keys.Esc)
        abort;
    end
end

%% Go!
for trial_i = 1:Param.Trial.Practice
    results.trial_num(trial_i) = trial_i;
    % results.set_size(trial_i) = mod(trial_index(trial_i),4)+3; %3,4,5,6
    % results.set_size(trial_i) = mod(trial_index(trial_i),5)+3; %3,4,5,6,7; 15
    
    remainder = mod(trial_index(trial_i),4)+3;
    if remainder == 3
        results.set_size(trial_i) = 1;
    elseif remainder == 4
        results.set_size(trial_i) = 3;
    elseif remainder == 5
        results.set_size(trial_i) = 5;
    elseif remainder == 6
        results.set_size(trial_i) = 7;
    end

    Theta = 5:10:355;
    Theta([9,18,27,36]) = [];
    if results.set_size(trial_i) == 1
        b = Theta(randperm(length(Theta), 7));
%     elseif results.set_size(trial_i) == 3 %120
%         b = Theta(randperm(length(Theta), 7));
%         while min(diff(sort(b))) < 100
%             b = Theta(randperm(length(Theta), 7));
%         end
        % elseif results.set_size(trial_i) == 4 %90
        %     b = Theta(randperm(length(Theta), 7));
        %     while min(diff(sort(b))) < 70
        %         b = Theta(randperm(length(Theta), 7));
        %     end
%     elseif results.set_size(trial_i) == 5 %70
%         b = Theta(randperm(length(Theta), 7));
%         while min(diff(sort(b))) < 50
%             b = Theta(randperm(length(Theta), 7));
%         end
        % elseif results.set_size(trial_i) == 6 %60
        %     b = Theta(randperm(length(Theta), 7));
        %     while min(diff(sort(b))) < 40
        %         b = Theta(randperm(length(Theta), 7));
        %     end
        %
    else 
        b = Theta(randperm(length(Theta), Param.Stimuli.LocationsNum));
        while abs(min(diff(sort(b)))) < 30
            b = Theta(randperm(length(Theta), Param.Stimuli.LocationsNum));
        end
    end

    % loc x/y
    temp1 = randperm(length(Param.Stimuli.Location_used));
    for loc_i = 1:Param.Stimuli.LocationsNum
        % relative to the centre
        Locations_temp(:,1)   = cosd(b) * Param.Stimuli.Eccentricity * Param.Settings.PixelPerDegree + Param.Settings.ScrnResolution(3)/2;
        Locations_temp(:,2)   = sind(b) * Param.Stimuli.Eccentricity * Param.Settings.PixelPerDegree + Param.Settings.ScrnResolution(4)/2;

        for ori_i = 1: length(b)
            results.(['ori' num2str(ori_i)])(trial_i) = b(ori_i);
        end

        results.(['loc' num2str(loc_i)])(trial_i) = Param.Stimuli.Location_used(temp1(loc_i));
        results.(['ovalcentre_x' num2str(loc_i)])(trial_i) = ((Locations_temp(results.(['loc' num2str(loc_i)])(trial_i), 1)  - Param.Settings.ScrnResolution(3)/2)) / Param.Settings.PixelPerDegree;
        results.(['ovalcentre_y' num2str(loc_i)])(trial_i) = ((Locations_temp(results.(['loc' num2str(loc_i)])(trial_i), 2)  - Param.Settings.ScrnResolution(4)/2)) / Param.Settings.PixelPerDegree;
    end

    % target is shown at Loc1
    sti_location_alltemp = reshape(table2array(results(trial_i, 10:23)), [],2);  %unit degree
    sti_location_alltemp_pixel = sti_location_alltemp * Param.Settings.PixelPerDegree + [Param.Settings.ScrnResolution(3)/2,Param.Settings.ScrnResolution(4)/2];
    tar_location_temp    = [results.ovalcentre_x1(trial_i),results.ovalcentre_y1(trial_i)];

    % define color
    color_temp1 = randperm(length(Param.Stimuli.Color_used),Param.Stimuli.LocationsNum);
    for color_i = 1:Param.Stimuli.LocationsNum
        results.(['color' num2str(color_i)])(trial_i) = color_temp1(color_i);
    end
    results.trial_onset(trial_i) = GetSecs;

    %% ITI
    %     Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);
    vbl = Screen('Flip',wnd);
    results.fix_onset_delay(trial_i) = vbl-results.trial_onset(trial_i);

    %% draw stim
    for oval_i = 1:results.set_size(trial_i)
        idx = find(b==results.(['ori' num2str(oval_i)])(trial_i));
        Screen('FillOval', wnd, Param.Stimuli.Color(results.(['color', num2str(oval_i)])(trial_i), :), [sti_location_alltemp_pixel(idx,1)-Param.Stimuli.OvalSize/2, sti_location_alltemp_pixel(idx,2)-Param.Stimuli.OvalSize/2, sti_location_alltemp_pixel(idx,1)+Param.Stimuli.OvalSize/2, sti_location_alltemp_pixel(idx,2)+Param.Stimuli.OvalSize/2]);
    end

    % for oval_i = 1:results.set_size(trial_i)
    %     Screen('BlendFunction', wnd, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    %     Screen('FillOval', wnd, Param.Stimuli.Color(results.(['color', num2str(oval_i)])(trial_i), :), [sti_location_alltemp_pixel(oval_i,1)-Param.Stimuli.OvalSize/2, sti_location_alltemp_pixel(oval_i,2)-Param.Stimuli.OvalSize/2, sti_location_alltemp_pixel(oval_i,1)+Param.Stimuli.OvalSize/2, sti_location_alltemp_pixel(oval_i,2)+Param.Stimuli.OvalSize/2]);
    % end

    % for tempi = 1:36
    %     Screen('DrawDots', wnd, Param.Stimuli.Locations(tempi,:), 3, white);
    % end

    % fixation
    %     Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);
    vbl = Screen('Flip',wnd,vbl+Param.Trial.ITI); %-Slack
    results.ITI(trial_i) = vbl-results.trial_onset(trial_i)-results.fix_onset_delay(trial_i);
    %% delay
    rhyme = randperm(length(Param.Trial.Delay),1);
    results.delay_dur(trial_i) = Param.Trial.Delay(rhyme);

    %     Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);
    %     vbl = Screen('Flip',wnd,vbl+Param.Trial.MaskDura);     %-Slack
    %     results.mask_dur(trial_i) = vbl-results.trial_onset(trial_i)-results.fix_onset_delay(trial_i)-results.ITI(trial_i)-results.stim_dur(trial_i);
    vbl = Screen('Flip',wnd,vbl+Param.Trial.StiDura);     %-Slack
    results.stim_dur(trial_i) = vbl-results.trial_onset(trial_i)-results.fix_onset_delay(trial_i)-results.ITI(trial_i);

    %% report
    results.color_test(trial_i) = results.color1(trial_i);
    Screen('BlendFunction', wnd, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    Screen(wnd,'FillOval', Param.Stimuli.Color(results.color_test(trial_i),:), [Param.Settings.ScrnResolution(3)/2-Param.Stimuli.OvalSize/2,Param.Settings.ScrnResolution(4)/2-Param.Stimuli.OvalSize/2,Param.Settings.ScrnResolution(3)/2+Param.Stimuli.OvalSize/2,Param.Settings.ScrnResolution(4)/2+Param.Stimuli.OvalSize/2]);
    %     Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);
%     num = results.color1(trial_i);
%     
%     results.ori_test(trial_i) = results.(['ori' num2str(num)])(trial_i);
% 
%     
    vbl = Screen('Flip',wnd,vbl+results.delay_dur(trial_i));%-Slack
    stimulus_onset = vbl;
    results.delay_dur(trial_i) = vbl-results.trial_onset(trial_i)-results.fix_onset_delay(trial_i)-results.ITI(trial_i)-results.stim_dur(trial_i);

    %% response
    %     ShowCursor('Arrow',wnd);

    WaitSetMouse(Param.Settings.ScrnResolution(3)/2, Param.Settings.ScrnResolution(4)/2,wnd); % set mouse at the center
    is_true = 0;
    while (is_true == 0)
        [mousex,mousey,buttons]=GetMouse(wnd);
        if any(buttons) && buttons(1)
            is_true = 1;
            results.reported_x(trial_i) = (mousex-Param.Settings.ScrnResolution(3)/2)/Param.Settings.PixelPerDegree;
            results.reported_y(trial_i) = (mousey-Param.Settings.ScrnResolution(4)/2)/Param.Settings.PixelPerDegree;

            rx = results.reported_x(trial_i);
            ry = results.reported_y(trial_i);
            rangle = cart2pol(rx, ry);
            results.ori_report(trial_i) = 90-(rad2deg(rangle));
%             results.ori_report(trial_i) =  90 - results.ori_report(trial_i);

            ovalcentre_x = tar_location_temp(1);
            ovalcentre_y = tar_location_temp(2);
            rangle_test = cart2pol(ovalcentre_x, ovalcentre_y);
            results.ori_test(trial_i) = 90-(rad2deg(rangle_test));
            
            results.delta_x(trial_i) = results.reported_x(trial_i) - ovalcentre_x;
            results.delta_y(trial_i) = results.reported_y(trial_i) - ovalcentre_y;
%             results.error(trial_i) = results.ori_report(trial_i)-results.ori_test(trial_i);
            results.error(trial_i) = anglediff(results.ori_report(trial_i),results.ori_test(trial_i));
            error = results.error(trial_i);
            DrawFormattedText(wnd,double(['error : ' num2str(error)]),'center','center',black);
            Screen('Flip',wnd);
            WaitSecs(1);

            RT_time = GetSecs;
            results.RT(trial_i) = RT_time - stimulus_onset;

        elseif buttons(2) % middle mouse
            abort;
        end    

        if ~is_true
%             for oval_i = 1:results.set_size(trial_i)
%                 idx = find(b==results.(['ori' num2str(oval_i)])(trial_i));
%                 Screen('FillOval', wnd, Param.Stimuli.Color(results.(['color', num2str(oval_i)])(trial_i), :), [sti_location_alltemp_pixel(idx,1)-Param.Stimuli.OvalSize/2, sti_location_alltemp_pixel(idx,2)-Param.Stimuli.OvalSize/2, sti_location_alltemp_pixel(idx,1)+Param.Stimuli.OvalSize/2, sti_location_alltemp_pixel(idx,2)+Param.Stimuli.OvalSize/2]);
%             end
%             Screen('Flip',wnd);

%             Screen('BlendFunction', wnd, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
            Screen('FrameOval',wnd,white,[Param.Settings.ScrnResolution(3)/2-Param.Stimuli.WheelOuter/2,Param.Settings.ScrnResolution(4)/2-Param.Stimuli.WheelOuter/2,Param.Settings.ScrnResolution(3)/2+Param.Stimuli.WheelOuter/2,Param.Settings.ScrnResolution(4)/2+Param.Stimuli.WheelOuter/2]);
            Screen('FrameOval', wnd, white, [Param.Settings.ScrnResolution(3)/2-Param.Stimuli.WheelInner/2,Param.Settings.ScrnResolution(4)/2-Param.Stimuli.WheelInner/2,Param.Settings.ScrnResolution(3)/2+Param.Stimuli.WheelInner/2,Param.Settings.ScrnResolution(4)/2+Param.Stimuli.WheelInner/2]);
            Screen(wnd,'FillOval', Param.Stimuli.Color(results.color_test(trial_i),:), [mousex-Param.Stimuli.OvalSize/2, mousey-Param.Stimuli.OvalSize/2, mousex+Param.Stimuli.OvalSize/2, mousey+Param.Stimuli.OvalSize/2]);
            Screen('Flip',wnd);
        end
    end

    results.trial_dur(trial_i) = GetSecs - results.trial_onset(trial_i);
    if results.error(trial_i) < min(diff(sort(b)))/2
        results.acc(trial_i) = 1;
    else
        results.acc(trial_i) = 0;
    end
end
WaitSecs(0.2);

    %% save the data
    resultsDir = [CurrDir '\Results\prac\' SubjID '\'];
    if ~isdir(resultsDir)
        mkdir(resultsDir);
    end
    
    cd(resultsDir);
    results_name = [SubjID '_results_sess' num2str(sess_num) '_run' num2str(run_num) '.mat'];
    save(results_name,'results','Param');
    cd(CurrDir);
    
    %% plot
    figure(1)
    bin1 =-180:45:180;
    subplot(2,2,1);
    hist(results.error(results.set_size==1),bin1);
    m1 = mean(abs(results.error(results.set_size==1)))
    subplot(2,2,2);
    hist(results.error(results.set_size==3),bin1);
    m2 = mean(abs(results.error(results.set_size==3)))
    subplot(2,2,3);
    hist(results.error(results.set_size==5),bin1);
    m3 = mean(abs(results.error(results.set_size==5)))
    subplot(2,2,4);
    hist(results.error(results.set_size==7),bin1);
    m4 = mean(abs(results.error(results.set_size==7)))

    figure(2)
    bin2 = 0:1:5;
    subplot(2,2,1);
    hist(results.RT(results.set_size==1),bin2);
    subplot(2,2,2);
    hist(results.RT(results.set_size==3),bin2);
    subplot(2,2,3);
    hist(results.RT(results.set_size==5),bin2);
    subplot(2,2,4);
    hist(results.RT(results.set_size==7),bin2);

    figure(3)
    plot([1,3,5,7],[m1 m2 m3 m4]);    
    %%
    reset_test_gamma;
    ShowCursor;
    Screen('CloseAll');
    delete *.asv