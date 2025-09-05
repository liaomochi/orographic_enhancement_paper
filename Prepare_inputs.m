% this scripts is for all world basins.
%% [Global paper below] read streamflow records from other mountain range in the world and create events text files


home_alps = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Alps/streamflow_gauges_daily/';
%home_andes = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Andes/streamflow_gauges_daily/';
%home_himalayas = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Himalayas/streamflow_gauges_daily/';

home_andes = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Andes/streamflow_daily_entire_Andes/';
home_himalayas = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Himalayas/streamflow_daily_entire_Himalayas/';

for mountain = [4 5]
    if mountain == 1
        home_moun = home_alps;mountname = 'Alps';
    end
    if mountain == 4
        home_moun = home_andes;mountname = 'Andes';
    end
    if mountain == 5
        home_moun = home_himalayas;mountname = 'Himalayas';
    end
    
    lastall = ls([home_moun '*.txt']);
    c2 = strsplit(lastall);
    
    parfor ii = 1:length(c2)-1
        
        warning('off')
        close all
        ffnm = c2{ii};
        discharge=readtable(ffnm);
        if size(discharge,2)<3
            % no records at all
            continue
        end
        
        discharge_v = table2array(discharge(:,3));
        obsnow=discharge_v;
        if length(discharge_v)<30
            continue
        end
        disp(ii)
        
        
        fid=fopen(ffnm,'r');
        n_lines = 10;
        outlet_info = cell(10,1);
        for jj = 1:15
            outlet_info(jj) = {fgetl(fid)};
        end
        fclose(fid);
        
        gauge_no = outlet_info{9}(26:end);
        river_nm = outlet_info{10}(26:end);
        area_km2 = outlet_info{15}(30:end);
        country_nm = outlet_info{12}(26:27);
        
        
        time_s = table2cell(discharge(:,1));
        record_len = length(time_s);
        plot_step = floor(record_len/12);
        xtk = 1:plot_step:record_len;
        [ystart,mstart,dstart] = ymd(time_s{1});[yend,mend,dend] = ymd(time_s{end});
        xlab = string;
        
        for il = 1:plot_step:record_len
            [ys,ms,ds] = ymd(time_s{il});
            ysshort = mod(ys,100);
            xlab = [xlab;strcat(num2str(ysshort,'%2.2d'),'/',num2str(ms,'%2.2d'))];
        end
        xlab=cellstr(xlab);
        %
        f = figure('visible','off');
        f.Position = [10 10 1300 500];
        plot(discharge_v)
        ylim([0 max(discharge_v)])
        
        xlabel([num2str(ystart) ' - ' num2str(yend) ' YY/MM'])
        title([ country_nm '  ' river_nm '  ' gauge_no '  ' area_km2 'km^2'])
        ylabel('Streamflow (m^3/s)')
        xlim([0 length(time_s)+0.05*record_len])
        set(gca,'XTick',xtk,'XTickLabel',xlab,'FontSize',20,'LineWidth',2.75,'FontWeight','bold')
        Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
        %[Left Bottom Right Top] spacing
        NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3)-0.01 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
        set(gca, 'Position', NewPos);% when use this set the figure visible off for better cropping figure
        saveas(gcf,['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/' mountname '/graphs/' country_nm '-' river_nm '-' gauge_no '-data_days' num2str(record_len)],'png')
        
        % ########################Event selection############################################################
        mkdir(['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/' mountname '/basins_events/'])
        mkdir(['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/' mountname '/basins_events/' gauge_no])
        
        
        ptl = 97.5;
        flowthre = prctile(discharge_v,ptl);% 99% will give 5-10 events per basin
        obsthre = flowthre;% when peak is above this value, it is considered as a flood event
        plot_step = floor(record_len/7);
        xtk = 1:plot_step:record_len;
        xlab = string;
        for il = 1:plot_step:record_len
            %kp=kp+1;
            [ys,ms,ds] = ymd(time_s{il});
            ysshort = mod(ys,100);
            xlab = [xlab;strcat(num2str(ysshort,'%2.2d'),'/',num2str(ms,'%2.2d'))];
        end
        xlab=cellstr(xlab);
        
        f = figure('visible','off');
        %         t = tiledlayout(1,1,'Padding','none');
        %         t.Units = 'inches';
        %         t.OuterPosition = [0.1 0.25 5.65 3.6];
        %         nexttile;
        f.Position = [5 5 700 400];
        sh = scatter(1:length(discharge_v),discharge_v,12,'b.');
        hold on
        plot([1 9999999],[obsthre obsthre],'r-.','LineWidth',2.5)
        hold off
        
        box on
        legend('Obs.',[num2str(obsthre,'%2.1f') 'm^3/s'],'interflow','± 5% obs.','Location','Northeast')
        ylabel('Streamflow (m^3/s)')
        %legend boxoff
        %yrf = prctile(obsnow,99.99);
        
        ylim([0 ceil(max(discharge_v)/30)*30])
        xlim([0 length(obsnow)])
        %title(['Basin' num2str(id,'%2.2d') ' Max=' num2str(obsmax,'%2.1f') 'm^3/s'])
        title([ country_nm '  ' river_nm '  ' gauge_no '  ' area_km2 'km^2'])
        xlabel([num2str(ystart) ' - ' num2str(yend) ' YY/MM'])
        set(gca,'XTick',xtk,'XTickLabel',xlab,'FontSize',20,'LineWidth',2.75,'FontWeight','bold')
        
        set(gca,'FontName','Times New Roman','FontSize',16,'LineWidth',2.5,'FontWeight','bold');
        %exportgraphics(gcf,['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/' mountname '/basins_events/' gauge_no '/flood_thre.jpg'],'Resolution',1000)
        % exportgraphics does not work with parfor for some reason
        Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
        %[Left Bottom Right Top] spacing
        NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3)-0.01 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
        set(gca, 'Position', NewPos);
        saveas(gcf,['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/' mountname '/basins_events/' gauge_no '/flood_thre.jpg'])
        
        % remove dir if falsely created some graphs
        %rmdir(['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/' mountname '/basins_events/' gauge_no '/events/'],'s')
        
        % identify events below
        mkdir(['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/' mountname '/basins_events/' gauge_no '/events'])
        mkdir(['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/' mountname '/basins_events/' gauge_no '/events/summer_only'])
        mkdir(['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/' mountname '/basins_events/' gauge_no '/events/winter_only'])
        
        % actually selected events
        mkdir(['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/' mountname '/basins_events/' gauge_no '/selected_events'])
        
        % check how many events
        
        events = [];events_rank = [];
        for i = 1:length(obsnow)-2
            if obsnow(i)<obsthre&&obsnow(i+1)>obsthre&&obsnow(i+2)>obsthre
                obs_period = obsnow(max((i-15),1):min((i+15),length(obsnow)));
                if max(obs_period)>=2*(min(obs_period))
                    events = [events;i];
                    events_rank = [events_rank;max(obs_period)];
                end
            end
        end
        [~,idx] = sort(events_rank,'descend');
        events = events(idx);
        evendate = time_s(events);% start with 2020-12-05 00:00:00
        
        time_s_str = string;
        for ijk = 1:record_len
            % change time_s datetime data format to char
            time_s_str(ijk,1) = string(char(time_s{ijk}));
        end
        
        % save an event list as txt file
        if isempty(events)
            continue
            % no qualified events, could be very big basin with annually
            % changes in hydrographs.
        end
        evchar = char;
        for ijk = 1:length(events)
            evchar(ijk,:) = char(evendate{ijk});
        end
        mths = evchar(:,6:7);mths = str2num(mths);summer_index = find(mths>=4&mths<10);winter_index = find(mths<4|mths>=10);
        evchar = strcat(evchar(:,1:4),evchar(:,6:7),evchar(:,9:10));
        evchar = str2num(evchar);
        evchar_summer = evchar(summer_index);
        evchar_winter = evchar(winter_index);
        filename = ['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/' mountname '/basins_events/' gauge_no '/events/summer_only/event_list_summer.txt'];  % better to use fullfile(path,name)
        fid = fopen(filename,'w');    % open file for writing (overwrite if necessary)
        fprintf(fid,'%d\n',evchar_summer);          % Write the char array, interpret newline as new line
        fclose(fid);
        filename = ['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/' mountname '/basins_events/' gauge_no '/events/winter_only/event_list_winter.txt'];  % better to use fullfile(path,name)
        fid = fopen(filename,'w');    % open file for writing (overwrite if necessary)
        fprintf(fid,'%d\n',evchar_winter);          % Write the char array, interpret newline as new line
        fclose(fid);
        
        % plot each event
        for kp = 1:length(evendate)
            close all
            datetime.setDefaultFormats('default','yyyy-MM-dd HH:mm')
            
            indx = find(time_s_str==string(char(evendate{kp})));% find the event index in the big matric of observations
            if indx<=16||indx>=(length(obsnow)-18)
                continue
            end
            obstep = 15*2;
            obsp = obsnow(indx-15:indx+15);
            obstime = char(time_s_str(indx-15:indx+15,:));
            
            tp = char(evendate{kp});
            tp = tp(1:10);
            monthindex = str2num(tp(6:7));
            
            f = figure('visible','off');
            %         t = tiledlayout(1,1,'Padding','none');
            %         t.Units = 'inches';
            %         t.OuterPosition = [0.1 0.25 5.65 3.6];
            %         nexttile;
            f.Position = [5 5 550 300];
            
            plot(1:length(obsp),obsp,'LineWidth',1.5);
            %sh.MarkerEdgeColor = [0,0,0];
            
            box on
            
            legend('Obs.','modify datum and soil depth','original (drain 6 days)','Location','Northeast')
            legend('Obs.','± 5% obs.','Location','Northeast')
            
            ylabel('Streamflow (m^3/s)')
            legend boxoff
            ylim([0 ceil(max(obsp)/30)*30])
            xlim([0 length(obsp)+1])
            
            title([ country_nm '  ' river_nm '  ' gauge_no '  ' area_km2 'km^2'])
            xlabel([obstime(1,1:4) '-MM/DD'])
            
            %             targetdate = [];
            %             for ik = 1:length(obstime)
            %                 if obstime(ik,1:10)==obstime(97,1:10)
            %                     targetdate = [targetdate;ik]; % this day will be labeled diff color
            %                 end
            %             end
            
            xticlabel = strcat(string(obstime(1:6:length(obsp),6:7)),'/',string(obstime(1:6:length(obsp),9:10)));
            % label color different color
            %         xtic = 1:72:288*2;% interval is 5 minutes
            %         nlabel = length(xtic);
            %         timetic = 1:24:192;% obs is actually every 15 minutes
            %         xticlabel = obstime(timetic,12:13);
            %         xticlcolor = cell(nlabel,1);
            %         for ik = 1:nlabel
            %             if timetic(ik)>min(targetdate)&&timetic(ik)<max(targetdate)
            %                 xticlcolor{ik} = sprintf('\\color[rgb]{%f, %f, %f}%s', [1 0 0], xticlabel(ik,:));
            %             else
            %                 xticlcolor{ik} = sprintf('\\color[rgb]{%f, %f, %f}%s', [0 0 0], xticlabel(ik,:));
            %             end
            %         end
            %        set(gca,'XTick',1:length(obsp),'XTickLabel',xticlcolor,'FontSize',12,'LineWidth',2.5,'FontWeight','bold')
            
            set(gca,'XTick',1:6:length(obsp),'XTickLabel',xticlabel,'FontSize',12,'LineWidth',2.5,'FontWeight','bold')
            
            hax = gca;
            hax.XAxis.MinorTickValues = linspace(1,length(obsp),length(obsp));
            hax.XMinorTick = 'on';
            
            set(gca,'FontName','Times New Roman','FontSize',16,'LineWidth',2.5,'FontWeight','bold');
            Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
            %[Left Bottom Right Top] spacing
            NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3)-0.01 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
            set(gca, 'Position', NewPos);
            %exportgraphics(gcf,['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/' mountname '/basins_events/' gauge_no '/events/event' num2str(kp,'%2.2d') '.png'])
            saveas(gcf,['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/' mountname '/basins_events/' gauge_no '/events/event' num2str(kp,'%2.2d') '.png'])
            
            if monthindex>=4&&monthindex<=9
                saveas(gcf,['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/' mountname '/basins_events/' gauge_no '/events/summer_only/Event' num2str(kp,'%2.2d') '.png'])
            else
                saveas(gcf,['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/' mountname '/basins_events/' gauge_no '/events/winter_only/Event' num2str(kp,'%2.2d') '.png'])
            end
        end
    end
    
end

%% Alps, Himalayas, Andes create table of basin outlets lat lon, basin area, records length

home_alps = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Alps/streamflow_gauges_daily/';
%home_andes = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Andes/streamflow_gauges_daily/';
%home_himalayas = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Himalayas/streamflow_gauges_daily/';

home_andes = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Andes/streamflow_daily_entire_Andes/';
home_himalayas = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Himalayas/streamflow_daily_entire_Himalayas/';

for mountain = [4 5]
    if mountain == 1
        home_moun = home_alps;mountname = 'Alps';
    end
    if mountain == 4
        home_moun = home_andes;mountname = 'Andes';
    end
    if mountain == 5
        home_moun = home_himalayas;mountname = 'Himalayas';
    end
    lastall = ls([home_moun '*.txt']);
    c2 = strsplit(lastall);
    streamgauge_latlon = [];
    for ii = 1:length(c2)-1
        
        warning('off')
        close all
        ffnm = c2{ii};
        discharge=readtable(ffnm);
        if size(discharge,2)<3
            % no records at all
            continue
        end
        
        discharge_v = table2array(discharge(:,3));
        obsnow=discharge_v;
        if length(discharge_v)<30
            continue
        end
        disp(ii)
        
        
        fid=fopen(ffnm,'r');
        n_lines = 10;
        outlet_info = cell(10,1);
        for jj = 1:15
            outlet_info(jj) = {fgetl(fid)};
        end
        fclose(fid);
        
        gauge_no = outlet_info{9}(26:end);
        river_nm = outlet_info{10}(26:end);
        area_km2 = outlet_info{15}(30:end);
        country_nm = outlet_info{12}(26:27);
        lat = outlet_info{13}(24:end);
        lon = outlet_info{14}(24:end);
        streamgauge_latlon = [streamgauge_latlon;[str2double(lat),str2double(lon),str2double(gauge_no),str2double(area_km2),length(discharge_v)]];
        
    end
    if mountain==1
        writematrix(streamgauge_latlon,[home_moun 'Alps_all.csv']) 
%     elseif mountain==2
%         writematrix(streamgauge_latlon,[home_moun 'Andes.csv']) 
%     elseif mountain==3
%         writematrix(streamgauge_latlon,[home_moun 'Himalayas.csv']) 
    elseif mountain==4
        writematrix(streamgauge_latlon,[home_moun 'Andes_all.csv']) 
    elseif mountain==5
        writematrix(streamgauge_latlon,[home_moun 'Himalayas_all.csv'])    
    end
end

%% systematically pick events based on criteria below
% has to be above a threshold
% has to be fast
% has to be a single peak

home_alps = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Alps/streamflow_gauges_daily/';
%home_andes = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Andes/streamflow_gauges_daily/';
%home_himalayas = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Himalayas/streamflow_gauges_daily/';

home_andes = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Andes/streamflow_daily_entire_Andes/';
home_himalayas = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Himalayas/streamflow_daily_entire_Himalayas/';

%all_gauge_events = [];
%all_gauge_events = {};

for mountain = 3%[1 2 3] 
    if mountain == 1
        home_moun = home_alps;mountname = 'Alps';
        mount_b_e = readmatrix(['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Alps/input_output/Alps_all.csv']);
    end
    if mountain == 2
        home_moun = home_andes;mountname = 'Andes';
        mount_b_e = readmatrix(['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Andes/input_output/Andes_all.csv']);
    end
    if mountain == 3
        home_moun = home_himalayas;mountname = 'Himalayas';
        mount_b_e = readmatrix(['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Himalayas/input_output/Himalayas_all.csv']);
    end
    
    lastall = ls([home_moun '*.txt']);
    c2 = strsplit(lastall);
    su = 0;
    for ii = 1:length(c2)-1
  
        su = su+1;
        
        warning('off')
        close all
        ffnm = c2{ii};
        discharge=readtable(ffnm);
        if size(discharge,2)<3
            % no records at all
            continue
        end
        
        discharge_v = table2array(discharge(:,3));
        obsnow=discharge_v;
        if length(discharge_v)<30
            continue
        end
        disp(ii)
        
        
        fid=fopen(ffnm,'r');
        n_lines = 10;
        outlet_info = cell(10,1);
        for jj = 1:15
            outlet_info(jj) = {fgetl(fid)};
        end
        fclose(fid);
        
        gauge_no = outlet_info{9}(26:end);
        g_id = str2double(gauge_no);
        river_nm = outlet_info{10}(26:end);
        area_km2 = outlet_info{15}(30:end);
        country_nm = outlet_info{12}(26:27);
        
        
        time_s = table2cell(discharge(:,1));
        record_len = length(time_s);
        plot_step = floor(record_len/12);
        xtk = 1:plot_step:record_len;
        [ystart,mstart,dstart] = ymd(time_s{1});[yend,mend,dend] = ymd(time_s{end});
        xlab = string;
        
        for il = 1:plot_step:record_len
            [ys,ms,ds] = ymd(time_s{il});
            ysshort = mod(ys,100);
            xlab = [xlab;strcat(num2str(ysshort,'%2.2d'),'/',num2str(ms,'%2.2d'))];
        end
        xlab=cellstr(xlab);
        %
        f = figure('visible','on');
        f.Position = [10 10 1300 500];
        plot(discharge_v)
        ylim([0 max(discharge_v)])
        
        xlabel([num2str(ystart) ' - ' num2str(yend) ' YY/MM'])
        title([ country_nm '  ' river_nm '  ' gauge_no '  ' area_km2 'km^2'])
        ylabel('Streamflow (m^3/s)')
        xlim([0 length(time_s)+0.05*record_len])
        set(gca,'XTick',xtk,'XTickLabel',xlab,'FontSize',20,'LineWidth',2.75,'FontWeight','bold')
        Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
        %[Left Bottom Right Top] spacing
        NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3)-0.01 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
        set(gca, 'Position', NewPos);% when use this set the figure visible off for better cropping figure
        %saveas(gcf,['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/' mountname '/graphs/' country_nm '-' river_nm '-' gauge_no '-data_days' num2str(record_len)],'png')
        
        % ########################Event selection############################################################
        mkdir(['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/' mountname '/basins_events/'])
        mkdir(['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/' mountname '/basins_events/' gauge_no])
        
        
        ptl = 97.5;
        flowthre = prctile(discharge_v,ptl);% 99% will give 5-10 events per basin
        obsthre = flowthre;% when peak is above this value, it is considered as a flood event
        plot_step = floor(record_len/7);
        xtk = 1:plot_step:record_len;
        xlab = string;
        for il = 1:plot_step:record_len
            %kp=kp+1;
            [ys,ms,ds] = ymd(time_s{il});
            ysshort = mod(ys,100);
            xlab = [xlab;strcat(num2str(ysshort,'%2.2d'),'/',num2str(ms,'%2.2d'))];
        end
        xlab=cellstr(xlab);
        
        f = figure('visible','off');
        %         t = tiledlayout(1,1,'Padding','none');
        %         t.Units = 'inches';
        %         t.OuterPosition = [0.1 0.25 5.65 3.6];
        %         nexttile;
        f.Position = [5 5 700 400];
        sh = scatter(1:length(discharge_v),discharge_v,12,'b.');
        hold on
        plot([1 9999999],[obsthre obsthre],'r-.','LineWidth',2.5)
        hold off
        
        box on
        legend('Obs.',[num2str(obsthre,'%2.1f') 'm^3/s'],'interflow','± 5% obs.','Location','Northeast')
        ylabel('Streamflow (m^3/s)')
        %legend boxoff
        %yrf = prctile(obsnow,99.99);
        
        ylim([0 ceil(max(discharge_v)/30)*30])
        xlim([0 length(obsnow)])
        %title(['Basin' num2str(id,'%2.2d') ' Max=' num2str(obsmax,'%2.1f') 'm^3/s'])
        title([ country_nm '  ' river_nm '  ' gauge_no '  ' area_km2 'km^2'])
        xlabel([num2str(ystart) ' - ' num2str(yend) ' YY/MM'])
        set(gca,'XTick',xtk,'XTickLabel',xlab,'FontSize',20,'LineWidth',2.75,'FontWeight','bold')
        
        set(gca,'FontName','Times New Roman','FontSize',16,'LineWidth',2.5,'FontWeight','bold');
        %exportgraphics(gcf,['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/' mountname '/basins_events/' gauge_no '/flood_thre.jpg'],'Resolution',1000)
        % exportgraphics does not work with parfor for some reason
        Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
        %[Left Bottom Right Top] spacing
        NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3)-0.01 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
        set(gca, 'Position', NewPos);
        %saveas(gcf,['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/' mountname '/basins_events/' gauge_no '/flood_thre.jpg'])
        
        % remove dir if falsely created some graphs
        %rmdir(['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/' mountname '/basins_events/' gauge_no '/events/'],'s')
        
%         % identify events below
%         mkdir(['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/' mountname '/basins_events/' gauge_no '/events'])
%         mkdir(['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/' mountname '/basins_events/' gauge_no '/events/summer_only'])
%         mkdir(['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/' mountname '/basins_events/' gauge_no '/events/winter_only'])
%         
        % actually selected events
        mkdir(['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/' mountname '/basins_events/' gauge_no '/selected_events'])
        
        % check how many events
        g_loc = find(mount_b_e(:,3)==g_id);
        g_lat = mount_b_e(g_loc,1);%get the gauge latitude
        
        events = [];events_rank = [];
        for i = 1:length(obsnow)-2
            % the 0.5*obsthre is for the v2 of the events, which are much
            % more than original WORLD_events, because a lot of events
            % quickly drop after its peak
            if obsnow(i)<obsthre&&obsnow(i+1)>obsthre&&obsnow(i+2)>(0.5*obsthre)
                obs_period = obsnow(max((i-15),1):min((i+15),length(obsnow)));
                allpeaks = findpeaks(obs_period);
                obs_peak = max(obs_period);
                peak_tot = length(find(allpeaks>0.667*obs_peak));
                
                
                if max(obs_period)>=2*(min(obs_period))&&peak_tot<=2% can't have more than 2 significant peaks
                    ddate = time_s(i);ddate = ddate{1};[ytt,mtt,dtt] = ymd(ddate);
                    
                    
                    if g_lat>-42&&g_lat<42
                        events = [events;i];
                        events_rank = [events_rank;max(obs_period)];
                        
                    end
                    
                    if g_lat>=42
                        if mtt>=4&&mtt<=9% outside 30N 30S only warm season is selected
                            events = [events;i];
                            events_rank = [events_rank;max(obs_period)];
                        end
                        
                    end
                    
                    if g_lat<=-42
                        if mtt<4||mtt>9% outside 30N 30S only warm season is selected
                            events = [events;i];
                            events_rank = [events_rank;max(obs_period)];
                        end
                        
                    end
                    
                end
                

            end
        end
        [~,idx] = sort(events_rank,'descend');
        events = events(idx);
        evendate = time_s(events);% start with 2020-12-05 00:00:00
        
        tms = string(gauge_no);
        beinfo = [tms;evendate]';
        all_gauge_events{mountain,ii} = beinfo;

        time_s_str = string;
        for ijk = 1:record_len
            % change time_s datetime data format to char
            time_s_str(ijk,1) = string(char(time_s{ijk}));
        end
        
        % save an event list as txt file
        if isempty(events)
            continue
            % no qualified events, could be very big basin with annually
            % changes in hydrographs.
        end
        evchar = char;
        for ijk = 1:length(events)
            evchar(ijk,:) = char(evendate{ijk});
        end
        mths = evchar(:,6:7);mths = str2num(mths);summer_index = find(mths>=4&mths<10);winter_index = find(mths<4|mths>=10);
        evchar = strcat(evchar(:,1:4),evchar(:,6:7),evchar(:,9:10));
        evchar = str2num(evchar);
        evchar_summer = evchar(summer_index);
        evchar_winter = evchar(winter_index);
%         filename = ['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/' mountname '/basins_events/' gauge_no '/events/summer_only/event_list_summer.txt'];  % better to use fullfile(path,name)
%         fid = fopen(filename,'w');    % open file for writing (overwrite if necessary)
%         fprintf(fid,'%d\n',evchar_summer);          % Write the char array, interpret newline as new line
%         fclose(fid);
%         filename = ['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/' mountname '/basins_events/' gauge_no '/events/winter_only/event_list_winter.txt'];  % better to use fullfile(path,name)
%         fid = fopen(filename,'w');    % open file for writing (overwrite if necessary)
%         fprintf(fid,'%d\n',evchar_winter);          % Write the char array, interpret newline as new line
%         fclose(fid);
        
        % plot each event
        for kp = 1:length(evendate)
            
            close all
            datetime.setDefaultFormats('default','yyyy-MM-dd HH:mm')
            
            indx = find(time_s_str==string(char(evendate{kp})));% find the event index in the big matric of observations
            if indx<=16||indx>=(length(obsnow)-18)% truncate 15 days before the event above threshold timing and after 
                continue
            end
            obstep = 15*2;
            obsp = obsnow(indx-15:indx+15);
            obstime = char(time_s_str(indx-15:indx+15,:));
            
            tp = char(evendate{kp});
            tp = tp(1:10);
            monthindex = str2num(tp(6:7));
            
            f = figure('visible','off');
            %         t = tiledlayout(1,1,'Padding','none');
            %         t.Units = 'inches';
            %         t.OuterPosition = [0.1 0.25 5.65 3.6];
            %         nexttile;
            f.Position = [5 5 550 300];
            
            plot(1:length(obsp),obsp,'LineWidth',1.5);
            %sh.MarkerEdgeColor = [0,0,0];
            
            box on
            
            legend('Obs.','modify datum and soil depth','original (drain 6 days)','Location','Northeast')
            legend('Obs.','± 5% obs.','Location','Northeast')
            
            ylabel('Streamflow (m^3/s)')
            legend boxoff
            ylim([0 ceil(max(obsp)/30)*30])
            xlim([0 length(obsp)+1])
            
            title([ country_nm '  ' river_nm '  ' gauge_no '  ' area_km2 'km^2'])
            xlabel([obstime(1,1:4) '-MM/DD'])
            
            %             targetdate = [];
            %             for ik = 1:length(obstime)
            %                 if obstime(ik,1:10)==obstime(97,1:10)
            %                     targetdate = [targetdate;ik]; % this day will be labeled diff color
            %                 end
            %             end
            
            xticlabel = strcat(string(obstime(1:6:length(obsp),6:7)),'/',string(obstime(1:6:length(obsp),9:10)));
            % label color different color
            %         xtic = 1:72:288*2;% interval is 5 minutes
            %         nlabel = length(xtic);
            %         timetic = 1:24:192;% obs is actually every 15 minutes
            %         xticlabel = obstime(timetic,12:13);
            %         xticlcolor = cell(nlabel,1);
            %         for ik = 1:nlabel
            %             if timetic(ik)>min(targetdate)&&timetic(ik)<max(targetdate)
            %                 xticlcolor{ik} = sprintf('\\color[rgb]{%f, %f, %f}%s', [1 0 0], xticlabel(ik,:));
            %             else
            %                 xticlcolor{ik} = sprintf('\\color[rgb]{%f, %f, %f}%s', [0 0 0], xticlabel(ik,:));
            %             end
            %         end
            %        set(gca,'XTick',1:length(obsp),'XTickLabel',xticlcolor,'FontSize',12,'LineWidth',2.5,'FontWeight','bold')
            
            set(gca,'XTick',1:6:length(obsp),'XTickLabel',xticlabel,'FontSize',12,'LineWidth',2.5,'FontWeight','bold')
            
            hax = gca;
            hax.XAxis.MinorTickValues = linspace(1,length(obsp),length(obsp));
            hax.XMinorTick = 'on';
            
            set(gca,'FontName','Times New Roman','FontSize',16,'LineWidth',2.5,'FontWeight','bold');
            Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
            %[Left Bottom Right Top] spacing
            NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3)-0.01 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
            set(gca, 'Position', NewPos);
            %exportgraphics(gcf,['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/' mountname '/basins_events/' gauge_no '/events/event' num2str(kp,'%2.2d') '.png'])
            
            saveas(gcf,['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/' mountname '/basins_events/' gauge_no '/selected_events/event' num2str(kp,'%2.2d') '.png'])
            
            if monthindex>=4&&monthindex<=9
                saveas(gcf,['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/' mountname '/basins_events/' gauge_no '/events/summer_only/Event' num2str(kp,'%2.2d') '.png'])
            else
                saveas(gcf,['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/' mountname '/basins_events/' gauge_no '/events/winter_only/Event' num2str(kp,'%2.2d') '.png'])
            end

        end
    end
    
end

%% regarding [all_gauge_events] variable
% now only select gauges that are actually located in the mountains but not
% in the plains. 
% load final gauges that are manually selected in ArcGIS pro. 

homedir_final = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/';
selected_gauge_events = {};
count_events = 0;
for mountain = 1%[1 2 3]
    count_events_sub = 0;
    if mountain == 1
        mountname = 'Alps';
        mount_gauges = readmatrix(['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Alps/input_output/Alps_final_basins.xlsx']);
        % 2025/05/01 note: row 1 to row 65 has been missing all the time.
        % modify this excel by adding a 'x' at the row 1 col 3, then this
        % variable is named mount_gauges_miss_alps which has 200+gauges
        % including west and east alps
    end
    if mountain == 2
        mountname = 'Andes';
        mount_gauges = readmatrix(['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Andes/input_output/Andes_final_basins.xlsx']);
    end
    if mountain == 3
        mountname = 'Himalayas';
        mount_gauges = readmatrix(['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Himalayas/input_output/Himalayas_final_basins.xlsx']);
    end
    
    k = 0;
    for i = 1:size(all_gauge_events,2)
        tmp = all_gauge_events(mountain,i);
        if ~isempty(tmp{1})
            
            tmp2 = tmp{1};
            g_id = str2double(tmp2(1));
            if ismember(g_id, mount_gauges)
                % if this gauge is inclued in the final gauges (manually selected)
                % then the events are passed to the variable selected_gauge_events
                % before passing the events, first do a control on data
                % availability, only events after 1951 should be included
                % for ERA5 land data (one year spinup). 
                if size(tmp2,2)>1
                    k = k+1;
                    events1951 = string;
                    events1951(1) = tmp2(1);
                    era5count = 0;
                    for eeid = 2:size(tmp2,2)% the first entry in tmp2 is the gauge id, the following are events
                        tmp3 = char(tmp2(eeid));
                        yeartmp = str2num(tmp3(1:4));
                        if yeartmp>=1951
                            era5count = era5count+1;
                            events1951(1+era5count) = tmp2(eeid);
                            count_events = count_events + 1;
                            count_events_sub = count_events_sub + 1;
                            
                        end
                    end
                    selected_gauge_events{mountain,k} = events1951;
                
                end
            end
        end
    end
    count_events_sub
end
count_events

%% Prepare DEM, Flow_direction,and some soil data for each basin
% get streamgauge locations in the matrix
home_alps = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Alps/input_output/';
home_andes = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Andes/input_output/';
home_himalayas = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Himalayas/input_output/';



for mountain = [1 2 3]
    Basins_sum = [];
    if mountain == 1
        home_moun = home_alps;mountname = 'Alps';
        fid = fopen([home_moun 'basins_' mountname '_all.asc']);
        out = textscan(fid,repmat('%f',1,1855),'headerlines',6,'collectoutput',1);
        fclose(fid);
        region_raster = out{1};
        
        fid=fopen([home_moun 'basins_' mountname '_all.asc'],'r');
        raster_info = cell(5,1);
        for ii = 1:5
            raster_info(ii) = {fgetl(fid)};
        end
        fclose(fid);
        col = strsplit(raster_info{1}); col = str2double(col{2});
        row = strsplit(raster_info{2}); row = str2double(row{2});
        xcor = strsplit(raster_info{3}); xcor = str2double(xcor{2});
        ycor = strsplit(raster_info{4}); ycor = str2double(ycor{2});
        reso = strsplit(raster_info{5}); reso = str2double(reso{2});
        
        fid = fopen([home_moun mountname '_demf_1000.asc']);
        out = textscan(fid,repmat('%f',1,1855),'headerlines',6,'collectoutput',1);
        fclose(fid);
        demf_raster = out{1};
        fid = fopen([home_moun mountname '_fdr_1000.asc']);
        out = textscan(fid,repmat('%f',1,1855),'headerlines',6,'collectoutput',1);
        fclose(fid);
        fdr = out{1};
        fid = fopen([home_moun mountname '_facc_1000.asc']);
        out = textscan(fid,repmat('%f',1,1855),'headerlines',6,'collectoutput',1);
        fclose(fid);
        facc = out{1};
        
        fc = imread([home_moun 'alps_fc_1000'],'tif');
        awc = imread([home_moun 'alps_awc_1000'],'tif');wilt=fc-awc;
        ks = imread([home_moun 'alps_ks_log_10'],'tif');ks = 10.^(ks)/86400/100;
        poro = imread([home_moun 'alps_poro_1000'],'tif');
    end
    if mountain == 2
        home_moun = home_andes;mountname = 'Andes';
        fid = fopen([home_moun 'basins_' mountname '_all.asc']);
        out = textscan(fid,repmat('%f',1,9000),'headerlines',6,'collectoutput',1);
        fclose(fid);
        region_raster = out{1};
        
        fid=fopen([home_moun 'basins_' mountname '_all.asc'],'r');
        raster_info = cell(5,1);
        for ii = 1:5
            raster_info(ii) = {fgetl(fid)};
        end
        fclose(fid);
        col = strsplit(raster_info{1}); col = str2double(col{2});
        row = strsplit(raster_info{2}); row = str2double(row{2});
        xcor = strsplit(raster_info{3}); xcor = str2double(xcor{2});
        ycor = strsplit(raster_info{4}); ycor = str2double(ycor{2});
        reso = strsplit(raster_info{5}); reso = str2double(reso{2});
        
        fid = fopen([home_moun mountname '_demf_900.asc']);
        out = textscan(fid,repmat('%f',1,9000),'headerlines',6,'collectoutput',1);
        fclose(fid);
        demf_raster = out{1};
        fid = fopen([home_moun mountname '_fdr_900.asc']);
        out = textscan(fid,repmat('%f',1,9000),'headerlines',6,'collectoutput',1);
        fclose(fid);
        fdr = out{1};
        fid = fopen([home_moun mountname '_facc_900.asc']);
        out = textscan(fid,repmat('%f',1,9000),'headerlines',6,'collectoutput',1);
        fclose(fid);
        facc = out{1};
        
        fc = imread([home_moun 'andes_fc_1000'],'tif');
        awc = imread([home_moun 'andes_awc_1000'],'tif');wilt=fc-awc;
        ks = imread([home_moun 'andes_ks_log_10'],'tif');ks = 10.^(ks)/86400/100;
        poro = imread([home_moun 'andes_poro_1000'],'tif');
        
    end
    if mountain == 3
        home_moun = home_himalayas;mountname = 'Himalayas';
        fid = fopen([home_moun 'basins_' mountname '_all.asc']);
        out = textscan(fid,repmat('%f',1,6000),'headerlines',6,'collectoutput',1);
        fclose(fid);
        region_raster = out{1};
        fid=fopen([home_moun 'basins_' mountname '_all.asc'],'r');
        raster_info = cell(5,1);
        for ii = 1:5
            raster_info(ii) = {fgetl(fid)};
        end
        fclose(fid);
        col = strsplit(raster_info{1}); col = str2double(col{2});
        row = strsplit(raster_info{2}); row = str2double(row{2});
        xcor = strsplit(raster_info{3}); xcor = str2double(xcor{2});
        ycor = strsplit(raster_info{4}); ycor = str2double(ycor{2});
        reso = strsplit(raster_info{5}); reso = str2double(reso{2});
        
        fid = fopen([home_moun mountname '_demf_900.asc']);
        out = textscan(fid,repmat('%f',1,6000),'headerlines',6,'collectoutput',1);
        fclose(fid);
        demf_raster = out{1};
        fid = fopen([home_moun mountname '_fdr_900.asc']);
        out = textscan(fid,repmat('%f',1,6000),'headerlines',6,'collectoutput',1);
        fclose(fid);
        fdr = out{1};
        fid = fopen([home_moun mountname '_facc_900.asc']);
        out = textscan(fid,repmat('%f',1,6000),'headerlines',6,'collectoutput',1);
        fclose(fid);
        facc = out{1};
        
        fc = imread([home_moun 'himalayas_fc_1000'],'tif');
        awc = imread([home_moun 'himalayas_awc_1000'],'tif');wilt=fc-awc;
        ks = imread([home_moun 'himalayas_ks_log_10'],'tif');ks = 10.^(ks)/86400/100;
        poro = imread([home_moun 'himalayas_poro_1000'],'tif');
        
    end
    
    outlet_latlon =  csvread([home_moun mountname '_final_outlets_sum.csv'],1,0);
    gauge_id = xlsread([home_moun mountname '_final_basins.xlsx']);
    lat_region = [];lon_region = [];
    for i = 1:(col+1)
        for j = 1:(row+1)
            lat_region(j,i) = ycor+(j-1)*reso;
            lon_region(j,i) = xcor+(i-1)*reso;
        end
    end
    lat_region = flipud(lat_region);
    
    gauges_all = [];
    problem_gauges = [];
    for basin_num = 1:length(gauge_id)
        close all
        
        lat_tmp = outlet_latlon(basin_num,1);
        lon_tmp = outlet_latlon(basin_num,2);
        
        stmp = 1;
        for iu = 1:col
            if stmp == 0
                break
            end
            for ju = 1:row
                if lat_region(ju,iu)>lat_tmp&&lat_region(ju+1,iu)<lat_tmp&&lon_region(ju,iu)<lon_tmp&&lon_region(ju,iu+1)>lon_tmp
                    gauge_loc = [ju,iu,gauge_id(basin_num,2)];
                    stmp = 0;
                    break
                end
            end
        end
        gauges_all = [gauges_all;gauge_loc];%stream gauge locations in the big region matrix
        
        i = gauge_loc(1,1);
        j = gauge_loc(1,2);
        neigh8 = [i-1,j;i+1,j;i-1,j-1;i,j-1;i+1,j-1;i-1,j+1;i,j+1;i+1,j+1];
        ups = [i,j];
        basinpts = ups;
        border = 0;
        
        for itrr = 1:70000
            if border == 1
                break
            end
            vps = ups;
            uss = size(vps,1);
            ups = [];
            
            for s = 1:uss
                i = vps(s,1);
                j = vps(s,2);
                % enforce below for border control
                if i == 1||i==row||j==1||j==col
                    border = 1;
                    problem_gauges = [problem_gauges;gauge_loc];
                    break
                end
                % enforce done
                
                if fdr(i-1,j)==4
                    ups = [ups;[i-1,j]];
                end
                if fdr(i+1,j)==64
                    ups = [ups;[i+1,j]];
                end
                if fdr(i-1,j-1)==2
                    ups = [ups;[i-1,j-1]];
                end
                if fdr(i,j-1)==1
                    ups = [ups;[i,j-1]];
                end
                if fdr(i+1,j-1)==128
                    ups = [ups;[i+1,j-1]];
                end
                if fdr(i-1,j+1)==8
                    ups = [ups;[i-1,j+1]];
                end
                if fdr(i,j+1)==16
                    ups = [ups;[i,j+1]];
                end
                if fdr(i+1,j+1)==32
                    ups = [ups;[i+1,j+1]];
                end
                
                basinpts = [basinpts;ups];
                
            end
            
            if border == 1
                break
            end
            
            if isempty(ups)
                break
            end
            
            if itrr>40000
                disp(b)
                disp('warning')
            end
        end
        
        if border == 1
            problem_gauges = [problem_gauges;gauge_loc];
            continue
        end
        
        rowmin = min(basinpts(:,1))-3;
        rowmax = max(basinpts(:,1))+3;
        colmin = min(basinpts(:,2))-3;
        colmax = max(basinpts(:,2))+3;
        
        
        basin_plot = zeros(rowmax-rowmin+1,colmax-colmin+1);
        demf = demf_raster(rowmin:rowmax,colmin:colmax);
        fdrf = fdr(rowmin:rowmax,colmin:colmax);
        
        for ik = 1:length(basinpts)
            basin_plot(basinpts(ik,1)-rowmin+1,basinpts(ik,2)-colmin+1) = 1;
        end
        
        ba = basin_plot';
        
        intt = find(ba>0.5);
        ind = find(ba<0.5);
        
        %save(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gauge_id(basin_num,2)) '_intt'],'intt')
        %save(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gauge_id(basin_num,2)) '_ind'],'ind')
        %
        
        figure
        imagesc(basin_plot)
        
        basin_area = length(find(basin_plot>0))*100*100*reso*reso;
        
        demfplot = demf;
        demfplot(basin_plot==0) = 0;
        A = jet(20);
        A(1,:) = [1 1 1];
        
        
        lat_b1 = lat_region(rowmin:rowmax,colmin);
        lon_b1 = lon_region(rowmin,colmin:colmax);
        lat_b = lat_region(rowmin:rowmax,colmin:colmax);
        lon_b = lon_region(rowmin:rowmax,colmin:colmax);
        save(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/latlon/Basin' num2str(gauge_id(basin_num,2)) '_lat'],'lat_b')
        save(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/latlon/Basin' num2str(gauge_id(basin_num,2)) '_lon'],'lon_b')
        
        soils_poro = poro(rowmin:rowmax,colmin:colmax);
        soils_fc = fc(rowmin:rowmax,colmin:colmax);
        soils_wilt = wilt(rowmin:rowmax,colmin:colmax);
        soils_ks = ks(rowmin:rowmax,colmin:colmax);
        %save(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/soils/Basin' num2str(gauge_id(basin_num,2)) '_poro'],'soils_poro')
        %save(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/soils/Basin' num2str(gauge_id(basin_num,2)) '_fc'],'soils_fc')
        %save(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/soils/Basin' num2str(gauge_id(basin_num,2)) '_wilt'],'soils_wilt')
        %save(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/soils/Basin' num2str(gauge_id(basin_num,2)) '_ks'],'soils_ks')

        
        
        f = figure('visible','off');
        
        
        f.Position = [5 5 500 400];
        
        imagesc(lon_b1,lat_b1,demfplot)
        title([mountname ' Basin ' num2str(gauge_id(basin_num,2))])
        colormap(A)
        colorbar
        xlabel(['Longitude E{\circ}  ' num2str(basin_area,'%1.0f') 'km^2'])
        caxis([min(demfplot(:))-10 ceil(max(demfplot(:))/10)*10])
        %         hold on
        %         p = nsidedpoly(3, 'Center', [lon_tmp, lat_tmp], 'SideLength', 0.05);
        %         plot(p);
        %         hold off
        
        ylabel(['Latitude N{\circ}'])
        set(gca,'YDir','normal')
        set(get(colorbar,'Title'),'string','m')
        set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',2.5,'FontWeight','bold');
        %exportgraphics(gcf,[home_moun 'basins_dem/Basin_' num2str(gauge_id(basin_num,2)) '.jpg'],'Resolution',600)
        Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
        %[Left Bottom Right Top] spacing
        NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3)-0.18 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
        set(gca, 'Position', NewPos);% when use this set the figure visible off for better cropping figure
        %saveas(gcf,[home_moun 'basins_dem/Basin_' num2str(gauge_id(basin_num,2)) '.png'])
        
        p128 = find(fdrf==128);
        p64 = find(fdrf==64);
        p32 = find(fdrf==32);
        p16 = find(fdrf==16);
        p8 = find(fdrf==8);
        p4 = find(fdrf==4);
        p2 = find(fdrf==2);
        p1 = find(fdrf==1);
        
        fdrf(p128) = 1;
        fdrf(p64) = 128;
        fdrf(p32) = 64;
        fdrf(p16) = 32;
        fdrf(p8) = 16;
        fdrf(p4) = 8;
        fdrf(p2) = 4;
        fdrf(p1) = 2;
        
        tmpdem = demf';
        tmpfdr = fdrf';
        Basins_sum(basin_num,:) = [rowmax-rowmin+1,colmax-colmin+1,gauge_loc(1,1)-rowmin+1,gauge_loc(1,2)-colmin+1,basin_area,gauge_id(basin_num,2)];
        
        
%         mkdir([home_moun 'fortran/Basin' num2str(gauge_id(basin_num,2))])
%         fidd = fopen([home_moun 'fortran/Basin' num2str(gauge_id(basin_num,2)) '/dem.bin'],'w');
%         fwrite(fidd,tmpdem,'single');
%         fclose(fidd);
%         fidd = fopen([home_moun 'fortran/Basin' num2str(gauge_id(basin_num,2)) '/fdr.bin'],'w');
%         fwrite(fidd,tmpfdr,'ubit16');
%         fclose(fidd);
        
    end
%     save([home_moun 'gauges_loc/Basins_dimention.mat'],'Basins_sum')
%     save([home_moun 'gauges_loc/Gauges_loc.mat'],'gauges_all')
%     save([home_moun 'gauges_loc/Gauges_loc_problem.mat'],'problem_gauges')
%     figure
%     t = tiledlayout(1,1,'Padding','none');
%     t.Units = 'inches';
%     t.OuterPosition = [0.25 0.25 4.5 3.4];
%     nexttile;
%     histogram(Basins_sum(:,5),0:100:6500)
%     xlabel('Drainage Area (km^2)')
%     ylabel('Count')
%     ylim([0 10])
%     set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',2.5,'FontWeight','bold');
%     exportgraphics(gcf,[home_moun 'Basins_size_distributions.jpg'],'Resolution',600)
%     
    
end

%% Prepare MPI variables for parallel computing for world basins
% Just check below
ffnm=['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Alps/input_output/fortran/Basin6559110/countfile.out'];
fid=fopen(ffnm);
faccAPL = fscanf(fid,'%d',[24,13]);
fclose(fid);
figure
imagesc(faccAPL')
% end of check

home_alps = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Alps/input_output/';
home_andes = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Andes/input_output/';
home_himalayas = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Himalayas/input_output/';

problemgauges = [];
for mountain = [1 2 3]
    
    if mountain == 1
        home_moun = home_alps;mtname = 'Alps';
        fid=fopen([home_moun 'basins_' mtname '_all.asc'],'r');
        raster_info = cell(5,1);
        for ii = 1:5
            raster_info(ii) = {fgetl(fid)};
        end
        fclose(fid);
        col = strsplit(raster_info{1}); col = str2double(col{2});
        row = strsplit(raster_info{2}); row = str2double(row{2});
        xcor = strsplit(raster_info{3}); xcor = str2double(xcor{2});
        ycor = strsplit(raster_info{4}); ycor = str2double(ycor{2});
        reso = strsplit(raster_info{5}); reso = str2double(reso{2});
    end
    if mountain == 2
        home_moun = home_andes;mtname = 'Andes';
        fid=fopen([home_moun 'basins_' mtname '_all.asc'],'r');
        raster_info = cell(5,1);
        for ii = 1:5
            raster_info(ii) = {fgetl(fid)};
        end
        fclose(fid);
        col = strsplit(raster_info{1}); col = str2double(col{2});
        row = strsplit(raster_info{2}); row = str2double(row{2});
        xcor = strsplit(raster_info{3}); xcor = str2double(xcor{2});
        ycor = strsplit(raster_info{4}); ycor = str2double(ycor{2});
        reso = strsplit(raster_info{5}); reso = str2double(reso{2});
    end
    if mountain == 3
        home_moun = home_himalayas;mtname = 'Himalayas';
        fid=fopen([home_moun 'basins_' mtname '_all.asc'],'r');
        raster_info = cell(5,1);
        for ii = 1:5
            raster_info(ii) = {fgetl(fid)};
        end
        fclose(fid);
        col = strsplit(raster_info{1}); col = str2double(col{2});
        row = strsplit(raster_info{2}); row = str2double(row{2});
        xcor = strsplit(raster_info{3}); xcor = str2double(xcor{2});
        ycor = strsplit(raster_info{4}); ycor = str2double(ycor{2});
        reso = strsplit(raster_info{5}); reso = str2double(reso{2});
    end
    
    lat_region = [];lon_region = [];
    for i = 1:(col+1)
        for j = 1:(row+1)
            lat_region(j,i) = ycor+(j-1)*reso;
            lon_region(j,i) = xcor+(i-1)*reso;
        end
    end
    lat_region = flipud(lat_region);
    gauges_all = load([home_moun 'gauges_loc/Gauges_loc.mat'],'gauges_all');
    gauges_all = gauges_all.gauges_all;
    
    tmp = load([home_moun 'gauges_loc/Basins_dimention.mat'],'Basins_sum');
    Basins_sum_region = tmp.Basins_sum;
    
    
    % use countfile from fortran DCHM, MPI delineation
    kou = 0;
    for bid = 1:size(Basins_sum_region,1)
        kou = kou+1;
        if mod(kou,10)==0
            close all
        end
        
        c = Basins_sum_region(bid,2);
        r = Basins_sum_region(bid,1);
        gid = Basins_sum_region(bid,6);
        %
        g_loc = find(gauges_all(:,3)==gid);
        g_row = gauges_all(g_loc,1);
        g_col = gauges_all(g_loc,2);
        
        
        ffnm=[home_moun 'fortran/Basin' num2str(gid) '/fdr.bin'];
        fid=fopen(ffnm);
        tmpdata=fread(fid,'ubit16');
        fdr = reshape(tmpdata,c,r);
        fclose(fid);
        % countfiles were in the same folders (by running Runfdr.sh in '~/IRC_LB/Other_mountains/derive_facc/)
        ffnm=[home_moun 'fortran/Basin' num2str(gid) '/countfile.out'];
        fid=fopen(ffnm);
        faccAPL = fscanf(fid,'%d',[c,r]);
        fclose(fid);
        figure
        imagesc(faccAPL')
        
        if length(faccAPL)<3
            problemgauges = [problemgauges;gid];
            continue
        end
        
        if size(fdr,1)~=size(faccAPL,1)||size(fdr,2)~=size(faccAPL,2)
            problemgauges = [problemgauges;gid];
            continue
        end
        
        flowmax = max(max(faccAPL))-1;
        B_plan = 1;
        P_plan = 1;
        fmean = flowmax/B_plan/P_plan*2;
        f_low = 1.25*fmean;
        f_max = 1.75*fmean;
        f_max = max(max(faccAPL))-10;% ensure B1 P1
        %f_max=1200;% 200 180 are good
        bpsum = cell(1,10);
        bpsum{1} = [Basins_sum_region(bid,4),Basins_sum_region(bid,3)];% first use column, second use row
        pt = 1;
        branchmap = zeros(Basins_sum_region(bid,2),Basins_sum_region(bid,1));
        reachesmap = zeros(Basins_sum_region(bid,2),Basins_sum_region(bid,1));
        critpt = zeros(Basins_sum_region(bid,2),Basins_sum_region(bid,1));
        for joints = 1:pt
            ifix = bpsum{1}(1);
            jfix = bpsum{1}(2);
            i = ifix;
            j = jfix;
            branchmap(ifix,jfix) = 1;
            reachesmap(ifix,jfix) = faccAPL(ifix,jfix);
            critpt(ifix,jfix) = 1;
            neigh8 = [i-1,j;i+1,j;i-1,j-1;i,j-1;i+1,j-1;i-1,j+1;i,j+1;i+1,j+1];
            ups = [i,j];
            flownet = [faccAPL(i,j)];
            for itrr = 1:70000
                vps = ups;
                uss = size(vps,1);
                ups = [];% now everytime only new added pixel are tracking upwads
                bc = 0;
                for s = 1:uss% loop through every pixel added last iteration
                    i = vps(s,1);
                    j = vps(s,2);
                    posbpt = [];
                    if fdr(i-1,j)==2
                        ups = [ups;[i-1,j]];
                        flownet = [flownet;faccAPL(i-1,j)];
                        if faccAPL(i-1,j)>(f_max)% if meets,then it is upstream pixel which is also qualifed for minimal flow accumulation condition
                            posbpt = [posbpt;[i-1,j]];
                        end
                    end
                    if fdr(i+1,j)==32
                        ups = [ups;[i+1,j]];
                        flownet = [flownet;faccAPL(i+1,j)];
                        if faccAPL(i+1,j)>(f_max)
                            posbpt = [posbpt;[i+1,j]];
                        end
                    end
                    if fdr(i-1,j-1)==4
                        ups = [ups;[i-1,j-1]];
                        flownet = [flownet;faccAPL(i-1,j-1)];
                        if faccAPL(i-1,j-1)>(f_max)
                            posbpt = [posbpt;[i-1,j-1]];
                        end
                    end
                    if fdr(i,j-1)==8
                        ups = [ups;[i,j-1]];
                        flownet = [flownet;faccAPL(i,j-1)];
                        if faccAPL(i,j-1)>(f_max)
                            posbpt = [posbpt;[i,j-1]];
                        end
                    end
                    if fdr(i+1,j-1)==16
                        ups = [ups;[i+1,j-1]];
                        flownet = [flownet;faccAPL(i+1,j-1)];
                        if faccAPL(i+1,j-1)>(f_max)
                            posbpt = [posbpt;[i+1,j-1]];
                        end
                    end
                    if fdr(i-1,j+1)==1
                        ups = [ups;[i-1,j+1]];
                        flownet = [flownet;faccAPL(i-1,j+1)];
                        if faccAPL(i-1,j+1)>(f_max)
                            posbpt = [posbpt;[i-1,j+1]];
                        end
                    end
                    if fdr(i,j+1)==128
                        ups = [ups;[i,j+1]];
                        flownet = [flownet;faccAPL(i,j+1)];
                        if faccAPL(i,j+1)>(f_max)
                            posbpt = [posbpt;[i,j+1]];
                        end
                    end
                    if fdr(i+1,j+1)==64
                        ups = [ups;[i+1,j+1]];
                        flownet = [flownet;faccAPL(i+1,j+1)];
                        if faccAPL(i+1,j+1)>(f_max)
                            posbpt = [posbpt;[i+1,j+1]];
                        end
                    end
                    
                    alpt = [];
                    branchnum = size(posbpt,1);
                    if branchnum>0% meaning the upstream pixel of the pixel is not qualified to be recognized(flow accu< a threshold)
                        for sk = 1:size(posbpt,1)
                            alpt(sk,:) = [posbpt(sk,1),posbpt(sk,2),faccAPL(posbpt(sk,1),posbpt(sk,2))];
                        end
                        alptsort = sortrows(alpt,3,'descend');
                        if branchnum == 1% in the same stream line
                            branchmap(alptsort(1,1),alptsort(1,2))=branchmap(i,j);
                            reachesmap(alptsort(1,1),alptsort(1,2))=reachesmap(i,j);
                            if faccAPL(alptsort(1,1),alptsort(1,2))<reachesmap(alptsort(1,1),alptsort(1,2))-f_max% when meets, we need a break, otherwise twoo big area for one processor.
                                reachesmap(alptsort(1,1),alptsort(1,2))=faccAPL(alptsort(1,1),alptsort(1,2));
                                branchmap(alptsort(1,1),alptsort(1,2))=branchmap(i,j)*10+1;%flow into its downstream e.g. 1111 flow into 111
                                critpt(alptsort(1,1),alptsort(1,2)) = 1;
                            end
                        else% now when go upstream, there are two or more options, meaning more branches.
                            for sk = 1:size(alptsort,1)
                                if sk<=(size(alptsort,1)-1)% for big branches, extend a couple of pixels back for MPI to work.MPI cant have multiple different processors with outlets leading to the same downstream pixel
                                    branchmap(alptsort(sk,1),alptsort(sk,2)) = branchmap(i,j);
                                    reachesmap(alptsort(sk,1),alptsort(sk,2)) = faccAPL(alptsort(sk,1),alptsort(sk,2));
                                else
                                    if ismember([branchmap(i,j)*10+sk],branchmap)% to avoid duplicate naming for sub branches
                                        branchmap(alptsort(sk,1),alptsort(sk,2)) = max(branchmap(:))+1;
                                        reachesmap(alptsort(sk,1),alptsort(sk,2))=faccAPL(alptsort(sk,1),alptsort(sk,2));
                                        critpt(alptsort(sk,1),alptsort(sk,2)) = 1;
                                    else
                                        branchmap(alptsort(sk,1),alptsort(sk,2)) = branchmap(i,j)*10+sk;
                                        reachesmap(alptsort(sk,1),alptsort(sk,2))=faccAPL(alptsort(sk,1),alptsort(sk,2));
                                        critpt(alptsort(sk,1),alptsort(sk,2)) = 1;
                                    end
                                end
                            end
                            
                        end
                        
                    end
                    
                end
                
                
                if isempty(ups)
                    break
                end
            end
        end
        itrr
        
        [outx,outy] = find(critpt>0);
        
        
        fdr_rot = fdr';
        fdrsubdeli = fdr_rot;% changing from fortran fdr to what subdeli can use fdr.
        p128 = find(fdr_rot==128);
        p64 = find(fdr_rot==64);
        p32 = find(fdr_rot==32);
        p16 = find(fdr_rot==16);
        p8 = find(fdr_rot==8);
        p4 = find(fdr_rot==4);
        p2 = find(fdr_rot==2);
        p1 = find(fdr_rot==1);
        
        fdrsubdeli(p128) = 32;
        fdrsubdeli(p64) = 64;
        fdrsubdeli(p32) = 128;
        fdrsubdeli(p16) = 1;
        fdrsubdeli(p8) = 2;
        fdrsubdeli(p4) = 4;
        fdrsubdeli(p2) = 8;
        fdrsubdeli(p1) = 16;
        
        spub = max(branchmap(:));
        max_digits = numel(num2str(spub));
        colm = 0;sall = [];rall = 0;
        allr = {};
        for digits = max_digits:-1:1
            colm = colm+1;
            s1 = unique(branchmap(find(branchmap>=(10^(digits-1))&branchmap<(2*10^(digits-1)))));
            for i = 1:length(s1)
                sall(i,colm) = s1(i);
                [lx,ly] = find(branchmap==s1(i)&critpt>0);
                rall = rall+1;
                allr{rall} = subdeli(ly,lx,fdrsubdeli,r,c);
            end
        end
        sall
        cellnum = length(find(sall(:)>0));
        psub = cell(cellnum,1);
        sall_nozero = sall(sall>0);
        outpointsr = [];
        
        for psu = 1:cellnum
            [lx,ly] = find(branchmap==sall_nozero(psu)&critpt>0);
            outpointsr = [outpointsr;[ly,lx]];
            if psu<=sum(sall(:,1)>0)% the first column get to be assigned to be part of the first processed
                psub{psu} = allr{psu};
            end
            if psu>sum(sall(:,1)>0)
                flo_up = floor(sall_nozero/10);
                [~,locts] = ismember(sall_nozero(psu),flo_up);
                if locts==0
                    psub{psu} = allr{psu};
                else
                    up_loc = find(flo_up==sall_nozero(psu));
                    up_num = length(up_loc);
                    tmpal = allr{psu};
                    for ik = 1:up_num
                        tmpal = setdiff(tmpal,allr{up_loc(ik)},'rows');
                    end
                    psub{psu} = tmpal;
                end
            end
        end
        A = jet(2);A(1,:) = [1 1 1];A(2,:) = [0 0 0];
        faccmap = faccAPL;faccmap(faccmap>5) = 500000;
        
        %     usge gauge locations (location in the big region) and basins
        %     information to locate the lat lon extent
        %     lat_b = lat_region((g_row-Basins_sum_region(bid,3)+1):(g_row-Basins_sum_region(bid,3)+Basins_sum_region(bid,1)),1:Basins_sum_region(bid,2));
        %     lon_b = lon_region(1:Basins_sum_region(bid,1),(g_col-Basins_sum_region(bid,4)+1):(g_col-Basins_sum_region(bid,4)+Basins_sum_region(bid,2)));
        %
        %     outx = lat_b(outx);
        %     outy = lon_b(outy);
        
        figure
        t = tiledlayout(1,1,'Padding','none');
        t.Units = 'inches';
        t.OuterPosition = [0.25 0.25 4.85 4.6];
        nexttile
        imagesc(faccmap')
        colormap(A)
        
        for i = 1:cellnum
            hold on
            cco = [rand(1,3)];
            scatter(psub{i}(:,2),psub{i}(:,1),6,cco,'filled');
            %scatter(outy,outx,6,cco,'filled');
        end
        hold on
        opu = scatter(outx,outy,25,'filled');
        opu.MarkerFaceColor = [1,0,0];
        
        hold off
        
        
        title([mtname ' GID: ' num2str(gid) ''])
        xlabel('Grid')
        ylabel('Grid')
        %xlabel(['Longitude E{\circ}'])
        %ylabel(['Latitude N{\circ}'])
        %set(gca,'YDir','normal')
        
        set(gca,'FontName','Calibri','FontSize',16,'LineWidth',3,'FontWeight','Bold');
        fig = gcf;
        fig.PaperUnits = 'inches';
        fig.PaperPosition = [0,0,4.6,4.2];
        exportgraphics(gcf,[home_moun 'basins_dem/MPI_Basin' num2str(gid) '.jpg'],'Resolution',300)
        
        Bmax = size(sall,2);
        pset = [];
        for i = 1:Bmax
            pset(i) = sum(sall(:,i)>0);
        end
        pset
        % write MPI vairables
        c = Basins_sum_region(bid,2);
        r = Basins_sum_region(bid,1);
        
        b=Bmax;
        p=max(pset);
        ntr = 5;
        homedir = ['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_fixed_900m/'];
        mkdir([homedir])
        mkdir([homedir 'ParalleliziedBasin_B' num2str(b) '_P' num2str(p) ''])
        
        ffnm=[home_moun 'fortran/Basin' num2str(gid) '/dem.bin'];
        fid=fopen(ffnm,'rb','ieee-le');
        tmpdata=fread(fid,inf,'single');
        tmp=reshape(tmpdata,c,r);
        dem = tmp';
        
        
        
        dpr = faccAPL';
        
        cou = 1;
        NLR = 7000;
        NCR = 7000;
        NODP = 99900;
        
        for bba = 1:length(pset)
            for partt = 1:pset(bba)
                
                tmp = psub{cou};
                % lateral points
                allpb = strcat('lateralrouting_B',num2str(bba),'_P',num2str(partt));
                xpp = sortrows(unique(dpr(tmp(:,3)),'rows'));
                lex = length(xpp);
                if (NLR-lex)<0
                    disp('error')
                end
                fxpp = [xpp;zeros(NLR-lex,1)];
                
                dlmwrite([homedir 'ParalleliziedBasin_B' num2str(b) '_P' num2str(p) '/' allpb '.ascii'],fxpp)
                
                % channel  pixel
                allpc = strcat('channelrouting_B',num2str(bba),'_P',num2str(partt));
                chapp = xpp(find(xpp>ntr));
                ley = length(chapp);
                if (NCR-lex)<0
                    disp('error2')
                end
                cpp = [chapp;zeros(NCR-ley,1)];
                
                dlmwrite([homedir 'ParalleliziedBasin_B' num2str(b) '_P' num2str(p) '/' allpc '.ascii'],cpp)
                
                dems = -1*ones(r,c);
                dems(tmp(:,3)) = dem(tmp(:,3));
                demf = dems';
                fidd = fopen([homedir 'ParalleliziedBasin_B' num2str(b) '_P' num2str(p) '/dem_B',num2str(bba),'_P',num2str(partt),'.bin'],'w');
                fwrite(fidd,demf,'single');
                fclose(fidd);
                
                allpe = strcat('dem_B',num2str(bba),'_P',num2str(partt));
                %varr.(allpe) = dems;
                dlmwrite([homedir 'ParalleliziedBasin_B' num2str(b) '_P' num2str(p) '/' allpe '.ascii'],demf)
                
                [sr,sc] = ind2sub([r,c],tmp(:,3));
                
                aft = sub2ind([c,r],sc,sr);
                
                % one dimensional point
                oned = sortrows(aft);
                lez = length(oned);
                if (NODP-lez)<0
                    disp('error3')
                end
                od = [oned;zeros(NODP-lez,1)];
                allpd = strcat('onedimensionalpoint_B',num2str(bba),'_P',num2str(partt));
                dlmwrite([homedir 'ParalleliziedBasin_B' num2str(b) '_P' num2str(p) '/' allpd '.ascii'],od)
                
                cou = cou+1;
            end
        end
        
        
        cou = 1;
        allo = [];
        for bba = 1:length(pset)
            outg = [];
            for partt = 1:pset(bba)
                outg = [outg,[outpointsr(cou,1) outpointsr(cou,2)]];
                cou = cou+1;
            end
            if partt<p
                su = p-partt;
                outg = [outg zeros(1,2*su)];
            end
            
            allo = [allo;outg];
        end
        
        dlmwrite([homedir '/ParalleliziedBasin_B' num2str(b) '_P' num2str(p) '/outlet_group.ascii'],allo)
        
        % write subbasin matrix
        
        cou = 1;
        suma = [];
        for bba = 1:length(pset)
            outgb = [];
            for partt = 1:pset(bba)
                
                tmp = psub{cou};
                
                outgb = [outgb,[min(tmp(:,1))-1 max(tmp(:,1))+1 min(tmp(:,2))-1 max(tmp(:,2))+1]];
                cou = cou+1;
            end
            if partt<p
                su = p-partt;
                outgb = [outgb zeros(1,4*su)];
            end
            
            suma = [suma;outgb];
        end
        dlmwrite([homedir 'ParalleliziedBasin_B' num2str(b) '_P' num2str(p) '/subbasin_matrix.ascii'],suma)
        matr = 30*ones(b,b*6);
        dlmwrite([homedir 'ParalleliziedBasin_B' num2str(b) '_P' num2str(p) '/subbasin_matrixunique.ascii'],matr)
        
        %}
        
        surflow = [5;3600;50;900;900;2000;2000;2000;2000;20;1];
        
        dlmwrite([homedir 'ParalleliziedBasin_B' num2str(b) '_P' num2str(p) '/surflow.ascii'],surflow)
        dlmwrite([homedir 'surflow.ascii'],surflow)
        
        %
        disp('Also copy fdr and dem and countfile to the fixed folder')
        % the countfile was copied in Runfdr.sh script from model outputs to /input/output/fortran folder,and then here copied again to fixed input folder 
        copyfile([home_moun 'fortran/Basin' num2str(gid) '/*'],homedir)
        %
    end
    
end

%% write additional soil-related inputs


home_alps = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Alps/input_output/';
home_andes = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Andes/input_output/';
home_himalayas = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Himalayas/input_output/';

for mountain = 1%[1 3 2]
    
    if mountain == 1
        home_moun = home_alps;mtname = 'Alps';
        
    end
    if mountain == 2
        home_moun = home_andes;mtname = 'Andes';
        
    end
    if mountain == 3
        home_moun = home_himalayas;mtname = 'Himalayas';
        
    end
    tmp = load([home_moun 'gauges_loc/Basins_dimention.mat'],'Basins_sum');
    Basins_sum_region = tmp.Basins_sum;
    kou = 1;
    for bid = 1:size(Basins_sum_region,1)
        kou = kou+1;
        if mod(kou,10)==0
            close all
        end
        
        cols = Basins_sum_region(bid,2);
        rows = Basins_sum_region(bid,1);
        gid = Basins_sum_region(bid,6);
        
        
        % first find the gd level
        %{
    filename = ['/shared/dondo/home/ml423/HMIOP/NoDAre/LA/obs/Basin' num2str(bid,'%2.2d') '_gd.txt'];
    fid = fopen(filename);
    tmp = textscan(fid,'%s','delimiter','\n');
    tmp1 = strfind(tmp{1}(1:50),'5s');% find where the headline ends
    
    for headline = 1:length(tmp1)
        if isempty(tmp1{headline})
        else
           headends = headline;
        end
    end
    
    filename = ['/shared/dondo/home/ml423/HMIOP/NoDAre/LA/obs/Basin' num2str(bid,'%2.2d') '_gd.txt'];
    fid = fopen(filename);
    formatspec=['%s%s%s%s%s%s%s%s%s'];%,repmat('%*s',1,20)];
    fobs = textscan(fid,formatspec,'HeaderLines',headends);
    depth_under_dem = fobs{6};
    dud = str2num(depth_under_dem{1})*0.3048;
        %}
        % end find the gd level
        %{
        
        tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_ind.mat'],'ind');
        ind = tmp.ind;
        tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_intt.mat'],'intt');
        intt = tmp.intt;
        
        
        
        
        ffnm=[home_moun 'fortran/Basin' num2str(gid) '/countfile.out'];
        fid=fopen(ffnm);
        facc=fscanf(fid,'%d',[cols,rows]);
        fclose(fid);
        
        ffnm=[home_moun 'fortran/Basin' num2str(gid) '/dem.bin'];
        fid=fopen(ffnm,'rb','ieee-le');
        tmpdata=fread(fid,inf,'single');
        fclose(fid);
        demb1 = reshape(tmpdata,cols,rows);
        
        totdep = 20.67;
        shadep = 1.75;
        
        l2min = 0.15;
        l2max = 0.65;
        l3min = 0.2;
        l3max = 0.9;
        
        d0 = 0.1*ones(cols,rows);
        
        datum_new = [];
        gd_new = [];
        
        dem0 = min(demb1(:));
        deminf = max(demb1(:));
        middem = median(demb1(intt));
        m20 = prctile(demb1(intt),20);
        m10 = prctile(demb1(intt),10);
        m5 = prctile(demb1(intt),5);
        m1 = prctile(demb1(intt),0.5);
        clear soild
        clear l2depth
        clear l3depth
        for i = 1:rows
            for j = 1:cols
                soild(j,i) = totdep-(totdep-shadep)*(demb1(j,i)-dem0)/(deminf-dem0);
                l2depth(j,i) = l2max-(l2max-l2min)*(demb1(j,i)-dem0)/(deminf-dem0);
                l3depth(j,i) = l3max-(l3max-l3min)*(demb1(j,i)-dem0)/(deminf-dem0);
            end
        end
        datum_new = demb1-soild;
        
        gd_new = zeros(cols,rows);
        
        gd_new(find(gd_new<datum_new))=datum_new(find(gd_new<datum_new));
        
        
        hihh = find(facc>5&demb1<m20);
        
        gd_new(hihh) = demb1(hihh)-d0(hihh)-l2depth(hihh)-l3depth(hihh)-0.0002;
        gd_new(find(gd_new<datum_new))=datum_new(find(gd_new<datum_new));
        
        
        at = gd_new-datum_new;
        
        avai = demb1-gd_new;
        baset = demb1-datum_new-0.1-l2depth-l3depth;
        min(baset(intt))
        at(ind) = 0;
        avai(ind) = 0;
        figure
        imagesc(at')
        title(['900m AT ' num2str(mean(at(intt)))])
        colorbar
        
        figure
        imagesc(avai')
        title(['900m available ' num2str(mean(avai(intt)))])
        colorbar
        
        figure
        imagesc(baset')
        title(['900m base layer thickness has to be >0'])
        colorbar
        %
        en2 = reshape(datum_new,rows*cols,1);
        en1 = reshape(gd_new,rows*cols,1);
        soii = [en1,en2];
        %dlmwrite(['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_fixed_900m/initialwatertable_basin.ascii'],soii,'precision',10)
        
        mean(l2depth(intt))
        mean(l3depth(intt))
        % soil depth and soil parameter
        
        dsmm = reshape(l2depth,rows*cols,1);
        dpm = reshape(l3depth,rows*cols,1);
        
        soii = [dsmm,dpm];
        %dlmwrite(['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_fixed_900m/soildepth.ascii'],soii,'precision',10)
        
        %}
        
        poro = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/soils/Basin' num2str(gid) '_poro'],'soils_poro');
        fc = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/soils/Basin' num2str(gid) '_fc'],'soils_fc');
        wilt = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/soils/Basin' num2str(gid) '_wilt'],'soils_wilt');
        ks = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/soils/Basin' num2str(gid) '_ks'],'soils_ks');
        poro = poro.soils_poro;
        fc = fc.soils_fc;
        wilt = wilt.soils_wilt;
        ks = ks.soils_ks;
        
        
        
        kh1 = ks';goodloc= find(kh1>0);badloc = find(kh1<=0|kh1>0.01);
        kh1(badloc) = min(kh1(goodloc));
        kh2 = ks';
        kh3 = ks';
        kh4 = ks';
        poro1 = poro';goodloc= find(poro1>0);badloc = find(poro1<=0|poro1>0.8);
        poro1(badloc) = min(poro1(goodloc));
        poro2 = poro';
        poro3 = poro';
        poro4 = poro';
        fcp1 = fc';goodloc= find(fcp1>0);badloc = find(fcp1<=0|fcp1>0.8);
        fcp1(badloc) = min(fcp1(goodloc));
        fcp2 = fc';
        fcp3 = fc';
        fcp4 = fc';
        wilt1 = wilt';goodloc= find(wilt1>0);badloc = find(wilt1<=0|wilt1>0.8);
        wilt1(badloc) = min(wilt1(goodloc));
        wilt2 = wilt';
        wilt3 = wilt';
        wilt4 = wilt';
        so1 = reshape(kh1,[cols*rows,1]);so2=so1;so3=so1;so4=so1;
        so5 = reshape(poro1,[cols*rows,1]);so6=so5;so7=so5;so8=so5;
        so9 = reshape(fcp1,[cols*rows,1]);so10=so9;so11=so9;so12=so9;
        so13 = reshape(wilt1,[cols*rows,1]);so14=so13;so15=so13;so16=so13;
        
        
        soi = [so1,so2,so3,so4,so5,so6,so7,so8,so9,so10,so11,so12,so13,so14,so15,so16];
        dlmwrite(['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_fixed_900m/soilparameters.ascii'],soi,'precision',10)
        
        %green ampt and surface roughness from literature
        
        sf = 0.13*ones(rows,cols);
        lamda = 0.25;
        n = 3+2./lamda;
        na = n*ones(rows,cols);
        
        sff = reshape(sf',rows*cols,1);
        naf = reshape(na',rows*cols,1);
        gam = [sff,naf];
        %dlmwrite(['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_fixed_900m/green-ampt.ascii'],gam,'precision',10)
        
        
        anmp = 0.05*ones(cols,rows);
        anmc = 0.03*ones(cols,rows);
        
        aaa = reshape(anmp,rows*cols,1);
        ccc = reshape(anmc,rows*cols,1);
        sal = [aaa,ccc];
        %dlmwrite(['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_fixed_900m/surface_roughness.ascii'],sal,'precision',10)
        
        % height zref from literature
        
        hei = 2*ones(rows,cols);
        diesh = 5.33*ones(rows,cols);
        rccmin = 200*ones(rows,cols);
        lttpe = 4*ones(rows,cols);
        
        va1 = reshape(hei,rows*cols,1);
        va2 = reshape(diesh,rows*cols,1);
        va3 = reshape(rccmin,rows*cols,1);
        va4 = reshape(lttpe,rows*cols,1);
        
        soii = [va1,va2,va3,va4];
        %dlmwrite(['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_fixed_900m/height_rcmin_lc.ascii'],soii,'precision',10)
        
        zref1 = 0.8*ones(rows,cols);
        zref2 = 0.8*ones(rows,cols);
        zref3 = 0.8*ones(rows,cols);
        
        z1 = reshape(zref1,rows*cols,1);
        z2 = reshape(zref2,rows*cols,1);
        z3 = reshape(zref3,rows*cols,1);
        
        soii = [z1,z2,z3];
        %dlmwrite(['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_fixed_900m/zref.ascii'],soii,'precision',10)
    end
end

%% calculate flow_distance to nearest stream pixel and output

home_alps = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Alps/input_output/';
home_andes = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Andes/input_output/';
home_himalayas = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Himalayas/input_output/';


for mountain = [1 3 2]
    
    if mountain == 1
        home_moun = home_alps;mtname = 'Alps';
        
    end
    if mountain == 2
        home_moun = home_andes;mtname = 'Andes';
        
    end
    if mountain == 3
        home_moun = home_himalayas;mtname = 'Himalayas';
        
    end
    tmp = load([home_moun 'gauges_loc/Basins_dimention.mat'],'Basins_sum');
    Basins_sum_region = tmp.Basins_sum;
    kou = 1;
    for bid = 1:size(Basins_sum_region,1)
        kou = kou+1;
        if mod(kou,10)==0
            close all
        end
        
        c = Basins_sum_region(bid,2);
        r = Basins_sum_region(bid,1);
        gid = Basins_sum_region(bid,6);
        
        
        tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid,'%2.2d') '_ind.mat'],'ind');
        ind = tmp.ind;
        tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid,'%2.2d') '_intt.mat'],'intt');
        intt = tmp.intt;
        
        ffnm=[home_moun 'fortran/Basin' num2str(gid) '/countfile.out'];
        fid=fopen(ffnm);
        facc=fscanf(fid,'%d',[c,r]);
        fclose(fid);
        facc(ind) = 100;
        faccu = facc';
        ffnm=[home_moun 'fortran/Basin' num2str(gid) '/fdr.bin'];
        fid=fopen(ffnm);
        tmpdata=fread(fid,'ubit16');
        fclose(fid);
        fdr = reshape(tmpdata,c,r);
        mask = zeros(c,r);
        mask(intt) = 1;
        fdrfull = fdr';
        
        
        
        figure
        imagesc(faccu)
        
        dis_outlet = zeros(c,r);
        dis_stream = 500*ones(c,r);
        dis_stream(find(facc>5)) = 0;
        dis_stream(ind) = nan;
        for i = 1:c
            for j = 1:r
                tmp1 = i;
                tmp2 = j;
                tmpz = fdr(tmp1,tmp2);
                if mask(i,j)==0
                    continue
                end
                distance = 0;
                for itr = 1:1000000
                    
                    if tmpz == 1
                        ik = tmp1+1;
                        jk = tmp2-1;
                    elseif  tmpz == 2
                        ik = tmp1+1;
                        jk = tmp2;
                    elseif  tmpz == 4
                        ik = tmp1+1;
                        jk = tmp2+1;
                    elseif  tmpz == 8
                        ik = tmp1;
                        jk = tmp2+1;
                    elseif  tmpz == 16
                        ik = tmp1-1;
                        jk = tmp2+1;
                    elseif  tmpz == 32
                        ik = tmp1-1;
                        jk = tmp2;
                    elseif  tmpz == 64
                        ik = tmp1-1;
                        jk = tmp2-1;
                    elseif  tmpz == 128
                        ik = tmp1;
                        jk = tmp2-1;
                    end
                    tmpz = fdr(ik,jk);
                    tmp1 = ik;
                    tmp2 = jk;
                    distance = distance+1;
                    if facc(ik,jk)>5
                        dis_stream(i,j) = min(dis_stream(i,j),distance);
                    end
                    
                    if mask(ik,jk)==0
                        dis_outlet(i,j) = distance;
                        break
                    end
                    
                    
                end
            end
        end
        
        save(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/distance/Basin' num2str(gid) 'dis_stream.mat'],'dis_stream')
        save(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/distance/Basin' num2str(gid) 'dis_outlet.mat'],'dis_outlet')
        
        figure
        imagesc(dis_outlet')
        colorbar
        figure
        t = tiledlayout(1,1,'Padding','none');
        t.Units = 'inches';
        t.OuterPosition = [0.25 0.25 3.5 2.8];
        nexttile;
        imagesc(dis_stream')
        colorbar
        
    end
    
end

%% prepare inputs (First row, northmost, first column, westmost)
% read emissivity from raw data 8day MODIS 0.05 degree 5km
for Edays = 1:8:365
    
    Emis_avg = zeros(3600,7200);
    for iyear = 2007:2013
        
        fdirect = '/shared/barros-srv-05/export/RAW-DATA/jt85/HMT-SE/Data/MODIS/';
        subdirect = 'LST-MOD11C2-5600m8day/';
        yearind = [num2str(iyear) '/'];
        filenm = ['MOD11C2.A' num2str(iyear) num2str(Edays,'%3.3d') '.*.hdf'];
        
        file = strtrim(ls([fdirect subdirect yearind filenm]));% strtrim remove extra 'newline' character
        E29 = hdfread(file, ...
            'MODIS_8DAY_0.05DEG_CMG_LST', 'Fields', 'Emis_29');
        E31 = hdfread(file, ...
            'MODIS_8DAY_0.05DEG_CMG_LST', 'Fields', 'Emis_31');
        E32 = hdfread(file, ...
            'MODIS_8DAY_0.05DEG_CMG_LST', 'Fields', 'Emis_32');
        
        Emissivity = 0.0139*(double(E29)*0.002+0.49)+0.4606*(double(E31)*0.002+0.49)+0.5256*(double(E32)*0.002+0.49);
        %Jin and Liang 2006 formula
        
        Emis_avg = Emis_avg+Emissivity;
    end
    Emis_avg = Emis_avg/7;
    mkdir(['/shared/dondo/home/ml423/world_mts/LULC/Emissivity/Day' num2str(Edays,'%3.3d') '/'])
    save(['/shared/dondo/home/ml423/world_mts/LULC/Emissivity/Day' num2str(Edays,'%3.3d') '/Emissivity_avg.mat'],'Emis_avg')
    Edays
    if Edays<100
       figure
       imagesc(Emis_avg)
       colorbar
    end
end    
    
%% Read LAI global monthly averaged LAI 0.25 degrees 25km
ncdisp('/shared/dondo/home/ml423/HMIOP/NoDAre/LA/LULC/LAI/LAI_mean_monthly_1981-2015.nc4')

LAI = ncread('/shared/dondo/home/ml423/HMIOP/NoDAre/LA/LULC/LAI/LAI_mean_monthly_1981-2015.nc4','LAI');
LAI_lat = ncread('/shared/dondo/home/ml423/HMIOP/NoDAre/LA/LULC/LAI/LAI_mean_monthly_1981-2015.nc4','lat');
LAI_time = ncread('/shared/dondo/home/ml423/HMIOP/NoDAre/LA/LULC/LAI/LAI_mean_monthly_1981-2015.nc4','time');

for monthind = 1:12
    LAI_avg = rot90(LAI(:,:,monthind));% convert to first row is northmost, first column is westmost
    
    mkdir(['/shared/dondo/home/ml423/world_mts/LULC/LAI/Month' num2str(monthind,'%2.2d') '/'])
    save(['/shared/dondo/home/ml423/world_mts/LULC/LAI/Month' num2str(monthind,'%2.2d') '/LAI_avg.mat'],'LAI_avg')
    
    CV_avg = 1-exp(-0.5*LAI_avg);
    mkdir(['/shared/dondo/home/ml423/world_mts/LULC/CV/Month' num2str(monthind,'%2.2d') '/'])
    save(['/shared/dondo/home/ml423/world_mts/LULC/CV/Month' num2str(monthind,'%2.2d') '/CV_avg.mat'],'CV_avg')

end

%% Read MODIS 0.05 Albedo daily data
fdirect = '/shared/dondo/home/ml423/HMIOP/NoDAre/LA/MODIS/';
subdirect = 'Albedo_5km/';
fils_all = ls([fdirect subdirect 'MCD*']);
cc = strsplit(fils_all);cc = cc';
yearday = [];
for i = 1:length(cc)-1
    yearday(i,1) = str2num(cc{i}(68:71));
    yearday(i,2) = str2num(cc{i}(72:74));
end
%{
for Adays = 1:365
    aloc = find(yearday(:,2)==Adays);
    yearloop = yearday(aloc,1);
    
    Albedo_avg = zeros(3600,7200);
    ko = 0;
    for j = 1:length(yearloop)
        iyear = yearloop(j);
        
        filenm = ['MCD43C3.A' num2str(iyear) num2str(Adays,'%3.3d') '.*.hdf'];
%         fid = fopen([fdirect subdirect filenm]);
%         if fid == -1
%            continue 
%         end

        file = strtrim(ls([fdirect subdirect filenm]));% strtrim remove extra 'newline' character
        file_exist = strfind(fils_all,file);
        if isempty(file_exist)% meaning no such file exists
            continue
        else
            ko = ko+1;
            file = strtrim(ls([fdirect subdirect filenm]));% strtrim remove extra 'newline' character
            Albedo = hdfread(file, ...
                'MCD_CMG_BRDF_0.05Deg', 'Fields', 'Albedo_WSA_shortwave');
        end
        Albedo = double(Albedo)*0.001;
        Albedo(Albedo>0.9) = 0.9;
        Albedo_avg = Albedo_avg+Albedo;
    end
    Albedo_avg = Albedo_avg/ko;
    mkdir(['/shared/dondo/home/ml423/world_mts/LULC/Albedo/Day' num2str(Adays,'%3.3d') '/'])
    save(['/shared/dondo/home/ml423/world_mts/LULC/Albedo/Day' num2str(Adays,'%3.3d') '/Albedo_avg.mat'],'Albedo_avg')
    Adays
    disp(ko)

    if mod(Adays,10)==0
       figure
       imagesc(Albedo_avg)
       colorbar
    end
end    
%}
% Create LULC inputs for each basin in the world
% Write data into basin input folders
%

t0 = datetime(1950,1,1,0,0,0);
te = datetime(2024,5,1,23,0,0);
clear days
tt = t0:days(1):te;
[yt,mt,dt] = ymd(tt);
ounm={'LAI','CV','Emissivity','Albedo',...
        'AirPressure_10m','AirTemp_10m','SpecHumi_10m','WindSpeed_10m','DSWR_surf','DLWR_surf',};
tdays = [];
mdays = [];mvec = [31,28,31,30,31,30,31,31,30,31,30,31];mvecn = [31,29,31,30,31,30,31,31,30,31,30,31];
ydays = [];
for i = 1950:2024
    if mod(i,4)==0
        tdays = [tdays;[1:366]'];
        ydays = [ydays;366];
    else
        tdays = [tdays;[1:365]'];
        ydays = [ydays;365];
    end
    
    for j = 1:12
        if mod(i,4)==0
            mdays = [mdays;j*ones(mvecn(j),1)];
        else
            mdays = [mdays;j*ones(mvec(j),1)];
        end
    end
end

mmdd = [];
mmonly = [];
for j = 1:12
    for k = 1:mvec(j)
        mmdd = [mmdd;j*100+k];
        mmonly = [mmonly;j];
    end
end


home_alps = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Alps/input_output/';
home_andes = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Andes/input_output/';
home_himalayas = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Himalayas/input_output/';

for mountain = 2%[1 2 3]
    
    if mountain == 1
        home_moun = home_alps;mtname = 'Alps';
    end
    if mountain == 2
        home_moun = home_andes;mtname = 'Andes';
    end
    if mountain == 3
        home_moun = home_himalayas;mtname = 'Himalayas';
    end
    tmp = load([home_moun 'gauges_loc/Basins_dimention.mat'],'Basins_sum');
    Basins_sum_region = tmp.Basins_sum;
    bnum = size(Basins_sum_region,1);
    for bid = 1:bnum
    % or parfor
        cols = Basins_sum_region(bid,2);
        rows = Basins_sum_region(bid,1);
        gid = Basins_sum_region(bid,6);
        
        tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_ind.mat'],'ind');
        ind = tmp.ind;
        tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_intt.mat'],'intt');
        intt = tmp.intt;
        
        outhomedir = ['/shared/dondo/home/ml423/world_mts/'];
        mkdir([outhomedir 'Basin' num2str(gid) '_ERA5_900m1hr'])
        mkdir([outhomedir 'Basin' num2str(gid) '_MODIS_900m1hr'])
        
        lat = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/latlon/Basin' num2str(gid) '_lat.mat'],'lat_b');
        basin_lat = lat.lat_b;
        lon = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/latlon/Basin' num2str(gid) '_lon.mat'],'lon_b');
        basin_lon = lon.lon_b;
        
        lat1d = basin_lat(:);lon1d = basin_lon(:);% one can adjust the lat lon here to get nearby lat lon if nan is present
        
        
        % grids are 3600x7200, 0.05 degree for emissivity albedo
        latE = (90-0.025):-0.05:-90;
        lonE = (-180+0.025):0.05:180;
        latEa = [90 latE(1:3599)];
        lonEa = [-180 lonE(1:7199)];
        loc_Emis = [];
        for k = 1:length(lat1d)
            tmp1 = find(latE<lat1d(k)&latEa>lat1d(k));
            tmp2 = find(lonE>lon1d(k)&lonEa<lon1d(k));
            loc_Emis(k) = (tmp2-1)*3600+tmp1;
        end
        
        % grids are 720x1440, 0.25 degree for LAI CV
            latE = (90-0.125):-0.25:-90;
            lonE = (-180+0.125):0.25:180;
            latEa = [90 latE(1:719)];
            lonEa = [-180 lonE(1:1439)];
            loc_LAI = [];
            for k = 1:length(lat1d)
                tmp1 = find(latE<lat1d(k)&latEa>lat1d(k));
                tmp2 = find(lonE>lon1d(k)&lonEa<lon1d(k));
                loc_LAI(k) = (tmp2-1)*720+tmp1;
            end

        
        disp(gid)
        %for i = 1:length(yt)
        for i = 1:365
            disp(mmdd(i))
            
            out_dir = [outhomedir 'Basin' num2str(gid) '_MODIS_900m1hr/' num2str(mmdd(i),'%4.4d') '/'];
            
            mkdir(out_dir)
            %{
            %now_day = tdays(i);
            
            for updays = 1:8
                file_day = i+updays-1;
                ffnm = (['/shared/dondo/home/ml423/world_mts/LULC/Emissivity/Day' num2str(file_day,'%3.3d') '/Emissivity_avg.mat']);
                fid = fopen(ffnm);
                if fid==-1
                    continue
                else
                    tmp = load(['/shared/dondo/home/ml423/world_mts/LULC/Emissivity/Day' num2str(file_day,'%3.3d') '/Emissivity_avg.mat'],'Emis_avg');
                    tmp = tmp.Emis_avg;
                    fclose(fid);
                    break
                end
            end
            
            input = tmp(loc_Emis);
            input = reshape(input,[size(basin_lat,1), size(basin_lat,2)]);
            data3D = [];
            
            for timi = 1:24
                data3D(:,:,timi)=input';
            end
            
            fnm=[out_dir,ounm{3}];
            predata=reshape(data3D,[size(basin_lat,2),size(basin_lat,1)*24]);
            fid = fopen(fnm,'wb');
            fwrite(fid,predata,'single');
            fclose(fid);
            
            fid=fopen([fnm,'_sta.txt'],'w');
            fprintf(fid,'Min    = %.4f\n',min(predata(:)));
            fprintf(fid,'Max    = %.4f\n',max(predata(:)));
            fprintf(fid,'Mean   = %.4f\n',mean(predata(:)));
            fprintf(fid,'Median = %.4f\n',median(predata(:)));
            fprintf(fid,'Std.   = %.4f\n',std(predata(:)));
            fclose(fid);
            
            if i==1
                figure
                imagesc(input)
                colorbar
                title(['Basin ' num2str(gid) 'Emissivity'])
                
            end
            
            %}
            % LAI
            
            
            tmp = load(['/shared/dondo/home/ml423/world_mts/LULC/LAI/Month' num2str(mmonly(i),'%2.2d') '/LAI_avg.mat'],'LAI_avg');
            tmp = tmp.LAI_avg;
          

            input = tmp(loc_LAI);
            input = reshape(input,[size(basin_lat,1), size(basin_lat,2)]);
            data3D = [];

            prob_loc = find(isnan(input)|input<0);% these 3 steps are new 10/26/2024 critical otherwise 0 flow or nan flow
            good_loc = find(input>0);
            input(prob_loc) = mean(input(good_loc));
            
            for timi = 1:24
                data3D(:,:,timi)=input';
            end
            
            fnm=[out_dir,ounm{1}];
            predata=reshape(data3D,[size(basin_lat,2),size(basin_lat,1)*24]);
            fid = fopen(fnm,'wb');
            fwrite(fid,predata,'single');
            fclose(fid);
            
            fid=fopen([fnm,'_sta.txt'],'w');
            fprintf(fid,'Min    = %.4f\n',min(predata(:)));
            fprintf(fid,'Max    = %.4f\n',max(predata(:)));
            fprintf(fid,'Mean   = %.4f\n',mean(predata(:)));
            fprintf(fid,'Median = %.4f\n',median(predata(:)));
            fprintf(fid,'Std.   = %.4f\n',std(predata(:)));
            fclose(fid);
            
            if i==1
                figure
                imagesc(input)
                colorbar
                title(['Basin ' num2str(gid) 'LAI'])
                
            end
            
            
            % CV
            
            tmp = load(['/shared/dondo/home/ml423/world_mts/LULC/CV/Month' num2str(mmonly(i),'%2.2d') '/CV_avg.mat'],'CV_avg');
            tmp = tmp.CV_avg;
         

            input = tmp(loc_LAI);
            input = reshape(input,[size(basin_lat,1), size(basin_lat,2)]);
            data3D = [];

            prob_loc = find(isnan(input)|input<0);
            good_loc = find(input>0);
            input(prob_loc) = mean(input(good_loc));
            
            for timi = 1:24
                data3D(:,:,timi)=input';
            end
            
            fnm=[out_dir,ounm{2}];
            predata=reshape(data3D,[size(basin_lat,2),size(basin_lat,1)*24]);
            fid = fopen(fnm,'wb');
            fwrite(fid,predata,'single');
            fclose(fid);
            
            fid=fopen([fnm,'_sta.txt'],'w');
            fprintf(fid,'Min    = %.4f\n',min(predata(:)));
            fprintf(fid,'Max    = %.4f\n',max(predata(:)));
            fprintf(fid,'Mean   = %.4f\n',mean(predata(:)));
            fprintf(fid,'Median = %.4f\n',median(predata(:)));
            fprintf(fid,'Std.   = %.4f\n',std(predata(:)));
            fclose(fid);
            
            if i==1
                figure
                imagesc(input)
                colorbar
                title(['Basin ' num2str(gid) 'CV'])
                
            end
            %
            %{
            % Albedo
            tmp = load(['/shared/dondo/home/ml423/world_mts/LULC/Albedo/Day' num2str(i,'%3.3d') '/Albedo_avg.mat'],'Albedo_avg');
            tmp = tmp.Albedo_avg;
    
            
            input = tmp(loc_Emis);
            input = reshape(input,[size(basin_lat,1), size(basin_lat,2)]);
            data3D = [];
            
            for timi = 1:24
                data3D(:,:,timi)=input';
            end
            
            fnm=[out_dir,ounm{4}];
            predata=reshape(data3D,[size(basin_lat,2),size(basin_lat,1)*24]);
            fid = fopen(fnm,'wb');
            fwrite(fid,predata,'single');
            fclose(fid);
            
            fid=fopen([fnm,'_sta.txt'],'w');
            fprintf(fid,'Min    = %.4f\n',min(predata(:)));
            fprintf(fid,'Max    = %.4f\n',max(predata(:)));
            fprintf(fid,'Mean   = %.4f\n',mean(predata(:)));
            fprintf(fid,'Median = %.4f\n',median(predata(:)));
            fprintf(fid,'Std.   = %.4f\n',std(predata(:)));
            fclose(fid);
            
            if i==1
                figure
                imagesc(input)
                colorbar
                title(['Basin ' num2str(gid) 'albedo'])
            end            
            %}

        end
        out_dircopy = [outhomedir 'Basin' num2str(gid) '_MODIS_900m1hr/' num2str(228,'%4.4d') '/'];
        out_dir = [outhomedir 'Basin' num2str(gid) '_MODIS_900m1hr/' num2str(229,'%4.4d') '/'];
        mkdir(out_dir)
        
        copyfile([out_dircopy '*'],out_dir)
        
        %{
    for v=1:4
        for k=1:days
            daily_dir=[data_dir,'/',int2str(doy+k-1),'/'];
            filename=[daily_dir,ounm{v},'.bin'];
            fid = fopen(filename,'rb','ieee-le');
            tmpdata = fread(fid,inf,'single');
            fclose(fid);
            data = reshape(tmpdata,ncol1km,nrow1km,tstep);
            predatatmp(:,:,(k-1)*tstep+1:k*tstep)=data(c1:c2,r1:r2,:);
        end
        %---------------------------------
        ttstep=tstep*days;data3D=NaN(col*4,row*4,(ttstep-houroff-hourini)*12);
        predata1hr=reshape(predatatmp,[col,row*ttstep]);
        
        predata30min=tohalfstep(predata1hr,row,col,ttstep);
        predata15min=tohalfstep(predata30min,row,col,ttstep*2);
        predata5min=tothirdstep(predata15min,row,col,ttstep*4);%1km5min
        
        data1km=reshape(predata5min,[col,row,ttstep*12]);
        for t=1+hourini*12:(ttstep-houroff)*12
            tmpp=data1km(:,:,t)';%row,col
            for i=1:row
                for j=1:col
                    tmp(4*i-3:4*i,4*j-3:4*j)=tmpp(i,j);
                end
            end
            data3D(:,:,t-hourini*12)=tmp';
        end
        
        fnm=[out_dir,ounm{v}];
        
        predata=reshape(data3D,[col*4,row*4*(ttstep-houroff-hourini)*12]);
        fid = fopen(fnm,'wb');
        fwrite(fid,predata,'single');
        fclose(fid);
        
        fid=fopen([fnm,'_sta.txt'],'w');
        fprintf(fid,'Min    = %.4f\n',min(predata(:)));
        fprintf(fid,'Max    = %.4f\n',max(predata(:)));
        fprintf(fid,'Mean   = %.4f\n',mean(predata(:)));
        fprintf(fid,'Median = %.4f\n',median(predata(:)));
        fprintf(fid,'Std.   = %.4f\n',std(predata(:)));
        fclose(fid);
    end
    
    
    for v=5:length(ounm)
        for k=1:days
            daily_dir=[data_dir,'/',int2str(doy+k-1),'/'];
            filename=[daily_dir,ounm{v},'.bin'];
            fid = fopen(filename,'rb','ieee-le');
            tmpdata = fread(fid,inf,'single');
            fclose(fid);
            data = reshape(tmpdata,ncol1km,nrow1km,tstep);
            predatatmp(:,:,(k-1)*tstep+1:k*tstep)=data(c1:c2,r1:r2,:);
        end
        %---------------------------------
        ttstep=tstep*days;data3D=NaN(col*4,row*4,(ttstep-houroff-hourini)*12);
        predata1hr=reshape(predatatmp,[col,row*ttstep]);
        
        predata30min=tohalfstep(predata1hr,row,col,ttstep);
        predata15min=tohalfstep(predata30min,row,col,ttstep*2);
        predata5min=tothirdstep(predata15min,row,col,ttstep*4);%1km5min
        
        data1km=reshape(predata5min,[col,row,ttstep*12]);
        for t=1+hourini*12:(ttstep-houroff)*12
            tmpp=data1km(:,:,t)';%row,col
            for i=1:row
                for j=1:col
                    tmp(4*i-3:4*i,4*j-3:4*j)=tmpp(i,j);
                end
            end
            data3D(:,:,t-hourini*12)=tmp';
        end
        
        fnm=[out_dir2,ounm{v}];
        
        predata=reshape(data3D,[col*4,row*4*(ttstep-houroff-hourini)*12]);
        fid = fopen(fnm,'wb');
        fwrite(fid,predata,'single');
        fclose(fid);
        
        fid=fopen([fnm,'_sta.txt'],'w');
        fprintf(fid,'Min    = %.4f\n',min(predata(:)));
        fprintf(fid,'Max    = %.4f\n',max(predata(:)));
        fprintf(fid,'Mean   = %.4f\n',mean(predata(:)));
        fprintf(fid,'Median = %.4f\n',median(predata(:)));
        fprintf(fid,'Std.   = %.4f\n',std(predata(:)));
        fclose(fid);
    end
        %}
    end
    %
    
end

%% read ERA5 and reformat data to basin level, day 1 to day 28

% when directly read grib, using nctoolbox-1.1.3 installed at /ml423/DA/
warning('off','all')
addpath('/shared/dondo/home/ml423/DA/nctoolbox-1.1.3')
setup_nctoolbox

lat = [90:-0.1:-90]';
lon = [-360:0.1:-0.1]';

ffnm=['/shared/dondo/home/ml423/ERA5/ERA5_data/data2/2019_01_Tdew.grib'];
gf = ncgeodataset(ffnm);
lat = gf.geovariable('lat');
lat = lat.data(:);
lon = gf.geovariable('lon');
lon = lon.data(:);
lon = lon-360;

home_alps = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Alps/input_output/';
home_andes = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Andes/input_output/';
home_himalayas = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Himalayas/input_output/';

varibs = {'T','P','DS','DL','Pres','wind_v','wind_u','Tdew'};%Downloaded ERA5 namings
fullname = {'2_metre_temperature_surface',
    'Total_precipitation_surface_1_Hour_Accumulation',
    'Surface_solar_radiation_surface_1_Hour_Accumulation',
    'Surface_thermal_radiation_surface_1_Hour_Accumulation',
    'Surface_pressure_surface',
    '10_metre_V_wind_component_surface',
    '10_metre_U_wind_component_surface',
    '2_metre_dewpoint_temperature_surface'};% ERA5 built-in default namings

ounmf = {'AirTemp_10m','Precipitation_ERA5_land_ori','DSWR_surf','DLWR_surf','AirPressure_10m','WindV_10m','WindU_10m','AirTempDew_10m'};


% based on selected events and basins, to create inputs, as this step takes
% a long time, so not every year, every month is systemeatically generated,
% but rather conditional on gauges and events selection
wge = load(['/shared/dondo/home/ml423/DA/Critical matrices/world_selected_gauge_events.mat'],'selected_gauge_events');
wge = wge.selected_gauge_events;
% for each gauge find the first and the last event.
gauge_events_first_last = [];
for i = 1:size(wge,1)
   for j = 1:size(wge,2)
       if length(wge{i,j})<=1% no events qualified for this gauge
       else
           gid = str2double(wge{i,j}(1));
           evets = (char(wge{i,j}(2:end)));
           dasum = [];
           for k = 1:length(wge{i,j})-1
               datetmp = evets(:,:,k);
               datenber = str2num(datetmp(1:4))*10000+str2num(datetmp(6:7))*100+str2num(datetmp(9:10));
               dasum = [dasum;datenber];
           end
           gauge_events_first_last =[gauge_events_first_last; [gid,min(dasum),max(dasum),i]];
       end
   end
end
%clear wge
%clear gauge_events_first_last

problem_month = [];
var = 3;
any_files_missing_left = [1980,6,3;
    1983,12,1;
    1984,3,2;
    1984,4,6;
    1984,4,8;
    1986,1,4;
    1990,12,1;
    1992,7,3;
    1992,7,4;
    1992,7,5;
    1994,3,6;
    1999,11,7;
    1999,12,3;
    1999,12,4;
    1999,12,5;
    1986,2,7;
    2017,2,7];
for lop = 16:size(any_files_missing_left,1)
    var = any_files_missing_left(lop,3);
for iyear = any_files_missing_left(lop,1)
    problem_day = [];
    for imonth = any_files_missing_left(lop,2)
        
        ffnm=['/shared/dondo/home/ml423/ERA5/ERA5_data/data2/' num2str(iyear) '_' num2str(imonth,'%2.2d') '_' varibs{var} '.grib'];
        fid = fopen(ffnm);
        
        if fid == -1
            continue
        end
        
        % representation gauge is 2969435, which is a gauge in himalayas, 
        % meaning if this basin has inputs written already, then no need to re write input for this month for all gauges in the world.
%         repgauge = '2969435';
%         ffnmfolder=['/shared/dondo/home/ml423/world_mts/Basin' repgauge '_ERA5_900m1hr/' num2str(iyear) num2str(imonth,'%2.2d') '28/'];
%         
%         if exist(ffnmfolder, 'dir')
%             
%             ffnmf=['/shared/dondo/home/ml423/world_mts/Basin' repgauge '_ERA5_900m1hr/' num2str(iyear) num2str(imonth,'%2.2d') '28/' ounmf{var}];
%             fidf = fopen(ffnmf);
%             if fidf ~= -1
%                 % the inputs for all basins for this month is written already
%                 disp([num2str(iyear) num2str(imonth,'%2.2d') ' is already finished'])
%                 fclose(fidf);
%                 continue
%             end
%         end
        
        
        
        gf = ncgeodataset(ffnm);
        vartmp = gf.geovariable(fullname{var});dime = size(vartmp);
        varmatrix = {};
        if dime(1)<672
            problem_month = [problem_month;[iyear,imonth,var]];
            continue
        end
        
        for vv =1:672
            varmatrix{vv} = squeeze(vartmp(vv,:,:));
        end
        
        %inum = iyear*10000+imonth*100;
        

        for mountain = [1 2 3]
            
            if mountain == 1
                home_moun = home_alps;mtname = 'Alps';
                
                repgauge = '6343531';
                ffnmfolder=['/shared/dondo/home/ml423/world_mts/Basin' repgauge '_ERA5_900m1hr/' num2str(iyear) num2str(imonth,'%2.2d') num2str(28,'%2.2d') '/'];
                
                fid = fopen([ffnmfolder ounmf{var}],'rb','ieee-le');
                if fid==-1
                    % file already exists, no need to rename it.
                    %problem_day = [problem_day;[iyear,imonth,1,var]];
                else
                    fclose(fid);
                    disp('Alps is good')
                    continue
                end
                
            end
            if mountain == 2
                home_moun = home_andes;mtname = 'Andes';
                repgauge = '3650630';
                ffnmfolder=['/shared/dondo/home/ml423/world_mts/Basin' repgauge '_ERA5_900m1hr/' num2str(iyear) num2str(imonth,'%2.2d') num2str(28,'%2.2d') '/'];
                
                fid = fopen([ffnmfolder ounmf{var}],'rb','ieee-le');
                if fid==-1
                    % file already exists, no need to rename it.
                    %problem_day = [problem_day;[iyear,imonth,1,var]];
                else
                    fclose(fid);
                    disp('Andes is good')
                    continue
                end
            end
            if mountain == 3
                home_moun = home_himalayas;mtname = 'Himalayas';
                repgauge = '2969435';
                ffnmfolder=['/shared/dondo/home/ml423/world_mts/Basin' repgauge '_ERA5_900m1hr/' num2str(iyear) num2str(imonth,'%2.2d') num2str(28,'%2.2d') '/'];
                
                fid = fopen([ffnmfolder ounmf{var}],'rb','ieee-le');
                if fid==-1
                    % file already exists, no need to rename it.
                    %problem_day = [problem_day;[iyear,imonth,1,var]];
                else
                    fclose(fid);
                    disp('Himalayas is good')
                    continue
                end
            end
            tmp = load([home_moun 'gauges_loc/Basins_dimention.mat'],'Basins_sum');
            Basins_sum_region = tmp.Basins_sum;
            bnum = size(Basins_sum_region,1);
            
            parfor bid = 1:bnum
                disp(iyear)
                disp(imonth)
                disp(mountain)
                
                cols = Basins_sum_region(bid,2);
                rows = Basins_sum_region(bid,1);
                gid = Basins_sum_region(bid,6);
                outhomedir = ['/shared/dondo/home/ml423/world_mts/'];
                mkdir([outhomedir 'Basin' num2str(gid) '_ERA5_900m1hr'])
                
%                 g_ind = find(gauge_events_first_last(:,1)==gid);
%                 
%                 if isempty(g_ind)% this gauge has no qualified events to look at
%                     continue
%                 elseif inum>gauge_events_first_last(g_ind,3)||inum<(gauge_events_first_last(g_ind,2)-10000)
%                     % meaning current data processing date is beyond the
%                     % selected events and not needed for more processing
%                     % for this gauge
%                     continue
%                 end
%                 
%                 %the part below is problematic, can cause matlab kill for
%                 %some reason
%                 
%                 file_exist = fopen([outhomedir 'Basin' num2str(gid) '_ERA5_900m1hr/' num2str(iyear) num2str(imonth,'%2.2d') num2str(28,'%2.2d') '/' ounmf{var}]);
%                 if file_exist~=-1% meaning file already processed, continue to next gauge
%                     fclose(file_exist);
%                     continue
%                 else% meaning file not generated yet, has to be generated
%                 end
        
                tic
                
                
                blat = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/latlon/Basin' num2str(gid) '_lat.mat'],'lat_b');
                basin_lat = blat.lat_b;
                blon = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/latlon/Basin' num2str(gid) '_lon.mat'],'lon_b');
                basin_lon = blon.lon_b;
                
                lat1d = basin_lat(:);lon1d = basin_lon(:);
                
                latE = (90-0.05):-0.1:-90.1;
                lonE = [(0-0.05):0.1:180 -179.95:0.1:-0.1];
                latEa = [90 latE(1:1800)];
                lonEa = [-0.1 lonE(1:3599)];
                loc_ERA5L = [];
                for k = 1:length(lat1d)
                    tmp1 = find(latE<lat1d(k)&latEa>lat1d(k));
                    tmp2 = find(lonE>lon1d(k)&lonEa<lon1d(k));
                    loc_ERA5L(k) = (tmp2-1)*1801+tmp1;
                end
                toc
                tic
                for iday = 1:28
           
                    outhomedir = ['/shared/dondo/home/ml423/world_mts/'];
                    outfiledir = [outhomedir 'Basin' num2str(gid) '_ERA5_900m1hr/' num2str(iyear) num2str(imonth,'%2.2d') num2str(iday,'%2.2d') '/'];
                    mkdir(outfiledir)
                    data3D = [];
                    for ihour = 1:24
                        %tmp = vartmp((iday-1)*24+ihour,:,:);
                        tmp = varmatrix{(iday-1)*24+ihour};
                        input = tmp(loc_ERA5L);
                        input = reshape(input,[size(basin_lat,1), size(basin_lat,2)]);
         
                        data3D(:,:,ihour)=input';
                        
                    end

                        fnm=[outfiledir,ounmf{var}];
                        predata=reshape(data3D,[size(basin_lat,2),size(basin_lat,1)*24]);
                        fid = fopen(fnm,'wb');
                        fwrite(fid,predata,'single');
                        fclose(fid);
                        
                        fid=fopen([fnm,'_sta.txt'],'w');
                        fprintf(fid,'Min    = %.4f\n',min(predata(:)));
                        fprintf(fid,'Max    = %.4f\n',max(predata(:)));
                        fprintf(fid,'Mean   = %.4f\n',mean(predata(:)));
                        fprintf(fid,'Median = %.4f\n',median(predata(:)));
                        fprintf(fid,'Std.   = %.4f\n',std(predata(:)));
                        fclose(fid);
                    
                end
                toc
            end
        end
    end
end
end
%% read ERA5 and reformat data to basin level, day 29 to day 31

% when directly read grib, using nctoolbox-1.1.3 installed at /ml423/DA/

addpath('/shared/dondo/home/ml423/DA/nctoolbox-1.1.3')
setup_nctoolbox

lat = [90:-0.1:-90]';
lon = [-360:0.1:-0.1]';

ffnm=['/shared/dondo/home/ml423/ERA5/ERA5_data/data2/2019_01_Tdew.grib'];
gf = ncgeodataset(ffnm);
lat = gf.geovariable('lat');
lat = lat.data(:);
lon = gf.geovariable('lon');
lon = lon.data(:);
lon = lon-360;

home_alps = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Alps/input_output/';
home_andes = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Andes/input_output/';
home_himalayas = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Himalayas/input_output/';

varibs = {'T2','P2','DS2','DL2','Pres2','wind_v2','wind_u2','Tdew2'};%Downloaded ERA5 namings
fullname = {'2_metre_temperature_surface',
    'Total_precipitation_surface_1_Hour_Accumulation',
    'Surface_solar_radiation_surface_1_Hour_Accumulation',
    'Surface_thermal_radiation_surface_1_Hour_Accumulation',
    'Surface_pressure_surface',
    '10_metre_V_wind_component_surface',
    '10_metre_U_wind_component_surface',
    '2_metre_dewpoint_temperature_surface'};% ERA5 built-in default namings

ounmf = {'AirTemp_10m','Precipitation_ERA5_land_ori','DSWR_surf','DLWR_surf','AirPressure_10m','WindV_10m','WindU_10m','AirTempDew_10m'};

var = 7;
% based on selected events and basins, to create inputs, as this step takes
% a long time, so not every year, every month is systemeatically generated,
% but rather conditional on gauges and events selection
wge = load(['/shared/dondo/home/ml423/DA/Critical matrices/world_selected_gauge_events.mat'],'selected_gauge_events');
wge = wge.selected_gauge_events;
% for each gauge find the first and the last event.
gauge_events_first_last = [];
for i = 1:size(wge,1)
   for j = 1:size(wge,2)
       if length(wge{i,j})<=1% no events qualified for this gauge
       else
           gid = str2double(wge{i,j}(1));
           evets = (char(wge{i,j}(2:end)));
           dasum = [];
           for k = 1:length(wge{i,j})-1
               datetmp = evets(:,:,k);
               datenber = str2num(datetmp(1:4))*10000+str2num(datetmp(6:7))*100+str2num(datetmp(9:10));
               dasum = [dasum;datenber];
           end
           gauge_events_first_last =[gauge_events_first_last; [gid,min(dasum),max(dasum),i]];
       end
   end
end
%clear wge
%clear gauge_events_first_last

imonthdays = [31 29 31 30 31 30 31 31 30 31 30 31];
resthours = 24*[3 1 3 2 3 2 3 3 2 3 2 3];
any_files_missing_left = [1986,2,7;
    2017,2,7];
for lop = 1:size(any_files_missing_left,1)
    var = any_files_missing_left(lop,3);
for iyear = any_files_missing_left(lop,1)
    %parpool('local',12)
    % if iyear==1992
    %     qi = 1;
    % else
    %     qi = 1;
    % end
    
    for imonth = any_files_missing_left(lop,2)
        
        if imonth==2&&mod(iyear,4)~=0
           continue 
        end
        
        
        
        ffnm=['/shared/dondo/home/ml423/ERA5/ERA5_data/data2/' num2str(iyear) '_' num2str(imonth,'%2.2d') '_' varibs{var} '.grib'];
        fid = fopen(ffnm);
        
        if fid == -1
            continue
        end
        
        % representation gauge is 2969435, which is a gauge in himalayas, 
        % meaning if this basin has inputs written already, then no need to re write input for this month for all gauges in the world.
%         repgauge = '2969435';
%         ffnmfolder=['/shared/dondo/home/ml423/world_mts/Basin' repgauge '_ERA5_900m1hr/' num2str(iyear) num2str(imonth,'%2.2d') num2str(imonthdays(imonth),'%2.2d') '/'];
%         
%         if exist(ffnmfolder, 'dir')
%             
%             ffnmf=['/shared/dondo/home/ml423/world_mts/Basin' repgauge '_ERA5_900m1hr/' num2str(iyear) num2str(imonth,'%2.2d') num2str(imonthdays(imonth),'%2.2d') '/' ounmf{var}];
%             fidf = fopen(ffnmf);
%             if fidf ~= -1
%                 % the inputs for all basins for this month is written already
%                 disp([num2str(iyear) num2str(imonth,'%2.2d') ' is already finished'])
%                 fclose(fidf);
%                 continue
%             end
%         end
        
        
        
        gf = ncgeodataset(ffnm);
        vartmp = gf.geovariable(fullname{var});dime = size(vartmp);
        varmatrix = {};
        if dime(1)<resthours(imonth)% partial data downloaded only, not the full month data
            problem_month = [problem_month;[iyear,imonth,var]];
            continue
        end
        for vv =1:resthours(imonth)
            varmatrix{vv} = squeeze(vartmp(vv,:,:));
        end
        
        %inum = iyear*10000+imonth*100;
        

        for mountain = [1 2 3]
            
            if mountain == 1
                home_moun = home_alps;mtname = 'Alps';
            end
            if mountain == 2
                home_moun = home_andes;mtname = 'Andes';
            end
            if mountain == 3
                home_moun = home_himalayas;mtname = 'Himalayas';
            end
            tmp = load([home_moun 'gauges_loc/Basins_dimention.mat'],'Basins_sum');
            Basins_sum_region = tmp.Basins_sum;
            bnum = size(Basins_sum_region,1);
            
            parfor bid = 1:bnum
                disp(iyear)
                disp(imonth)
                disp(bid)
                
                cols = Basins_sum_region(bid,2);
                rows = Basins_sum_region(bid,1);
                gid = Basins_sum_region(bid,6);
                outhomedir = ['/shared/dondo/home/ml423/world_mts/'];
                mkdir([outhomedir 'Basin' num2str(gid) '_ERA5_900m1hr'])
                
%                 g_ind = find(gauge_events_first_last(:,1)==gid);
%                 
%                 if isempty(g_ind)% this gauge has no qualified events to look at
%                     continue
%                 elseif inum>gauge_events_first_last(g_ind,3)||inum<(gauge_events_first_last(g_ind,2)-10000)
%                     % meaning current data processing date is beyond the
%                     % selected events and not needed for more processing
%                     % for this gauge
%                     continue
%                 end
%                 
%                 %the part below is problematic, can cause matlab kill for
%                 %some reason
%                 
%                 file_exist = fopen([outhomedir 'Basin' num2str(gid) '_ERA5_900m1hr/' num2str(iyear) num2str(imonth,'%2.2d') num2str(28,'%2.2d') '/' ounmf{var}]);
%                 if file_exist~=-1% meaning file already processed, continue to next gauge
%                     fclose(file_exist);
%                     continue
%                 else% meaning file not generated yet, has to be generated
%                 end
        
                tic
                
                
                blat = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/latlon/Basin' num2str(gid) '_lat.mat'],'lat_b');
                basin_lat = blat.lat_b;
                blon = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/latlon/Basin' num2str(gid) '_lon.mat'],'lon_b');
                basin_lon = blon.lon_b;
                
                lat1d = basin_lat(:);lon1d = basin_lon(:);
                
                latE = (90-0.05):-0.1:-90.1;
                lonE = [(0-0.05):0.1:180 -179.95:0.1:-0.1];
                latEa = [90 latE(1:1800)];
                lonEa = [-0.1 lonE(1:3599)];
                loc_ERA5L = [];
                for k = 1:length(lat1d)
                    tmp1 = find(latE<lat1d(k)&latEa>lat1d(k));
                    tmp2 = find(lonE>lon1d(k)&lonEa<lon1d(k));
                    loc_ERA5L(k) = (tmp2-1)*1801+tmp1;
                end
                toc
                tic
                for iday = 29:(imonthdays(imonth))
           
                    outhomedir = ['/shared/dondo/home/ml423/world_mts/'];
                    outfiledir = [outhomedir 'Basin' num2str(gid) '_ERA5_900m1hr/' num2str(iyear) num2str(imonth,'%2.2d') num2str(iday,'%2.2d') '/'];
                    mkdir(outfiledir)
                    data3D = [];
                    for ihour = 1:24
                        %tmp = vartmp((iday-1)*24+ihour,:,:);
                        tmp = varmatrix{(iday-29)*24+ihour};
                        input = tmp(loc_ERA5L);
                        input = reshape(input,[size(basin_lat,1), size(basin_lat,2)]);
         
                        data3D(:,:,ihour)=input';
                        
                    end

                        fnm=[outfiledir,ounmf{var}];
                        predata=reshape(data3D,[size(basin_lat,2),size(basin_lat,1)*24]);
                        fid = fopen(fnm,'wb');
                        fwrite(fid,predata,'single');
                        fclose(fid);
                        
                        fid=fopen([fnm,'_sta.txt'],'w');
                        fprintf(fid,'Min    = %.4f\n',min(predata(:)));
                        fprintf(fid,'Max    = %.4f\n',max(predata(:)));
                        fprintf(fid,'Mean   = %.4f\n',mean(predata(:)));
                        fprintf(fid,'Median = %.4f\n',median(predata(:)));
                        fprintf(fid,'Std.   = %.4f\n',std(predata(:)));
                        fclose(fid);
                    
                end
                toc
            end
        end
    end
end
save(['/shared/dondo/home/ml423/world_mts/a_problem/problem_month.mat'],'problem_month')
end    
%% Create Wind and specific humidity, convert all variables units used by DCHM
home_alps = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Alps/input_output/';
home_andes = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Andes/input_output/';
home_himalayas = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Himalayas/input_output/';

varibs = {'T','P','DS','DL','Pres','wind_v','wind_u','Tdew'};%Downloaded ERA5 namings
fullname = {'2_metre_temperature_surface',
    'Total_precipitation_surface_1_Hour_Accumulation',
    'Surface_solar_radiation_surface_1_Hour_Accumulation',
    'Surface_thermal_radiation_surface_1_Hour_Accumulation',
    'Surface_pressure_surface',
    '10_metre_V_wind_component_surface',
    '10_metre_U_wind_component_surface',
    '2_metre_dewpoint_temperature_surface'};% ERA5 built-in default namings

ounmf = {'AirTemp_10m','Precipitation_ERA5_land_ori','DSWR_surf','DLWR_surf','AirPressure_10m','WindV_10m','WindU_10m','AirTempDew_10m'};

imonthdays_leap = [31 29 31 30 31 30 31 31 30 31 30 31];
imonthdays_normal = [31 28 31 30 31 30 31 31 30 31 30 31];
resthours = 24*[3 1 3 2 3 2 3 3 2 3 2 3];
% first let's double check every variable grib file is downloaded
% then check if every grib is converted to fortran readable already?
problem_month = [];
var = 8;
for iyear = 2021:-1:1973
    %parpool('local',12)
    if mod(iyear,4)==0
        imonthdays = imonthdays_leap;
    else
        imonthdays = imonthdays_normal;
    end
    for imonth = 1:12
        
        % representation gauge is 2969435, which is a gauge in himalayas,
        % meaning if this basin has inputs written already, then no need to re write input for this month for all gauges in the world.
%         repgauge = '2969435';
%         ffnmfolder=['/shared/dondo/home/ml423/world_mts/Basin' repgauge '_ERA5_900m1hr/' num2str(iyear) num2str(imonth,'%2.2d') num2str(imonthdays(imonth),'%2.2d') '/'];
%         
%         
%         ffnm=['/shared/dondo/home/ml423/ERA5/ERA5_data/data2/' num2str(iyear) '_' num2str(imonth,'%2.2d') '_' varibs{var} '.grib'];
%         fid = fopen(ffnm);
%         
%         if fid == -1
%             problem_month = [problem_month;[iyear,imonth,var]];
%             disp('no grib 28')
%             continue
%         end
%         
%         gf = ncgeodataset(ffnm);
%         vartmp = gf.geovariable(fullname{var});dime = size(vartmp);
%         varmatrix = {};
%         if dime(1)<672
%             problem_month = [problem_month;[iyear,imonth,var]];
%             disp('grib cut off 28')
%             continue
%         end
%         
%         
%         ffnm=['/shared/dondo/home/ml423/ERA5/ERA5_data/data2/' num2str(iyear) '_' num2str(imonth,'%2.2d') '_' varibs{var} '2.grib'];
%         fid = fopen(ffnm);
%         
%         if fid == -1
%             problem_month = [problem_month;[iyear,imonth,var]];
%             disp('no grib 31')
%             continue
%         end
%         
%         gf = ncgeodataset(ffnm);
%         vartmp = gf.geovariable(fullname{var});dime = size(vartmp);
%         varmatrix = {};
%         if dime(1)<resthours
%             problem_month = [problem_month;[iyear,imonth,var]];
%             disp('grib cut off 31')
%             continue
%         end
        
        
        
        
        
        %         if exist(ffnmfolder, 'dir')
        %
        %             ffnmf=['/shared/dondo/home/ml423/world_mts/Basin' repgauge '_ERA5_900m1hr/' num2str(iyear) num2str(imonth,'%2.2d') num2str(imonthdays(imonth),'%2.2d') '/WindSpeed_10m'];
        %             fidfw = fopen(ffnmf);
        %             ffnmf=['/shared/dondo/home/ml423/world_mts/Basin' repgauge '_ERA5_900m1hr/' num2str(iyear) num2str(imonth,'%2.2d') num2str(imonthdays(imonth),'%2.2d') '/SpecHumi_10m'];
        %             fidfs = fopen(ffnmf);
        %             if fidfw ~= -1 && fidfs ~= -1
        %                 % the inputs for all basins for this month is written already
        %                 disp([num2str(iyear) num2str(imonth,'%2.2d') ' is already finished'])
        %                 fclose(fidfw);
        %                 fclose(fidfs);
        %                 continue
        %             end
        %         end
        
        
        for mountain = [2]
            
            if mountain == 1
                home_moun = home_alps;mtname = 'Alps';
            end
            if mountain == 2
                home_moun = home_andes;mtname = 'Andes';
            end
            if mountain == 3
                home_moun = home_himalayas;mtname = 'Himalayas';
            end
            tmp = load([home_moun 'gauges_loc/Basins_dimention.mat'],'Basins_sum');
            Basins_sum_region = tmp.Basins_sum;
            bnum = size(Basins_sum_region,1);
            
            parfor bid = 1:bnum
                disp(iyear)
                disp(imonth)
                disp(bid)
                
                cols = Basins_sum_region(bid,2);
                rows = Basins_sum_region(bid,1);
                gid = Basins_sum_region(bid,6);
                for iday = 1:imonthdays(imonth)
                   
                    currentday = datetime(iyear,imonth,iday);
                    nextday = currentday+days(1);
                    [ny,nm,nd] = ymd(nextday);
                    outhomedir = ['/shared/dondo/home/ml423/world_mts/'];
                    outfiledir = [outhomedir 'Basin' num2str(gid) '_ERA5_900m1hr/' num2str(iyear) num2str(imonth,'%2.2d') num2str(iday,'%2.2d') '/'];
                    outfilenextdir = [outhomedir 'Basin' num2str(gid) '_ERA5_900m1hr/' num2str(ny) num2str(nm,'%2.2d') num2str(nd,'%2.2d') '/'];
                    fid = fopen([outfiledir 'WindV_10m'],'rb','ieee-le');
                    tmpdata=fread(fid,inf,'single');fclose(fid);
                    windv=reshape(tmpdata,cols,rows,24);
                    fid = fopen([outfiledir 'WindU_10m'],'rb','ieee-le');
                    tmpdata=fread(fid,inf,'single');fclose(fid);
                    windu=reshape(tmpdata,cols,rows,24);
                    
                    
                    fid = fopen([outfiledir 'AirTempDew_10m'],'rb','ieee-le');
                    tmpdata=fread(fid,inf,'single');fclose(fid);
                    tdew=reshape(tmpdata,cols,rows,24);
                    fid = fopen([outfiledir 'AirTemp_10m'],'rb','ieee-le');
                    tmpdata=fread(fid,inf,'single');fclose(fid);
                    tair=reshape(tmpdata,cols,rows,24);
                    fid = fopen([outfiledir 'AirPressure_10m'],'rb','ieee-le');
                    tmpdata=fread(fid,inf,'single');fclose(fid);
                    presr=reshape(tmpdata,cols,rows,24);
                    fid = fopen([outfiledir 'DLWR_surf'],'rb','ieee-le');
                    tmpdata=fread(fid,inf,'single');fclose(fid);
                    dlwr=reshape(tmpdata,cols,rows,24);
                    fid = fopen([outfiledir 'DSWR_surf'],'rb','ieee-le');
                    tmpdata=fread(fid,inf,'single');fclose(fid);
                    dswr=reshape(tmpdata,cols,rows,24);
                    fid = fopen([outfiledir 'Precipitation_ERA5_land_ori'],'rb','ieee-le');
                    tmpdata=fread(fid,inf,'single');fclose(fid);
                    preci=reshape(tmpdata,cols,rows,24);
                    
                    %SH = 0.622*611./presr.*exp(17.27*(tdew-273.15)./(237.3+tdew-273.15));
                    %presrmb = presr/100;
                    %wind = sqrt(windv^2+windu^2);
                    %precip = preci.*1000;
                    
                    
                    new_wind = [];
                    new_p = [];
                    new_sh = [];
                    new_dlwr = [];
                    new_dswr = [];
                    
                    
                    for ihour = 1:24
                        windui = windu(:,:,ihour);
                        windvi = windv(:,:,ihour);
                        presri = presr(:,:,ihour);
                        tdewi = tdew(:,:,ihour);
                        SH = 0.622*611./presri.*exp(17.27*(tdewi-273.15)./(237.3+tdewi-273.15));
                        presrimb = presri/100;
                        wind = sqrt(windvi.^2+windui.^2);
                        if ihour == 1
                            precipi = preci(:,:,ihour+1);
                            dlwri = dlwr(:,:,ihour+1);
                            dswri = dswr(:,:,ihour+1);
                        end
                        
                        if ihour > 1 && ihour <= 23
                            precipi = preci(:,:,ihour+1)-preci(:,:,ihour);
                            dlwri = dlwr(:,:,ihour+1)-dlwr(:,:,ihour);
                            dswri = dswr(:,:,ihour+1)-dswr(:,:,ihour);
                        end
                        if ihour == 24
                            fid = fopen([outfilenextdir 'Precipitation_ERA5_land_ori'],'rb','ieee-le');
                            tmpdata=fread(fid,inf,'single');fclose(fid);
                            precj=reshape(tmpdata,cols,rows,24);
                            precipi = precj(:,:,1) - preci(:,:,ihour);
                            
                            fid = fopen([outfilenextdir 'DLWR_surf'],'rb','ieee-le');
                            tmpdata=fread(fid,inf,'single');fclose(fid);
                            dlwrj=reshape(tmpdata,cols,rows,24);
                            fid = fopen([outfilenextdir 'DSWR_surf'],'rb','ieee-le');
                            tmpdata=fread(fid,inf,'single');fclose(fid);
                            dswrj=reshape(tmpdata,cols,rows,24);
                            dlwri = dlwrj(:,:,1) - dlwr(:,:,ihour);
                            dswri = dswrj(:,:,1) - dswr(:,:,ihour);
                            
                            
                        end
                        precipi = squeeze(precipi)*1000;
                        dlwri = squeeze(dlwri)/3600;
                        dswri = squeeze(dswri)/3600;
                        
                        new_wind(:,:,ihour) = wind;
                        new_p(:,:,ihour) = precipi;
                        new_sh(:,:,ihour) = SH;
                        new_dlwr(:,:,ihour) = dlwri;
                        new_dswr(:,:,ihour) = dswri;
                        
                        
                    end
                    
                    
                    fnm=[outfiledir,'WindSpeed_10m'];
                    predata=reshape(new_wind,[cols,rows*24]);
                    fid = fopen(fnm,'wb');
                    fwrite(fid,predata,'single');
                    fclose(fid);
                    
                    fid=fopen([fnm,'_sta.txt'],'w');
                    fprintf(fid,'Min    = %.4f\n',min(predata(:)));
                    fprintf(fid,'Max    = %.4f\n',max(predata(:)));
                    fprintf(fid,'Mean   = %.4f\n',mean(predata(:)));
                    fprintf(fid,'Median = %.4f\n',median(predata(:)));
                    fprintf(fid,'Std.   = %.4f\n',std(predata(:)));
                    fclose(fid);
                    
                    fnm=[outfiledir,'SpecHumi_10m'];
                    predata=reshape(new_sh,[cols,rows*24]);
                    fid = fopen(fnm,'wb');
                    fwrite(fid,predata,'single');
                    fclose(fid);
                    
                    fid=fopen([fnm,'_sta.txt'],'w');
                    fprintf(fid,'Min    = %.4f\n',min(predata(:)));
                    fprintf(fid,'Max    = %.4f\n',max(predata(:)));
                    fprintf(fid,'Mean   = %.4f\n',mean(predata(:)));
                    fprintf(fid,'Median = %.4f\n',median(predata(:)));
                    fprintf(fid,'Std.   = %.4f\n',std(predata(:)));
                    fclose(fid);
                    
                    fnm=[outfiledir,'Precipitation_ERA5_land_oriunit'];
                    predata=reshape(new_p,[cols,rows*24]);
                    fid = fopen(fnm,'wb');
                    fwrite(fid,predata,'single');
                    fclose(fid);
                    
                    fid=fopen([fnm,'_sta.txt'],'w');
                    fprintf(fid,'Min    = %.4f\n',min(predata(:)));
                    fprintf(fid,'Max    = %.4f\n',max(predata(:)));
                    fprintf(fid,'Mean   = %.4f\n',mean(predata(:)));
                    fprintf(fid,'Median = %.4f\n',median(predata(:)));
                    fprintf(fid,'Std.   = %.4f\n',std(predata(:)));
                    fclose(fid);
                    
                    fnm=[outfiledir,'DSWR_surfunit'];
                    predata=reshape(new_dswr,[cols,rows*24]);
                    fid = fopen(fnm,'wb');
                    fwrite(fid,predata,'single');
                    fclose(fid);
                    
                    fid=fopen([fnm,'_sta.txt'],'w');
                    fprintf(fid,'Min    = %.4f\n',min(predata(:)));
                    fprintf(fid,'Max    = %.4f\n',max(predata(:)));
                    fprintf(fid,'Mean   = %.4f\n',mean(predata(:)));
                    fprintf(fid,'Median = %.4f\n',median(predata(:)));
                    fprintf(fid,'Std.   = %.4f\n',std(predata(:)));
                    fclose(fid);
                    
                    fnm=[outfiledir,'DLWR_surfunit'];
                    predata=reshape(new_dlwr,[cols,rows*24]);
                    fid = fopen(fnm,'wb');
                    fwrite(fid,predata,'single');
                    fclose(fid);
                    
                    fid=fopen([fnm,'_sta.txt'],'w');
                    fprintf(fid,'Min    = %.4f\n',min(predata(:)));
                    fprintf(fid,'Max    = %.4f\n',max(predata(:)));
                    fprintf(fid,'Mean   = %.4f\n',mean(predata(:)));
                    fprintf(fid,'Median = %.4f\n',median(predata(:)));
                    fprintf(fid,'Std.   = %.4f\n',std(predata(:)));
                    fclose(fid);
                    
                end
            end
        end
        
    end
    
end




%% (key) assemble key basin information basins_info_sum 
% based on selected events and basins, to create inputs, as this step takes
% a long time, so not every year, every month is systemeatically generated,
% but rather conditional on gauges and events selection
wge = load(['/shared/dondo/home/ml423/DA/Critical matrices/world_selected_gauge_events_v3.mat'],'selected_gauge_events');
wge = wge.selected_gauge_events;
selected_gauge_events = wge;
% for each gauge find the first and the last event.
gauge_events_first_last = [];
for i = 1:size(wge,1)
    for j = 1:size(wge,2)
        if length(wge{i,j})<=1% no events qualified for this gauge
        else
            gid = str2double(wge{i,j}(1));
            evets = (char(wge{i,j}(2:end)));
            dasum = [];
            for k = 1:length(wge{i,j})-1
                datetmp = evets(:,:,k);
                datenber = str2num(datetmp(1:4))*10000+str2num(datetmp(6:7))*100+str2num(datetmp(9:10));
                dasum = [dasum;datenber];
            end
            gauge_events_first_last =[gauge_events_first_last; [gid,min(dasum),max(dasum),i]];
        end
    end
end

home_alps = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Alps/streamflow_gauges_daily/';
home_andes = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Andes/streamflow_daily_entire_Andes/';
home_himalayas = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Himalayas/streamflow_daily_entire_Himalayas/';

home_alpsdim = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Alps/input_output/';
home_andesdim = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Andes/input_output/';
home_himalayasdim = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Himalayas/input_output/';

for k = 1:length(gauge_events_first_last)
    mountain = gauge_events_first_last(k,4);
    gid = gauge_events_first_last(k,1);
    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_ind.mat'],'ind');
    ind = tmp.ind;
    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_intt.mat'],'intt');
    intt = tmp.intt;

    if mountain == 1
        home_moun = home_alps;mountname = 'Alps';home_mounin = home_alpsdim;
    end
    if mountain == 2
        home_moun = home_andes;mountname = 'Andes';home_mounin = home_andesdim;
    end
    if mountain == 3
        home_moun = home_himalayas;mountname = 'Himalayas';home_mounin = home_himalayasdim;
    end

    ffnm=[home_mounin 'fortran/Basin' num2str(gid) '/dem.bin'];
    fid=fopen(ffnm,'rb','ieee-le');
    tmpdata=fread(fid,inf,'single');
    fclose(fid);
    

    basin_info = load([home_mounin 'gauges_loc/Basins_dimention.mat'],'Basins_sum');
    binfo = basin_info.Basins_sum;
    tind = find(binfo(:,6)==gid);

    demb1 = reshape(tmpdata,binfo(tind,2),binfo(tind,1));
    demb1(ind) = nan;


    gauge_events_first_last(k,5) = binfo(tind,5);
    gauge_events_first_last(k,6) = binfo(tind,1);
    gauge_events_first_last(k,7) = binfo(tind,2);
    gauge_events_first_last(k,8) = binfo(tind,3);
    gauge_events_first_last(k,9) = binfo(tind,4);
    gauge_events_first_last(k,10) = mean(demb1(intt));%mean elevation
    gauge_events_first_last(k,11) = max(demb1(intt))-min(demb1(intt));%relief

end

%% Write Qi, obs as in fortran fixed (Qi), and single day, allrecords (obs) format
% based on selected events and basins, to create inputs, as this step takes
% a long time, so not every year, every month is systemeatically generated,
% but rather conditional on gauges and events selection
wge = load(['/shared/dondo/home/ml423/DA/Critical matrices/world_selected_gauge_events_v3.mat'],'selected_gauge_events');
wge = wge.selected_gauge_events;
selected_gauge_events = wge;
% for each gauge find the first and the last event.
% gauge_events_first_last = [];
% for i = 1:size(wge,1)
%    for j = 1:size(wge,2)
%        if length(wge{i,j})<=1% no events qualified for this gauge
%        else
%            gid = str2double(wge{i,j}(1));
%            evets = (char(wge{i,j}(2:end)));
%            dasum = [];
%            for k = 1:length(wge{i,j})-1
%                datetmp = evets(:,:,k);
%                datenber = str2num(datetmp(1:4))*10000+str2num(datetmp(6:7))*100+str2num(datetmp(9:10));
%                dasum = [dasum;datenber];
%            end
%            gauge_events_first_last =[gauge_events_first_last; [gid,min(dasum),max(dasum),i]];
%        end
%    end
% end
% or just extend the selected_gauge_events for 30 days before and after

home_alps = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Alps/streamflow_gauges_daily/';
home_andes = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Andes/streamflow_daily_entire_Andes/';
home_himalayas = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Himalayas/streamflow_daily_entire_Himalayas/';

home_alpsdim = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Alps/input_output/';
home_andesdim = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Andes/input_output/';
home_himalayasdim = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Himalayas/input_output/';


for mountain = 1
    if mountain == 1
        home_moun = home_alps;mountname = 'Alps';home_mounin = home_alpsdim;
    end
    if mountain == 2
        home_moun = home_andes;mountname = 'Andes';home_mounin = home_andesdim;
    end
    if mountain == 3
        home_moun = home_himalayas;mountname = 'Himalayas';home_mounin = home_himalayasdim;
    end
    
    for gauge = 1:length(selected_gauge_events)
        disp(mountain)
        disp(gauge)
        
        
        qc = selected_gauge_events{mountain,gauge};
        if length(qc)<2
            continue
        end
        
        gid = str2num(qc(1));

        if ismember(gid,mount_gauges_miss_alps)&&~ismember(gid,WORLD_events)
            disp('processing missing')
            disp(gid)
        else
            continue
        end
        
        odir2 = ['/shared/dondo/home/ml423/world_mts/obs_events/'];
        %if isfile([odir2 'obs_' num2str(gid,'%7.0d') '.mat'])
            
           % continue
        if 1==0
        else
            
            
            disp(gid)
            disp(mountain)
            
            outdir = ['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_fixed_900m/Qi/'];
            mkdir(outdir)
            
            
            
            basin_info = load([home_mounin 'gauges_loc/Basins_dimention.mat'],'Basins_sum');
            binfo = basin_info.Basins_sum;
            tind = find(binfo(:,6)==gid);
            
            charac = binfo(tind,[1 2 6 3 4]);
            charac(3) = charac(1)*charac(2);
            fileID = fopen(['/shared/dondo/home/ml423/world_mts/Basins_events_chars/Basin' num2str(gid,'%7.0d') '.txt'],'w');
            fprintf(fileID,'%d\n',charac);
            fclose(fileID);
            
            % write obs for plotting purpose and IRC purpose
            %
            lastall = ls([home_moun num2str(gid) '*.txt']);
            c2 = strsplit(lastall);
            ii = 1;
            
            warning('off')
            close all
            ffnm = c2{ii};
            discharge=readtable(ffnm);
            if size(discharge,2)<3
                % no records at all
                continue
            end
            
            discharge_v = table2array(discharge(:,3));
            obsnowall=discharge_v;
            if length(discharge_v)<30
                continue
            end
            disp(ii)
            
            
            fid=fopen(ffnm,'r');
            n_lines = 10;
            outlet_info = cell(10,1);
            for jj = 1:15
                outlet_info(jj) = {fgetl(fid)};
            end
            fclose(fid);
            
            gauge_no = outlet_info{9}(26:end);
            g_id = str2double(gauge_no);
            
            time_s = table2cell(discharge(:,1));
            record_len = length(time_s);
            
            
            odir = ['/shared/dondo/home/ml423/world_mts/obs_events/days/obs_' num2str(g_id,'%7.0d') '/'];
            mkdir(odir)
            odir2 = ['/shared/dondo/home/ml423/world_mts/obs_events/'];
            mkdir(odir2)
            obs = [];
            
            for il = 1:record_len
                [ys,ms,ds] = ymd(time_s{il});
                datstr = ys*10000+ms*100+ds;
                obsnow = obsnowall(il);
                save([odir num2str(datstr) '.mat'],'obsnow')
                obs(il,1) = datstr;
                obs(il,2) = obsnow;
            end
            save([odir2 'obs_' num2str(g_id,'%7.0d') '.mat'],'obs')
            
            %}



        
        lastall = ls([home_moun num2str(gid) '*.txt']);
        c2 = strsplit(lastall);
        
        
        warning('off')
        close all
        ffnm = c2{1};
        discharge=readtable(ffnm);
        if size(discharge,2)<3
            % no records at all
            continue
        end
        
        discharge_v = table2array(discharge(:,3));
        obsnow=discharge_v;
        if length(discharge_v)<30
            continue
        end
       % clear days
        discharge_t = table2array(discharge(:,1));

        parfor event = 2:length(qc)
            
            datm = datetime(qc(event));
            
            disp(mountname)
            disp(gid)
            disp(datm)
            
            datmi = datm-days(20);
            datme = datm+days(20);
            
            t1loc = find(discharge_t==datmi);
            enddateloc = min([length(discharge_t),t1loc+41]);
            % should end either when there is no discharge obs, or 41 day
            % boundary
            for tloc = t1loc:enddateloc
                Qi = discharge_v(tloc);
                [yo,mo,do] = ymd(discharge_t(tloc));
                daystr = [num2str(yo) num2str(mo,'%2.2d') num2str(do,'%2.2d')];
                % write iniQi in the fixed directory
                fid = fopen([outdir 'iniQi_' daystr '.ascii'],'w');
                fprintf(fid,'%f',Qi);
                fclose(fid);
                
            end


        end
            %}
        end
    end
end


        

%% plot a few inputs for graphs


ffnm=['/shared/dondo/home/ml423/ERA5/ERA5_data/data2/2021_01_wind_u.grib'];
gf = ncgeodataset(ffnm);

mountain = 1;

if mountain == 1
    home_moun = home_alps;mtname = 'Alps';
end
if mountain == 2
    home_moun = home_andes;mtname = 'Andes';
end
if mountain == 3
    home_moun = home_himalayas;mtname = 'Himalayas';
end
tmp = load([home_moun 'gauges_loc/Basins_dimention.mat'],'Basins_sum');
Basins_sum_region = tmp.Basins_sum;

bid = 1;
cols = Basins_sum_region(bid,2);
rows = Basins_sum_region(bid,1);
gid = Basins_sum_region(bid,6);
outhomedir = ['/shared/dondo/home/ml423/world_mts/'];
ffnm = ([outhomedir 'Basin' num2str(gid) '_ERA5_900m1hr/19880101/AirTemp_10m']);
fid=fopen(ffnm,'rb','ieee-le');
tmpdata=fread(fid,inf,'single');
airtemp1=reshape(tmpdata,cols,rows,24);
air = airtemp1(:,:,1);
tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_ind.mat'],'ind');
        ind = tmp.ind;
        tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_intt.mat'],'intt');
        intt = tmp.intt;
        air(ind) = nan;
figure
imagesc(air')

ffnm=[home_moun 'fortran/Basin' num2str(gid) '/dem.bin'];
        fid=fopen(ffnm,'rb','ieee-le');
        tmpdata=fread(fid,inf,'single');
        fclose(fid);
        demb1 = reshape(tmpdata,cols,rows);
        demb1(ind) = nan;
        figure
        imagesc(demb1')
gid = 6125361;        
ffnm = ([outhomedir 'Basin' num2str(gid) '_ERA5_900m1hr/19780402/Precipitation_ERA5_land_oriunit']);
fid=fopen(ffnm,'rb','ieee-le');
tmpdata=fread(fid,inf,'single');
airtemp1=reshape(tmpdata,cols,rows,24);
air = airtemp1(:,:,1);
%% Plot streamflow observations as an example
wge = load(['/shared/dondo/home/ml423/DA/Critical matrices/world_selected_gauge_events_allinfo_v3.mat'],'gauge_events_first_last');
ge_info = wge.gauge_events_first_last;clear wge;

tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/WORLD_events_v3.mat'],'WORLD_events');
WORLD_events = tmp.WORLD_events;clear tmp;


for g_shell_index = [299]
    close all
    gid = WORLD_events(g_shell_index,1);%6559101;%6935540;%;%6235535;%6559100;%6125361;

    tmp = load(['/shared/dondo/home/ml423/world_mts/obs_events/obs_' num2str(gid,'%7.0d') '.mat'],'obs');
    obs = tmp.obs;
    figure
    set(gcf,'outerposition', [10 100 1500 350]);
    plot(obs(:,2))
    ylim([-10 1.2*max(obs(:,2))])
    title([num2str(g_shell_index)])
    saveas(gcf,[picdir, 'simple_obs' num2str(g_shell_index) '_Basin' num2str(gid,'%7.0d') ],'png')

end
%% Plot spinup period streamflow from DCHM

wge = load(['/shared/dondo/home/ml423/DA/Critical matrices/world_selected_gauge_events_allinfo_v3.mat'],'gauge_events_first_last');
ge_info = wge.gauge_events_first_last;clear wge;
%flow_vol_r = [];
plotfig = 1;

tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/WORLD_events_v3.mat'],'WORLD_events');
WORLD_events = tmp.WORLD_events;clear tmp;


for g_shell_index = [250]%147%1:450%[29 31 34 45 46 54 60 61 66 67 68 69 70]%[297:566]%99%[234:566]%[1:22 24:105 107:160]% problem 23 106
    %close all
    preset_gauge = WORLD_events(g_shell_index,1);%6559101;%6935540;%;%6235535;%6559100;%6125361;
    preset_gauge
    suffix = 'S2';
    intdays = 30;
    spinitr = [3];% spinup iteration to be plotted
    smooth_simu = 0;%[-1 is original simu resolution hourly, 0 is pick 00:00, 1 is average of 00:00-24:00]
    
    galoc = find(ge_info(:,1)==preset_gauge);mtindex = ge_info(galoc,4);
    stdate = ge_info(galoc,2);etdate = ge_info(galoc,3);
    if stdate<19730911
        stdate = 19730911;
    end
    fyear = floor(stdate/10000);fmonth = floor((stdate-fyear*10000)/100);fday=mod(stdate,100);
    startdate = datetime(fyear,fmonth,fday)-days(9);[siy,sim,sid] = ymd(startdate);
    lyear = floor(etdate/10000);lmonth = floor((etdate-lyear*10000)/100);lday=mod(etdate,100);
    enddate = datetime(lyear,lmonth,lday)+days(9);[siy2,sim2,sid2] = ymd(enddate);
    
    spin_info = [siy,sim,sid;...
        siy+1,4,30];   % iteration period
    spin_info = [siy+2,sim,sid;...
                 siy2,sim2,sid2];   % full spinup results
    %spin_info = [1996,7,20;1996,9,19];
    
    
    home_alpsdim = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Alps/input_output/';
    home_andesdim = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Andes/input_output/';
    home_himalayasdim = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Himalayas/input_output/';
    
    
    flowtype = {'streamflowarea','interflowarea','overlandflowarea','baseflowarea','groundwaterflow'};
    ftp = 1;yenforce=0;ymmax=300;
    
    totdays = datetime(spin_info(1,1),spin_info(1,2),spin_info(1,3)):days(1):datetime(spin_info(2,1),spin_info(2,2),spin_info(2,3));
    
    plotdays_itv = ceil(length(totdays)/9)-2;% then we have about 10 ticks per graph
    
    mountain = mtindex;complete_simu = 1;
    
    if mountain == 1
        mountname = 'Alps';home_mounin = home_alpsdim;
    end
    if mountain == 2
        mountname = 'Andes';home_mounin = home_andesdim;
    end
    if mountain == 3
        mountname = 'Himalayas';home_mounin = home_himalayasdim;
    end
    
    basin_info = load([home_mounin 'gauges_loc/Basins_dimention.mat'],'Basins_sum');
    binfo = basin_info.Basins_sum;
    gloc = find(binfo(:,6)==preset_gauge);
    gauge = gloc;
    
    
    disp(mountain)
    disp(binfo(gauge,6))
    
    r = binfo(gauge,1);
    c = binfo(gauge,2);
    rg = binfo(gauge,3);
    cg = binfo(gauge,4);
    gid = binfo(gauge,6);
    area = binfo(gauge,5);
    
    
    
    
    spin_start_year = spin_info(1,1);spin_end_year = spin_info(2,1);
    spin_start_month = spin_info(1,2);spin_end_month = spin_info(2,2);
    spin_start_day = spin_info(1,3);spin_end_day = spin_info(2,3);
    nm = spin_end_year*10000+spin_end_month*100+spin_end_day;
    nsm = spin_start_year*10000+spin_start_month*100+spin_start_day;
    
    
    
    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid,'%7.0d') '_intt.mat'],'intt');
    intt = tmp.intt;
    
    
    
    datetime.setDefaultFormats('default','yyyy-MM-dd hh:mm:ss')
    
    home_dir2 = ['/shared/dondo/home/ml423/world_mts_runs/figs/'];
    picdir=[home_dir2,'Spinup_images/'];mkdir(picdir);
    clear days
    spinstart = datetime(spin_start_year,spin_start_month,spin_start_day);
    spinend = datetime(spin_end_year,spin_end_month,spin_end_day);
    spindays = spinstart:days(1):spinend;
    [ys,ms,ds] = ymd(spindays);
    
    
    
    startstep = 1;
    showstep = length(spindays)*24;
    xt=1:plotdays_itv*24:showstep;k=1;
    for i=xt
        label = ms(max(floor(i/24),1))*100+ds(max(floor(i/24),1));label=num2str(label,'%4.4d');% get the date for each tick
        xs{k}=label;
        k=k+1;
    end
    xstr=['Time ' num2str(spin_start_year) '-' num2str(spin_end_year) '/MMDD'];
    
    if plotfig == 1
        % calculate heytograph
        rainall = [];
        ii = 0;
        for i = 1:intdays:length(spindays)
            ii = ii+1;
            daystr = strcat(num2str(ys(i)),num2str(ms(i),'%2.2d'),num2str(ds(i),'%2.2d'));
            ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid,'%7.0d') '_ERA5_900m1hr/' daystr '/Precipitation_ERA5_land_oriunit'];
            fid=fopen(ffnm,'rb','ieee-le');
            tmpdata=fread(fid,inf,'single');
            tmp=reshape(tmpdata,c,r,24);
            fclose(fid);
            for j  = 1:24
                tmp1 = tmp(:,:,j);
                rainall((ii-1)*24+j) = mean(tmp1(intt));
            end
        end
        % calculate temp
        tempall = [];
        ii = 0;
        for i = 1:intdays:length(spindays)
            ii = ii+1;
            daystr = strcat(num2str(ys(i)),num2str(ms(i),'%2.2d'),num2str(ds(i),'%2.2d'));
            ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid,'%7.0d') '_ERA5_900m1hr/' daystr '/AirTemp_10m'];
            fid=fopen(ffnm,'rb','ieee-le');
            tmpdata=fread(fid,inf,'single');
            tmp=reshape(tmpdata,c,r,24);
            fclose(fid);
            for j  = 1:24
                tmp1 = tmp(:,:,j);
                tempall((ii-1)*24+j) = mean(tmp1(intt));
            end
        end
        rainallcold = rainall;
        rainallwarm = rainall;
        rainallcold(find(tempall>=273.15)) = nan;
        rainallwarm(find(tempall<273.15)) = nan;
    end
    
    % plotting flows
    for spinup = spinitr

        s_flow = [];
        obs_flow = [];
        for i = 1:intdays:length(spindays)
            disp(i)
            daystr = strcat(num2str(ys(i)),num2str(ms(i),'%2.2d'),num2str(ds(i),'%2.2d'));
            ffnm=['/shared/dondo/home/ml423/world_mts_outputs/Basin' num2str(gid,'%7.0d') '_output_900m1hr_' suffix '_spinup/' daystr '/' flowtype{ftp} '_' num2str(spinup)];
            fid=fopen(ffnm,'rb','ieee-le');
            if fid == -1
                complete_simu = 0;
                break
            else
                complete_simu = 1;
            end
            
            tmpdata=fread(fid,inf,'single');
            fclose(fid);
            tmp=reshape(tmpdata,c,r,24);
            flow = [];
            for k = 1:24
                flow(k) = tmp(cg,rg,k);
            
            end
            if max(flow>6000)%||min(flow<=0.01)
                disp(daystr)
                break
            end
            s_flow = [s_flow;flow'];
            
            % obs
            tmp = load(['/shared/dondo/home/ml423/world_mts/obs_events/days/obs_' num2str(gid,'%7.0d') '/' daystr '.mat'],'obsnow');
            obs = tmp.obsnow;
            obs_flow = [obs_flow;obs];
        end
            
        
        if complete_simu == 0
            break
        end
        
        obstep = length(obs_flow);
        
        simu = s_flow(1:24:obstep*24);
        smoo_flow = [];
        for smoo = 1:floor(length(s_flow)/24)
            smoo_flow(smoo,1) = mean(s_flow(((smoo-1)*24+1):(smoo*24)));
        end
        tmp0 = 0;
        tmp1 = 0;
        
   
        n_loc = find(obs_flow>=0);
        
        obsnowx = obs_flow(n_loc);
        simuna = simu(n_loc);
        tmp2 = mean(obsnowx);
        for i = 1:length(obsnowx)
            tmp0 = tmp0 + (simuna(i)-obsnowx(i))^2;
            tmp1 = tmp1 + (tmp2-obsnowx(i))^2;
        end
        NSE = 1-tmp0./tmp1;
        
        rrr = corrcoef(simuna',obsnowx);
        rstar = rrr(1,2);
        obsstd = std(obsnowx);
        simustd = std(simuna);
        obsmean = mean(obsnowx);
        simumean = mean(simuna);
        KGE = 1-sqrt((rstar-1)^2+(simustd/obsstd-1)^2+(simumean/obsmean-1)^2);
        %flow_vol_r = [flow_vol_r;[gid,mountain,area,sum(simuna)/sum(obsnowx),KGE]];
        
        if plotfig == 1
            
            ymx=ceil(max(max(obs_flow),max(simu))/30)*30;pmax=ceil(1.5*max(rainall)/15)*15;
            if yenforce==1
                ymx = ymmax;
            end
            
            figure('Color',[1 1 1]);
            set(gcf,'outerposition', [10 100 1500 550]);
            
            %-----------------------------------------
            
            subplot('Position',[0.095 0.18 0.82 0.78]);
            
            [AX,H1,H2] = plotyy(1:24*intdays:obstep*24*intdays,obs_flow,1:intdays:length(rainall)*intdays,rainallwarm,'line','stem'); hold on;
            if smooth_simu == 0
                H4=plot(1:24*intdays:size(simu,1)*24*intdays,simu(:));
            elseif smooth_simu == -1
                H4=plot(1:size(s_flow,1),s_flow(:));
            else
                H4=plot(1:24:size(smoo_flow,1)*24,smoo_flow(:));
            end
            
            hold(AX(2));
            H3 = stem(AX(2),1:intdays:length(rainall)*intdays,rainallcold);
            set(H3,'Marker','none')
            set(H1,'LineStyle','-','Linewidth',5,'Marker','d','MarkerSize',0.5,'Color','k');
            set(H2,'LineStyle','-','Linewidth',1.5,'Marker','.','MarkerSize',0.001,'Color',[0,0,0.7]);
            set(H3,'LineStyle','-','Linewidth',2.5,'Marker','.','MarkerSize',0.001,'Color',[1,0,0]);
            set(H4,'LineStyle','-','Linewidth',2.75,'Marker','d','MarkerSize',0.5,'Color','g');
            axis(AX(1),[startstep showstep 0 ymx]);axis(AX(2),[startstep showstep 0 pmax]);
            
            set(AX(1),'XTick',xt,'XTickLabel',xs,'Fontsize',24,'FontWeight','bold');
            set(AX(2),'XTick',xt,'XTickLabel',xs,'Fontsize',24,'FontWeight','bold');
            set(AX(1),'YTick',[0:ymx/5:ymx],'YTickLabel',[0:ymx/5:ymx],'YColor','k','Fontsize',24);
            set(AX(2),'YDir','reverse','YTick',[0:pmax/5:pmax],'YColor','k','Fontsize',24);
            set(get(AX(1),'Ylabel'),'String','Streamflow(m^3/s)','Fontsize',24,'FontWeight','bold');
            set(get(AX(2),'Ylabel'),'String','Precip(mm/hr)','Fontsize',24,'FontWeight','bold');
            grid on;
            
            P=legend([H2 H1 H4],'Precip.','Obs.',['Sim. ' num2str(spinup) ' Basin' num2str(gid,'%7.0d') ' Area ' num2str(area,'%2.0f') 'km^2 ' suffix]);
            set(P,'Color','None','Location','NorthOutside','Orientation','horizontal','Fontsize',22,'linewidth',2);
            legend boxoff
            xlabel(xstr,'Fontsize',24,'FontWeight','bold');
            
            titstr=['NSE = ' num2str(NSE,'%2.2f')];
            xminmax=get(gca,'xlim');xmin=xminmax(1);xmax=xminmax(2);
            %text(xmin+(xmax-xmin)/100,ymx/2,titstr,'Fontsize',24,'FontWeight','bold');
            set(gca,'linewidth',3);
            
            set(gcf,'PaperPositionMode','auto');
            hold off
            hold(AX(2),'off');
            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0, 0, 14.3, 5.3];
            %saveas(gcf,[picdir, 'HRRR_F01_Basin' num2str(bid,'%2.2d') '_spinup-' flowtype{ftp} '-' num2str(spinup) '-' suffix],'png')
            
            saveas(gcf,[picdir, 'tmpS2spin' num2str(g_shell_index) '_Basin' num2str(gid,'%7.0d') '_' num2str(nsm) '_' num2str(nm) '_spinup-' flowtype{ftp} '-' num2str(spinup) '-' suffix 'wTemp'],'png')
        end
    end
    
end
%% plot spinup when determining coeflateral coefficients

wge = load(['/shared/dondo/home/ml423/DA/Critical matrices/world_selected_gauge_events_allinfo_v3.mat'],'gauge_events_first_last');
ge_info = wge.gauge_events_first_last;clear wge;
%flow_vol_r = [];
plotfig = 1;

tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/WORLD_events_v3.mat'],'WORLD_events');
WORLD_events = tmp.WORLD_events;clear tmp;

%ratio_changes = [];
for cali = [1]% calibration of coeflateral in the loop in the shell script, 1 is original 1000 300 30 3, and 2 is after one multiplication
    %close all
    
    disp(cali)
    for g_shell_index = [1:522]%[1:60 101:487]%99%[234:566]%[1:22 24:105 107:160]% problem 23 106
        close all
        preset_gauge = WORLD_events(g_shell_index,1);
       
        suffix = 'S9';
        spinitr = [1];% spinup iteration to be plotted
        smooth_simu = 0;%[-1 is original simu resolution hourly, 0 is pick 00:00, 1 is average of 00:00-24:00]
        tmp = load(['/shared/dondo/home/ml423/world_mts/obs_events/obs_' num2str(preset_gauge) '.mat'],'obs');
        obsgid = tmp.obs; clear tmp

        galoc = find(ge_info(:,1)==preset_gauge);mtindex = ge_info(galoc,4);
        stdate = ge_info(galoc,2);etdate = ge_info(galoc,3);

        if stdate<19730911
            final1 = 19730911;
        else
            final1 = stdate;
        end
        
        obsloc = find(obsgid(:,1)==final1);
        
        for obsloop = obsloc:length(obsgid)
            if obsgid(obsloop,2)>0
                final1 = obsgid(obsloop,1);
                break
            end
        end
        
        fyear = floor(final1/10000);fmonth = floor((final1-fyear*10000)/100);fday=mod(final1,100);
        startdate = datetime(fyear,fmonth,fday)-days(9);[siy,sim,sid] = ymd((startdate));
        lyear = floor(etdate/10000);lmonth = floor((etdate-lyear*10000)/100);lday=mod(etdate,100);
        enddate = datetime(lyear,lmonth,lday)+days(9);[siy2,sim2,sid2] = ymd(enddate);
        
        spin_info = [siy,sim,sid;...
            siy+1,4,30];   % iteration period
        %spin_info = [siy,sim,sid;...
        %             siy2,sim2,sid2];   % full spinup results
        
        
        home_alpsdim = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Alps/input_output/';
        home_andesdim = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Andes/input_output/';
        home_himalayasdim = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Himalayas/input_output/';
        
        
        flowtype = {'streamflowarea','interflowarea','overlandflowarea','baseflowarea','groundwaterflow'};
        ftp = 1;yenforce=0;ymmax=300;
        
        totdays = datetime(spin_info(1,1),spin_info(1,2),spin_info(1,3)):days(1):datetime(spin_info(2,1),spin_info(2,2),spin_info(2,3));
        
        plotdays_itv = ceil(length(totdays)/9)-2;% then we have about 10 ticks per graph
        
        mountain = mtindex;complete_simu = 1;
        
        if mountain == 1
            mountname = 'Alps';home_mounin = home_alpsdim;
        end
        if mountain == 2
            mountname = 'Andes';home_mounin = home_andesdim;
        end
        if mountain == 3
            mountname = 'Himalayas';home_mounin = home_himalayasdim;
        end
        
        basin_info = load([home_mounin 'gauges_loc/Basins_dimention.mat'],'Basins_sum');
        binfo = basin_info.Basins_sum;
        gloc = find(binfo(:,6)==preset_gauge);
        gauge = gloc;
        
        
        disp(mountain)
        disp(binfo(gauge,6))
        
        r = binfo(gauge,1);
        c = binfo(gauge,2);
        rg = binfo(gauge,3);
        cg = binfo(gauge,4);
        gid = binfo(gauge,6);
        area = binfo(gauge,5);
        
        
        
        
        spin_start_year = spin_info(1,1);spin_end_year = spin_info(2,1);
        spin_start_month = spin_info(1,2);spin_end_month = spin_info(2,2);
        spin_start_day = spin_info(1,3);spin_end_day = spin_info(2,3);
        nm = spin_end_year*10000+spin_end_month*100+spin_end_day;
        nsm = spin_start_year*10000+spin_start_month*100+spin_start_day;
        
        
        
        tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid,'%7.0d') '_intt.mat'],'intt');
        intt = tmp.intt;
        
        
        
        datetime.setDefaultFormats('default','yyyy-MM-dd hh:mm:ss')
        
        home_dir2 = ['/shared/dondo/home/ml423/world_mts_runs/figs/'];
        picdir=[home_dir2,'Spinup_images/'];mkdir(picdir);
        clear days
        spinstart = datetime(spin_start_year,spin_start_month,spin_start_day);
        spinend = datetime(spin_end_year,spin_end_month,spin_end_day);
        spindays = spinstart:days(1):spinend;
        [ys,ms,ds] = ymd(spindays);
        
        
        
        startstep = 1;
        showstep = length(spindays)*24;
        xt=1:plotdays_itv*24:showstep;k=1;
        for i=xt
            label = ms(max(floor(i/24),1))*100+ds(max(floor(i/24),1));label=num2str(label,'%4.4d');% get the date for each tick
            xs{k}=label;
            k=k+1;
        end
        xstr=['Time ' num2str(spin_start_year) '-' num2str(spin_end_year) '/MMDD'];
        
        if plotfig == 1
            % calculate heytograph
            rainall = [];
            for i = 1:length(spindays)
                daystr = strcat(num2str(ys(i)),num2str(ms(i),'%2.2d'),num2str(ds(i),'%2.2d'));
                ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid,'%7.0d') '_ERA5_900m1hr/' daystr '/Precipitation_ERA5_land_oriunit'];
                fid=fopen(ffnm,'rb','ieee-le');
                tmpdata=fread(fid,inf,'single');
                tmp=reshape(tmpdata,c,r,24);
                fclose(fid);
                for j  = 1:24
                    tmp1 = tmp(:,:,j);
                    rainall((i-1)*24+j) = mean(tmp1(intt));
                end
            end
            % calculate temp
            tempall = [];
            for i = 1:length(spindays)
                daystr = strcat(num2str(ys(i)),num2str(ms(i),'%2.2d'),num2str(ds(i),'%2.2d'));
                ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid,'%7.0d') '_ERA5_900m1hr/' daystr '/AirTemp_10m'];
                fid=fopen(ffnm,'rb','ieee-le');
                tmpdata=fread(fid,inf,'single');
                tmp=reshape(tmpdata,c,r,24);
                fclose(fid);
                for j  = 1:24
                    tmp1 = tmp(:,:,j);
                    tempall((i-1)*24+j) = mean(tmp1(intt));
                end
            end
            rainallcold = rainall;
            rainallwarm = rainall;
            rainallcold(find(tempall>=273.15)) = nan;
            rainallwarm(find(tempall<273.15)) = nan;
        end
        
        % plotting flows
        for spinup = spinitr
            s_flow = [];
            obs_flow = [];
            for i = 1:length(spindays)
                daystr = strcat(num2str(ys(i)),num2str(ms(i),'%2.2d'),num2str(ds(i),'%2.2d'));
                ffnm=['/shared/dondo/home/ml423/world_mts_outputs/Basin' num2str(gid,'%7.0d') '_output_900m1hr_' suffix '_spinupsum/Basin' num2str(gid,'%7.0d') '_output_900m1hr_' suffix '_spinup' num2str(cali) '/' daystr '/' flowtype{ftp} '_' num2str(spinup)];
                fid=fopen(ffnm,'rb','ieee-le');
                if fid == -1
                    complete_simu = 0;
                    complete_simu
                    disp(ffnm)
                    break
                else
                    complete_simu = 1;
                end
                
                tmpdata=fread(fid,inf,'single');
                fclose(fid);
                tmp=reshape(tmpdata,c,r,24);
                flow = [];
                for k = 1:24
                    flow(k) = tmp(cg,rg,k);
                end
                s_flow = [s_flow;flow'];
                
                % obs
                tmp = load(['/shared/dondo/home/ml423/world_mts/obs_events/days/obs_' num2str(gid,'%7.0d') '/' daystr '.mat'],'obsnow');
                obs = tmp.obsnow;
                obs_flow = [obs_flow;obs];
            end
            
            if complete_simu == 0
                break
            end
            
            obstep = length(obs_flow);
            
            simu = s_flow(1:24:obstep*24);
            smoo_flow = [];
            for smoo = 1:floor(length(s_flow)/24)
                smoo_flow(smoo,1) = mean(s_flow(((smoo-1)*24+1):(smoo*24)));
            end
            tmp0 = 0;
            tmp1 = 0;
            
            %n_ind = ~isnan(obs_flow);
            n_loc = find(obs_flow>0);
            
            obsnowx = obs_flow(n_loc);
            simuna = simu(n_loc);
            tmp2 = mean(obsnowx);
            for i = 1:length(obsnowx)
                tmp0 = tmp0 + (simuna(i)-obsnowx(i))^2;
                tmp1 = tmp1 + (tmp2-obsnowx(i))^2;
            end
            NSE = 1-tmp0./tmp1;
            
            [obsmax,obsmloc] = max(obsnowx);
            lbd = obsmloc-5;
            rbd = obsmloc+5;
            if lbd<1
                lbd=1;
            end
            if rbd>length(simuna)
                rbd=length(simuna);
            end
            
            flowpeak_ratio = (max(simuna(lbd:rbd))-obsmax)/obsmax;
            nanmean(obsnowx(ceil(0.33*length(simuna)):end))
            nanmean(simuna(ceil(0.33*length(simuna)):end))
            ratio_changes(g_shell_index,cali) = flowpeak_ratio;
            
            rrr = corrcoef(simuna',obsnowx);
            rstar = rrr(1,2);
            obsstd = std(obsnowx);
            simustd = std(simuna);
            obsmean = mean(obsnowx);
            simumean = mean(simuna);
            KGE = 1-sqrt((rstar-1)^2+(simustd/obsstd-1)^2+(simumean/obsmean-1)^2);
            flow_vol_r(g_shell_index,:) = [gid,mountain,area,sum(simuna)/sum(obsnowx),KGE];
            
            if plotfig == 1
                
                ymx=ceil(max(max(obs_flow),max(simu))/30)*30;pmax=ceil(1.5*max(rainall)/15)*15;
                if yenforce==1
                    ymx = ymmax;
                end
                
                figure('Color',[1 1 1]);
                set(gcf,'outerposition', [10 100 1500 550]);
                
                %-----------------------------------------
                
                subplot('Position',[0.095 0.18 0.82 0.78]);
                
                [AX,H1,H2] = plotyy(1:24:obstep*24,obs_flow,1:length(rainall),rainallwarm,'line','stem'); hold on;
                if smooth_simu == 0
                    H4=plot(1:24:size(simu,1)*24,simu(:));
                elseif smooth_simu == -1
                    H4=plot(1:size(s_flow,1),s_flow(:));
                else
                    H4=plot(1:24:size(smoo_flow,1)*24,smoo_flow(:));
                end
                
                hold(AX(2));
                H3 = stem(AX(2),1:length(rainall),rainallcold);
                set(H3,'Marker','none')
                set(H1,'LineStyle','-','Linewidth',5,'Marker','d','MarkerSize',0.5,'Color','k');
                set(H2,'LineStyle','-','Linewidth',1.5,'Marker','.','MarkerSize',0.001,'Color',[0,0,0.7]);
                set(H3,'LineStyle','-','Linewidth',2.5,'Marker','.','MarkerSize',0.001,'Color',[1,0,0]);
                set(H4,'LineStyle','-','Linewidth',2.75,'Marker','d','MarkerSize',0.5,'Color','g');
                axis(AX(1),[startstep showstep 0 ymx]);axis(AX(2),[startstep showstep 0 pmax]);
                
                set(AX(1),'XTick',xt,'XTickLabel',xs,'Fontsize',24,'FontWeight','bold');
                set(AX(2),'XTick',xt,'XTickLabel',xs,'Fontsize',24,'FontWeight','bold');
                set(AX(1),'YTick',[0:ymx/5:ymx],'YTickLabel',[0:ymx/5:ymx],'YColor','k','Fontsize',24);
                set(AX(2),'YDir','reverse','YTick',[0:pmax/5:pmax],'YColor','k','Fontsize',24);
                set(get(AX(1),'Ylabel'),'String','Streamflow(m^3/s)','Fontsize',24,'FontWeight','bold');
                set(get(AX(2),'Ylabel'),'String','Precip(mm/hr)','Fontsize',24,'FontWeight','bold');
                grid on;
                
                P=legend([H2 H1 H4],'Precip.','Obs.',['Sim. ' num2str(spinup) ' Basin' num2str(gid,'%7.0d') ' Area ' num2str(area,'%2.0f') 'km^2 ' suffix ' ' num2str(flowpeak_ratio,'%0.2f')]);
                set(P,'Color','None','Location','NorthOutside','Orientation','horizontal','Fontsize',22,'linewidth',2);
                legend boxoff
                xlabel(xstr,'Fontsize',24,'FontWeight','bold');
                
                titstr=['NSE = ' num2str(NSE,'%2.2f')];
                xminmax=get(gca,'xlim');xmin=xminmax(1);xmax=xminmax(2);
                %text(xmin+(xmax-xmin)/100,ymx/2,titstr,'Fontsize',24,'FontWeight','bold');
                set(gca,'linewidth',3);
                
                set(gcf,'PaperPositionMode','auto');
                hold off
                hold(AX(2),'off');
                fig = gcf;
                fig.PaperUnits = 'inches';
                fig.PaperPosition = [0, 0, 14.3, 5.3];
                %saveas(gcf,[picdir, 'HRRR_F01_Basin' num2str(bid,'%2.2d') '_spinup-' flowtype{ftp} '-' num2str(spinup) '-' suffix],'png')
                
                saveas(gcf,[picdir, 'Basin' num2str(gid,'%7.0d') '_' num2str(nsm) '_' num2str(nm) '_spinup-' flowtype{ftp} '-' num2str(spinup) '-' suffix '_cali' num2str(cali) 'wTemp'],'png')
            end
        end
        
    end
    
end


%% plot groundwater  
% to make sure groundwater is not below the datum
wge = load(['/shared/dondo/home/ml423/DA/Critical matrices/world_selected_gauge_events_allinfo.mat'],'gauge_events_first_last');
ge_info = wge.gauge_events_first_last;clear wge;
gauge1 = 6246155;%6342935;%6139240;
gauge2 = 6342935;
g1loc = find(ge_info(:,1)==gauge1);
r = ge_info(g1loc,6);
c = ge_info(g1loc,7);


tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gauge1) '_ind.mat'],'ind');
ind = tmp.ind;
tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gauge1) '_intt.mat'],'intt');
intt = tmp.intt;

fixedstr = ['/shared/dondo/home/ml423/world_mts/Basin' num2str(gauge1) '_fixed_900m/'];
ffnm=[fixedstr 'soildepth.ascii'];
soi = dlmread(ffnm);
dsm = reshape(soi(:,1),c,r);
dp = reshape(soi(:,2),c,r);
ds0 = 0.1*ones(c,r);
ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gauge1) '_fixed_900m/soilparameters.ascii'];
soi = dlmread(ffnm);
khh1 = reshape(soi(:,1),c,r);
figure
imagesc(khh1')
wcs1 = reshape(soi(:,5),c,r);wcs1 = wcs1*0.1;
wcs2 = reshape(soi(:,6),c,r);wcs2 = wcs2.*dsm;
wcs3 = reshape(soi(:,7),c,r);wcs3 = wcs3.*dp;
wcs4 = reshape(soi(:,8),c,r);

wfc1 = reshape(soi(:,9),c,r);wfc1 = wfc1*0.1;
wfc2 = reshape(soi(:,10),c,r);wfc2 = wfc2.*dsm;
wfc3 = reshape(soi(:,11),c,r);wfc3 = wfc3.*dp;
wfc4 = reshape(soi(:,11),c,r);
wwp1 = reshape(soi(:,13),c,r);
wwp2 = reshape(soi(:,14),c,r);
wwp3 = reshape(soi(:,15),c,r);
wwp4 = reshape(soi(:,16),c,r);

figure
imagesc(wfc1')
figure
imagesc(wwp1')

ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gauge1) '_fixed_900m/initialwatertable_basin.ascii'];
soi = dlmread(ffnm);
gd = reshape(soi(:,1),c,r);
datum = reshape(soi(:,2),c,r);

ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gauge1) '_fixed_900m/dem.bin'];
fid=fopen(ffnm,'rb','ieee-le');
tmpdata=fread(fid,inf,'single');
fclose(fid);
dem = reshape(tmpdata,c,r);
dem(ind) = nan;
demf = dem';
figure
imagesc(dem')


aquifer_thick = gd-datum;
unsaturated_thick = dem-gd;
aquifer_thick(ind) = nan;
unsaturated_thick(ind) = nan;
figure
imagesc(aquifer_thick')
title('aquifer thick')
colorbar
figure
imagesc(unsaturated_thick')
title('unsaturated thick')
colorbar
min(aquifer_thick(intt))
min(unsaturated_thick(intt))

%% plot each input and output of the DCHM to show examples 
gauge1 = 6343570;
g1loc = find(ge_info(:,1)==gauge1);
r = ge_info(g1loc,6);
c = ge_info(g1loc,7);
ffnm=['/shared/dondo/home/ml423/world_mts_outputs/Basin' num2str(gauge1) '_output_900m1hr_S1_spinup/19730906/streamflowarea_1'];
fid=fopen(ffnm,'rb','ieee-le');
tmpdata=fread(fid,inf,'single');
tmp=reshape(tmpdata,c,r,24);
fclose(fid);
for j = 1:9
tmpplot = tmp(:,:,j);
figure
imagesc(tmpplot')
title(['flow ' num2str(j) ])
colorbar
end

%


ffnm=['/shared/dondo/home/ml423/world_mts_outputs/Basin' num2str(gauge1) '_output_900m1hr_S1_spinup/19730907/laststep_Qu1'];
fid=fopen(ffnm,'rb','ieee-le');
tmpdata=fread(fid,inf,'single');
tmp=reshape(tmpdata,c,r);
fclose(fid);
figure
imagesc(tmp')
colorbar
title('0907 Qu1')
caxis([0 10])
ffnm=['/shared/dondo/home/ml423/world_mts_outputs/Basin' num2str(gauge1) '_output_900m1hr_S1_spinup/19730907/laststep_Qd1'];
fid=fopen(ffnm,'rb','ieee-le');
tmpdata=fread(fid,inf,'single');
tmp=reshape(tmpdata,c,r);
fclose(fid);
figure
imagesc(tmp')
colorbar
title('0907 Qd1')
caxis([0 10])
ffnm=['/shared/dondo/home/ml423/world_mts_outputs/Basin' num2str(gauge1) '_output_900m1hr_S1_spinup/19730908/laststep_Qu1'];
fid=fopen(ffnm,'rb','ieee-le');
tmpdata=fread(fid,inf,'single');
tmp=reshape(tmpdata,c,r);
fclose(fid);
figure
imagesc(tmp')
colorbar
title('0908 Qu1')
caxis([0 10])
ffnm=['/shared/dondo/home/ml423/world_mts_outputs/Basin' num2str(gauge1) '_output_900m1hr_S1_spinup/19730908/laststep_Qd1'];
fid=fopen(ffnm,'rb','ieee-le');
tmpdata=fread(fid,inf,'single');
tmp=reshape(tmpdata,c,r);
fclose(fid);
figure
imagesc(tmp')
colorbar
title('0908 Qd1')
caxis([0 10])







ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gauge1) '_ERA5_900m1hr/19730908/Precipitation_ERA5_land_oriunit'];
fid=fopen(ffnm,'rb','ieee-le');
tmpdata=fread(fid,inf,'single');
tmp=reshape(tmpdata,c,r,24);
fclose(fid);
for j = 1:9
tmpplot = tmp(:,:,j);
figure
imagesc(tmpplot')
title(['rain ' num2str(j) ])
end

ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gauge1) '_ERA5_900m1hr/19730908/DSWR_surfunit'];
fid=fopen(ffnm,'rb','ieee-le');
tmpdata=fread(fid,inf,'single');
tmp=reshape(tmpdata,c,r,24);
fclose(fid);
for j = 1:9
tmpplot = tmp(:,:,j);
figure
imagesc(tmpplot')
colorbar
title(['DSWR ' num2str(j) ])
end

ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gauge1) '_ERA5_900m1hr/19730908/DLWR_surfunit'];
fid=fopen(ffnm,'rb','ieee-le');
tmpdata=fread(fid,inf,'single');
tmp=reshape(tmpdata,c,r,24);
fclose(fid);
for j = 1:9
tmpplot = tmp(:,:,j);
figure
imagesc(tmpplot')
colorbar
title(['DLWR ' num2str(j) ])
end

ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gauge1) '_ERA5_900m1hr/19730908/AirPressure_10munit'];
fid=fopen(ffnm,'rb','ieee-le');
tmpdata=fread(fid,inf,'single');
tmp=reshape(tmpdata,c,r,24);
fclose(fid);
for j = 1:9
tmpplot = tmp(:,:,j);
figure
imagesc(tmpplot')
colorbar
title(['press ' num2str(j) ])
end

ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gauge1) '_ERA5_900m1hr/19730908/AirTemp_10m'];
fid=fopen(ffnm,'rb','ieee-le');
tmpdata=fread(fid,inf,'single');
tmp=reshape(tmpdata,c,r,24);
fclose(fid);
for j = 1:9
tmpplot = tmp(:,:,j);
figure
imagesc(tmpplot')
colorbar
title(['at ' num2str(j) ])
end

ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gauge1) '_ERA5_900m1hr/19730908/WindSpeed_10m'];
fid=fopen(ffnm,'rb','ieee-le');
tmpdata=fread(fid,inf,'single');
tmp=reshape(tmpdata,c,r,24);
fclose(fid);
for j = 1:9
tmpplot = tmp(:,:,j);
figure
imagesc(tmpplot')
colorbar
title(['wind ' num2str(j) ])
end

ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gauge1) '_ERA5_900m1hr/19730908/SpecHumi_10m'];
fid=fopen(ffnm,'rb','ieee-le');
tmpdata=fread(fid,inf,'single');
tmp=reshape(tmpdata,c,r,24);
fclose(fid);
for j = 1:9
tmpplot = tmp(:,:,j);
figure
imagesc(tmpplot')
colorbar
title(['sh ' num2str(j) ])
end


ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gauge1) '_MODIS_900m1hr/0908/Albedo'];
fid=fopen(ffnm,'rb','ieee-le');
tmpdata=fread(fid,inf,'single');
tmp=reshape(tmpdata,c,r,24);
fclose(fid);
for j = 1:9
tmpplot = tmp(:,:,j);
figure
imagesc(tmpplot')
title(['Albedo ' num2str(j) ])
end


ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gauge1) '_MODIS_900m1hr/0908/LAI'];
fid=fopen(ffnm,'rb','ieee-le');
tmpdata=fread(fid,inf,'single');
tmp=reshape(tmpdata,c,r,24);
fclose(fid);
for j = 1:9
tmpplot = tmp(:,:,j);
figure
imagesc(tmpplot')
title(['LAI ' num2str(j) ])
end

%% [data quality control for all basins all events all inputs] Forcing
% dicard data that has NaN, -1, -999, -9999, 9999, +inf, -inf in the data
varnm = {'AirTemp_10m','Precipitation_ERA5_land_oriunit','DSWR_surfunit','DLWR_surfunit','WindSpeed_10m','AirPressure_10munit','SpecHumi_10m'};
probmm = [];
wge = load(['/shared/dondo/home/ml423/DA/Critical matrices/world_selected_gauge_events_allinfo_v3.mat'],'gauge_events_first_last');
ge_info = wge.gauge_events_first_last;clear wge;
for bid = 134:465%472:522%:length(ge_info)
    gid = WORLD_events(bid,1);
    i = find(ge_info(:,1)==gid);
    
    ymd1 = ge_info(i,2);
    ymd2 = ge_info(i,3);
    r = ge_info(i,6);
    c = ge_info(i,7);

    if ymd1<19730911
        ymd1=19730911;
    end
    if ymd2<19740901
        disp('record not match study period')
        continue
    end
    if (ymd2-ymd1)<10000
        disp('record too short')
        continue
    end
    
    ymd1y = floor(ymd1/10000);
    ymd1d = mod(ymd1,100);
    ymd1m = floor((ymd1-ymd1y*10000)/100);
    ymd2y = floor(ymd2/10000);
    ymd2d = mod(ymd2,100);
    ymd2m = floor((ymd2-ymd2y*10000)/100);
    
    tt =  (datetime(ymd1y,ymd1m,ymd1d)-days(10)):days(1):(datetime(ymd2y,ymd2m,ymd2d)+days(9));
    tt = datetime(1973,1,1):days(1):datetime(2021,12,30);% to be general
    
    [yyl,mml,ddl] = ymd(tt);
    
    parfor dt = 1:length(yyl)
        
        dtm = yyl(dt)*10000+mml(dt)*100+ddl(dt);
        out_dir=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid,'%7.0d') '_ERA5_900m1hr/' num2str(dtm) '/'];
        if exist(out_dir,'dir')
        else
            %probmm = [probmm;[bid,dtm]];
            continue
        end
        
        for j = 1:7
            
            fnm = [out_dir,varnm{j}];
            staa = cell(4,1);
            fid = fopen([fnm,'_sta.txt'],'r');
            for ii = 1:3
                staa(ii) = {fgetl(fid)};
            end
            fclose(fid);
            if staa{3}==-1
                probmm = [probmm;[bid,dtm]];
               continue
               % meaning the txt file is not there,for what ever reason
            end
            meannum = staa{3}(10:12);
            meannum=str2double(meannum);
            
            if meannum>=0
                
            else
                probmm = [probmm;[bid,dtm,j]];
                %{
                disp(gid)
                disp(dtm)
                disp('readjust NaN')
                %probmm = [probmm;[gid,dtm,j]];
                ffnm=[out_dir varnm{j}];
                fid=fopen(ffnm,'rb','ieee-le');
                tmpdata=fread(fid,inf,'single');
                tmp=reshape(tmpdata,c,r,24);
                tmp_update = [];
                
                for k = 1:24
                    tmpp = tmp(:,:,k);
                    badloc = find(isnan(tmpp));
                    if isempty(badloc)
                    else
                        goodloc = find(~isnan(tmpp));
                        tmpp(badloc) = mean(tmpp(goodloc));
                    end
                    tmp_update(:,:,k)=tmpp;
                    
                end
                fnm=[out_dir varnm{j}];
                predata=reshape(tmp_update,[c,r*24]);
                fid = fopen(fnm,'wb');
                fwrite(fid,predata,'single');
                fclose(fid);
                
                fid=fopen([fnm,'_sta.txt'],'w');
                fprintf(fid,'Min    = %.4f\n',min(predata(:)));
                fprintf(fid,'Max    = %.4f\n',max(predata(:)));
                fprintf(fid,'Mean   = %.4f\n',mean(predata(:)));
                fprintf(fid,'Median = %.4f\n',median(predata(:)));
                fprintf(fid,'Std.   = %.4f\n',std(predata(:)));
                fclose(fid);
                %}
            end
            
        end
    end
    disp(i)
end

% probtmp = probmm;
%{
for ii = 1:length(probtmp)
    gid = WORLD_events(probtmp(ii,1),1);
    dtm = probtmp(ii,2);
    if mod(dtm,100)==1
        dtmtmp =dtm+1;
    else
        dtmtmp = dtm-1;
    end
    out_dir=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid,'%7.0d') '_ERA5_900m1hr/' num2str(dtmtmp) '/'];
    out_dirnew=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid,'%7.0d') '_ERA5_900m1hr/' num2str(dtm) '/'];

    i = find(ge_info(:,1)==gid);

    ymd1 = ge_info(i,2);
    ymd2 = ge_info(i,3);
    r = ge_info(i,6);
    c = ge_info(i,7);
    ffnm=[out_dir varnm{3}];
    fid=fopen(ffnm,'rb','ieee-le');
    tmpdata=fread(fid,inf,'single');
    tmp=reshape(tmpdata,c,r,24);

    fnm=[out_dirnew varnm{3}];
    predata=reshape(tmp,[c,r*24]);
    fid = fopen(fnm,'wb');
    fwrite(fid,predata,'single');
    fclose(fid);

    fid=fopen([fnm,'_sta.txt'],'w');
    fprintf(fid,'Min    = %.4f\n',min(predata(:)));
    fprintf(fid,'Max    = %.4f\n',max(predata(:)));
    fprintf(fid,'Mean   = %.4f\n',mean(predata(:)));
    fprintf(fid,'Median = %.4f\n',median(predata(:)));
    fprintf(fid,'Std.   = %.4f\n',std(predata(:)));
    fclose(fid);

end
%}

% fix 65 soil parameters



%% [data quality control for all basins all events all inputs] LULC
varnm = {'LAI','Emissivity','Albedo','CV'};
probmmlulc = [];
mday = [31,29,31,30,31,30,31,31,30,31,30,31];
mstr = [];
for mmo = 1:12
    
    for j = 1:mday(mmo)
        mstr = [mstr;mmo*100+j]; 
    end
end

for bid = 65%466:522%:length(ge_info)
    gid = WORLD_events(bid,1);
    i = find(ge_info(:,1)==gid);
 
    ymd1 = ge_info(i,2);
    ymd2 = ge_info(i,3);
    r = ge_info(i,6);
    c = ge_info(i,7);
    

    
    if ymd1<19730911
        ymd1=19730911;
    end
    if ymd2<19740901
        disp('record not match study period')
        continue
    end
    if (ymd2-ymd1)<10000
        disp('record too short')
        continue
    end
    
    
    
    parfor dt = 1:366
        
        dtm = mstr(dt);
        out_dir=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid,'%7.0d') '_MODIS_900m1hr/' num2str(dtm,'%4.4d') '/'];
        for j = 1:4

            fnm = [out_dir,varnm{j}];
            staa = cell(4,1);
            fid = fopen([fnm,'_sta.txt'],'r');
            for ii = 1:3
                staa(ii) = {fgetl(fid)};
            end
            fclose(fid);
           
            meannum = staa{3}(10:12);
            meannum=str2double(meannum);
            
            if meannum>=0
                
            else
                probmmlulc = [probmmlulc;[gid,dtm,j]];
                
                %{
                disp('readjust NaN')
                %probmm = [probmm;[gid,dtm,j]];
                ffnm=[out_dir varnm{j}];
                fid=fopen(ffnm,'rb','ieee-le');
                tmpdata=fread(fid,inf,'single');
                tmp=reshape(tmpdata,c,r,24);
                tmp_update = [];
                
                for k = 1:24
                    tmpp = tmp(:,:,k);
                    badloc = find(isnan(tmpp));
                    if isempty(badloc)
                    else
                        goodloc = find(~isnan(tmpp));
                        tmpp(badloc) = mean(tmpp(goodloc));
                    end
                    tmp_update(:,:,k)=tmpp;
                    
                end
                fnm=[out_dir varnm{j}];
                predata=reshape(tmp_update,[c,r*24]);
                fid = fopen(fnm,'wb');
                fwrite(fid,predata,'single');
                fclose(fid);
                
                fid=fopen([fnm,'_sta.txt'],'w');
                fprintf(fid,'Min    = %.4f\n',min(predata(:)));
                fprintf(fid,'Max    = %.4f\n',max(predata(:)));
                fprintf(fid,'Mean   = %.4f\n',mean(predata(:)));
                fprintf(fid,'Median = %.4f\n',median(predata(:)));
                fprintf(fid,'Std.   = %.4f\n',std(predata(:)));
                fclose(fid);
                %}
            end
            
        end
    end
    
end

%% [key based on selected_gauge_events] determine the final events

% 1. Streamflow gauge observations has to be greater than 0
% 2. Events have to be after 19730901
% 3. Basin size has to be less than 4000km2/1200km2 for the paper

tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/world_selected_gauge_events_v3.mat'],'selected_gauge_events');
selected_gauge_events = tmp.selected_gauge_events;clear tmp
tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/world_selected_gauge_events_allinfo_v3.mat'],'gauge_events_first_last');
ge_info = tmp.gauge_events_first_last;% this variable is just a reformat of the selected_gauge_events variable changing from cell variable to matrix variable
WORLD_events = nan(length(ge_info),300);
enid = 0;
biid = 0;
for mount = 1
    for bas = 1:size(selected_gauge_events,2)
        disp(bas)
        disp(mount)
        if size(selected_gauge_events{mount,bas},2)>=2
            % means this basin actually has events
            entry = selected_gauge_events{mount,bas};
            g_id = str2double(entry(1));

            if ismember(g_id,WORLD_eventsold)
                continue
            end
            tmp = load(['/shared/dondo/home/ml423/world_mts/obs_events/obs_' num2str(g_id,'%7.0d') '.mat'],'obs');
            obsgid = tmp.obs; clear tmp
            gloc = find(ge_info(:,1)==g_id);

            if ge_info(gloc,5)<=4000
                % basins bigger than 4000km2 are abandoned, later changed
                % to 1200 for the paper
                biid = biid+1;
                enid = 0;
                for eid = 1:size(selected_gauge_events{mount,bas},2)-1
                    tmpe = char(entry(1+eid));
                    yth = str2num(tmpe(1:4));
                    mth = str2num(tmpe(6:7));
                    dth = str2num(tmpe(9:10));
                    thresdate = datetime(yth,mth,dth);edate = yth*10000+mth*100+dth;
                    eventbd1 = thresdate-days(15);[bd1y,bd1m,bd1d] = ymd(eventbd1);bd1str = bd1y*10000+bd1m*100+bd1d;
                    eventbd2 = thresdate+days(15);[bd2y,bd2m,bd2d] = ymd(eventbd2);bd2str = bd2y*10000+bd2m*100+bd2d;
                    bd1loc = find(obsgid(:,1)==bd1str);bd2loc = find(obsgid(:,1)==bd2str);
                    obsperiod = obsgid(bd1loc:bd2loc,2);
                    if bd1str>19730901% start date has to be after 19730901.
                        if bd2str<20211231% ERA5 is processed till 20211231
                            if min(obsperiod)>-0.0001 % can not have NaN flow observations
                                enid = enid+1;
                                WORLD_events(biid,1) = g_id;
                                WORLD_events(biid,1+enid) = edate;

                            end

                        end
                    end


                end

            end

        end
    end
end

noeventloc = find(isnan(WORLD_events(:,2)));% This means this basin has no event

WORLD_events(noeventloc,:) = [];

%% write out 30day stream obs, rainfall into the IRC folders for calculating rainfall uncertainty


tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/world_selected_gauge_events_allinfo_v3.mat'],'gauge_events_first_last');
ge_info = tmp.gauge_events_first_last;clear tmp


for i = 1:size(WORLD_events,1)
    gid = WORLD_events(i,1);

    g_loc = find(ge_info(:,1)==gid);
    r = ge_info(g_loc,6);
    c = ge_info(g_loc,7);
    % g_loc = find(Basins_sum_region(:,6)==gid);
    % r = Basins_sum_region(g_loc,1);
    % c = Basins_sum_region(g_loc,2);
    tmp = load(['/shared/dondo/home/ml423/world_mts/obs_events/obs_' num2str(gid) '.mat']);
    tmp = tmp.obs;
    obsgid = tmp;
    disp(i)
    for j = 1:size(WORLD_events,2)-1
        event = WORLD_events(i,1+j);
%         if isfile(['/shared/dondo/home/ml423/world_mts/obs_events/events/obs_' num2str(gid) '_event' num2str(event) '.mat'])
%             continue
%         else
%             disp(i)
%             disp(j)
%         end
        
        if ~isnan(event)
            yth = floor(event/10000);
            dth = mod(event,100);
            mth = floor((event-yth*10000)/100);
            
            thresdate = datetime(yth,mth,dth);edate = yth*10000+mth*100+dth;
            eventbd1 = thresdate-days(15);[bd1y,bd1m,bd1d] = ymd(eventbd1);bd1str = bd1y*10000+bd1m*100+bd1d;
            eventbd2 = thresdate+days(15);[bd2y,bd2m,bd2d] = ymd(eventbd2);bd2str = bd2y*10000+bd2m*100+bd2d;
            % This is a 31 day window defined above
            bd1loc = find(obsgid(:,1)==bd1str);bd2loc = find(obsgid(:,1)==bd2str);
            
            obsnow = obsgid(bd1loc:bd2loc,2);
            save(['/shared/dondo/home/ml423/world_mts/obs_events/events/obs_' num2str(gid) '_event' num2str(event)],'obsnow')
            
            durati = (thresdate-days(15)):days(1):(thresdate+days(15));
            [dury,durm,durd] = ymd(durati);
            p_event = {};
            p_count = 0;
            for kk = 1:length(dury)
                p_count = p_count+1;
               pda = dury(kk)*10000+durm(kk)*100+durd(kk); 
               ffnm = ['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_ERA5_900m1hr/' num2str(pda) '/Precipitation_ERA5_land_oriunit'];
               fid=fopen(ffnm,'rb','ieee-le');
               tmpdata=fread(fid,inf,'single');fclose(fid);
               pra=reshape(tmpdata,c,r,24);
               for kk2 = 1:24
                p_event{(p_count-1)*24+kk2} = pra(:,:,kk2);
               end
            end
            rain = p_event(1:720);% get the 30 days rainfall
            
            save(['/shared/dondo/home/ml423/world_mts/Precip/Event/ERA5_Basin' num2str(gid) '_' num2str(event) '.mat'],'rain')
            
            ounm = {'Precipitation_ERA5_land_oriunit'};
            
            out_dirp = ['/shared/dondo/home/ml423/world_mts/Precip/Event/Basin' num2str(gid) '_input_ERA5_900m_FIX/' num2str(event) '/'];
            mkdir(out_dirp)
            dataf = [];
            for j = 1:720
                datas = rain{j}; % can add stuff
                dataf(:,:,j) = datas;
            end
            dataf(dataf<0) = 0;
            fnm = [out_dirp,ounm{1}];
            raindata = reshape(dataf,[c,r,720]);
            fid = fopen(fnm,'wb');
            fwrite(fid,raindata,'single');
            fclose(fid);
            fid = fopen([fnm,'_sta.txt'],'w');
            fprintf(fid,'Min    = %.4f\n',min(raindata(:)));
            fprintf(fid,'Max    = %.4f\n',max(raindata(:)));
            fprintf(fid,'Mean    = %.4f\n',mean(raindata(:)));
            fprintf(fid,'Std    = %.4f\n',std(raindata(:)));
            fclose(fid);

        else
            break
        end
    end
end 
%% write out 30day MODIS, other Forcings into the IRC folders for calculating rainfall uncertainty

tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/world_selected_gauge_events_allinfo_v3.mat'],'gauge_events_first_last');
ge_info = tmp.gauge_events_first_last;clear tmp


for i = 500%466:size(WORLD_events,1)
    gid = WORLD_events(i,1);

    g_loc = find(ge_info(:,1)==gid);
    r = ge_info(g_loc,6);
    c = ge_info(g_loc,7);
    % g_loc = find(Basins_sum_region(:,6)==gid);
    % r = Basins_sum_region(g_loc,1);
    % c = Basins_sum_region(g_loc,2);

    disp(i)
    %parfor j = 1:size(WORLD_events,2)-1
    for j = 29
        event = WORLD_events(i,1+j);
        % if isfile(['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_input_MODIS_900m_tmp/' num2str(event) '/Albedo'])
        %     continue
        % else
        %     disp(i)
        %     disp(j)
        % end
        
        if ~isnan(event)
            yth = floor(event/10000);
            dth = mod(event,100);
            mth = floor((event-yth*10000)/100);

            thresdate = datetime(yth,mth,dth);edate = yth*10000+mth*100+dth;

            durati = (thresdate-days(15)):days(1):(thresdate+days(15));
            [dury,durm,durd] = ymd(durati);
            ounm = {'AirTemp_10m','DLWR_surfunit','DSWR_surfunit','WindSpeed_10m','SpecHumi_10m','AirPressure_10munit'};

            for var = 3%1:6
                p_event = {};
                p_count = 0;
                for kk = 1:length(dury)
                    p_count = p_count+1;
                    pda = dury(kk)*10000+durm(kk)*100+durd(kk);
                    ffnm = ['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_ERA5_900m1hr/' num2str(pda) '/' ounm{var}];
                    fid=fopen(ffnm,'rb','ieee-le');
                    tmpdata=fread(fid,inf,'single');fclose(fid);
                    pra=reshape(tmpdata,c,r,24);
                    for kk2 = 1:24
                        p_event{(p_count-1)*24+kk2} = pra(:,:,kk2);
                    end
                end
                rain = p_event(1:720);% get the 30 days rainfall

                out_dirp = ['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_input_ERA5_900m_tmp/' num2str(event) '/'];
                mkdir(out_dirp)
                dataf = [];
                for jk = 1:720
                    datas = rain{jk}; % can add stuff
                    dataf(:,:,jk) = datas;
                end
                dataf(dataf<0) = 0;
                fnm = [out_dirp,ounm{var}];
                raindata = reshape(dataf,[c,r,720]);
                fid = fopen(fnm,'wb');
                fwrite(fid,raindata,'single');
                fclose(fid);
                fid = fopen([fnm,'_sta.txt'],'w');
                fprintf(fid,'Min    = %.4f\n',min(raindata(:)));
                fprintf(fid,'Max    = %.4f\n',max(raindata(:)));
                fprintf(fid,'Mean    = %.4f\n',mean(raindata(:)));
                fprintf(fid,'Std    = %.4f\n',std(raindata(:)));
                fclose(fid);
                %
            end
            %{

            ounmla = {'Emissivity','LAI','CV','Albedo'};

            for var = 1:4
                p_event = {};
                p_count = 0;
                for kk = 1:length(dury)
                    p_count = p_count+1;
                    pda = durm(kk)*100+durd(kk);
                    ffnm = ['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_MODIS_900m1hr/' num2str(pda,'%4.4d') '/' ounmla{var}];
                    fid=fopen(ffnm,'rb','ieee-le');
                    tmpdata=fread(fid,inf,'single');fclose(fid);
                    pra=reshape(tmpdata,c,r,24);
                    for kk2 = 1:24
                        p_event{(p_count-1)*24+kk2} = pra(:,:,kk2);
                    end
                end
                rain = p_event(1:720);% get the 30 days rainfall

                out_dirp = ['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_input_MODIS_900m_tmp/' num2str(event) '/'];
                mkdir(out_dirp)
                dataf = [];
                for jk = 1:720
                    datas = rain{jk}; % can add stuff
                    dataf(:,:,jk) = datas;
                end
                dataf(dataf<0) = 0;
                fnm = [out_dirp,ounmla{var}];
                raindata = reshape(dataf,[c,r,720]);
                fid = fopen(fnm,'wb');
                fwrite(fid,raindata,'single');
                fclose(fid);
                fid = fopen([fnm,'_sta.txt'],'w');
                fprintf(fid,'Min    = %.4f\n',min(raindata(:)));
                fprintf(fid,'Max    = %.4f\n',max(raindata(:)));
                fprintf(fid,'Mean    = %.4f\n',mean(raindata(:)));
                fprintf(fid,'Std    = %.4f\n',std(raindata(:)));
                fclose(fid);
                
            end
            %}


        else
            %break
        end
    end
end
