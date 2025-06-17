wge = load(['/shared/dondo/home/ml423/DA/Critical matrices/world_selected_gauge_events_allinfo_v3.mat'],'gauge_events_first_last');
ge_info = wge.gauge_events_first_last;clear wge;
%flow_vol_r = [];
plotfig = 1;

tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/WORLD_events_v3.mat'],'WORLD_events');
WORLD_events = tmp.WORLD_events;




fid = fopen('basin_code.txt','r');
fom = '%d';
bid = fscanf(fid,fom);
fclose(fid);
fid = fopen('cloop.txt','r');
fom = '%d';
loopitr = fscanf(fid,fom);
fclose(fid);

for g_shell_index = bid%[1:22 24:105 107:160]% problem 23 106
    close all
    preset_gauge = WORLD_events(g_shell_index,1);%6559101;%6935540;%;%6235535;%6559100;%6125361;
    preset_gauge
    suffix = 'S9';
    spinitr = [1];% spinup iteration to be plotted
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
            s_flow = [s_flow;flow'];
            
            % obs
            obsfile = ['/shared/dondo/home/ml423/world_mts/obs_events/days/obs_' num2str(gid,'%7.0d') '/' daystr '.mat'];
            if isfile(obsfile)
                tmp = load(obsfile,'obsnow');
                obs = tmp.obsnow;
                obs_flow = [obs_flow;obs];
            else
                obs_flow = [obs_flow;nan];
            end
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
        
        n_ind = ~isnan(obs_flow);
        n_loc = find(n_ind>0);
        
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
        flow_vol_r = sum(simuna)/sum(obsnowx);

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
        
        % if flow_vol_r>0.8&&flow_vol_r<1.2
        %    changecoeflater = 1;
        % elseif flow_vol_r>=1.2
        %     changecoeflater = 2; 
        % elseif flow_vol_r<=0.8
        %    changecoeflater = 0; 
        % end

        if flowpeak_ratio>-0.33&&flowpeak_ratio<0.33
           changecoeflater = 1;
        elseif flowpeak_ratio>=0.33
            changecoeflater = 2; 
        elseif flowpeak_ratio<=-0.33
           changecoeflater = 0; 
        end
        
        fid = fopen('update_coef.txt','w');
        fprintf(fid,'%d',changecoeflater);
        fclose(fid);
        fid = fopen(['flow_pvolume_ratio_itr' num2str(loopitr) '.txt'],'w');
        fprintf(fid,'%2.3f',flowpeak_ratio);
        fclose(fid);
        
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

            saveas(gcf,[picdir, 'Basin' num2str(gid,'%7.0d') '_' num2str(nsm) '_' num2str(nm) '_spinup-' flowtype{ftp} '-' num2str(spinup) '-' suffix 'wTemp_lateral_iter_' num2str(loopitr)],'png')
        end
    end
    
end
quit