% this script is to plot results from IRC
%% map APGD data to each basins
working_dir = '/shared/dondo/home/ml423/DA/alps_rg_based_qpe/';
ncdisp([working_dir 'RapdD_al05.etrs.laea_197908.nc'])

ncfile = [working_dir 'RapdD_al05.etrs.laea_199409.nc'];
lon = ncread(ncfile,'lon') ; 
lat = ncread(ncfile,'lat') ; 
rapdd = ncread(ncfile,'RapdD') ; 

figure
imagesc(rapdd(:,:,21))
colorbar
% for each basin, find the nearest mapping, nearest neighbor
home_alps = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Alps/input_output/';
home_andes = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Andes/input_output/';
home_himalayas = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Himalayas/input_output/';

home_moun = home_alps;mtname = 'Alps';
tmp = load([home_moun 'gauges_loc/Basins_dimention.mat'],'Basins_sum');
Basins_sum_region = tmp.Basins_sum;
bnum = size(Basins_sum_region,1);

bidapgdloc = cell(522,1);
for bid = 1:bnum
    cols = Basins_sum_region(bid,2);
    rows = Basins_sum_region(bid,1);
    gid = Basins_sum_region(bid,6);

    gbid = find(WORLD_events==gid);

    if isempty(gbid)
        continue
    end
    outdir = ['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_APGD_5km1day'];
    mkdir(outdir)

    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/latlon/Basin' num2str(gid) '_lat.mat'],'lat_b');
    basin_lat = tmp.lat_b;
    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/latlon/Basin' num2str(gid) '_lon.mat'],'lon_b');
    basin_lon = tmp.lon_b;

    lat1d = basin_lat(:);lon1d = basin_lon(:);% one can adjust the lat lon here to get nearby lat lon if nan is present

    if max(lat1d)<min(lat(:))||max(lon1d)<min(lon(:))
        continue
    end

    loc_apgd = [];
    for k = 1:length(lat1d)
        lsquare = (lat-lat1d(k)).^2+(lon-lon1d(k)).^2;
        [~,nnloc] = min(lsquare(:));
        loc_apgd(k) = nnloc;
        if lat1d(k)<min(lat(:))||lon1d(k)<min(lon(:))
            loc_apgd(k) = 1;% 1 indicates out of apgd area
        end

    end
    bidapgdloc{gbid} = loc_apgd;
end



imonthdays_leap = [31 29 31 30 31 30 31 31 30 31 30 31];
imonthdays_normal = [31 28 31 30 31 30 31 31 30 31 30 31];

for iyear = 2019:-1:1971
    %parpool('local',12)
    if mod(iyear,4)==0
        imonthdays = imonthdays_leap;
    else
        imonthdays = imonthdays_normal;
    end
    for imonth = 1:12
        disp(iyear)
        disp(imonth)

        home_moun = home_alps;mtname = 'Alps';
        tmp = load([home_moun 'gauges_loc/Basins_dimention.mat'],'Basins_sum');
        Basins_sum_region = tmp.Basins_sum;
        bnum = size(Basins_sum_region,1);

        ncfile = [working_dir 'RapdD_al05.etrs.laea_' num2str(iyear) num2str(imonth,'%2.2d') '.nc'];

        rapdd = ncread(ncfile,'RapdD') ;

        parfor bid = 1:bnum

            cols = Basins_sum_region(bid,2);
            rows = Basins_sum_region(bid,1);
            gid = Basins_sum_region(bid,6);

            gbid = find(WORLD_events==gid);

            if isempty(gbid)||isempty(bidapgdloc{gbid})
                continue
            end


            for iday = 1:imonthdays(imonth)

                outdir = ['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_APGD_5km1day'];

                outfiledir = [outdir  '/' num2str(iyear) num2str(imonth,'%2.2d') num2str(iday,'%2.2d') '/'];
                mkdir(outfiledir)

                tmpr = rapdd(:,:,iday);
                basinr = tmpr(bidapgdloc{gbid});
                outbd = find(bidapgdloc{gbid}==1);
                basinr(outbd) = nan;
                basinr2d=reshape(basinr,rows,cols);
                bfinal =  basinr2d';% transpose to fortran style in dchm


                fnm=[outfiledir,'APGD_QPE'];
                predata=reshape(bfinal,[cols,rows]);
                fid = fopen(fnm,'wb');
                fwrite(fid,predata,'single');
                fclose(fid);

                fid=fopen([fnm,'_sta.txt'],'w');
                fprintf(fid,'Min    = %.4f\n',min(predata(:)));
                fprintf(fid,'Max    = %.4f\n',max(predata(:)));
                fprintf(fid,'Mean   = %.4f\n',mean(predata(:)));
                fprintf(fid,'Median = %.4f\n',nanmedian(predata(:)));
                fprintf(fid,'Std.   = %.4f\n',std(predata(:)));
                fclose(fid);

            end
        end
    end

end
%% calculate basin averaged APGD rain
tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/world_selected_gauge_events_allinfo_v3.mat'],'gauge_events_first_last');
ge_info = tmp.gauge_events_first_last;clear tmp
apgd_qpe = cell(522,150);
apgd_mean = nan(522,150);

for i = [1:133 466:522]
    %466:size(WORLD_events,1)
    gid = WORLD_events(i,1);

    g_loc = find(ge_info(:,1)==gid);
    r = ge_info(g_loc,6);
    c = ge_info(g_loc,7);
    disp(i)
    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_ind.mat'],'ind');
    ind = tmp.ind;
    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_intt.mat'],'intt');
    intt = tmp.intt;

    %parfor j = 1:size(WORLD_events,2)-1
    ixx = find(~isnan(WORLD_events(i,:)));
    etot = length(ixx)-1;
    if etot>150
        etot = 150;
    end
    for j = 1:etot
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
            ounm = {'APGD_QPE'};

            paccu = 0;
            pqpe = zeros(c,r);
            for kk = 1:30
                %length(dury)
                pda = dury(kk)*10000+durm(kk)*100+durd(kk);
                ffnm = ['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_APGD_5km1day/' num2str(pda) '/' ounm{1}];
                fid=fopen(ffnm,'rb','ieee-le');
                if fid ==-1
                    apgd_qpe{i,j} = nan(c,r);
                    apgd_mean(i,j) = nan;
                    break
                end
                tmpdata=fread(fid,inf,'single');fclose(fid);
                pra=reshape(tmpdata,c,r);
                pavg = nanmean(pra(intt));
                paccu = paccu+pavg;

                pqpe = pqpe+pra;
                apgd_qpe{i,j} = pqpe;
                apgd_mean(i,j) = paccu;
            end
            

        end
    end

end
%% prepare events info for using by shell scripts


for bid = 1:522

    gid = WORLD_events(bid,1);
    gloc = find(gauge_events_first_last(:,1)==gid);
    area = gauge_events_first_last(gloc,5);

    ixx = find(~isnan(WORLD_events(bid,:)));
    etot = length(ixx)-1;
    if etot == 0
        binfo = [gid;etot];
        ffnm = ['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/basins_events/'];
        fid = fopen([ffnm 'Basin_index_' num2str(bid,'%3.3d') '.txt'],'w');
        fprintf(fid,'%d\n',binfo);
        fclose(fid);
    else
        binfo = [gid;etot];
        ffnm = ['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/basins_events/'];
        fid = fopen([ffnm 'Basin_index_' num2str(bid,'%3.3d') '.txt'],'w');
        fprintf(fid,'%d\n',binfo);
        fclose(fid);

        for eid = 1:etot
            %[1:15]%[1 2 3 4 6 8]%:5%1:etot%[5 6 7]%[1 2 3 4 5 6 7 8 9 10 11]
            close all
            if ~isnan(event_quality(bid,eid+1))
                event =  WORLD_events(bid,eid+1);

                tmp = load(['/shared/dondo/home/ml423/world_mts/obs_events/events/obs_' num2str(gid) '_event' num2str(event) '.mat']);
                tmp = tmp.obsnow;
                obsnow = tmp;clear tmp

                obsnan=find(isnan(obsnow));
                if isempty(obsnan)
                else
                    for na = 1:length(obsnan)
                        obsnow(obsnan(na)) = 0.5*(obsnow(obsnan(na)-1)+obsnow(obsnan(na)+1));
                    end
                end

                obs2 = interp1(1:31,obsnow,[1:0.33333:31]);
                obsnow = obs2(1:90);

                ffnm = ['/shared/dondo/home/ml423/world_mts_runs/Basins_events_chars/'];
                fid = fopen([ffnm 'Basin' num2str(gid) '_event' num2str(event) '.txt'],'r');
                if fid == -1
                    
                    eventinfo = [0;0;0;0];
                    ffnm = ['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/basins_events/'];
                    fid = fopen([ffnm 'Basin_' num2str(gid) '_event_' num2str(eid) '.txt'],'w');
                    fprintf(fid,'%d\n',eventinfo);
                    fclose(fid);

                    continue
                end
                windowtimes = fscanf(fid,'%d\n',[4,1]);
                fclose(fid);
                if length(windowtimes)<4
                    eventinfo = [0;0;0;0];
                    ffnm = ['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/basins_events/'];
                    fid = fopen([ffnm 'Basin_' num2str(gid) '_event_' num2str(eid) '.txt'],'w');
                    fprintf(fid,'%d\n',eventinfo);
                    fclose(fid);
                    continue
                end

                tmp = load(['/shared/dondo/home/ml423/world_mts/Precip/Event/ERA5_Basin' num2str(gid) '_' num2str(event) '.mat']);
                era5 = tmp.rain;
                rmax = zeros(1,720);
                for te = 1:720
                    rmax(te) = max(era5{te}(intt));
                end
                launchp = find(rmax>0.5);% when rainfall >1mm/hr launch particles
                ptmaxthre = 0.5*max(rmax);% used as a threshold to separate events
                la2 = find(launchp<(windowtimes(4)*8));
                la1 = find(launchp<(windowtimes(3)*8));
                accu0 = 0;
                if isempty(la1)
                    eventinfo = [0;0;0;0];
                    ffnm = ['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/basins_events/'];
                    fid = fopen([ffnm 'Basin_' num2str(gid) '_event_' num2str(eid) '.txt'],'w');
                    fprintf(fid,'%d\n',eventinfo);
                    fclose(fid);
                    continue
                end
                for tmpla = launchp(la1(end)):-1:1
                    if ismember(tmpla,launchp)
                        accu0 = 0;
                        rainbegintime = tmpla;
                    else
                        accu0 = accu0+1;
                    end

                    if area>1000
                        noraininterval = 12;% for big basins, it is likely there is always rain in the basin, therefore using a small norain limit to distinguish different rain events.
                    else
                        noraininterval = 24;
                    end
                    if accu0>noraininterval
                        % consecutive x hours of basin max rain less than 0.5mm/hr,
                        % this should be counted as a different rain cell, should not
                        % be tracked!
                        % one exception is if another major rain is 48-24 ahead, should be
                        % considered as well
                        rsmallest = max([(tmpla-noraininterval),1]);

                        reconsider_period = rsmallest:tmpla;
                        %reconsider_period = (tmpla-noraininterval):tmpla;
                        secondary_rain = find(rmax(reconsider_period)>0.5);

                        if isempty(secondary_rain)
                            eventinfo = [0;0;0;0];
                            ffnm = ['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/basins_events/'];
                            fid = fopen([ffnm 'Basin_' num2str(gid) '_event_' num2str(eid) '.txt'],'w');
                            fprintf(fid,'%d\n',eventinfo);
                            fclose(fid);
                            continue
                        end

                        if max(rmax(reconsider_period))>ptmaxthre
                            % adjust the rainbegintime here
                            rainbegintime = reconsider_period(secondary_rain(1));
                        end
                        smalltmpla = max([floor(tmpla/8),2]);
                        if obsnow(smalltmpla)>obsnow(smalltmpla-1)
                            % meaning obsflow is still rising.
                            if ~isempty(secondary_rain)
                                rainbegintime = reconsider_period(secondary_rain(1));
                            end
                        end

                        break
                    end
                end

                la0 = find(launchp == rainbegintime);

                refinedlaunchp = launchp(la0:la2(end));


                evnt = WORLD_events(bid,1+eid);
                ymd2y = floor(evnt/10000);
                ymd2d = mod(evnt,100);
                ymd2m = floor((evnt-ymd2y*10000)/100);
                icdate = datetime(ymd2y,ymd2m,ymd2d)-days(16);% careful here! the intial condition is the file from 16 days ago not 15
                [icy,icm,icd] = ymd(icdate);
                ic = icy*10000+icm*100+icd;

                eventinfo = [1;event;ic;refinedlaunchp(1)];% quality,current event, IC date, trueIC location
            else
                eventinfo = [0;0;0;0];

            end
       
            ffnm = ['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/basins_events/'];
            fid = fopen([ffnm 'Basin_' num2str(gid) '_event_' num2str(eid) '.txt'],'w');
            fprintf(fid,'%d\n',eventinfo);
            fclose(fid);
        end
        disp(bid);
    end
end


% for i = 1:150
%     [close figure i]
% end
%{
r00 = [];
rop = [];ropold = [];
r00rg = [];
r00p = [];% use old IRC results where 0.0999 is
ropp = [];
r00kge = [];
ropkge = [];
r00pn = [];%use  new IRC results where 0.000999 is used for IRC basethick
roppn = [];
r00kgen = [];
ropkgen = [];
for bid = [507:522]%[50 79 91 93 519 92 94 89 110 472 481 484 490 491 494 104]
    ixx = find(~isnan(WORLD_events(bid,:)));
    etot = length(ixx)-1;
    if etot == 0
        continue
    elseif etot>149
        etot=149;
    end
    ffnm = ['/shared/dondo/home/ml423/world_mts_runs/IRC1/'];
    fid = fopen([ffnm 'Basin' num2str(WORLD_events(bid,1)) '_' num2str(bid) '.txt'],'r');
    lastupdate = fscanf(fid,'%d');
    fclose(fid);
    for eid = 1:etot%lastupdate-1
        r00 = [r00,rain_orim(bid,eid)];
        ropold = [ropold,rain_optmold(bid,eid)];
        rop = [rop,rain_optm(bid,eid)];
        loc = find(rain_compareold(:,1)==bid&rain_compareold(:,2)==eid);
        locn = find(rain_compare(:,1)==bid&rain_compare(:,2)==eid);
        if ~isempty(loc)
            r00rg = [r00rg;rain_compareold(loc,6)];
            r00p = [r00p;rain_compareold(loc,7)];
            ropp = [ropp;rain_compareold(loc,8)];
            r00kge = [r00kge;rain_compareold(loc,9)];
            ropkge = [ropkge;rain_compareold(loc,10)];

            r00pn = [r00pn;rain_compare(locn,7)];
            roppn = [roppn;rain_compare(locn,8)];
            r00kgen = [r00kgen;rain_compare(locn,9)];
            ropkgen = [ropkgen;rain_compare(locn,10)];
            if max(rain_compare(locn,8))>500
                disp([bid,eid,max(rain_compare(locn,8))])
            end
        end
    end
end
figure
scatter(r00,ropold,'ro')
hold on
scatter(r00,rop,'b+')
hold on
plot([0 500],[0 500],'r-')

figure
scatter(r00rg,r00p,'r+')
hold on
scatter(r00rg,ropp,'b.')
hold on
plot([0 500],[0 500],'k-')

figure
scatter(r00rg,r00pn,'r+')
hold on
scatter(r00rg,roppn,'b.')
hold on
plot([0 500],[0 500],'k-')


figure
scatter(r00kgen,ropkgen,'r.')
hold on
plot([-1 1],[-1 1],'k-')
xlim([-2 1])
%}

lat = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/latlon/Basin' num2str(WORLD_events(91,1)) '_lat.mat'],'lat_b');
basin_lat = lat.lat_b;
lon = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/latlon/Basin' num2str(WORLD_events(91,1)) '_lon.mat'],'lon_b');
basin_lon = lon.lon_b;
figure
imagesc(basin_lat)
figure
imagesc(basin_lon)


figure
scatter(rain_orim(519,1:62),rain_optmold(519,1:62))
hold on
scatter(rain_orim(519,1:62),rain_optm(519,1:62),12,kge_optm(519,1:62),'filled')
hold on 
plot([0 500],[0 500],'r')
clim([0.2 1])

%% event quality, make sure events statisfy 'flash floods' and eliminate duplicate events from events selection

% event_flash = [];

for bid = [1:465]%[1:133 466:522]%tmpjjj = 1:length(bidt)% [1:133 466:522]%1:length(reprob)
    disp(bid)
    %bid = reprob(prob,1);%[74 59 517 60 25 121 69 512 502 107 519 82 521 70 56]
    %bid = bidt(tmpjjj);
    %[1:133 466:522]
    %tmpii = 1:length(listrg) %[50 79 91 93 519 92 94 104 89 110 472 481 484 490 491 494]%[1:133 466:522]%[268:300]%37:100%26:100%[ 355 359 86 22 360] %98:190
    %bid = listrg(tmpii);
    %close all
    % kge_optm(bid,1:149) = nan;
    % kge_orim(bid,1:149) = nan;
    % rain_orim(bid,1:149) = nan;
    % rain_optm(bid,1:149) = nan;
    % rainsd_orim(bid,1:149) = nan;
    % rainsd_optm(bid,1:149) = nan;


    gid = WORLD_events(bid,1);

    gloc = find(Basins_sum(:,1)==gid);

    r = Basins_sum(gloc,6);
    c = Basins_sum(gloc,7);
    rg = Basins_sum(gloc,8);
    cg = Basins_sum(gloc,9);
    area = Basins_sum(gloc,5);

    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_ind.mat'],'ind');
    ind = tmp.ind;
    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_intt.mat'],'intt');
    intt = tmp.intt;


    ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_fixed_900m/dem.bin'];
    fid=fopen(ffnm,'rb','ieee-le');
    tmpdata=fread(fid,inf,'single');
    fclose(fid);
    dem = reshape(tmpdata,c,r);
    dem(ind) = nan;
    demf = dem';
    colo = slanCM('terrain',80);
    colo(1,:) = [0.8 0.8 0.8];
    figure
    imagesc(demf)
    colormap(colo)

    ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_fixed_900m/countfile.out'];
    fid=fopen(ffnm);
    facc=fscanf(fid,'%d',[c,r]);
    fclose(fid);
    facc(ind) = nan;
    [streamx,streamy] = find(facc>5);


    ICn = 4000;
    ICn = 4140;
    ICC = 0;% IF ICC is part of the algo

    fileind = -100000*ICn;
    KGEsu = [];
    evesu = [];
    budgetalldbkc = [];
    budgetallirc = [];

    ixx = find(~isnan(WORLD_events(bid,:)));
    etot = length(ixx)-1;

    % check if this basin has rgs

    if etot == 0
        continue
    end


    for eid = 1:etot%eidt(tmpjjj)%1:etot
        %reprob(prob,2)%1:etot%[1:15]%[1 2 3 4 6 8]%:5%1:etot%[5 6 7]%[1 2 3 4 5 6 7 8 9 10 11]
        close all
        %rglists = find(rgr(:,1)==bid&rgr(:,2)==eid);
        %rgeve = rgr(rglists,5:7);

        if isnan(event_quality(bid,eid+1))
            disp('Event quality is bad, peak too early or too late')
            continue
        else
            rainevent = num2str(WORLD_events(bid,eid+1));
            yth = floor(WORLD_events(bid,eid+1)/10000);
            dth = mod(WORLD_events(bid,eid+1),100);
            mth = floor((WORLD_events(bid,eid+1)-yth*10000)/100);

            thresdate = datetime(yth,mth,dth);

            durati = (thresdate-days(15)):days(1):(thresdate+days(15));
            [dury,durm,durd] = ymd(durati);
            tst1 = num2str(durm(1)*100+durd(1),'%4.4d');
            tst2 = num2str(durm(30)*100+durd(30),'%4.4d');

            xtlabels ={};
            for tl = 1:30
                xtlabels{tl} =  [ num2str(durm(tl),'%2.2d') '/' num2str(durd(tl),'%2.2d') ];
            end



            ffnm=['/shared/dondo/home/ml423/world_mts_runs/IRC1/Basin' num2str(gid) 'outputs/refxtmp' num2str(str2num(rainevent)*100+fileind+0) 'fullresults/' rainevent '_SW/' rainevent '/streamflowarea'];
            fid=fopen(ffnm,'rb','ieee-le');
            if fid == -1
                continue
            else
                fclose(fid);

                ffnm = ['/shared/dondo/home/ml423/world_mts_runs/Basins_events_chars/'];
                fid = fopen([ffnm 'Basin' num2str(gid) '_event' rainevent '.txt'],'r');
                if fid == -1
                    fclose(fid);
                    continue
                end
                windowtimes = fscanf(fid,'%d\n',[4,1]);
                fclose(fid);

                %close all


                if 1==1

                    %fh = APL_frontcut(bid,eid);

                    tstr = [num2str(dury(1)) '-' num2str(durm(1),'%2.2d') '/' num2str(durd(1),'%2.2d') '-' num2str(durm(30),'%2.2d') '/' num2str(durd(30),'%2.2d')];


                    obstep = 90;

                    tmp = load(['/shared/dondo/home/ml423/world_mts/obs_events/events/obs_' num2str(gid) '_event' rainevent '.mat']);
                    tmp = tmp.obsnow;
                    obsnow = tmp;clear tmp
                    obsnan=find(isnan(obsnow));
                    if isempty(obsnan)
                    else
                        for na = 1:length(obsnan)
                            obsnow(obsnan(na)) = 0.5*(obsnow(obsnan(na)-1)+obsnow(obsnan(na)+1));
                        end
                    end
                    obs2 = interp1(1:31,obsnow,[1:0.33333:31]);
                    obsnow = obs2(1:90);



                    tmp = load(['/shared/dondo/home/ml423/world_mts/Precip/Event/ERA5_Basin' num2str(gid) '_' rainevent '.mat']);
                    s4dbkc = tmp.rain;

                    MCp = zeros(c,r);
                    pts = [];
                    dbkcr = 0;rmax = [];
                    for i = 1:720
                        MCp = MCp + s4dbkc{i};
                        pts(i) = mean(s4dbkc{i}(intt));
                        dbkcr = dbkcr+pts(i);
                        rmax(i) = max(s4dbkc{i}(intt));
                    end
                    MCp(ind) = nan;
                    ori_mean = nanmean(MCp(:));
                    ama = jet(40);ama(1,:) = [1 1 1];
                    figure
                    imagesc(MCp')
                    colorbar
                    colormap(ama)
                    caxis([0.2*min(MCp(:)) 1.1*max(MCp(:))])
                    % hold on
                    % scatter(streamx,streamy,'k.')
                    % scatter(rgeve(:,2),rgeve(:,1),rgeve(:,3)*20,'kx')
                    % hold off




                    % for the paper

                    launchp = find(rmax>0.5);% when rainfall >0.5mm/hr launch particles
                    ptmaxthre = 0.5*max(rmax);% used as a threshold to separate events
                    la2 = find(launchp<(windowtimes(4)*8));
                    la1 = find(launchp<(windowtimes(3)*8));
                    accu0 = 0;
                    if isempty(la1)
                        continue
                    end
                    for tmpla = launchp(la1(end)):-1:1
                        if ismember(tmpla,launchp)
                            accu0 = 0;
                            rainbegintime = tmpla;
                        else
                            accu0 = accu0+1;
                        end

                        if area>1000
                            noraininterval = 12;% for big basins, it is likely there is always rain in the basin, therefore using a small norain limit to distinguish different rain events.
                        else
                            noraininterval = 24;
                        end
                        if accu0>noraininterval
                            % consecutive x hours of basin max rain less than 0.5mm/hr,
                            % this should be counted as a different rain cell, should not
                            % be tracked!
                            % one exception is if another major rain is 48-24 ahead, should be
                            % considered as well
                            rsmallest = max([(tmpla-noraininterval),1]);

                            reconsider_period = rsmallest:tmpla;
                            %reconsider_period = (tmpla-noraininterval):tmpla;
                            secondary_rain = find(rmax(reconsider_period)>0.5);

                            if max(rmax(reconsider_period))>ptmaxthre
                                % adjust the rainbegintime here
                                rainbegintime = reconsider_period(secondary_rain(1));
                            end
                            smalltmpla = max([floor(tmpla/8),2]);
                            if obsnow(smalltmpla)>obsnow(smalltmpla-1)
                                % meaning obsflow is still rising.
                                if ~isempty(secondary_rain)
                                    rainbegintime = reconsider_period(secondary_rain(1));
                                end
                            end

                            break
                        end
                    end

                    la0 = find(launchp == rainbegintime);

                    refinedlaunchp = launchp(la0:la2(end));
                    stat_left_bd_hydro = windowtimes(1);% systematic defined
                    stat_left_bd_rain = floor(rainbegintime/8);%target rainevent adjusted
                    stat_left_bd = max([stat_left_bd_hydro,stat_left_bd_rain]);

                    upcap = min([90 (windowtimes(4)+6)]);

                    obsref = obsnow(stat_left_bd);obstmp = obsref;
                    [maxq,maxloc] = max(obsnow(stat_left_bd:windowtimes(4)));
                    maxloc = maxloc-1+stat_left_bd;
                    if isempty(maxloc)
                        event_quality_window = [event_quality_window;[bid,eid]];
                        event_quality_windowall = [event_quality_windowall;[bid,eid]];
                        continue
                    end

                    for ileft = stat_left_bd:-1:1
                        if obsnow(ileft)<=1.1*obstmp% allows alittle bit fluctuations
                            obstmp = obsnow(ileft);
                            dips = ileft;
                        else
                            break
                        end
                    end
                    % if (stat_left_bd-dips)>=1*(maxloc-stat_left_bd)&&obsref>(0.67*obstmp+0.33*maxq)
                    %     if rain_prct(bid,eid)>0.5
                    %         event_quality_window = [event_quality_window;[bid,eid]];
                    %     end
                    %     event_quality_windowall = [event_quality_windowall;[bid,eid]];
                    % end
                    %if (stat_left_bd-dips)>=15
                    tmplo = (2*stat_left_bd-maxloc);
                    if tmplo<=1
                       tmplo = 1;
                    end
                    if sum(obsnow(tmplo:stat_left_bd))>0.33*sum(obsnow(stat_left_bd:maxloc))&&obstmp>0.33*maxq
                         event_flash = [event_flash;[bid,eid]];
                    end

                    if upcap<=stat_left_bd
                        continue
                    end

                    shj = [1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7];
                    tww = 3;
                    KGEwxy = [];
                    no3sum = [];
                    ftot = [];
                    %{
                    f = figure('visible','on');
                    f.Position = [5 5 700 350];

                    pos1 = [0.13 0.8 0.67 0.15];
                    %subplot('Position',pos1)
                    pos2 = [0.15 0.32 0.78 0.42];

                    subplot('Position',pos2)
                    % hjh = plot(2:3:obstep*3,obsdatav1(:,6),'k',1:288,no3,'g',1:288,no2,'r',1:288,no5,'m');
                    % set(hjh,'LineWidth',2.5)
                    % legend('Obs.','IMERG_D','STIV_D','STIV_{DBK}','Location','Northeast')
                    %yyaxis left
                    xj = 2:8:obstep*8;
                    xj = xj';
                    bd1 = 1.05*obsnow;
                    bd2 = 0.95*obsnow;
                    sh = scatter(2:8:obstep*8,obsnow,45,'g','filled');
                    sh.MarkerEdgeColor = [0,1,0];
                    hold on
                    ylim([0 ceil(max(obsnow)/20)*20])
                    title('W00')

                    hjb = plot([1+8*stat_left_bd 1+8*stat_left_bd],[0 1000000],'k-.');
                    hold on
                    hjb2 = plot([1+8*(windowtimes(4)+6) 1+8*(windowtimes(4)+6)],[0 1000000],'r-.');

                    % ymax is based on max of y but rounded to
                    % nearest interval defined by max of y

                    ylabel('Streamflow (m^3/s)')

                    box on


                    xlabel([tstr ' ' num2str(bid) ' ' num2str(eid)])
                    %set(gca,'XTick',1:3:288,'XTickLabel',[])
                    xlim([0 719])
                    set(gca,'XTick',0:144:720,'XTickLabel',xtlabels(1:6:30),'FontSize',12,'LineWidth',2.5,'FontWeight','bold')


                    hax = gca;
                    hax.XAxis.MinorTickValues = linspace(0,720,31);
                    hax.XMinorTick = 'on';


                    set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',2.5,'FontWeight','bold');
                    Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
                    %[Left Bottom Right Top] spacing
                    NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3)-0.01 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
                    set(gca, 'Position', NewPos);
                    %}

                end
            end
        end
    end
end
%
event_same = [];
% repititve events
for bid = [134:465]%tmpjjj = 1:length(bidt)% [1:133 466:522]%1:length(reprob)
    disp(bid)

    gid = WORLD_events(bid,1);

    gloc = find(Basins_sum(:,1)==gid);


    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_ind.mat'],'ind');
    ind = tmp.ind;
    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_intt.mat'],'intt');
    intt = tmp.intt;


    ixx = find(~isnan(WORLD_events(bid,:)));
    etot = length(ixx)-1;

    % check if this basin has rgs

    if etot == 0
        continue
    end


    for eid = 1:etot-1%eidt(tmpjjj)%1:etot
        %reprob(prob,2)%1:etot%[1:15]%[1 2 3 4 6 8]%:5%1:etot%[5 6 7]%[1 2 3 4 5 6 7 8 9 10 11]
        close all
        %rglists = find(rgr(:,1)==bid&rgr(:,2)==eid);
        %rgeve = rgr(rglists,5:7);

        if isnan(event_quality(bid,eid+1))
            disp('Event quality is bad, peak too early or too late')
            continue
        else
            rainevent = num2str(WORLD_events(bid,eid+1));
            yth = floor(WORLD_events(bid,eid+1)/10000);
            dth = mod(WORLD_events(bid,eid+1),100);
            mth = floor((WORLD_events(bid,eid+1)-yth*10000)/100);

            cud = datetime(yth,mth,dth);

            rainevent_next = num2str(WORLD_events(bid,eid+2));
            yth2 = floor(WORLD_events(bid,eid+2)/10000);
            dth2 = mod(WORLD_events(bid,eid+2),100);
            mth2 = floor((WORLD_events(bid,eid+2)-yth2*10000)/100);

            cud2 = datetime(yth2,mth2,dth2);

            if (cud2-cud)<days(8)&&(cud2-cud)>days(0)
                event_same = [event_same;[bid,eid+1]];

            end
        end
    end
end

for i = 1:length(event_flash)
    eqcontrol(event_flash(i,1),event_flash(i,2)) = nan;

end
for i = 1:length(event_same)
    eqcontrol(event_same(i,1),event_same(i,2)) = nan;

end


%% write out IRC window boundary, and event quality
%write IRC window points.

tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/world_selected_gauge_events_allinfo_v3.mat'],'gauge_events_first_last');
Basins_sum = tmp.gauge_events_first_last;clear tmp;
tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/WORLD_events_v3.mat'],'WORLD_events');
WORLD_events = tmp.WORLD_events;clear tmp;

%event_quality = nan(size(WORLD_events,1),size(WORLD_events,2));


for bid = 466:length(WORLD_events)
    close all
    gid = WORLD_events(bid,1);
    % gloc = find(Basins_sum(:,1)==gid);
    % area = Basins_sum(gloc,5);
    event_quality(bid,1:300) = nan;
    ixx = find(~isnan(WORLD_events(bid,:)));
    etot = length(ixx)-1;
    event_quality(bid,1) = gid;
    for eid = 1:etot%[5 6 7]%[1 2 3 4 5 6 7 8 9 10 11]
        close all
        riseqc = 1;
        rainevent = num2str(WORLD_events(bid,eid+1));
        
        tmp = load(['/shared/dondo/home/ml423/world_mts/obs_events/events/obs_' num2str(gid) '_event' rainevent '.mat']);
        tmp = tmp.obsnow;
        obsnow = tmp;clear tmp
         
        
        if isempty(obsnow)
            continue
        end
        
        [~,maxloc] = max(obsnow);
        if maxloc<5||maxloc>26% if the peak too early or too late within the 30 day window, discard it because not useful for IRC
            continue
        end
        
        p66 = prctile([obsnow(1:maxloc)],66.6);
        if p66>0.5*(obsnow(1)+obsnow(maxloc))
           continue 
        end
        
        
%         method 1 below for v2 quality, not good, use the above
%         riseaccu = 0;
%         for ik = maxloc:-1:2
%             if obsnow(ik)>(obsnow(ik-1)+0.1)% to prevent flat pre rising, therefore 0.1
%                 riseaccu = riseaccu+1;
%             else
%                 riseaccu = 0;
%             end
%             if riseaccu>=10
%                riseqc = 0;
%                break
%             end
%         end
%         
%         if riseqc == 0
%            continue 
%         end
        

        obsnan=find(isnan(obsnow));
        if isempty(obsnan)
        else
            for na = 1:length(obsnan)
                obsnow(obsnan(na)) = 0.5*(obsnow(obsnan(na)-1)+obsnow(obsnan(na)+1));
            end
        end
        
        obs2 = interp1(1:31,obsnow,[1:0.33333:31]);
        obsnow = obs2(1:90);
        
        [~,p3] = max(obsnow);
        for k = (p3-1):-1:1
            if obsnow(k)>obsnow(k+1)
               p2 = k+1;
               break
            end
        end
        p1 = p2-6;
        if p1<1
           p1=1; 
        end
        
        ptmp = min([p3+6,90]);
        kre = (obsnow(p3)-obsnow(ptmp))/(ptmp-p3);
        basef = mean(obsnow(82:90));
        
        for p4 = (p3+1):90
           extendvalue = obsnow(p3)-(p4-p3)*kre;
           if extendvalue<=basef
              break 
           end
        end
        
        figure
        plot(obsnow)
        hold on
        scatter([p1 p2 p3 p4],[obsnow(p1) obsnow(p2) obsnow(p3) obsnow(p4)],40,'ro','filled')
        hold off
        title(['Basin ' num2str(gid) ' event ' rainevent ])
        set(gca,'XTick',0:10:90,'XTickLabel',0:10:90,'FontSize',12,'LineWidth',2.5,'FontWeight','bold')
        saveas(gcf,['/shared/dondo/home/ml423/world_mts_runs/Basins_events_chars/Basin' num2str(gid) '_event' rainevent '_fig'],'jpg')
        
        windows = [p1;p2;p3;p4];
        ffnm = ['/shared/dondo/home/ml423/world_mts_runs/Basins_events_chars/'];
        fid = fopen([ffnm 'Basin' num2str(gid) '_event' rainevent '.txt'],'w');
        fprintf(fid,'%d\n',windows);
        fclose(fid);
        
        if p3<15||p3>75
            
        else
            event_quality(bid,eid+1) = 1;
        end
        
    end
end

%
% yth = floor(WORLD_events(bid,eid+1)/10000);
% dth = mod(WORLD_events(bid,eid+1),100);
% mth = floor((WORLD_events(bid,eid+1)-yth*10000)/100);

%             ffnm = ['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LA/Basins_events_chars/'];
%             fid = fopen([ffnm 'Basin' num2str(bd,'%2.2d') '_event' num2str(event_ind,'%2.2d') '.txt'],'w');
%             fprintf(fid,'%d\n',windows);
%             fclose(fid);

%% calculate/plot IRC results, with statistics 
clear tmpp
clear tmp
clear predatatmp
clear data3D
addpath('/shared/dondo/home/ml423/HMIOP/NoDAre/New era/')
addpath('/shared/dondo/home/ml423/DA/Critical matrices/')
addpath('/shared/dondo/home/ml423/DA/slanCM/')

% KGE_ori = nan(size(APL_events,1),size(APL_events,2));
% KGE_new = nan(size(WORLD_events,1),size(WORLD_events,2));
% NSE_ori = nan(size(WORLD_events,1),size(WORLD_events,2));
% NSE_new = nan(size(WORLD_events,1),size(WORLD_events,2));
% EPV_ori = nan(size(WORLD_events,1),size(WORLD_events,2));
% EPV_new = nan(size(WORLD_events,1),size(WORLD_events,2));
% EPT_ori = nan(size(WORLD_events,1),size(WORLD_events,2));
% EPT_new = nan(size(WORLD_events,1),size(WORLD_events,2));
% EV_ori = nan(size(WORLD_events,1),size(WORLD_events,2));
% EV_new = nan(size(WORLD_events,1),size(WORLD_events,2));
%best_perform = nan(size(WORLD_events,1),size(WORLD_events,2));
tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/world_selected_gauge_events_allinfo_v3.mat'],'gauge_events_first_last');
Basins_sum = tmp.gauge_events_first_last;clear tmp;
tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/WORLD_events_v3.mat'],'WORLD_events');
WORLD_events = tmp.WORLD_events;clear tmp;

flowtype = {'streamflowarea','interflowarea','overlandflowarea','baseflowarea'};
ftp = 1;

% kge_opt = [];kge_ori = [];
% rain00_sum = cell(522,149);
% rain34_sum = cell(522,149);
% eval_ind = {};%cell(522,149);
% simu_flood = {};%cell(522,149);
% cvflood = nan(522,149);maxflood = nan(522,149);meanflood = nan(522,149);
% 
% EVtot = nan(522,149);
% kge_optm(1:522,1:149) = nan;
% kge_orim(1:522,1:149) = nan;
% rain_orim(1:522,1:149) = nan;
% rain_optm(1:522,1:149) = nan;
% rainsd_orim(1:522,1:149) = nan;
% rainsd_optm(1:522,1:149) = nan;

%
%listrg = unique(rgs_loc(:,2));
% bidt = bidbig(:,1);
% eidt = bidbig(:,3);


for bid = [1:133 466:522]%1:length(reprob)
    %bid = reprob(prob,1);%[74 59 517 60 25 121 69 512 502 107 519 82 521 70 56]
    %bid = bidt(tmpjjj);
    %[1:133 466:522]
    %tmpii = 1:length(listrg) %[50 79 91 93 519 92 94 104 89 110 472 481 484 490 491 494]%[1:133 466:522]%[268:300]%37:100%26:100%[ 355 359 86 22 360] %98:190
    %bid = listrg(tmpii);
    %close all



    gid = WORLD_events(bid,1);

    gloc = find(Basins_sum(:,1)==gid);

    r = Basins_sum(gloc,6);
    c = Basins_sum(gloc,7);
    rg = Basins_sum(gloc,8);
    cg = Basins_sum(gloc,9);
    area = Basins_sum(gloc,5);

    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_ind.mat'],'ind');
    ind = tmp.ind;
    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_intt.mat'],'intt');
    intt = tmp.intt;


    ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_fixed_900m/dem.bin'];
    fid=fopen(ffnm,'rb','ieee-le');
    tmpdata=fread(fid,inf,'single');
    fclose(fid);
    dem = reshape(tmpdata,c,r);
    dem(ind) = nan;
    demf = dem';
    colo = slanCM('terrain',80);
    colo(1,:) = [0.8 0.8 0.8];
    figure
    imagesc(demf)
    colormap(colo)

    ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_fixed_900m/countfile.out'];
    fid=fopen(ffnm);
    facc=fscanf(fid,'%d',[c,r]);
    fclose(fid);
    facc(ind) = nan;
    [streamx,streamy] = find(facc>5);


    ICn = 4000;
    ICn = 4140;
    ICC = 0;% IF ICC is part of the algo

    fileind = -100000*ICn;
    KGEsu = [];
    evesu = [];
    budgetalldbkc = [];
    budgetallirc = [];

    ixx = find(~isnan(WORLD_events(bid,:)));
    etot = length(ixx)-1;

    % check if this basin has rgs

    if etot == 0
        continue
    end

    % bidrg = find(cell2mat(rgs_loc_sum(:,2))==bid);
    % rgr = [];
    % if ~isempty(bidrg)
    %     % this basin has >=1 raingauges
    %     for irg = 1:length(bidrg)
    %         for eid = 1:etot
    %             rgrecord = rgs_loc_sum(bidrg(irg),eid+6);
    %             if isempty(rgrecord{1})
    %                 %this rg has no recods for this event
    %                 rgr = [rgr;[bid,eid,cell2mat(rgs_loc_sum(bidrg(irg),1:4)),0.2]];
    % 
    %             else
    %                 %this rg has records for this event
    %                 rgr = [rgr;[bid,eid,cell2mat(rgs_loc_sum(bidrg(irg),1:4)),1]];
    %             end
    %         end
    %     end
    % end

    for eid = 1:etot%eidt(tmpjjj)%1:etot
        %reprob(prob,2)%1:etot%[1:15]%[1 2 3 4 6 8]%:5%1:etot%[5 6 7]%[1 2 3 4 5 6 7 8 9 10 11]
        close all
        %rglists = find(rgr(:,1)==bid&rgr(:,2)==eid);
        %rgeve = rgr(rglists,5:7);

        if isnan(event_quality(bid,eid+1))
            disp('Event quality is bad, peak too early or too late')
            continue
        else
            rainevent = num2str(WORLD_events(bid,eid+1));
            yth = floor(WORLD_events(bid,eid+1)/10000);
            dth = mod(WORLD_events(bid,eid+1),100);
            mth = floor((WORLD_events(bid,eid+1)-yth*10000)/100);

            thresdate = datetime(yth,mth,dth);

            durati = (thresdate-days(15)):days(1):(thresdate+days(15));
            [dury,durm,durd] = ymd(durati);
            tst1 = num2str(durm(1)*100+durd(1),'%4.4d');
            tst2 = num2str(durm(30)*100+durd(30),'%4.4d');

            xtlabels ={};
            for tl = 1:30
                xtlabels{tl} =  [ num2str(durm(tl),'%2.2d') '/' num2str(durd(tl),'%2.2d') ];
            end



            ffnm=['/shared/dondo/home/ml423/world_mts_runs/IRC1/Basin' num2str(gid) 'outputs/refxtmp' num2str(str2num(rainevent)*100+fileind+0) 'fullresults/' rainevent '_SW/' rainevent '/streamflowarea'];
            fid=fopen(ffnm,'rb','ieee-le');
            if fid == -1
                continue
            else
                fclose(fid);

                ffnm = ['/shared/dondo/home/ml423/world_mts_runs/Basins_events_chars/'];
                fid = fopen([ffnm 'Basin' num2str(gid) '_event' rainevent '.txt'],'r');
                if fid == -1
                    fclose(fid);
                    continue
                end
                windowtimes = fscanf(fid,'%d\n',[4,1]);
                fclose(fid);

                %close all


                if 1==1

                    %fh = APL_frontcut(bid,eid);

                    tstr = [num2str(dury(1)) '-' num2str(durm(1),'%2.2d') '/' num2str(durd(1),'%2.2d') '-' num2str(durm(30),'%2.2d') '/' num2str(durd(30),'%2.2d')];


                    obstep = 90;

                    tmp = load(['/shared/dondo/home/ml423/world_mts/obs_events/events/obs_' num2str(gid) '_event' rainevent '.mat']);
                    tmp = tmp.obsnow;
                    obsnow = tmp;clear tmp
                    obsnan=find(isnan(obsnow));
                    if isempty(obsnan)
                    else
                        for na = 1:length(obsnan)
                            obsnow(obsnan(na)) = 0.5*(obsnow(obsnan(na)-1)+obsnow(obsnan(na)+1));
                        end
                    end
                    obs2 = interp1(1:31,obsnow,[1:0.33333:31]);
                    obsnow = obs2(1:90);



                    tmp = load(['/shared/dondo/home/ml423/world_mts/Precip/Event/ERA5_Basin' num2str(gid) '_' rainevent '.mat']);
                    s4dbkc = tmp.rain;

                    MCp = zeros(c,r);
                    pts = [];
                    dbkcr = 0;rmax = [];
                    for i = 1:720
                        MCp = MCp + s4dbkc{i};
                        pts(i) = mean(s4dbkc{i}(intt));
                        dbkcr = dbkcr+pts(i);
                        rmax(i) = max(s4dbkc{i}(intt));
                    end
                    MCp(ind) = nan;
                    ori_mean = nanmean(MCp(:));
                    ama = jet(40);ama(1,:) = [1 1 1];
                    % figure
                    % imagesc(MCp')
                    % colorbar
                    % colormap(ama)
                    % caxis([0.2*min(MCp(:)) 1.1*max(MCp(:))])
                    % hold on
                    % scatter(streamx,streamy,'k.')
                    % scatter(rgeve(:,2),rgeve(:,1),rgeve(:,3)*20,'kx')
                    % hold off




                    % for Nature Science paper, consider the major events only

                    launchp = find(rmax>0.5);% when rainfall >1mm/hr launch particles
                    ptmaxthre = 0.5*max(rmax);% used as a threshold to separate events
                    la2 = find(launchp<(windowtimes(4)*8));
                    la1 = find(launchp<(windowtimes(3)*8));
                    accu0 = 0;
                    if isempty(la1)
                        continue
                    end
                    for tmpla = launchp(la1(end)):-1:1
                        if ismember(tmpla,launchp)
                            accu0 = 0;
                            rainbegintime = tmpla;
                        else
                            accu0 = accu0+1;
                        end

                        if area>1000
                            noraininterval = 12;% for big basins, it is likely there is always rain in the basin, therefore using a small norain limit to distinguish different rain events.
                        else
                            noraininterval = 24;
                        end
                        if accu0>noraininterval
                            % consecutive x hours of basin max rain less than 0.5mm/hr,
                            % this should be counted as a different rain cell, should not
                            % be tracked!
                            % one exception is if another major rain is 48-24 ahead, should be
                            % considered as well
                            rsmallest = max([(tmpla-noraininterval),1]);

                            reconsider_period = rsmallest:tmpla;
                            %reconsider_period = (tmpla-noraininterval):tmpla;
                            secondary_rain = find(rmax(reconsider_period)>0.5);

                            if max(rmax(reconsider_period))>ptmaxthre
                                % adjust the rainbegintime here
                                rainbegintime = reconsider_period(secondary_rain(1));
                            end
                            smalltmpla = max([floor(tmpla/8),2]);
                            if obsnow(smalltmpla)>obsnow(smalltmpla-1)
                                % meaning obsflow is still rising.
                                if ~isempty(secondary_rain)
                                    rainbegintime = reconsider_period(secondary_rain(1));
                                end
                            end

                            break
                        end
                    end

                    la0 = find(launchp == rainbegintime);

                    refinedlaunchp = launchp(la0:la2(end));
                    stat_left_bd_hydro = windowtimes(1);% systematic defined
                    stat_left_bd_rain = floor(rainbegintime/8);%target rainevent adjusted
                    stat_left_bd = max([stat_left_bd_hydro,stat_left_bd_rain]);

                    upcap = min([90 (windowtimes(4)+6)]);

                    if upcap<=stat_left_bd
                        continue
                    end

                    shj = [1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7];
                    tww = 3;
                    KGEwxy = [];
                    no3sum = [];
                    ftot = [];
                    for p = 0:15

                        %itranum = 0;

                        if p==0
                            filenam = str2num(rainevent)*100+fileind;
                            wnum = mod(filenam,10);
                            itranum = 0;
                        else
                            wnum = mod(filenam,10);
                            itranum = shj(p);

                            if mod(p,tww)==1%window 2
                                zpo = windowtimes(1);
                            elseif mod(p,tww)==2%window 3
                                zpo = windowtimes(2);
                            elseif mod(p,tww)==0%window 4
                                zpo = windowtimes(3);
                            end

                        end


                        ffnmx=['/shared/dondo/home/ml423/world_mts_runs/IRC1/Basin' num2str(gid) 'outputs/refxtmp' num2str(filenam) 'fullresults/' rainevent '_SW/' rainevent '/' flowtype{ftp}];
                        fid=fopen(ffnmx,'rb','ieee-le');

                        if fid==-1
                            %    realIRCstat = [realIRCstat;[raininfob1(rit,1),nan,nan,nan,nan,nan]];
                            %disp('error')
                            %ffnm

                            % meaning the files are not there, so move
                            % on to next iteration. THis can be because
                            % the user delete intermediate files to
                            % save space.
                            if p==0
                                filenam = filenam+2;
                            else
                                if mod(p,tww)==1
                                    filenam = filenam+1;
                                elseif mod(p,tww)==2
                                    filenam = filenam+1;
                                elseif mod(p,tww)==0
                                    filenam = filenam+8;
                                end
                            end
                            no3sum = [no3sum;nan(1,720)];

                            if p<9
                                continue
                            else
                                break
                            end
                        end
                        tmpdata=fread(fid,inf,'single');
                        fclose(fid);
                        NODAc2o5=reshape(tmpdata,c,r,720);
                        flownoda5 = [];
                        for k = 1:720
                            flownoda5(k) = abs(NODAc2o5(cg,rg,k));
                        end

                        no3 = flownoda5;
                        no3sum = [no3sum;no3];

                        if p==0
                            nostivdbkc = no3;
                        end


                        ffnm=['/shared/dondo/home/ml423/world_mts_runs/IRC1/Basin' num2str(gid) 'outputs/refxtmp' num2str(filenam) 'fullresults/' rainevent '_SW/' rainevent '/precipitation'];
                        fid=fopen(ffnm,'rb','ieee-le');
                        tmpdata=fread(fid,inf,'single');
                        fclose(fid);
                        NODAc2o5=reshape(tmpdata,c,r,720);
                        flownoda5 = [];
                        meanr250 = 0;
                        tmpsum = zeros(c,r);
                        for k = 1:720
                            tmp = NODAc2o5(:,:,k);
                            tmpsum = tmpsum+tmp;
                            flownoda5(k) = mean(tmp(intt));
                            meanr250 = meanr250+flownoda5(k);
                        end
                        no2 = flownoda5;%GREEN
                        rainsd = std(tmpsum(intt));

                        cal_window = (stat_left_bd-1)*8+1:(upcap-1)*8+1;

                        if p==9

                            % figure
                            % imagesc(tmpsum')
                            % colorbar
                            % colormap(ama)
                            % hold on
                            % scatter(rgeve(:,2),rgeve(:,1),rgeve(:,3)*30,'kx')
                            % hold off
                            % caxis([0.2*min(MCp(:)) 1.1*max(MCp(:))])
                            % figure
                            % heatmap(tmpsum')
                            % colorbar
                            % colormap(ama)
                            
                            %caxis([0.2*min(MCp(:)) 1.1*max(MCp(:))])


                        end

                        if p==0
                            tstmp = find(flownoda5>0.1);% only consider rain>0.1mm/h in the original field
                            % ircicc rain>0.1 hours are a little bit
                            % more than original rain>0.1 hours have to
                            % be consistent for the number of hours for
                            % error modeling
                            rts1 = intersect(tstmp,cal_window);
                            rain00_sum{bid,eid} = NODAc2o5(:,:,rts1);
                            eval_ind{bid,eid} = rts1;
                            simu_flood{bid,eid} = no3(rts1);
                            if isempty(rts1)
                                cvflood(bid,eid) = nan;
                                maxflood(bid,eid) = nan;
                                meanflood(bid,eid) = nan;
                            else
                                cvflood(bid,eid) = mean(no3(rts1))/std(no3(rts1));
                                maxflood(bid,eid) = max(no3(rts1));
                                meanflood(bid,eid) = mean(no3(rts1));
                            end
                        end
                        %
                        if p<8
                            rain34_sum{bid,eid} = [];
                        end
                        %disp(p)
                        if p==8
                            rain33 = NODAc2o5(:,:,rts1);
                        end
                        if p==9
                            rain34_sum{bid,eid} = NODAc2o5(:,:,rts1);
                        end


                        % mass balance check fo W00 and W34
                        %{
                        filenamp = str2num(rainevent)*100+fileind+2;% the very first P
                        cc = ls(['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Basin' num2str(bid,'%2.2d') 'outputs/refxtmp' num2str(filenamp) '/xtmpbb1*']);
                        cc2 = strsplit(cc);
                        tmpC  = struct2cell(load([cc2{end-1}]));
                        P00 = tmpC{1};
                        filenamp = str2num(rainevent)*100+fileind+24;% the very first P
                        cc = ls(['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Basin' num2str(bid,'%2.2d') 'outputs/refxtmp' num2str(filenamp) '/xtmpbb1*']);
                        cc2 = strsplit(cc);
                        tmpC  = struct2cell(load([cc2{1}]));
                        P34 = tmpC{1};
                        P00tot = 0;
                        P34tot = 0;
                        for k = 1:288
                            tmp = P00{k};
                            P00tot = P00tot+mean(tmp(intt));
                            tmp = P34{k};
                            P34tot = P34tot+mean(tmp(intt));
                        end
                        % mass balance quick calculations
                        runoff250 = sum(no3)*300/length(intt)/250/250*1000;
                        
                        
                        
                        fixedstr = ['/shared/dondo/home/ml423/HMIOP/NoDAre/LB/APL_Basin' num2str(bid,'%2.2d') '_fixed_250m/'];
                        datastr = ['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Basin' num2str(bid,'%2.2d') 'outputs/refxtmp' num2str(filenam) 'fullresults/Basin' num2str(bid,'%2.2d') '_output_250m5min_' rainevent '_SW/' rainevent '00/'];
                        
                        lastnote = 'c';
                        if p==0
                            laststr = ['/shared/dondo/home/ml423/HMIOP/NoDAre/LB/Basin' num2str(bid,'%2.2d') '_output_250m5min_' rainevent '_' lastnote '/'];
                            [smdelta250,gddelta250] = waterbalance(0,fixedstr,c,r,288,288,laststr,datastr,'',intt,250);
                        elseif p==9
                            if ICC==1% if ICC is part of the algo
                                laststr = ['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/ICall/' num2str(ICn) '/Basin' num2str(bid,'%2.2d') '_output_250m5min_' rainevent lastnote '/'];
                            else
                                laststr = ['/shared/dondo/home/ml423/HMIOP/NoDAre/LB/Basin' num2str(bid,'%2.2d') '_output_250m5min_' rainevent '_' lastnote '/'];
                            end
                            [smdelta250,gddelta250] = waterbalance(0,fixedstr,c,r,288,288,laststr,datastr,'',intt,250);
                        end
                        
                        if p==0
                            rainDBKC250 = dbkcr;% in mm
                            riverDBKC250 = runoff250;% in mm
                            deltasDBKC250 = smdelta250/250/250/length(intt)*1000; %in mm
                            gdDBKC250 = gddelta250/250/250/length(intt)*1000; %in mm
                            budgetDBKC250 = rainDBKC250-riverDBKC250-(deltasDBKC250+gdDBKC250);
                            wbterms_DBKC = [rainDBKC250,riverDBKC250,deltasDBKC250,gdDBKC250,budgetDBKC250];
                            budgetalldbkc = [budgetalldbkc;[str2num(rainevent),wbterms_DBKC]];
                        elseif p==9
                            rainIRC250 = P34tot;% in mm
                            riverIRC250 = runoff250;% in mm
                            deltasIRC250 = smdelta250/250/250/length(intt)*1000; %in mm
                            gdIRC250 = gddelta250/250/250/length(intt)*1000; %in mm
                            budgetIRC250 = rainIRC250-riverIRC250-(deltasIRC250+gdIRC250);
                            wbterms_IRC = [rainIRC250,riverIRC250,deltasIRC250,gdIRC250,budgetIRC250];
                            budgetallirc = [budgetallirc;[str2num(rainevent),wbterms_IRC]];
                        end
                        
                        
                        %}


                        %
                        simu = no3(1:8:end);


                        simustat = simu(stat_left_bd:upcap);
                        obsnowstat = obsnow(stat_left_bd:upcap);
                        tmp = 0;
                        tmp1 = 0;
                        tmp2 = mean(obsnowstat);
                        for i = 1:length(obsnowstat)
                            tmp = tmp + (simustat(i)-obsnowstat(i))^2;
                            tmp1 = tmp1 + (tmp2-obsnowstat(i))^2;
                        end
                        NSE = 1-tmp./tmp1;

                        rr = corrcoef(simustat',obsnowstat);
                        rstar = rr(1,2);
                        obsstd = std(obsnowstat);
                        simustd = std(simustat);
                        obsmean = mean(obsnowstat);
                        simumean = mean(simustat);
                        KGE = 1-sqrt((rstar-1)^2+(simustd/obsstd-1)^2+(simumean/obsmean-1)^2);
                        KGEwxy = [KGEwxy;KGE];
                        [vsimu,tsimu]=max(simustat);
                        [vobs,tobs]=max(obsnowstat);
                        EPV = vsimu-vobs;
                        EPT = (tsimu-tobs)*1;
                        simvol = 0;
                        obsvol = 0;
                        for vlo=2:length(obsnowstat)-1
                            obsvol = obsvol+obsnowstat(vlo)*1;
                            simvol = simvol+simustat(vlo)*1;
                        end
                        obsvol = obsvol+0.5*(obsnowstat(1)+obsnowstat(end))*1;
                        simvol = simvol+0.5*(simustat(1)+simustat(end))*1;

                        EV = (simvol-obsvol)./obsvol;
                        % KGE is period defined, while ther other stats
                        % use the whole series
                        if p==0
                            EVtot(bid,eid) = EV;
                        end

                        if p == 0
                            KGEori = KGE;
                            tmp = get_statistics(WORLD_events(bid,eid+1),no3,obsnow,8,1);
                            % KGE_ori(bid,eid+1) = tmp(3);
                            % NSE_ori(bid,eid+1) = tmp(2);
                            % EPV_ori(bid,eid+1) = tmp(4);
                            % EPT_ori(bid,eid+1) = tmp(5);
                            % EV_ori(bid,eid+1) = tmp(6);
                            KGE_base = KGE;
                            p_base = p;
                            NSE_base = tmp(2);
                            EPV_base = tmp(4);
                            EPT_base = tmp(5);
                            EV_base = tmp(6);
                            %best_perform(bid,eid+1) = 0;
                            fopt = filenam;
                            fopt_last = filenam;
                            rain_ori = meanr250;
                            rain_base = meanr250;
                            rainsd_ori = rainsd;
                            rainsd_base = rainsd;
                        end
                        if p==8
                            % W33
                            tmp = get_statistics(WORLD_events(bid,eid+1),no3,obsnow,8,1);
                            KGE_base = KGE;
                            rain_base = meanr250;%20250606 added
                            p_base = p;
                            rain_opt = meanr250
                            rainsd_opt = rainsd;
                            NSE_base = tmp(2);
                            EPV_base = tmp(4);
                            EPT_base = tmp(5);
                            EV_base = tmp(6);
                            KGE33 = KGE;
                        end


                        if p >= 8
                            tmp = get_statistics(WORLD_events(bid,eid+1),no3,obsnow,8,1);
                            if KGE>KGE_base

                                KGEopt = KGE;
                                opt_p = p;
                                p_base = p;
                                KGE_base = KGE;
                                NSEopt = tmp(2);
                                EPVopt = tmp(4);
                                EPTopt = tmp(5);
                                EVopt = tmp(6);
                                %best_perform(bid,eid+1) = itranum*10+wnum;
                                fopt_last = fopt;
                                fopt = filenam;
                                rain_opt = meanr250;
                                rain_base = meanr250;
                                rainsd_opt = rainsd;
                                rainsd_base = rainsd;
                            else
                                KGEopt = KGE_base;
                                opt_p = p_base;
                                NSEopt = NSE_base;
                                EPVopt = EPV_base;
                                EPTopt = EPT_base;
                                EVopt = EV_base;
                                rain_opt = rain_base;
                                rainsd_opt = rainsd_base;
                            end
                            if KGE_base==KGEori&&p>=9
                                %basically w33 is missing,
                                KGEopt = KGE;
                                rain_opt = meanr250;
                            end

                        end
                        %}

                        %{
                        % disp(['NSE is ' num2str(NSE)])
                        % disp(['KGE is ' num2str(KGE)])
                        % disp(['EPV is ' num2str(EPV)])
                        % disp(['EPT is ' num2str(EPT)])
                        %

                        if p==0||p==9
                            %                             figure
                            %                             t = tiledlayout(1,1,'Padding','none');
                            %                             t.Units = 'inches';
                            %                             t.OuterPosition = [0.25 0.25 4.8 3.2];
                            %                             nexttile;

                            f = figure('visible','on');
                            f.Position = [5 5 700 350];

                            pos1 = [0.13 0.8 0.67 0.15];
                            %subplot('Position',pos1)
                            pos2 = [0.15 0.32 0.78 0.42];

                            subplot('Position',pos2)
                            % hjh = plot(2:3:obstep*3,obsdatav1(:,6),'k',1:288,no3,'g',1:288,no2,'r',1:288,no5,'m');
                            % set(hjh,'LineWidth',2.5)
                            % legend('Obs.','IMERG_D','STIV_D','STIV_{DBK}','Location','Northeast')
                            %yyaxis left
                            xj = 2:8:obstep*8;
                            xj = xj';
                            bd1 = 1.05*obsnow;
                            bd2 = 0.95*obsnow;
                            sh = scatter(2:8:obstep*8,obsnow,45,'g','filled');
                            sh.MarkerEdgeColor = [0,1,0];
                            hold on
                            hjh = plot(1:8:720,no3(1:8:720),'m-.');
                            max(no3)
                            hold on
                            hjb = plot([1+8*stat_left_bd 1+8*stat_left_bd],[0 1000000],'k-.');
                            hold on
                            hjb2 = plot([1+8*(windowtimes(4)+6) 1+8*(windowtimes(4)+6)],[0 1000000],'k-.');

                            % ymax is based on max of y but rounded to
                            % nearest interval defined by max of y
                            if p==0
                                ymm = ceil(max(max(obsnow),max(no3))/(max(obsnow)/4))*(max(obsnow)/4);
                            end
                            tf = text(16,0.85*ymm,['' num2str(KGEori,'%2.2f')]);
                            tf.FontSize = 14;
                            tf.FontWeight = 'bold';

                            tf = text(16,0.72*ymm,['' num2str(KGE,'%2.2f')],'Color','m');
                            tf.FontSize = 14;
                            tf.FontWeight = 'bold';
                            set(hjh,'LineWidth',2.5)
                            %             hold on
                            %             patch([xj; flipud(xj)]', [bd1; flipud(bd2)]',  'k', 'FaceAlpha',0.2, 'EdgeColor','none')
                            %             hold off
                            ylim([0 ymm])
                            ylabel('Streamflow (m^3/s)')

                            box on

                            title([num2str(gid) ' IRC W' num2str(itranum) num2str(wnum) ' | ' num2str(area,'%0.0f') ' km^2'])


                            xlabel([tstr ' ' num2str(bid) ' ' num2str(eid)])
                            %set(gca,'XTick',1:3:288,'XTickLabel',[])
                            xlim([0 719])
                            set(gca,'XTick',0:144:720,'XTickLabel',xtlabels(1:6:30),'FontSize',12,'LineWidth',2.5,'FontWeight','bold')


                            hax = gca;
                            hax.XAxis.MinorTickValues = linspace(0,720,31);
                            hax.XMinorTick = 'on';


                            set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',2.5,'FontWeight','bold');
                            Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
                            %[Left Bottom Right Top] spacing
                            NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3)-0.01 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
                            set(gca, 'Position', NewPos);
                            if p==0
                                saveas(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/IRC_images/Bid_' num2str(bid,'%3.3d') '_Eid' num2str(eid,'%3.3d') 'Basin' num2str(gid) '_' rainevent '_' num2str(ICn) '_W' num2str(itranum) num2str(wnum) 'check_prctall.jpg'])
                            end
                            if p==9
                                chanr = (rain_opt-rain_ori)/rain_ori;
                                if chanr>0.5

                                    saveas(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/IRC_images/Bid_' num2str(bid,'%3.3d') '_Eid' num2str(eid,'%3.3d') 'Basin' num2str(gid) '_' rainevent '_' num2str(ICn) '_W' num2str(itranum) num2str(wnum) 'check_prctover05.jpg'])
                            
                                end
                            end

                            if p==0
                                %exportgraphics(gcf,['/shared/dondo/home/ml423/HMIOP/NoDAre/LB/pics/Basin' num2str(bid,'%2.2d') '_event' num2str(eid,'%2.2d') '_XYIRC.jpg'],'Resolution',600)
                                figure
                                bar(no2,1)
                                xlim([0 719])
                                set(gca,'XTick',0:144:720,'XTickLabel',xtlabels(1:6:30),'FontSize',12,'LineWidth',2.5,'FontWeight','bold')
                                hax = gca;
                                hax.XAxis.MinorTickValues = linspace(0,720,31);
                                hax.XMinorTick = 'on';
                                xlabel([tstr ' '])
                            end

                        end

                        %}

                        %
                        if p==0
                            filenam = filenam+2;
                        else
                            if mod(p,tww)==1
                                filenam = filenam+1;
                            elseif mod(p,tww)==2
                                filenam = filenam+1;
                            elseif mod(p,tww)==0
                                filenam = filenam+8;
                            end
                        end
                        ftot = [ftot;filenam];
                        %}

                    end
                    %
                    % done plot optimum Wxy
                    if isempty(rain34_sum{bid,eid})
                        KGEopt = nan;
                    end
                    %kge_opt = [kge_opt;[bid,eid,KGEopt]];
                    %kge_ori = [kge_ori;[bid,eid,KGEori]];
                    kge_optm(bid,eid) = KGEopt;
                    kge_orim(bid,eid) = KGEori;
                    rain_orim(bid,eid) = rain_ori;
                    rain_optm(bid,eid) = rain_opt;
                    rainsd_orim(bid,eid) = rainsd_ori;
                    rainsd_optm(bid,eid) = rainsd_opt;
                    %}

                    %{
                        
                        f = figure('visible','on');
                                f.Position = [5 95 700 350];
                                
                                pos1 = [0.13 0.8 0.67 0.15];
                                %subplot('Position',pos1)
                                pos2 = [0.15 0.12 0.78 0.42];
                        
                        subplot('Position',pos2)
                        % hjh = plot(2:3:obstep*3,obsdatav1(:,6),'k',1:288,no3,'g',1:288,no2,'r',1:288,no5,'m');
                        % set(hjh,'LineWidth',2.5)
                        % legend('Obs.','IMERG_D','STIV_D','STIV_{DBK}','Location','Northeast')
                        %yyaxis left
                        xj = 2:8:obstep*8;
                        xj = xj';
                        bd1 = 1.05*obsnow;
                        bd2 = 0.95*obsnow;
                        sh = scatter(2:8:obstep*8,obsnow,45,'g','filled');
                        sh.MarkerEdgeColor = [0,1,0];
                        hold on
                        % if size(no3sum,1)<5
                        %     opt_p = 0;% meaning a proper IRC is not done.
                        % end                       
                        mp = min([opt_p+1 10]);

                        hjh = plot(1:8:720,no3sum(mp,1:8:720),'m-.');
                        hold on
                        hjb = plot([1+8*stat_left_bd 1+8*stat_left_bd],[0 1000000],'k-.');
                        hold on
                        hjb2 = plot([1+8*(windowtimes(4)+6) 1+8*(windowtimes(4)+6)],[0 1000000],'k-.');
                        
                        ymm = ceil(max(max(obsnow),max(no3sum(opt_p+1,1:8:720)))/15)*15;
                        
                        
                        hold on
                        hjjh = plot(1:720,nostivdbkc,'k:');
                        set(hjjh,'LineWidth',2.5)
                        ymm = max([ymm,ceil(max(nostivdbkc)/15)*15]);
                        hold off
                        
                        tf = text(16,0.85*ymm,['' num2str(KGEori,'%2.2f')]);
                        tf.FontSize = 14;
                        tf.FontWeight = 'bold';
                        
                        tf = text(16,0.72*ymm,['' num2str(KGEopt,'%2.2f')],'Color','m');
                        tf.FontSize = 14;
                        tf.FontWeight = 'bold';
                        set(hjh,'LineWidth',2.5)
                        %             hold on
                        %             patch([xj; flipud(xj)]', [bd1; flipud(bd2)]',  'k', 'FaceAlpha',0.2, 'EdgeColor','none')
                        %             hold off
                        ylim([0 ymm])
                        ylabel('Streamflow (m^3/s)')
                        
                        box on
                        
                      
                        title([num2str(gid) ' ERA5_{D}^{IRC*} | ' num2str(area,'%0.0f') ' km^2'])

                        xlabel([tstr ' '])
                        %set(gca,'XTick',1:3:288,'XTickLabel',[])
                        xlim([0 719])
                        
                        %set(gca,'XTick',1:36:288,'XTickLabel',['00';'03';'06';'09';'12';'15';'18';'21'],'FontSize',12,'LineWidth',2.5,'FontWeight','bold')
                        %set(gca,'YTick',0:0.5:3,'YTickLabel',0:0.5:3,'FontSize',12,'LineWidth',2.5,'FontWeight','bold')
                        set(gca,'XTick',0:144:720,'XTickLabel',xtlabels(1:6:30),'FontSize',12,'LineWidth',2.5,'FontWeight','bold')
                        
                        
                        hax = gca;
                        hax.XAxis.MinorTickValues = linspace(0,720,31);
                        hax.XMinorTick = 'on';

                        set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',2.5,'FontWeight','bold');
                        set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',2.5,'FontWeight','bold');
                                Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
                                %[Left Bottom Right Top] spacing
                                NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3)-0.01 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
                                set(gca, 'Position', NewPos);
                        
                        %saveas(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/IRC_images/Basin' num2str(gid) '_' rainevent '_' num2str(ICn) '_Woptimum.jpg'])
    
                        
                        %
                        
                        figure
                        %                     t = tiledlayout(1,1,'Padding','none');
                        %                     t.Units = 'inches';
                        %                     t.OuterPosition = [0.15 0.25 4.8 3.8];
                        %                     nexttile;
                        
                        pos2 = [0.11 0.32 0.60 0.45];
                        
                        subplot('Position',pos2)
                        plot(KGEwxy,'-o','MarkerEdgeColor','k','MarkerFaceColor',[0 0 0], 'LineWidth',2)
                        ylabel('KGE')
                        xlabel([tstr ' '])
                        xlim([0 11])
                        xtla = {'W11','W12','W13','W14','W22','W23','W24','W32','W33','W34','W42','W43','W44','W52','W53','W54'};
                        set(gca,'XTick',1:3:15,'XTickLabel',xtla(1:3:15),'FontSize',12,'LineWidth',2.5,'FontWeight','bold')
                        hax = gca;
                        hax.XAxis.MinorTickValues = linspace(1,18,18);
                        hax.XMinorTick = 'on';
                        set(gca,'FontName','Times New Roman','FontSize',16,'LineWidth',3,'FontWeight','bold');
                        %exportgraphics(gcf,['/shared/dondo/home/ml423/HMIOP/NoDAre/LB/pics/Basin' num2str(bid,'%2.2d') '_event' num2str(eid,'%2.2d') '_KGE.jpg'],'Resolution',500)

                    %}

                    %{
                        % plot opt KGE corresponded rain
                        ffnm=['/shared/dondo/home/ml423/world_mts_runs/IRC1/Basin' num2str(gid) 'outputs/refxtmp' num2str(fopt) 'fullresults/' rainevent '_SW/' rainevent '/precipitation'];
                        fid=fopen(ffnm,'rb','ieee-le');
                        tmpdata=fread(fid,inf,'single');
                        fclose(fid);
                        NODAc2o5=reshape(tmpdata,c,r,720);
                        rtot = zeros(c,r);
                        for k = 1:720
                            tmp = NODAc2o5(:,:,k);
                            rtot = rtot+tmp;
                        end
                        rtot(ind) = nan;
                        ama = jet(40);ama(1,:) = [1 1 1];
                        figure
                        imagesc(rtot')
                        colorbar
                        colormap(ama)
                        title(['DCHM rain 42 ' num2str(max(rtot(:)),'%2.1f') ' ' num2str(nanmean(rtot(:)),'%2.1f') 'start' num2str(ori_mean,'%2.1f')])
                        caxis([0.2*min(MCp(:)) 1.1*max(MCp(:))])
                    %}
                    %{
                        
%                         if ~exist(['/shared/dondo/home/ml423/world_mts_runs/IRC1/Basin' num2str(gid) 'outputs/refxtmp' num2str(fopt_last) '/']);
%                             disp(gid)
%                             disp(rainevent)
%                             continue
%                         else
%                             
%                         end
                        
                        cc = ls(['/shared/dondo/home/ml423/world_mts_runs/IRC1/Basin' num2str(gid) 'outputs/refxtmp' num2str(fopt) '/xtmpbb1*']);
                        cc2 = strsplit(cc);
                        tmpC  = struct2cell(load([cc2{1}]));
                        inb2 = tmpC{1};
                        ICC = zeros(c,r);p_ts = [];
                        for i = 1:720
                            tmp = inb2{i};
                            %tmp(tmp>3) = 3;
                            p_ts(i) = max(tmp(intt));
                            ICC = ICC + tmp;
                        end
                        ICC(ind) = nan;
                        figure
                        imagesc(ICC')
                        colorbar
                        colormap(ama)
                        
                        title(['DCHM Input xtmpp 34 ' num2str(max(ICC(:))) ' ' num2str(nanmean(ICC(:)),'%2.1f') 'start' num2str(ori_mean,'%2.1f')])
                        caxis([0.2*min(MCp(:)) 1.1*max(MCp(:))])
                        hold on
                        scatter(streamx,streamy,20,'ko','filled')
                    %}
                end
            end
        end
        [bid,eid]
    end


    %bid
end



%% plot soil moisture outputs as examples 
tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/world_selected_gauge_events_allinfo_v3.mat'],'gauge_events_first_last');
Basins_sum = tmp.gauge_events_first_last;clear tmp;

for gloc = 104%79
    gid = WORLD_events(gloc,1);
    gind = find(Basins_sum(:,1)==gid);
    ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_fixed_900m/soildepth.ascii'];
    soi = dlmread(ffnm);
    dsm = reshape(soi(:,1),Basins_sum(gind,7),Basins_sum(gind,6));
    dp = reshape(soi(:,2),Basins_sum(gind,7),Basins_sum(gind,6));
    d0 = 0.1*ones(Basins_sum(gind,7),Basins_sum(gind,6));
    
    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_ind.mat'],'ind');
    ind = tmp.ind;
    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_intt.mat'],'intt');
    intt = tmp.intt;
    
    IRCICC = 4140;
    for eid = 2
        
        
        eve = WORLD_events(gloc,eid+1);
        
        
        
        % ffnm=['/shared/dondo/home/ml423/world_mts_outputs/Basin' num2str(gid) '_output_900m1hr/' num2str(eve) '_SW/midout/laststep_ws'];
        % fid=fopen(ffnm,'rb','ieee-le');
        % tmpdata=fread(fid,inf,'single');
        % wsori=reshape(tmpdata,Basins_sum(gind,7),Basins_sum(gind,6));
        % ffnm=['/shared/dondo/home/ml423/world_mts_outputs/Basin' num2str(gid) '_output_900m1hr/' num2str(eve) '_SW/midout/laststep_wsm'];
        % fid=fopen(ffnm,'rb','ieee-le');
        % tmpdata=fread(fid,inf,'single');
        % wsmori=reshape(tmpdata,Basins_sum(gind,7),Basins_sum(gind,6));
        % ffnm=['/shared/dondo/home/ml423/world_mts_outputs/Basin' num2str(gid) '_output_900m1hr/' num2str(eve) '_SW/midout/laststep_wsd'];
        % fid=fopen(ffnm,'rb','ieee-le');
        % tmpdata=fread(fid,inf,'single');
        % wsdori=reshape(tmpdata,Basins_sum(gind,7),Basins_sum(gind,6));
        % % spinup
        ffnm=['/shared/dondo/home/ml423/world_mts_outputs/Basin' num2str(gid) '_output_900m1hr_S2_spinup/' num2str(20100831) '/laststep_ws'];
        fid=fopen(ffnm,'rb','ieee-le');
        tmpdata=fread(fid,inf,'single');
        wsorispin=reshape(tmpdata,Basins_sum(gind,7),Basins_sum(gind,6));
        ffnm=['/shared/dondo/home/ml423/world_mts_outputs/Basin' num2str(gid) '_output_900m1hr_S2_spinup/' num2str(20100831) '/laststep_wsb'];
        fid=fopen(ffnm,'rb','ieee-le');
        tmpdata=fread(fid,inf,'single');
        wsborispin=reshape(tmpdata,Basins_sum(gind,7),Basins_sum(gind,6));
        
        ffnm=['/shared/dondo/home/ml423/world_mts_outputs/Basin' num2str(gid) '_output_900m1hr_IC/' num2str(eve) 'c/laststep_ws'];
        fid=fopen(ffnm,'rb','ieee-le');
        tmpdata=fread(fid,inf,'single');
        wsspincopy=reshape(tmpdata,Basins_sum(gind,7),Basins_sum(gind,6));


        ffnm=['/shared/dondo/home/ml423/world_mts_runs/IRC1/Basin6559100outputs/refxtmp1596091600fullresults/20100916_SW/20100916/soilmoisture1'];
        fid=fopen(ffnm,'rb','ieee-le');
        tmpdata=fread(fid,inf,'single');
        ws00=reshape(tmpdata,Basins_sum(gind,7),Basins_sum(gind,6),720);
        ws00tmp = ws00(:,:,1)/10;

        ffnm=['/shared/dondo/home/ml423/world_mts_runs/IRC1/Basin6559100outputs/refxtmp1596091600fullresults/20100916_SW/20100916/soilmoisture4'];
        fid=fopen(ffnm,'rb','ieee-le');
        tmpdata=fread(fid,inf,'single');
        wsb00=reshape(tmpdata,Basins_sum(gind,7),Basins_sum(gind,6),720);
        wsb00tmp = wsb00(:,:,1);
        ffnm=['/shared/dondo/home/ml423/world_mts_runs/IRC1/Basin6559100outputs/refxtmp1596091600fullresults/20100916_SW/20100916/baselayerthick'];
        fid=fopen(ffnm,'rb','ieee-le');
        tmpdata=fread(fid,inf,'single');
        wsbd00=reshape(tmpdata,Basins_sum(gind,7),Basins_sum(gind,6),720);
        wsbd00tmp = wsbd00(:,:,1);
        wsball = wsb00tmp.*wsbd00tmp;

        ffnm=['/shared/dondo/home/ml423/world_mts_runs/IRC1/Basin6559100outputs/refxtmp1596091600fullresults/20100916_SW/20100916/aquiferthick'];
        fid=fopen(ffnm,'rb','ieee-le');
        tmpdata=fread(fid,inf,'single');
        aqt00=reshape(tmpdata,Basins_sum(gind,7),Basins_sum(gind,6),720);
        aqt00tmp = aqt00(:,:,10);
        % intermediate
        % #$^$%#^%#%@#$%@#$%@#$%@#$%@#$%@#$%@#$%@#$%@#$%@#$%@#$%@#$%@#$%
        ffnm=['/shared/dondo/home/ml423/world_mts_outputs/Basin' num2str(gid) '_output_900m1hr_IC/' num2str(eve) 'c/TrueIC/laststep_ws'];
        fid=fopen(ffnm,'rb','ieee-le');
        tmpdata=fread(fid,inf,'single');
        wsint=reshape(tmpdata,Basins_sum(gind,7),Basins_sum(gind,6));
        ffnm=['/shared/dondo/home/ml423/world_mts_outputs/Basin' num2str(gid) '_output_900m1hr_IC/' num2str(eve) 'c/TrueIC/laststep_wsm'];
        fid=fopen(ffnm,'rb','ieee-le');
        tmpdata=fread(fid,inf,'single');
        wsmint=reshape(tmpdata,Basins_sum(gind,7),Basins_sum(gind,6));
        ffnm=['/shared/dondo/home/ml423/world_mts_outputs/Basin' num2str(gid) '_output_900m1hr_IC/' num2str(eve) 'c/TrueIC/laststep_wsd'];
        fid=fopen(ffnm,'rb','ieee-le');
        tmpdata=fread(fid,inf,'single');
        wsdint=reshape(tmpdata,Basins_sum(gind,7),Basins_sum(gind,6));
        
        % figure
        % imagesc(wsori')
        % colorbar
        % figure
        % imagesc(wsint')
        % colorbar
        figure
        imagesc(wsorispin')
        colorbar
        title('ws spin')
        figure
        imagesc(wsborispin')
        colorbar
        title('wsb spin')

        figure
        imagesc(wsspincopy')
        colorbar
        title('ws spin copied')
        figure
        imagesc(ws00tmp')
        colorbar
        title('ws progress')
        clim([0.03 0.059])
        wsb00tmp(wsb00tmp<0)=0;   
        figure
        imagesc(wsb00tmp')
        colorbar
        title('wsb progress')
        figure
        imagesc(wsball')
        colorbar
        title('wsb w depth progress')


        figure
        imagesc(aqt00tmp')
        colorbar
        title('aqt progress')

        figure
        imagesc(wsmori')
        colorbar
        figure
        imagesc(wsmint')
        colorbar
        % output
        ffnm=['/shared/dondo/home/ml423/world_mts_runs/IRC1/ICall/Basin' num2str(gid) '/' num2str(IRCICC) '/' num2str(eve) '/TrueIC/laststep_ws'];
        fid=fopen(ffnm,'rb','ieee-le');
        if fid<0
            continue
        end
        tmpdata=fread(fid,inf,'single');
        ws=reshape(tmpdata,Basins_sum(gind,7),Basins_sum(gind,6));
        ffnm=['/shared/dondo/home/ml423/world_mts_runs/IRC1/ICall/Basin' num2str(gid) '/' num2str(IRCICC) '/' num2str(eve) '/TrueIC/laststep_wsm'];
        fid=fopen(ffnm,'rb','ieee-le');
        tmpdata=fread(fid,inf,'single');
        wsm=reshape(tmpdata,Basins_sum(gind,7),Basins_sum(gind,6));
        ffnm=['/shared/dondo/home/ml423/world_mts_runs/IRC1/ICall/Basin' num2str(gid) '/' num2str(IRCICC) '/' num2str(eve) '/TrueIC/laststep_wsd'];
        fid=fopen(ffnm,'rb','ieee-le');
        tmpdata=fread(fid,inf,'single');
        wsd=reshape(tmpdata,Basins_sum(gind,7),Basins_sum(gind,6));
        
        wsori(ind) = nan;wsmori(ind) = nan;wsdori(ind) = nan;ws(ind) = nan;wsm(ind) = nan;wsd(ind) = nan;
        wsdiff = ws-wsori;wsdiff_norm = wsdiff./d0;
        wsmdiff = wsm-wsmori;wsmdiff_norm = wsmdiff./dsm;
        wsddiff = wsd-wsdori;wsddiff_norm = wsddiff./dp;
        
        maxdiff = max([prctile([wsdiff_norm(intt),wsmdiff_norm(intt),wsddiff_norm(intt)],99)]);
        mindiff = min([prctile([wsdiff_norm(intt),wsmdiff_norm(intt),wsddiff_norm(intt)],1)]);
        
        A = flipud(jet(200));
        A(1,:) = [1 1 1];
        B = flipud(redblue(40));
        B(1,:) = [0.5 0.5 0.5];
        
        for layer = 1
            if layer == 1
                wstmp = wsori; wsnew = ws; wdiff = wsdiff;
                
            end
            if layer == 2
                wstmp = wsmori; wsnew = wsm; wdiff = wsmdiff;
            end
            if layer == 3
                wstmp = wsdori; wsnew = wsd; wdiff = wsddiff;
            end
            cx = max([max(wstmp(:)),max(wsnew(:))]);
            if layer == 1
                
                cx = 0.055;
            end
            
            figure
%             t = tiledlayout(1,1,'Padding','none');
%             t.Units = 'inches';
%             t.OuterPosition = [0.25 0.25 3.4 2.8];
%             nexttile;
            imagesc(wstmp')
            caxis([-0.001 cx])
            hold on
            %text(40,77,rainevent,'FontSize',16)
            hold off
            colormap(A)
            title(['Layer ' num2str(layer) ' ori'])
            axis off
            colorbar
            
            %exportgraphics(gcf,['/shared/dondo/home/ml423/HMIOP/NoDAre/LB/pics/sm/CCB-' rainevent '-original IC layer ' num2str(layer) '.jpg'],'Resolution',800)
            
            figure
%             t = tiledlayout(1,1,'Padding','none');
%             t.Units = 'inches';
%             t.OuterPosition = [0.25 0.25 3.4 2.8];
%             nexttile;
            imagesc(wsnew')
            caxis([-0.001 cx])
            hold on
            %text(40,77,rainevent,'FontSize',16)
            hold off
            colormap(A)
            title(['Layer ' num2str(layer) ' new ' num2str(max(wsnew(:)))])
            axis off
            colorbar
            %saveas(gcf,['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_v2/pic/CCB-' rainevent '-new IC layer ' num2str(layer) ' ' num2str(indic)],'png')
            %exportgraphics(gcf,['/shared/dondo/home/ml423/HMIOP/NoDAre/LB/pics/sm/CCB-' rainevent '-new IC layer ' num2str(layer) ' ' num2str(indic) '.jpg'],'Resolution',800)
            
            
            %maxdiff = max(wdiff(:));% comment out so diff maps have the
            %same plot limit accross layers
            %mindiff = min(wdiff(:));
            
            
            figure
%             t = tiledlayout(1,1,'Padding','none');
%             t.Units = 'inches';
%             t.OuterPosition = [0.25 0.25 3.4 2.8];
%             nexttile;
            imagesc(wdiff')
            %caxis([0.01 0.355])
            if maxdiff<=0.0001
                CO = B(1:20,:);
                colormap(CO)
                cm = max(abs([maxdiff,mindiff]));
                caxis([-1.2*cm 0])
            elseif mindiff>=-0.0001
                
                CO = B(21:40,:);
                colormap(CO)
                cm = max(abs([maxdiff,mindiff]));
                caxis([-0.1*cm 1.2*cm])
            else
                cm = max(abs([maxdiff,mindiff]));
                colormap(B)
                caxis([-1.2*cm 1.2*cm])
            end
            hold on
            %text(40,77,rainevent,'FontSize',16)
            hold off
            title(['Layer ' num2str(layer) ' Diff.'])
            axis off
            colorbar
            %saveas(gcf,['/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_v2/pic/CCB-' rainevent '-diff IC layer ' num2str(layer) ' ' num2str(indic)],'png')
            %exportgraphics(gcf,['/shared/dondo/home/ml423/HMIOP/NoDAre/LB/pics/sm/CCB-' rainevent '-diff IC layer ' num2str(layer) ' ' num2str(indic) '.jpg'],'Resolution',800)
            
        end
    end
    
end


%% plot number of events figures Figure1

Mb = shaperead('world-administrative-boundaries');
% the above shapefile is downloaded from https://public.opendatasoft.com/explore/dataset/world-administrative-boundaries/export/?flg=en-us

figure
m_proj('miller','lat',[-45 55],'lon',[-120 130]);
m_coast('patch',[1 1 1],'edgecolor','k');
m_grid('box','on','linestyle','-','gridcolor','w','linewidth',3,'fontsize',18,'fontweight','bold');

for k=1:length(Mb)
     m_line(Mb(k).X(:),Mb(k).Y(:),'Color','k'); 
end


hold on
m_scatter(gauge_event_loc(:,3),gauge_event_loc(:,2),10,gauge_event_loc(:,4),'filled')

% hold on
% m_pcolor(lon,lat,T_p_sum);
%imagesc(lon,lat,T_p_sum)
%colormap(A)
%colormap(flipud(m_colmap('Blues')))
colormap(jet(30))
xlabel(['Longitude '])
ylabel(['Latitude '])
set(gca,'YDir','normal')
title(['Number of Events'])
caxis([0 100])
colorbar
set(get(colorbar,'Title'),'string','')
set(gca,'FontSize',14)
% set(gca,'YDir','normal')

set(gca,'FontName','Times New Roman','FontSize',16,'LineWidth',2,'FontWeight','bold');

x0=10;
y0=10;
width=1550;
height=900;
set(gcf,'position',[x0,y0,width,height])
%saveas(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/number_of_events'],'png')%
exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/number_of_events_country.jpg'],'Resolution',1000)
                    

%% Alps results dependency on dem and relief
figure

t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.25 0.25 5.0 3.4];
nexttile;
scatter(gauge_event_loc(:,8),gauge_event_loc(:,9),14,'ko','filled')
%title('QPE% - mean DEM')
xlabel('\Delta QPE')
ylabel('DEM (m)')
xlim([-0.5 2.5])
box on
set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',2.5,'FontWeight','bold');
%exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/QPE_dem.jpg'],'Resolution',800)
 
           

figure
scatter(gauge_event_loc(:,8),gauge_event_loc(:,10))
title('QPE% - Relief')


figure
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.25 0.25 5.0 3.4];
nexttile;
scatter(gauge_event_loc(:,8),gauge_event_loc(:,5),14,'ko','filled')
hold on
scatter(gauge_event_loc(:,8),gauge_event_loc(:,6),14,'go','filled')
%title('QPE% - KGE')
xlabel('\Delta QPE')
ylabel('KGE')
ylim([-2.5 1])
xlim([-0.5 2.5])
box on
set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',2.5,'FontWeight','bold');
%exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/QPE_kge.jpg'],'Resolution',800)
 

figure
scatter(gauge_event_loc(:,8),gauge_event_loc(:,11))
title('QPE% - basin area')


figure

t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.25 0.25 5.0 3.4];
nexttile;
scatter(gauge_event_loc(:,8),gauge_event_loc(:,7),14,'ko','filled')
%title('QPE% - mean DEM')
xlabel('\Delta QPE')
ylabel('QPE ori (mm)')
xlim([-0.5 2.5])
box on
set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',2.5,'FontWeight','bold');


%% plot full events but conditional on time interval every 15 years

timestamps = gauge_event_locsum(:,2);
figure
histogram(timestamps,[19740000:10000:20240000])
timeyear = floor(timestamps/10000);
yr1 = find(timeyear>=1975&timeyear<=1989);
yr2 = find(timeyear>=1990&timeyear<=2004);
yr3 = find(timeyear>=2005&timeyear<=2019);
%yr1 = find(timeyear>=1974&timeyear<=1998);
%yr2 = find(timeyear>=1999&timeyear<=2023);
figure
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.25 0.25 5.0 3.4];
nexttile;
scatter(gauge_event_locsum(yr1,6),gauge_event_locsum(yr1,7),14,'bo','filled')
hold on
scatter(gauge_event_locsum(yr2,6),gauge_event_locsum(yr2,7),14,'ko','filled')
hold on
scatter(gauge_event_locsum(yr3,6),gauge_event_locsum(yr3,7),14,'ro','filled')

%title('QPE% - mean DEM')
xlabel('\Delta QPE')
ylabel('DEM (m)')
xlim([-0.5 3])
box on
set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',2.5,'FontWeight','bold');
%exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/fullQPE_dem_years.jpg'],'Resolution',800)
 
f = figure
% t = tiledlayout(1,1,'Padding','none');
% t.Units = 'inches';
% t.OuterPosition = [0.25 0.25 5.2 3.2];
f.Position = [50 50 500 300];
h1 = cdfplot(gauge_event_locsum(yr1,6));
hold on
h2 = cdfplot(gauge_event_locsum(yr2,6));
hold on
h3 = cdfplot(gauge_event_locsum(yr3,6));
set(h1,'Color','b','LineWidth',2)
set(h2,'Color','k','LineWidth',2)
set(h3,'Color','r','LineWidth',2)
xlabel('\Delta QPE')
title('Alps')
ylabel('CDF')
xlim([-1 3])
set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',2.5,'FontWeight','bold');
hax = gca;
hax.XAxis.MinorTickValues = linspace(-10,10,81);
hax.XMinorTick = 'on';
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
%[Left Bottom Right Top] spacing
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);% when use this set the figure visible off for better cropping figure

%exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/Alps_cdfdeltaqpe_more.jpg'],'Resolution',600)
 



figure
h1 = cdfplot(gauge_event_locsum(yr1,5));
hold on
h2 = cdfplot(gauge_event_locsum(yr2,5));
hold on
h3 = cdfplot(gauge_event_locsum(yr3,5));
set(h1,'Color','b')
set(h2,'Color','k')
set(h3,'Color','r')
nanmedian(gauge_event_locsum(yr1,5))
nanmedian(gauge_event_locsum(yr2,5))
nanmedian(gauge_event_locsum(yr3,5))




% plot full events but conditional on time every 5 years (tried 10, 12 not work, 5 not work neither)
timestamps = gauge_event_locsum(:,2);
figure
histogram(timestamps)
timeyear = floor(timestamps/10000);
yr1 = find(timeyear>=1975&timeyear<=1979);
yr2 = find(timeyear>=1980&timeyear<=1984);
yr3 = find(timeyear>=1985&timeyear<=1990);
yr4 = find(timeyear>=1995&timeyear<=2000);
yr5 = find(timeyear>=2001&timeyear<=2005);

figure
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.25 0.25 5.0 3.4];
nexttile;
scatter(gauge_event_locsum(yr1,6),gauge_event_locsum(yr1,7),14,'bo','filled')
hold on
scatter(gauge_event_locsum(yr2,6),gauge_event_locsum(yr2,7),14,'ko','filled')
hold on
scatter(gauge_event_locsum(yr3,6),gauge_event_locsum(yr3,7),14,'ro','filled')
hold on
scatter(gauge_event_locsum(yr4,6),gauge_event_locsum(yr4,7),14,'mo','filled')
%hold on
%scatter(gauge_event_locsum(yr5,6),gauge_event_locsum(yr5,7),14,'co','filled')
%title('QPE% - mean DEM')
xlabel('\Delta QPE')
ylabel('DEM (m)')
xlim([-0.5 3])
box on
set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',2.5,'FontWeight','bold');
%exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/fullQPE_dem_years.jpg'],'Resolution',800)
 
figure
h1 = cdfplot(gauge_event_locsum(yr1,6));
hold on
h2 = cdfplot(gauge_event_locsum(yr2,6));
hold on
h3 = cdfplot(gauge_event_locsum(yr3,6));
hold on
h4 = cdfplot(gauge_event_locsum(yr4,6));
hold on
h5 = cdfplot(gauge_event_locsum(yr5,6));
set(h1,'Color','b')
set(h2,'Color','k')
set(h3,'Color','r')
set(h4,'Color','m')
set(h5,'Color','c')
xlim([-1 3])

figure
cdfplot(gauge_event_locsum(yr3,5))
hold on
cdfplot(gauge_event_locsum(yr2,5))
hold on
cdfplot(gauge_event_locsum(yr1,5))
nanmedian(gauge_event_locsum(yr1,6))
nanmedian(gauge_event_locsum(yr2,6))
nanmedian(gauge_event_locsum(yr3,6))
nanmedian(gauge_event_locsum(yr4,5))



addpath('/shared/dondo/home/ml423/DA/Critical matrices/m_map')

Mb = shaperead('world-administrative-boundaries');
% the above shapefile is downloaded from https://public.opendatasoft.com/explore/dataset/world-administrative-boundaries/export/?flg=en-us

timestamps = gauge_event_locsum(:,2);
basins_actual = unique(gauge_event_locsum(:,1));
figure
histogram(timestamps,[19740000:10000:20240000])
delta_change = [];
for i = 1:length(basins_actual)
    loc = find(gauge_event_locsum(:,1)==basins_actual(i));
    eoc = gauge_event_locsum(loc,2);
    deltap = gauge_event_locsum(loc,5);
    area = gauge_event_locsum(loc,11);
    [~,eind] = sort(eoc);
    raints = deltap(eind);
    [hv,pv,sv] = Mann_Kendall(raints,0.05);


    eoc_year = floor(eoc/10000);
    a = min(eoc_year);b = max(eoc_year);

    int1 = a+floor((b-a)/3)-1;
    int2 = a+2*floor((b-a)/3)-1;

    period1 = find(eoc_year>=a&eoc_year<=int1);
    period2 = find(eoc_year>int1&eoc_year<=int2);
    period3 = find(eoc_year>int2&eoc_year<=b);

    finalchange = (nanmean(deltap(period3))-nanmean(deltap(period1)))/nanmean(deltap(period1));
    delta_change(i,1) = basins_actual(i);
    delta_change(i,2) = gauge_event_loc(find(gauge_event_loc(:,1)==basins_actual(i)),2);
    delta_change(i,3) = gauge_event_loc(find(gauge_event_loc(:,1)==basins_actual(i)),3);
    delta_change(i,4) = a;
    delta_change(i,5) = b;
    delta_change(i,6) = length(period1);
    delta_change(i,7) = length(period3);
    delta_change(i,8) = nanmean(deltap(period1));
    delta_change(i,9) = nanmean(deltap(period3));
    delta_change(i,10) = finalchange;
    delta_change(i,11) = area(1);
    delta_change(i,12) = hv;
    delta_change(i,13) = pv;
    delta_change(i,14) = sv;

    if finalchange>0.25
        delta_change(i,15) = 2;
    elseif finalchange>0.125
        delta_change(i,15) = 1;
    elseif finalchange<-0.25
        delta_change(i,15) = -2;
    elseif finalchange<-0.125
        delta_change(i,15) = -1;
    end

    if delta_change(i,8)<150
        delta_change(i,16) = 1;
    elseif delta_change(i,8)>300
        delta_change(i,16) = 3;
    else
        delta_change(i,16) = 2;
    end


end


figure

m_proj('miller','lat',[alpslat1 alpslat2],'lon',[alpslon1 alpslon2]);
%m_coast('patch',[1 1 1],'edgecolor','k');
m_grid('box','on','linestyle','-','gridcolor','w','linewidth',3,'fontsize',18,'fontweight','bold');

for k=1:length(Mb)
     m_line(Mb(k).X(:),Mb(k).Y(:),'Color','k','LineStyle',':'); 
end
hold on 
m_elev('contour',[500:500:5000],'edgecolor','k','linewidth',1.5,'ShowText','on');
hold on
hm = m_scatter(delta_change(:,3),delta_change(:,2),delta_change(:,16)*20,delta_change(:,15),'filled');


% hold on
% m_pcolor(lon,lat,T_p_sum);
%imagesc(lon,lat,T_p_sum)
%colormap(A)
%colormap(flipud(m_colmap('Blues')))
A = redblue(5);
A(3,:) = [0 0 0];
colormap(A)


xlabel(['Longitude '])
ylabel(['Latitude '])
set(gca,'YDir','normal')
title(['changes in rain'])
caxis([-2 2])
colorbar
set(get(colorbar,'Title'),'string','')
set(gca,'FontSize',14)
% set(gca,'YDir','normal')

set(gca,'FontName','Times New Roman','FontSize',16,'LineWidth',2,'FontWeight','bold');

x0=10;
y0=10;
width=1550;
height=900;
set(gcf,'position',[x0,y0,width,height])
%saveas(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/number_of_events'],'png')%
%exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/Changes_in_the_rain_more.jpg'],'Resolution',1000)
  

% plot this in the Critical matrices/hgt_dem_data folder
figure
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.25 0.25 5.1 4.2];
nexttile;

readhgt(alpslat1:alpslat2,alpslon1:alpslon2)
%title('Extreme rainfall changes')
title('')
hold on 
for k=1:length(Mb)
     line(Mb(k).X(:),Mb(k).Y(:),'Color','y','LineWidth',1,'LineStyle',':'); 
end
hold on
sc = scatter(delta_change(:,3),delta_change(:,2),delta_change(:,16)*8,delta_change(:,15),'filled');
sc.MarkerEdgeColor = 'k';
sc.LineWidth = 1;
ylim([alpslat1 alpslat2])
xlim([alpslon1 alpslon2])
A = redblue(5);
A(3,:) = [1 1 1];
colormap(A)
clim([-2 2])
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
%[Left Bottom Right Top] spacing
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3)-0.15 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);
exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/alps_rainfall_change_mean_more.jpg'],'Resolution',900)


%% read hgt dem file and plot Andes graphs (median KGE, and QPE changes)
% use the readhgt.m file in /DA/Critical matrices/hgt_dem_data/ folder so
% that it is gray scale. 

alpsdem = '/shared/dondo/home/ml423/DA/Critical matrices/hgt_dem_data/alps/';
andesdem = '/shared/dondo/home/ml423/DA/Critical matrices/hgt_dem_data/andes/';
himalayasdem = '/shared/dondo/home/ml423/DA/Critical matrices/hgt_dem_data/himalayas/';
addpath('/shared/dondo/home/ml423/DA/Critical matrices/')

andes1lat1 = -5;
andes1lat2 = 20;
andes1lon1 = -85;
andes1lon2 = -72;
andes2lat1 = -35;
andes2lat2 = -3;
andes2lon1 = -57;
andes2lon2 = -33;


andeslat1 = andes2lat1;
andeslat2 = andes2lat2;
andeslon1 = andes2lon1;
andeslon2 = andes2lon2;
% plot preIRC kge
figure
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.25 0.25 5.1 4.2];
nexttile;

readhgt(andeslat1:andeslat2,andeslon1:andeslon2,'srtm3')
title('The Andes pre-IRC median KGE')
hold on 
for k=1:length(Mb)
     line(Mb(k).X(:),Mb(k).Y(:),'Color','y','LineWidth',1); 
end
hold on
sc = scatter(gauge_event_loc(:,3),gauge_event_loc(:,2),26,gauge_event_loc(:,5),'filled');
sc.MarkerEdgeColor = 'k';
sc.LineWidth = 1;
ylim([andeslat1 andeslat2])
xlim([andeslon1 andeslon2])
colormap(jet(30))
clim([-1 1])
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
%[Left Bottom Right Top] spacing
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3)-0.15 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);
exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/andes_prekge.jpg'],'Resolution',700)


% plot post IRC kge
figure
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.25 0.25 5.1 4.2];
nexttile;
readhgt(andeslat1:andeslat2,andeslon1:andeslon2,'srtm3')
title('The Andes post-IRC median KGE')
hold on 
for k=1:length(Mb)
     line(Mb(k).X(:),Mb(k).Y(:),'Color','y','LineWidth',1); 
end
hold on
sc = scatter(gauge_event_loc(:,3),gauge_event_loc(:,2),26,gauge_event_loc(:,6),'filled');
sc.MarkerEdgeColor = 'k';
sc.LineWidth = 1;
ylim([andeslat1 andeslat2])
xlim([andeslon1 andeslon2])
colormap(jet(30))
clim([-1 1])
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
%[Left Bottom Right Top] spacing
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3)-0.15 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);
exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/andes_postkge.jpg'],'Resolution',700)

% QPE change percentage
figure
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.25 0.25 5.1 4.2];
nexttile;
ax1 = axes;
readhgt(andeslat1:andeslat2,andeslon1:andeslon2,'srtm3')
title('The Andes median QPE change (%)')
hold on 
for k=1:length(Mb)
     line(Mb(k).X(:),Mb(k).Y(:),'Color','y','LineWidth',1); 
end
hold on
sc = scatter(gauge_event_loc(:,3),gauge_event_loc(:,2),0.1*gauge_event_loc(:,7),gauge_event_loc(:,8),'filled');
sc.MarkerEdgeColor = 'k';
sc.LineWidth = 1;
ylim([andeslat1 andeslat2])
xlim([andeslon1 andeslon2])
colormap(redblue(30))
clim([-1 1])
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
%[Left Bottom Right Top] spacing
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3)-0.15 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);
%exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/andes_rain_change.jpg'],'Resolution',700)
%% Andes results dependency on dem and relief
figure

t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.25 0.25 5.0 3.4];
nexttile;
scatter(gauge_event_loc(:,8),gauge_event_loc(:,9),14,'ko','filled')
%title('QPE% - mean DEM')
xlabel('\Delta QPE')
ylabel('DEM (m)')
xlim([-0.5 2.5])
box on
set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',2.5,'FontWeight','bold');
%exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/QPE_dem.jpg'],'Resolution',800)
 
           

figure
scatter(gauge_event_loc(:,8),gauge_event_loc(:,10))
title('QPE% - Relief')


figure
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.25 0.25 5.0 3.4];
nexttile;
scatter(gauge_event_loc(:,8),gauge_event_loc(:,5),14,'ko','filled')
hold on
scatter(gauge_event_loc(:,8),gauge_event_loc(:,6),14,'go','filled')
%title('QPE% - KGE')
xlabel('\Delta QPE')
ylabel('KGE')
ylim([-2.5 1])
xlim([-0.5 2.5])
box on
set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',2.5,'FontWeight','bold');
%exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/QPE_kge.jpg'],'Resolution',800)
 

figure
scatter(gauge_event_loc(:,8),gauge_event_loc(:,11))
title('QPE% - basin area')


figure

t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.25 0.25 5.0 3.4];
nexttile;
scatter(gauge_event_loc(:,8),gauge_event_loc(:,7),14,'ko','filled')
%title('QPE% - mean DEM')
xlabel('\Delta QPE')
ylabel('QPE ori (mm)')
xlim([-0.5 2.5])
box on
set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',2.5,'FontWeight','bold');

%% plot Andes full scatter results, not only median values.
% get the alps results 
gauge_event_locsum = gauge_event_locsumall(find(gauge_event_locsumall(:,1)<6000000),:);

figure
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.25 0.25 5.0 3.4];
nexttile;
scatter(gauge_event_locsum(:,6),gauge_event_locsum(:,9),14,'ko','filled')
%title('QPE% - mean DEM')
xlabel('\Delta QPE')
ylabel('DEM (m)')
xlim([-0.5 3])
box on
set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',2.5,'FontWeight','bold');
exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/fullQPE_dem_andes.jpg'],'Resolution',800)

figure
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.25 0.25 5.0 3.4];
nexttile;
scatter(gauge_event_locsum(:,9),gauge_event_locsum(:,8)./gauge_event_locsum(:,7),14,'ko','filled')
%title('QPE% - mean DEM')
xlabel('DEM')
ylabel('QPE opt std')
%xlim([-0.5 3])
box on
set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',2.5,'FontWeight','bold');

figure
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.25 0.25 5.0 3.4];
nexttile;
scatter(gauge_event_locsum(:,6),gauge_event_locsum(:,5),14,'ko','filled')
%title('QPE% - mean DEM')
xlabel('\Delta QPE')
ylabel('QPE ori (mm)')
xlim([-0.5 3])
box on
set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',2.5,'FontWeight','bold');

figure
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.25 0.25 5.0 3.4];
nexttile;
scatter(gauge_event_locsum(:,6),gauge_event_locsum(:,3),14,'ko','filled')
hold on
scatter(gauge_event_locsum(:,6),gauge_event_locsum(:,4),14,'go','filled')
%title('QPE% - KGE')
xlabel('\Delta QPE')
ylabel('KGE')
ylim([-2.5 1])
xlim([-0.5 3])
box on
set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',2.5,'FontWeight','bold');
exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/fullQPE_kge_andes.jpg'],'Resolution',800)
   


% plot full events but conditional on time every 15 years
timestamps = gauge_event_locsum(:,2);
figure
histogram(timestamps,[19740000:10000:20240000])
timeyear = floor(timestamps/10000);
yr1 = find(timeyear>=1975&timeyear<=1984);
yr2 = find(timeyear>=1985&timeyear<=1994);
yr3 = find(timeyear>=1995&timeyear<=2004);
%yr1 = find(timeyear>=1974&timeyear<=1998);
%yr2 = find(timeyear>=1999&timeyear<=2023);
figure
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.25 0.25 5.0 3.4];
nexttile;
scatter(gauge_event_locsum(yr1,6),gauge_event_locsum(yr1,7),14,'bo','filled')
hold on
scatter(gauge_event_locsum(yr2,6),gauge_event_locsum(yr2,7),14,'ko','filled')
hold on
scatter(gauge_event_locsum(yr3,6),gauge_event_locsum(yr3,7),14,'ro','filled')

%title('QPE% - mean DEM')
xlabel('\Delta QPE')
ylabel('DEM (m)')
xlim([-0.5 3])
box on
set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',2.5,'FontWeight','bold');
%exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/fullQPE_dem_years.jpg'],'Resolution',800)
 
f = figure
% t = tiledlayout(1,1,'Padding','none');
% t.Units = 'inches';
% t.OuterPosition = [0.25 0.25 5.2 3.2];
f.Position = [50 50 500 300];
h1 = cdfplot(gauge_event_locsum(yr1,6));
hold on
h2 = cdfplot(gauge_event_locsum(yr2,6));
hold on
h3 = cdfplot(gauge_event_locsum(yr3,6));
set(h1,'Color','b','LineWidth',2)
set(h2,'Color','k','LineWidth',2)
set(h3,'Color','r','LineWidth',2)
xlabel('\Delta QPE')
title('Andes')
ylabel('CDF')
xlim([-1 3])
set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',2.5,'FontWeight','bold');
hax = gca;
hax.XAxis.MinorTickValues = linspace(-10,10,81);
hax.XMinorTick = 'on';
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
%[Left Bottom Right Top] spacing
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);% when use this set the figure visible off for better cropping figure

exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/Andes_cdfdeltaqpe.jpg'],'Resolution',600)
 

figure
h1 = cdfplot(gauge_event_locsum(yr1,5));
hold on
h2 = cdfplot(gauge_event_locsum(yr2,5));
hold on
h3 = cdfplot(gauge_event_locsum(yr3,5));
set(h1,'Color','b')
set(h2,'Color','k')
set(h3,'Color','r')
nanmedian(gauge_event_locsum(yr1,5))
nanmedian(gauge_event_locsum(yr2,5))
nanmedian(gauge_event_locsum(yr3,5))



% plot full events but conditional on time every 5 years (tried 10, 12 not work, 5 not work neither)
timestamps = gauge_event_locsum(:,2);
figure
histogram(timestamps)
timeyear = floor(timestamps/10000);
yr1 = find(timeyear>=1975&timeyear<=1979);
yr2 = find(timeyear>=1980&timeyear<=1984);
yr3 = find(timeyear>=1985&timeyear<=1990);
yr4 = find(timeyear>=1995&timeyear<=2000);
yr5 = find(timeyear>=2001&timeyear<=2005);

figure
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.25 0.25 5.0 3.4];
nexttile;
scatter(gauge_event_locsum(yr1,6),gauge_event_locsum(yr1,7),14,'bo','filled')
hold on
scatter(gauge_event_locsum(yr2,6),gauge_event_locsum(yr2,7),14,'ko','filled')
hold on
scatter(gauge_event_locsum(yr3,6),gauge_event_locsum(yr3,7),14,'ro','filled')
hold on
scatter(gauge_event_locsum(yr4,6),gauge_event_locsum(yr4,7),14,'mo','filled')
%hold on
%scatter(gauge_event_locsum(yr5,6),gauge_event_locsum(yr5,7),14,'co','filled')
%title('QPE% - mean DEM')
xlabel('\Delta QPE')
ylabel('DEM (m)')
xlim([-0.5 3])
box on
set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',2.5,'FontWeight','bold');
%exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/fullQPE_dem_years.jpg'],'Resolution',800)
 
figure
h1 = cdfplot(gauge_event_locsum(yr1,6));
hold on
h2 = cdfplot(gauge_event_locsum(yr2,6));
hold on
h3 = cdfplot(gauge_event_locsum(yr3,6));
hold on
h4 = cdfplot(gauge_event_locsum(yr4,6));
hold on
h5 = cdfplot(gauge_event_locsum(yr5,6));
set(h1,'Color','b')
set(h2,'Color','k')
set(h3,'Color','r')
set(h4,'Color','m')
set(h5,'Color','c')
xlim([-1 3])

figure
cdfplot(gauge_event_locsum(yr3,5))
hold on
cdfplot(gauge_event_locsum(yr2,5))
hold on
cdfplot(gauge_event_locsum(yr1,5))
nanmedian(gauge_event_locsum(yr1,6))
nanmedian(gauge_event_locsum(yr2,6))
nanmedian(gauge_event_locsum(yr3,6))
nanmedian(gauge_event_locsum(yr4,5))


addpath('/shared/dondo/home/ml423/DA/Critical matrices/m_map')

Mb = shaperead('world-administrative-boundaries');
% the above shapefile is downloaded from https://public.opendatasoft.com/explore/dataset/world-administrative-boundaries/export/?flg=en-us

timestamps = gauge_event_locsum(:,2);
basins_actual = unique(gauge_event_locsum(:,1));
figure
histogram(timestamps,[19740000:10000:20240000])
delta_change = [];
for i = 1:length(basins_actual)
    loc = find(gauge_event_locsum(:,1)==basins_actual(i));
    eoc = gauge_event_locsum(loc,2);
    deltap = gauge_event_locsum(loc,5);
    area = gauge_event_locsum(loc,11);
    [~,eind] = sort(eoc);
    raints = deltap(eind);
    [hv,pv,sv] = Mann_Kendall(raints,0.05);


    eoc_year = floor(eoc/10000);
    a = min(eoc_year);b = max(eoc_year);

    int1 = a+floor((b-a)/3)-1;
    int2 = a+2*floor((b-a)/3)-1;

    period1 = find(eoc_year>=a&eoc_year<=int1);
    period2 = find(eoc_year>int1&eoc_year<=int2);
    period3 = find(eoc_year>int2&eoc_year<=b);

    finalchange = (nanmean(deltap(period3))-nanmean(deltap(period1)))/nanmean(deltap(period1));
    delta_change(i,1) = basins_actual(i);
    delta_change(i,2) = gauge_event_loc(find(gauge_event_loc(:,1)==basins_actual(i)),2);
    delta_change(i,3) = gauge_event_loc(find(gauge_event_loc(:,1)==basins_actual(i)),3);
    delta_change(i,4) = a;
    delta_change(i,5) = b;
    delta_change(i,6) = length(period1);
    delta_change(i,7) = length(period3);
    delta_change(i,8) = nanmean(deltap(period1));
    delta_change(i,9) = nanmean(deltap(period3));
    delta_change(i,10) = finalchange;
    delta_change(i,11) = area(1);
    delta_change(i,12) = hv;
    delta_change(i,13) = pv;
    delta_change(i,14) = sv;

    if finalchange>0.25
        delta_change(i,15) = 2;
    elseif finalchange>0.125
        delta_change(i,15) = 1;
    elseif finalchange<-0.25
        delta_change(i,15) = -2;
    elseif finalchange<-0.125
        delta_change(i,15) = -1;
    end

    if delta_change(i,8)<150
        delta_change(i,16) = 1;
    elseif delta_change(i,8)>300
        delta_change(i,16) = 3;
    else
        delta_change(i,16) = 2;
    end


end


figure

m_proj('miller','lat',[andeslat1 andeslat2],'lon',[andeslon1 andeslon2]);
%m_coast('patch',[1 1 1],'edgecolor','k');
m_grid('box','on','linestyle','-','gridcolor','w','linewidth',3,'fontsize',18,'fontweight','bold');

for k=1:length(Mb)
     m_line(Mb(k).X(:),Mb(k).Y(:),'Color','k','LineStyle',':'); 
end
hold on 
m_elev('contour',[500:500:5000],'edgecolor','k','linewidth',1.5,'ShowText','on');
hold on
hm = m_scatter(delta_change(:,3),delta_change(:,2),delta_change(:,16)*20,delta_change(:,15),'filled');


% hold on
% m_pcolor(lon,lat,T_p_sum);
%imagesc(lon,lat,T_p_sum)
%colormap(A)
%colormap(flipud(m_colmap('Blues')))
A = redblue(5);
A(3,:) = [0 0 0];
colormap(A)


xlabel(['Longitude '])
ylabel(['Latitude '])
set(gca,'YDir','normal')
title(['changes in rain'])
caxis([-2 2])
colorbar
set(get(colorbar,'Title'),'string','')
set(gca,'FontSize',14)
% set(gca,'YDir','normal')

set(gca,'FontName','Times New Roman','FontSize',16,'LineWidth',2,'FontWeight','bold');

x0=10;
y0=10;
width=1550;
height=900;
set(gcf,'position',[x0,y0,width,height])
%saveas(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/number_of_events'],'png')%
%exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/Changes_in_the_rain.jpg'],'Resolution',1000)
  

% plot this in the Critical matrices/hgt_dem_data folder
figure
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.25 0.25 5.1 4.2];
nexttile;


readhgt(andeslat1:andeslat2,andeslon1:andeslon2,'srtm3')
%title('Extreme rainfall changes')
title('')
hold on 
for k=1:length(Mb)
     line(Mb(k).X(:),Mb(k).Y(:),'Color','y','LineWidth',1,'LineStyle',':'); 
end
hold on
sc = scatter(delta_change(:,3),delta_change(:,2),delta_change(:,16)*8,delta_change(:,15),'filled');
sc.MarkerEdgeColor = 'k';
sc.LineWidth = 1;
ylim([andeslat1 andeslat2])
xlim([andeslon1 andeslon2])
A = redblue(5);
A(3,:) = [1 1 1];
colormap(A)
clim([-2 2])
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
%[Left Bottom Right Top] spacing
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3)-0.15 1-Tight(2)-Tight(4)-0.05]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);
exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/andes_rainfall_change_mean.jpg'],'Resolution',900)



%% read hgt dem file and plot Appalahians
% use the readhgt.m file in /DA/Critical matrices/hgt_dem_data/ folder so
% that it is gray scale. 

alpsdem = '/shared/dondo/home/ml423/DA/Critical matrices/hgt_dem_data/alps/';
andesdem = '/shared/dondo/home/ml423/DA/Critical matrices/hgt_dem_data/andes/';
himalayasdem = '/shared/dondo/home/ml423/DA/Critical matrices/hgt_dem_data/himalayas/';
addpath('/shared/dondo/home/ml423/DA/Critical matrices/hgt_dem_data/')
bdbox = [34.2 44.8;-85 -70];

figure
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.25 0.25 5.1 4.6];
nexttile;

readhgt(bdbox(1,1):bdbox(1,2),bdbox(2,1):bdbox(2,2))
set(gca,'XTick',35:2:44,'XTickLabel',[35:2:44])
title('The Appalachians pre-IRC median KGE')
hold on
sc = scatter(gauge_event_loc(:,3),gauge_event_loc(:,2),26,gauge_event_loc(:,5),'filled');
sc.MarkerEdgeColor = 'k';
sc.LineWidth = 1;
ylim([bdbox(1,1) bdbox(1,2)])
xlim([bdbox(2,1) bdbox(2,2)])
colormap(jet(30))
clim([-1 1])
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
%[Left Bottom Right Top] spacing
NewPos = [Tight(1)+0.1 Tight(2)+0.15 1-Tight(1)-Tight(3)-0.25 1-Tight(2)-Tight(4)-0.15]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);
%exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/appalachians_prekge.jpg'],'Resolution',700)

figure
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.25 0.25 5.1 4.6];
nexttile;

readhgt(bdbox(1,1):bdbox(1,2),bdbox(2,1):bdbox(2,2))
title('The Appalachians post-IRC median KGE')
hold on
sc = scatter(gauge_event_loc(:,3),gauge_event_loc(:,2),26,gauge_event_loc(:,6),'filled');
sc.MarkerEdgeColor = 'k';
sc.LineWidth = 1;
ylim([bdbox(1,1) bdbox(1,2)])
xlim([bdbox(2,1) bdbox(2,2)])
colormap(jet(30))
clim([-1 1])
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
%[Left Bottom Right Top] spacing
NewPos = [Tight(1)+0.1 Tight(2)+0.15 1-Tight(1)-Tight(3)-0.25 1-Tight(2)-Tight(4)-0.15]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);
%exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/appalachians_postkge.jpg'],'Resolution',700)



%% [Results section below] append stats to plot global graphs.
addpath('/shared/dondo/home/ml423/DA/ERA5_read/m_map/')
A = jet(30);
A(1,:) = [1 1 1];

tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/world_selected_gauge_events_allinfo_v3.mat'],'gauge_events_first_last');
Basins_sum = tmp.gauge_events_first_last;clear tmp;
tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/WORLD_events_v3.mat'],'WORLD_events');
WORLD_events = tmp.WORLD_events;clear tmp;


home_alps = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Alps/input_output/';
home_andes = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Andes/input_output/';
home_himalayas = '/shared/dondo/home/ml423/HMIOP/NoDAre/IRC_LB/Other_mountains/Himalayas/input_output/';


outlet_all = [];
gaugeid_all = [];
for mountain = [1:3]
    
    if mountain == 1
        home_moun = home_alps;mountname = 'Alps';
    end
    if mountain == 2
        home_moun = home_andes;mountname = 'Andes';
    end
    if mountain == 3
        home_moun = home_himalayas;mountname = 'Himalayas';
    end
    
    outlet_latlon =  csvread([home_moun mountname '_final_outlets_sum.csv'],1,0);
    gauge_id = xlsread([home_moun mountname '_final_basins.xlsx']);
    outlet_all = [outlet_all;outlet_latlon];
    gaugeid_all = [gaugeid_all;gauge_id];
end
outlet_all = [outlet_all,gaugeid_all(:,2)];
gauge_event_loc = [];
for i = 1:size(WORLD_events,1)
    gauge_event_loc(i,1) = WORLD_events(i,1);
    gauge_event_loc(i,2) = outlet_all(find(gaugeid_all(:,2)==WORLD_events(i,1)),1);
    gauge_event_loc(i,3) = outlet_all(find(gaugeid_all(:,2)==WORLD_events(i,1)),2);
    gauge_event_loc(i,4) = length(find(WORLD_events(i,:)>1))-1;
end

for i = 1:length(outlet_all)
    tmp = find(gauge_events_first_last(:,1)==outlet_all(i,4));
    if ~isempty(tmp)
        outlet_all(i,5) = gauge_events_first_last(tmp,5);
    end
end

outletsgis = outlet_all(find(outlet_all(:,5)>0&outlet_all(:,5)<1000),:);
%csvwrite('outletpaper.csv',outletsgis)

% APL outlets

addpath('/shared/dondo/home/ml423/DA/Critical matrices/')

APL_LAT = imread('APL_LAT','tif');
APL_LON = imread('APL_LON','tif');

APL_gauges = load(['APL_gauges.mat'],'APL_gauges');% outlet locations in the tif images
APL_gauges = APL_gauges.APL_gauges;

for i = 1:length(APL_gauges)
    latloninfo(i,:) = [APL_LAT(APL_gauges(i,1),APL_gauges(i,2)),APL_LON(APL_gauges(i,1),APL_gauges(i,2))];
end
for i = 1:length(APL_gauges)
    if ~ismember(i,[7 12])
        gauge_event_loc = [gauge_event_loc;[i,latloninfo(i,1),latloninfo(i,2),length(find(APL_events(i,:)>0))]];
    end
end
%% append stats plot KGE values pre post-IRC,QPE change, DEM, relief, basin area in Alps
gauge_event_locsum1 = [];gauge_event_locsum2 = [];gauge_event_locsum3 = [];gauge_event_locsum4 = [];gauge_event_locsum5 = [];gauge_event_locsum6 = [];gauge_event_locsum7 = [];gauge_event_locsum8 = [];gauge_event_locsum9 = [];gauge_event_locsum10 = [];gauge_event_locsum11 = [];
gauge_event_loc(:,5) = nan;gauge_event_loc(:,6) = nan;
gauge_event_loc(:,7) = nan;gauge_event_loc(:,8) = nan;
gauge_event_loc(:,9) = nan;gauge_event_loc(:,10) = nan;%9 is dem and 10 is relief
kge_orim(kge_orim==0) = nan;kge_optm(kge_optm==0) = nan;
rain_orim(rain_orim==0) = nan;rain_optm(rain_optm==0) = nan;
%eqcontrol = event_quality(:,2:150);
%kge_orim(isnan(eqcontrol)) = nan;kge_optm(isnan(eqcontrol)) = nan;
%rain_orim(isnan(eqcontrol)) = nan;rain_optm(isnan(eqcontrol)) = nan;
eqcontrol = event_quality(:,2:150);
for i = 1:length(event_quality_window)
    eqcontrol(event_quality_window(i,1),event_quality_window(i,2)) = nan;
end
for i = 1:length(event_flash)
    eqcontrol(event_flash(i,1),event_flash(i,2)) = nan;
end
for i = 1:length(event_same)
    eqcontrol(event_same(i,1),event_same(i,2)) = nan;
end
for i = 1:length(reprob_150)
    eqcontrol(reprob_150(i,1),reprob_150(i,2)) = nan;
end


rain_prct = (rain_optm-rain_orim)./rain_orim;
% eqc control for basin 130 131 are all nan


for i = [1:129 132:133 466:480 482:504 506:522]
    % get these KGE values from world IRC
    loconan = find(~isnan(eqcontrol(i,:)));    
    gauge_event_loc(find(gauge_event_loc(:,1)==WORLD_events(i,1)),5)=nanmedian(kge_orim(i,loconan));
    gauge_event_loc(find(gauge_event_loc(:,1)==WORLD_events(i,1)),6)=nanmedian(kge_optm(i,loconan));
    gauge_event_loc(find(gauge_event_loc(:,1)==WORLD_events(i,1)),7)=nanmedian(rain_orim(i,loconan));
    gauge_event_loc(find(gauge_event_loc(:,1)==WORLD_events(i,1)),8)=nanmedian(rain_prct(i,loconan));
    
    loco = find(~isnan(kge_orim(i,:))&~isnan(eqcontrol(i,:)));
    gauge_event_locsum1 = [gauge_event_locsum1;kge_orim(i,loco)'];
    gauge_event_locsum2 = [gauge_event_locsum2;kge_optm(i,loco)'];
    gauge_event_locsum3 = [gauge_event_locsum3;rain_orim(i,loco)'];
    gauge_event_locsum4 = [gauge_event_locsum4;rain_prct(i,loco)'];
    gauge_event_locsum5 = [gauge_event_locsum5;rainsd_orim(i,loco)'];
    gauge_event_locsum6 = [gauge_event_locsum6;rainsd_optm(i,loco)'];
    gauge_event_locsum7 = [gauge_event_locsum7;ones(length(loco),1).*WORLD_events(i,1)];
    gauge_event_locsum8 = [gauge_event_locsum8;ones(length(loco),1).*gauge_events_first_last(find(gauge_events_first_last(:,1)==WORLD_events(i,1)),10)];
    gauge_event_locsum9 = [gauge_event_locsum9;ones(length(loco),1).*gauge_events_first_last(find(gauge_events_first_last(:,1)==WORLD_events(i,1)),11)];
    gauge_event_locsum10 = [gauge_event_locsum10;ones(length(loco),1).*gauge_events_first_last(find(gauge_events_first_last(:,1)==WORLD_events(i,1)),5)];
    % store event number
    gauge_event_locsum11 = [gauge_event_locsum11;WORLD_events(i,1+loco)'];
end
gauge_event_locsumall = [gauge_event_locsum7,gauge_event_locsum11,gauge_event_locsum1,gauge_event_locsum2,gauge_event_locsum3,gauge_event_locsum4,gauge_event_locsum5,gauge_event_locsum6,gauge_event_locsum8,gauge_event_locsum9,gauge_event_locsum10];
clear gauge_event_locsum1 gauge_event_locsum2 gauge_event_locsum3 gauge_event_locsum4 gauge_event_locsum5 gauge_event_locsum6 gauge_event_locsum7 gauge_event_locsum8 
% for i = 1:size(KGE_ori,1)
%     % get these KGE values from IRCICC Paper IC IRC run
%     gauge_event_loc(find(gauge_event_loc(:,1)==i),5)=nanmedian(KGE_ori(i,:));
%     gauge_event_loc(find(gauge_event_loc(:,1)==i),6)=nanmedian(KGE_new(i,:));
% end
% dem and relief
for i = 1:size(gauge_event_loc,1)
    locc = find(gauge_events_first_last(:,1)==gauge_event_loc(i,1));
    if ~isempty(locc)
        gauge_event_loc(i,9) = gauge_events_first_last(find(gauge_events_first_last(:,1)==gauge_event_loc(i,1)),10);
        gauge_event_loc(i,10) = gauge_events_first_last(find(gauge_events_first_last(:,1)==gauge_event_loc(i,1)),11);
        gauge_event_loc(i,11) = gauge_events_first_last(find(gauge_events_first_last(:,1)==gauge_event_loc(i,1)),5);%basin area
    else
        %disp(gauge_event_loc(i,1))
    end
end


figure

m_proj('miller','lat',[39 52],'lon',[-10 20]);
m_coast('patch',[1 1 1],'edgecolor','k');
m_grid('box','on','linestyle','-','gridcolor','w','linewidth',2);
%m_elev('contourf',[500:500:5000]);

hold on
m_scatter(gauge_event_loc(:,3),gauge_event_loc(:,2),18,gauge_event_loc(:,5),'filled')

% hold on
% m_pcolor(lon,lat,T_p_sum);
%imagesc(lon,lat,T_p_sum)
%colormap(A)
%colormap(flipud(m_colmap('Blues')))
colormap(jet(30))
%colormap([m_colmap('blues',32);m_colmap('gland',128)]); 
xlabel(['Longitude '])
ylabel(['Latitude '])
set(gca,'YDir','normal')
title(['Number of Events'])
caxis([-1 1])
colorbar
set(get(colorbar,'Title'),'string','')
set(gca,'FontSize',14)
% set(gca,'YDir','normal')

set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',2,'FontWeight','bold');

x0=10;
y0=10;
width=1050;
height=600;
set(gcf,'position',[x0,y0,width,height])
%saveas(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/Alps_preIRC_kge'],'png')%
%exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/number_of_events.jpg'],'Resolution',800)
 

figure

m_proj('miller','lat',[39 52],'lon',[-10 20]);
m_coast('patch',[1 1 1],'edgecolor','k');
m_grid('box','on','linestyle','-','gridcolor','w','linewidth',2);
%m_elev('contourf',[500:500:5000]);

hold on
m_scatter(gauge_event_loc(:,3),gauge_event_loc(:,2),18,gauge_event_loc(:,6),'filled')

% hold on
% m_pcolor(lon,lat,T_p_sum);
%imagesc(lon,lat,T_p_sum)
%colormap(A)
%colormap(flipud(m_colmap('Blues')))
colormap(redblue(30))
%colormap([m_colmap('blues',32);m_colmap('gland',128)]); 
xlabel(['Longitude '])
ylabel(['Latitude '])
set(gca,'YDir','normal')
title(['Number of Events'])
caxis([-1 1])
colorbar
set(get(colorbar,'Title'),'string','')
set(gca,'FontSize',14)
% set(gca,'YDir','normal')

set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',2,'FontWeight','bold');

x0=10;
y0=10;
width=1050;
height=600;
set(gcf,'position',[x0,y0,width,height])
%saveas(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/Alps_postIRC_kgeredblue'],'png')%
%exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/number_of_events.jpg'],'Resolution',800)
 
% find problem events 
be_prob = gauge_event_locsumall(find(gauge_event_locsumall(:,6)>2),[1 2 6]);

for i = 1:length(be_prob)
    bid = find(WORLD_events(:,1)==be_prob(i,1));
    be_prob(i,4) = bid;% basin number
    be_prob(i,5) = find(WORLD_events(bid,:)==be_prob(i,2))-1;% event number
end

%% [Figure1] read hgt dem file and plot alps graphs (median KGE, and QPE changes)
% use the readhgt.m file in /DA/Critical matrices/hgt_dem_data/ folder so
% that it is gray scale. 
% gauge_event_locsumall 7000row
%basinID, event, oriKGE, optKGE, rainori, rainchange, rainorisd,rainoptsd,
%DEM relief area

% gauge_event_loc 465row
%basinID, LAT, LON, #of events, oriKGE, optKGE,rainori,rainchange,dem,relief,area 


alpsdem = '/shared/dondo/home/ml423/DA/Critical matrices/hgt_dem_data/alps/';
andesdem = '/shared/dondo/home/ml423/DA/Critical matrices/hgt_dem_data/andes/';
himalayasdem = '/shared/dondo/home/ml423/DA/Critical matrices/hgt_dem_data/himalayas/';
addpath('/shared/dondo/home/ml423/DA/Critical matrices/')
Mb = shaperead('world-administrative-boundaries');

alpslat1 = 42;
alpslat2 = 49;
alpslon1 = -3;
alpslon2 = 17;
alpslat1 = 43.5;
alpslat2 = 49;
alpslon1 = 4;
alpslon2 = 17;

% ori ERA5 amounts

figure
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.25 0.25 5.1 4.2];
nexttile;
readhgt(alpslat1:alpslat2,alpslon1:alpslon2)
alpha(.2)
title(['Median QPE (mm)'],'FontWeight','normal')
hold on 
for k=1:length(Mb)
     line(Mb(k).X(:),Mb(k).Y(:),'Color','k','LineWidth',0.5,'LineStyle',':'); 
end
hold on
sc = scatter(gauge_event_loc(:,3),gauge_event_loc(:,2),26,gauge_event_loc(:,7),'filled');
sc.MarkerEdgeColor = 'k';
sc.LineWidth = 1;
ylim([alpslat1 alpslat2])
xlim([alpslon1 alpslon2])
colormap(jet(20))
clim([80 280])
cb = colorbar;
a = cb.Position;
set(cb,'Position',[a(1)+0.03 a(2)-0.06 0.02 0.65])
set(gca,'FontName','Times New Roman','FontSize',14)%,'LineWidth',2.5);%,'FontWeight','bold');
box on  
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
%[Left Bottom Right Top] spacing
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3)-0.15 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);
                          
exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/alps_QPE_mean.jpg'],'Resolution',500)

% correction percentage
figure
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.25 0.25 5.1 4.2];
nexttile;

readhgt(alpslat1:alpslat2,alpslon1:alpslon2)
alpha(.2)
title('Median   QPE (%)','FontWeight','normal')

hold on 
for k=1:length(Mb)
     line(Mb(k).X(:),Mb(k).Y(:),'Color','k','LineWidth',0.5,'LineStyle',':'); 
end
hold on
sc = scatter(gauge_event_loc(:,3),gauge_event_loc(:,2),26,100*gauge_event_loc(:,8),'filled');
sc.MarkerEdgeColor = 'k';
sc.LineWidth = 1;
ylim([alpslat1 alpslat2])
xlim([alpslon1 alpslon2])
rb = redblue(10);
colormap(rb(4:9,:))
clim([-50 100])
cb = colorbar;
a = cb.Position;
set(cb,'Position',[a(1)+0.03 a(2)-0.06 0.02 0.65])
set(gca,'FontName','Times New Roman','FontSize',14)%,'LineWidth',2.5);%,'FontWeight','bold');
box on  
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
%[Left Bottom Right Top] spacing
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3)-0.15 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);
                          
exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/alps_QPE_percentage.jpg'],'Resolution',500)

figure
%legend only
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.25 0.25 5.1 4.2];
nexttile;
peaks(20)
rb = redblue(120);
colormap(rb(50:120,:))
clim([-0.5 3])
% colormap(jet(256))
% clim([80 280])
ccc = gray(256);
ccc = [ones(40,3);ccc(41:256,:)];
colormap(ccc)
clim([0 2900])
colorbar
%exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/legend_elev.jpg'],'Resolution',500)

% plot post IRC kge
figure
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.25 0.25 5.1 4.2];
nexttile;
ax1 = axes;
readhgt(alpslat1:alpslat2,alpslon1:alpslon2)
title('The Alps post-IRC median KGE')
hold on 
for k=1:length(Mb)
     line(Mb(k).X(:),Mb(k).Y(:),'Color','y','LineWidth',1); 
end
hold on
sc = scatter(gauge_event_loc(:,3),gauge_event_loc(:,2),26,gauge_event_loc(:,6),'filled');
sc.MarkerEdgeColor = 'k';
sc.LineWidth = 1;
ylim([alpslat1 alpslat2])
xlim([alpslon1 alpslon2])
colormap(jet(30))
%clim([-1 1])
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
%[Left Bottom Right Top] spacing
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3)-0.15 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);
%exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/alps_postkge_more.jpg'],'Resolution',700)

% QPE change percentage
figure
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.25 0.25 5.1 4.2];
nexttile;
ax1 = axes;
readhgt(alpslat1:alpslat2,alpslon1:alpslon2)
title('The Alps median QPE change (%)')
hold on 
for k=1:length(Mb)
     line(Mb(k).X(:),Mb(k).Y(:),'Color','y','LineWidth',1); 
end
hold on
sc = scatter(gauge_event_loc(:,3),gauge_event_loc(:,2),0.1*gauge_event_loc(:,7),gauge_event_loc(:,8),'filled');
sc.MarkerEdgeColor = 'k';
sc.LineWidth = 1;
ylim([alpslat1 alpslat2])
xlim([alpslon1 alpslon2])
colormap(redblue(30))
%clim([-1 1])
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
%[Left Bottom Right Top] spacing
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3)-0.15 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);
%exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/alps_rain_change.jpg'],'Resolution',700)


%% [Figure1] plot 15 15 15 oginal QPE and QPE histogram

% gauge_event_locsumall 7000row
%basinID, event, oriKGE, optKGE, rainori, rainchange, rainorisd,rainoptsd,
%DEM relief area

% gauge_event_loc 465row
%basinID, LAT, LON, #of events, oriKGE, optKGE,rainori,rainchange,dem,relief,area 

gauge_15 = [];

for i = 1:length(WORLD_events)
    gid = WORLD_events(i,1);
    gauge_15(i,1) = gid;

    gauge_15(i,2) = gauge_event_loc(i,2);

    gauge_15(i,3) = gauge_event_loc(i,3);

    esum = find(gauge_event_locsumall==gid);
    if ~isempty(esum)
        eind = find(gauge_event_locsumall(:,1)==gid&gauge_event_locsumall(:,2)>=19750000&gauge_event_locsumall(:,2)<=19890000);
        if ~isempty(eind)
            gauge_15(i,4) = mean(gauge_event_locsumall(eind,5));
            gauge_15(i,7) = mean(gauge_event_locsumall(eind,6));
        else
            gauge_15(i,4) = nan;
            gauge_15(i,7) = nan;
        end

        eind = find(gauge_event_locsumall(:,1)==gid&gauge_event_locsumall(:,2)>=19900000&gauge_event_locsumall(:,2)<=20040000);
        if ~isempty(eind)
            gauge_15(i,5) = mean(gauge_event_locsumall(eind,5));
            gauge_15(i,8) = mean(gauge_event_locsumall(eind,6));
        else
            gauge_15(i,5) = nan;
            gauge_15(i,8) = nan;
        end

        eind = find(gauge_event_locsumall(:,1)==gid&gauge_event_locsumall(:,2)>=20050000&gauge_event_locsumall(:,2)<=20190000);
        if ~isempty(eind)
            gauge_15(i,6) = mean(gauge_event_locsumall(eind,5));
            gauge_15(i,9) = mean(gauge_event_locsumall(eind,6));
        else
            gauge_15(i,6) = nan;
            gauge_15(i,9) = nan;
        end

    else
        gauge_15(i,4) = nan;
        gauge_15(i,5) = nan;
        gauge_15(i,6) = nan;
        gauge_15(i,7) = nan;
        gauge_15(i,8) = nan;
        gauge_15(i,9) = nan;
    end

end

alpslat1 = 43;
alpslat2 = 49;
alpslon1 = 4;
alpslon2 = 17;
% plot preIRC kge
figure
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.25 0.25 5.1 4.2];
nexttile;

readhgt(alpslat1:alpslat2,alpslon1:alpslon2)
title('The Alps pre-IRC median KGE')
hold on 
for k=1:length(Mb)
     line(Mb(k).X(:),Mb(k).Y(:),'Color','y','LineWidth',1); 
end
hold on
sc = scatter(gauge_15(:,3),gauge_15(:,2),26,gauge_15(:,6),'filled');
sc.MarkerEdgeColor = 'k';
sc.LineWidth = 1;
ylim([alpslat1 alpslat2])
xlim([alpslon1 alpslon2])
colormap(jet(30))
clim([100 300])
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
%[Left Bottom Right Top] spacing
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3)-0.15 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);
%exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/alps_15_3.jpg'],'Resolution',700)


figure
histogram(gauge_15(:,4),[60:10:350])
hold on
histogram(gauge_15(:,5),[60:10:350])
hold on
histogram(gauge_15(:,6),[60:10:350])

figure
cdfplot(gauge_15(:,4))
hold on
cdfplot(gauge_15(:,5))
hold on
cdfplot(gauge_15(:,6))

figure
cdfplot(gauge_15(:,7))
hold on
cdfplot(gauge_15(:,8))
hold on
cdfplot(gauge_15(:,9))
xlim([0 3])

% entire region based analysis
loc1 = find(gauge_event_locsumall(:,2)>=19750000&gauge_event_locsumall(:,2)<=19890000);
loc2 = find(gauge_event_locsumall(:,2)>=19900000&gauge_event_locsumall(:,2)<=20040000);
loc3 = find(gauge_event_locsumall(:,2)>=20050000&gauge_event_locsumall(:,2)<=20190000);

figure
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.25 0.25 5.5 3.9];
nexttile;
h1 = cdfplot(gauge_event_locsumall(loc1,6)*100);
hold on
h2 = cdfplot(gauge_event_locsumall(loc2,6)*100);
hold on
h3 = cdfplot(gauge_event_locsumall(loc3,6)*100);
set(h1,'Color','b','LineWidth',2)
set(h2,'Color','k','LineWidth',2)
set(h3,'Color','r','LineWidth',2)

xlabel('\Delta QPE (%)')
title(' ')
ylabel('CDF')
xlim([-100 100])
ylim([0 1.02])
set(gca,'FontName','Times New Roman','FontSize',16)%,'LineWidth',2.5,'FontWeight','bold');
hax = gca;
hax.XAxis.MinorTickValues = linspace(-100,100,21);
hax.XMinorTick = 'on';
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
%[Left Bottom Right Top] spacing
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);% when use this set the figure visible off for better cropping figure

exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/alps_QPE_percentage_cdf.jpg'],'Resolution',600)
 
figure
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.25 0.25 5.5 2.5];
nexttile;
h1 = cdfplot(gauge_event_locsumall(loc1,5));
hold on
h2 = cdfplot(gauge_event_locsumall(loc2,5));
hold on
h3 = cdfplot(gauge_event_locsumall(loc3,5));
set(h1,'Color','b','LineWidth',2)
set(h2,'Color','k','LineWidth',2)
set(h3,'Color','r','LineWidth',2)

xlabel('QPE (mm)')
title(' ')
ylabel('CDF')
xlim([0 400])
ylim([0 1.02])
set(gca,'FontName','Times New Roman','FontSize',16)%,'LineWidth',2.5,'FontWeight','bold');
hax = gca;
hax.XAxis.MinorTickValues = linspace(-100,100,21);
hax.XMinorTick = 'on';
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
%[Left Bottom Right Top] spacing
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);% when use this set the figure visible off for better cropping figure

exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/alps_QPE_3interv_cdf.jpg'],'Resolution',600)
 


for iy = 1:50
    loc = find(gauge_event_locsumall(:,2)>=(19740000+10000*iy)&gauge_event_locsumall(:,2)<=(19750000+iy*10000));
    extyear(iy) = nanmedian(gauge_event_locsumall(loc,5));
end
%
figure
set(gcf,'outerposition', [100 100 800 400]);
% t = tiledlayout(1,1,'Padding','none');
% t.Units = 'inches';
% t.OuterPosition = [0.15 0.15 5.9 3.9];
% t.Units = 'points';
% t.OuterPosition = [0.15 0.15 800 450];

nexttile;

bar(1975:2024,extyear,'FaceColor',[34 71 91]/100)
hold on
hp = plot([1974 2022],[nanmean(extyear) nanmean(extyear)],'Color','k');%[0 45 69]/100);
hold off
xlabel('Year')
title('')

ylabel('Median QPE (mm)')
set(gca,'FontName','Times New Roman','FontSize',16)%,'LineWidth',2.5,'FontWeight','bold');
xticks([1975:10:1990 1992 1997 2005 2010:10:2025])
xticklabels([1975:10:1990 1992 1997 2005 2010:10:2025])

xlim([1974 2022])
%set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
%[Left Bottom Right Top] spacing
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);% when use this set the figure visible off for better cropping figure

exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/alps_QPE_time_series.jpg'],'Resolution',600)
 


for iy = 1:50
    loc = find(gauge_event_locsumall(:,2)>=(19740000+10000*iy)&gauge_event_locsumall(:,2)<=(19750000+iy*10000));
    tmp = [];
    for j = 1:length(loc)
        lontmp = gauge_event_loc(find(gauge_event_loc(:,1) == gauge_event_locsumall(loc(j),1)),3);
        if lontmp>9
            tmp = [tmp;gauge_event_locsumall(loc(j),5)];
        end
    end
    extyeareast(iy) = nanmedian(tmp);
end
%
figure
set(gcf,'outerposition', [100 100 800 400]);
% t = tiledlayout(1,1,'Padding','none');
% t.Units = 'inches';
% t.OuterPosition = [0.15 0.15 5.9 3.9];
% t.Units = 'points';
% t.OuterPosition = [0.15 0.15 800 450];

nexttile;

bar(1975:2024,extyeareast,'FaceColor',[34 71 91]/100)
hold on
hp = plot([1974 2022],[nanmean(extyeareast) nanmean(extyeareast)],'Color','k');%[0 45 69]/100);
hold off
xlabel('Year')
title('')

ylabel('Median QPE (mm)')

%% [Figure2] plot alps full conditional on DEM band, not only median values.
% get the alps results 
gauge_event_locsum = gauge_event_locsumall(find(gauge_event_locsumall(:,1)>6000000),:);

figure
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.25 0.25 5.0 3.4];
nexttile;
scatter(gauge_event_locsum(:,6),gauge_event_locsum(:,9),14,'ko','filled')
%title('QPE% - mean DEM')
xlabel('\Delta QPE')
ylabel('DEM (m)')
xlim([-0.5 3])
box on
set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',2.5,'FontWeight','bold');
%exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/fullQPE_dem_more.jpg'],'Resolution',800)
%
gauge_event_locsum = gauge_event_locsumall(find(gauge_event_locsumall(:,1)>6000000),:);
% years bracket below
%gauge_event_locsum = gauge_event_locsumall(find(gauge_event_locsumall(:,1)>6000000&gauge_event_locsumall(:,2)>=19750000&gauge_event_locsumall(:,2)<=19890000),:);
%gauge_event_locsum = gauge_event_locsumall(find(gauge_event_locsumall(:,1)>6000000&gauge_event_locsumall(:,2)>=19900000&gauge_event_locsumall(:,2)<=20040000),:);
%gauge_event_locsum = gauge_event_locsumall(find(gauge_event_locsumall(:,1)>6000000&gauge_event_locsumall(:,2)>=20050000&gauge_event_locsumall(:,2)<=20190000),:);

numband = 15;
demint = 200;
demgroups = nan(numband,5000);% QPE change
demgroups2 = nan(numband,5000);% QPE mm
demgroups3 = nan(numband,5000);% ori KGE
demgroups4 = nan(numband,5000);% opt KGE
eventsum = nan(numband,1);
eventsum_year = nan(numband,3);
for demband = 1:numband
    demloc = find(gauge_event_locsum(:,9)>=(200+(demband-1)*demint)&gauge_event_locsum(:,9)<(200+demband*demint));
    if ~isempty(demloc)
        demgroups(demband,1:length(demloc)) = gauge_event_locsum(demloc,6)';
        demgroups2(demband,1:length(demloc)) = gauge_event_locsum(demloc,5)';
        demgroups3(demband,1:length(demloc)) = gauge_event_locsum(demloc,3)';
        demgroups4(demband,1:length(demloc)) = gauge_event_locsum(demloc,4)';
        
    end
    eventsum(demband) = length(find(demgroups(demband,:)>-99));
    eventsum_year(demband,1) = length(find(gauge_event_locsum(:,9)>=(200+(demband-1)*demint)&gauge_event_locsum(:,9)<(200+demband*demint)&gauge_event_locsum(:,2)>=19750000&gauge_event_locsum(:,2)<=19890000));
    eventsum_year(demband,2) = length(find(gauge_event_locsum(:,9)>=(200+(demband-1)*demint)&gauge_event_locsum(:,9)<(200+demband*demint)&gauge_event_locsum(:,2)>=19900000&gauge_event_locsum(:,2)<=20040000));
    eventsum_year(demband,3) = length(find(gauge_event_locsum(:,9)>=(200+(demband-1)*demint)&gauge_event_locsum(:,9)<(200+demband*demint)&gauge_event_locsum(:,2)>=20050000&gauge_event_locsum(:,2)<=20190000));

end
% 
bidbigcell = {};
ibinfocell = {};

for demband = 1:numband-1
    demloc = find(gauge_event_locsum(:,9)>=(200+(demband-1)*demint)&gauge_event_locsum(:,9)<(200+demband*demint));

    deltabig = (gauge_event_locsum(demloc,1));
    bidbig = [];
    for i = 1:length(deltabig)
        bidbig(i,1) = find(WORLD_events(:,1)==deltabig(i));
        bidbig(i,2) = (gauge_event_locsum(demloc(i),6));
        bidbig(i,3) = find(WORLD_events(bidbig(i,1),:)==gauge_event_locsum(demloc(i),2))-1;
    end
    ib = unique(bidbig(:,1));ibinfo = [];
    for j = 1:length(ib)
        ibinfo(j,:) = [ib(j),median(bidbig(find(bidbig(:,1)==ib(j)),2))];
    end
    bidbigcell{demband} = bidbig;
    ibinfocell{demband} = ibinfo;
end

xts = 300:200:3200;
%DEM - QPE change below
figure
%set(gcf,'outerposition', [100 100 700 500]);
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.15 0.15 3.9 5.9];
nexttile;
boxplot(demgroups',xts,'Symbol','k.')
%xlim([200 3200])
ylim([-1 2])

xlabel('DEM (m)')
ylabel('\Delta QPE')
view([90 -90])
set(gca,'FontName','Times New Roman','FontSize',14)
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
%[Left Bottom Right Top] spacing
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);% when use this set the figure visible off for better cropping figure
%set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/alps_DEM_boxplot.jpg'],'Resolution',600)



% dem ~ QPE below
figure
%set(gcf,'outerposition', [100 100 700 500]);
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.15 0.15 3.9 5.9];
nexttile;
boxplot(demgroups2',xts,'Symbol','k.')
%xlim([200 3200])
ylim([0 500])
xlabel('DEM (m)')
ylabel('QPE (mm)')
view([90 -90])
set(gca,'FontName','Times New Roman','FontSize',14)
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
%[Left Bottom Right Top] spacing
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);% when use this set the figure visible off for better cropping figure
%set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/alps_DEM_boxplot2.jpg'],'Resolution',600)

% number of events below
figure
%set(gcf,'outerposition', [100 100 500 800]);
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.15 0.15 3.9 5.9];
nexttile;
hb = bar(xts,eventsum,1);
set(hb,'FaceColor',[82 94 95]/100)
xlim([200 3200])
ylim([0 1500])
xlabel('DEM (m)')
ylabel('Number of Events')
view([90 -90])
set(gca,'FontName','Times New Roman','FontSize',14)
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
%[Left Bottom Right Top] spacing
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);% when use this set the figure visible off for better cropping figure
%set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/alps_DEM_bar_events.jpg'],'Resolution',600)
  
% oriKGE vs DEM
figure
%set(gcf,'outerposition', [100 100 700 500]);
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.15 0.15 3.9 5.9];
nexttile;
boxplot(demgroups3',xts,'Symbol','k.')
%xlim([200 3200])
ylim([-2.5 1])
xlabel('DEM (m)')
ylabel('Original KGE')
view([90 -90])
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
%[Left Bottom Right Top] spacing
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);% when use this set the figure visible off for better cropping figure
%set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
%exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/alps_DEM_oriKGE.jpg'],'Resolution',600)

% opt KGE vs DEM
figure
%set(gcf,'outerposition', [100 100 700 500]);
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.15 0.15 3.9 5.9];
nexttile;
boxplot(demgroups4',xts,'Symbol','k.')
%xlim([200 3200])
ylim([-2.5 1])
xlabel('DEM (m)')
ylabel('Optimal KGE')
view([90 -90])
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
%[Left Bottom Right Top] spacing
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);% when use this set the figure visible off for better cropping figure
%set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
%exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/alps_DEM_optKGE.jpg'],'Resolution',600)

% combine into one box chart
trial1 = demgroups3';
trial2 = demgroups4';

% These grouping matrices label the columns:
grp1 = repmat(1:length(xts),size(trial1,1),1);
grp2 = repmat(1:length(xts),size(trial2,1),1);
% These color matrices label the matrix id:
clr1 = repmat(1,size(trial1));
clr2 = repmat(2,size(trial2));
% Combine the above matrices into one for x, y, and c:
x = [grp1;grp2];
y = [trial1;trial2];
c = [clr1;clr2];
% Convert those matrices to vectors:
x = x(:);
y = y(:);
c = c(:);
% Multiply x by 2 so that they're spread out:
x = x*2;
% Make the boxchart, 
figure
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.15 0.15 3.6 6.5];
nexttile;
b = boxchart(x(:),y(:),'GroupByColor',c(:),'MarkerStyle','.','JitterOutliers','on');
% Set the x ticks and labels, and add a legend
colormap(jet(10))
b(1).BoxFaceColor = [1 0 0];b(1).MarkerColor = [1 0 0];
b(2).BoxFaceColor = [0 0 1];b(2).MarkerColor = [0 0 1];
xlim([1 31])
xticks(2:2:2*length(xts));
xticklabels(xts)
xlabel('DEM (m)')
ylim([-3 1])
hax = gca;
hax.YAxis.MinorTickValues = linspace(-10,1,45);
hax.YMinorTick = 'on';
ylabel('KGE')
box on
view([90 -90])
set(gca,'FontName','Times New Roman','FontSize',12)
legend(["KGE original" "KGE optimal"],'Location','NorthOutside')
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
%[Left Bottom Right Top] spacing
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);% when use this set the figure visible off for better cropping figure
%set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/alps_DEM_KGE_all.jpg'],'Resolution',600)


%% read Alps raingauge measurement data

alps_rg = ['/shared/dondo/home/ml423/DA/alps_rg/Metadata/'];
lastall = ls([alps_rg '*.txt']);
args = strsplit(lastall);
p_allrg = [];

for j = 1:length(args)-1
    ffnm = args{j};
    rgslist=readtable([ffnm]);

    pdata = table2array(rgslist(:,9));
    pexist = [];
    for i = 1:length(pdata)
        tmp = pdata{i};
        if length(tmp)==0
            % P data does not exist for this station
        else
            pexist = [pexist;i];
        end
    end

    psum = rgslist(pexist,[1 3 4 10 11]);
    p_allrg = [p_allrg;psum];
end

edate = table2array(p_allrg(:,5));
latall = table2array(p_allrg(:,3));
lonall = table2array(p_allrg(:,2));
qualify_time_rg = find(edate>19730000);
rgs_loc = [];
% rg index, bid, rowind, colind,dem
for bid = [1:133 466:522]
    gid = WORLD_events(bid,1);
    gloc = find(gauge_events_first_last(:,1)==gid);
    r = gauge_events_first_last(gloc,6);
    c = gauge_events_first_last(gloc,7);
    area = gauge_events_first_last(gloc,5);
    lat = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/latlon/Basin' num2str(gid) '_lat.mat'],'lat_b');
    basin_lat = lat.lat_b;
    lon = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/latlon/Basin' num2str(gid) '_lon.mat'],'lon_b');
    basin_lon = lon.lon_b;
    tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_ind.mat'],'ind');
    ind = tmp.ind;
    blatf = basin_lat';
    blonf = basin_lon';
    blatf(ind) = nan;
    blonf(ind) = nan;
    basin_lat = blatf';
    basin_lon = blonf';

    latmin = min(basin_lat(:));
    lonmin = min(basin_lon(:));
    latmax = max(basin_lat(:));
    lonmax = max(basin_lon(:));
    for iq = 1:length(qualify_time_rg)
        tmplon = table2array(p_allrg(qualify_time_rg(iq),2));
        tmplat = table2array(p_allrg(qualify_time_rg(iq),3));
        if tmplon>lonmax||tmplon<lonmin||tmplat>latmax||tmplat<latmin
            % this rain station is outside basin box, not considered
        else
            for i = 1:r-2
                for j = 1:c-2
                    if tmplat<(basin_lat(i,j)+basin_lat(i+1,j))/2&&tmplat>(basin_lat(i+1,j)+basin_lat(i+2,j))/2&&tmplon>(basin_lon(i,j)+basin_lon(i,j+1))/2&&tmplon<(basin_lon(i,j+1)+basin_lon(i,j+2))/2
                        ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_fixed_900m/dem.bin'];
                        fid=fopen(ffnm,'rb','ieee-le');
                        tmpdata=fread(fid,inf,'single');
                        fclose(fid);
                        dem = reshape(tmpdata,c,r);
                        demf = dem';
                        ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_fixed_900m/countfile.out'];
                        fid=fopen(ffnm);
                        facc=fscanf(fid,'%d',[c,r]);
                        faccf = facc';
                        fclose(fid);


                        rgs_loc = [rgs_loc;[qualify_time_rg(iq),bid,i+1,j+1,demf(i+1,j+1),faccf(i+1,j+1),area]]; 
                    end
                end
            end
        end
    end
    bid
end

[Urg,ia,ic] = unique(rgs_loc(:,1));

% because some basins include smaller basins, therefore, one raingauge can
% be in multiple basins.

% reformat raingauge data based on rainfall events we have.
rgs_loc_sum = {};
for i = 1:length(rgs_loc)
    rgs_loc_sum{i,1} = rgs_loc(i,1);rgs_loc_sum{i,2} = rgs_loc(i,2);rgs_loc_sum{i,3} = rgs_loc(i,3);rgs_loc_sum{i,4} = rgs_loc(i,4);rgs_loc_sum{i,5} = rgs_loc(i,5);
    rgs_loc_sum{i,6} = rgs_loc(i,6);
    rgname = char(table2cell(p_allrg(rgs_loc(i),1)));
    rgnametable = p_allrg(rgs_loc(i),1);

    for j = 1:length(args)-1
        ffnm = args{j};
        rgslist=readtable([ffnm]);
        if ismember(rgnametable,rgslist(:,1))
            %this raingauge is in this folder/organization
            [~,foldernm,~] = fileparts(ffnm);
            foldernm = foldernm(1:end-5);
            if ismember(foldernm,['eHYD','Geosphere','MeteoSwiss','MeteoAM','ARPAV'])
                continue
                % eHYD and xxx folders have no daily P data
            end
            % if ~exist(['/shared/dondo/home/ml423/DA/alps_rg/' foldernm ],'dir')
            %     % meaning this folder does not exist in the uploaded data
            %     % on zenodo.
            %     continue
            % end
            tmprg = ls(['/shared/dondo/home/ml423/DA/alps_rg/' foldernm '/' rgname '*.txt']);
            rg_ts = readtable([tmprg(1:end-1)]);% the .txt extension has some issues here. 
            rg_dates = table2array(rg_ts(:,1));
            etot = length(find(WORLD_events(rgs_loc(i,2),:)>0))-1;
            for eve = 1:etot
                event = WORLD_events(rgs_loc(i,2),eve+1);
                yth = floor(event/10000);
                dth = mod(event,100);
                mth = floor((event-yth*10000)/100);
                thresdate = datetime(yth,mth,dth);
                if ~ismember(thresdate,rg_dates)
                    continue
                end


                durati = (thresdate-days(15)):days(1):(thresdate+days(15));
                [dury,durm,durd] = ymd(durati);
                ko = 0;
                rg_ts_values = [];
                for kk = 1:length(dury)
                    timetmp = datetime(dury(kk),durm(kk),durd(kk));
                    ko = ko+1;
                    if ismember(timetmp,rg_dates)
                        rg_ts_values(ko) = table2array(rg_ts(find(rg_dates==timetmp),13));
                    else
                        rg_ts_values(ko) = nan;
                    end
                end
                rg_ts_values = rg_ts_values(1:30);
                rgs_loc_sum{i,6+eve} = rg_ts_values;
            end

        end

    end
    disp(i)
end

%% [Figure4] raingauge verification Alps

% use eval_ind to get raingauge measurements for the IRC period
rain_compare = [];
% bid, eid, raingauge, rain00(basin avg), rainopt(basin avg), raingauge
% values, rain00(rg loc), rain34(rg loc), elevation, basin area
for bid = [1:129 132:133 466:504 507:522]

    ixx = find(~isnan(WORLD_events(bid,:)));
    etot = length(ixx)-1;
    if etot == 0
        continue
    end
    rgs_inbasin = find(rgs_loc(:,2)==bid);
    if isempty(rgs_inbasin)
        continue
        % no raingauges in current basin
    end
    if etot>146
        etot = 146;%this is to consider the raingauge observations in the alps, each basin not more than 146 events
    end
    for eid = 1:etot
        if isnan(event_quality(bid,eid+1))
            continue
            % event quality is no good per IRC policy
        end
        if isnan(eqcontrol(bid,eid))
            continue
        end
        for ii = 1:length(rgs_inbasin)
            rgv = rgs_loc_sum{rgs_inbasin(ii),eid+6};
            rowind = rgs_loc(rgs_inbasin(ii),3);
            colind = rgs_loc(rgs_inbasin(ii),4);
            elev = rgs_loc(rgs_inbasin(ii),5);
            area = rgs_loc(rgs_inbasin(ii),7);
            if isempty(rgv)
                continue
            end
            if isempty(rain00_sum{bid,eid})||isempty(rain34_sum{bid,eid})
                continue
            end

            ircts = eval_ind{bid,eid};
            % exact rg location pixel rainfall
            tmp = rain00_sum{bid,eid};
            tmp2 = tmp(colind,rowind,:);
            rgloc00 = sum(tmp2);
            tmp = rain34_sum{bid,eid};
            tmp2 = tmp(colind,rowind,:);
            rgloc34 = sum(tmp2);

            % exact rg with a neighbor 3x3, and get the nan median rainfall
            % tmpnn = rain00_sum{bid,eid};
            % tmpmm = rain34_sum{bid,eid};
            % allnn = [];
            % allmm = [];
            % for ni = -1:1
            %     for nj = -1:1
            %         tmpnn2 = tmpnn(colind+ni,rowind+nj,:);
            %         allnn = [allnn;sum(tmpnn2)];
            %         tmpmm2 = tmpmm(colind+ni,rowind+nj,:);
            %         allmm = [allmm;sum(tmpmm2)];
            %     end
            % end
            % rgloc00 = nanmedian(allnn);
            % rgloc34 = nanmedian(allmm);

            
            ircts = unique(round(ircts/24));% raingauge measurements is daily scale
            rgtot = sum(rgv(ircts));% raingauge event total rainfall mm
            rain_compare = [rain_compare;[bid,eid,rgs_loc(rgs_inbasin(ii),1),rain_orim(bid,eid),rain_optm(bid,eid),rgtot,rgloc00,rgloc34,kge_orim(bid,eid),kge_optm(bid,eid),elev,area,rgs_loc(rgs_inbasin(ii),6)]];

        end

    end
    disp(bid)

end


%% [Figure3&4] rain gauge point scale comparison, raingauge QPE watershed scale comparison, error diagnosis in the Alps
% characteristics of some rainfal thresholds, see common features perhaps
prob_bid = rain_compare(find(rain_compare(:,8)>370),1);
prob_bid(:,2) = rain_compare(find(rain_compare(:,8)>370),8);
prob_bids = unique(prob_bid(:,1));
for i = 1:length(prob_bids)
    tmp = find(prob_bid==prob_bids(i));
    prob_bids(i,2) = length(tmp);
    prob_bids(i,3) = length(find(WORLD_events(prob_bids(i,1),:)>0));
    prob_bids(i,4:200) = nan;
    prob_bids(i,4:3+length(tmp)) = prob_bid(find(prob_bid(:,1)==prob_bids(i)),2)';
end
figure
histogram(prob_bid(:,1))


interloc = find(rain_compare(:,1)==519&rain_compare(:,2)<62);

prob_bid = rain_compare(find(rain_compare(:,9)<0),1);
prob_bids = unique(prob_bid);
for i = 1:length(prob_bids)
    prob_bids(i,2) = length(find(prob_bid==prob_bids(i)));
    prob_bids(i,3) = length(find(WORLD_events(prob_bids(i,1),:)>0));
end

% find post IRC >800mm
prob_bid = rain_compare(find(rain_compare(:,8)>500),1);
prob_bid(:,2) = rain_compare(find(rain_compare(:,8)>500),8);
prob_bids = unique(prob_bid(:,1));
for i = 1:length(prob_bids)
    tmp = find(prob_bid==prob_bids(i));
    prob_bids(i,2) = length(tmp);
    prob_bids(i,3) = length(find(WORLD_events(prob_bids(i,1),:)>0));
    prob_bids(i,4:200) = nan;
    prob_bids(i,4:3+length(tmp)) = prob_bid(find(prob_bid(:,1)==prob_bids(i)),2)';
end
% rg >300

prob_bid = rain_compare(find(rain_compare(:,6)>500),1);
prob_bid(:,2) = rain_compare(find(rain_compare(:,6)>500),8);
prob_bids = unique(prob_bid(:,1));
for i = 1:size(prob_bids,1)
    tmp = find(prob_bid==prob_bids(i));
    prob_bids(i,2) = length(tmp);
    prob_bids(i,3) = length(find(WORLD_events(prob_bids(i,1),:)>0));
    prob_bids(i,4:200) = nan;
    prob_bids(i,4:3+length(tmp)) = prob_bid(find(prob_bid(:,1)==prob_bids(i)),2)';
end

figure
scatter(rain_compare(:,4),rain_compare(:,5),20,rain_compare(:,10),'filled')
clim([-2 1])
figure
scatter(rain_compare(:,6),rain_compare(:,7),'filled')
hold on
scatter(rain_compare(:,6),rain_compare(:,8),'filled')

figure
histogram(rain_compare(:,9),[-0.1:0.05:1])
hold on
histogram(rain_compare(:,10),[-0.1:0.05:1])


figure
scatter(rain_compare(:,6),rain_compare(:,4),'filled')
hold on
scatter(rain_compare(:,6),rain_compare(:,5),'filled')
hold on
plot([0 500],[0 500],'r')
xlim([0 700])
ylim([0 700])

nanjj = [];
for jj = 1:length(rain_compare)
    if isnan(rain_compare(jj,6))||isnan(rain_compare(jj,7))||isnan(rain_compare(jj,8))
        nanjj = [nanjj;jj];
    end
end


% plot rg point scale comparison

figure
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.15 0.15 4.8 4.2];
nexttile;
scatter(rain_compare(:,6),rain_compare(:,7),12,'ro','filled')
hold on
scatter(rain_compare(:,6),rain_compare(:,8),12,'bo','filled')
hold on
plot([0 650],[0 650],'k')
xlim([0 650])
ylim([0 650])
xlabel('RG (mm)')
ylabel('QPE (mm)')
box on
set(gca,'FontName','Times New Roman','FontSize',14)
exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/RG_comparison.jpg'],'Resolution',1000)



figure
cdfplot(rain_compare(:,6))
hold on
cdfplot(rain_compare(:,7))
hold on
cdfplot(rain_compare(:,8))
xlim([0 300])

rain_compare(find(rain_compare(:,8)>(rain_compare(:,6)+150)),1:2)


numband = 15;
demint = 200;
demgroupsrg = nan(numband,5000);% rg 
demgroupsrg2 = nan(numband,5000);% ERA5 ori at rg location
demgroupsrg3 = nan(numband,5000);% ERA5 Opt at rg location


for demband = 1:numband
    demloc = find(rain_compare(:,11)>=(200+(demband-1)*demint)&rain_compare(:,11)<(200+demband*demint));
    if ~isempty(demloc)
        demgroupsrg(demband,1:length(demloc)) = rain_compare(demloc,6)';
        demgroupsrg2(demband,1:length(demloc)) = rain_compare(demloc,7)';
        demgroupsrg3(demband,1:length(demloc)) = rain_compare(demloc,8)';

    end
end

% combine into one box chart
trial1 = demgroupsrg';
trial2 = demgroupsrg2';
trial3 = demgroupsrg3';

% These grouping matrices label the columns:
grp1 = repmat(1:length(xts),size(trial1,1),1);
grp2 = repmat(1:length(xts),size(trial2,1),1);
grp3 = repmat(1:length(xts),size(trial3,1),1);
% These color matrices label the matrix id:
clr1 = repmat(1,size(trial1));
clr2 = repmat(2,size(trial2));
clr3 = repmat(3,size(trial3));
% Combine the above matrices into one for x, y, and c:
x = [grp1;grp2;grp3];
y = [trial1;trial2;trial3];
c = [clr1;clr2;clr3];
% Convert those matrices to vectors:
x = x(:);
y = y(:);
c = c(:);
% Multiply x by 2 so that they're spread out:
x = x*2;
% Make the boxchart, 
figure
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.15 0.15 3.9 5.9];
nexttile;
b = boxchart(x(:),y(:),'GroupByColor',c(:),'MarkerStyle','.','JitterOutliers','on');
% Set the x ticks and labels, and add a legend
xlim([1 31])
xticks(2:2:2*length(xts));
xticklabels(xts)
xlabel('DEM (m)')
ylim([0 500])
hax = gca;
hax.YAxis.MinorTickValues = linspace(-10,1,45);
hax.YMinorTick = 'on';
ylabel('KGE')
box on
view([90 -90])
legend(["KGE original" "KGE optimal"],'Location','NorthOutside')
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
%[Left Bottom Right Top] spacing
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);% when use this set the figure visible off for better cropping figure
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

% exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/alps_DEM_KGE_all.jpg'],'Resolution',600)


tmp = rain_compare;
% divide into percentile and compare
rg25 = prctile(tmp(:,6),25);
rg50 = prctile(tmp(:,6),50);
rg75 = prctile(tmp(:,6),75);
rg90 = prctile(tmp(:,6),90);
rg95 = prctile(tmp(:,6),95);
rg99 = prctile(tmp(:,6),99);
rg100 = prctile(tmp(:,6),100);

locr1 = find(tmp(:,6)>=0&tmp(:,6)<rg25);
locr2 = find(tmp(:,6)>=0&tmp(:,6)<rg50);
locr3 = find(tmp(:,6)>=rg50&tmp(:,6)<rg75);
locr4 = find(tmp(:,6)>=rg75&tmp(:,6)<rg90);
locr5 = find(tmp(:,6)>=rg90&tmp(:,6)<rg95);
locr6 = find(tmp(:,6)>=rg95&tmp(:,6)<rg99);
locr7 = find(tmp(:,6)>=rg99&tmp(:,6)<rg100);

rgrain = [];
era5rain = [];
era5ircrain = [];
era5rainarea = [];
era5ircrainarea = [];
kgeori = [];
kgeopt = [];
for i = 2:7
    locrtmp = eval(['locr' num2str(i)]);
    rgrain(i-1) = nanmean(tmp(locrtmp,6));
    era5rain(i-1) = nanmean(tmp(locrtmp,7));
    era5ircrain(i-1) = nanmean(tmp(locrtmp,8));
    tmp1 = tmp(locrtmp,4);tmp1(tmp1<-1000) = [];
    tmp2 = tmp(locrtmp,5);tmp2(tmp2<-1000) = [];
    era5rainarea(i-1) = nanmean(tmp1);
    era5ircrainarea(i-1) = nanmean(tmp2);

    tmp1 = rain_compare(locrtmp,9);tmp1(tmp1<-1000) = [];
    tmp2 = rain_compare(locrtmp,10);tmp2(tmp2<-1000) = [];
    kgeori(i-1) = nanmean(tmp1);
    kgeopt(i-1) = nanmean(tmp2);

end
rgcall = [rgrain',era5rain',era5ircrain'];
kgeall = [kgeori',kgeopt'];

figure
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.25 0.25 4.95 4.2];
nexttile;
b = bar(rgcall);
b(1).FaceColor = [0 0 0 ];
b(2).FaceColor = [1 0 0 ];
b(3).FaceColor = [0 0 1];
xlabel('Extreme rainfall group')
ylabel('Group median rainfall (mm)')
box on
ylim([0 530])
set(gca,'XTick',1:6, 'XTickLabel',{'50%','75%','90%','95%','99%','100%'})
set(gca,'FontName','Times New Roman','FontSize',14)
legend(["RG" "ERA5L" "ERA5L-IRC"],'Location','Northwest')
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
%[Left Bottom Right Top] spacing
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);% when use this set the figure visible off for better cropping figure
%set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/alps_rg_groups.jpg'],'Resolution',600)




figure

t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.25 0.25 4.95 3.8];
nexttile;
b = bar(kgeall);
b(1).EdgeColor = [1 0 0];
b(2).EdgeColor = [0 0 1];
b(1).FaceColor = 'none';
b(2).FaceColor = 'none';
xlabel('Extreme rainfall group')
ylabel('Group median KGE')
box on
ylim([0 1])
hax = gca;
hax.YAxis.MinorTickValues = linspace(0,1,11);
hax.YMinorTick = 'on';
set(gca,'XTick',1:6, 'XTickLabel',{'50%','75%','90%','95%','99%','100%'})
set(gca,'FontName','Times New Roman','FontSize',14)
legend(["ERA5L" "ERA5L-IRC"],'Location','Northwest')
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
%[Left Bottom Right Top] spacing
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);% when use this set the figure visible off for better cropping figure
%set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/alps_rg_groups_kge.jpg'],'Resolution',600)



figure
scatter(rain_compare(:,6),rain_compare(:,8),10,rain_compare(:,11),'filled')
hold on
plot([0 500],[0 500],'r')
xlim([0 500])
ylim([0 500])

% area scale qpe compare
rori = [];
ropt = [];
rref = [];
for bid = [1:129 132 133 466:522]
    disp(bid)
    for eid = 1:149
        if ~isnan(WORLD_events(bid,1+eid))
            if ~isnan(eqcontrol(bid,eid))
            if ~isnan(rain_orim(bid,eid))
                if ~isnan(apgd_mean(bid,eid))
                %if ~ismember([bid,eid],reprob,'rows')
                    rori = [rori;rain_orim(bid,eid)];
                    ropt = [ropt;rain_optm(bid,eid)];
                    rref = [rref;apgd_mean(bid,eid)];

                %end
                end
            end
            end
        end
    end
end
figure
scatter(rori,ropt,'filled')
hold on
plot([0 600],[0 600],'r')
hold on
scatter(rori,rref,'filled')


% plot RG QPE watershed scale comparison

figure

t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.15 0.15 4.8 4.2];
nexttile;
scatter(rref,rori,12,'ro','filled')
hold on
scatter(rref,ropt,12,'bo','filled')
hold on
plot([0 1000],[0 1000],'K')

xlim([0 1050])
ylim([0 1050])
xlabel('APGD QPE (mm)')
ylabel('QPE (mm)')
box on
set(gca,'FontName','Times New Roman','FontSize',14)
exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/RG_QPE_comparison.jpg'],'Resolution',1000)



figure
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.25 0.25 4.95 4.2];
nexttile;
cd1 = cdfplot(rref);set(cd1,'Color','k','LineWidth',2.5);
hold on
cd2 = cdfplot(rori);set(cd2,'Color','r','LineWidth',2);
hold on
cd3 = cdfplot(ropt);set(cd3,'Color','b','LineWidth',2);
xlim([0 450])
xlabel('Event total rainfall (mm)')
ylabel('CDF')
title(' ')
box on
legend(["APGD" "ERA5L" "ERA5L-IRC"],'Location','Northwest')
set(gca,'FontName','Times New Roman','FontSize',14)
exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/RG_QPE_CDF_comparison.jpg'],'Resolution',600)


figure
t = tiledlayout(1,1,'Padding','none');
t.Units = 'inches';
t.OuterPosition = [0.25 0.25 4.95 3.6];
nexttile;
cd1 = cdfplot(kge_orim(:));set(cd1,'Color','r','LineWidth',2);
hold on
cd2 = cdfplot(kge_optm(:));set(cd2,'Color','b','LineWidth',2);
hold off
xlim([-3 1])
ylim([0 1.02])
xlabel('KGE')
ylabel('CDF')
title(' ')
box on
legend(["ERA5L" "ERA5L-IRC"],'Location','Northwest')
set(gca,'FontName','Times New Roman','FontSize',14)
exportgraphics(gcf,['/shared/dondo/home/ml423/world_mts_runs/figs/synthesis/KGE_CDF_comparison.jpg'],'Resolution',600)




[bidarea,eidarea] = find(rain_optm>(rain_orim+150));
tmp = [bidarea,eidarea];
for i = 1:length(tmp)
    tmp(i,3) = rain_optm(bidarea(i),eidarea(i))-rain_orim(bidarea(i),eidarea(i));
end
beprob150new = sortrows(tmp);

eliminated_cases = setdiff(beprob_old,beprob,'rows');

% 519 2, and 521 10 are good examples of Era5L underestimating
for i = 1:length(beprob150new)
bb = beprob150new(i,1);
ee = beprob150new(i,2);
if ~isempty(eidsum{bb})
if ismember(ee,eidsum{bb})
disp(i)
end
end
end


[bidarea,eidarea] = find(rain_optm>(apgd_mean+200));
tmp = [bidarea,eidarea];
for i = 1:length(tmp)
    tmp(i,3) = rain_optm(bidarea(i),eidarea(i))-apgd_mean(bidarea(i),eidarea(i));
    tmp(i,4) = rain_orim(bidarea(i),eidarea(i));
    tmp(i,5) = rain_optm(bidarea(i),eidarea(i));
    tmp(i,6) = apgd_mean(bidarea(i),eidarea(i));
end
beprob200apgd = sortrows(tmp);




rg50 = prctile(rref,50);
rg75 = prctile(rref,75);
rg90 = prctile(rref,90);
rg95 = prctile(rref,95);
rg99 = prctile(rref,99);
rg100 = prctile(rref,100);

locr2 = find(rref>=0&rref<rg50);
locr3 = find(rref>=rg50&rref<rg75);
locr4 = find(rref>=rg75&rref<rg90);
locr5 = find(rref>=rg90&rref<rg95);
locr6 = find(rref>=rg95&rref<rg99);
locr7 = find(rref>=rg99&rref<rg100);
rgrainqpe = [];
era5rainqpe = [];
era5ircrainqpe = [];

for i = 2:7
    locrtmp = eval(['locr' num2str(i)]);
    rgrainqpe(i-1) = nanmean(rref(locrtmp));
    era5rainqpe(i-1) = nanmean(rori(locrtmp));
    era5ircrainqpe(i-1) = nanmean(ropt(locrtmp));
end
rgcallqpe = [rgrainqpe',era5rainqpe',era5ircrainqpe'];

figure
bar(rgcallqpe)








