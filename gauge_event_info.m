% check gauge and events here
tooshort = 0;
partialdone = 0;

fid = fopen('basin_code.txt','r');
fom = '%d';
bid = fscanf(fid,fom);
fclose(fid);

tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/world_selected_gauge_events_allinfo_v3.mat'],'gauge_events_first_last');
gauge_all_info = tmp.gauge_events_first_last;

tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/WORLD_events_v3.mat'],'WORLD_events');
WORLD_events = tmp.WORLD_events;
tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/WORLD_events_qc_v3.mat'],'event_quality');
event_quality = tmp.event_quality;
gid = WORLD_events(bid,1);

tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/world_selected_gauge_events_allinfo.mat'],'gauge_events_first_last');
gauge_all_info_old = tmp.gauge_events_first_last;
gtmp = find(gauge_all_info_old(:,1)==gid);
if isempty(gtmp)
    original_last = 999999;
else
    original_last = gauge_all_info_old(gtmp,3);
end


goodevents = find(event_quality(bid,:)<2);
if length(goodevents)<=1
    % only 0 or 1 event, this basin should be abandoned.
    tooshort=1;
    fid = fopen('too_short.txt','w');
    fprintf(fid,'%d',tooshort);
    fclose(fid);
    quit
else
    eventslist = WORLD_events(bid,goodevents);
end


tmp = load(['/shared/dondo/home/ml423/world_mts/obs_events/obs_' num2str(gid) '.mat'],'obs');
obsgid = tmp.obs; clear tmp

gloc = find(gauge_all_info(:,1)==gid);
area = gauge_all_info(gloc,5);

% if area is less than 1000, then just do the all available flow records
% for spinup
if area<1000
    firstevent = gauge_all_info(gloc,2);
    lastevent = gauge_all_info(gloc,3);
else
    % if area is too big, then just cover the interested events and do the
    % spinup, otherwise the spinup takes too much time.
    firstevent = min(eventslist);
    lastevent = max(eventslist);
end


% this study is 50 year study 1974 to 2024
% comparing 19730911
if firstevent<19730911
    final1 = 19730911;
else
    final1 = firstevent; 
end

obsloc = find(obsgid(:,1)==final1);
% find the first non-999 record after firstevent, and this should rather be
% the firstdate to do the spinup. 
for obsloop = obsloc:length(obsgid)
    if obsgid(obsloop,2)>0
        final1 = obsgid(obsloop,1);
        break
    end
end

fyear = floor(final1/10000);fmonth = floor((final1-fyear*10000)/100);fday=mod(final1,100);
startdate = datetime(fyear,fmonth,fday)-days(10);[syear,smonth,sday] = ymd(startdate);
lyear = floor(lastevent/10000);lmonth = floor((lastevent-lyear*10000)/100);lday=mod(lastevent,100);
enddate = datetime(lyear,lmonth,lday)+days(10);[eyear,emonth,eday] = ymd(enddate);
% simulation ends 10 days after the first event

fid = fopen('spinsufix.txt','r');
fom = '%s';
spinsuffix = fscanf(fid,fom);
fclose(fid);

ffnm = ['/shared/dondo/home/ml423/world_mts_outputs/Basin' num2str(gid) '_output_900m1hr_' spinsuffix '_spinup/' num2str(lastevent) '/totalrunoff'];
ffnmoriginallast = ['/shared/dondo/home/ml423/world_mts_outputs/Basin' num2str(gid) '_output_900m1hr_' spinsuffix '_spinup/' num2str(original_last) '/totalrunoff'];

if isfile(ffnm)
    tooshort = 1;% meaning the spinup outputs already completely done
else
    if isfile(ffnmoriginallast)
        partialdone = 1;% partially done from version1 events
    end
end

avaidays = datetime(fyear,fmonth,fday):days(1):datetime(lyear,lmonth,lday);

if length(avaidays)<365
    tooshort = 1;
end

if area>2000% when area > 4000km2, discard 2000.20250613 update
    tooshort = 1;
end

fid = fopen('too_short.txt','w');
fprintf(fid,'%d',tooshort);
fclose(fid);

spinstartyear = syear;
spinstartmonth = smonth;
spinstartday = sday;

sintyear = syear+1;
sintmonth = 4;
sintday = 30;

if partialdone == 0
    sintyear2 = syear+1;
    sintmonth2 = 5;
    sintday2 = 1;
    startfromwarmup = 1;
else
    sintyear2 = floor(original_last/10000);
    sintmonth2 = floor((original_last-sintyear2*10000)/100);
    sintday2 = mod(original_last,100);
    startfromwarmup = 0;
end

spinendyear = eyear;
spinendmonth = emonth;
spinendday = eday;

fid = fopen('sfy.txt','w');
fprintf(fid,'%d',spinstartyear);
fclose(fid);
fid = fopen('sfm.txt','w');
fprintf(fid,'%d',spinstartmonth);
fclose(fid);
fid = fopen('sfd.txt','w');
fprintf(fid,'%d',spinstartday);
fclose(fid);
fid = fopen('sly.txt','w');
fprintf(fid,'%d',sintyear);
fclose(fid);
fid = fopen('slm.txt','w');
fprintf(fid,'%d',sintmonth);
fclose(fid);
fid = fopen('sld.txt','w');
fprintf(fid,'%d',sintday);
fclose(fid);

fid = fopen('sfy2.txt','w');
fprintf(fid,'%d',sintyear2);
fclose(fid);
fid = fopen('sfm2.txt','w');
fprintf(fid,'%d',sintmonth2);
fclose(fid);
fid = fopen('sfd2.txt','w');
fprintf(fid,'%d',sintday2);
fclose(fid);
fid = fopen('sly2.txt','w');
fprintf(fid,'%d',spinendyear);
fclose(fid);
fid = fopen('slm2.txt','w');
fprintf(fid,'%d',spinendmonth);
fclose(fid);
fid = fopen('sld2.txt','w');
fprintf(fid,'%d',spinendday);
fclose(fid);

fid = fopen('gauge_number.txt','w');
fprintf(fid,'%d',gid);
fclose(fid);
fid = fopen('startfrom0.txt','w');
fprintf(fid,'%d',startfromwarmup);
fclose(fid);


quit


