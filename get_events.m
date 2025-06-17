%rainf
% This script is designed to check if the plan is working

fid = fopen('basin_code.txt','r');
fom = '%d';
bid = fscanf(fid,fom);
fclose(fid);

fid = fopen('event_code.txt','r');
fom = '%d';
eid = fscanf(fid,fom);
fclose(fid);

tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/WORLD_events_v3.mat'],'WORLD_events');
WORLD_events = tmp.WORLD_events;

tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/WORLD_events_qc_v3.mat'],'event_quality');
event_quality = tmp.event_quality;
tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/WORLD_events_qc_control.mat'],'eqcontrol');
eqcontrol = tmp.eqcontrol;

evnt = WORLD_events(bid,1+eid);gid = WORLD_events(bid,1);
ymd2y = floor(evnt/10000);
ymd2d = mod(evnt,100);
ymd2m = floor((evnt-ymd2y*10000)/100);
icdate = datetime(ymd2y,ymd2m,ymd2d)-days(16);% careful here! the intial condition is the file from 16 days ago not 15
[icy,icm,icd] = ymd(icdate);
ic = icy*10000+icm*100+icd;

fid = fopen('datnm.txt','w');
fprintf(fid,'%d',evnt);
fclose(fid);
fid = fopen('icdaten.txt','w');
fprintf(fid,'%d',ic);
fclose(fid);



if isnan(event_quality(bid,1+eid))
    eq = 0;
else
    eq = 1;
end



IRC=4400;
IRC=4140;

flnm = evnt*100-IRC*100000+24;
ffnm = ['/shared/dondo/home/ml423/world_mts_runs/IRC1/Basin' num2str(gid) 'outputs/refxtmp' num2str(flnm) 'fullresults/' num2str(evnt) '_SW/' num2str(evnt) '/totalrunoff'];
if isfile(ffnm)
    eq = 0;% meaning the IRC outputs already completely done
    eq = 1; %force 2025 05 16
end

if isnan(eqcontrol(bid,eid))
    eq = 0;
end

fid = fopen('event_qc.txt','w');
fprintf(fid,'%d',eq);
fclose(fid);


quit