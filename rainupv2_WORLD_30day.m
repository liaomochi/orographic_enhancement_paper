%rainf
% This script is designed to check if the plan is working
fid = fopen('current_basin.txt','r');
fom = '%d';
basinid = fscanf(fid,fom);
fclose(fid);

fid = fopen('intervalmulti.txt','r');
fom = '%d';
rainmulti = fscanf(fid,fom);
fclose(fid);

fid = fopen('utchh.txt','r');
fom = '%d';
utch = fscanf(fid,fom);
fclose(fid);

tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/world_selected_gauge_events_allinfo_v3.mat'],'gauge_events_first_last');
Basins_sum = tmp.gauge_events_first_last;clear tmp;
tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/WORLD_events_v3.mat'],'WORLD_events');
WORLD_events = tmp.WORLD_events;

gid = WORLD_events(basinid,1);

gloc = find(Basins_sum(:,1)==gid);


tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_ind.mat'],'ind');
ind = tmp.ind;
tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_intt.mat'],'intt');
intt = tmp.intt;

tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/distance/Basin' num2str(gid) 'dis_stream.mat'],'dis_stream');
d_stream = tmp.dis_stream;
tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/distance/Basin' num2str(gid) 'dis_outlet.mat'],'dis_outlet');
d_outlet = tmp.dis_outlet;        


fid = fopen('datnm.txt','r');
fom = '%d';
event = fscanf(fid,fom);
fclose(fid);


nmm = num2str(event);
tww = 3;

fid = fopen('bl.txt','r');
fom = '%d';
p = fscanf(fid,fom);
fclose(fid);

fid = fopen('pre_rise_pt.txt','r');
fom = '%d';
rains_bd = fscanf(fid,fom);
fclose(fid);

fid = fopen('filenum.txt','r');
fom = '%d';
filnm = fscanf(fid,fom);
fclose(fid);

fid = fopen('filu.txt','r');
fom = '%d';
filnmc = fscanf(fid,fom);
fclose(fid);

wron = load(['/shared/dondo/home/ml423/world_mts_runs/IRC1_sup/wrontest5MMORIw' num2str(ceil(p/tww)) num2str(mod(filnm,10)) 'cc' num2str(filnm) '_' num2str(gid) '_' num2str(utch,'%2.2d') '.mat']);
wron = wron.wron;

% pbb2 = load(['/shared/dondo/home/ml423/world_mts_runs/IRC1_sup/pbb2test5MMORIw' num2str(ceil(p/tww)) num2str(mod(filnm,10)) 'cc' num2str(filnm) '_' num2str(gid) '_' num2str(utch,'%2.2d') '.mat'],'pbb2');
% pbb2 = pbb2.pbb2;


workingwindow = 3;
modd = 2;

fid = fopen('loop.txt','r');% for world basins loop is not 0-96 for 15 minutes in a day, but 0-90 for 8hours in 30 day
fom = '%d';
ck = fscanf(fid,fom);
fclose(fid);

if ck==1
    ck=2;% border adjust
end

    pao=1;


% strr = ['corrtest5MMORIw' num2str(ceil(p/tww)) num2str(mod(filnm,10)) 'cc' nmm 'narrowwide'];
% corp = load(['/shared/dondo/home/ml423/HMIOP/NoDAre/New era/' strr '.mat']);
% corr = corp.corr;


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

% now obs is 8-hourly resolution for 30 day period.


%ounm = {'Precipitation_IMERGEnewb111'};
%ounm = {'Precipitation_St4dbknewb111'};
 ounm = {'Precipitation_ERA5_land_oriunit'};

rainfallb1 = load('rainfallb1.mat');
rainfallb1 = rainfallb1.rainfallb1;

raints = [];
for i = 1:720
    tmp = rainfallb1{i};
    raints(i) = mean(tmp(intt));
end
rsum = sum(raints);


fid = fopen('loop2.txt','r');
fom = '%d';
itr = fscanf(fid,fom);
fclose(fid);

% fid = fopen('loop3.txt','r');
% fom = '%d';
% lp3 = fscanf(fid,fom);
% fclose(fid);



xtmpbb1 = load(['xtmpbb1.mat']);
xtmpbb1 = xtmpbb1.xtmpbb1;
% mode 1 all 10 replicates
% obs2mat = load('obs2mat.mat');
% obs2mat = obs2mat.obs2mat;
% mode 2 coarse hydrograph
% obsdata = load('b2agg10.mat');
% obsdata = obsdata.obsdatav4;
% original observation


% obsnow = obs2mat(:,lp3);
% disp(['$$$$$$$$$$$ loop3 value is ' num2str(lp3)])


ak = 4+(ck-1)*8;%maximum is 720 for ak, and 90 for ck, ck is ranging from 1-90 for 8-hourly in 30 day

ffnm=['/shared/dondo/home/ml423/world_mts_outputs/Basin' num2str(gid) '_output_900m1hr/' nmm '_SW/' nmm '/streamflowarea'];
fid=fopen(ffnm,'rb','ieee-le');
tmpdata=fread(fid,inf,'single');
NODAc2o1=reshape(tmpdata,Basins_sum(gloc,7),Basins_sum(gloc,6),720);
fclose(fid);

for k = 1:720
    flownoda3(k) = NODAc2o1(Basins_sum(gloc,9),Basins_sum(gloc,8),k);
end

no3 = flownoda3;
simu = no3(1:8:720);

percx = (simu(ck)-obsnow(ck))./obsnow(ck);
perc = abs(simu(ck)-obsnow(ck))./obsnow(ck);

dif11 = simu(ck-1)-obsnow(ck-1);
dif22 = simu(ck)-obsnow(ck);

difv = (dif11+dif22)*28800/2;% 8 hour is 28800 seconds
disp(['$$$ water diff after sm change is ' num2str(difv)])

tmp = 0;
tmp1 = 0;

tmp2 = mean(obsnow);
for i = 1:ck
    tmp = tmp + (simu(i)-obsnow(i))^2;
    tmp1 = tmp1 + (tmp2-obsnow(i))^2;
end
NSE = 1-tmp./tmp1;

tmp = 0;
tmp1 = 0;
tmp2 = mean(obsnow);
for i = 1:length(simu)
    tmp = tmp + (simu(i)-obsnow(i))^2;
    tmp1 = tmp1 + (tmp2-obsnow(i))^2;
end
NSEall = 1-tmp./tmp1;

pertmp = perc;
fid = fopen(['pernn' num2str(itr) '.txt'],'w');
fprintf(fid,'%f',pertmp);
fclose(fid);

if itr>=1
    fid = fopen('qualflag.txt','r');
    fom = '%d';
    qf = fscanf(fid,fom);
    fclose(fid);
end

if itr>=1&&qf==100
    
    percall = [];
    for pkj = 0:itr
        fid = fopen(['pernn' num2str(pkj) '.txt'],'r');
        fom = '%f';
        percxt = fscanf(fid,fom);
        fclose(fid);
        percall = [percall;percxt];
    end
    disp(percall)
    mino = percall(end-1)-percall(end);
    if abs(mino)<0.01
        stopit = 1;
        fid = fopen('stopit.txt','w');
        fprintf(fid,'%d',stopit);
        fclose(fid);
    else
        stopit = 0;
        fid = fopen('stopit.txt','w');
        fprintf(fid,'%d',stopit);
        fclose(fid);
        
    end
    
    %         if NSE>0.99999
    %             stopit = 1;
    %             fid = fopen('stopit.txt','w');
    %             fprintf(fid,'%d',stopit);
    %             fclose(fid);
    %         end
    
end



%if perc < 0.25
%    bdy = min(ak,24);
%else
bdy = min((ak-24),250);%24 is 2:00am, rainfall starts at 2:00am
% change to 20240521
bdy = min((ak-12),250);
% change to 20240805
bdy = ak-2;

% change below 11/02/2024
bdy = ak-rains_bd;
if bdy<8
    bdy=8;
end
%donechange
disp(['bdy  rain_starts and current modifying time are are' num2str(bdy) ' - ' num2str(rains_bd) ' - ' num2str(ak)])
disp(['Event is ' num2str(event)])
cdy = 2;
cdy = 1;%11/06/2024

if modd == 2

    dif1 = simu(ck-1)-obsnow(ck-1);
    dif2 = simu(ck)-obsnow(ck);

    if obsnow(ck)<0||isnan(obsnow(ck))
        quit
    end

    % rdif3 = 0.5*(2*dif1+1/3*(dif2-dif1))*9600;
    % rdif2 = 0.5*(dif1+dif2)*9600;
    % rdif1 = 0.5*(2*dif2+1/3*(dif1-dif2))*9600;
    %

    % rdif = [rdif1;rdif2;rdif3];
    rdif = rainmulti*(dif1+dif2)*8*3600/2;%cubic meters of water difference for each ck tick, which is 8 hours in this study.

    disp(['$$$$$$$$$$$ ck value is ' num2str(ck)])
    disp(['$$$$$$$$$$$ simu value is ' num2str(simu(ck))])
    disp(['$$$$$$$$$$$ percx value is ' num2str(percx)])
    disp(['$$$$$$$$$$$ water diff value is ' num2str(rdif)])

    %end
    % mod 1
    %rainshare = rdif.*(pbb2(1:bdy)./sum(pbb2(1:bdy)));
    % mod 2
    %rainshare = rdif.*ones(bdy,1)./bdy;
    % one can consider to add iterative processes
    if perc > 0.05

        rainshare = rdif*ones((bdy-cdy+1),1)./(bdy-cdy+1);%rainshare has to share again across timesteps between each ck tick
        %rainshare = rdif(part).*(pbb2(cdy:bdy)./sum(pbb2(cdy:bdy)));
        %rainshare = rdif.*[4/5*ones(1,8-cdy+1)./(8-cdy+1),1/5*ones(1,bdy-8)./(bdy-8)];
        disp(['============== NSE upto now is ' num2str(NSE)])
        disp(['============== NSE  entire  is ' num2str(NSEall)])
        rainshare_sub = rainshare/8/rainmulti;%/8;%3640
        if ak<=4
            inckmin = 0;
        else
            inckmin = -4*rainmulti;
        end

        if ak>=716
            inckmax = 0;
        else
            inckmax = 4*rainmulti;
        end

        for interck = inckmin:inckmax% in total 8 hours around the ck tick

            for jh = cdy:bdy
                if (ak-(jh-1)+interck)<=(rains_bd-20)

                    disp('corrections toooo early')
                    disp(['ak,ck,jh,interck, rainstart are ' num2str(ak) ' ' num2str(ck) ' ' num2str(jh) ' ' num2str(interck) ' ' num2str(rains_bd)])
                end

                if (ak-(jh-1)+interck)<=0
                    break
                end
                tmpp1 = find(wron{ak+interck,ak-(jh-1)+interck}(:,:)>0);
                [l1,l2] = find(wron{ak+interck,ak-(jh-1)+interck}(:,:)>0);

                tmpp2 = find(rainfallb1{ak-(jh-1)+interck}(tmpp1)>0);
                gro = length(tmpp2);

                % add new below for inter step interpolation of particle
                % 11/06 2024
                % locations
                tmpp1up = find(wron{ak+interck,ak-(jh-2)+interck}(:,:)>0);

                %disp(['$$$$$$$$$$$ total pixel is ' num2str(gro) ' while ' num2str(length(tmpp1up)) ' is with interpolation'])

                if length(tmpp1up)>0&&length(tmpp1)>0
                    d_str_max = max(d_stream(tmpp1up));
                    d_str_min = min(d_stream(tmpp1));

                    d_out_max = max(d_outlet(tmpp1up));
                    d_out_min = min(d_outlet(tmpp1));
                    d_out_median = median(d_outlet([tmpp1up;tmpp1]));


                    if d_out_max<=6
                        continue
                    else
                        d_out_min = 5;
                    end
%                    tmpp1modified = find(d_outlet>=d_out_min&d_outlet<=d_out_max&d_stream>=d_str_min&d_stream<=d_str_max);
                    tmpp1modified = find(d_outlet>=d_out_median&d_outlet<=d_out_max&d_stream==d_str_max);%3000
                    if length(tmpp1modified)<=1%3100 meanning tracked pixel too limited! have to think
                        continue
                    end
                    [l1mod,l2mod] = ind2sub([Basins_sum(gloc,7),Basins_sum(gloc,6)],tmpp1modified);
                    tmpp2modified = find(rainfallb1{ak-(jh-1)+interck}(tmpp1modified)>0);
                    gromod = length(tmpp2modified);
                    % the 4 variables above should replace previous vars.
                    tmpp1 = tmpp1modified;
                    tmpp2 = tmpp2modified;
                    l1 = l1mod;
                    l2 = l2mod;
                    gro = gromod;
                %end
                % finished new changes % 2700 above

             
                %disp(['$$$$$$$$$$$ total pixel is ' num2str(gro) ' while ' num2str(gromod) ' is with interpolation'])

                sumrain = sum(rainfallb1{ak-(jh-1)+interck}(tmpp1(tmpp2)));
                if sumrain>0
                    doublecorr = zeros(Basins_sum(gloc,7),Basins_sum(gloc,6));
                    % if saturonly <= 3
                    for ju = 1:gro
                        correction = rainshare_sub(jh-cdy+1)/900/900*1000*rainfallb1{ak-(jh-1)+interck}(tmpp1(tmpp2(ju)))./sumrain;
                        %correction = rainshare_sub(jh-cdy+1)/900/900*1000/gro;%IRC 2900
                        % if correction>2
                        %     correction=2;
                        % end
                        % if correction<-2
                        %     correction=-2;%3200
                        % end


                        % if ju==1
                        %     disp(['$$$$$$$$$$$ each pixel correction is ' num2str(correction) ' mm'])
                        % end

                        for kl = -1:1
                            for kl2 = -1:1
                                working_row = l1(tmpp2(ju))+kl;
                                working_col = l2(tmpp2(ju))+kl2;
                                
                                if working_row>0&&working_row<=Basins_sum(gloc,7)&&working_col>0&&working_col<=Basins_sum(gloc,6)

                                    if doublecorr(working_row,working_col)==1

                                    else
                                        xtmpbb1{ak-(jh-1)+interck}(working_row,working_col) = xtmpbb1{ak-(jh-1)+interck}(working_row,working_col)-correction;%2.25
                                        doublecorr(working_row,working_col)=1;
                                    end
                                %else
                                %    disp(['IRC current pixel outside of boundary cxr ' num2str(working_row) ' ' num2str(working_col) ' w kl kl2 ' num2str(kl) ' ' num2str(kl2)])
                                end
                            end
                        end
                    end
                end


                end
                % above end is 2800

                xtmpbb1{ak-(jh-1)+interck}(find(xtmpbb1{ak-(jh-1)+interck}(:,:)<0)) = 0;

            end

        end


        pasn = 0;
        fid = fopen('pasn.txt','w');
        fprintf(fid,'%d',pasn);
        fclose(fid);


    else
        pasn = 1;
        fid = fopen('pasn.txt','w');
        fprintf(fid,'%d',pasn);
        fclose(fid);

    end



    save(['/shared/dondo/home/ml423/world_mts_runs/IRC1/Basin' num2str(gid) 'outputs/refxtmp' num2str(filnmc) '/no3_' num2str(ck) '_' num2str(itr)],'no3')
    save(['xtmpbb1.mat'],'xtmpbb1')
    save(['/shared/dondo/home/ml423/world_mts_runs/IRC1/Basin' num2str(gid) 'outputs/refxtmp' num2str(filnmc) '/xtmpbb1_' num2str(ck) '_' num2str(itr)],'xtmpbb1')
    %    save(['/shared/dondo/home/ml423/HMIOP/NoDAre/New era/' strr '.mat'],'corr')%
    %    commented out 20240420 no use

end


for i = 1  %:46

    out_dir = ['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_input_ERA5_900m_tmp/' nmm '/'];
    dataf = [];
    for j = 1:720
        datas = xtmpbb1{j+(i-1)*720}; % can add stuff
        dataf(:,:,j) = datas;
    end
    dataf(dataf<0) = 0;
    fnm = [out_dir,ounm{1}];
    raindata = reshape(dataf,[Basins_sum(gloc,7),Basins_sum(gloc,6),720]);
    %    mrains = [mrains;mean(raindata(:))];
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

clear all;
%

quit