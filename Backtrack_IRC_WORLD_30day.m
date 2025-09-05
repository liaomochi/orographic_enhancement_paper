
% load current basin, event date
fid = fopen('current_basin.txt','r');
fom = '%d';
basinid = fscanf(fid,fom);
fclose(fid);

fid = fopen('utchh.txt','r');
fom = '%d';
utch = fscanf(fid,fom);
fclose(fid);

fid = fopen('datnm.txt','r');
fom = '%d';
event = fscanf(fid,fom);
fclose(fid);


tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/world_selected_gauge_events_allinfo_v3.mat'],'gauge_events_first_last');
Basins_sum = tmp.gauge_events_first_last;clear tmp;

tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/WORLD_events_v3.mat'],'WORLD_events');
WORLD_events = tmp.WORLD_events;

gid = WORLD_events(basinid,1);

gloc = find(Basins_sum(:,1)==gid);

area = Basins_sum(gloc,5);

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
% GRDC observation interpolation
obs2 = interp1(1:31,obsnow,[1:0.33333:31]);
obsnow = obs2(1:90);

% load basin boundary
tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_ind.mat'],'ind');
ind = tmp.ind;
tmp = load(['/shared/dondo/home/ml423/DA/Critical matrices/World_Basins/inttind/Basin' num2str(gid) '_intt.mat'],'intt');
intt = tmp.intt;

nmm = num2str(event);
pmm = num2str(event);
rst = 2;%1. original, 2. dbk

fid = fopen('bl.txt','r');
fom = '%d';
p = fscanf(fid,fom);
fclose(fid);

tww = 3;


% p  = total windows corrected , = 2 IF W12 AND W13, = 3 IF W12,W13,W14

% windows are to be changed using the shell script
w2bd1 = 17;
w2bd = 28;
w3bd1 = 28;
w3bd = 43;
w4bd1 = 43;
w4bd = 59;

fid = fopen('inflection_pt.txt','w');
fprintf(fid,'%d',w4bd);
fclose(fid);
fid = fopen('rising_pt.txt','w');
fprintf(fid,'%d',w2bd);
fclose(fid);


% classic 4-window IRC boundary settings
if p==0%window 2 for p==0
    fid = fopen('filenum.txt','r');
    fom = '%d';
    filnm = fscanf(fid,fom);
    fclose(fid);
    
    filnm = filnm+2;
    disp(filnm)
    bd = w2bd;
    fid = fopen('bd.txt','w');
    fprintf(fid,'%d',bd);
    fclose(fid);
    
    bd1 = w2bd1;
    fid = fopen('bd1.txt','w');
    fprintf(fid,'%d',bd1);
    fclose(fid);
    
else
    
    if mod(p,tww)==1%window 3
        fid = fopen('filenum.txt','r');
        fom = '%d';
        filnm = fscanf(fid,fom);
        fclose(fid);
        
        filnm = filnm+1;
        
        bd = w3bd;
        fid = fopen('bd.txt','w');
        fprintf(fid,'%d',bd);
        fclose(fid);
        
        bd1 = w3bd1;
        fid = fopen('bd1.txt','w');
        fprintf(fid,'%d',bd1);
        fclose(fid);
       
    elseif mod(p,tww)==2%window 4
        fid = fopen('filenum.txt','r');
        fom = '%d';
        filnm = fscanf(fid,fom);
        fclose(fid);
        
        filnm = filnm+1;
        
        bd = w4bd;
        fid = fopen('bd.txt','w');
        fprintf(fid,'%d',bd);
        fclose(fid);
        
        bd1 = w4bd1;
        fid = fopen('bd1.txt','w');
        fprintf(fid,'%d',bd1);
        fclose(fid);
        
    elseif mod(p,tww)==0%window2
        fid = fopen('filenum.txt','r');
        fom = '%d';
        filnm = fscanf(fid,fom);
        fclose(fid);
        
        filnm = filnm+8;
        
        bd = w2bd;
        fid = fopen('bd.txt','w');
        fprintf(fid,'%d',bd);
        fclose(fid);
        
        bd1 = w2bd1;
        fid = fopen('bd1.txt','w');
        fprintf(fid,'%d',bd1);
        fclose(fid);
        
    end
end

fid = fopen('filu.txt','w');
fprintf(fid,'%d',filnm);
fclose(fid);

%=========================== above is name changing
%parpool('local',12)




fid = fopen('filenum.txt','r');
fom = '%d';
filnm = fscanf(fid,fom);
fclose(fid);

% flash floods use surface soil velocity only, along with streamflow
% velocity
flowch = {'q9m','q9m2','vlc','overlandflow','streamflow','interflow','landint1','landint2','landint3','soilmoisture1','fluxdown','fluxup','q1rec','q9m','vlc','flowdepth'};
L2=0;
if L2==1
    moee = [2 3];
else
    moee = [1 3];
end


% % during W00, window 2 is adjusted to account for when rainfall starts
if p>0
    fid = fopen('first_sttime.txt','r');
    fom = '%d';
    w2bd1 = fscanf(fid,fom);
    fclose(fid);
end


for modee = moee


    
    disp(filnm)
    disp(['mod = ' num2str(mod(filnm,10))])
    disp(['p = ' num2str(p)])
    disp(['filnm using lower to get higher ' num2str(filnm)])
    if mod(filnm,10)==2%window2
        ffnm=['/shared/dondo/home/ml423/world_mts_runs/IRC1/Basin' num2str(gid) 'outputs/refxtmp' num2str(filnm) 'fullresults/iteration' num2str(w2bd1) '-0/' pmm '_SW/' nmm '/' flowch{modee} 'area'];
        % below ffnm is to reduce output size, no iteration outputs anymore
        % 03 06 2025
        ffnm=['/shared/dondo/home/ml423/world_mts_outputs/Basin' num2str(gid) '_output_900m1hr/' pmm '_SW/' nmm '/' flowch{modee} 'area'];
        
        xtmp = load(['/shared/dondo/home/ml423/world_mts_runs/IRC1/Basin' num2str(gid) 'outputs/refxtmp' num2str(filnm) '/xtmpbb1_' num2str(w2bd1) '_0.mat']);
        xtmp = xtmp.xtmpbb1;
        rainfallb1 = xtmp;
         
    elseif mod(filnm,10)==3%window3
        ffnm=['/shared/dondo/home/ml423/world_mts_runs/IRC1/Basin' num2str(gid) 'outputs/refxtmp' num2str(filnm) 'fullresults/iteration' num2str(w3bd1) '-0/' pmm '_SW/' nmm '/' flowch{modee} 'area'];
        % below ffnm is to reduce output size, no iteration outputs anymore
        % 03 06 2025
        ffnm=['/shared/dondo/home/ml423/world_mts_outputs/Basin' num2str(gid) '_output_900m1hr/' pmm '_SW/' nmm '/' flowch{modee} 'area'];
        
        xtmp = load(['/shared/dondo/home/ml423/world_mts_runs/IRC1/Basin' num2str(gid) 'outputs/refxtmp' num2str(filnm) '/xtmpbb1_' num2str(w3bd1) '_0.mat']);
        xtmp = xtmp.xtmpbb1;
        rainfallb1 = xtmp;
        

    elseif mod(filnm,10)==4%window4
        ffnm=['/shared/dondo/home/ml423/world_mts_runs/IRC1/Basin' num2str(gid) 'outputs/refxtmp' num2str(filnm) 'fullresults/iteration' num2str(w4bd1) '-0/' pmm '_SW/' nmm '/' flowch{modee} 'area'];
        % below ffnm is to reduce output size, no iteration outputs anymore
        % 03 06 2025
        ffnm=['/shared/dondo/home/ml423/world_mts_outputs/Basin' num2str(gid) '_output_900m1hr/' pmm '_SW/' nmm '/' flowch{modee} 'area'];
        
        xtmp = load(['/shared/dondo/home/ml423/world_mts_runs/IRC1/Basin' num2str(gid) 'outputs/refxtmp' num2str(filnm) '/xtmpbb1_' num2str(w4bd1) '_0.mat']);
        xtmp = xtmp.xtmpbb1;
        rainfallb1 = xtmp;
         
    elseif mod(filnm,10)==0
        
        ffnm=['/shared/dondo/home/ml423/world_mts_outputs/Basin' num2str(gid) '_output_900m1hr/' pmm '_SW/' nmm '/' flowch{modee} 'area'];
        if rst==2
            tmp = load(['/shared/dondo/home/ml423/world_mts/Precip/Event/ERA5_Basin' num2str(gid) '_' nmm '.mat']);
            era5 = tmp.rain;
            xtmp = era5;
        end
       

        rainfallb1 = xtmp; % original rainfall inputs

    end
    
    save('rainfallb1.mat','rainfallb1')
    disp(ffnm)
    fid=fopen(ffnm,'rb','ieee-le');
    tmpdata=fread(fid,inf,'single');
    overl=reshape(tmpdata,Basins_sum(gloc,7),Basins_sum(gloc,6),720);
    fclose(fid)
    

    
    if modee ==2
        q9m2 = overl;
    elseif modee ==3
        vlc = overl;
    elseif modee ==1
        q9m = overl;
    end
    
end
clear qint3
%clear rainfallb1
maxr = zeros(1,720);
meanr = zeros(1,720);
for te = 1:720
    maxr(te) = max(rainfallb1{te}(intt));
    meanr(te) = mean(rainfallb1{te}(intt));
end
% launchp = find(maxr>0.1);% when rainfall >0.1mm/hr launch particles
% coarselaunch = launchp(1:2:end);

% for the paper, find when rainfall truly begins by defining a thresold of
% 0.5

launchp = find(maxr>0.5);% when rainfall >1mm/hr launch particles
ptmaxthre = 0.5*max(maxr);% used as a threshold to separate events

coarselaunch = launchp(1:2:end);


la2 = find(launchp<(w4bd*8));
la1 = find(launchp<(w3bd*8));
accu0 = 0;
for tmpla = launchp(la1(end)):-1:1
    if ismember(tmpla,launchp)% max p >0.5mm/hr, basic requirement satisfied
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
        % consecutive 24 hours of basin max rain less than 0.5mm/hr,
        % this should be counted as a different rain cell, should not
        % be tracked!

        % one exception is if another major rain is 48-24 ahead, should be
        % considered as well
        rsmallest = max([(tmpla-noraininterval),1]);

        reconsider_period = rsmallest:tmpla;
        secondary_rain = find(maxr(reconsider_period)>0.5);

        if max(maxr(reconsider_period))>ptmaxthre
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


% save the actual rainfall start timings
la0 = find(launchp == rainbegintime);

refinedlaunchp = launchp(la0:la2(end));
disp('launch p below')
launchp
disp('refined launch p below')
refinedlaunchp

if p==0% if is added afternoon 11/04/2024
fid = fopen('pre_rise_pt.txt','w');
fprintf(fid,'%d',refinedlaunchp(1));
fclose(fid);
end
% convert to whole hours
r_start_time = floor(refinedlaunchp(1)/8);

if mod(p,tww)==0%window2, only modifying window 2 to dynamically account for rainfall starts
%if p==0%window2, only modifying window 2 to dynamically account for rainfall starts
        
        if r_start_time>w2bd
            bd = w2bd;
            fid = fopen('bd.txt','w');
            fprintf(fid,'%d',bd);
            fclose(fid);

            bd1 = w2bd-1;
            fid = fopen('bd1.txt','w');
            fprintf(fid,'%d',bd1);
            fclose(fid);

            fid = fopen('first_sttime.txt','w');
            fprintf(fid,'%d',bd1);
            fclose(fid);


            disp('option A rain_start time')
            disp(['bd bd1 are ' num2str(bd) ' ' num2str(bd1)])

            fid = fopen(['Basin_' num2str(gid) '_event' num2str(event) 'window_problem.txt'],'w');
            fprintf(fid,'%d',-9999);
            fclose(fid);
            % rainfall is later than rising point
            % above might be intermittent rainfall with multiple 24 hour
            % <0.5mm/hr rainfall, therefore leading to falsely 'target' rain starting timing
            % latter than rising point
        % above is just to renew the tracking, and just do 1 time step
        elseif r_start_time>w2bd1
            bd = w2bd;
            fid = fopen('bd.txt','w');
            fprintf(fid,'%d',bd);
            fclose(fid);

            bd1 = r_start_time;
            fid = fopen('bd1.txt','w');
            fprintf(fid,'%d',bd1);
            fclose(fid);

            fid = fopen('first_sttime.txt','w');
            fprintf(fid,'%d',bd1);
            fclose(fid);


            disp('option B rain_start time')
            disp(['bd bd1 are ' num2str(bd) ' ' num2str(bd1)])
           
        else
            % rainfall starts way earlier than rising point, so, just use
            % rising points
            bd = w2bd;
            fid = fopen('bd.txt','w');
            fprintf(fid,'%d',bd);
            fclose(fid);

            bd1 = w2bd1;
            fid = fopen('bd1.txt','w');
            fprintf(fid,'%d',bd1);
            fclose(fid);

            fid = fopen('first_sttime.txt','w');
            fprintf(fid,'%d',bd1);
            fclose(fid);

            disp('option C rain_start time')
            disp(['bd bd1 are ' num2str(bd) ' ' num2str(bd1)])
           
        end
end
% above: there is no point running water difference checking before target
% rainfall event starts, therefore reset the window 2 (pre-rising) boundary
% above

% donewith critical change here 11/02/2024


% Simplify the codes
% addpath('C:\Users\ml423\Downloads\m_map')
% tracking is below

for iuy = 4% this is an option, for the paper, it is default 4
    close all

    tic
    c = Basins_sum(gloc,7); r = Basins_sum(gloc,6);basinpixels = c*r;

    wron = cell(720,720);
    for zh = 1:720
        for ju = 1:zh
            wron{zh,ju} = sparse(c,r);
        end
    end
    cal_done = cell(720);
    for t0 = 1:720
        cal_done{t0} = zeros(basinpixels,1);
    end

    ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_fixed_900m/countfile.out'];
    fid=fopen(ffnm);
    facc = fscanf(fid,'%d',[Basins_sum(gloc,7),Basins_sum(gloc,6)]);
    fclose(fid);
    facc(ind) = nan;

    ffnm=['/shared/dondo/home/ml423/world_mts/Basin' num2str(gid) '_fixed_900m/fdr.bin'];
    fid=fopen(ffnm);
    tmpdata=fread(fid,'ubit16');
    fclose(fid);
    fdr = reshape(tmpdata,Basins_sum(gloc,7),Basins_sum(gloc,6));
    fdrfull = fdr;

    q9m(vlc>0) = vlc(vlc>0);

    fac1 = find(facc==1);
    itarget = fac1(1);

    pbb2m = cell(length(fac1),1);
    % getting trajectories for each pixel in the basin 

    for fac1loop = 1:length(fac1)
        itarget = fac1(fac1loop);

        tmp1 = mod(itarget,c);
        tmp2 = ceil(itarget/c);
        tmpz = fdrfull(tmp1,tmp2);
        itarget_trajectory = [];

        for iloop = 1:9999

            if tmpz == 1
                ik = tmp1+1;
                jk = tmp2-1;
                dista = 900*sqrt(2);
            elseif  tmpz == 2
                ik = tmp1+1;
                jk = tmp2;
                dista = 900;
            elseif  tmpz == 4
                ik = tmp1+1;
                jk = tmp2+1;
                dista = 900*sqrt(2);
            elseif  tmpz == 8
                ik = tmp1;
                jk = tmp2+1;
                dista = 900;
            elseif  tmpz == 16
                ik = tmp1-1;
                jk = tmp2+1;
                dista = 900*sqrt(2);
            elseif  tmpz == 32
                ik = tmp1-1;
                jk = tmp2;
                dista = 900;
            elseif  tmpz == 64
                ik = tmp1-1;
                jk = tmp2-1;
                dista = 900*sqrt(2);
            elseif  tmpz == 128
                ik = tmp1;
                jk = tmp2-1;
                dista = 900;
            end

            itarget_trajectory = [itarget_trajectory;[(tmp2-1)*c+tmp1,dista,tmp1,tmp2]];

            tmp1 = ik;
            tmp2 = jk;
            tmpz = fdrfull(ik,jk);
            loctmp = (tmp2-1)*c+tmp1;
            if ismember(loctmp,ind)
                break
            end
        end
        %outlet = itarget_trajectory(end,1);
        %

        endingtrack = min([(w4bd+2)*8,720]);

        launchp_first = refinedlaunchp(1);% the  timestep when particles are launched.

        T_all = [];
        for itra = 1:size(itarget_trajectory,1)
            itarget_series = ((launchp_first-1)*basinpixels+itarget_trajectory(itra,1)):basinpixels:(endingtrack*basinpixels);
            T = itarget_trajectory(itra,2)./abs(q9m(itarget_series))/3600;
            % next pixels caculate
            T_all = [T_all;T];
        end
        %
        [totlen,totsteps] = size(T_all);

        arrival = nan(totlen,totsteps);
        % save the time each particle is launched, spent in the basin, and
        % the time it exits the basin. 

        for t0 = 1:totsteps
            for i0 = 1:totlen
                if cal_done{launchp_first+t0-1}(itarget_trajectory(i0,1)) == 1
                    continue
                else
                    cal_done{launchp_first+t0-1}(itarget_trajectory(i0,1)) = 1;
                    prevaccu = 0;
                    t = t0;
                    i = i0;
                    multipj = 1;
                    for itr = 1:9999
                        % disp(itr)
                        % disp(multipj)
                        % disp(t)
                        accu = prevaccu+multipj*T_all(i,t);
                        if accu>1
                            adj = (1-prevaccu)/T_all(i,t);
                            multipj = multipj-adj;
                            t = t+1;
                            prevaccu = 0;
                        else
                            i = i+1;
                            multipj = 1;
                            prevaccu=accu;
                        end

                        if i == totlen || i>totlen
                            arrival(i0,t0) = launchp_first+(t-1);
                            wron{launchp_first+t-1,launchp_first+t0-1}(itarget_trajectory(i0,1)) = 1;
                            break
                        end
                        if t == totsteps || t>totsteps
                            arrival(i0,t0) = 9999;
                            break
                        end
                    end
                end
            end
        end
        % pbb2m{fac1loop} = arrival - ones(totlen,1).*(launchp_first:1:(totsteps+launchp_first-1))+0.5;
    end
    % pbb2sum = [];
    % for i = 1:length(pbb2m)
    %     pbb2sum = [pbb2sum;pbb2m{i}(pbb2m{i}>0&pbb2m{i}<720)];
    % end

    save(['/shared/dondo/home/ml423/world_mts_runs/IRC1_sup/wrontest5MMORIw' num2str(ceil(p/tww)) num2str(mod(filnm,10)) 'cc' num2str(filnm) '_' num2str(gid) '_' num2str(utch,'%2.2d') '.mat'],'wron')%,'-v7.3')
       
    % 
    % figure
    % h = histogram(pbb2sum,[0:1:720],'Normalization','probability');
    % pbb2 = h.Values;
    % %title(['PDF of travel time 5/15 event b1'])
    % xlabel('Hours')
    % %  set(gca,'XTick',[])
    % xlim([0 200])
    % %  set(gca,'xticklabel',[])cc
    % ylabel('Probability');
    % set(gca,'XTick',0:48:200,'XTickLabel',0:48:200,'FontSize',12,'LineWidth',2.5,'FontWeight','bold')
    % 
    % set(gca,'FontSize',14)
    % 
    % set(gca,'FontName','Times New Roman','FontSize',12,'LineWidth',2.5,'FontWeight','bold');
    % fig = gcf;
    % fig.PaperUnits = 'inches';
    % fig.PaperPosition = [0, 0, 5, 3.1];
    % 
    % if iuy == 4
    %     save(['/shared/dondo/home/ml423/world_mts_runs/IRC1_sup/pbb2test5MMORIw' num2str(ceil(p/tww)) num2str(mod(filnm,10)) 'cc' num2str(filnm) '_' num2str(gid) '_' num2str(utch,'%2.2d') '.mat'],'pbb2')
    %     saveas(gcf,['/shared/dondo/home/ml423/world_mts_runs/IRC1_sup/travel time pdf test5MMORIw' num2str(ceil(p/tww)) num2str(mod(filnm,10)) 'cc' nmm '-' num2str(filnm)],'png')
    % end


    toc

%     pcou3 = {};
%     pcou3 = pcou2{1};
%     for i = 1:length(pcou2)-1
%         pcou3(i+1:end,:) = pcou3(i+1:end,:)+pcou2{i+1};
%     end
%     if iuy == 4
%         save(['pcou3test5MMORIw' num2str(ceil(p/tww)) num2str(mod(filnm,10)) 'cc' nmm '-' num2str(filnm) '.mat'],'pcou3')
%     end
    % if iuy == 5
    %     save('pcou3riverb1l2smw00trackall.mat','pcou3')
    % end
end

quit
    
