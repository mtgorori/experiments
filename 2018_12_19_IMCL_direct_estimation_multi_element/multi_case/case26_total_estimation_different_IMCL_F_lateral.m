%%%%%%%%%%%%%%%%%%%%
% �ΏہFcase1
% �g�ʒx���v���t�@�C���ɂ�艹������
% ����x���v���t�@�C���Ǝ����x���v���t�@�C���̑��ݑ���
% ��p���ĉ����𐄒肷��
% �}��   IMCL0%
%%%%%%%%%%%%%%%%%%%%
%% �����ݒ�i���ʁj
clear
load("H:/data/kwave/config/t_pos_2board.mat");
load("H:\data\kwave\medium\2019_01_24_realisticScatter_for_5MHz\case26_IMCL0.0.mat")
load("H:\data\kwave\result\2019_01_24_realisticScatter_for_5MHz\case26_IMCL0.0\rfdata.mat")
load("H:\data\kwave\result\2019_01_24_realisticScatter_for_5MHz\case26_IMCL0.0\kgrid.mat")
load("H:\experiments\2018_12_19_IMCL_direct_estimation_multi_element\multi_case\condition\2019_01_25_case26.mat")
% �����l
v_fat        = 1450;%[m/s]
v_muscle = 1580;%[m/s]

% IMCL�����i�����j
IMCL_rate                  = linspace(0,20,11);%[%]
num_IMCL                  = length(IMCL_rate);
v_muscle_with_IMCL = v_fat * IMCL_rate/100 + v_muscle*(1-IMCL_rate/100);%��������[m/s]

% �f�q�z�u
t_facing_distance      = 0.04;%[m]
[~,num_receiver,~]  = size(rfdata);
num_transmitter        = num_receiver/2;
num_receiver             = num_receiver/2;
element_pitch           = abs(t_pos(1,1) - t_pos(1,2));
minimum_elementNum = 20;
lateral_range_max = 4.0*1e-3;
lateral_range_min = -4.0*1e-3;
num_lateral = round((lateral_range_max -lateral_range_min) / element_pitch)+1;%0.0~4.8 mm�܂�
lateral_focus_point = linspace(lateral_range_min,lateral_range_max,num_lateral);

% �T���ʒu
assumed_depth          = 19e-3:-kgrid.dx*2:0;
[~,ind_end_point_depth] = min(abs(assumed_depth - 10e-3));
assumed_distance      = 20e-3 - assumed_depth;
num_assumed_depth = length(assumed_depth);
ind_assumed_depth   = zeros(num_assumed_depth,1);% kgrid��ł͋��E�ʒu�͂ǂ̃C���f�b�N�X�ŕ\����邩�����Ƃ߂�D
assumed_point         = zeros(2,num_assumed_depth,num_lateral);
for i = 1:num_assumed_depth
    ind_assumed_depth(i) = find(single(kgrid.x_vec) == single(assumed_depth(i)));
end
for i = 1:num_lateral
    assumed_point(1,:,i) = lateral_focus_point(i);% ���M�t�H�[�J�X�_
    assumed_point(2,:,i) = kgrid.x_vec(ind_assumed_depth);
end

% �x���v���t�@�C��
distance_from_assumed_point = zeros(1,num_transmitter);
distance_round_trip                  = zeros(1,num_transmitter);
delay_time_assumed                = zeros(1,num_transmitter);
window_left = 128;%���֑�
window_right = 301;
% ����l
assumed_IMCL_rate                  = linspace(1,20,20);%[%]
num_assumed_IMCL                  = length(assumed_IMCL_rate);
v_muscle_with_assumed_IMCL = v_fat * assumed_IMCL_rate/100 + v_muscle*(1-assumed_IMCL_rate/100);%��������[m/s]
assumed_SOS = linspace(v_muscle,v_muscle_with_assumed_IMCL(end),num_assumed_IMCL);
num_assumed_SOS  = length(assumed_SOS);
correlation                = zeros(num_assumed_SOS,ind_end_point_depth,num_lateral);
estimated_velocity   = zeros(1,num_IMCL);
estimated_IMCL       = zeros(1,num_IMCL);
best_lateral              = zeros(1,num_IMCL);

% ����l
correct_velocity = zeros(1,num_IMCL);
for i = 1:num_IMCL
    correct_velocity(:,i) = v_muscle_with_IMCL(i);
end


%% �������菈����
% ����x���v���t�@�C���Ǝ����x���v���t�@�C���̑��ݑ��ւ����߂�
for nn = 1:num_IMCL
    
    loadpath = sprintf('H:/data/kwave/result/2019_01_24_realisticScatter_for_5MHz/case26_IMCL%0.1f',IMCL_rate(nn));
    load([loadpath,'/rfdata.mat'])
    load([loadpath,'/kgrid.mat'])
    load([loadpath,'/sourse_wave.mat'])
    [num_sample,~,~] = size(rfdata);
    
    
    % rf�f�[�^���`%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % �T���v������4�{�ɂ��邽�߂ɃX�v���C�����
    interp_sourcewave                     = interp1(source_wave,linspace(1,num_sample,num_sample*4)','spline');
    [~,offset_interp_sourcewave]   = max(interp_sourcewave);% ���M�g�`�̍ő�l�����_�D�x���Ȑ��ɎU�z�}���d�ˍ��킹�邱�ƂɎg���D
    rfdata_echo_only                       = zeros(num_sample,num_receiver,num_transmitter);
    
    for ii = 1:num_transmitter
        for jj = 1:num_receiver
            delay_transmitted_wave           = round(((abs(t_pos(1,jj) - t_pos(1,ii)))/v_muscle)/(kgrid.dt));
            rfdata_echo_only(:,jj,ii)            = rfdata(:,jj,ii);
            rfdata_echo_only(1:delay_transmitted_wave+50,jj,ii) = 0;
        end
    end
    
    % ��������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for mm = 1:num_lateral
        for  kk = 1:ind_end_point_depth
            % �쓮�f�q�̑I��
            target_element = find_target_element(assumed_distance,assumed_point,lateral_focus_point,num_receiver,t_pos,minimum_elementNum,kk,mm);
            
            for ll = 1:num_assumed_SOS
                % ���M�t�H�[�J�X
                [focused_rfdata,distance_from_assumed_point] = transmit_focusing(num_transmitter,assumed_point,t_pos,kk,ll,mm,assumed_SOS,kgrid,rfdata_echo_only,target_element,num_sample,num_receiver);

                % �x���v���t�@�C���̉���
                distance_round_trip   = distance_from_assumed_point + min(distance_from_assumed_point);%[m]
                delay_time_assumed = round(distance_round_trip / assumed_SOS(ll) / (kgrid.dt/4));%[sample]
                
                for ii = 1:length(target_element)
                    tmp_delay = delay_time_assumed(target_element(ii));
                    source_wave2cat = interp_sourcewave;
                    source_wave2cat(num_sample*4-tmp_delay+1:end,1)=NaN;
                    source_wave2cat(isnan(source_wave2cat)) = [];
                    delay_sourcewave = cat(1,zeros(tmp_delay,1),source_wave2cat);
                    % ���ȑ��֌W���̍ő�l���e�M���ŋ��߂�irfdata, source wave�j
                    hilb_focused_rfdata = hilbert(focused_rfdata(:,target_element(ii)));
                    hilb_delay_sourcewave = hilbert(delay_sourcewave);
                    delay_sourcewave_trim = hilb_delay_sourcewave(tmp_delay+window_left:tmp_delay+window_right);
                    focused_rfdata_trim = hilb_focused_rfdata(tmp_delay+window_left:tmp_delay+window_right);
                    auto_correlation_rfdata            =  focused_rfdata_trim' * focused_rfdata_trim;
                    auto_correlation_source_wave = delay_sourcewave_trim' * delay_sourcewave_trim;
                    tmp_xcorr = (focused_rfdata_trim' * delay_sourcewave_trim) / sqrt(auto_correlation_rfdata * auto_correlation_source_wave);
                    tmp_xcorr = abs(tmp_xcorr) * sign(real(tmp_xcorr));
                    % ���ݑ��ւ̐ώZ
                    correlation(ll,kk,mm) = correlation(ll,kk,mm) + tmp_xcorr;%...
                end
                
                correlation(ll,kk,mm) = correlation(ll,kk,mm)/length(target_element);
                
            end
            
            dispname = sprintf('estimated depth # is %d, lateral # is %0.1f, IMCL # is %d, target element # is %d',kk,mm,nn,length(target_element));
            disp(dispname) %#ok<DSPS>
            
        end
    end
    
    [~,ind_estimate_l] = max(max(max(correlation)));
    [~,ind_estimate_d] = max(max(correlation(:,:,ind_estimate_l)));
    [~,ind_estimate_v]  = max(correlation(:,ind_estimate_d,ind_estimate_l));
    estimated_velocity(1,nn) = assumed_SOS(1,ind_estimate_v);
    estimated_IMCL(1,nn) = 100*((v_muscle-estimated_velocity(1,nn))/(v_muscle-v_fat));
    best_lateral(1,nn) = lateral_focus_point(1,ind_estimate_l);
    % delay sourcewave �� rfdata ���d�˂ĕ\�����邽�߂̏����D%%%%%%
    target_element = find_target_element(assumed_distance,assumed_point,lateral_focus_point,num_receiver,t_pos,minimum_elementNum,ind_estimate_d,ind_estimate_l);
    [focused_rfdata,distance_from_assumed_point] = transmit_focusing(num_transmitter,assumed_point,t_pos,ind_estimate_d,ind_estimate_v,ind_estimate_l,assumed_SOS,kgrid,rfdata_echo_only,target_element,num_sample,num_receiver);
    distance_round_trip = distance_from_assumed_point + min(distance_from_assumed_point);%[m]
    delay_time_assumed = round(distance_round_trip / assumed_SOS(ind_estimate_v) / (kgrid.dt/4));%[sample]
    tmp_delay = delay_time_assumed(round(median(target_element)));
    source_wave2cat = interp_sourcewave;
    source_wave2cat(num_sample*4-tmp_delay+1:end,1)=NaN;
    source_wave2cat(isnan(source_wave2cat)) = [];
    delay_sourcewave = cat(1,zeros(tmp_delay,1),source_wave2cat);
    t_array = interp1(kgrid.t_array,linspace(1,num_sample,num_sample*4)','spline');
    % �摜�ۑ�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cat_dst_path = sprintf('IMCL%d%%',IMCL_rate(nn));
    dst_path1 = [dst_path,cat_dst_path];
    
    if ~exist(dst_path1, 'dir')
        mkdir(dst_path1);
    end
    
    figure;
    imagesc(assumed_IMCL_rate,20-assumed_depth(1:ind_end_point_depth)*1e3,correlation(:,:,ind_estimate_l)');
    xlabel('IMCL content[%]')
    ylabel('depth[mm]')
    titlename = sprintf('IMCL: %0.1f %%, lateral: %0.1f mm',IMCL_rate(nn),lateral_focus_point(ind_estimate_l));
    title({'Correlation of delay curve';titlename})
    colorbar
    savefilename = sprintf('/correlation');
    savefig([dst_path1,savefilename,'.fig'])
    exportfig([dst_path1,savefilename],'png',[300,200])
    close gcf
    
    figure;
    imagesc(focused_rfdata);
    hold on
    scatter(target_element,delay_time_assumed(target_element)+offset_interp_sourcewave,'r.')
    caxis([min(min(focused_rfdata))/5 max(max(focused_rfdata))/5])
    ylim([min(delay_time_assumed(target_element)) max(delay_time_assumed(target_element)+2*offset_interp_sourcewave)])
    xlabel('element number')
    ylabel('time[sample]')
    axis square
    title(titlename)
    colorbar
    savefilename = sprintf('/solved_delay_curve');
    savefig([dst_path1,savefilename,'.fig'])
    exportfig([dst_path1,savefilename],'png',[350,300])
    close gcf
    
    figure;
    plot(t_array(tmp_delay+window_left:tmp_delay+window_right),...
        focused_rfdata(tmp_delay+window_left:tmp_delay+window_right,round(median(target_element)))...
        /max(abs(focused_rfdata(tmp_delay+window_left:tmp_delay+window_right,round(median(target_element))))))
    hold on
    plot(t_array(tmp_delay+window_left:tmp_delay+window_right),...
        delay_sourcewave(tmp_delay+window_left:tmp_delay+window_right)...
        /max(abs(delay_sourcewave(tmp_delay+window_left:tmp_delay+window_right))))
    xlabel('����[sample]')
    ylabel('����[a.u.]')
    xlim([t_array(tmp_delay+window_left) t_array(tmp_delay+window_right)])
    title(titlename)
    legend('���M��','����x���M��','Location','bestoutside')
    savefilename = sprintf('/singal_match');
    savefig([dst_path1,savefilename,'.fig'])
    exportfig([dst_path1,savefilename],'png',[500,200])
    close gcf
    
    % �ϐ��ۑ�
    savefilename = '/result';
    save([dst_path1,savefilename],'correlation','focused_rfdata','ind_estimate_d','ind_estimate_v','ind_estimate_l');
end

%% �ۑ���
if ~exist(dst_path2, 'dir')
    mkdir(dst_path2);
end

figure;
plot(correct_velocity(1,:),estimated_velocity(1,:),'LineWidth',1);
hold on
plot(correct_velocity(1,:),correct_velocity(1,:),'k--','LineWidth',0.25);
hold off
xlabel('correct velocity [m/s]')
ylabel('estimated velocity [m/s]')
xlim([v_fat*0.2+v_muscle*0.8 v_muscle])
ylim([v_fat*0.2+v_muscle*0.8 v_muscle])
savefilename = sprintf('/velocity');
savefig([dst_path2,savefilename,'.fig'])
exportfig([dst_path2,savefilename],'png',[300,300])
close gcf

figure;
plot(IMCL_rate(1,:),estimated_IMCL(1,:),'LineWidth',1);
hold on
plot(IMCL_rate(1,:),IMCL_rate(1,:),'k--','LineWidth',0.25);
hold off
xlabel('correct IMCL content [%]')
ylabel('estimated IMCL content [%]')
xlim([0 20])
ylim([0 20])
savefilename = sprintf('/IMCL');
savefig([dst_path2,savefilename,'.fig'])
exportfig([dst_path2,savefilename],'png',[300,300])
close gcf


if ~exist(dst_path3, 'dir')
    mkdir(dst_path3);
end
clear focused_rfdata
clear rfdata
clear rfdata_echo_only
savefilename = sprintf('all_result');
save([dst_path3,savefilename]);

%% �֐���
function target_element = find_target_element(assumed_distance,assumed_point,lateral_focus_point,num_receiver,t_pos,minimum_elementNum,kk,mm)
target_element = find((-assumed_distance(1,kk) + assumed_point(1,1,mm)<=t_pos(1,1:num_receiver))&...
    ((t_pos(1,1:num_receiver)<=assumed_distance(1,kk)+ assumed_point(1,1,mm))));
if length(target_element) < minimum_elementNum
    [~,ind_central_element] = min(abs(t_pos(1,:)-lateral_focus_point(mm)));
    target_element = ceil(ind_central_element-minimum_elementNum/2+1):ceil(ind_central_element+minimum_elementNum/2);
end
end


function [focused_rfdata,distance_from_assumed_point] = transmit_focusing(num_transmitter,assumed_point,t_pos,kk,ll,mm,assumed_SOS,kgrid,rfdata_echo_only,target_element,num_sample,num_receiver)
distance_from_assumed_point = zeros(1,num_transmitter);

for ii = 1:num_transmitter
    distance_from_assumed_point(1,ii) = norm(assumed_point(:,kk,mm)-t_pos(:,ii));%[m]
end

delay_length_focusing = distance_from_assumed_point - min(distance_from_assumed_point);
delay_time_focusing    = round(delay_length_focusing/ assumed_SOS(ll) / (kgrid.dt));
focused_rfdata        = transmit_focus(rfdata_echo_only,target_element,...
    num_sample,delay_time_focusing,num_receiver);
interp_rfdata = zeros(num_sample*4,num_receiver);

for jj = 1:num_receiver
    interp_rfdata(:,jj) = interp1(focused_rfdata(:,jj),linspace(1,num_sample,num_sample*4),'spline');
end

focused_rfdata = interp_rfdata;
end