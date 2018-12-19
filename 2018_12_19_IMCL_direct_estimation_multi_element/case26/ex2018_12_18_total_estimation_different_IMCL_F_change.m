%%%%%%%%%%%%%%%%%%%%
% �ΏہFcase26
% �g�ʒx���v���t�@�C���ɂ�艹������
% ����x���v���t�@�C���Ǝ����x���v���t�@�C���̑��ݑ���
% ��p���ĉ����𐄒肷��
% �}��   IMCL0%~20%(21���)
%%%%%%%%%%%%%%%%%%%%

%% �����ݒ�i���ʁj
clear
dst_path = sprintf('H:/result/2018_12_06_IMCL_direct_estimation');
load("H:/data/kwave/config/t_pos_2board.mat");
load("H:\data\kwave\medium\2018_09_28_realisticScatter_variousIMCL\corrected\case26_IMCL0.0_pure.mat")
load("H:\data\kwave\result\2018_11_11_case26_variousIMCL\case26_IMCL1.0_pure\rfdata.mat")
load("H:\data\kwave\result\2018_11_11_case26_variousIMCL\case26_IMCL1.0_pure\kgrid.mat")

% �����l
v_fat        = 1450;%[m/s]
v_muscle = 1580;%[m/s]

% IMCL����
IMCL_rate                  = linspace(1,20,20);%[%]
num_IMCL                  = length(IMCL_rate);
v_muscle_with_IMCL = v_fat * IMCL_rate/100 + v_muscle*(1-IMCL_rate/100);%��������[m/s]

% �T���ʒu
assumed_depth          = 19e-3:-kgrid.dx:0;
assumed_distance = 20e-3 - assumed_depth;
num_assumed_depth = length(assumed_depth);
ind_assumed_depth   = zeros(num_assumed_depth,1);% kgrid��ł͋��E�ʒu�͂ǂ̃C���f�b�N�X�ŕ\����邩�����Ƃ߂�D
for i = 1:num_assumed_depth
    ind_assumed_depth(i) = find(single(kgrid.x_vec) == single(assumed_depth(i)));
end
assumed_point         = zeros(2,num_assumed_depth);
assumed_point(1,:) = t_pos(1,51);%���M�t�H�[�J�X�_
assumed_point(2,:) = kgrid.x_vec(ind_assumed_depth);

% �f�q�z�u
t_facing_distance      = 0.04;%[m]
[~,num_receiver,~]  = size(rfdata);
num_transmitter        = num_receiver/2;
num_receiver             = num_receiver/2; 
element_pitch           = abs(t_pos(1,1) - t_pos(1,2));
minimum_elementNum = 20;

% �x���v���t�@�C��
distance_from_assumed_point = zeros(1,num_transmitter);
distance_round_trip                  = zeros(1,num_transmitter);
delay_time_assumed                = zeros(1,num_transmitter);

% ����l
assumed_SOS           = v_muscle:-1:v_muscle_with_IMCL(end);
num_assumed_SOS  = length(v_muscle_with_IMCL(end):v_muscle);
correlation                = zeros(num_assumed_depth,num_assumed_SOS);
estimated_velocity   = zeros(1,num_IMCL);
estimated_IMCL       = zeros(1,num_IMCL);

% ����l
correct_velocity = zeros(1,num_IMCL);
for i = 1:num_IMCL
    correct_velocity(:,i) = v_muscle_with_IMCL(i);
end

%% �������菈����
% ����x���v���t�@�C���Ǝ����x���v���t�@�C���̑��ݑ��ւ����߂�
    for nn = 1:num_IMCL
        
        loadpath = sprintf('H:/data/kwave/result/2018_11_11_case26_variousIMCL/case26_IMCL%0.1f_pure',IMCL_rate(nn));
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
        for  kk = 1:num_assumed_depth
            % �쓮�f�q�̑I��
                target_element = find_target_element(assumed_distance,num_receiver,t_pos,minimum_elementNum,kk);
            
            for ll = 1:num_assumed_SOS
                
                % ���M�t�H�[�J�X
                [focused_rfdata,distance_from_assumed_point] = transmit_focusing(num_transmitter,assumed_point,t_pos,kk,ll,assumed_SOS,kgrid,rfdata_echo_only,target_element,num_sample,num_receiver);
                
                % ���ȑ��֌W���̍ő�l���e�M���ŋ��߂�irfdata, source wave�j
                auto_correlation_rfdata            = diag(focused_rfdata.' * focused_rfdata);
                auto_correlation_source_wave = interp_sourcewave.' * interp_sourcewave;
                auto_correlation_source_wave = repmat(auto_correlation_source_wave,num_receiver,1);
                
                % �x���v���t�@�C���̉���
                distance_round_trip   = distance_from_assumed_point + min(distance_from_assumed_point);%[m]
                delay_time_assumed = round(distance_round_trip / assumed_SOS(ll) / (kgrid.dt/4));%[sample]
                delay_sourcewave     = zeros(num_sample*4, num_receiver);
                
                for ii = 1:length(target_element)
                    source_wave2cat = interp_sourcewave;
                    source_wave2cat(num_sample*4-delay_time_assumed(target_element(ii))+1:end,1)=NaN;
                    source_wave2cat(isnan(source_wave2cat)) = [];
                    delay_sourcewave(:,target_element(ii)) = cat(1,zeros(delay_time_assumed(target_element(ii)),1),source_wave2cat);
                    % ���ݑ��ւ̐ώZ
                    correlation(kk,ll) = correlation(kk,ll) + (focused_rfdata(:,target_element(ii)).'*delay_sourcewave(:,target_element(ii))...
                        /sqrt(auto_correlation_rfdata(target_element(ii),1)*auto_correlation_source_wave(target_element(ii),1)));
                end
                
                correlation(kk,ll) = correlation(kk,ll)/length(target_element);
                
            end
            
            dispname = sprintf('estimated depth # is %d, IMCL # is %d, target element # is %d',kk,nn,length(target_element));
            disp(dispname) %#ok<DSPS>
            
        end
        
        [~,ind_estimate_v] = max(max(correlation));
        [~,ind_estimate_d] = max(correlation(:,ind_estimate_v));
        estimated_velocity(1,nn) = assumed_SOS(1,ind_estimate_v);
        estimated_IMCL(1,nn) = 100*((v_muscle-estimated_velocity(1,nn))/(v_muscle-v_fat));
        
        % delay sourcewave �� rfdata ���d�˂ĕ\�����邽�߂̏����D%%%%%%
        target_element = find_target_element(assumed_distance,num_receiver,t_pos,minimum_elementNum,ind_estimate_d);
        [focused_rfdata,distance_from_assumed_point] = transmit_focusing(num_transmitter,assumed_point,t_pos,ind_estimate_d,ind_estimate_v,assumed_SOS,kgrid,rfdata_echo_only,target_element,num_sample,num_receiver);
        distance_round_trip = distance_from_assumed_point + min(distance_from_assumed_point);%[m]
        delay_time_assumed = round(distance_round_trip / assumed_SOS(ind_estimate_v) / (kgrid.dt/4));%[sample]
        
                
        % �摜�ۑ�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dst_path = sprintf('H:/result/2018_12_19_IMCL_direct_estimation_multi_element/case26/2018_12_19_focused/IMCL%d%%',IMCL_rate(nn));
        if ~exist(dst_path, 'dir')
            mkdir(dst_path);
        end
        
        figure;
        imagesc(IMCL_rate,20-assumed_depth*1e3,correlation);
        xlabel('IMCL content[%]')
        ylabel('depth[mm]')
        titlename = sprintf('IMCL: %d %%',IMCL_rate(nn)*1e3);
        title({'Correlation of delay curve';titlename})
        colorbar
        savefilename = sprintf('/correlation');
        savefig([dst_path,savefilename,'.fig'])
        exportfig([dst_path,savefilename],'png',[300,200])
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
        savefig([dst_path,savefilename,'.fig'])
        exportfig([dst_path,savefilename],'png',[350,300])
        close gcf
        
        % �ϐ��ۑ�
        savefilename = '/result';
        save([dst_path,savefilename],'correlation','focused_rfdata','ind_estimate_d','ind_estimate_v');
    end

%% �ۑ���
dst_path2 = sprintf('H:/result/2018_12_19_IMCL_direct_estimation_multi_element/case26/2018_12_19_focused/figure');
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
    

dst_path3 = sprintf('H:/result/2018_12_19_IMCL_direct_estimation_multi_element/case26/2018_12_19_focused/');
if ~exist(dst_path3, 'dir')
    mkdir(dst_path3);
end
savefilename = sprintf('2018_12_19_all_result');
save([dst_path3,savefilename],'estimated_IMCL','estimated_velocity',...
    'v_muscle_with_IMCL','target_element','correct_velocity','IMCL_rate','t_pos');

%% �֐���

function target_element = find_target_element(assumed_distance,num_receiver,t_pos,minimum_elementNum,kk)
       target_element = find((-assumed_distance(1,kk)<=t_pos(1,1:num_receiver))&((t_pos(1,1:num_receiver)<=assumed_distance(1,kk))));
        if length(target_element) < minimum_elementNum
            target_element = ceil((num_receiver-minimum_elementNum)/2+1):ceil((num_receiver+minimum_elementNum)/2) ;
       end
end

function [focused_rfdata,distance_from_assumed_point] = transmit_focusing(num_transmitter,assumed_point,t_pos,kk,ll,assumed_SOS,kgrid,rfdata_echo_only,target_element,num_sample,num_receiver)
       distance_from_assumed_point = zeros(1,num_transmitter);

    for ii = 1:num_transmitter
        distance_from_assumed_point(1,ii) = norm(assumed_point(:,kk)-t_pos(:,ii));%[m]
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