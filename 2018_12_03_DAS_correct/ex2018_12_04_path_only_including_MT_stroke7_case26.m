%%%%%%%%%%%%%%%%%%%%%%%%%%%%��������������%%%%%
% �ΏہFcase26
% �œ_�����ʒu�Œ�Fy=0
% ���E�ʂ����o���Q�Ɠ_�̑Ó����]��
% ch 50�����̃s�N�Z���Ƀt�H�[�J�X��������D
% �Q�Ɠ_�F�������Z�ɂ�����Q�Ɠ_
% ��r����_�F�ech�ł�rf�M���ő�l(��Βl�͂Ƃ�Ȃ�)
%�댟�o�������邪�C���ˋ��x���l�����Č댟�o�Ɣ��ʂł���悤�ɂ���D
% �t�H�[�J�X�[�x���P�O���b�h���ω������āC�J��������M�f�[�^�𐶐�����D
% F�l���Œ肷��D�ŋߐڋ����������ɐݒ肷��H
% x-axis: 20 mm�̂Ƃ��ɑS�f�q���g���Ƃ����O���݂���D
% �Ώۂɂ���}���f�[�^�Fcase26��p����D
% IMCL������0 %�ɌŒ肷��D
% kgrid.x_vec ��0�ƂȂ�̂�251�Ԗڂ̗v�f�D
% �O���b�h����0.1 mm
% �f�q�ԃs�b�`��0.4 mm
% ����āC�ŋߐڐ[�x(�ŋߐڋ�����0.4 mm)
% �[�x���Ƃׂ̍����f�q�����́Cfloor()��p����D
% �������Z�̍ۂɎQ�Ƃ��鉹����1580 m/s�Ƃ���D
% �������Z�̑O�ɎQ�Ɠ_�Ǝ�Mch���Ƃ̐U��(�q���x���g�ϊ����Βl)�ő�l�Ƃ̋�����
% ���v��]���֐��Ƃ��Ĕ}���̋ώ�����]�����邱�Ƃ������ɍs���D
% update:rf-data�̔z��T�C�Y���傫���Ȃ��Ă����̂ŁC�ϐ��敪���ו������Ċe�ϐ��̌Ăяo�����x���グ��D[2018/11/05]
% update:for�����ƂɃt�H���_���쐬����悤�ɂ���D�t�H���_���Ƃ�rf�f�[�^��ۑ�����D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% �����ݒ�ɕK�v�ȃf�[�^�̌Ăяo��
clear;
load("H:/data/kwave/config/t_pos_2board.mat");
pathname = sprintf('H:/data/kwave/result/2018_11_11_case26_variousIMCL/case26_IMCL%0.1f_pure',1.0);
cd(pathname);
load('rfdata.mat');
load('kgrid.mat');

%% �����ݒ�
v_fat = 1450;%[m/s]
v_muscle = 1580;%[m/s]
rate_IMCL = linspace(1,20,20);
v_reference = zeros(1,length(rate_IMCL));
t_facing_distance = 0.04;%[m]
[~,num_receiver,num_transmitter] = size(rfdata);
num_echo_receiver = num_transmitter;
num_rate_IMCL = length(rate_IMCL);
num_medium  = num_rate_IMCL;
reference_point = zeros(1,num_echo_receiver);
reference_point_lowerlimit = zeros(1,num_echo_receiver);%�ώ����]���̂��߂�RF�f�[�^�}�X�L���O�Ɏg��
reference_point_upperlimit = zeros(1,num_echo_receiver);%�ώ����]���̂��߂�RF�f�[�^�}�X�L���O�Ɏg��
point_max_in_mask = zeros(1,num_echo_receiver);%�}�X�N�������RF�f�[�^�ŐU���ő�̃T���v���_���
distance_from_focal_point_all = zeros(1,num_echo_receiver);

num_depth = (t_pos(2,1)-t_pos(2,101))/kgrid.dx/2 - 3;%'3'�Ƃ���̂́C�ŋߐڋ�����0.4 mm�ł��邱�Ƃ��l�����Ă���D
focal_depth = zeros(1,num_depth);
focal_point = zeros(2,num_depth);
element_pitch = abs(t_pos(1,1) - t_pos(1,2));
lateral_range_max = 4.8*1e-3;
lateral_range_min = 0;
num_lateral = round((lateral_range_max -lateral_range_min) / element_pitch)+1;%0.0~4.8 mm�܂�

ind_max_signal = zeros(num_rate_IMCL,num_medium,num_lateral);
max_signal = zeros(num_rate_IMCL,num_medium,num_lateral);
ind_search_depth = zeros(num_lateral,num_medium);%�T���͈͂����ڂ̐���Ō��߂�D���̍ۂ̒T���͈͂̒����l
displacement = zeros(num_echo_receiver,num_rate_IMCL,num_medium,num_lateral);
displacement_ref = zeros(num_echo_receiver,num_rate_IMCL,num_medium,num_lateral);
error_wavefront = zeros(num_rate_IMCL,num_medium,num_lateral);
focal_signal_total = zeros(num_depth,num_rate_IMCL,num_medium,num_lateral);
focal_signal = zeros(num_rate_IMCL,num_medium,num_lateral);
detected_boundary = zeros(num_rate_IMCL,num_medium,num_lateral);

%% ������
for mm = 1:num_lateral
    for ii = 9:num_depth
        focal_depth(1,ii) = single((ii+3)*kgrid.dx);
        focal_point(2,ii) = single(t_pos(2,1)-focal_depth(1,ii));
        focal_point(1,ii) = single((mm-1) * element_pitch);
    end
    for ll = 1:num_medium
        %% �Q�Ɣz��̌Ăяo��
        pathname = sprintf('H:/data/kwave/result/2018_11_11_case26_variousIMCL/case26_IMCL%0.1f_pure',ll);
        cd(pathname);
        load('rfdata.mat');
        load('kgrid.mat');
        [num_sample,~,~] = size(rfdata);
        %% ���E�ʒu�̌���I�肷��D(��␔�Fnum_IMCL x num_medium)
        for kk = 1:num_rate_IMCL
            %% ���E�ʒu�̃A�^��������(full depth)
            if kk == 1
                %% ���ˋ��x�v���t�@�C�������߂�D
                v_reference(1,kk) = v_muscle*(1-rate_IMCL(1,kk)/100) + v_fat*(rate_IMCL(1,kk)/100);
                for ii = 9:num_depth
                    target_element = find((-focal_depth(1,ii)/2+focal_point(1,ii)<=t_pos(1,1:100)&(t_pos(1,1:100)<=focal_depth(1,ii)/2+focal_point(1,ii))));
                    %��M�p�̎Q�Ɠ_�Z�o
                    for jj = 1:num_echo_receiver
                        distance_from_focal_point_all(1,jj) = norm(t_pos(:,jj) - focal_point(:,ii));
                        delay_time_all = round(((distance_from_focal_point_all - focal_depth(1,ii))/v_reference(1,kk))/kgrid.dt);%[sample]
                        reference_point(1,jj) = round(delay_time_all(1,jj)+(2*focal_depth(1,ii)/v_reference(1,kk))/kgrid.dt+25-1);
                    end
                    %���M�r�[���t�H�[�~���O�i���ʁj
                    focused_rfdata = transmit_focus(rfdata,target_element,num_sample,delay_time_all,num_echo_receiver);
                    %��M�r�[���t�H�[�J�V���O
                    focal_signal_total(ii,kk,ll,mm) =  receive_focus(focused_rfdata,target_element,reference_point);
                end
                %% ���ˋ��x���ő�̏œ_�ʒu��T���D
                [max_tmp,ind_tmp] = findpeaks(abs(focal_signal_total(9:end,kk,ll,mm)),'Npeaks',1,'SortStr','descend');
                [min_ind_tmp,ind_ind_tmp] = min(ind_tmp);
                ind_max_signal(kk,ll,mm) = min_ind_tmp+8;
                max_signal(kk,ll,mm) = focal_signal_total(ind_max_signal(kk,ll,mm),kk,ll,mm);
                ind_search_depth(mm,ll) = ind_max_signal(kk,ll,mm);
                ii = ind_search_depth(mm,ll);
                %% �Q�Ɠ_�Ɣg�ʂƂ̌`�󑊊ւ����߂�D�i���ݑ��֖@�j
                target_element = find((-focal_depth(1,ii)/2+focal_point(1,ii)<=t_pos(1,1:100)&(t_pos(1,1:100)<=focal_depth(1,ii)/2+focal_point(1,ii))));
                %��M�p�̎Q�Ɠ_�Z�o
                for jj = 1:num_echo_receiver
                    distance_from_focal_point_all(1,jj) = norm(t_pos(:,jj) - focal_point(:,ii));
                    delay_time_all = round(((distance_from_focal_point_all - focal_depth(1,ii))/v_reference(1,kk))/kgrid.dt);%[sample]
                    reference_point(1,jj) = round(delay_time_all(1,jj)+(2*focal_depth(1,ii)/v_reference(1,kk))/kgrid.dt+25-1);
                    %25��focal_amplitude���ő�ɂ���I�t�Z�b�g�D
                end
                max_delay = abs(diff(reference_point)) + 2;
                %���M�r�[���t�H�[�~���O�i���ʁj
                focused_rfdata = transmit_focus(rfdata,target_element,num_sample,delay_time_all,num_echo_receiver);
                %RF�f�[�^�}�X�L���O�i�ώ����]���̂���)
                min_reference = min(reference_point) - 5;
                max_reference = min(reference_point) + 40 - 1;
                mask_rfdata = zeros(size(focused_rfdata));
                mask_rfdata(min_reference:max_reference,:) = 1;
                focused_rfdata_masked = mask_rfdata .* focused_rfdata;
                for jj = 1:length(target_element)-1
                    %�g�ʌ`��]��
                    displacement(target_element(jj),kk,ll,mm)  = - finddelay(focused_rfdata(:,target_element(jj+1)),...
                        focused_rfdata_masked(:,target_element(jj)),max_delay(target_element(jj)));
                end
                displacement_ref(target_element(1:end-1),kk,ll,mm) = diff(reference_point(target_element));
                error_wavefront(kk,ll,mm) = sum(abs(displacement(target_element(1:end-1),kk,ll,mm) - ...
                    displacement_ref(target_element(1:end-1),kk,ll,mm)))/(length(target_element)-1);
                %% �o�̓f�[�^�̕ۑ�
                dst_path = sprintf('H:/result/2018_12_03_DAS_correct/2018_12_04_case26/lateral%0.1fmm/true%d_assumption%d',...
                    focal_point(1,9)*1e3,rate_IMCL(1,ll),rate_IMCL(1,kk));
                if ~exist(dst_path, 'dir')
                    mkdir(dst_path);
                end
                displacement_tmp = displacement(:,kk,ll,mm);
                save([dst_path,'\pre_rfdata.mat'],'focused_rfdata','focused_rfdata_masked','displacement_tmp','reference_point');
                figure;
                imagesc(focused_rfdata_masked);
                hold on
                scatter(min(target_element):max(target_element),reference_point(min(target_element):max(target_element)),'red');
                hold off
                xlabel('receiver[ch]');
                ylabel('time[sample]');
                axis square;
                axis tight;
                xlim([min(target_element) max(target_element)])
                ylim([min_reference max_reference])
                colormap(bone);
                colorbar;
                caxis([min(min(focused_rfdata_masked)) max(max(focused_rfdata_masked))])
                savefig([dst_path,'\pre_rfdata_masked.fig'])
                exportfig([dst_path,'\pre_rfdata_masked'],'png',[400,400])
                close gcf
                figure;
                imagesc(focused_rfdata);
                hold on
                scatter(min(target_element):max(target_element),reference_point(min(target_element):max(target_element)),2,'red');
                scatter(min(target_element):max(target_element),point_max_in_mask(min(target_element):max(target_element)),2,'blue','filled');
                hold off
                xlabel('receiver[ch]');
                ylabel('time[sample]');
                axis square;
                axis tight;
                colormap(bone);
                colorbar;
                caxis([min(min(focused_rfdata_masked)) max(max(focused_rfdata_masked))])
                savefig([dst_path,'\pre_rfdata_whole.fig'])
                exportfig([dst_path,'\pre_rfdata_whole'],'png',[400,400])
                close gcf
                figure;
                plot(focal_depth*1e3,focal_signal_total(:,kk,ll,mm));
                hold on
                scatter(focal_depth(ind_max_signal(kk,ll,mm))*1e3,max_signal(kk,ll,mm),'red','filled')
                hold off
                xlabel('�œ_�[��[mm]')
                ylabel('�M�����x[au]')
                savefig([dst_path,'\pre_focal_signal.fig'])
                exportfig([dst_path,'\pre_focal_signal'],'png',[400,400])
                close gcf
                fprintf('prepareing... lateral number is %d, medium number is %d, and estimation number is %d\n',mm,ll,kk);
            else
            %% �����i���ċ��E�ʒu�����߂�D(limited depth)
                %% ���ˋ��x�v���t�@�C�������߂�D
                v_reference(1,kk) = v_muscle*(1-rate_IMCL(1,kk)/100) + v_fat*(rate_IMCL(1,kk)/100);
                for ii = ind_search_depth(mm,ll)-5:ind_search_depth(mm,ll)+5
                    target_element = find((-focal_depth(1,ii)/2+focal_point(1,ii)<=t_pos(1,1:100)&(t_pos(1,1:100)<=focal_depth(1,ii)/2+focal_point(1,ii))));
                    %��M�p�̎Q�Ɠ_�Z�o
                    for jj = 1:num_echo_receiver
                        distance_from_focal_point_all(1,jj) = norm(t_pos(:,jj) - focal_point(:,ii));
                        delay_time_all = round(((distance_from_focal_point_all - focal_depth(1,ii))/v_reference(1,kk))/kgrid.dt);%[sample]
                        reference_point(1,jj) = round(delay_time_all(1,jj)+(2*focal_depth(1,ii)/v_reference(1,kk))/kgrid.dt+25-1);
                    end
                    %���M�r�[���t�H�[�~���O�i���ʁj
                    focused_rfdata = transmit_focus(rfdata,target_element,num_sample,delay_time_all,num_echo_receiver);
                    %��M�r�[���t�H�[�J�V���O
                    focal_signal_total(ii,kk,ll,mm) =  receive_focus(focused_rfdata,target_element,reference_point);
                end
                %% ���ˋ��x���ő�̏œ_�ʒu��T���D
                [max_tmp,ind_tmp] = findpeaks(abs(focal_signal_total(9:end,kk,ll,mm)),'Npeaks',1,'SortStr','descend');
                [min_ind_tmp,ind_ind_tmp] = min(ind_tmp);
                ind_max_signal(kk,ll,mm) = min_ind_tmp+8;
                max_signal(kk,ll,mm) = focal_signal_total(ind_max_signal(kk,ll,mm),kk,ll,mm);
                ii = ind_max_signal(kk,ll,mm);
                %% �Q�Ɠ_�Ɣg�ʂƂ̌`�󑊊ւ����߂�D�i���ݑ��֖@�j
                target_element = find((-focal_depth(1,ii)/2+focal_point(1,ii)<=t_pos(1,1:100)&(t_pos(1,1:100)<=focal_depth(1,ii)/2+focal_point(1,ii))));
                %��M�p�̎Q�Ɠ_�Z�o
                for jj = 1:num_echo_receiver
                    distance_from_focal_point_all(1,jj) = norm(t_pos(:,jj) - focal_point(:,ii));
                    delay_time_all = round(((distance_from_focal_point_all - focal_depth(1,ii))/v_reference(1,kk))/kgrid.dt);%[sample]
                    reference_point(1,jj) = round(delay_time_all(1,jj)+(2*focal_depth(1,ii)/v_reference(1,kk))/kgrid.dt+25-1);
                    %25��focal_amplitude���ő�ɂ���I�t�Z�b�g�D
                end
                max_delay = abs(diff(reference_point)) + 2;
                %���M�r�[���t�H�[�~���O�i���ʁj
                focused_rfdata = transmit_focus(rfdata,target_element,num_sample,delay_time_all,num_echo_receiver);
                %RF�f�[�^�}�X�L���O�i�ώ����]���̂���)
                min_reference = min(reference_point) - 5;
                max_reference = min(reference_point) + 40 - 1;
                mask_rfdata = zeros(size(focused_rfdata));
                mask_rfdata(min_reference:max_reference,:) = 1;
                focused_rfdata_masked = mask_rfdata .* focused_rfdata;
                for jj = 1:length(target_element)-1
                    %�g�ʌ`��]��
                    displacement(target_element(jj),kk,ll,mm)  = - finddelay(focused_rfdata(:,target_element(jj+1)),...
                        focused_rfdata_masked(:,target_element(jj)),max_delay(target_element(jj)));
                end
                displacement_ref(target_element(1:end-1),kk,ll,mm) = diff(reference_point(target_element));
                error_wavefront(kk,ll,mm) = sum(abs(displacement(target_element(1:end-1),kk,ll,mm) - ...
                    displacement_ref(target_element(1:end-1),kk,ll,mm)))/(length(target_element)-1);
                %% �o�̓f�[�^�̕ۑ�
                dst_path = sprintf('H:/result/2018_12_03_DAS_correct/2018_12_04_case26/lateral%0.1fmm/true%d_assumption%d',...
                    focal_point(1,9)*1e3,rate_IMCL(1,ll),rate_IMCL(1,kk));
                if ~exist(dst_path, 'dir')
                    mkdir(dst_path);
                end
                displacement_tmp = displacement(:,kk,ll,mm);
                save([dst_path,'\pre_rfdata.mat'],'focused_rfdata','focused_rfdata_masked','displacement_tmp','reference_point');
                figure;
                imagesc(focused_rfdata_masked);
                hold on
                scatter(min(target_element):max(target_element),reference_point(min(target_element):max(target_element)),'red');
                hold off
                xlabel('receiver[ch]');
                ylabel('time[sample]');
                axis square;
                axis tight;
                xlim([min(target_element) max(target_element)])
                ylim([min_reference max_reference])
                colormap(bone);
                colorbar;
                caxis([min(min(focused_rfdata_masked)) max(max(focused_rfdata_masked))])
                savefig([dst_path,'\pre_rfdata_masked.fig'])
                exportfig([dst_path,'\pre_rfdata_masked'],'png',[400,400])
                close gcf
                figure;
                imagesc(focused_rfdata);
                hold on
                scatter(min(target_element):max(target_element),reference_point(min(target_element):max(target_element)),2,'red');
                scatter(min(target_element):max(target_element),point_max_in_mask(min(target_element):max(target_element)),2,'blue','filled');
                hold off
                xlabel('receiver[ch]');
                ylabel('time[sample]');
                axis square;
                axis tight;
                colormap(bone);
                colorbar;
                caxis([min(min(focused_rfdata_masked)) max(max(focused_rfdata_masked))])
                savefig([dst_path,'\pre_rfdata_whole.fig'])
                exportfig([dst_path,'\pre_rfdata_whole'],'png',[400,400])
                close gcf
                figure;
                plot(focal_depth*1e3,focal_signal_total(:,kk,ll,mm));
                hold on
                scatter(focal_depth(ind_max_signal(kk,ll,mm))*1e3,max_signal(kk,ll,mm),'red','filled')
                hold off
                xlabel('�œ_�[��[mm]')
                ylabel('�M�����x[au]')
                savefig([dst_path,'\pre_focal_signal.fig'])
                exportfig([dst_path,'\pre_focal_signal'],'png',[400,400])
                close gcf
                fprintf('prepareing...lateral number is %d, medium number is %d, and estimation number is %d\n',mm,ll,kk);
            end
        end
        %% ���E�ʒu����肷��D(���E�ʒu���Fnum_medium)
        [~,ind_ind_max_signal] = max(max_signal(:,ll,mm));
        detected_boundary(:,ll,mm) =  ind_max_signal(ind_ind_max_signal,ll,mm);
        %% �e�}��(num_medium)�ɑ΂��ē���̏œ_�ʒu�ł�RF�f�[�^����g�`�`�󓙂�]������D
        for kk = 1:num_rate_IMCL
                %% ���ˋ��x�v���t�@�C�������߂�D
            v_reference(1,kk) = v_muscle*(1-rate_IMCL(1,kk)/100) + v_fat*(rate_IMCL(1,kk)/100);
            for ii = detected_boundary(1,ll,mm)
                target_element = find((-focal_depth(1,ii)/2+focal_point(1,ii)<=t_pos(1,1:100)&(t_pos(1,1:100)<=focal_depth(1,ii)/2+focal_point(1,ii))));
                %��M�p�̎Q�Ɠ_�Z�o
                for jj = 1:num_echo_receiver
                    distance_from_focal_point_all(1,jj) = norm(t_pos(:,jj) - focal_point(:,ii));
                    delay_time_all = round(((distance_from_focal_point_all - focal_depth(1,ii))/v_reference(1,kk))/kgrid.dt);%[sample]
                    reference_point(1,jj) = round(delay_time_all(1,jj)+(2*focal_depth(1,ii)/v_reference(1,kk))/kgrid.dt+25-1);
                end
                %���M�r�[���t�H�[�~���O�i���ʁj
                focused_rfdata = transmit_focus(rfdata,target_element,num_sample,delay_time_all,num_echo_receiver);
                %��M�r�[���t�H�[�J�V���O
                focal_signal(kk,ll,mm) =  receive_focus(focused_rfdata,target_element,reference_point);
            end
                %% �Q�Ɠ_�Ɣg�ʂƂ̌`�󑊊ւ����߂�D�i���ݑ��֖@�j
            %RF�f�[�^�}�X�L���O�i�ώ����]���̂���)
            min_reference = min(reference_point) - 5;
            max_reference = min(reference_point) + 40 - 1;
            mask_rfdata = zeros(size(focused_rfdata));
            mask_rfdata(min_reference:max_reference,:) = 1;
            focused_rfdata_masked = mask_rfdata .* focused_rfdata;
            for jj = 1:length(target_element)-1
                %�g�ʌ`��]��
                displacement(target_element(jj),kk,ll,mm)  = - finddelay(focused_rfdata(:,target_element(jj+1)),...
                    focused_rfdata_masked(:,target_element(jj)),max_delay(target_element(jj)));
            end
            displacement_ref(target_element(1:end-1),kk,ll,mm) = diff(reference_point(target_element));
            error_wavefront(kk,ll,mm) = sum(abs(displacement(target_element(1:end-1),kk,ll,mm) - ...
                displacement_ref(target_element(1:end-1),kk,ll,mm)))/(length(target_element)-1);
            dst_path = sprintf('H:/result/2018_12_03_DAS_correct/2018_12_04_case26/lateral%0.1fmm/true%d_assumption%d',...
                focal_point(1,9)*1e3,rate_IMCL(1,ll),rate_IMCL(1,kk));
            if ~exist(dst_path, 'dir')
                mkdir(dst_path);
            end
            displacement_tmp = displacement(:,kk,ll,mm);
            save([dst_path,'\rfdata.mat'],'focused_rfdata','focused_rfdata_masked','displacement_tmp','reference_point');
            figure;
            imagesc(focused_rfdata_masked);
            hold on
            scatter(min(target_element):max(target_element),reference_point(min(target_element):max(target_element)),'red');
            hold off
            xlabel('receiver[ch]');
            ylabel('time[sample]');
            axis square;
            axis tight;
            xlim([min(target_element) max(target_element)])
            ylim([min_reference max_reference])
            colormap(bone);
            colorbar;
            caxis([min(min(focused_rfdata_masked)) max(max(focused_rfdata_masked))])
            savefig([dst_path,'\rfdata_masked.fig'])
            exportfig([dst_path,'\rfdata_masked'],'png',[400,400])
            close gcf
            figure;
            imagesc(focused_rfdata);
            hold on
            scatter(min(target_element):max(target_element),reference_point(min(target_element):max(target_element)),2,'red');
            scatter(min(target_element):max(target_element),point_max_in_mask(min(target_element):max(target_element)),2,'blue','filled');
            hold off
            xlabel('receiver[ch]');
            ylabel('time[sample]');
            axis square;
            axis tight;
            colormap(bone);
            colorbar;
            caxis([min(min(focused_rfdata_masked)) max(max(focused_rfdata_masked))])
            savefig([dst_path,'\rfdata_whole.fig'])
            exportfig([dst_path,'\rfdata_whole'],'png',[400,400])
            close gcf
            fprintf('detected... lateral number is %d, medium number is %d, and estimation number is %d\n',mm,ll,kk);
        end
    end
    %% �o�̓f�[�^�̕ۑ�
    max_signal_normalized = zeros(num_rate_IMCL,num_medium,num_lateral);
    for ll = 1:num_medium
        tmp2 = focal_signal(:,ll,mm);
        max_signal_normalized(:,ll,mm) = tmp2/max(tmp2);
    end
    
    dst_path = sprintf('H:/result/2018_12_03_DAS_correct/2018_12_04_case26/lateral%0.1fmm/',focal_point(1,9)*1e3);
    save([dst_path,'\result.mat'],'ind_max_signal','focal_depth','displacement',...
        'reference_point','focal_signal_total','focal_signal','max_signal','max_signal_normalized');
    
    
    figure;
    tmp = zeros(num_rate_IMCL,num_medium);
    for ii = 1:num_rate_IMCL
        tmp(ii,:) = focal_depth(1,ind_max_signal(ii,:,mm))*1e3;
    end
    heatmap(tmp);
    title('boundary depth[mm]')
    xlabel('correct IMCL percentage [%]');
    ylabel('predicted IMCL percentage [%]');
    exportfig([dst_path,'\pre_boundary_eng'],'png',[300,250]);
    close gcf
    
    figure;
    tmp = zeros(num_rate_IMCL,num_medium);
    for ii = 1:num_rate_IMCL
        tmp(ii,:) = focal_depth(1,ind_max_signal(ii,:,mm))*1e3;
    end
    heatmap(tmp);
    title('�}�����E�[�x[mm]')
    xlabel('����IMCL��L�� [%]');
    ylabel('�\��IMCL��L�� [%]');
    exportfig([dst_path,'\pre_boundary'],'png',[300,250]);
    close gcf
    
        figure;
    tmp = zeros(num_rate_IMCL,num_medium);
    for ii = 1:num_rate_IMCL
        tmp(ii,:) = focal_depth(1,detected_boundary(ii,:,mm))*1e3;
    end
    heatmap(tmp);
    title('boundary depth[mm]')
    xlabel('correct IMCL percentage [%]');
    ylabel('predicted IMCL percentage [%]');
    exportfig([dst_path,'\detected_boundary_eng'],'png',[300,250]);
    close gcf
    
    figure;
    tmp = zeros(num_rate_IMCL,num_medium);
    for ii = 1:num_rate_IMCL
        tmp(ii,:) = focal_depth(1,detected_boundary(ii,:,mm))*1e3;
    end
    heatmap(tmp);
    title('�}�����E�[�x[mm]')
    xlabel('����IMCL��L�� [%]');
    ylabel('�\��IMCL��L�� [%]');
    exportfig([dst_path,'\detected_boundary'],'png',[300,250]);
    close gcf
    
    figure;
    heatmap(error_wavefront(:,:,mm));
    title('wavefront distortion');
    xlabel('correct IMCL percentage [%]');
    ylabel('predicted IMCL percentage [%]')
    exportfig([dst_path,'\detected_wavefront_eng'],'png',[300,250]);
    close gcf
    
    figure;
    heatmap(error_wavefront(:,:,mm));
    title('�g�ʘc�ݎw�W');
    xlabel('����IMCL��L�� [%]');
    ylabel('�\��IMCL��L�� [%]');
    exportfig([dst_path,'\detected_wavefront]]'],'png',[300,250]);
    close gcf
    
    figure;
    heatmap(max_signal(:,:,mm));
    title('�M�����x[au]');
    xlabel('����IMCL��L�� [%]');
    ylabel('�\��IMCL��L�� [%]');
    exportfig([dst_path,'\pre_intensity'],'png',[300,250]);
    close gcf
    
    figure;
    heatmap(max_signal(:,:,mm));
    title('Intensity[au]');
    xlabel('correct IMCL percentage [%]');
    ylabel('predicted IMCL percentage [%]');
    exportfig([dst_path,'\pre_intensity_eng'],'png',[300,250]);
    close gcf
    
    figure;
    heatmap(focal_signal(:,:,mm));
    title('�M�����x[au]');
    xlabel('����IMCL��L�� [%]');
    ylabel('�\��IMCL��L�� [%]');
    exportfig([dst_path,'\detected_intensity'],'png',[300,250]);
    close gcf
    
    figure;
    heatmap(focal_signal(:,:,mm));
    title('Intensity[au]');
    xlabel('correct IMCL percentage [%]');
    ylabel('predicted IMCL percentage [%]');
    exportfig([dst_path,'\detected_intensity_eng'],'png',[300,250]);
    close gcf    
end
%% ���ׂĂ̕ϐ���ۑ�
savefilename = sprintf('H:/result/2018_12_03_DAS_correct/2018_12_04_case26/total_result.mat');
save(savefilename);