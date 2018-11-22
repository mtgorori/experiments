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

clear;
load("H:/data/kwave/config/t_pos_2board.mat");
pathname = sprintf('H:/data/kwave/result/2018_11_11_case26_variousIMCL/case26_IMCL%0.1f',1.0);
cd(pathname);
load('rfdata.mat');
load('kgrid.mat');
% load("H:/data/kwave/result/2018_10_21_point_medium/point_mudium5MHz/rfdata.mat");
% load("H:/data/kwave/result/2018_10_21_point_medium/point_mudium5MHz/medium.mat");
% load("H:/data/kwave/result/2018_10_21_point_medium/point_mudium5MHz/kgrid.mat");

% �����ݒ�
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

num_depth = (t_pos(2,1)-t_pos(2,101))/kgrid.dx - 3;%'3'�Ƃ���̂́C�ŋߐڋ�����0.4 mm�ł��邱�Ƃ��l�����Ă���D
focal_depth = zeros(1,num_depth);
focal_point = zeros(2,num_depth);
for ii = 1:num_depth
    focal_depth(1,ii) = (ii+3)*kgrid.dx;
    focal_point(2,ii) = t_pos(2,1)-focal_depth(1,ii);
    focal_point(1,ii) = 0;
end


ind_max_signal = zeros(num_rate_IMCL,num_medium);
max_signal = zeros(num_rate_IMCL,num_medium);
homogeneity_percel = zeros(num_rate_IMCL,num_echo_receiver);
homogeneity_total = zeros(num_rate_IMCL,num_medium);
focal_signal_total = zeros(num_depth,num_rate_IMCL,num_medium);
for ll = 1
    pathname = sprintf('H:/data/kwave/result/2018_11_07_layer_medium/Layer_medium_boundary_7.9mm_ IMCL%d%%',ll);
    cd(pathname);
    load('rfdata.mat');
    load('kgrid.mat');
    [num_sample,~,~] = size(rfdata);
    for kk = 1
        %% ���ˋ��x�v���t�@�C�������߂�D
        v_reference(1,kk) = v_muscle*(1-rate_IMCL(1,kk)/100) + v_fat*(rate_IMCL(1,kk)/100);
        for ii = 1:num_depth
            focused_rfdata = zeros(num_sample,num_echo_receiver);
            target_element = find((-focal_depth(1,ii)/2<=t_pos(1,1:100)&(t_pos(1,1:100)<=focal_depth(1,ii)/2)));
            %��M�p�̎Q�Ɠ_�Z�o
            for jj = 1:num_echo_receiver
                distance_from_focal_point_all(1,jj) = norm(t_pos(:,jj) - focal_point(:,ii));
                delay_time_all = round(((distance_from_focal_point_all - focal_depth(1,ii))/v_reference(1,kk))/kgrid.dt);%[sample]
                reference_point(1,jj) = round(delay_time_all(1,jj)+(2*focal_depth(1,ii)/v_reference(1,kk))/kgrid.dt+25-1);
            end
            %���M�r�[���t�H�[�~���O�i���ʁj
            distance_from_focal_point = distance_from_focal_point_all(target_element);
            % �x������
            delay_time = delay_time_all(1,target_element);%[sample]
            for jj = 1:length(target_element)
                read_range_rfdata = length(delay_time+1:num_sample);
                focused_rfdata(1:read_range_rfdata,:) = focused_rfdata(1:read_range_rfdata,:)...
                    +  rfdata(delay_time+1:num_sample,1:100,target_element(jj));%�������Z
            end
            for jj = 1:length(target_element)
                %��M�r�[���t�H�[�~���O�i�������Z�̂��߁j
                focal_signal_total(ii,kk,ll) = focal_signal_total(ii,kk,ll)+ ...
                    focused_rfdata(reference_point(1,target_element(1,jj)),target_element(1,jj))/length(target_element);
            end
        end
        %% ���ˋ��x���ő�̏œ_�ʒu��T���D
        [max_tmp,ind_tmp] = findpeaks(focal_signal_total(:,kk,ll),'Npeaks',1,'SortStr','descend');
        ind_max_signal(kk,ll) = ind_tmp;
        max_signal(kk,ll) = max_tmp;
        ii = ind_tmp;
        %% �Q�Ɠ_�Ɣg�ʂƂ̌덷��]���D
        focused_rfdata = zeros(num_sample,num_echo_receiver);
        target_element = find((-focal_depth(1,ii)/2<=t_pos(1,1:100)&(t_pos(1,1:100)<=focal_depth(1,ii)/2)));
        %��M�p�̎Q�Ɠ_�Z�o
        for jj = 1:num_echo_receiver
            distance_from_focal_point_all(1,jj) = norm(t_pos(:,jj) - focal_point(:,ii));
            delay_time_all = round(((distance_from_focal_point_all - focal_depth(1,ii))/v_reference(1,kk))/kgrid.dt);%[sample]
            reference_point(1,jj) = round(delay_time_all(1,jj)+(2*focal_depth(1,ii)/v_reference(1,kk))/kgrid.dt+25-1);
            %25��focal_amplitude���ő�ɂ���I�t�Z�b�g�D
            reference_point_lowerlimit(1,jj) ...
                = round(delay_time_all(1,jj)*(v_reference(1,kk)/v_muscle)+(2*focal_depth(1,ii)/v_muscle)/kgrid.dt+25-1-1);
            reference_point_upperlimit(1,jj) ...
                = round(delay_time_all(1,jj)*(v_reference(1,kk)/v_fat)+(2*focal_depth(1,ii)/v_fat)/kgrid.dt+25-1);
            %�ǂ�Ȃɒx�����Ă��������B���Ă����͈͓̔��ɏœ_�ʒu����̃G�R�[�p���X�������Ă���ł��낤����E����
        end
        %���M�r�[���t�H�[�~���O�i���ʁj
        distance_from_focal_point = distance_from_focal_point_all(1,target_element);
        % �x������
        delay_time = delay_time_all(target_element);%[sample]
        read_range_rfdata = num_sample-delay_time;
        focused_rfdata(1:read_range_rfdata,:) = sum(rfdata(delay_time+1:num_sample,1:100,target_element),3);%�������Z
        focused_rfdata_masked = focused_rfdata;
        %RF�f�[�^�}�X�L���O�i�ώ����]���̂��߁j
        for jj = 1:num_echo_receiver
            focused_rfdata_masked(1:reference_point_lowerlimit,:) = 0;
            focused_rfdata_masked(reference_point_upperlimit(1,jj):end,jj) = 0;
        end
        [~,point_max_in_mask] = max(focused_rfdata_masked,[],1);
        for jj = 1:length(target_element)
            %�ώ����]���w�W
            homogeneity_percel(kk,jj) = abs(point_max_in_mask(1,target_element(jj)) - reference_point(1,target_element(jj)));
        end
        homogeneity_total(kk,ll) = sum(homogeneity_percel(kk,:))/length(target_element);
        dst_path = sprintf('H:/result/2018_11_07_IMCL_estimation_principle_verifiacation/2018_11_11_maxrfsignal_detect_boundary_case26/true%d_assumption%d',...
            rate_IMCL(1,ll),rate_IMCL(1,kk));
        if ~exist(dst_path, 'dir')
            mkdir(dst_path);
        end
        save([dst_path,'\rfdata.mat'],'focused_rfdata','focused_rfdata_masked');
    end
end
% cd 'H:/result/2018_11_07_IMCL_estimation_principle_verifiacation'
% save('2018_11_10_with_rfsignal_homogeneity_SingleLayer_boundary_7.9mm_varied_focaldepth_20x20_ver2.mat','ind_tmp','ind_max_signal','focal_depth','homogeneity_total','homogeneity_percel',...
%     'reference_point','focused_rfdata','focused_rfdata_masked','focal_signal_total','max_signal');