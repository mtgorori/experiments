%%%%%%%%%%%%%%%%%%%%%%%%%%%%��������������%%%%%
% ���E�ʂ����o���Q�Ɠ_�̑Ó����]��
% ch 50�����̃s�N�Z���Ƀt�H�[�J�X��������D
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
load("H:/data/kwave/config/t_pos_2board.mat");
load("H:/data/kwave/medium/2018_11_07_layer_medium/reference.mat")
pathname = sprintf('H:/data/kwave/result/2018_11_07_layer_medium/Layer_medium_boundary_7.9mm_ IMCL%d%%',1);
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
point_maxAmp_in_mask = zeros(1,num_echo_receiver);%�}�X�N�������RF�f�[�^�ŐU���ő�̃T���v���_���
distance_from_focal_point_all = zeros(1,num_echo_receiver);

num_depth = (t_pos(2,1)-t_pos(2,101))/kgrid.dx - 3;%'3'�Ƃ���̂́C�ŋߐڋ�����0.4 mm�ł��邱�Ƃ��l�����Ă���D
focal_depth = zeros(1,num_depth);
focal_point = zeros(2,num_depth);
for ii = 1:num_depth
    focal_depth(1,ii) = (ii+3)*kgrid.dx;
    focal_point(2,ii) = t_pos(2,1)-focal_depth(1,ii);
    focal_point(1,ii) = 0;
end

focal_signal_total = zeros(num_rate_IMCL,num_depth);
ind_max_signal = zeros(num_rate_IMCL,num_medium);
homogeneity_percel = zeros(num_rate_IMCL,num_echo_receiver);
homogeneity_total = zeros(num_rate_IMCL,num_medium);

for ll = 1:num_medium
    pathname = sprintf('H:/data/kwave/result/2018_11_07_layer_medium/Layer_medium_boundary_7.9mm_ IMCL%d%%',ll);
    cd(pathname);
    load('rfdata.mat');
    load('kgrid.mat');
    [num_sample,~,~] = size(rfdata);
    focused_rfdata = zeros(num_sample,num_echo_receiver);
    focused_rfdata_masked = zeros(num_sample,num_echo_receiver);
    focused_rfdata_amp = zeros(num_sample,num_echo_receiver);
    focused_rfdata_amp_masked = zeros(num_sample,num_echo_receiver);
    for kk = 1:num_rate_IMCL
        %% ���ˋ��x�v���t�@�C�������߂�D
        v_reference(1,kk) = v_muscle*(1-rate_IMCL(1,kk)/100) + v_fat*(rate_IMCL(1,kk)/100);
        for ii = 1:num_depth
            target_element = find((-focal_depth(1,ii)/2<=t_pos(1,1:100)&(t_pos(1,1:100)<=focal_depth(1,ii)/2)));
            distance_from_focal_point = zeros(1,length(target_element));
            %��M�p�̎Q�Ɠ_�Z�o
            for jj = 1:num_echo_receiver
                distance_from_focal_point_all(1,jj) = norm(t_pos(:,jj) - focal_point(:,ii));
                delay_time_all = round(((distance_from_focal_point_all - focal_depth(1,ii))/v_reference(1,kk))/kgrid.dt);%[sample]
                reference_point(1,jj) = round(delay_time_all(1,jj)+1+(2*focal_depth(1,ii)/v_reference(1,kk))/kgrid.dt+25);
            end
            %���M�r�[���t�H�[�~���O�i���ʁj
            for jj = 1:length(target_element)
                distance_from_focal_point(1,jj) = distance_from_focal_point_all(1,target_element(jj));
                % �x������
                delay_time = delay_time_all(1,target_element(jj));%[sample]
                read_range_rfdata = length(delay_time+1:num_sample);
                focused_rfdata(1:read_range_rfdata,:) = focused_rfdata(1:read_range_rfdata,:)...
                    +  rfdata(delay_time+1:num_sample,1:100,target_element(jj));%�������Z
            end
            for jj = 1:length(target_element)
                %��M�r�[���t�H�[�~���O�i�������Z�̂��߁j
                focal_signal_total(kk,ii) = focal_signal_total(kk,ii)+ ...
                    focused_rfdata(reference_point(1,target_element(1,jj)),target_element(1,jj))/length(target_element);
            end
        end
        %% ���ˋ��x���ő�̏œ_�ʒu��T���D
        [~,ind_max_signal(kk,ll)] = max(abs(focal_signal_total(kk,:)));
        ii = ind_max_signal(kk,ll);
        %% �Q�Ɠ_�Ɣg�ʂƂ̌덷��]���D
        target_element = find((-focal_depth(1,ii)/2<=t_pos(1,1:100)&(t_pos(1,1:100)<=focal_depth(1,ii)/2)));
        distance_from_focal_point = zeros(1,length(target_element));
        %��M�p�̎Q�Ɠ_�Z�o
        for jj = 1:num_echo_receiver
            distance_from_focal_point_all(1,jj) = norm(t_pos(:,jj) - focal_point(:,ii));
            delay_time_all = round(((distance_from_focal_point_all - focal_depth(1,ii))/v_reference(1,kk))/kgrid.dt);%[sample]
            reference_point(1,jj) = round(delay_time_all(1,jj)+1+(2*focal_depth(1,ii)/v_reference(1,kk))/kgrid.dt+25);
            %25��focal_amplitude���ő�ɂ���I�t�Z�b�g�D
            reference_point_lowerlimit(1,jj) ...
                = round(delay_time_all(1,jj)*(v_reference(1,kk)/v_muscle)+1+(2*focal_depth(1,ii)/v_muscle)/kgrid.dt+25-1);
            reference_point_upperlimit(1,jj) ...
                = round(delay_time_all(1,jj)*(v_reference(1,kk)/v_fat)+1+(2*focal_depth(1,ii)/v_fat)/kgrid.dt+25);
            %�ǂ�Ȃɒx�����Ă��������B���Ă����͈͓̔��ɏœ_�ʒu����̃G�R�[�p���X�������Ă���ł��낤����E����
        end
        %���M�r�[���t�H�[�~���O�i���ʁj
        for jj = 1:length(target_element)
            distance_from_focal_point(1,jj) = distance_from_focal_point_all(1,target_element(jj));
            % �x������
            delay_time = delay_time_all(1,target_element(jj));%[sample]
            read_range_rfdata = length(delay_time+1:num_sample);
            focused_rfdata(1:read_range_rfdata,:) = focused_rfdata(1:read_range_rfdata,:)...
                +  rfdata(delay_time+1:num_sample,1:100,target_element(jj));%�������Z
        end
        hilb_rfdata = hilbert(focused_rfdata);
        focused_rfdata_amp = abs(hilb_rfdata);
        focused_rfdata_amp_masked = focused_rfdata_amp;
        %RF�f�[�^�}�X�L���O�i�ώ����]���̂��߁j
        for jj = 1:num_echo_receiver
            focused_rfdata_amp_masked(1:reference_point_lowerlimit(1,jj),jj) = 0;
            focused_rfdata_amp_masked(reference_point_upperlimit(1,jj):end,jj) = 0;
        end
        [~,point_maxAmp_in_mask] = max(focused_rfdata_amp_masked,[],1);
        for jj = 1:length(target_element)
            %�ώ����]���w�W
            homogeneity_percel(kk,jj) = abs(point_maxAmp_in_mask(1,target_element(jj)) - reference_point(1,target_element(jj)));
        end
        homogeneity_total(kk,ll) = sum(homogeneity_percel(kk,:))/length(target_element);
    end
end