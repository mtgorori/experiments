%%%%%%%%%%%%%%%%%%%%%%%%%%%%��������������%%%%%
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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
load("H:\data\kwave\config\t_pos_2board.mat");
load("H:\data\kwave\result\2018_10_21_highFrequency_variousIMCL\case26_IMCL4.0_5MHz\rfdata.mat");
load("H:\data\kwave\result\2018_10_21_highFrequency_variousIMCL\case26_IMCL4.0_5MHz\medium.mat");
load("H:\data\kwave\result\2018_10_21_highFrequency_variousIMCL\case26_IMCL4.0_5MHz\kgrid.mat");
load("H:\data\kwave\result\2018_10_21_highFrequency_variousIMCL\case26_IMCL4.0_5MHz\sourse_wave.mat");
% load("H:\data\kwave\result\2018_10_21_point_medium\point_mudium5MHz\rfdata.mat");
% load("H:\data\kwave\result\2018_10_21_point_medium\point_mudium5MHz\medium.mat");
% load("H:\data\kwave\result\2018_10_21_point_medium\point_mudium5MHz\kgrid.mat");
% load("H:\data\kwave\result\2018_10_21_point_medium\point_mudium5MHz\sourse_wave.mat");

% �����ݒ�
v_fat = 1450;%[m/s]
v_muscle = 1580;%[m/s]
rate_IMCL = [1,2,3,4,5,6,7,8,9,10];
v_reference = zeros(1,length(rate_IMCL));
t_facing_distance = 0.04;%[m]
num_depth = (t_pos(2,1)-t_pos(2,101))/kgrid.dx - 3;%'3'�Ƃ���̂́C�ŋߐڋ�����0.4 mm�ł��邱�Ƃ��l�����Ă���D
focal_depth = zeros(1,num_depth);
[num_sample,num_receiver,num_transmitter] = size(rfdata);
num_echo_receiver = num_transmitter;
reference_point = zeros(num_echo_receiver,num_depth,length(rate_IMCL));
reference_point_lowerlimit = zeros(num_echo_receiver,num_depth,length(rate_IMCL));%�ώ����]���̂��߂�RF�f�[�^�}�X�L���O�Ɏg��
reference_point_upperlimit = zeros(num_echo_receiver,num_depth,length(rate_IMCL));%�ώ����]���̂��߂�RF�f�[�^�}�X�L���O�Ɏg��
point_maxAmp_in_mask = zeros(num_echo_receiver,num_depth,length(rate_IMCL));%�}�X�N�������RF�f�[�^�ŐU���ő�̃T���v���_���
distance_from_focal_point_all = zeros(1,num_echo_receiver);
focal_signal_total = zeros(length(rate_IMCL),num_depth);
focal_signal_total2 = zeros(length(rate_IMCL),num_depth);
focal_amp = zeros(length(rate_IMCL),num_depth,num_echo_receiver);
% focal_phase = zeros(length(rate_IMCL),num_depth,num_echo_receiver);
homogeneity_percel = zeros(length(rate_IMCL),num_depth,num_echo_receiver);
homogeneity_total = zeros(length(rate_IMCL),num_depth);
focused_rfdata = zeros(num_sample,num_echo_receiver,num_depth,length(rate_IMCL));
focused_rfdata_masked = zeros(num_sample,num_echo_receiver,num_depth,length(rate_IMCL));
focused_rfdata_amp = zeros(num_sample,num_echo_receiver,num_depth,length(rate_IMCL));
focused_rfdata_amp_masked = zeros(num_sample,num_echo_receiver,num_depth,length(rate_IMCL));

for kk = 1:length(rate_IMCL)
    v_reference(1,kk) = v_muscle*(1-rate_IMCL(kk)/100) + v_fat*(rate_IMCL(kk)/100);
    for ii = 76
        focal_depth(1,ii) = (ii+3)*kgrid.dx;
        focal_point = [0;t_pos(2,1)-focal_depth(1,ii)];
        target_element = find((-focal_depth(1,ii)/2<=t_pos(1,1:100)&(t_pos(1,1:100)<=focal_depth(1,ii)/2)));
        distance_from_focal_point = zeros(1,length(target_element));
        %��M�p�̎Q�Ɠ_�Z�o
        for jj = 1:num_echo_receiver
            distance_from_focal_point_all(1,jj) = norm(t_pos(:,jj) - focal_point);
            delay_time_all = round(((distance_from_focal_point_all - focal_depth(1,ii))/v_reference(1,kk))/kgrid.dt);%[sample]
            reference_point(jj,ii,kk) = round(delay_time_all(1,jj)+1+(2*focal_depth(1,ii)/v_reference(1,kk))/kgrid.dt+25);
            %25��focal_amplitude���ő�ɂ���I�t�Z�b�g�D
            reference_point_lowerlimit(jj,ii,kk) ...
                = round(delay_time_all(1,jj)*(v_reference(1,kk)/v_muscle)+1+(2*focal_depth(1,ii)/v_muscle)/kgrid.dt+25-1);
            reference_point_upperlimit(jj,ii,kk) ...
                = round(delay_time_all(1,jj)*(v_reference(1,kk)/v_fat)+1+(2*focal_depth(1,ii)/v_fat)/kgrid.dt+25);
            %�ǂ�Ȃɒx�����Ă��������B���Ă����͈͓̔��ɏœ_�ʒu����̃G�R�[�p���X�������Ă���ł��낤����E����
        end
        %���M�r�[���t�H�[�~���O�i���ʁj
        for jj = 1:length(target_element)
            distance_from_focal_point(1,jj) = distance_from_focal_point_all(1,target_element(jj));
            % �x������
            delay_time = delay_time_all(1,target_element(jj));%[sample]
            read_range_rfdata = length(delay_time+1:num_sample);
            focused_rfdata(1:read_range_rfdata,:,ii,kk) = focused_rfdata(1:read_range_rfdata,:,ii,kk)...
                +  rfdata(delay_time+1:num_sample,1:100,target_element(jj));%�������Z
        end
        hilb_rfdata = hilbert(focused_rfdata(:,:,ii,kk));
        focused_rfdata_amp(:,:,ii,kk) = abs(hilb_rfdata);
        %         focused_rfdata_phase(:,:,ii) = atan(imag(hilb_rfdata)./real(hilb_rfdata));
        focused_rfdata_amp_masked(:,:,ii,kk) = focused_rfdata_amp(:,:,ii,kk);
        focused_rfdata_masked(:,:,ii,kk) = focused_rfdata(:,:,ii,kk);
        %RF�f�[�^�}�X�L���O�i�ώ����]���̂��߁j
        for jj = 1:num_echo_receiver
            focused_rfdata_amp_masked(1:reference_point_lowerlimit(jj,ii),jj,ii,kk) = 0;
            focused_rfdata_amp_masked(reference_point_upperlimit(jj,ii):end,jj,ii,kk) = 0;
            focused_rfdata_masked(1:reference_point_lowerlimit(jj,ii),jj,ii,kk) = 0;
            focused_rfdata_masked(reference_point_upperlimit(jj,ii):end,jj,ii,kk) = 0;
        end
        [~,point_maxAmp_in_mask(:,ii,kk)] = max(focused_rfdata_amp_masked(:,:,ii,kk),[],1);
        for jj = 1:length(target_element)
            %��M�r�[���t�H�[�~���O�i�������Z�̂��߁j
            focal_signal_total(kk,ii) = focal_signal_total(kk,ii)+ ...
                focused_rfdata_amp(reference_point(target_element(1,jj),ii),target_element(1,jj),ii,kk)/length(target_element);
            focal_signal_total2(kk,ii) = focal_signal_total2(kk,ii)+ ...
                focused_rfdata(reference_point(target_element(1,jj),ii),target_element(1,jj),ii,kk)/length(target_element);
            focal_amp(kk,ii,target_element(jj)) = ...
                focused_rfdata_amp(reference_point(target_element(1,jj),ii),target_element(1,jj),ii,kk)/length(target_element);
            %�ώ����]���w�W
            homogeneity_percel(kk,ii,jj) = abs(point_maxAmp_in_mask(jj,ii,kk) - reference_point(jj,ii));
        end
        homogeneity_total(kk,ii) = sum(homogeneity_percel(kk,ii,:))/length(target_element);
    end
end