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
t_facing_distance = 0.04;%[m]
rate_IMCL = [0, 2, 4, 6, 8];
num_depth = (t_pos(2,1)-t_pos(2,101))/kgrid.dx - 3;%'3'�Ƃ���̂́C�ŋߐڋ�����0.4 mm�ł��邱�Ƃ��l�����Ă���D
[num_sample,num_receiver,num_transmitter] = size(rfdata);
num_echo_receiver = num_transmitter;
reference_point = zeros(2,num_echo_receiver,num_depth);
reference_point_base = linspace(1,100,100)';
reference_point_base = repmat(reference_point_base,1,num_depth);
reference_point(1,:,:) = reference_point_base;
distance_from_focal_point_all = zeros(1,num_echo_receiver);
focal_signal = zeros(length(rate_IMCL),num_depth);
focal_amp = zeros(length(rate_IMCL),num_depth,num_echo_receiver);
focal_phase = zeros(length(rate_IMCL),num_depth,num_echo_receiver);
for kk = 1:length(rate_IMCL)
    v_reference = v_muscle*(1-rate_IMCL(kk)/100) + v_fat*(rate_IMCL(kk)/100);
    focused_rfdata = zeros(num_sample,num_echo_receiver,num_depth);
    focused_rfdata_amp = zeros(num_sample,num_echo_receiver,num_depth);
    focused_rfdata_phase = zeros(num_sample,num_echo_receiver,num_depth);
    for ii = 1:num_depth
        focal_depth = (ii+3)*kgrid.dx;
        focal_point = [0;t_pos(2,1)-focal_depth];
        target_element = find((-focal_depth/2<=t_pos(1,1:100)&(t_pos(1,1:100)<=focal_depth/2)));
        distance_from_focal_point = zeros(1,length(target_element));
        for jj = 1:num_echo_receiver
            distance_from_focal_point_all(1,jj) = norm(t_pos(:,jj) - focal_point);
            delay_time_all = round(((distance_from_focal_point_all - focal_depth)/v_reference)/kgrid.dt);%[sample]
            reference_point(2,jj,ii) = round(delay_time_all(1,jj)+1+(2*focal_depth/v_reference)/kgrid.dt+25);
            %25��focal_amplitude���ő�ɂ���I�t�Z�b�g�D
        end
        for jj = 1:length(target_element)
            distance_from_focal_point(1,jj) = distance_from_focal_point_all(1,target_element(jj));
            % �x������
            delay_time = delay_time_all(1,target_element(jj));%[sample]
            read_range_rfdata = length(delay_time+1:num_sample);
            focused_rfdata(1:read_range_rfdata,:,ii) = focused_rfdata(1:read_range_rfdata,:,ii)...
                +  rfdata(delay_time+1:num_sample,1:100,target_element(jj));%�������Z
        end
        hilb_rfdata = hilbert(focused_rfdata(:,:,ii));
        focused_rfdata_amp(:,:,ii) = abs(hilb_rfdata);
        focused_rfdata_phase(:,:,ii) = atan(imag(hilb_rfdata)./real(hilb_rfdata));
        for jj = 1:length(target_element)
            tmp = focused_rfdata(reference_point(2,target_element(1,jj),ii),target_element(1,jj),ii)/length(target_element);
            focal_signal(kk,ii) = focal_signal(kk,ii)+ tmp;
            focal_amp(kk,ii,target_element(jj)) = focused_rfdata_amp(reference_point(2,target_element(1,jj),ii),target_element(1,jj),ii)/length(target_element);
            focal_phase(kk,ii,target_element(jj)) = focused_rfdata_phase(reference_point(2,target_element(1,jj),ii),target_element(1,jj),ii);
        end
        focal_signal(kk,ii) = abs(hilbert(focal_signal(kk,ii)));
    end
end