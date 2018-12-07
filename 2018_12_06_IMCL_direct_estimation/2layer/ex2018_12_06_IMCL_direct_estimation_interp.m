%%%%%%%%%%%%%%%%%%%%
% �ΏہF�񑊔}���C���E�����F15.0mm
% �ݒ艹���F1580[m/s]����������
% �œ_�����ʒu�Œ�Fy=0
% ���E�ʒu�F���m
% IMCL������0 %�ɌŒ肷��D
% ��M�g�ʂ̋ȗ��ɂ�艹������
% �f�q�Ԏ�M�g���ԍ��͎��M���̐��K�����ݑ��ւ�p����
%%%%%%%%%%%%%%%%%%%%

%% �����ݒ�
clear
load("H:/data/kwave/config/t_pos_2board.mat");
load("H:\data\kwave\result\2018_11_07_layer_medium\Layer_medium_boundary_15.0mm_ IMCL0%\sourse_wave.mat")
load("H:\data\kwave\result\2018_11_07_layer_medium\Layer_medium_boundary_15.0mm_ IMCL0%\rfdata.mat")
load("H:\data\kwave\medium\2018_11_07_layer_medium\Layer_medium_boundary_15.0mm_ IMCL0%.mat")
load("H:\data\kwave\result\2018_11_07_layer_medium\Layer_medium_boundary_15.0mm_ IMCL0%\kgrid.mat")
v_fat = 1450;%[m/s]
v_muscle = 1580;%[m/s]
t_facing_distance = 0.04;%[m]
[num_sample,num_receiver,num_transmitter] = size(rfdata);
num_echo_receiver = num_transmitter;
reference_point = zeros(1,num_echo_receiver);
reference_point_lowerlimit = zeros(1,num_echo_receiver);%�ώ����]���̂��߂�RF�f�[�^�}�X�L���O�Ɏg��
reference_point_upperlimit = zeros(1,num_echo_receiver);%�ώ����]���̂��߂�RF�f�[�^�}�X�L���O�Ɏg��
distance_from_focal_point_all = zeros(1,num_echo_receiver);
num_depth = (t_pos(2,1)-t_pos(2,101))/kgrid.dx/2 - 3;%'3'�Ƃ���̂́C�ŋߐڋ�����0.4 mm�ł��邱�Ƃ��l�����Ă���D
focal_depth = medium_focal_depth(14);
focal_point = [0 ;kgrid.x_vec(ind_focal_point(14))];
element_pitch = abs(t_pos(1,1) - t_pos(1,2));
lateral_range_max = 0;
lateral_range_min = 0;
num_lateral = round((lateral_range_max -lateral_range_min) / element_pitch)+1;%0.0~4.8 mm�܂�
ind_max_signal = 0;
max_signal = 0;
displacement = zeros(1,num_echo_receiver);
displacement_ref = zeros(1,num_echo_receiver);
error_wavefront = 0;
% focal_signal_total = zeros(1,num_depth);
detected_boundary = 0;
estimated_velocity  = 0;

%% �t�H�[�J�V���O
%�쓮�f�q�̌���
v_reference = v_muscle;
target_element = find((-focal_depth/2<=t_pos(1,1:100)&(t_pos(1,1:100)<=focal_depth/2)));
%��M�p�̎Q�Ɠ_�Z�o
for jj = 1:num_echo_receiver
    distance_from_focal_point_all(1,jj) = norm(t_pos(:,jj) - focal_point);
    delay_time_all = floor(((distance_from_focal_point_all - focal_depth)/v_reference)/kgrid.dt);%[sample]
    reference_point(1,jj) = floor(delay_time_all(1,jj)+(2*focal_depth/v_reference)/kgrid.dt+25-1);
end
%���M�r�[���t�H�[�~���O�i���ʁj
focused_rfdata = transmit_focus(rfdata,target_element,num_sample,delay_time_all,num_echo_receiver);
%��M�r�[���t�H�[�J�V���O
focal_signal_total =  receive_focus(focused_rfdata,target_element,reference_point);
interp_rfdata = zeros(num_sample*4,num_echo_receiver);
for jj = 1:num_echo_receiver
    interp_rfdata(:,jj) = interp(focused_rfdata(:,jj),4);
end
%% �ȗ��v�Z
%RF�f�[�^�}�X�L���O(�v�Z�̈�)
min_mask = (min(reference_point) - 10)*4;
max_mask = (min(reference_point) + 100 - 1)*4;
mask_rfdata = zeros(size(interp_rfdata));
mask_rfdata(min_mask:max_mask,:) = 1;
focused_rfdata_mask = mask_rfdata .* interp_rfdata;
min_reference = (min(reference_point) - 30)*4;
max_reference = (min(reference_point) + 500 - 1)*4;
reference_rfdata = zeros(size(interp_rfdata));
reference_rfdata(min_reference:max_reference,:) = 1;
focused_rfdata_reference = reference_rfdata .* interp_rfdata;
% �g�ʌ`�󐄒�
delay_profile = zeros(1,length(target_element));
for jj = 1:length(target_element)-1
    [acor,lag] = xcorr(focused_rfdata_mask(:,target_element(jj+1)),focused_rfdata_reference(:,target_element(jj)),'coeff');
    [~,I] = max(abs(acor));
    displacement(1,target_element(jj)+1) = lag(I);
    delay_profile(1,jj+1) = sum(displacement(1,target_element(1):target_element(jj)+1));
end
poly_delay_profile_fitted = polyfit(t_pos(1,target_element),(delay_profile*kgrid.dt/4).^2,2); 
displacement_ref(1,target_element(2:end)) = diff(reference_point(target_element));
estimated_velocity = 1/sqrt(-poly_delay_profile_fitted(1))/4.107;

%% �ۑ���