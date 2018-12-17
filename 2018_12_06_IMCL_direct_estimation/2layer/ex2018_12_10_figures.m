%%%%%%%%%%%%%%
% Œ‹‰Ê‚Ìfigure‚ğì¬‚·‚éD2018/12/10
%%%%%%%%%%%%%%
load("H:\result\2018_12_06_IMCL_direct_estimation\2layer\2018_12_07_multi_layer_variable.mat")
%%%%%%%%%%%%%%%%%%%%
% ‘ÎÛFƒƒCƒ„ƒ^[ƒQƒbƒgC‹«ŠEŒú‚³F2mm~19mm
% İ’è‰¹‘¬F1580[m/s]³‰ğ‰¹‘¬
% Å“_…•½ˆÊ’uŒÅ’èFy=0
% ‹«ŠEˆÊ’uFŠù’m
% IMCLŠ„‡‚ğ0 %‚ÉŒÅ’è‚·‚éD
% ”g–Ê’x‰„ƒvƒƒtƒ@ƒCƒ‹‚É‚æ‚è‰¹‘¬„’è
% ‘fqŠÔóM”gŠÔ·‚ÍÀM†‚Ì³‹K‰»‘ŠŒİ‘ŠŠÖ‚ğ—p‚¢‚é
%%%%%%%%%%%%%%%%%%%%

%% ‹«ŠE[‚³‚Æ‰¹‘¬„’è¸“x‚ÌŠÖŒW

% ‘S‘Ì}
figure;
plot(focal_depth*1e3,estimated_velocity);
hold on
plot(focal_depth*1e3,correct_velocity,'--');
hold off
xlabel('boundary depth[mm]');
ylabel('estimated velocity');
legend('estimated','correct');
savefilename = sprintf('/estimation_velocity');
savefig([dst_path,savefilename,'.fig']);
exportfig([dst_path,savefilename],'png',[400,300]);
% 4mmˆÈ~‚ÌŠg‘å}
figure;
plot(focal_depth(3:end)*1e3,estimated_velocity(3:end));
hold on
plot(focal_depth(3:end)*1e3,correct_velocity(3:end),'--');
hold off
xlabel('boundary depth[mm]');
ylabel('estimated velocity');
legend('estimated','correct')
savefilename = sprintf('/estimation_velocity_detail');
savefig([dst_path,savefilename,'.fig']);
exportfig([dst_path,savefilename],'png',[400,300]);
% [“x‚²‚Æ‚Ì„’èŒë·‚ğ–_ƒOƒ‰ƒt‚Å¦‚·
figure;
bar(focal_depth*1e3,estimated_velocity - correct_velocity);
xlabel('boundary depth[mm]');
ylabel('error of estimation[m/s] ');
savefilename = sprintf('/error_of_estimation');
savefig([dst_path,savefilename,'.fig']);
exportfig([dst_path,savefilename],'png',[400,300]);
% 4mmˆÈ~‚ÌŠg‘å}
figure;
bar(focal_depth(3:end)*1e3,estimated_velocity(3:end) - correct_velocity(3:end));
xlabel('boundary depth[mm]');
ylabel('error of estimation[m/s] ');
savefilename = sprintf('/error_of_estimation_detail');
savefig([dst_path,savefilename,'.fig']);
exportfig([dst_path,savefilename],'png',[400,300]);