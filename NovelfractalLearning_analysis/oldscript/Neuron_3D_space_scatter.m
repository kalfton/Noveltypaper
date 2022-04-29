
% map the electrode ID on the array
StatisticalThreshold = 0.01;
NovelExcited=find([Neuronlist_good(:).pred_nov_vs_fam]>0 & [Neuronlist_good(:).P_pred_nov_vs_fam]<=StatisticalThreshold)';
NovelInhibited=find([Neuronlist_good(:).pred_nov_vs_fam]<0 & [Neuronlist_good(:).P_pred_nov_vs_fam]<=StatisticalThreshold)';
NotNoveltySelective=find([[Neuronlist_good(:).P_pred_nov_vs_fam]>=StatisticalThreshold])';
NoveltySelective=find([[Neuronlist_good(:).P_pred_nov_vs_fam]<StatisticalThreshold])';

% electrodeML_AP is a 128*3 array where first column is AP, second is ML, and
% 3rd is electrode ID (For validation)
electrodeML = [1:6,repmat(0:7,1,14),1:6, nan, nan, nan, nan]';%M- L+
temp = ones(8,1)*(14:-1:1);
electrodeAP = [15*ones(1,6), temp(:)', 0*ones(1,6), nan, nan, nan, nan]';%P- A+
electrodeML_AP = [electrodeML,electrodeAP,(1:128)'];
clear electrodeML electrodeAP temp;


electrodeID = [Neuronlist_good(:).electrodeID];
% electrodedepth_origin = [Neuronlist_good(:).electrodeDepth]';
% electrodeML_origin = [electrodeML_AP(electrodeID,1)];
% electrodeAP_origin = [electrodeML_AP(electrodeID,2)];
find(cellfun(@isempty,{Neuronlist_good(:).electrodeLocation}))
electrode3D = vertcat(Neuronlist_good(:).electrodeLocation);
electrodeML = electrode3D(:,1);
electrodeAP = electrode3D(:,2);
electrodedepth = electrode3D(:,3); %+0.2*(rand(size(electrodeID))'-0.5);

%CCA analysis?
electrode3D_origin = [electrodedepth_origin,electrodeML_origin,electrodeAP_origin];
tempind = all(~isnan(electrode3D_origin),2);
[A,B,r,U,V] = canoncorr(electrode3D(tempind,:), electrode3D_origin(tempind,:))
figure;
scatter3(U(:,1),U(:,2), U(:,3),5,'filled','b');
xlim([-4,4]);
ylim([-4,4]);
zlim([-4,4]);
figure;
scatter3(V(:,1),V(:,2), V(:,3),5,'filled','b');
xlim([-4,4]);
ylim([-4,4]);
zlim([-4,4]);

figure;hold on;
scatter3(electrodeML(NotNoveltySelective),electrodeAP(NotNoveltySelective),-electrodedepth(NotNoveltySelective),5,[0.7,0.7,0.7],'filled');
scatter3(electrodeML(NovelExcited),electrodeAP(NovelExcited),-electrodedepth(NovelExcited),20,'filled','r'); 
scatter3(electrodeML(NovelInhibited),electrodeAP(NovelInhibited),-electrodedepth(NovelInhibited),20,'filled','b');

legend('Other', 'Nov Excited', 'Nov Inhibited');

xlabel('Medial-Lateral');
ylabel('Anterior- Posterior');
zlabel('Depth');

xlim([-0.5,15.5]);
ylim([-0.5,15.5]);
zlim([-50,-0.5]);

saveas(gcf, fullfile(plotpath,['Plot_3D_all_session' '.fig']));
