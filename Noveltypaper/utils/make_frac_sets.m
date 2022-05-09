% Sets of fractals:
function frac_sets = make_frac_sets(Generaltask, fracsplitlogical)
% this function is to make some common fractal sets.

% successful fractal & violation fractal in type 2&3
% in trial type 2&3, pick the same number of normal fractal as
% comparison
frac_sets.successfulfrac = find(~isnan(Generaltask.Fractals(5,:)));
frac_sets.successfulfrac_split = find(~isnan(Generaltask.Fractals(5,:))& fracsplitlogical);
numberb=2;

temp_ind = intersect(frac_sets.successfulfrac_split,find(ismember(Generaltask.Fractals(1,:),[6401])));
violation_frac_2nd_6401 = temp_ind(Generaltask.Fractals(1,temp_ind-1) ~= Generaltask.Fractals(1,temp_ind)-1);
temp_ind = setdiff(temp_ind,violation_frac_2nd_6401);
normal_frac_2nd_contrast_6401 = temp_ind(randperm(numel(temp_ind)));
normal_frac_2nd_contrast_6401 = sort(normal_frac_2nd_contrast_6401(1:min(numberb*numel(violation_frac_2nd_6401),numel(normal_frac_2nd_contrast_6401))));

temp_ind = intersect(frac_sets.successfulfrac_split,find(ismember(Generaltask.Fractals(1,:),[6402])));
violation_frac_3rd_6402 = temp_ind((Generaltask.Fractals(1,temp_ind-1) ~= Generaltask.Fractals(1,temp_ind)-1) & (Generaltask.Fractals(1,temp_ind-2) ~= Generaltask.Fractals(1,temp_ind)-2)); % 3nd violation fractal
violation_frac_3rd_2_6402 = temp_ind(Generaltask.Fractals(1,temp_ind-1) ~= Generaltask.Fractals(1,temp_ind)-1 & Generaltask.Fractals(1,temp_ind-2) == Generaltask.Fractals(1,temp_ind)-2); % 3rd fractal following the 2nd violation fractal
temp_ind = setdiff(temp_ind,violation_frac_3rd_6402);
normal_frac_3rd_contrast_6402 = temp_ind(randperm(numel(temp_ind)));
try
    normal_frac_3rd_contrast_6402 = sort(normal_frac_3rd_contrast_6402(1:numberb*numel(violation_frac_3rd_6402)));
catch %some session happened to have a lot of violated fractals and the normal fractal is not enough to match them, here is the solution
    normal_frac_3rd_contrast_6402 = sort(normal_frac_3rd_contrast_6402(1:floor(length(normal_frac_3rd_contrast_6402))));
end

%%%%%%%%%%%%%%%%%%%%

temp_ind = intersect(frac_sets.successfulfrac_split,find(ismember(Generaltask.Fractals(1,:),[6404])));
violation_frac_2nd_6404 = temp_ind(Generaltask.Fractals(1,temp_ind-1) ~= Generaltask.Fractals(1,temp_ind)-1);
temp_ind = setdiff(temp_ind,violation_frac_2nd_6404);
normal_frac_2nd_contrast_6404 = temp_ind(randperm(numel(temp_ind)));
normal_frac_2nd_contrast_6404 = sort(normal_frac_2nd_contrast_6404(1:numberb*numel(violation_frac_2nd_6404)));

temp_ind = intersect(frac_sets.successfulfrac_split,find(ismember(Generaltask.Fractals(1,:),[6405])));
violation_frac_3rd_6405 = temp_ind((Generaltask.Fractals(1,temp_ind-1) ~= Generaltask.Fractals(1,temp_ind)-1) & (Generaltask.Fractals(1,temp_ind-2) ~= Generaltask.Fractals(1,temp_ind)-2)); % 3nd violation fractal
violation_frac_3rd_2_6405 = temp_ind(Generaltask.Fractals(1,temp_ind-1) ~= Generaltask.Fractals(1,temp_ind)-1 & Generaltask.Fractals(1,temp_ind-2) == Generaltask.Fractals(1,temp_ind)-2); % 3rd fractal following the 2nd violation fractal
temp_ind = setdiff(temp_ind,violation_frac_3rd_6405);
normal_frac_3rd_contrast_6405 = temp_ind(randperm(numel(temp_ind)));
normal_frac_3rd_contrast_6405 = sort(normal_frac_3rd_contrast_6405(1:numberb*numel(violation_frac_3rd_6405)));

violation_frac_2nd = [violation_frac_2nd_6401,violation_frac_2nd_6404];
violation_frac_3rd = [violation_frac_3rd_6402,violation_frac_3rd_6405];
violation_frac_3rd_2 = [violation_frac_3rd_2_6402,violation_frac_3rd_2_6405];
normal_frac_2nd_contrast = [normal_frac_2nd_contrast_6401,normal_frac_2nd_contrast_6404];
normal_frac_3rd_contrast = [normal_frac_3rd_contrast_6402,normal_frac_3rd_contrast_6405];


frac_sets.successfulfrac_noviol = setdiff(frac_sets.successfulfrac_split,[violation_frac_2nd, violation_frac_3rd, violation_frac_3rd_2,normal_frac_3rd_contrast,normal_frac_2nd_contrast],'sorted');

frac_sets.violation_frac_2nd_6401 = violation_frac_2nd_6401;
frac_sets.violation_frac_2nd_6404 = violation_frac_2nd_6404;
frac_sets.violation_frac_3rd_6402 = violation_frac_3rd_6402;
frac_sets.violation_frac_3rd_6405 = violation_frac_3rd_6405;
frac_sets.normal_frac_2nd_contrast_6401 = normal_frac_2nd_contrast_6401;
frac_sets.normal_frac_2nd_contrast_6404 = normal_frac_2nd_contrast_6404;
frac_sets.normal_frac_3rd_contrast_6402 = normal_frac_3rd_contrast_6402;
frac_sets.normal_frac_3rd_contrast_6405 = normal_frac_3rd_contrast_6405;

frac_sets.violation_frac_2nd = violation_frac_2nd;
frac_sets.violation_frac_3rd = violation_frac_3rd;
frac_sets.normal_frac_2nd_contrast = normal_frac_2nd_contrast;
frac_sets.normal_frac_3rd_contrast = normal_frac_3rd_contrast;


end