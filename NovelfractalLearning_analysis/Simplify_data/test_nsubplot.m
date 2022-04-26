tic
figure
nsubplot(5, 5, 1, 1);
plot(1:10, sin(1:10))
toc 

tic
figure
nsubplot(10, 10, 1:2, 1:2);
plot(1:10, sin(1:10))
toc

tic
figure
nsubplot(169, 169, 1:20, 1:20);
plot(1:10, sin(1:10))
toc

running_time = [];
for i = 1: 169
    close all;
    tic
    figure
    nsubplot(169, 169, 1:i, 1:i);
    plot(1:10, sin(1:10))
%     if i==150
%         pause(3)
%     end
    t = toc;
    running_time(end+1) = t;
end



