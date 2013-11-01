function [] = stt_ver6(NSimulations)

close all;
clc;

val=randi(15,1);
for tt=1:NSimulations
    file=strcat('sims_results1/lifetime_',num2str(tt),'.txt');
    y=load(file);
    res(tt)=mean(y);
    if tt==val
        figure
        hist(y,20)
    end
end

s=sprintf('The mean lifetime for all simulations is: %f seconds (%f min)...\n',mean(res),mean(res)/60)
set(gca,'fontweight','b','fontsize',16);
xlabel('Time (s)','fontweight','b','fontsize',16);
ylabel('# Extinct MTs','fontweight','b','fontsize',16);
end