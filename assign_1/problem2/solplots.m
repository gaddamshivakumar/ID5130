clear all
close all
clc

%% 
array = importdata("output_n25.csv");
xpoints = array(:,1);
ana = array(:,2);
serial = array(:,3);
parallel = array(:,4);


figure; hold on;
set(gca,'FontSize',12)
box on
% xlim([xstart,xend])
ylim([-6 8])
pbaspect([1 1 1]);
xlabel('x')
ylabel('df/dx')
legend('Location','northeast')
set(gcf,'units','pixels','position',[100 100 450 450]);
title("problem-2(a): analytical vs. LU decomp.")
plot(xpoints,ana,'.k-','MarkerSize',20,'LineWidth',2,'DisplayName','analytical')
plot(xpoints,serial,'.-','MarkerSize',10,'LineWidth',1,'DisplayName','serial-LU')

%% 

tarray = importdata("time.csv");
p2time = tarray(1,2:end-1);
p4time = tarray(2,2:end-1);
p8time = tarray(3,2:end-1);
threads = [2,4,8];

figure;
boxchart(tarray(:,2:end-1)'*10)
hold on
plot(mean(tarray(:,2:end-1),2)*10,'-o','LineWidth',1.5)
hold off
legend(["time data","time Mean"])
set(gca,'FontSize',12)
box on
% xlim([xstart,xend])
% ylim([-6 8])
pbaspect([1 1 1]);
xlabel('number of threads')
ylabel('time (microseconds)')
legend('Location','northeast')
set(gcf,'units','pixels','position',[1700 100 450 450]);
title("problem-2(b): time.")
xticklabels(threads)










