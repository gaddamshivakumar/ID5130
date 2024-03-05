clear all
close all
clc
run("output.m")

%%
[Lia,Locb] = ismember(0.5,data.ypoints)

%%
figure; hold on;
set(gca,'FontSize',12)
box on
% xlim([xstart,xend])
ylim([0 0.9])
pbaspect([1 1 1]);
xlabel('x (y=0.5)')
ylabel('\phi')
legend('Location','northeast')
set(gcf,'units','pixels','position',[100 100 450 450]);
title("3(a): analytical vs. Gauss-Seidel serial")
plot(data.xpoints,data.phi_ana((data.npoints+1-16),:),'.k-','MarkerSize',20,'LineWidth',2,'DisplayName','analytical')
plot(data.xpoints,data.phi_gs_serial((data.npoints+1-16),:),'.-','MarkerSize',10,'LineWidth',1,'DisplayName','serial-GS')

%%
figure; hold on;
set(gca,'FontSize',12)
box on
% xlim([xstart,xend])
ylim([0 0.9])
pbaspect([1 1 1]);
xlabel('x (y=0.5)')
ylabel('\phi')
legend('Location','northeast')
set(gcf,'units','pixels','position',[100 100 450 450]);
title("3(a): analytical vs. Gauss-Seidel serial")
plot(data.xpoints,data.phi_gs_serial((data.npoints+1-16),:),'.k-','MarkerSize',20,'LineWidth',2,'DisplayName','serial-GS')
plot(data.xpoints,data.phi_gs_diag((data.npoints+1-16),:),'.-','MarkerSize',10,'LineWidth',1,'DisplayName','parallel-diag')

%%
figure; hold on;
set(gca,'FontSize',12)
box on
% xlim([xstart,xend])
ylim([0 0.9])
pbaspect([1 1 1]);
xlabel('x (y=0.5)')
ylabel('\phi')
legend('Location','northeast')
set(gcf,'units','pixels','position',[100 100 450 450]);
title("3(a): analytical vs. Gauss-Seidel serial")
plot(data.xpoints,data.phi_gs_serial((data.npoints+1-16),:),'.k-','MarkerSize',20,'LineWidth',2,'DisplayName','serial-GS')
plot(data.xpoints,data.phi_gs_rb((data.npoints+1-16),:),'.-','MarkerSize',10,'LineWidth',1,'DisplayName','parallel-redblack')

%%
tarray = importdata("comptime.csv");
algorithm = tarray(:,1);
time = mean(tarray(:,2:end-1)')./tarray(:,end)';
d1 = time(1:3:end)*1e9
d01 = time(2:3:end)*1e6
d005 =time(3:3:end)/10
gridsize = [1 2 3]

%%
figure; hold on;
set(gca,'FontSize',12)
box on
xlim([0.5,3.5])
ylim([3 10])
pbaspect([1 1 1]);
xlabel('grid size \Delta')
ylabel('time')
legend('Location','northeast')
set(gcf,'units','pixels','position',[100 100 450 450]);
title("3(c): time vs. grid size")
xticks([1 2 3])
xticklabels(["0.1(ns)","0.01(\mu s)","0.005(\times 10s)"])
plot(gridsize,[d1(1),d01(1),d005(1)],'.','MarkerSize',25,'LineWidth',1,'DisplayName','serial-GS')
plot(gridsize,[d1(2),d01(2),d005(2)],'.','MarkerSize',15,'LineWidth',1,'DisplayName','parallel-diag')
plot(gridsize,[d1(3),d01(3),d005(3)],'.','MarkerSize',10,'LineWidth',1,'DisplayName','parallel-redblack')

%%
tparray = [2 2 4 8 16;
1 82.733 97.885 119.612 197.312;
2 81.36 81.465 87.520 128.096]

xthreads = tparray(1,2:end)
time_diag = tparray(2,2:end)/2
time_redblack = tparray(3,2:end)/2

%%
figure; hold on;
set(gca,'FontSize',12)
box on
xlim([0.5,4.5])
% ylim([3 10])
pbaspect([1 1 1]);
xlabel('no. of threads')
ylabel('time (sec)')
legend('Location','northwest')
set(gcf,'units','pixels','position',[100 100 450 450]);
title("3(d): time vs. threads")
xticklabels(["2", "4", "6", "8"])
plot(1:4,time_diag,'.-','MarkerSize',15,'LineWidth',1,'DisplayName','parallel-diag')
plot(1:4,time_redblack,'.-','MarkerSize',15,'LineWidth',1,'DisplayName','parallel-redblack')
















