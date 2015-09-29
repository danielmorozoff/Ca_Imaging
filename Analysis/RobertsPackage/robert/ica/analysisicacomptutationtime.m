%% this is the script used in the keyboard of b1 to identify computation time
% currentx=128;
% currenty=128;
% xrad = 5;
% yrad = 5;
% tic
% FindCells;
% toc


%%

cpuavgtime1 = [1.88 2.22 4.86 11.50 47.35 171.7 382];
cpustdtime1 = [1.11 1.35 2.59 1.62 4.37 27.2 17];
px1 = (2*[5 10 15 20 25 30 35]+1).^2;
cpuavgtime2 = [1.06 1.39 1.77 5.95 23.4 61.0 130 309];
cpustdtime2 = [0.72 0.35 0.58 0.68 6.7 7.6 2.08 9.33];
px2 = (2*[5 10 15 20 25 30 35 40]+1).^2;


cpumultiple1 = px1(2:end)/px1(2);
cpumultiple2 = px2(2:end)/px1(2);
expavgtime1 = cpuavgtime1(1)*cpumultiple1;
expavgtime2 = cpuavgtime2(2)*cpumultiple2;
expstdtime1 = cpustdtime1(1)*cpumultiple1;
expstdtime2 = cpustdtime2(2)*cpumultiple2;



figure;
errorbar(px1,cpuavgtime1,cpustdtime1,'b')
hold on
errorbar(px1(2:end),expavgtime1,expstdtime1,'r');
hold off
xlabel('pixels')
ylabel('cpu time')
set(gca,'XTick',[0 2500 5000 7500],'YTick',[0 200 400])
set(gcf,'Units','centimeters','Position',[5 5 8 6])
box off

figure;
errorbar(px2,cpuavgtime2,cpustdtime2)
hold on
errorbar(px2(2:end),expavgtime2,expstdtime2,'r');
hold off
xlabel('pixels')
ylabel('cpu time')
set(gca,'XTick',[0 2500 5000 7500],'YTick',[0 200 400])
set(gcf,'Units','centimeters','Position',[5 5 8 6])
box off