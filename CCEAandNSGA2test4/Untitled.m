dat = load('outdata.txt');
y1 = dat(:,1)';
y2 = dat(:,2)';
plot(y1,y2,'.');
xlabel('f1'),ylabel('f2');
title("内部演化次数为10，外部演化次数为10");