clear all;
load CoRDF222.dat;
load CoRDF82000.dat;
x = CoRDF222;
y = CoRDF82000;


plot(x(:,1)/10,x(:,2),'b',y(:,1)/10,y(:,2),'r','LineWidth',1.5);
xlim([0.2,2]);
xlabel('distance (nm)');
ylabel('g(r)');
legend('random distribution','after 80000 MC');
set(gca,'FontSize',12);