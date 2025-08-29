function C_plot(nodes,Dati)
figure()
x = [Dati.domain(1),Dati.domain(1)];
y = [0,Dati.T];
hold on 
plot(x,y,'b')
x = [Dati.domain(1),Dati.domain(2)];
y = [Dati.T,Dati.T];
plot(x,y,'b')
x=[Dati.domain(2),Dati.domain(2)];
y = [0,Dati.T];
plot(x,y,'b')
x = [Dati.domain(1),Dati.domain(1)];
y = [0,0];
plot(x,y,'b')
for i=2:size(nodes,1)
    plot(nodes(1,:),nodes(i,:),'r')
end
end