%% Função de equações e derivadas

function dx= eqsp(Kn, Cn, q, tvec)

q1= q(1); qp1= q(4);
q2= q(2); qp2= q(5);
q3= q(3); qp3= q(6);

qpp1= -Cn(1,1)*qp1-Kn(1,1)*q1;          % derivada segunda de q
qpp2= -Cn(2,2)*qp2-Kn(2,2)*q2;
qpp3= -Cn(3,3)*qp3-Kn(3,3)*q3;

dx= [qp1;qp2;qp3;qpp1;qpp2;qpp3];

end