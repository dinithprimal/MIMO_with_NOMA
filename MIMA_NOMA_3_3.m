clc; clear variables; close all;

%Distances
d1 = 500; d2 = 200;

%Power allocation coefficients
a1 = 0.8; a2 = 0.2;

N = 5*10^5;
eta = 4;%Path loss exponent

%Rayleigh fading channels
h11 = sqrt(d1^-eta)*(randn(N,1) + 1i*randn(N,1))/sqrt(2);
h12 = sqrt(d1^-eta)*(randn(N,1) + 1i*randn(N,1))/sqrt(2);
h13 = sqrt(d1^-eta)*(randn(N,1) + 1i*randn(N,1))/sqrt(2);
h21 = sqrt(d2^-eta)*(randn(N,1) + 1i*randn(N,1))/sqrt(2);
h22 = sqrt(d2^-eta)*(randn(N,1) + 1i*randn(N,1))/sqrt(2);
h23 = sqrt(d2^-eta)*(randn(N,1) + 1i*randn(N,1))/sqrt(2);

% 2*2 MIMO
h1_2 = h11+h12;
h2_2 = h21+h22;

% 3*3 MIMO
h1_3 = h11+h12+h13;
h2_3 = h21+h22+h23;

%Channel gains 2*2 MIMO
g1_2 = (abs(h1_2)).^2;
g2_2 = (abs(h2_2)).^2;

%Channel gains 3*3 MIMO
g1_3 = (abs(h1_3)).^2;
g2_3 = (abs(h2_3)).^2;

%Transmit power
Pt = -20:5:40; %in dBm
pt = (10^-3)*db2pow(Pt); %linear scale

BW = 10^6;  %bandwidth
%Noise power
No = -174 + 10*log10(BW);   %in dBm
no = (10^-3)*db2pow(No);    %in linear scale

%Target rates
R1 = 1; % Far user target rate
R2 = 4; % Near user target rate

p1n_2 = zeros(1,length(pt));
p2n_2 = zeros(1,length(pt));
p1n_3 = zeros(1,length(pt));
p2n_3 = zeros(1,length(pt));
p1o_2 = zeros(1,length(pt));
p2o_2 = zeros(1,length(pt));
p1o_3 = zeros(1,length(pt));
p2o_3 = zeros(1,length(pt));

for u = 1:length(pt)
   
   %Achievable rates for MIMO-NOMA 2*2
   R1n_2 = log2(1 + pt(u)*a1.*g1_2./(pt(u)*a2.*g1_2 + no));
   R12n_2 = log2(1 + pt(u)*a1.*g2_2./(pt(u)*a2.*g2_2 + no));
   R2n_2 = log2(1 + pt(u)*a2.*g2_2/no);
   
   %Achievable rates for MIMO-NOMA 3*3
   R1n_3 = log2(1 + pt(u)*a1.*g1_3./(pt(u)*a2.*g1_3 + no));
   R12n_3 = log2(1 + pt(u)*a1.*g2_3./(pt(u)*a2.*g2_3 + no));
   R2n_3 = log2(1 + pt(u)*a2.*g2_3/no);
   
   %Achievable rates for MIMO-OMA 2*2
   R1o_2 = 0.5*log2(1 + pt(u)*g1_2/no);
   R2o_2 = 0.5*log2(1 + pt(u)*g2_2/no);
   
   %Achievable rates for MIMO-OMA 2*2
   R1o_3 = 0.5*log2(1 + pt(u)*g1_3/no);
   R2o_3 = 0.5*log2(1 + pt(u)*g2_3/no);
    
   %Sum rates
   Rn_2(u) = mean(R1n_2+R2n_2);    %MIMO-NOMA 2*2
   Rn_3(u) = mean(R1n_3+R2n_3);    %MIMO-NOMA 2*2
   Ro_2(u) = mean(R1o_2+R2o_2);    %MIMO-OMA 2*2
   Ro_3(u) = mean(R1o_3+R2o_3);    %MIMO-OMA 2*2
   
   %Invidual user rates
   R1n_2_av(u) = mean(R1n_2);   %MIMO-NOMA 2*2 USER 1 (FAR)
   R2n_2_av(u) = mean(R2n_2);   %MIMO-NOMA 2*2 USER 2 (NEAR)
   
   R1n_3_av(u) = mean(R1n_3);   %MIMO-NOMA 3*3 USER 1 (FAR)
   R2n_3_av(u) = mean(R2n_3);   %MIMO-NOMA 3*3 USER 2 (NEAR)
   
   R1o_2_av(u) = mean(R1o_2);   %MIMO-OMA 2*2 USER 1 (FAR)
   R2o_2_av(u) = mean(R2o_2);   %MIMO-OMA 2*2 USER 2 (NEAR)
   
   R1o_3_av(u) = mean(R1o_3);   %MIMO-OMA 3*3 USER 1 (FAR)
   R2o_3_av(u) = mean(R2o_3);   %MIMO-OMA 3*3 USER 2 (NEAR)
   
   %Outage calculation
   for k = 1:N
       %MIMO-NOMA 2*2 USER 1 (FAR)
       if R1n_2(k) < R1
           p1n_2(u) = p1n_2(u)+1;
       end
       %MIMO-NOMA 2*2 USER 2 (NEAR)
       if (R12n_2(k)<R1)||((R12n_2(k)>R1)&&(R2n_2(k) < R2))
           p2n_2(u) = p2n_2(u)+1;
       end
       
       %MIMO-NOMA 3*3 USER 1 (FAR)
       if R1n_3(k) < R1
           p1n_3(u) = p1n_3(u)+1;
       end
       %MIMO-NOMA 3*3 USER 2 (NEAR)
       if (R12n_3(k)<R1)||((R12n_3(k)>R1)&&(R2n_3(k) < R2))
           p2n_3(u) = p2n_3(u)+1;
       end
       
       %MIMO-OMA 2*2 USER 1 (FAR)
       if R1o_2(k) < R1
           p1o_2(u) = p1o_2(u)+1;
       end
       %MIMO-OMA 2*2 USER 2 (NEAR)
       if R2o_2(k) < R2
           p2o_2(u) = p2o_2(u)+1;
       end
       
       %MIMO-OMA 3*3 USER 1 (FAR)
       if R1o_3(k) < R1
           p1o_3(u) = p1o_3(u)+1;
       end
       %MIMO-OMA 3*3 USER 2 (NEAR)
       if R2o_3(k) < R2
           p2o_3(u) = p2o_3(u)+1;
       end
   end
end
figure;
plot(Pt, Rn_2, '-sb', 'linewidth',1.5); hold on; grid on;
plot(Pt, Rn_3, '-sc', 'linewidth',1.5); hold on; grid on;
plot(Pt, Ro_2, '-*r','linewidth',1.5);
plot(Pt, Ro_3, '-*g','linewidth',1.5);
legend('MIMO-NOMA 2\times2','MIMO-NOMA 3\times3','MIMO-OMA 2\times2','MIMO-OMA 3\times3');
xlabel('Transmit power (dBm)');
ylabel('Achievable sum rates (bps/Hz)');
str = '$$ Bandwidth = 1MHz , \it(NOMA - \alpha_{1} = '+string(a1)+' (Far User = '+string(d1)+'), \alpha_{2} = '+string(a2)+' (Near User = '+string(d2)+')) $$';
text(1.1,0.5,str,'Interpreter','latex')
title('Sum rate comparison');

figure;
semilogy(Pt,p1n_2/N,'-*b','linewidth',1.5); hold on; grid on;
semilogy(Pt,p2n_2/N,'-ob','linewidth',1.5);
semilogy(Pt,p1n_3/N,'-*g','linewidth',1.5); hold on; grid on;
semilogy(Pt,p2n_3/N,'-og','linewidth',1.5);
semilogy(Pt,p1o_2/N,'-*r','linewidth',1.5); hold on; grid on;
semilogy(Pt,p2o_2/N,'-or','linewidth',1.5);
semilogy(Pt,p1o_3/N,'-*c','linewidth',1.5); hold on; grid on;
semilogy(Pt,p2o_3/N,'-oc','linewidth',1.5);
legend('MIMO-NOMA 2\times2 far','MIMO-NOMA 2\times2 near','MIMO-NOMA 3\times3 far','MIMO-NOMA 3\times3 near','MIMO-OMA 2\times2 far', 'MIMO-OMA 2\times2 near','MIMO-OMA 3\times3 far', 'MIMO-OMA 3\times3 near')
xlabel('Transmit power (dBm)');
ylabel('Outage probability');
str = '$$ Bandwidth = 1MHz , \it(NOMA - \alpha_{1} = '+string(a1)+' (Far User = '+string(d1)+'), \alpha_{2} = '+string(a2)+' (Near User = '+string(d2)+')) $$';
text(1.1,0.5,str,'Interpreter','latex')
title('Outage comparison');

figure;
plot(Pt,R1n_2_av,'-*b','linewidth',1.5); hold on; grid on;
plot(Pt,R2n_2_av,'-ob','linewidth',1.5); 
plot(Pt,R1n_3_av,'-*g','linewidth',1.5); hold on; grid on;
plot(Pt,R2n_3_av,'-og','linewidth',1.5); 
plot(Pt,R1o_2_av,'-*r','linewidth',1.5); 
plot(Pt,R2o_2_av,'-or','linewidth',1.5); 
plot(Pt,R1o_3_av,'-*c','linewidth',1.5); 
plot(Pt,R2o_3_av,'-oc','linewidth',1.5); 
legend('MIMO-NOMA 2\times2 far','MIMO-NOMA 2\times2 near','MIMO-NOMA 3\times3 far','MIMO-NOMA 3\times3 near','MIMO-OMA 2\times2 far', 'MIMO-OMA 2\times2 near','MIMO-OMA 3\times3 far', 'MIMO-OMA 3\times3 near')
xlabel('Transmit power (dBm)');
ylabel('Achievable rates (bps/Hz)');
str = '$$ Bandwidth = 1MHz , \it(NOMA - \alpha_{1} = '+string(a1)+' (Far User = '+string(d1)+'), \alpha_{2} = '+string(a2)+' (Near User = '+string(d2)+')) $$';
text(1.1,0.5,str,'Interpreter','latex')
title('Individual user rates')