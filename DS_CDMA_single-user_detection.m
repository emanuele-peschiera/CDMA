% This script simulates a baseband synchronous CDMA system with spreading 
% factor Q and number of users K. The spreading codes are randomly
% generated. A simple single-user matched filter (i.e., correlator
% receiver) is used to reconstruct the data. The simulated bit error rates
% are compared to the ideal BPSK and to the model in [1], p. 60.

% [1] L. Hanzo et al., "Single and multi-carrier DS-CDMA: multi-user 
% detection, space-time spreading, synchronisation, networking, and 
% standards," Wiley-IEEE Press, 2004.

clc
clear
close all

%% Parameters initializationQ       = 80;        	% spreading factor
K       = 4;         	% number of users

EbN0    = 0:10;      	% energy per bit to noise power spectral density ratio in dB

% Adapt the iterations to the noise level
Nmax    = 1e8;                                          % maximum number of iterations
Nerr    = 300;                                          % minimum number of errors
Nit     = Nerr./berawgn(EbN0,'psk',2,'nondiff');        % number of iterations

%% Spreading codes generation
Codes = zeros(Q,K);   	% spreading matrix
for k = 1:K
    Codes(:,k) = randi([0 1],1,Q)*2-1;
end

%% System simulation
BER = zeros(K,length(EbN0));

for no = 1:length(EbN0)
    
    fprintf(['Simulating\tEb/N0\t=\t' num2str(EbN0(no)) '\tdB ...\n'])
    err = zeros(K,1);
    
    it = 0;
    
    while it<Nmax && it<Nit(no)

        % Generate the users bits and BPSK modulate
        d = randi([0 1],K,1)*2-1;

        % Apply spreading, sum and normalize the power
        s = Codes*d;

        % Introduce AWGN
        r = awgn(s,10*log10(2)+EbN0(no)-10*log10(Q)+10*log10(K),'measured');

        % Apply despreading and integrate
        y = Codes'*r;
        
        % Demodulate
        x = (sign(y)+1)*0.5;

        % Find the errors
        err = err+double((d+1)*0.5~=x);

        it = it+1;
        
    end
    
    % Evaluate the BER for the current Eb/N0
    BER(:,no) = err./it;
    
end

%% Results
EbN0_ov = EbN0(1):1e-2:EbN0(end);   % oversampled Eb/N0
% Compute analytical BER of BPSK
BER_bpsk = qfunc( (1 ./ (2.*(10.^(EbN0_ov./10))) ).^-0.5 );

% Compute the analytical BER of synchronous CDMA [1]
BER_cdma = zeros(K,length(EbN0_ov));
MAI = zeros(1,K);   % multiple access interference term
for k1 = 1:K
    for k2 = 1:K
        if k1 ~= k2
            MAI(k1) = MAI(k1) + (sum(Codes(:,k1).*Codes(:,k2))/Q)^2;
        end
    end
end
for k = 1:K
    BER_cdma(k,:) = qfunc( ( (1 ./ (2.*(10.^(EbN0_ov./10)))) + MAI(k) ).^-0.5 );
end

% Visualize BER curves for the first user
figure
semilogy(EbN0_ov,BER_bpsk,'k','LineWidth',1.5)
hold on
semilogy(EbN0_ov,BER_cdma(1,:),'-.k','LineWidth',1.5)
semilogy(EbN0,BER(1,:),'--^k','MarkerFaceColor','w','MarkerSize',7,'LineWidth',1.5)
grid on
box on
xlabel('$E_{\mathrm{b}}/N_0$ (dB)','interpreter','latex')
ylabel('BER','interpreter','latex')
xlim([EbN0(1) EbN0(end)])
ylim([1e-6 2e-1])
set(gca,'TickLabelInterpreter','latex')
legend({'BPSK (analytical)','DS-CDMA (analytical)','DS-CDMA (simulation)'},...
    'interpreter','latex','location','southwest')
title(['$Q = ' num2str(Q) ', K = ' num2str(K) '$'],'interpreter','latex')
set(gca,'FontSize',15)
