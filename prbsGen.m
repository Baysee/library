function S_OS=prbsGen(L,OSR)


%% Digital modulation parameters

% P = 7;                                                  % PRBS order
%  L = 2^P-1;                                              % PRBS length

digimod_format = 'QAM';                                 % Digital modulation format
M = 2;                                                 % Number of symbols (modulation order)

K = log2(M);                                            % Number of bits per symbol
N = L*K;                                                % Number of bits
B = idinput([N,2],'prbs',[1,1],[0,1]).'; B = B(1,:);    % Bit sequence
D = bi2de(reshape(B,L,K)).';                            % Symbol indices


%% Symbol sequence

switch digimod_format
    case 'SP-PAM'           % Signle-polarity PAM
        S = D;
    case 'DP-PAM'           % Dual-polarity PAM
        S = pammod(D,M);
    case 'PSK'              % PSK
        S = pskmod(D,M);
    case 'DPSK'             % Differential PSK
        S = dpskmod(D,M);
    case 'QAM'              % QAM
        S = qammod(D,M);
end

S_I = real(S);              % In-phase channel
S_Q = imag(S);              % Quadrature channel

S_M = abs(S);               % Magnitude
S_P = angle(S);             % Phase


%% Oversampling

% OSR = 1;                                   % Oversampling rate

S_OS = S(ones(1,OSR),:); S_OS = S_OS(:).';  % Oversampled symbol sequence

S_OS_I = real(S_OS);                        % Oversampled In-phase channel
S_OS_Q = imag(S_OS);                        % Oversampled Quadrature channel

S_OS_M = abs(S_OS);                         % Oversampled Magnitude
S_OS_P = angle(S_OS);                       % Oversampled Phase


%% Results
% 
% figure
% subplot(2,1,1)
% bar(S_OS_I)
% subplot(2,1,2)
% bar(S_OS_Q)
% 
% figure
% subplot(2,1,1)
% bar(S_OS_M)
% subplot(2,1,2)
% bar(S_OS_P/pi)
% 
% figure
% plot(S_I,S_Q,'o')
% axis square


end