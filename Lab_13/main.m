DATA_NUM = 1200
DFT_SIZE = 2048                %N
USED_SUB_CARRIER = DATA_NUM/4  %M
UP_FACTOR = 4
DOWN_FACTOR = 4
AWGN_DB = 15
CP_NUM = 12
CP_SIG_LEN = CP_NUM + DFT_SIZE

load('./filter/practice_filter')
data = randi([0 1],DATA_NUM,1)
% ------------------Tx------------------------
% transform data into QAM data
modulated_data = QAM_MAPPER(data)

% Make QAM data modulate on subcarrier
sub_data = subcarrier_mapping(modulated_data,DFT_SIZE,USED_SUB_CARRIER)

% idft
idft_data = idft(sub_data,DFT_SIZE)

% cyclic prefix
cp_data = cyclic_prefix(idft_data,CP_NUM)

% DAC
up_data = DAC(cp_data,UP_FACTOR)

% filter the signal
% TXfilter_data = filter(practice_filter,up_data)

% AWGN Channel
% channel_data = add_awgn_noise(TXfilter_data,AWGN_DB)
% ------------------Rx------------------------

% filter in reciever
% RXfilter_data = filter(practice_filter,channel_data)

% down sampling data
down_data = ADC(up_data,DOWN_FACTOR)


% remove cyclic prefix
decp_data = cp_remove(down_data,CP_SIG_LEN,CP_NUM)

% FFT Block
dft_data = dft(decp_data,DFT_SIZE)

% subcarrier demapping
desub_data = subcarrier_demapping(dft_data,DFT_SIZE,USED_SUB_CARRIER)

% demutiplexor
demux_data = demultiplexing(desub_data,1)

% QAM demodulate
QAM_demod_data = QAM_demod(demux_data)
