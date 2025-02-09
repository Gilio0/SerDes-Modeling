clear;
close all;
clc;

%generation data stream
num_bits =100000 ;
binary_data_stream = randi( [0,1] , 1,num_bits );
binary_data_stream2 = randi( [0,1] , 1,num_bits );
%% parameter
Eb=1;
%snr_db= -5 : 1 : 20 ;
snr_db=10;         %no noise
snr_linear = 10 .^( snr_db / 10 );
N0 = Eb ./snr_linear;
mu=0.0001;
%% pam4 mapping
modulated_data=PAM4_mapping(binary_data_stream);
modulated_data2=PAM4_mapping(binary_data_stream2);
pam4_table = [(-3), (-1), (1) , (3)];
%% filter taps
FFE_taps=21;
DFE_taps=2;

%% channel response
h_channel=[1.2 0.7 -0.5];
channel_length=length(h_channel);
channel_out=filter(h_channel,1,modulated_data);
channel_out2=filter(h_channel,1,modulated_data2);
%%adding AWGN noise
Es=( sum( ( real(modulated_data) .^ 2 ) )  ) / (num_bits/2) ;
E_avg=Es/2;
noise_inphase = randn(length(modulated_data),1) .* sqrt( (N0/2)*E_avg  )   ;
channel_out_noise=channel_out+noise_inphase;

Es=( sum( ( real(modulated_data2) .^ 2 ) )  ) / (num_bits/2) ;
E_avg=Es/2;
noise_inphase = randn(length(modulated_data2),1) .* sqrt( (N0/2)*E_avg  )   ;
channel_out_noise2=channel_out2+noise_inphase;

%% outputs
num_data=length(modulated_data);
FFE_output=zeros(1,num_data);
DFE_output=zeros(1,num_data);
desicicon_in=zeros(1,num_data);
des_out = zeros(1,num_data);
error = zeros(1,num_data);
% ffe_history=zeros(num_data,FFE_taps);
% dfe_history=zeros(num_data,DFE_taps);

%% training seq
seq_length=50000;
Forword_Buffer=zeros(1,FFE_taps);Feedback_Buffer=zeros(1,DFE_taps);
Forward_taps=zeros(FFE_taps,num_data+1);Forward_taps(10,:)=1;
Feedback_taps=zeros(DFE_taps,num_data+1);
for n=1:1:seq_length
    Forword_Buffer = [channel_out_noise(n), Forword_Buffer(1:end-1)];
    FFE_output(n)= Forward_taps(:,n)' * Forword_Buffer' ;
    DFE_output(n)= Feedback_taps(:,n)'* Feedback_Buffer' ;
   
    desicicon_in(n)=FFE_output(n)+DFE_output(n);
    
   
    error(n)= modulated_data(n)- desicicon_in(n);
    
    Forward_taps(:,n+1) = Forward_taps(:,n) + (mu*error(n)) .* Forword_Buffer';
    Feedback_taps(:,n+1) = Feedback_taps(:,n) + (mu*error(n)) .* Feedback_Buffer';
    
    Feedback_Buffer = [modulated_data(n), Feedback_Buffer(1:end-1)];
end
%%exact
%w_mmse=MMSE_exact(channel_out_noise,modulated_data,21);


%% data
ffe_taps=Forward_taps(:,seq_length);
dfe_taps=Feedback_taps(:,seq_length);
des_out2 = zeros(1,num_data);
% ffe_history=zeros(FFE_taps,num_data+1);
% dfe_history=zeros(DFE_taps,num_data+1);
for n=1:1:num_data
    Forword_Buffer = [channel_out_noise2(n), Forword_Buffer(1:end-1)];
    FFE_output(n)= ffe_taps' * Forword_Buffer' ;
    DFE_output(n)= dfe_taps'* Feedback_Buffer' ;
   
    desicicon_in(n)=FFE_output(n)+DFE_output(n);
    
    distance = abs( desicicon_in( n ) - pam4_table );
    [min_value , index ] = min(distance);
    des_out2(n)= pam4_table(index);
    
    error(n)= des_out2(n)- desicicon_in(n);
    
    ffe_taps = ffe_taps + (mu*error(n)) .* Forword_Buffer';
    dfe_taps = dfe_taps + (mu*error(n)) .* Feedback_Buffer';
    
    Feedback_Buffer = [des_out2(n), Feedback_Buffer(1:end-1)];
end
%%
num_error=sum(des_out2 ~=modulated_data2' );
%disp(num_error);
fprintf('number of error = %d',num_error);
%%
figure('name','Equalizers taps' );
x_axis=1:num_data+1;
subplot(2,1,1);
plot(x_axis,Forward_taps);
title('FFE taps');
ylabel('taps weight' );
xlabel('iterations');
grid on ;

subplot(2,1,2);
plot(x_axis,Feedback_taps);
title('DFE taps');
ylabel('taps weight' );
xlabel('iterations');
grid on ;

function pam4_mapping = PAM4_mapping( pam4_mapping )

sympol_data_binary = reshape( pam4_mapping , [ 2 length(pam4_mapping)/2 ]  );
sympol_data_binary = sympol_data_binary' ;
sympol_in_decimal = bin2dec( num2str(sympol_data_binary ) );
pam4_table = [(-3), (-1), (1) , (3)];
transmiited_out = zeros(length(pam4_mapping)/2, 1);
for i=1:size(sympol_in_decimal)
    
transmiited_out( i ) = pam4_table(sympol_in_decimal(i) +1);
end
pam4_mapping=transmiited_out;

end