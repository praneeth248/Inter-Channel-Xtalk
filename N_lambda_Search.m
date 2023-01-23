clc;
clear all;
close all;

wavelength = input('Enter the first wavelength of the band : '); % Example: 1570e-9 (1570 nm)
quality_factor = input('Enter the value of Quality Factor: '); % Example: 7000 
FSR = input('Enter the value of FSR : '); % Example: 20e-9 (For a value of 20 nm)
WG_length = input('Enter the length of Waveguide (cm) : '); % Example: 4 ( If the length of the waveguide is 4 cms) 
WG_loss_dB_cm = input('Enter the Waveguide Loss (dB/cm) : '); % Example: 1.5 ( If the propagation loss is 1.5 dB/cm)
MAOP_WL = input('Enter the value of Maximum Signal Power (per Wavelength) (dBm) : '); % Example: 20 ( For 20 dBm)
MAOP_WG = input('Enter the value of Maximum Signal Power (per Waveguide) (dBm) : '); % Example: 20 ( For 20 dBm)
Coupling_loss_dB = input('Enter the value of Coupling Loss (dB) : '); % Example: 4 ( For a loss value of 4dB)

Vmod = 3;
Vdd = 1.2;
Cmod = 2e-15;
Cref = 50e-15;
c = 3*(10^8);
q0 = 0.04;
g = 0;
h = 0;
m = 0;
d = 0;
modulator = cell(1:1);
detector = cell(1:1);


for Q = quality_factor % If searching for different values of Q, erase the quality_factor here and provide the range. Example Q = 5000:250:10000
    m = 0;
    d = d+1;
for Rb = 1:5 % Provide the range of bit-rates. For each combination of bit-rate and quality factor, an optimal value of N-lambda is calculated
for n = 2:150 % This is the array of N-lambda values. You can change the maximum value from 150 if you want to see whether much higher values are supported 

spac(Rb,n) = FSR/(n+1); % Spacing between the channels
    
for j = 0:n-1
    wl(j+1)=wavelength-j*(spac(Rb,n)); 
end

for j = 1:n
    freq(j) = (c/wl(j));
end

for i = 1:n
     
    for j = 1:n
        FWHM(j) = (freq(j)/Q); % 3-dB Bandwidth 
    end
    
     for j = 1:n
        fd(i,j) = (freq(i)-freq(j)); % Frequency Spacing between the channels
     end
    
    for j = 1:n
        Ei(j) = FWHM(j)/(2*((Rb)*1e9));
    end
    
    for j = 1:n
        Q_intrinsic(Rb,j) = 2*pi*4.2/(100*wl(j)); %Instrinsic Quality Factor
    end
    
    for j = 1:n
        ring_drop_IL(Rb,j) = -10*log10(1-(Q/Q_intrinsic(Rb,j)));
    end
    
    for j = 1:n
        eps(i,j) = fd(i,j);
    end
    
    for j = 1:n
        pp(i,j) = -5*log10((((2*eps(i,j)/FWHM(i))^2)+q0)/(((2*eps(i,j)/FWHM(i))^2)+1));
    end
    
    sum = 0;
    for j = 1:n
        if(i ~= j)
            modulator{Rb,n}.pp_mux(Rb,i) = sum+pp(i,j); %Modulator Crosstalk Penalty
            sum = modulator{Rb,n}.pp_mux(Rb,i);
        end
    end
   
    for j = 1:n
        beta(i,j) = (2*(fd(i,j)))/(FWHM(i));
    end
    
    for j = 1:n
        beta_square(i,j) = beta(i,j)^2;
    end
    
     for j = 1:n
        gamma(i,j) = (1/(1+((beta_square(i,j)))))-((1/(2*pi*Ei(i))).*real((1-exp(-2*pi*Ei(i)*(1-(1j*beta(i,j)))))./((1-beta_square(i,j)-1j*2*beta(i,j)))));
     end
    
    gamma_sqrt = sqrt(gamma);
    
    sum1 = 0;
    for j = 1:n
        p1(Rb,i) = sum1+(gamma(i,j));
        sum1 = p1(Rb,i);
    end
    
    sum_beta_square = 1./(1+beta_square);
    
    sum2 = 0;
    for j = 1:n
        p2(Rb,i) = sum2+(sum_beta_square(i,j));
        sum2 = p2(Rb,i);
    end
    
    pp_demux1(Rb,i) = -5*log10((0.8^2)*(p1(Rb,i))*(p2(Rb,i)));
    
    sum3 = 0;
    for j = 1:n
        if(i ~= j)
            p3(Rb,i) = sum3+(gamma(i,j));
            sum3 = p3(Rb,i);
        end
    end
    
    pp_demux2(Rb,i) = -10*log10(1-(2*p3(Rb,i)));
    
    detector{Rb,n}.pp_demux(Rb,i) = pp_demux1(Rb,i)+pp_demux2(Rb,i)+ring_drop_IL(Rb,i); % Detector Crosstalk Penalty
end


S(Rb) = (-0.0008*((Rb)^2))+(0.8279*(Rb))-30.547; % This equation calculates detector sensitivity based on bit-rate

Per_wavelength_Budget(Rb,n) = MAOP_WL - S(Rb); % This is used to calculate per-wavelength budget

Per_waveguide_Budget(Rb,n) = MAOP_WG - S(Rb); % This is used to calculate per-waveguide budget


Insertion_loss_dB = 0; % If you have any other losses other than propagation and coupling, Please update this before running the code
WG_propogation_loss_dB = WG_length*WG_loss_dB_cm; % Total Propagation Loss in the link

if( imag((detector{Rb,n}.pp_demux(Rb,1:n))) == 0 )
h = n; 
Total_loss(Rb,n) = (Insertion_loss_dB+WG_propogation_loss_dB+Coupling_loss_dB+max(modulator{Rb,n}.pp_mux(Rb,1:n))+max(abs(detector{Rb,n}.pp_demux(Rb,1:n)))); % This calcultes total loss in the link
error_function(Rb,n) = (Per_wavelength_Budget(Rb,n)-Total_loss(Rb,n)); % This error-function (per-wavelength) is based on checking if our total loss is within the per-wavelength budget or not
error_function_1(Rb,n) = (Per_waveguide_Budget(Rb,n)-Total_loss(Rb,n)-10*log10(n)); % This error-function (per-waveguide) is based on checking if our total loss is within the per-waveguide budget or not
end
end

% In order to find the optimal value of N-lambda for a given quality
% factor, bit rate, loss and other paramters, the value of error-function should be a minimum
% positive value. In the upcoming sections of code, we filter out all the
% positive values of error function and we extract the least positive value
% from the filtered out positive values. The N-lambda corresponding to
% minimum positive value of error-function is our Optimal N-lambda.


%% In this section of Code, We filter out positive values of Error function
for k = 2:h
    if((Total_loss(Rb,k)) <= Per_wavelength_Budget(Rb,k))
    if((Total_loss(Rb,k)) <= Per_waveguide_Budget(Rb,k))
        if((error_function(Rb,k))> 0)
            if((error_function_1(Rb,k))> 0)
                loss(Rb,k) = Total_loss(Rb,k);
                ef(Rb,k) = error_function(Rb,k); 
                ef_1(Rb,k) = error_function_1(Rb,k);
                g = k;
            end
        end
    end
    end
end

%% In this section of code, we extract the N-lambda corresponding to minimum positive value of error-function and calculate aggregate bandwidth and energy per bit for the link
if(Total_loss(Rb,2:g) <= Per_wavelength_Budget(Rb,2:g))
if(Total_loss(Rb,2:g) <= Per_waveguide_Budget(Rb,2:g))
if(error_function(Rb,2:g)>=0)
if(error_function_1(Rb,2:g)>=0)
    
m = m+1;
err_fun(m) = min(ef(Rb,2:g));
tot_loss(m) = max(loss(Rb,2:g));
N_lambda(m) = g;
bit_rate(m) = Rb;
aggregate_bandwidth(m) = bit_rate(m)*N_lambda(m);

slope = (1.4*(1e-23))*((Vmod)/(2*(Vdd^2)))*(Cmod/Cref);
constant = (8.4*(1e-14))+((0.25)*((Cmod*(Vmod^2))-(Cref*((2*Vdd)^2))));
Dynamic_Energy(m) = (slope*(Rb*1e9))+constant;


slope = (1.4*(1e-23))*((Vmod)/(2*(Vdd^2)))*((1e-15)/Cref);
constant = (8.4*(1e-14))+((0.25)*((1e-15*(Vmod^2))-(Cref*((2*Vdd)^2))));
Heating_power(m) = (slope*(Rb*1e9))+constant;

Laser_power(m) = 100*((10^(MAOP_WG/10))*1e-3)/(aggregate_bandwidth(m)*(1e9)*10); % Laser power

Receiver_power(m) = 1e-13;

Eb(m) = (Laser_power(m))+(Heating_power(m))+(Dynamic_Energy(m))+(Receiver_power(m)); % Total energy per bit for the link in (pJ/bit)


Nlambda(d,Rb) = N_lambda(m);  % NLambda correspoding to each quality factor (row-wise) and bit-rate (column-wise)
aggban(d,Rb) = aggregate_bandwidth(m);  % Aggregate Bandwidth correspoding to each quality factor (row-wise) and bit-rate (column-wise)
Enpb(d,Rb) = Eb(m);  % Energy Per Bit correspoding to each quality factor (row-wise) and bit-rate (column-wise)

end
end
end
end
end
end
