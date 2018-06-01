clear all; close all;
load('sparam_waveguide_data.mat')
a1abs=abs(a2).^2;
a2abs=abs(a2).^2;
b1abs=abs(b1).^2;
b2abs=abs(b2).^2;

figure()
subplot(2,1,1)
hold on
plot(abs(a1).^2)
plot(abs(a2).^2)
plot(abs(b1).^2,'--')
plot(abs(b2).^2,'--')

S12 = abs(a2./(a1)).^2;
plot(abs(b2./(a1)).^2,'linewidth',2)
hold off

legend('a1','a2','b1','b2')

subplot(2,1,2)
fftS12 = abs(fftshift(fft(S12))).^2;
plot(fftS12)