clc;
w=2.2745d0*sqrt(2d0);
sigma=w/(2d0*sqrt(2d0*log(2d0)));

b=1/(2d0*(sigma^2));
rho_1=zeros(2001,2001);
near_field_x=zeros(1,2001);
rho_1p5=zeros(2001,2001);
coherence_dist=zeros(1,2001);
initial_coh_dist=zeros(1,2001);
probability=zeros(1,2001);


for i=1:2001
 for j=1:2001
  x=i*10d0/2000d0;
  y=j*10d0/2000d0;
  rho_1(i,j)=exp(-b*((x-5.005)^2+(y-5.005)^2));
 end
end

for i=1:2001
 near_field_x(i)=-5d0+(i-1)*.005d0;
end

%fwhm1p67=16.48; .4628
%fwhm1p67=189.1; %Attempt to obtain density width of 557nm (525nm)
%fwhm1p67=210.2; %Attempt to obtain density of 512nm (lcoh=570nm),new si
fwhm1p67=215.44; %Attempt to obtain density of 524.1nm (lcoh=583.1229nm), old si
%fwhm1p67=220.5; %Attempt to obtain density of 535.7nm (lcoh=595.5127m), gold


beta=fwhm1p67/(2d0*sqrt(log(2d0)));

for i=1:2001
  for j=1:2001
   rho_1p5(i,j)=rho_1(i,j)*exp(-((i-j)/beta)^2);
  end
end


for i=1:2001
 j=2002-i;
 coherence_dist(i)=rho_1p5(i,j);
 initial_coh_dist(i)=rho_1(i,j);
 probability(i)=rho_1(i,i);
end


plot(near_field_x,coherence_dist,'.-k')
hold on
plot(near_field_x,initial_coh_dist,'.-b')
plot(near_field_x,probability,'.-m')
hold off


f = fit(transpose(near_field_x),transpose(coherence_dist),'gauss1')

clc;
coeffvals = coeffvalues(f);
width=2*sqrt(log(2))*coeffvals(3)


%%
%operator which acts as a grating to the wavefunction distribution
grating_operator=zeros(1,numel(near_field_x));

for i=1:numel(near_field_x)
        if mod(near_field_x(i)*100,20)>= 0
            if mod(near_field_x(i)*100,20)< 10
                grating_operator(i)=1;
            end
        end
end

%%
%send coherent_function_through_grating
clc;
%first normalize rho_1
rho_1=rho_1/trace(rho_1);
for i=1:2001 
 wavefunction_1(i)=sqrt(rho_1(i,i))*grating_operator(i);
end
plot(near_field_x,wavefunction_1)


%%
clc;
wavefunction_1extended(1:1000)=0;
wavefunction_1extended(1001:3001)=wavefunction_1;
wavefunction_1extended(3002:4001)=0;
plot(wavefunction_1extended)
%%

farfield_x=[-2000:2000]*.72; %far field domain to obtain a ff peak_periodicity of 72um

Y=fft(wavefunction_1extended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result_coherent=abs(farfield_pattern)/max(abs(farfield_pattern));
%%
plot(farfield_x,fourier_result_coherent,'.-k')
hold on
%xlim([68 76])
hold off
%%
clc;
% 0th order peak
format long
f = fit(transpose(farfield_x(1960:2040)),transpose(fourier_result_coherent(1960:2040)),'gauss1')
coeffvals = coeffvalues(f);
width=2*sqrt(log(2))*coeffvals(3);
.1*72/width % = 2.573496999326959
%%
clc;
% 1st order peak
format long
f = fit(transpose(farfield_x(2060:2140)),transpose(fourier_result_coherent(2060:2140)),'gauss1')
coeffvals = coeffvalues(f);
width=2*sqrt(log(2))*coeffvals(3);
.1*72/width % = 2.576545361683533

%%
clc;
2.2745/2.573496999326959%ratio of diagonal width to coherence length
2.2745/2.576545361683533
%%
clc;
1/(2*3.14159265)
%%
clc;

72*.1/(384.7587-367.3198)
