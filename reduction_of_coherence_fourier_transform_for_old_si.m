
plot(near_field_x,coherence_dist,'.-k')
hold on
plot(near_field_x,initial_coh_dist,'.-b')
plot(near_field_x,probability,'.-m')
hold off
%%

clc;

for i=1:2001
 near_field_x(i)=-5d0+(i-1)*.005d0;
end

w=2.2745d0;
%wdblprime=.5;% produces lcoh=.557;
%wdblprime=.5120; %produces lcoh = .570
wdblprime=.5241; %produces lcoh = 583.1229, old si
%wdblprime=.5357; %produces lcoh = 595.5127, gold
sigmadbl=wdblprime/(2d0*sqrt(2d0*log(2d0)));
delx=.1;

wprime=sqrt((w*w)-(wdblprime*wdblprime));
sigmaprime=wprime/(2d0*sqrt(2d0*log(2d0)));


a0=exp(-(0).^2/(2*sigmaprime*sigmaprime));
a1=exp(-(delx).^2/(2*sigmaprime*sigmaprime));
a2=exp(-(2*delx).^2/(2*sigmaprime*sigmaprime));
a3=exp(-(3*delx).^2/(2*sigmaprime*sigmaprime));
a4=exp(-(4*delx).^2/(2*sigmaprime*sigmaprime));
a5=exp(-(5*delx).^2/(2*sigmaprime*sigmaprime));
a6=exp(-(6*delx).^2/(2*sigmaprime*sigmaprime));
a7=exp(-(7*delx).^2/(2*sigmaprime*sigmaprime));
a8=exp(-(8*delx).^2/(2*sigmaprime*sigmaprime));
a9=exp(-(9*delx).^2/(2*sigmaprime*sigmaprime));
a10=exp(-(10*delx).^2/(2*sigmaprime*sigmaprime));
a11=exp(-(11*delx).^2/(2*sigmaprime*sigmaprime));
a12=exp(-(12*delx).^2/(2*sigmaprime*sigmaprime));
a13=exp(-(13*delx).^2/(2*sigmaprime*sigmaprime));
a14=exp(-(14*delx).^2/(2*sigmaprime*sigmaprime));
a15=exp(-(15*delx).^2/(2*sigmaprime*sigmaprime));
a16=exp(-(16*delx).^2/(2*sigmaprime*sigmaprime));
a17=exp(-(17*delx).^2/(2*sigmaprime*sigmaprime));
a18=exp(-(18*delx).^2/(2*sigmaprime*sigmaprime));
a19=exp(-(19*delx).^2/(2*sigmaprime*sigmaprime));
a20=exp(-(20*delx).^2/(2*sigmaprime*sigmaprime));
a21=exp(-(21*delx).^2/(2*sigmaprime*sigmaprime));
a22=exp(-(22*delx).^2/(2*sigmaprime*sigmaprime));
a23=exp(-(23*delx).^2/(2*sigmaprime*sigmaprime));
a24=exp(-(24*delx).^2/(2*sigmaprime*sigmaprime));
a25=exp(-(25*delx).^2/(2*sigmaprime*sigmaprime));
a26=exp(-(26*delx).^2/(2*sigmaprime*sigmaprime));
a27=exp(-(27*delx).^2/(2*sigmaprime*sigmaprime));
a28=exp(-(28*delx).^2/(2*sigmaprime*sigmaprime));
a29=exp(-(29*delx).^2/(2*sigmaprime*sigmaprime));
a30=exp(-(30*delx).^2/(2*sigmaprime*sigmaprime));

rho_1p5_prob0=a0*exp(-(near_field_x).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob1a=a1*exp(-(near_field_x-delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob1b=a1*exp(-(near_field_x+delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob2a=a2*exp(-(near_field_x-2*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob2b=a2*exp(-(near_field_x+2*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob3a=a3*exp(-(near_field_x-3*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob3b=a3*exp(-(near_field_x+3*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob4a=a4*exp(-(near_field_x-4*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob4b=a4*exp(-(near_field_x+4*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob5a=a5*exp(-(near_field_x-5*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob5b=a5*exp(-(near_field_x+5*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob6a=a6*exp(-(near_field_x-6*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob6b=a6*exp(-(near_field_x+6*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob7a=a7*exp(-(near_field_x-7*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob7b=a7*exp(-(near_field_x+7*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob8a=a8*exp(-(near_field_x-8*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob8b=a8*exp(-(near_field_x+8*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob9a=a9*exp(-(near_field_x-9*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob9b=a9*exp(-(near_field_x+9*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob10a=a10*exp(-(near_field_x-10*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob10b=a10*exp(-(near_field_x+10*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob11a=a11*exp(-(near_field_x-11*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob11b=a11*exp(-(near_field_x+11*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob12a=a12*exp(-(near_field_x-12*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob12b=a12*exp(-(near_field_x+12*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob13a=a13*exp(-(near_field_x-13*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob13b=a13*exp(-(near_field_x+13*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob14a=a14*exp(-(near_field_x-14*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob14b=a14*exp(-(near_field_x+14*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob15a=a15*exp(-(near_field_x-15*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob15b=a15*exp(-(near_field_x+15*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob16a=a16*exp(-(near_field_x-16*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob16b=a16*exp(-(near_field_x+16*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob17a=a17*exp(-(near_field_x-17*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob17b=a17*exp(-(near_field_x+17*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob18a=a18*exp(-(near_field_x-18*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob18b=a18*exp(-(near_field_x+18*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob19a=a19*exp(-(near_field_x-19*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob19b=a19*exp(-(near_field_x+19*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob20a=a20*exp(-(near_field_x-20*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob20b=a20*exp(-(near_field_x+20*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob21a=a21*exp(-(near_field_x-21*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob21b=a21*exp(-(near_field_x+21*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob22a=a22*exp(-(near_field_x-22*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob22b=a22*exp(-(near_field_x+22*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob23a=a23*exp(-(near_field_x-23*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob23b=a23*exp(-(near_field_x+23*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob24a=a24*exp(-(near_field_x-24*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob24b=a24*exp(-(near_field_x+24*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob25a=a25*exp(-(near_field_x-25*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob25b=a25*exp(-(near_field_x+25*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob26a=a26*exp(-(near_field_x-26*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob26b=a26*exp(-(near_field_x+26*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob27a=a27*exp(-(near_field_x-27*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob27b=a27*exp(-(near_field_x+27*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob28a=a28*exp(-(near_field_x-28*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob28b=a28*exp(-(near_field_x+28*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob29a=a29*exp(-(near_field_x-29*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob29b=a29*exp(-(near_field_x+29*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob30a=a30*exp(-(near_field_x-30*delx).^2/(2*sigmadbl*sigmadbl));
rho_1p5_prob30b=a30*exp(-(near_field_x+30*delx).^2/(2*sigmadbl*sigmadbl));

rho_1p5_total=rho_1p5_prob0+rho_1p5_prob1a+rho_1p5_prob1b+rho_1p5_prob2a+rho_1p5_prob2b...
                           +rho_1p5_prob3a+rho_1p5_prob3b+rho_1p5_prob4a+rho_1p5_prob4b...
                           +rho_1p5_prob5a+rho_1p5_prob5b+rho_1p5_prob6a+rho_1p5_prob6b...
                           +rho_1p5_prob7a+rho_1p5_prob7b+rho_1p5_prob8a+rho_1p5_prob8b...
                           +rho_1p5_prob9a+rho_1p5_prob9b+rho_1p5_prob10a+rho_1p5_prob10b...
                           +rho_1p5_prob11a+rho_1p5_prob11b+rho_1p5_prob12a+rho_1p5_prob12b...
                           +rho_1p5_prob13a+rho_1p5_prob13b+rho_1p5_prob14a+rho_1p5_prob14b...
                           +rho_1p5_prob15a+rho_1p5_prob15b+rho_1p5_prob16a+rho_1p5_prob16b...
                           +rho_1p5_prob17a+rho_1p5_prob17b+rho_1p5_prob18a+rho_1p5_prob18b...
                           +rho_1p5_prob19a+rho_1p5_prob19b+rho_1p5_prob20a+rho_1p5_prob20b...
                           +rho_1p5_prob21a+rho_1p5_prob21b+rho_1p5_prob22a+rho_1p5_prob22b...
                           +rho_1p5_prob23a+rho_1p5_prob23b+rho_1p5_prob24a+rho_1p5_prob24b...
                           +rho_1p5_prob25a+rho_1p5_prob25b+rho_1p5_prob26a+rho_1p5_prob26b...
                           +rho_1p5_prob27a+rho_1p5_prob27b+rho_1p5_prob28a+rho_1p5_prob28b...
                           +rho_1p5_prob29a+rho_1p5_prob29b+rho_1p5_prob30a+rho_1p5_prob30b;
                       
plot(near_field_x,probability,'m')
hold on
plot(near_field_x,rho_1p5_prob2a,'g')
plot(near_field_x,rho_1p5_prob15a,'g')
%plot(near_field_x,rho_1p5_prob1,'b')
plot(near_field_x,rho_1p5_total./max(rho_1p5_total),'k')
hold off

%%

for i=1:2001 
 wavefun_1p5_0(i)=sqrt(rho_1p5_prob0(i))*grating_operator(i);
 wavefun_1p5_1a(i)=sqrt(rho_1p5_prob1a(i))*grating_operator(i);
 wavefun_1p5_1b(i)=sqrt(rho_1p5_prob1b(i))*grating_operator(i);
 wavefun_1p5_2a(i)=sqrt(rho_1p5_prob2a(i))*grating_operator(i);
 wavefun_1p5_2b(i)=sqrt(rho_1p5_prob2b(i))*grating_operator(i);
 wavefun_1p5_3a(i)=sqrt(rho_1p5_prob3a(i))*grating_operator(i);
 wavefun_1p5_3b(i)=sqrt(rho_1p5_prob3b(i))*grating_operator(i);
 wavefun_1p5_4a(i)=sqrt(rho_1p5_prob4a(i))*grating_operator(i);
 wavefun_1p5_4b(i)=sqrt(rho_1p5_prob4b(i))*grating_operator(i);
 wavefun_1p5_5a(i)=sqrt(rho_1p5_prob5a(i))*grating_operator(i);
 wavefun_1p5_5b(i)=sqrt(rho_1p5_prob5b(i))*grating_operator(i);
 wavefun_1p5_6a(i)=sqrt(rho_1p5_prob6a(i))*grating_operator(i);
 wavefun_1p5_6b(i)=sqrt(rho_1p5_prob6b(i))*grating_operator(i);
 wavefun_1p5_7a(i)=sqrt(rho_1p5_prob7a(i))*grating_operator(i);
 wavefun_1p5_7b(i)=sqrt(rho_1p5_prob7b(i))*grating_operator(i);
 wavefun_1p5_8a(i)=sqrt(rho_1p5_prob8a(i))*grating_operator(i);
 wavefun_1p5_8b(i)=sqrt(rho_1p5_prob8b(i))*grating_operator(i);
 wavefun_1p5_9a(i)=sqrt(rho_1p5_prob9a(i))*grating_operator(i);
 wavefun_1p5_9b(i)=sqrt(rho_1p5_prob9b(i))*grating_operator(i);
 wavefun_1p5_10a(i)=sqrt(rho_1p5_prob10a(i))*grating_operator(i);
 wavefun_1p5_10b(i)=sqrt(rho_1p5_prob10b(i))*grating_operator(i);
 wavefun_1p5_11a(i)=sqrt(rho_1p5_prob11a(i))*grating_operator(i);
 wavefun_1p5_11b(i)=sqrt(rho_1p5_prob11b(i))*grating_operator(i);
 wavefun_1p5_12a(i)=sqrt(rho_1p5_prob12a(i))*grating_operator(i);
 wavefun_1p5_12b(i)=sqrt(rho_1p5_prob12b(i))*grating_operator(i);
 wavefun_1p5_13a(i)=sqrt(rho_1p5_prob13a(i))*grating_operator(i);
 wavefun_1p5_13b(i)=sqrt(rho_1p5_prob13b(i))*grating_operator(i);
 wavefun_1p5_14a(i)=sqrt(rho_1p5_prob14a(i))*grating_operator(i);
 wavefun_1p5_14b(i)=sqrt(rho_1p5_prob14b(i))*grating_operator(i);
 wavefun_1p5_15a(i)=sqrt(rho_1p5_prob15a(i))*grating_operator(i);
 wavefun_1p5_15b(i)=sqrt(rho_1p5_prob15b(i))*grating_operator(i);
 wavefun_1p5_16a(i)=sqrt(rho_1p5_prob16a(i))*grating_operator(i);
 wavefun_1p5_16b(i)=sqrt(rho_1p5_prob16b(i))*grating_operator(i);
 wavefun_1p5_17a(i)=sqrt(rho_1p5_prob17a(i))*grating_operator(i);
 wavefun_1p5_17b(i)=sqrt(rho_1p5_prob17b(i))*grating_operator(i);
 wavefun_1p5_18a(i)=sqrt(rho_1p5_prob18a(i))*grating_operator(i);
 wavefun_1p5_18b(i)=sqrt(rho_1p5_prob18b(i))*grating_operator(i);
 wavefun_1p5_19a(i)=sqrt(rho_1p5_prob19a(i))*grating_operator(i);
 wavefun_1p5_19b(i)=sqrt(rho_1p5_prob19b(i))*grating_operator(i);
 wavefun_1p5_20a(i)=sqrt(rho_1p5_prob20a(i))*grating_operator(i);
 wavefun_1p5_20b(i)=sqrt(rho_1p5_prob20b(i))*grating_operator(i); 
 wavefun_1p5_21a(i)=sqrt(rho_1p5_prob21a(i))*grating_operator(i);
 wavefun_1p5_21b(i)=sqrt(rho_1p5_prob21b(i))*grating_operator(i);
 wavefun_1p5_22a(i)=sqrt(rho_1p5_prob22a(i))*grating_operator(i);
 wavefun_1p5_22b(i)=sqrt(rho_1p5_prob22b(i))*grating_operator(i);
 wavefun_1p5_23a(i)=sqrt(rho_1p5_prob23a(i))*grating_operator(i);
 wavefun_1p5_23b(i)=sqrt(rho_1p5_prob23b(i))*grating_operator(i);
 wavefun_1p5_24a(i)=sqrt(rho_1p5_prob24a(i))*grating_operator(i);
 wavefun_1p5_24b(i)=sqrt(rho_1p5_prob24b(i))*grating_operator(i);
 wavefun_1p5_25a(i)=sqrt(rho_1p5_prob25a(i))*grating_operator(i);
 wavefun_1p5_25b(i)=sqrt(rho_1p5_prob25b(i))*grating_operator(i);
 wavefun_1p5_26a(i)=sqrt(rho_1p5_prob26a(i))*grating_operator(i);
 wavefun_1p5_26b(i)=sqrt(rho_1p5_prob26b(i))*grating_operator(i);
 wavefun_1p5_27a(i)=sqrt(rho_1p5_prob27a(i))*grating_operator(i);
 wavefun_1p5_27b(i)=sqrt(rho_1p5_prob27b(i))*grating_operator(i);
 wavefun_1p5_28a(i)=sqrt(rho_1p5_prob28a(i))*grating_operator(i);
 wavefun_1p5_28b(i)=sqrt(rho_1p5_prob28b(i))*grating_operator(i);
 wavefun_1p5_29a(i)=sqrt(rho_1p5_prob29a(i))*grating_operator(i);
 wavefun_1p5_29b(i)=sqrt(rho_1p5_prob29b(i))*grating_operator(i);
 wavefun_1p5_30a(i)=sqrt(rho_1p5_prob30a(i))*grating_operator(i);
 wavefun_1p5_30b(i)=sqrt(rho_1p5_prob30b(i))*grating_operator(i);

end
%%
wavefun_1p5_0extended(1:1000)=0;
wavefun_1p5_0extended(1001:3001)=wavefun_1p5_0;
wavefun_1p5_0extended(3002:4001)=0;

wavefun_1p5_1aextended(1:1000)=0;
wavefun_1p5_1aextended(1001:3001)=wavefun_1p5_1a;
wavefun_1p5_1aextended(3002:4001)=0;
wavefun_1p5_1bextended(1:1000)=0;
wavefun_1p5_1bextended(1001:3001)=wavefun_1p5_1b;
wavefun_1p5_1bextended(3002:4001)=0;
wavefun_1p5_2aextended(1:1000)=0;
wavefun_1p5_2aextended(1001:3001)=wavefun_1p5_2a;
wavefun_1p5_2aextended(3002:4001)=0;
wavefun_1p5_2bextended(1:1000)=0;
wavefun_1p5_2bextended(1001:3001)=wavefun_1p5_2b;
wavefun_1p5_2bextended(3002:4001)=0;
wavefun_1p5_3aextended(1:1000)=0;
wavefun_1p5_3aextended(1001:3001)=wavefun_1p5_3a;
wavefun_1p5_3aextended(3002:4001)=0;
wavefun_1p5_3bextended(1:1000)=0;
wavefun_1p5_3bextended(1001:3001)=wavefun_1p5_3b;
wavefun_1p5_3bextended(3002:4001)=0;
wavefun_1p5_4aextended(1:1000)=0;
wavefun_1p5_4aextended(1001:3001)=wavefun_1p5_4a;
wavefun_1p5_4aextended(3002:4001)=0;
wavefun_1p5_4bextended(1:1000)=0;
wavefun_1p5_4bextended(1001:3001)=wavefun_1p5_4b;
wavefun_1p5_4bextended(3002:4001)=0;
wavefun_1p5_5aextended(1:1000)=0;
wavefun_1p5_5aextended(1001:3001)=wavefun_1p5_5a;
wavefun_1p5_5aextended(3002:4001)=0;
wavefun_1p5_5bextended(1:1000)=0;
wavefun_1p5_5bextended(1001:3001)=wavefun_1p5_5b;
wavefun_1p5_5bextended(3002:4001)=0;

wavefun_1p5_6aextended(1:1000)=0;
wavefun_1p5_6aextended(1001:3001)=wavefun_1p5_6a;
wavefun_1p5_6aextended(3002:4001)=0;
wavefun_1p5_6bextended(1:1000)=0;
wavefun_1p5_6bextended(1001:3001)=wavefun_1p5_6b;
wavefun_1p5_6bextended(3002:4001)=0;
wavefun_1p5_7aextended(1:1000)=0;
wavefun_1p5_7aextended(1001:3001)=wavefun_1p5_7a;
wavefun_1p5_7aextended(3002:4001)=0;
wavefun_1p5_7bextended(1:1000)=0;
wavefun_1p5_7bextended(1001:3001)=wavefun_1p5_7b;
wavefun_1p5_7bextended(3002:4001)=0;
wavefun_1p5_8aextended(1:1000)=0;
wavefun_1p5_8aextended(1001:3001)=wavefun_1p5_8a;
wavefun_1p5_8aextended(3002:4001)=0;
wavefun_1p5_8bextended(1:1000)=0;
wavefun_1p5_8bextended(1001:3001)=wavefun_1p5_8b;
wavefun_1p5_8bextended(3002:4001)=0;
wavefun_1p5_9aextended(1:1000)=0;
wavefun_1p5_9aextended(1001:3001)=wavefun_1p5_9a;
wavefun_1p5_9aextended(3002:4001)=0;
wavefun_1p5_9bextended(1:1000)=0;
wavefun_1p5_9bextended(1001:3001)=wavefun_1p5_9b;
wavefun_1p5_9bextended(3002:4001)=0;
wavefun_1p5_10aextended(1:1000)=0;
wavefun_1p5_10aextended(1001:3001)=wavefun_1p5_10a;
wavefun_1p5_10aextended(3002:4001)=0;
wavefun_1p5_10bextended(1:1000)=0;
wavefun_1p5_10bextended(1001:3001)=wavefun_1p5_10b;
wavefun_1p5_10bextended(3002:4001)=0;

wavefun_1p5_11aextended(1:1000)=0;
wavefun_1p5_11aextended(1001:3001)=wavefun_1p5_11a;
wavefun_1p5_11aextended(3002:4001)=0;
wavefun_1p5_11bextended(1:1000)=0;
wavefun_1p5_11bextended(1001:3001)=wavefun_1p5_11b;
wavefun_1p5_11bextended(3002:4001)=0;
wavefun_1p5_12aextended(1:1000)=0;
wavefun_1p5_12aextended(1001:3001)=wavefun_1p5_12a;
wavefun_1p5_12aextended(3002:4001)=0;
wavefun_1p5_12bextended(1:1000)=0;
wavefun_1p5_12bextended(1001:3001)=wavefun_1p5_12b;
wavefun_1p5_12bextended(3002:4001)=0;
wavefun_1p5_13aextended(1:1000)=0;
wavefun_1p5_13aextended(1001:3001)=wavefun_1p5_13a;
wavefun_1p5_13aextended(3002:4001)=0;
wavefun_1p5_13bextended(1:1000)=0;
wavefun_1p5_13bextended(1001:3001)=wavefun_1p5_13b;
wavefun_1p5_13bextended(3002:4001)=0;
wavefun_1p5_14aextended(1:1000)=0;
wavefun_1p5_14aextended(1001:3001)=wavefun_1p5_14a;
wavefun_1p5_14aextended(3002:4001)=0;
wavefun_1p5_14bextended(1:1000)=0;
wavefun_1p5_14bextended(1001:3001)=wavefun_1p5_14b;
wavefun_1p5_14bextended(3002:4001)=0;
wavefun_1p5_15aextended(1:1000)=0;
wavefun_1p5_15aextended(1001:3001)=wavefun_1p5_15a;
wavefun_1p5_15aextended(3002:4001)=0;
wavefun_1p5_15bextended(1:1000)=0;
wavefun_1p5_15bextended(1001:3001)=wavefun_1p5_15b;
wavefun_1p5_15bextended(3002:4001)=0;

wavefun_1p5_16aextended(1:1000)=0;
wavefun_1p5_16aextended(1001:3001)=wavefun_1p5_16a;
wavefun_1p5_16aextended(3002:4001)=0;
wavefun_1p5_16bextended(1:1000)=0;
wavefun_1p5_16bextended(1001:3001)=wavefun_1p5_16b;
wavefun_1p5_16bextended(3002:4001)=0;
wavefun_1p5_17aextended(1:1000)=0;
wavefun_1p5_17aextended(1001:3001)=wavefun_1p5_17a;
wavefun_1p5_17aextended(3002:4001)=0;
wavefun_1p5_17bextended(1:1000)=0;
wavefun_1p5_17bextended(1001:3001)=wavefun_1p5_17b;
wavefun_1p5_17bextended(3002:4001)=0;
wavefun_1p5_18aextended(1:1000)=0;
wavefun_1p5_18aextended(1001:3001)=wavefun_1p5_18a;
wavefun_1p5_18aextended(3002:4001)=0;
wavefun_1p5_18bextended(1:1000)=0;
wavefun_1p5_18bextended(1001:3001)=wavefun_1p5_18b;
wavefun_1p5_18bextended(3002:4001)=0;
wavefun_1p5_19aextended(1:1000)=0;
wavefun_1p5_19aextended(1001:3001)=wavefun_1p5_19a;
wavefun_1p5_19aextended(3002:4001)=0;
wavefun_1p5_19bextended(1:1000)=0;
wavefun_1p5_19bextended(1001:3001)=wavefun_1p5_19b;
wavefun_1p5_19bextended(3002:4001)=0;
wavefun_1p5_20aextended(1:1000)=0;
wavefun_1p5_20aextended(1001:3001)=wavefun_1p5_20a;
wavefun_1p5_20aextended(3002:4001)=0;
wavefun_1p5_20bextended(1:1000)=0;
wavefun_1p5_20bextended(1001:3001)=wavefun_1p5_20b;
wavefun_1p5_20bextended(3002:4001)=0;

wavefun_1p5_21aextended(1:1000)=0;
wavefun_1p5_21aextended(1001:3001)=wavefun_1p5_21a;
wavefun_1p5_21aextended(3002:4001)=0;
wavefun_1p5_21bextended(1:1000)=0;
wavefun_1p5_21bextended(1001:3001)=wavefun_1p5_21b;
wavefun_1p5_21bextended(3002:4001)=0;
wavefun_1p5_22aextended(1:1000)=0;
wavefun_1p5_22aextended(1001:3001)=wavefun_1p5_22a;
wavefun_1p5_22aextended(3002:4001)=0;
wavefun_1p5_22bextended(1:1000)=0;
wavefun_1p5_22bextended(1001:3001)=wavefun_1p5_22b;
wavefun_1p5_22bextended(3002:4001)=0;
wavefun_1p5_23aextended(1:1000)=0;
wavefun_1p5_23aextended(1001:3001)=wavefun_1p5_23a;
wavefun_1p5_23aextended(3002:4001)=0;
wavefun_1p5_23bextended(1:1000)=0;
wavefun_1p5_23bextended(1001:3001)=wavefun_1p5_23b;
wavefun_1p5_23bextended(3002:4001)=0;
wavefun_1p5_24aextended(1:1000)=0;
wavefun_1p5_24aextended(1001:3001)=wavefun_1p5_24a;
wavefun_1p5_24aextended(3002:4001)=0;
wavefun_1p5_24bextended(1:1000)=0;
wavefun_1p5_24bextended(1001:3001)=wavefun_1p5_24b;
wavefun_1p5_24bextended(3002:4001)=0;
wavefun_1p5_25aextended(1:1000)=0;
wavefun_1p5_25aextended(1001:3001)=wavefun_1p5_25a;
wavefun_1p5_25aextended(3002:4001)=0;
wavefun_1p5_25bextended(1:1000)=0;
wavefun_1p5_25bextended(1001:3001)=wavefun_1p5_25b;
wavefun_1p5_25bextended(3002:4001)=0;

wavefun_1p5_26aextended(1:1000)=0;
wavefun_1p5_26aextended(1001:3001)=wavefun_1p5_26a;
wavefun_1p5_26aextended(3002:4001)=0;
wavefun_1p5_26bextended(1:1000)=0;
wavefun_1p5_26bextended(1001:3001)=wavefun_1p5_26b;
wavefun_1p5_26bextended(3002:4001)=0;
wavefun_1p5_27aextended(1:1000)=0;
wavefun_1p5_27aextended(1001:3001)=wavefun_1p5_27a;
wavefun_1p5_27aextended(3002:4001)=0;
wavefun_1p5_27bextended(1:1000)=0;
wavefun_1p5_27bextended(1001:3001)=wavefun_1p5_27b;
wavefun_1p5_27bextended(3002:4001)=0;
wavefun_1p5_28aextended(1:1000)=0;
wavefun_1p5_28aextended(1001:3001)=wavefun_1p5_28a;
wavefun_1p5_28aextended(3002:4001)=0;
wavefun_1p5_28bextended(1:1000)=0;
wavefun_1p5_28bextended(1001:3001)=wavefun_1p5_28b;
wavefun_1p5_28bextended(3002:4001)=0;
wavefun_1p5_29aextended(1:1000)=0;
wavefun_1p5_29aextended(1001:3001)=wavefun_1p5_29a;
wavefun_1p5_29aextended(3002:4001)=0;
wavefun_1p5_29bextended(1:1000)=0;
wavefun_1p5_29bextended(1001:3001)=wavefun_1p5_29b;
wavefun_1p5_29bextended(3002:4001)=0;
wavefun_1p5_30aextended(1:1000)=0;
wavefun_1p5_30aextended(1001:3001)=wavefun_1p5_30a;
wavefun_1p5_30aextended(3002:4001)=0;
wavefun_1p5_30bextended(1:1000)=0;
wavefun_1p5_30bextended(1001:3001)=wavefun_1p5_30b;
wavefun_1p5_30bextended(3002:4001)=0;

%%

farfield_x=[-2000:2000]*.72; %far field domain to obtain a ff peak_periodicity of 72um

Y=fft(wavefun_1p5_0extended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_0=a0*fourier_result;

Y=fft(wavefun_1p5_1aextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_1a=a1*interp1(farfield_x-1*delx,fourier_result,farfield_x);
Y=fft(wavefun_1p5_1bextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_1b=a1*interp1(farfield_x+1*delx,fourier_result,farfield_x);

Y=fft(wavefun_1p5_2aextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_2a=a2*interp1(farfield_x-2*delx,fourier_result,farfield_x);
Y=fft(wavefun_1p5_2bextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_2b=a2*interp1(farfield_x+2*delx,fourier_result,farfield_x);

Y=fft(wavefun_1p5_3aextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_3a=a3*interp1(farfield_x-3*delx,fourier_result,farfield_x);
Y=fft(wavefun_1p5_3bextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_3b=a3*interp1(farfield_x+3*delx,fourier_result,farfield_x);

Y=fft(wavefun_1p5_4aextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_4a=a4*interp1(farfield_x-4*delx,fourier_result,farfield_x);
Y=fft(wavefun_1p5_4bextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_4b=a4*interp1(farfield_x+4*delx,fourier_result,farfield_x);

Y=fft(wavefun_1p5_5aextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_5a=a5*interp1(farfield_x-5*delx,fourier_result,farfield_x);
Y=fft(wavefun_1p5_5bextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_5b=a5*interp1(farfield_x+5*delx,fourier_result,farfield_x);

Y=fft(wavefun_1p5_6aextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_6a=a6*interp1(farfield_x-6*delx,fourier_result,farfield_x);
Y=fft(wavefun_1p5_6bextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_6b=a6*interp1(farfield_x+6*delx,fourier_result,farfield_x);

Y=fft(wavefun_1p5_7aextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_7a=a7*interp1(farfield_x-7*delx,fourier_result,farfield_x);
Y=fft(wavefun_1p5_7bextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_7b=a7*interp1(farfield_x+7*delx,fourier_result,farfield_x);

Y=fft(wavefun_1p5_8aextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_8a=a8*interp1(farfield_x-8*delx,fourier_result,farfield_x);
Y=fft(wavefun_1p5_8bextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_8b=a8*interp1(farfield_x+8*delx,fourier_result,farfield_x);

Y=fft(wavefun_1p5_9aextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_9a=a9*interp1(farfield_x-9*delx,fourier_result,farfield_x);
Y=fft(wavefun_1p5_9bextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_9b=a9*interp1(farfield_x+9*delx,fourier_result,farfield_x);

Y=fft(wavefun_1p5_10aextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_10a=a10*interp1(farfield_x-10*delx,fourier_result,farfield_x);
Y=fft(wavefun_1p5_10bextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_10b=a10*interp1(farfield_x+10*delx,fourier_result,farfield_x);

Y=fft(wavefun_1p5_11aextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_11a=a11*interp1(farfield_x-11*delx,fourier_result,farfield_x);
Y=fft(wavefun_1p5_11bextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_11b=a11*interp1(farfield_x+11*delx,fourier_result,farfield_x);

Y=fft(wavefun_1p5_12aextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_12a=a12*interp1(farfield_x-12*delx,fourier_result,farfield_x);
Y=fft(wavefun_1p5_12bextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_12b=a12*interp1(farfield_x+12*delx,fourier_result,farfield_x);

Y=fft(wavefun_1p5_13aextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_13a=a13*interp1(farfield_x-13*delx,fourier_result,farfield_x);
Y=fft(wavefun_1p5_13bextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_13b=a13*interp1(farfield_x+13*delx,fourier_result,farfield_x);

Y=fft(wavefun_1p5_14aextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_14a=a14*interp1(farfield_x-14*delx,fourier_result,farfield_x);
Y=fft(wavefun_1p5_14bextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_14b=a14*interp1(farfield_x+14*delx,fourier_result,farfield_x);

Y=fft(wavefun_1p5_15aextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_15a=a15*interp1(farfield_x-15*delx,fourier_result,farfield_x);
Y=fft(wavefun_1p5_15bextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_15b=a15*interp1(farfield_x+15*delx,fourier_result,farfield_x);

Y=fft(wavefun_1p5_16aextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_16a=a16*interp1(farfield_x-16*delx,fourier_result,farfield_x);
Y=fft(wavefun_1p5_16bextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_16b=a16*interp1(farfield_x+16*delx,fourier_result,farfield_x);

Y=fft(wavefun_1p5_17aextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_17a=a17*interp1(farfield_x-17*delx,fourier_result,farfield_x);
Y=fft(wavefun_1p5_17bextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_17b=a17*interp1(farfield_x+17*delx,fourier_result,farfield_x);

Y=fft(wavefun_1p5_18aextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_18a=a18*interp1(farfield_x-18*delx,fourier_result,farfield_x);
Y=fft(wavefun_1p5_18bextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_18b=a18*interp1(farfield_x+18*delx,fourier_result,farfield_x);

Y=fft(wavefun_1p5_19aextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_19a=a19*interp1(farfield_x-19*delx,fourier_result,farfield_x);
Y=fft(wavefun_1p5_19bextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_19b=a19*interp1(farfield_x+19*delx,fourier_result,farfield_x);

Y=fft(wavefun_1p5_20aextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_20a=a20*interp1(farfield_x-20*delx,fourier_result,farfield_x);
Y=fft(wavefun_1p5_20bextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_20b=a20*interp1(farfield_x+20*delx,fourier_result,farfield_x);


Y=fft(wavefun_1p5_21aextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_21a=a21*interp1(farfield_x-21*delx,fourier_result,farfield_x);
Y=fft(wavefun_1p5_21bextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_21b=a21*interp1(farfield_x+21*delx,fourier_result,farfield_x);

Y=fft(wavefun_1p5_22aextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_22a=a22*interp1(farfield_x-22*delx,fourier_result,farfield_x);
Y=fft(wavefun_1p5_22bextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_22b=a22*interp1(farfield_x+22*delx,fourier_result,farfield_x);

Y=fft(wavefun_1p5_23aextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_23a=a23*interp1(farfield_x-23*delx,fourier_result,farfield_x);
Y=fft(wavefun_1p5_23bextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_23b=a23*interp1(farfield_x+23*delx,fourier_result,farfield_x);

Y=fft(wavefun_1p5_24aextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_24a=a24*interp1(farfield_x-24*delx,fourier_result,farfield_x);
Y=fft(wavefun_1p5_24bextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_24b=a24*interp1(farfield_x+24*delx,fourier_result,farfield_x);

Y=fft(wavefun_1p5_25aextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_25a=a25*interp1(farfield_x-25*delx,fourier_result,farfield_x);
Y=fft(wavefun_1p5_25bextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_25b=a25*interp1(farfield_x+25*delx,fourier_result,farfield_x);

Y=fft(wavefun_1p5_26aextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_26a=a26*interp1(farfield_x-26*delx,fourier_result,farfield_x);
Y=fft(wavefun_1p5_26bextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_26b=a26*interp1(farfield_x+26*delx,fourier_result,farfield_x);

Y=fft(wavefun_1p5_27aextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_27a=a27*interp1(farfield_x-27*delx,fourier_result,farfield_x);
Y=fft(wavefun_1p5_27bextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_27b=a27*interp1(farfield_x+27*delx,fourier_result,farfield_x);

Y=fft(wavefun_1p5_28aextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_28a=a28*interp1(farfield_x-28*delx,fourier_result,farfield_x);
Y=fft(wavefun_1p5_28bextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_28b=a28*interp1(farfield_x+28*delx,fourier_result,farfield_x);

Y=fft(wavefun_1p5_29aextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_29a=a29*interp1(farfield_x-29*delx,fourier_result,farfield_x);
Y=fft(wavefun_1p5_29bextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_29b=a29*interp1(farfield_x+29*delx,fourier_result,farfield_x);
 
Y=fft(wavefun_1p5_30aextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_30a=a30*interp1(farfield_x-30*delx,fourier_result,farfield_x);
Y=fft(wavefun_1p5_30bextended);
farfieldintensity = real( Y .* conj( Y ) );
farfield_pattern=fftshift( farfieldintensity );
fourier_result=abs(farfield_pattern)/max(abs(farfield_pattern));
scaledresult_1p5_30b=a30*interp1(farfield_x+30*delx,fourier_result,farfield_x);

%%

totalscaledresult=scaledresult_1p5_0+scaledresult_1p5_1a+scaledresult_1p5_1b+scaledresult_1p5_2a+scaledresult_1p5_2b+...
                                     scaledresult_1p5_3a+scaledresult_1p5_3b+scaledresult_1p5_4a+scaledresult_1p5_4b+...
                                     scaledresult_1p5_5a+scaledresult_1p5_5b+scaledresult_1p5_6a+scaledresult_1p5_6b+...
                                     scaledresult_1p5_7a+scaledresult_1p5_7b+scaledresult_1p5_8a+scaledresult_1p5_8b+...
                                     scaledresult_1p5_9a+scaledresult_1p5_9b+scaledresult_1p5_10a+scaledresult_1p5_10b+...
                                     scaledresult_1p5_11a+scaledresult_1p5_11b+scaledresult_1p5_12a+scaledresult_1p5_12b+...
                                     scaledresult_1p5_13a+scaledresult_1p5_13b+scaledresult_1p5_14a+scaledresult_1p5_14b+...                                     
                                     scaledresult_1p5_15a+scaledresult_1p5_15b+scaledresult_1p5_16a+scaledresult_1p5_16b+...
                                     scaledresult_1p5_17a+scaledresult_1p5_17b+scaledresult_1p5_18a+scaledresult_1p5_18b+...
                                     scaledresult_1p5_19a+scaledresult_1p5_19b+scaledresult_1p5_20a+scaledresult_1p5_20b+...
                                     scaledresult_1p5_21a+scaledresult_1p5_21b+scaledresult_1p5_22a+scaledresult_1p5_22b+...                                     
                                     scaledresult_1p5_23a+scaledresult_1p5_23b+scaledresult_1p5_24a+scaledresult_1p5_24b+...
                                     scaledresult_1p5_25a+scaledresult_1p5_25b+scaledresult_1p5_26a+scaledresult_1p5_26b+...                                     
                                     scaledresult_1p5_27a+scaledresult_1p5_27b+scaledresult_1p5_28a+scaledresult_1p5_28b+...
                                     scaledresult_1p5_29a+scaledresult_1p5_29b+scaledresult_1p5_30a+scaledresult_1p5_30b;

%%
%plot(farfield_x,fourier_result_coherent,'k')
hold on
%plot(farfield_x,fourier_result_1p5_0,'b')
plot(farfield_x,totalscaledresult./max(totalscaledresult),'m')
xlim([-30 30])
hold off
%%

plot(totalscaledresult,'m')
hold on
xlim([1950 2050])
hold off


%%

format long
f = fit(transpose(farfield_x(1950:2050)),transpose(totalscaledresult(1950:2050)),'gauss1');
coeffvals = coeffvalues(f);
width=2*sqrt(log(2))*coeffvals(3);
length=.1*72/width % = 0.618384553064904

