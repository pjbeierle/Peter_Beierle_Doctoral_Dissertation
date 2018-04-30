xslope=[266.5,292.75,320.5];
yslope=[287,286,282.5];

cftool(xslope,yslope)
%%
clc;

i004_nearsurf_Au=zeros(480,640);
for i=1:640
    for j=1:480
        gogo=(j-1)*640+i;
        i004_nearsurf_Au(j,i)=intensity(gogo);
    end
end
%%
contour(i004_nearsurf_Au)
%%
slope=-.084;
%%
i004_nearsurf_Au(173:174,53)=0;
%intmatrix(171:172,52)=0;
%intmatrix(131,216)=0;
%intmatrix(134,216)=0;
contour(log(i004_nearsurf_Au))
%%

%%
%Nice_Contour_Generator

Xscale = [1:640]*((72e-6)/27.1913)-7.73325e-4;
Yscale = 1.2e-6*[1:480]-3.45e-4;
centery=((1.2e-6*(317.5))-3.81e-4)-((1.2e-6*(240))-3.81e-4)
[mxx,myy] = meshgrid(Yscale,Xscale);
contour(myy,mxx,transpose((i004_nearsurf_Au)),1000)
%shading interp;
box on;
%view([0,90])
hold on
set(gca,'xtick',[])
set(gca,'ytick',[])
 
ylim([1.2e-6*(246)-3.45e-4 1.2e-6*(320)-3.45e-4])
%ylim([(1.2e-6)-3.81e-4 ((1.2e-6)*480)-3.81e-4]+centery)
xlim([-200e-6 200e-6])
hold off
 

%%
ylineout_004nearsurf_Au=zeros(1,480);
for i=1:480 
    for j=277:308
ylineout_004nearsurf_Au(i)=ylineout_004nearsurf_Au(i)+i004_nearsurf_Au(i,j);
    end
end
%%

clc;
figure
plot(1.2e-6*[1:480]-3.45e-4,ylineout_004nearsurf_Au./max(ylineout_004nearsurf_Au),'.k','markers',15) %246-320
hold on
plot(simulation_domain+6e-5,revised_nearsurf_simulation_Gold,'LineWidth',3)
xlim([1.2e-6*(246)-3.45e-4 1.2e-6*(320)-3.45e-4])
xlabel('Y (in microns)','FontSize',40)
ylabel('#e^{-}','FontSize',40)
title('Vertical Electron Lineout','FontSize',40)
legend('fit','experimental data')
%
hold off
%%
%creating lineouts in x direction

newdomain=zeros(640*480,2);

newdomain(:,1)=VarName1;
newdomain(:,2)=VarName2;
F = scatteredInterpolant(newdomain,intensity,'natural','nearest')

%%

clc;
xxdom=1:.1:640;
yydom=1:.1:480;
%red lines

for k=246
line1=slope.*(xxdom-xslope(2))+k-2;
line2=slope.*(xxdom-xslope(2))+k+2;

lineout_temp=zeros(1,numel(xxdom));
for i=1:numel(xxdom)
   for j=1:numel(yydom)
       if yydom(j)>=line1(i)
           if yydom(j)<=line2(i)
               lineout_temp(i)=lineout_temp(i)+F(yydom(j),xxdom(i));
           end
       end
   end
end
lineout_im004_Au(k,:)=lineout_temp;


end
beep on
beep

%%
clc;
100*delta/(1.9*(2*sqrt(2*log(2))))
3.1*(2*sqrt(2*log(2)))
%%
clc;
k=246;
%guessing section
%plot(lineout(345,:)./max(lineout(345,:)))
xdata=1:.1:640;
xdata2=1.00003:.10003:.10018*6391;

enamp= 0.95%.99%0.9509;
alpha=0.046%.055%0.0532;
v=289.4%339.1864;
c1=291.2%336.7331;
delta=27.30%25.7624;
%delta=27.4
sigma=2.1%1.9589;
Scl=4 %4.4265 %Scaling factor determined on 11/30/16 from "beam only" fit
a1=0.9%.9861; %Amplitude factor determined on 11/30/16 from "beam only" fit
bckamp=0.05;%0.0113;
bckwidth=0.0005;%0.0005;
bckmid=290.;%341.2022%350.0158;


%x=[enampout_007_1(k),alphaout_007_1(k),vout_007_1(k),c1out_007_1(k),dout_007_1(k),wout_007_1(k)/(2*sqrt(2*log(2))),bckamp1out_007_1(k),bckwidth1out_007_1(k),bckmid1out_007_1(k),a1_1(k)]
x=[enamp,alpha,v,c1,delta,sigma,bckamp,bckwidth,bckmid,a1];
guess=x(1)*((sin(x(2)*(xdata2-x(3)))./(x(2)*(xdata2-x(3)))).^2).*...
       ( x(10)*exp(-((xdata2-x(4)).^2)/(2*x(6)^2))+(1-x(10))*exp(-((xdata2-x(4)).^2)/(2*Scl*x(6)^2))...
        +x(10)*exp(-((xdata2-x(4)+x(5)).^2)/(2*x(6)^2))+(1-x(10))*exp(-((xdata2-x(4)+x(5)).^2)/(2*Scl*x(6)^2))...
        +x(10)*exp(-((xdata2-x(4)-x(5)).^2)/(2*x(6)^2))+(1-x(10))*exp(-((xdata2-x(4)-x(5)).^2)/(2*Scl*x(6)^2))...
        +x(10)*exp(-((xdata2-x(4)+2*x(5)).^2)/(2*x(6)^2))+(1-x(10))*exp(-((xdata2-x(4)+2*x(5)).^2)/(2*Scl*x(6)^2))...
        +x(10)*exp(-((xdata2-x(4)-2*x(5)).^2)/(2*x(6)^2))+(1-x(10))*exp(-((xdata2-x(4)-2*x(5)).^2)/(2*Scl*x(6)^2))...
        +x(10)*exp(-((xdata2-x(4)+3*x(5)).^2)/(2*x(6)^2))+(1-x(10))*exp(-((xdata2-x(4)+3*x(5)).^2)/(2*Scl*x(6)^2))...
        +x(10)*exp(-((xdata2-x(4)-3*x(5)).^2)/(2*x(6)^2))+(1-x(10))*exp(-((xdata2-x(4)-3*x(5)).^2)/(2*Scl*x(6)^2))...
        )+x(7)*exp(-x(8)*(xdata2-x(9)).^2);


plot(xdata2,guess)
hold on
plot(xdata,lineout_im004_Au(k,:)./max(lineout_im004_Au(k,:)),'r')

plot(xdata,x(7)*exp(-x(8)*(xdata2-x(9)).^2),'k')
xlim([235 355])
%xlim([290 320])
%xlim([1 640])
ylim([0 1.02])
hold off



%%
clc;
clock
%%

clc;
begintime=clock
xxdom=1:.1:640;
yydom=1:.1:480;
%alpha=alphaout1(321);
%k=256;
%x=[enampout_007_1(k),alphaout_007_1(k),vout_007_1(k),c1out_007_1(k),dout_007_1(k),wout_007_1(k)/(2*sqrt(2*log(2))),bckamp1out_007_1(k),bckwidth1out_007_1(k),bckmid1out_007_1(k),a1_1(k)]
x=[enamp,alpha,v,c1,delta,sigma,bckamp,bckwidth,bckmid,a1];
for k=246:320 %maxima to ~5% of maxima
%  k=322-running
%line1=slope.*(xxdom-xslope(2))+k-2;
%line2=slope.*(xxdom-xslope(2))+k+2;

%lineout_temp=zeros(1,numel(xxdom));
%for i=1:numel(xxdom)
%   for j=1:numel(yydom)
%       if yydom(j)>=line1(i)
%           if yydom(j)<=line2(i)
%               lineout_temp(i)=lineout_temp(i)+F(yydom(j),xxdom(i));
%           end
%       end
%   end
%end
%lineout_im004_Au(k,:)=lineout_temp;

ydata=lineout_im004_Au(k,:)./max(lineout_im004_Au(k,:));


%x=[enamp,alpha,v,c1,delta,sigma,bckamp,bckwidth,bckmid,bckamp2,bckwidth2,bckmid2];

%lower and upper bounds force background to not fit peaks, or become
%negative
if k<316
lb=[-inf,-inf,-inf,-inf,-inf,-inf,0,0,0,.50];
ub=[inf,inf,inf,inf,inf,inf,inf,0.0007,inf,.95];
end
if k>=316
lb=[-inf,-inf,-inf,-inf,-inf,-inf,0,0,0,.01];
ub=[inf,inf,inf,inf,inf,inf,inf,0.0007,inf,.5];
end

%lb=[-inf,-inf,-inf,-inf,-inf,-inf,0,0,0,.50];
%ub=[inf,inf,inf,inf,inf,inf,inf,0.0007,inf,.95];

f = @(x,xdata2)x(1)*((sin(x(2)*(xdata2-x(3)))./(x(2)*(xdata2-x(3)))).^2).*...
       ( x(10)*exp(-((xdata2-x(4)).^2)/(2*x(6)^2))+(1-x(10))*exp(-((xdata2-x(4)).^2)/(2*Scl*x(6)^2))...
        +x(10)*exp(-((xdata2-x(4)+x(5)).^2)/(2*x(6)^2))+(1-x(10))*exp(-((xdata2-x(4)+x(5)).^2)/(2*Scl*x(6)^2))...
        +x(10)*exp(-((xdata2-x(4)-x(5)).^2)/(2*x(6)^2))+(1-x(10))*exp(-((xdata2-x(4)-x(5)).^2)/(2*Scl*x(6)^2))...
        +x(10)*exp(-((xdata2-x(4)+2*x(5)).^2)/(2*x(6)^2))+(1-x(10))*exp(-((xdata2-x(4)+2*x(5)).^2)/(2*Scl*x(6)^2))...
        +x(10)*exp(-((xdata2-x(4)-2*x(5)).^2)/(2*x(6)^2))+(1-x(10))*exp(-((xdata2-x(4)-2*x(5)).^2)/(2*Scl*x(6)^2))...
        +x(10)*exp(-((xdata2-x(4)+3*x(5)).^2)/(2*x(6)^2))+(1-x(10))*exp(-((xdata2-x(4)+3*x(5)).^2)/(2*Scl*x(6)^2))...
        +x(10)*exp(-((xdata2-x(4)-3*x(5)).^2)/(2*x(6)^2))+(1-x(10))*exp(-((xdata2-x(4)-3*x(5)).^2)/(2*Scl*x(6)^2))...
        )+x(7)*exp(-x(8)*(xdata2-x(9)).^2);
    
x = lsqcurvefit(f,x,xdata2,ydata,lb,ub);
x2=[x(4),x(5),x(6)];

f2 = @(x2,xdata2)x(1)*((sin(x(2)*(xdata2-x(3)))./(x(2)*(xdata2-x(3)))).^2).*...
       ( x(10)*exp(-((xdata2-x2(1)).^2)/(2*x2(3)^2))+(1-x(10))*exp(-((xdata2-x2(1)).^2)/(2*Scl*x2(3)^2))...
        +x(10)*exp(-((xdata2-x2(1)+x2(2)).^2)/(2*x2(3)^2))+(1-x(10))*exp(-((xdata2-x2(1)+x2(2)).^2)/(2*Scl*x2(3)^2))...
        +x(10)*exp(-((xdata2-x2(1)-x2(2)).^2)/(2*x2(3)^2))+(1-x(10))*exp(-((xdata2-x2(1)-x2(2)).^2)/(2*Scl*x2(3)^2))...
        +x(10)*exp(-((xdata2-x2(1)+2*x2(2)).^2)/(2*x2(3)^2))+(1-x(10))*exp(-((xdata2-x2(1)+2*x2(2)).^2)/(2*Scl*x2(3)^2))...
        +x(10)*exp(-((xdata2-x2(1)-2*x2(2)).^2)/(2*x2(3)^2))+(1-x(10))*exp(-((xdata2-x2(1)-2*x2(2)).^2)/(2*Scl*x2(3)^2))...
        +x(10)*exp(-((xdata2-x2(1)+3*x2(2)).^2)/(2*x2(3)^2))+(1-x(10))*exp(-((xdata2-x2(1)+3*x2(2)).^2)/(2*Scl*x2(3)^2))...
        +x(10)*exp(-((xdata2-x2(1)-3*x2(2)).^2)/(2*x2(3)^2))+(1-x(10))*exp(-((xdata2-x2(1)-3*x2(2)).^2)/(2*Scl*x2(3)^2))...
        )+x(7)*exp(-x(8)*(xdata2-x(9)).^2);
    
x2 = lsqcurvefit(f2,x2,xdata2,ydata);

result_004_Au_02(k,:)=x(1)*((sin(x(2)*(xdata2-x(3)))./(x(2)*(xdata2-x(3)))).^2).*...
       ( x(10)*exp(-((xdata2-x2(1)).^2)/(2*x2(3)^2))+(1-x(10))*exp(-((xdata2-x2(1)).^2)/(2*Scl*x2(3)^2))...
        +x(10)*exp(-((xdata2-x2(1)+x2(2)).^2)/(2*x2(3)^2))+(1-x(10))*exp(-((xdata2-x2(1)+x2(2)).^2)/(2*Scl*x2(3)^2))...
        +x(10)*exp(-((xdata2-x2(1)-x2(2)).^2)/(2*x2(3)^2))+(1-x(10))*exp(-((xdata2-x2(1)-x2(2)).^2)/(2*Scl*x2(3)^2))...
        +x(10)*exp(-((xdata2-x2(1)+2*x2(2)).^2)/(2*x2(3)^2))+(1-x(10))*exp(-((xdata2-x2(1)+2*x2(2)).^2)/(2*Scl*x2(3)^2))...
        +x(10)*exp(-((xdata2-x2(1)-2*x2(2)).^2)/(2*x2(3)^2))+(1-x(10))*exp(-((xdata2-x2(1)-2*x2(2)).^2)/(2*Scl*x2(3)^2))...
        +x(10)*exp(-((xdata2-x2(1)+3*x2(2)).^2)/(2*x2(3)^2))+(1-x(10))*exp(-((xdata2-x2(1)+3*x2(2)).^2)/(2*Scl*x2(3)^2))...
        +x(10)*exp(-((xdata2-x2(1)-3*x2(2)).^2)/(2*x2(3)^2))+(1-x(10))*exp(-((xdata2-x2(1)-3*x2(2)).^2)/(2*Scl*x2(3)^2))...
        )+x(7)*exp(-x(8)*(xdata2-x(9)).^2);

%chisqr calculation
%chisqr=0;
%norm_result=result_14(k,:)./sum(result_14(k,:));
%norm_ydata=ydata./sum(ydata);

%for i=1:numel(result_14(k,:))
%    chisqr=chisqr+((result(k,i)-ydata(i))^2)/result(k,i);
%     chisqr=chisqr+((norm_result(i)-norm_ydata(i))^2)/norm_result(i);
%end

%nu=numel(ydata)-(8); %9=free parameters/covariates 
%reduced_chisqr(k)=chisqr/nu;

%reset x (initial guess for each k)
x(4)=x2(1);
x(5)=x2(2);
x(6)=x2(3);

dout_004_Au_2(k)=x(5);
wout_004_Au_2(k)=2*sqrt(2*log(2))*x(6);

enampout_004_Au_2(k)=x(1);
alphaout_004_Au_2(k)=x(2);
vout_004_Au_2(k)=x(3);
c1out_004_Au_2(k)=x(4);
bckamp1out_004_Au_2(k)=x(7);
bckwidth1out_004_Au_2(k)=x(8);
bckmid1out_004_Au_2(k)=x(9);
a1_004_Au_2(k)=x(10);

k
end
beep on
beep
endtime=clock
x(2)

alpha
alpha=x(2)


%%
begintime-endtime
%%

%dout_007_1(256:334)*100./wout_007_1(256:334);
dout_004_Au_2(256)
wout_004_Au_2(256)
%%
(8.2152-5.2795)/5.2795
(25.7624-22.9957)/25.7624
%%
x(6)
wout(362)
%%
%figure %275-289
k=286; %246:320

x_lineout_data_004=result_004_Au_02(k,:);
xfit_004=lineout_im004_Au(k,:)./max(lineout_im004_Au(k,:));

plot(xdata*((72e-6)/27.1913)-7.73325e-4,x_lineout_data_004,'LineWidth',3)
hold on
plot(xdata*((72e-6)/27.1913)-7.73325e-4,xfit_004,'.k','markers',15)
xlim([-200e-6 200e-6])
ylim([0 1.02])
xlabel('Y (in microns)','FontSize',40)
ylabel('#e^{-}','FontSize',40)
title('Horizontal Electron Lineout','FontSize',40)
legend('fit','experimental data')
hold off

%%
%275-289
%290-308


plot(1.2e-6*[246:320]-3.45e-4,lcoh_im_004_Au_2,'.-k','markers',15)%564.125
hold on
plot(simulation_domain+6e-5,real_coherence_length_zurek_gold,'g','LineWidth',3)
plot(simulation_domain+6e-5,real_coherence_length_mach_gold,'LineWidth',3)
plot(simulation_domain+6e-5,real_coherence_length_scheel_gold,'r','LineWidth',3)
plot(simulation_domain+6e-5,real_coherence_length_mach_gold_dblcheck,'m','LineWidth',3)
plot(simulation_domain+6e-5,real_coherence_length_levinson_gold,'k','LineWidth',3)


xlim([1.2e-6*(246)-3.45e-4 1.2e-6*(320)-3.45e-4])
xlabel('Y (in pixels)','FontSize',40)
ylabel('L_{t} (in nm)','FontSize',40)
title('Gold','FontSize',40)
legend('experiment','Zurek','Machnikowski','Scheel & Buhmann','Machnikowski dblcheck')
ylim([0 680])
hold off

%%
mean(lcoh_im_004_Au_2)
%%
plot(wout_004_Au_2(246:320),'.k')
%%
enamp= 0.7669%.99%0.9509;
alpha=0.048%.055%0.0532;
v=337.1342%339.1864;
c1=336.9912%336.7331;
delta=26.0218%25.7624;
%delta=27.4
sigma=1.9%1.9589;
Scl=4 %4.4265 %Scaling factor determined on 11/30/16 from "beam only" fit
a1=0.5698%.9861; %Amplitude factor determined on 11/30/16 from "beam only" fit
bckamp=0.2034;%0.0113;
bckwidth=0.0005;%0.0005;
bckmid=340.5628;%341.2022%350.0158;
%%
%plot(283:364,alphaout_006_1(283:364),'.-k')
%plot(283:364,bckmid1out_006_1(283:364),'.-k')
%plot(283:364,bckwidth1out_037_1(283:364),'.-k')
plot(246:320,bckamp1out_004_Au_2(246:320),'.-k')
plot(246:320,a1_004_Au_2(246:320),'.-k')
hold on
%xlim([333.9 334.1])
hold off

%%
clc;
plot(1:86,lcoh_im038_01,'.-k')
hold on
plot([1:74]+11,lcoh_im034_01,'.-b')
hold off

%%

plot(1.3e-6*[1:480]-4.173e-4,ylineout_i007_nearsurf_h4p435_1000_36000./max(ylineout_i007_nearsurf_h4p435_1000_36000),'.-b')
hold on
plot(1.3e-6*[307:335]-4.173e-4,ylineout_i007_nearsurf_h4p435_1000_36000(307:335)./max(ylineout_i007_nearsurf_h4p435_1000_36000(307:335)),'.-m')
plot(1.3e-6*[300:345]-4.173e-4,lcoh_im038_01/500,'.-r')
plot(simulation_domain,1.03*edist_no_surf_no_exclusion./max(edist_no_surf_no_exclusion),'k')
plot(simulation_domain,1.03*edist_nosurf_exludedpaths./max(edist_nosurf_exludedpaths),'g')
xlim([-.5e-4 .5e-4])
ylim([0 1.1])
title('Crude Fit, Verticle Electron Distribution of 0th order')
xlabel('Position (in meters)')
ylabel('Intensity (normalized to max=1)')
legend('Experiment','Classical Simulation')
hold off

%%

plot(1.3e-6*[307:335]-4.18e-4,ylineout_i007_nearsurf_h4p435_1000_36000(307:335)./max(ylineout_i007_nearsurf_h4p435_1000_36000(307:335)),'.-b')
hold on
plot(simulation_domain,1.03*edist_no_surf_no_exclusion./max(edist_no_surf_no_exclusion),'k')
xlim([-1e-4 1e-4])
ylim([0 1.1])
title('Crude Fit, Verticle Electron Distribution of 0th order')
xlabel('Position (in meters)')
ylabel('Intensity (normalized to max=1)')
legend('Experiment','Classical Simulation')
hold off

%%

mean(dout_004_Au_2(246:320))


