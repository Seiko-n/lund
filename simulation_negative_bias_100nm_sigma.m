%% Quantum well capacitance
n1 = 1;% the number of sub bands populated
n2 = 2;% the number of sub bands populated
m = 9.109*10^-31; % [kg] electron mass
ms = 0.041*m; % effective mass of electron in InGaAs
q = 1.602*10^-19; %[C] elementary charge
h = 1.054*10^-34; % [m^2kg/s] reduced plank constant
Cq1 = (q^2*ms)/(n1*pi*h^2)*10^-4; % [F/cm2] 1st quantum well capacitance 
Cq2 = (q^2*ms)/(n2*pi*h^2)*10^-4;% [F/cm2] 2nd quantum well capacitance 

%% Subband calculation
h_ch = 7*10^-9;%[m] hight of channel
w_ch = 100*10^-9;%[m] width of channel
E1 = (h^2/(2*ms))*(n1^2*pi^2/h_ch^2 + n1^2*pi^2/w_ch^2)/q; %[eV] 1st subband of 2D quantum well
E2 = (h^2/(2*ms))*(n2^2*pi^2/h_ch^2 + n2^2*pi^2/w_ch^2)/q; %[eV] 2nd subband of 2D quantum well

%% Oxide capacitance 
tox = 5*10^-9; % oxide thickness [m]
ep = 8.854*10^-12; %[Fm^-1] % permittivity of vacuum
epox = 18; % permittivity of oxide. assume HfO2
Cox = (ep*epox/tox)*10^-4; % oxide capacitance [F/cm2]


%% plot experiment data
%H2 = [0.0084,0.0148,0.018,0.0456,0.03,0.0529,0.0505,0.0681,0.0694,0.0949];
%Vov2 = [200,400,500,600,700,800,900,1000,1100,1200];

H3 = [0.0086,0.0143,0.0202,0.0193,0.0225,0.0284,0.0208,0.0239,0.026,0.0363,0.0360,0.0579,0.0582,0.0788,0.0665];
Vov3 = [100,150,200,250,300,350,400,450,500,550,600,650,700,750,800];

%H4 = [0.0304,0.0235,0.045,0.0350,0.0416,0.0324,0.033,0.0639,0.0769];
%Vov4 = [400,500,600,700,800,900,1000,1100,1200];

%H5 = [0.012,0.0234,0.026,0.038,0.0386,0.048,0.0683,0.0985];
%Vov5 = [500,600,700,800,900,1000,1100,1200];

%H6 = [0.0171,0.0164,0.0325,0.0294,0.02,0.038,0.0514,0.049,0.0573,0.0901];
%Vov6 = [300,400,500,600,700,800,900,1000,1100,1200];

plot((Vov3/1000),H3,'o','LineWidth',2);
%plot((Vov2/1000)-0.2,H2,'d',(Vov3/1000)-0.2,H3,'d',(Vov4/1000)-0.2,H4,'d',(Vov5/1000)-0.2,H5,'d',(Vov6/1000)-0.2,H6,'d');
hold on


%% spline
%sp_H = csaps((Vov3/1000),H3,0.999);
%nVov = 0.1:0.001:1;
%nH = ppval(sp_H,nVov);
%plot(nVov,nH,'--','LineWidth',2);
%hold on

%% Integration from negative bias
%tic;
%[Qit,Cit,Qqu,Cqu,Vox,Vov,Vg,dV] = mos_calculation(0.95-E1,0.17,-1-E1,0.5,0.4*10^14,-0.37,2,0.2,Cox,Cq1);
%toc;


%%
k =1;
S = 1000;
h = waitbar(0,'Please wait...');
for Dit = 0.5:0.05:0.8
    for mu1 = 0.9:0.05:1.1
        for mu2 = -0.9:-0.05:-1.1
            for sigma1 = 0.15:0.025:0.25
                for sigma2 = 0.5:0.025:0.6
                    tic;
                    [Qit,Cit,Qqu,Cqu,Vox,Vov,Vg,dV] = mos_calculation(mu1-E1,sigma1,mu2-E1,sigma2,Dit*10^14,0.001,1,0,Cox,Cq1);
                    sp_Hs = csaps(Vov,dV,0.999);
                    nHs = ppval(sp_Hs,Vov3/1000);
                    %plot(nVov,nHs,'--','LineWidth',2);
                    %hold on
                    E = sum(log((abs(nHs-H3)).^2));
                    if E < S;
                        S = E;
                        m1 = mu1;
                        m2 = mu2;
                        s1 = sigma1;
                        s2 = sigma2;
                        %fs = f_start;
                        D = Dit*10^14;
                    end
                    toc;
                end
            end
        end
    end
    waitbar(Dit/1)
end
close(h)
%toc;


%% Error estimation
                    %[Qit,Cit,Qqu,Cqu,Vox,Vov,Vg,dV] = mos_calculation(0.95-E1,0.17,-1-E1,0.5,0.7*10^14,0.001,1,0,Cox,Cq1);
                    %sp_Hs = csaps(Vov,dV,0.999);
                    %nHs = ppval(sp_Hs,nVov);
                    %plot(nVov,nHs,'--','LineWidth',2);
                    %set(gca, 'yscale','log')
                    %set(gca, 'xscale','log')
                    %axis([0.1,1,0.001,1])
                    %hold on
             
%% Plot Vg vs hysteresis
%plot(Vov,dV,'or','LineWidth',2);
%plot(Vov,H,'o','MarkerSize', 8,'MarkerFaceColor','r','MarkerEdgeColor','k');
%set(gca, 'yscale','log')
%set(gca, 'xscale','log')
%axis([-0.5,1,0.001,1])
%xlabel('V_{ov} = V_{g} - V_{th} (V)');
%ylabel('\Deltahysteresis(V)');
%set(gca,'FontName','Helvetica','FontSize',16); 
%legend('experiment data1(W = 25nm)','experiment data2(W = 25nm)','experiment data3(W = 25nm)','experiment data4(W = 25nm)','experiment data5(W = 25nm)','simulation','Location','SouthEast')
%legend('simulation','Location','SouthEast')
%set(gca, ...
%  'TickDir'     , 'out'     , ...
%  'TickLength'  , [.02 .02] , ...
%  'XMinorTick'  , 'on'      , ...
%  'YMinorTick'  , 'on'      , ...
%  'YGrid'       , 'on'      , ...
%  'XGrid'       , 'on'      , ...
%  'XColor'      , [.3 .3 .3], ...
%  'YColor'      , [.3 .3 .3], ...
%  'LineWidth'   , 2         );
%hold on
%
%saveas(gcf,'fitting_W=100nm-device_with_HfO2_distribution.png');




