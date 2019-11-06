
clear

%range of shear velocities to model
ustar1 = linspace(0.05,0.75,200);

%range of equilibrium sizes to calculate
D1 = (200:100:17000).*10^-6;

%parameters and constants for abrasion rate calcs
kv = 90*10^4;  %dimensionless tuning parameter
rho_s = 2800; %[kg/m^3] density of aragonite
rho_w = 1025; %[kg/m^3] density of seawater
R = (rho_s - rho_w)/rho_w; %submerged specific density
g = 9.81; %[m/s^2]
nu = 0.937*10^-6; %[m^2/s] kinematic viscosity of water
young = 20*10^9; %[Pa] young's modulus
strength = 1*10^6; %[Pa] tensile strength
tauc = 0.043; %Critical Shields number.  
Stc = 9; %critical Stokes threshold for viscous damping
A1 = 0.36; %dimensionless constant for transport calcs
intermittency = .01;
%in general larger bedforms = less frequent movement = smaller value for
%intermittency in the range (0 1]
H = 5;  %[m] water depth

%parameters and constants for precipitation rate calcs
Omega1 = 1:24;
n = 2.26; %empirical values from Zhong & Mucci 1989
k = 10^1.11; %[umol/m^2/hr]
%preallocate space for the calculate Deq's
D_fOmega = zeros(length(ustar1),length(Omega1));

%calculate settling velocities and drag coeficients
CSF = 1;  %1 is for spheres, 0.8 is for natural
PS = 6;  %6 is for spheres, 3.5 is for natural
Dstar = (R.*g.*D1.^3)./(nu.^2);
X = log10(Dstar);
R1 = -3.76715+1.92944.*X - 0.09815.*(X.^2) - 0.00575.*(X.^3) + 0.00056.*(X.^4);
R2 = log10(1-((1-CSF)./0.85))-(((1-CSF).^2.3).*tanh(X-4.6)) + 0.3.*(0.5-CSF).*((1-CSF).^2).*(X-4.6);
R3 = (0.65-((CSF./2.83).*tanh(X-4.6))).^(1+((3.5-PS)./2.5));
Wstar = R3.*10.^(R2+R1);
ws1 = (R.*g.*nu.*Wstar).^(1./3); %[m/s]
cdrag1 = (4/3).*(R.*g.*D1)./(ws1.^2);

for counter1 = 1:length(Omega1)
    Omega = Omega1(counter1);
    Precip_rate = k*(Omega - 1)^n; %[umol/m^2/hr]
    SA_ssa = pi().*D1.^2.*23; %specific surface area estimate
    Precip_rate_vol = Precip_rate.*10^-6.*10^-3.*100.0869./rho_s.*SA_ssa;

    for counter2 = 1:length(ustar1)
        ustar = ustar1(counter2);
        Rabrasion_mat = zeros(1,length(D1));
        
        for m = 1:length(D1)
            D = D1(m);
            ws = ws1(m);
            cdrag = cdrag1(m);
            tau = ustar^2/(R*g*D);
            tstage = tau/tauc;
            susp_abrasion_calculations_giantooids
            Rabrasion_mat(m) = Rabrasion;
        end
        
        dVA = 4*pi().*(D1./2).^2.*Rabrasion_mat; %[m^3/yr]

        differences = (Precip_rate_vol - dVA.*intermittency);

        [mm,index] = min(differences>0);
        D_fOmega(counter2,counter1) = D1(index);

    end
end

figure
contourf(ustar1,Omega1,10^3*transpose(D_fOmega),1:14,'ShowText','on')
xlabel('u_* (m/s)')
ylabel('\Omega_a_r')
c = colorbar;
c.Label.String = 'equilibrium ooid size (mm)';
xlim([0.15 0.75])
caxis([1 14])
