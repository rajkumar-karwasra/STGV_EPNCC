clear
close all
addpath('NoiseEstimation')
addpath('utils')
addpath('RTV')

para.alpha = 0.10;%α
para.beta = 0.01;%β
para.delta = 1000;%δ
para.maxiter = 20;%K
para.vareps = 1e-2;%ɛ
para.debug = true;

for j = 1
    img1 = im2double(imread("dataset\5.png"));
    figure;imshow(img1);title("Input");
    hsvI=rgb2hsv(img1);
    gamma=log10(0.3)/log10(mean2(hsvI(:,:,3))); 
    if(gamma>1) 
        gamma=1; 
    end
    hsvI(:,:,3)=hsvI(:,:,3).^gamma; 
    E0=hsv2rgb(hsvI);  
    S0=tsmooth(E0,0.015,3);
 
    S00=S0;
    
    nSig=NoiseEstimation(E0,1);
    if(nSig>0.06)
        para.delta=0.5;
    else
        para.delta=1000;
    end
    
    %Target noise histogram statistics
    var1=0.1;
    N_target=var1^0.5*randn(size(img1));
    [C0,edges0] = histcounts(N_target,2000,'Normalization','cdf');
    [C,ia,ic]=unique(C0);
    C0=[0,C];
    C0(1,end)=1;
    edges01=edges0(2:end);
    edges0=[edges0(1),edges01(ia)];
    rng('default') 
    if(nSig>0.1)
        N_add=0.00001^0.5*randn(size(img1)); 
    else
        N_add =0.000001^0.5*randn(size(img1)); %for SDSD
    end

    maxiter1=3; 
    N2=zeros(size(img1));
    for iter=1:maxiter1
        hsvI=rgb2hsv(S0);
        hsvI(:,:,3)=hsvI(:,:,3).^(1/gamma);
        S_low=hsv2rgb(hsvI);
        if(iter==1)
            S_low0=S_low;
        end
        N1=img1-S_low;
        N1=N1+N_add;
        
        %Noise transformation
        for ch=1:size(img1,3)
           N2(:,:,ch)=Noise_trans(N1(:,:,ch),C0,edges0);
        end

        Nimg2=S0+N2;
        Nimg2(Nimg2<0)=0;
        Nimg2(Nimg2>1)=1;
        
        if(iter==maxiter1)

%% Adaptive ERGF Parameters

nSig = NoiseEstimation(Nimg2,5);

if nSig < 0.03

    sigma_s = 1;
    sigma_r = 0.005;
    iteration = 1;

elseif nSig < 0.06

    sigma_s = 1;
    sigma_r = 0.007;
    iteration = 1;

else

    sigma_s = 1;
    sigma_r = 0.009;
    iteration = 1;

end

GaussianPrecision = 0.02;

T_map = RollingGuidanceFilter( ...
            im2double(Nimg2), ...
            sigma_s, ...
            sigma_r, ...
            iteration, ...
            GaussianPrecision);

T_map = max(min(T_map,1),0);

break;


else
    T_map = Nimg2;

        end


%% Texture Transformation (Edge-Aware Adaptive)

T2 = T_map - S0;

% Prevent division by zero
X = T2 ./ (N2 + 1e-6);

% Compute edge strength from current structure
grayS = rgb2gray(S0);

[Gx,Gy] = imgradientxy(grayS);

edgeMap = sqrt(Gx.^2 + Gy.^2);

edgeMap = edgeMap ./ (max(edgeMap(:)) + eps);

% Adaptive weight
w = exp(-2*edgeMap);

% Expand to RGB
w = repmat(w,[1 1 3]);

% Adaptive texture update
T1 = w .* N1 .* X;

% Limit excessive updates
idx = abs(X)>1;

T1(idx)=T2(idx);


S0 = S0 + T1;

S0 = max(min(S0,1),0);

    end
 
    %STGV
    S_map=S00;
    [I,R,N_target]=STGV(img1,para,S_map,(T_map-S_map));
    I_en=I.^(1/2.2);
    E1=R.*I_en;  
    
    %NCC
    E2=NCC_Enhancement(img1);
    figure; imshow(E1); title('STGV');
    figure; imshow(E2); title('NCC');
    %Integration
    S1=S_estimate(E1,0.03,3);
    S2=S_estimate1(E2,0.015,3,T_map);
    F=S2+E1-S1;
    figure; imshow(F); title('Final');

end