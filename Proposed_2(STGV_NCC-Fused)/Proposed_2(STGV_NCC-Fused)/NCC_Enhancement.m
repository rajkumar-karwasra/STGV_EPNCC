function E2 = NCC_Enhancement(img1)

% ==========================================================
% Edge-Preserving Adaptive Neighborhood Color Correction (EP-ANCC)
% Conservative modification of original NCC
% ==========================================================

avg = mean2(im2double(img1));

IMAX = max(double(img1(:)));

input_img = double(img1)/IMAX;

%% ---------------------------------------------------------
% Gamma (close to original)
%% ---------------------------------------------------------

if avg < 0.25
    gam = 0.45;
elseif avg > 0.50
    gam = 0.75;
else
    gam = 0.60;
end

input_img = input_img.^gam;

%% ---------------------------------------------------------
% RGB Channels
%% ---------------------------------------------------------

R_channel = input_img(:,:,1);
G_channel = input_img(:,:,2);
B_channel = input_img(:,:,3);

%% ---------------------------------------------------------
% ORIGINAL NCC
%% ---------------------------------------------------------

NCC = 0.299*R_channel + ...
      0.587*G_channel + ...
      0.114*B_channel ...
      - max(input_img,[],3) ...
      + min(input_img,[],3);

%% ---------------------------------------------------------
% Edge-preserving filtering
%% ---------------------------------------------------------

sigmaNCC = std(NCC(:));

smooth = 3e-4 + 1e-3*sigmaNCC;

NCC = imguidedfilter(NCC,...
    'NeighborhoodSize',[21 21],...
    'DegreeOfSmoothing',smooth);

%% ---------------------------------------------------------
% Atmospheric Light (unchanged)
%% ---------------------------------------------------------

atmsphrclight = atmLight(input_img,NCC);

%% ---------------------------------------------------------
% Adaptive transmission
%% ---------------------------------------------------------

localStd = stdfilt(NCC,ones(5));

localStd = mat2gray(localStd);

%omega = 0.22 + 0.06*localStd;   % range ≈ 0.22–0.28
lambda_e = 0.02;
omega = 0.22 + lambda_e*localStd;

lsar = 1 - omega.*NCC;

lsar = imguidedfilter(lsar,...
    'NeighborhoodSize',[15 15],...
    'DegreeOfSmoothing',0.0005);

lsar = max(lsar,0.15);

%% ---------------------------------------------------------
% Recovery
%% ---------------------------------------------------------

t0 = max(0.10,mean(lsar(:))*0.5);

R_out = atmsphrclight(1) + ...
    (R_channel-atmsphrclight(1))./max(lsar,t0);

G_out = atmsphrclight(2) + ...
    (G_channel-atmsphrclight(2))./max(lsar,t0);

B_out = atmsphrclight(3) + ...
    (B_channel-atmsphrclight(3))./max(lsar,t0);

E2 = cat(3,R_out,G_out,B_out);

%% ---------------------------------------------------------
% Clip
%% ---------------------------------------------------------

E2 = min(max(E2,0),1);

end