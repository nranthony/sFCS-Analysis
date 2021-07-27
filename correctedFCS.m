function [cor, pwave] = correctedFCS(decay, macrot, microt, cx)
%CORRECTEDFCS Time symmetric autocorrelation function
%   decay --> histogram of photon counts binned into tcspc channels
%   macro --> macrotime values; integer value of which laser pulse the
%   count occured afterwards
%   micro --> for each macro item, a microtime integer value of which tcspc
%   bin the count was recorded
%   cx --> log spaced abcissa for acf; units of macro laser pulses

%   cor --> returned corrrected autocorrelation function; same length as
%   input cx
%   pwave --> currently unknown; maybe error in acf

% Resvised version of Procedure.ipf for use with Matlab, adapted from:
% Correction of the afterpulsing effect in fluorescence correlation spectroscopy using time symmetry analysis
% Kunihiko Ishii 1,2 and Tahei Tahara 1,2
% 1 Molecular Spectroscopy Laboratory, RIKEN, 2-1 Hirosawa, Wako, Saitama 351-0198, Japan
% 2 RIKEN Center for Advanced Photonics, RIKEN, 2-1 Hirosawa, Wako, Saitama 351-0198, Japan
% 2015 Optical Society of America
% 14 Dec 2015 | Vol. 23, No. 25 | DOI:10.1364/OE.23.032387 | OPTICS EXPRESS 32387


    macrot_dims = size(macrot);
    microt_dims = size(microt);
    tcspc_chns_dims = size(decay);
    cx_dims = size(cx);
    % check 2D arrays, that the second dim is 1, and that macro and micro
    % have the same number of elements
    if ( numel(macrot_dims) ~= 2 || numel(microt_dims) ~= 2 )
        fprintf("Microtime and Macrotime arrays should both be dimensions N x 1, with N photon events.");
    end
    if ( macrot_dims(2) ~= 1 || microt_dims(2) ~= 1 )
        fprintf("Microtime and Macrotime arrays should both be dimensions N x 1, with N photon events.");
    end
    if not( isequal( macrot_dims(1), microt_dims(1) ) )
        fprintf("Microtime and Macrotime arrays should be the same length with N photon events.");
    end
    % check deacy data is correct dimensions	
    if ( numel(tcspc_chns_dims) ~= 2 || tcspc_chns_dims(2) ~= 1 )
        fprintf("TCSPC decay array should be dimensions N x 1, with N TCSPC channels.");
    end
    % check cx data is correct dimensions
    if ( numel(cx_dims) ~= 2 || cx_dims(2) ~= 1 )
        fprintf("Abscissa array should be dimensions N x 1.");
    end    
    
    n = macrot_dims(1); % number of photon counts -> should equal sumdecay below
    %maxt = double(macrot(end)); % last element in macrotime, i.e. the largest timepoint of experiment
    maxt = double(macrot(end)-macrot(1));
    tcspc_chns = tcspc_chns_dims(1); % number of TCSPC channels
    nt = cx_dims(1); % nt for number of tau; number of points in ACF
    cor = zeros(nt-1, 1); % correlation output; note, one element shorter than abscissa
    pwave = zeros(nt-1, 1); % probability wave...  ??  unknown at this point
    
    % bmat = zeros(tcspc_chns, nt-1); % b matrix for more detailed error -
    % to be revisited
    
    decay = double(decay); % change to floating point for ease in below; more accurate for weighting
    sumdecay = sum(decay); % total number of photon counts
    decay_mask = repmat(decay, 1); % copy of decay for masking
    decay_mask(decay>0) = 1; % create mask of all non zero TCSPC channels
    
    ave = sumdecay / sum(decay_mask); % average number of photons in all TCSPC channels
    decay0sum = decay - decay_mask .* ave; % center around mean, sum = 0
    inv_sqrt_decay = ones(tcspc_chns,1); % create empty array of 1's for below
    inv_sqrt_decay = inv_sqrt_decay ./ sqrt(decay); % create inv relative sigma error for multiplying for speed
    inv_sqrt_decay(isinf(inv_sqrt_decay)) = 0; % set inf's to 0
    siga = decay0sum .* inv_sqrt_decay; % sigma weighted
    siga_norm = norm(siga)^2; % ||sig||^2  the sum of squared relative sigma errors
    
	parfor i = 1:nt-1  %  loop over each abcissa; each big delta T
		
        tau = double(cx(i));
		dtau = double(cx(i+1) - cx(i));
        % cast all values to correct type for c function
        [b1, b2] = calcSymmetricACF(uint32(tau), ...
                                    uint32(dtau), ...
                                    uint16(microt), ...
                                    uint64(macrot), ...
                                    uint16(tcspc_chns));
        sumb1 = sum(b1);
        b = double(b1 - b2);
        
        sigb = b .* inv_sqrt_decay; % sigma weighted        
        c = dot(siga, sigb) / siga_norm;
        
        tmp = sumb1 - c * sumdecay;        
        tmp = tmp .* maxt ./ ( n .* ( maxt - tau - 0.5 .* dtau ) );
        tmp = tmp .* maxt ./ (n .* dtau);

        cor(i) = tmp;  %  chop last cx value when plotting
        
		pwave(i) = c / dtau .* maxt ./ ( maxt - tau - 0.5 .* dtau);

	end
	

end

