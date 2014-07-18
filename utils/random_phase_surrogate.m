function [z_surrogate, Z_amps] = random_phase_surrogate(z, varargin)
% [z_surrogate, Z_amps] = random_phase_surrogate(z, varargin)
%
% Creates a shuffled (surrogate) version of the time series z that as the
% same auto-correlation (time-dynamics).
%
% In:
%   z - univariate time series
% 
% Out:
%   z_surrogate - time series with the same amplitude (power) spectrum as z
%                   but random phases.
%

opt= propertylist2struct(varargin{:});
opt= set_defaults(opt, ...
                  'z_amps',[]...
                  );
z = z(:)';
% compute the amplitude spectrum of z (if not already given)
n = length(z);
if isempty(opt.z_amps)
    Z = fft(z);
    Z_amps = abs(Z); % amplitude spectrum
else
    Z_amps = opt.z_amps;
end

% construct random phases
rand_phases = rand(1,floor(n/2))*2*pi;
start = length(rand_phases);
if mod(n,2) == 0
    start = start-1;
end
rand_phases = [0, rand_phases, -rand_phases(start:-1:1)];
% put amps and phases together for complex Fourier spectrum
Z_surrogate = Z_amps .* exp(1i*rand_phases);
% project the complex spectrum back to the time domain
z_surrogate = real(ifft(Z_surrogate));
