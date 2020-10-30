function [A,W] = ssd(X, freq, Fs, varargin)
% The function implimenting Spatio-Spectral Decomposition

% INPUT
%       * X [time x channel]: multi-channel signalto be decomposed
%        freq -  3 x 2 matrix with the cut-off frequencies. 
%             First row: cut-off frequencies for band-pass of the to be extracted 
%             oscillations.
%             Second row: cut-off frequencies for the lowest and highest 
%             frequencies defining flanking intervals.
%             Third row: cut-off frequencies for the band-stop filtering of 
%             the central frequency process.
%       * Fs: sampling frequency
% Inputs in varargin (those marked with ** are mandatory):
%       ** 'method' (str): 'fft' (default) or 'time'
%       * 'KeepN' (int>0): number of desired components. defauls is not
%           to exclude any component
%       * 'freqmask' (str): the frequency mask to be applied to the frequency
%        binds. default is the butterworth filter. accepted values are
%        'kaiser', 'butter'
%       *(*) 'maskparams: the freqmask parameters. In case butter window is
%       selected, this argument is mandatory. accepted values are:
%           - kaiser window: a structure with fields winlen and shape.
%           default values are winlen=2*Fs, shape=4.
%           - butter window: strucure containing the filter coefficients.
%           it has fields b,a, b_f,a_f, b_s,a_s, which would be computed as
%           [b,a]=butter(2, freq(1,:)/(Fs/2));
%           [b_f,a_f]=butter(2, freq2(2,:)/(Fs/2));
%           [b_s,a_s]=butter(2, freq(3,:)/(Fs/2),'stop');

% OUTPUTs:
%       * A [channel x KeepN]: the activation pattern
%       * W: The spaitial filter
%
%
%Example:
% [A_fft, W_fft] = ssd(sig_elec, FiltFreq, Fs, 'method', 'fft', ...
%        'KeepN', SourceNum, 'freqmask', 'kaiser', 'maskparams', kaiser_params);
%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Reference (please cite):
% Nikulin VV, Nolte G, Curio G. A novel method for reliable and fast extraction
% of neuronal EEG/MEG oscillations on the basis of spatio-spectral decomposition.
% NeuroImage, 2011, 55: 1528-1535.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% * This functon is basically an extension to the ssd function by Haufe et
% al, 2014. please check the bbci toolbox (https://github.com/bbci/bbci_public)
% ~~~~~~~~
% * FFT part by (c) Mina Jamshidi Idaji
% https://github.com/minajamshidi
% (minajamshidi91@gmail.com)
%% versions
% V1 (of the fft part): Aug 2019
% TODO: hanle epoch indices
%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% * This code is supposed to perform good!  Please report any problem to 
% the corresponding author of the paper

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS 
% OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
% AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER 
% OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
% DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER 
% IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
% OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%
%% default values
method = 'fft';
freqmask = 'butter';
KeepN = NaN;
%% input check ------------------------------------
% ------------------------------------------------

% X 
if size(X,1)<size(X,2)
    warning('this function accepts raw data as time*channel arrangement. It seems the raw data is not this format. Therefore it is transposed')
    X = X';
end

if (rem(length(varargin),2) == 1)
    error('Optional parameters should always go by pairs');
else
    for j = 1:2:(length(varargin)-1)
        if ~ischar (varargin{j})
            error (['Unknown type of optional parameter name (parameter' ...
                ' names must be strings).']);
        end
        switch lower (varargin{j})
            case 'method'
                method = varargin{j+1};
            case 'keepn'
                KeepN = varargin{j+1};
            case 'freqmask'
                freqmask = varargin{j+1};
            case 'maskparams'
                maskparams = varargin{j+1};
        end
    end
end
if ~exist('Fs')
    error('No input for sampling frequency.')   
end
%% 
if strcmp(freqmask,'butter')
    if ~exist('maskparams')
        [b,a]=butter(2, freq(1,:)/(Fs/2));
        [b_f,a_f]=butter(2, freq(2,:)/(Fs/2));
        [b_s,a_s]=butter(2, freq(3,:)/(Fs/2),'stop');
    else      
        try
            b = maskparams.b; a = maskparams.a;
            b_f = maskparams.b_f; a_f = maskparams.a_f;
            b_s = maskparams.b_s; a_s = maskparams.a_s;
        catch
            error('invalid mask parameters')
        end
    end
elseif strcmp(freqmask,'kaiser')
    if ~exist('maskparams')
        shape = 4;
        winlen = 2*Fs;
    else
        try
            shape = maskparams.shape;
            winlen = maskparams.winlen;
        catch
            error('invalid mask parameters')
        end
    end
end
%%
if strcmp(method, 'time')
    [b,a]=butter(2, freq(1,:)/(Fs/2));
    [b_f,a_f]=butter(2, freq(2,:)/(Fs/2));
    [b_s,a_s]=butter(2, freq(3,:)/(Fs/2),'stop');
    
    
    % Covariance matrix for the center frequencies (signal)
    X_s = filtfilt(b,a,X);
    C_s = cov(X_s,1);
    
    % Covariance matrix for the flanking frequencies (noise)
    X_tmp = filtfilt(b_f,a_f,X);
    X_tmp = filtfilt(b_s,a_s,X_tmp);
    C_n = cov(X_tmp,1);
    clear X_tmp
elseif strcmp(method, 'fft')
    % fft info
    T = size(X,1);
    nfft = length(X);
    f = (0:nfft-1)*2/nfft*Fs/2;
    f = f-Fs/2;
    SIG = fft(X,nfft);
    if strcmp(freqmask, 'kaiser')
        M = winlen;
        w = kaiser(M,shape);
        W = fft(w,nfft);
        Hs = (abs(f)>=freq(1,1) & abs(f)<=freq(1,2));
        hs = conv(Hs,W,'same')';     
        Hn = (abs(f)>=freq(2,1) & abs(f)<=freq(3,1)) | (abs(f)>=freq(3,2) & abs(f)<=freq(2,2));
        hn = conv(Hn,W,'same')';
        Xs = SIG.*hs;
        Xn = SIG.*hn;
    elseif strcmp(freqmask, 'butter')
        hs = fft(impz(b,a),nfft);
        hn1 = fft(impz(b_f,a_f),nfft);
        hn2 = fft(impz(b_s,a_s),nfft);
        Xs = SIG.*hs;
        Xn = SIG.*hn1.*hn2;
    else
        error('invalid frequency mask specified. Accepted values are kaiser and butter.')
    end
      
    C_n = real(Xn'*Xn)/(T*nfft);
    C_s = real(Xs'*Xs)/(T*nfft);
    
else
    error('invalid method specified. Accepted values are fft and time.')
end
        

%% Generalized eigenvalue decomposition
% dim-reduction of X does not have full rank
C = C_s;
[V, D] = eig(C);
[ev_sorted, sort_idx] = sort(diag(D), 'descend');
V = V(:,sort_idx);
% compute an estimate of the rank of the data
tol = ev_sorted(1) * 10^-6;
r = sum(ev_sorted > tol);
if r < size(X,2)
%     fprintf('SSD: Input data does not have full rank. Only %d components can be computed.\n',r);
    M = V(:,1:r) * diag(ev_sorted(1:r).^-0.5);
else
    M = eye(size(X,2));
end


C_s_r = M' * C_s * M;
C_n_r = M' * C_n * M;
[W,D]= eig(C_s_r,C_s_r+C_n_r);
[~, sort_idx] = sort(diag(D), 'descend');
W = W(:,sort_idx);

if ~isnan(KeepN)
    W = W(:,1:KeepN);
end

W = M * W;
% A is the matrix with the patterns (in columns)
A = C * W / (W'* C * W);
A = A./sqrt(sum(A.^2));

end

