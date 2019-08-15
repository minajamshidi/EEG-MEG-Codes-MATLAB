function  [bad_cmp_active,bad_cmp_passive] = select_components(mixmat,icasig,Fs,chanlocs,comp_var,varargin)
% This is a function to plot the mixing patterns and spectrum of components of a
% source separation method, and select the bad components.

% ACTIVE and PASSIVE bad components: 
% - Active bad components are those which you wish to clean.
% - Passive bad components are those which you wanna save, but not clean.

%   * For selecting active (passive) bad components--> left (right) click 
% on the mouse at the area of the mixing pattern of the component you wanna
% select.


% for instructions how to use the function, look at section "Functions"
% -------------------------------------------------------------------------
% INPUTS:
% ~~~~~~~
%   - mixmat: mixing matrix -> channel*#component
%   - icasig = the components -> #component*time
%   - Fs = sampling rate
%   - chanlocs = channel locations -> from EEG-lab structure
%   - comp_var = variance explained by each component -> #component*1
% -------------------------------------------------------------------------
% OUTPUTS:
% ~~~~~~~~
%   - bad_cmp_active = indeces of the active bad components
%   - bad_cmp_passive = indeces of the passive bad components
% -------------------------------------------
% Functions:
% ~~~~~~~~~~
%   mouse click:
%       - left click
%       - right click
%   keyboard press:
%       - esc: finish the inspection process
%       - a,p  : undo a specific active (a) or passive (p) component selection.
%       u will be asked to enter the index of the component u selected by mistake.
%       - z  : undo the very last selection.
%           - Note: the memory of the function is 1-back. meaning that u
%           cannot press z mutiple times before selecting a new component,
%           or u cannot press z a% Mina Jamshidi @ MPI CBS, Leipzigfter modifying your selections with a or
%           p.
%       - s: u will be asked to give the component number(s) u wish to see the spectrum of the sensor space projection. for example if u give [1,2,3] as the inoput, u will see the spec of the sensor-space signal constructed from those components.
%       - t: u will be asked for a component number that u want to see the time course. I added this mainly for heart.
%       - f: u will see the spectrum of the good components in
%       sensor-space.
%       -c : spectrum of a specific channel
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% VERSIONS:
% ~~~~~~~~~
% V0 - 21/01/2019
% Mina Jamshidi Idaji (minajamshidi91@gmail.com)
% - Exceptions handled:
%   - selecting a component which is selected before.
% -------------------------------------------------------------------------
% V1 - 01/02/2019 - \Mina
% - pressing z: undoing the very last selection.
% - right and left click, for selecting active and passive bad components
% - pressing a: deleting a specific active component selected before.
% - pressing p: deleting a specific passive component selected before.
% 
% - Exceptions handled:
%   - pressing a or p, inputing a component index, which is not selected
%   before.
%   -pressing a or p, without any input
%   - pressing a or p, not giving an integer - exception i and j. but the
%   error is handled as the input is not selected before.
%   - pressing z, but not components selected.
%   - pressing z multiple times. 
%   - selecting the same value for passive and active
% -------------------------------------------------------------------------
% xxxxxxxxxx TODO xxxxxxx
% - Exceptions to be handled:
%   - click on the wrong place!  
%   - if the number of input figure handles are not the equal to the figures needed 
% *** MULTIPLE windows
% specmax and specmin
%  - any error happens while waiting for key press
% - this function is good for 32 cannels (or less), it should be modified for more
% channels
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS 
% OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
% AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER 
% OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
% DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER 
% IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
% OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
% DAMAGE.
%% varargin

FIG = [];
if (rem(length(varargin),2) == 1)
    error('Optional parameters should always go by pairs');
else
    for j = 1:2:(length(varargin)-1)
        if ~ischar (varargin{j})
            error (['Unknown type of optional parameter name (parameter' ...
                ' names must be strings).']);
        end
        switch lower (varargin{j})
            case 'fig'
                FIG = varargin{j+1}; 
        end
    end
end
%% spectrum calculations
[ChN,~] = size(icasig);
[X_spec,F] = pwelch(icasig',1*Fs,0,Fs,Fs);
X_spec = X_spec(F<=45,:);
F = F(F<=45);
%% plotting
plot_per_fig_1 = 5;
plot_per_fig_2 = 7;
plot_per_fig = plot_per_fig_1*plot_per_fig_2;
N_fig = ceil(ChN/plot_per_fig);
% FIG = cell(N_fig,1);
fprintf('Plotting... Just wait a few seconds!\n');

% modify this later - now works for 32 channel
if isempty(FIG)
    fig = figure('Name',sprintf('fastICA Results - plot %d',1), 'units','normalized','outerposition',[0 0 1 1]);
else
    fig = FIG{1};
    clf(fig,'reset'); pause(0.5);
    set(0, 'CurrentFigure', fig); clf(fig,'reset')
end
for n_fig = 1:N_fig
   
    
    n_subplot = (n_fig-1)*plot_per_fig+1 : n_fig*plot_per_fig;
    n_subplot(n_subplot>ChN) = [];
    for k_component = n_subplot
        k_sp = mod(k_component,plot_per_fig);
        k_sp = k_sp + (k_sp==0)*plot_per_fig;
        
        subplot(plot_per_fig_1, plot_per_fig_2*2, 2*(k_sp-1)+1)
        topoplot(mixmat(:,k_component),chanlocs);
        title(sprintf('var=%.2f', comp_var(k_component)))
        
        subplot(plot_per_fig_1,plot_per_fig_2*2,2*k_sp)
        plot(F,10*log10(X_spec(:,k_component)),'linewidth',2.5);
        specmax = max(10*log10(X_spec(:,k_component)));
        specmin = min(10*log10(X_spec(:,k_component)));
        hold on;
        plot([10,10],[specmin,specmax],'r--'); hold on;
        plot([20,20],[specmin,specmax],'g--');
        grid on,
        title(num2str(k_component))
    end
    drawnow;
    %FIG{n_fig} = fig; % \Mina. commented. but kept for the future!;)
    
end
fprintf('Enjoy cleaning!\n');
%% selection
bad_cmp_active = NaN(ChN,1);
bad_cmp_passive = NaN(ChN,1);
n_bc_active = 0;
n_bc_passive = 0;
last_selection_type = NaN;
flag_z_inactivate = 0;
set(0, 'CurrentFigure', fig);
while 1 == 1
    
    wait1 = waitforbuttonpress;
    switch wait1
        
        % Keyboard press --------------------------------------------------
        case 1 
            key = get(gcf,'currentcharacter');
            
            % the Esc key -------------------------------------------------
            if key==27 
                break
                
            % press z, undo the previous action ---------------------------  
            elseif key=='z'
                if ~flag_z_inactivate
                    if strcmp(last_selection_type,'active')
                        if n_bc_active==0
                            fprintf('Warning: you asked to undo! BUT there is no components selected\n');
                            Selected_Cmp_Info_Print(bad_cmp_active,bad_cmp_passive,n_bc_active,n_bc_passive);
                        else
                            bad_cmp_active(n_bc_active) = NaN;
                            n_bc_active = n_bc_active - 1;
                            flag_z_inactivate = 1;
                            fprintf('Undo to active compenents applied!\n');
                            Selected_Cmp_Info_Print(bad_cmp_active,bad_cmp_passive,n_bc_active,n_bc_passive);
                        end
                    elseif strcmp(last_selection_type,'passive')
                        if n_bc_passive==0
                            fprintf('Warning: you asked to undo! BUT there is no components selected\n');
                        else
                            bad_cmp_passive(n_bc_passive) = NaN;
                            n_bc_passive = n_bc_passive - 1;
                            flag_z_inactivate = 1;
                            fprintf('Undo to passive compenents applied!\n');
                            Selected_Cmp_Info_Print(bad_cmp_active,bad_cmp_passive,n_bc_active,n_bc_passive);
                        end
                    elseif isnan(last_selection_type)
                        fprintf('Warning: No component is already selected at all!!\n');
                        Selected_Cmp_Info_Print(bad_cmp_active,bad_cmp_passive,n_bc_active,n_bc_passive);
                    end
                else
                    fprintf('Ohhhh... My memory is only 1-back!!\n');
                    fprintf('Use p or a to delete an unwanted selected component!!\n');
                end
                
                
               
            % press a, delete from active selection -----------------------
            elseif key=='a' 
                in_del = input('which active selected component do u wanna delete?','s');
                in_del = str2double(in_del);
                if isnan(in_del)
                   fprintf('Warning: you asked to undo an active component, but u entered an invalid index!!\n');
                    Selected_Cmp_Info_Print(bad_cmp_active,bad_cmp_passive,n_bc_active,n_bc_passive);
                elseif ~numel(in_del)
                    fprintf('Warning: you asked to undo an active component, but u did not specify any!!\n');
                    Selected_Cmp_Info_Print(bad_cmp_active,bad_cmp_passive,n_bc_active,n_bc_passive);
                elseif ~ismember(in_del,bad_cmp_active)
                    fprintf('Warning: you asked to undo a component which is not selected!!\n');
                    Selected_Cmp_Info_Print(bad_cmp_active,bad_cmp_passive,n_bc_active,n_bc_passive);
                else
                    bad_cmp_active(bad_cmp_active==in_del) = [];
                    bad_cmp_active = [bad_cmp_active;NaN];
                    n_bc_active = n_bc_active - 1;
                    flag_z_inactivate = 1;
                    fprintf('An active component deleted!\n');
                    Selected_Cmp_Info_Print(bad_cmp_active,bad_cmp_passive,n_bc_active,n_bc_passive)
                end
                
           
            elseif key=='p' 
                in_del = input('which passive selected component do u wanna delete?','s');
                in_del = str2double(in_del);
                if isnan(in_del)
                   fprintf('Warning: you asked to undo a passive component, but u entered an invalid index!!\n');
                   Selected_Cmp_Info_Print(bad_cmp_active,bad_cmp_passive,n_bc_active,n_bc_passive);
                elseif ~numel(in_del)
                    fprintf('Warning: you asked to undo a passive component, but u did not specify any!!\n');
                    Selected_Cmp_Info_Print(bad_cmp_active,bad_cmp_passive,n_bc_active,n_bc_passive);
                elseif ~ismember(in_del,bad_cmp_passive)
                    fprintf('Warning: you asked to undo a component which is not selected!!\n');
                    Selected_Cmp_Info_Print(bad_cmp_active,bad_cmp_passive,n_bc_active,n_bc_passive);
                else
                    bad_cmp_passive(bad_cmp_passive==in_del) = [];
                    bad_cmp_passive = [bad_cmp_passive;NaN];
                    n_bc_passive = n_bc_passive - 1;
                    flag_z_inactivate = 1;
                    fprintf('A passive component deleted!\n');
                    Selected_Cmp_Info_Print(bad_cmp_active,bad_cmp_passive,n_bc_active,n_bc_passive)
                end
                
            % press s, plot the spec of sensor space sig of the selected component ----------------------       
            elseif key=='s'
                comp_n_spec = input('sensor-space spectrum of which component?');
                if ~isempty(comp_n_spec)
                    sig_n_sensor = mixmat(:,comp_n_spec)*icasig(comp_n_spec,:);
                    figure, plot_spec(sig_n_sensor',Fs,45,'freqmark',[10,20]);
                end
                set(0, 'CurrentFigure', fig);
                
                
            % press t, plot time course of the selected component ----------------------       
            elseif key=='t'
                comp_n_spec = input('time-course of which component?');
                figure('Name',sprintf('time course of component %d',comp_n_spec)),
                if ~isempty(comp_n_spec)
                    L = size(icasig,2);
                    TT = (0:L-1)/Fs;
                    plot(TT,icasig(comp_n_spec,:));
                    xlabel('time(s)');
                end
                set(0, 'CurrentFigure', fig);
                
            % press f, plot time spec of remiained components ----------------------       
            elseif key=='f'
                comp_n_spec = 1:size(icasig,1);
                IND_del = [bad_cmp_active;bad_cmp_passive];
                IND_del(isnan(IND_del)) = [];
                comp_n_spec(IND_del) = [];
                sig_n_sensor = mixmat(:,comp_n_spec)*icasig(comp_n_spec,:);
                figure,
                set(gcf,'position',[10,10,1200,500]);
                subplot(2,7,[1,10]),
                plot_spec(sig_n_sensor',Fs,45,'freqmark',[10,20]);
                fprintf('--var rejected = %.2f (total), %.2f(active), %.2f(passive)\n',sum(comp_var(IND_del)),sum(comp_var(bad_cmp_active(~isnan(bad_cmp_active)))),sum(comp_var(bad_cmp_passive(~isnan(bad_cmp_passive)))));
                
                X = mixmat*icasig; X = X';
                [b_alphat,a_alpha] = butter(2,[8,12]/Fs*2);
                [b_beta,a_beta] = butter(2,[16,24]/Fs*2);
                
                sig_beta = filtfilt(b_beta,a_beta,double(sig_n_sensor)');
                p_sig_beta = var(sig_beta);
                sig_alpha = filtfilt(b_alphat,a_alpha,double(sig_n_sensor)');
                p_sig_alpha = var(sig_alpha);
                X_beta = filtfilt(b_beta,a_beta,double(X));
                p_X_beta= var(X_beta);
                X_alpha = filtfilt(b_alphat,a_alpha,double(X));
                p_X_alpha = var(X_alpha);
                
                subplot(2,7,[4,12]),
                stem(p_sig_alpha./p_X_alpha), title('clean/dirty alpha power')
                hold on, plot(1:size(X,2),0.95*ones(1,size(X,2)),'r--')
                hold on, plot(1:size(X,2),0.8*ones(1,size(X,2)),'g--')
                subplot(2,7,[6,14]),
                stem(p_sig_beta./p_X_beta), title('clean/dirty beta power')
                hold on, plot(1:size(X,2),0.95*ones(1,size(X,2)),'r--')
                hold on, plot(1:size(X,2),0.8*ones(1,size(X,2)),'g--')

                set(0, 'CurrentFigure', fig);
            % press c, plot spec of a specific chanenl ----------------------       
            elseif key=='c'
                chan_n_spec = input('spectrum of which channel in the sensor-space?');
                comp_n_spec = 1:size(icasig,1);
                IND_del = [bad_cmp_active;bad_cmp_passive];
                IND_del(isnan(IND_del)) = [];
                comp_n_spec(IND_del) = [];
                sig_n_sensor = mixmat(:,comp_n_spec)*icasig(comp_n_spec,:);
                figure, plot_spec(sig_n_sensor(chan_n_spec,:)',Fs,45,'freqmark',[10,20]);
                set(0, 'CurrentFigure', fig);
            end
            
        % Mouse click -----------------------------------------------------
        case 0 % mouse click
            f = gcf;
            axesClicked = gca;
            allAxes = findobj(f.Children,'Type','axes');
            numClicked = find(axesClicked==allAxes);
            nn_temp = length(n_subplot):-1:1;
            n_selected = n_subplot(nn_temp(numClicked/2));
            which_key = get(f,'SelectionType');
            
            % left click -------------------------------------------------- 
            if strcmp(which_key,'normal') 
                if ismember(n_selected,bad_cmp_passive)
                    fprintf('Warning: this component is already selected as passive bad component!\n');
                    fprintf('If you wish to select it as active bad component, first delete it from passive bad components by pressing p.\n');
                    Selected_Cmp_Info_Print(bad_cmp_active,bad_cmp_passive,n_bc_active,n_bc_passive);
                elseif ~ismember(n_selected,bad_cmp_active)
                    last_selection_type = 'active';
                    n_bc_active = n_bc_active + 1;                    
                    bad_cmp_active(n_bc_active) = n_selected;
                    flag_z_inactivate = 0;
                    Selected_Cmp_Info_Print(bad_cmp_active,bad_cmp_passive,n_bc_active,n_bc_passive)
                else
                    fprintf('Warning: this active component is selected before!\n');
                    Selected_Cmp_Info_Print(bad_cmp_active,bad_cmp_passive,n_bc_active,n_bc_passive);
                end
            
            % right click -------------------------------------------------- 
            elseif strcmp(which_key,'alt') 
                if ismember(n_selected,bad_cmp_active)
                    fprintf('Warning: this component is already selected as active bad component!\n');
                    fprintf('If you wish to select it as passive bad component, first delete it from active bad components by pressing a.\n');
                    Selected_Cmp_Info_Print(bad_cmp_active,bad_cmp_passive,n_bc_active,n_bc_passive);
                elseif ~ismember(n_selected,bad_cmp_passive)
                    last_selection_type = 'passive';
                    n_bc_passive = n_bc_passive + 1;
                    flag_z_inactivate = 0;
                    bad_cmp_passive(n_bc_passive) = n_selected;
                    Selected_Cmp_Info_Print(bad_cmp_active,bad_cmp_passive,n_bc_active,n_bc_passive)
                else
                    fprintf('Warning: this passive component is selected before!\n');
                    Selected_Cmp_Info_Print(bad_cmp_active,bad_cmp_passive,n_bc_active,n_bc_passive);
                end
                
            end
            
    end
end


bad_cmp_active(n_bc_active+1:end) = [];
bad_cmp_passive(n_bc_passive+1:end) = [];

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function Selected_Cmp_Info_Print(bad_cmp_active,bad_cmp_passive,n_bc_active,n_bc_passive)

if n_bc_passive==0
    bad_cmp_passive_str = '-';
else
    bad_cmp_passive_str = sprintf('%d, ',bad_cmp_passive(1:n_bc_passive));
end
if n_bc_active==0
    bad_cmp_active_str = '-';
else
    bad_cmp_active_str = sprintf('%d, ',bad_cmp_active(1:n_bc_active));
end

fprintf('Active bad components until now are %s\n',bad_cmp_active_str);
fprintf('Passive bad components until now are %s\n',bad_cmp_passive_str);
fprintf('----------------\n');
end


