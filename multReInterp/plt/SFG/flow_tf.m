function varargout = flow_tf(flwfilename)

% Signal flowgraph (SFG) analysis for discrete-time, continuous-time and 
% purely symbolic systems. 
%
% Usage:
% output = flow_tf(flwfilename)
% 
% The SFG is defined in a textual nodelist file, where the syntax has five
% main types
%
% * type 1: comments
% * Comment begins with an asterisk
%
% * type 2: Signal flow
% * the signal flow is presented as:
% innode outnode gain
% * gain can be symbolic or numeric (e.g. c*z/(z-1))
%
% * type 3: indicating transfer function nodes
% * transfer function command '.tf', if you need transfer function:
% .tf fromnode tonode
%
% * type 4: pre-initialization and 
% * type 5: post-initialization
% * '.pre' and '.post' commands are in MATLAB syntax.
%
% * '.pre'-commands: initialization for symbols e.g.
% .pre syms z c
%
% * '.post' command is optional. These help further analysis.
%  Three examples:
% .post z=tf('z',Ts);c=1; *  the complex LTI-object frequency-variable 
% .post s=tf('s')     * (continuous time)
% .post syms w;z=exp(j*w); * or just symbolic re-evaluation
% 
% OUTPUT H (There are '.tf' command(s) in the nodelist file)
%  Transfer functions in H.
%  A struct array including symbolic structure H.sym and optionally H.tf
%  as a result of .post command.
%  If there are more than one ".tf " -lines, H.sym{1} is the 
%  first-defined transfer function in the nodelist file.
%
% OUTPUT y (No '.tf' command in the nodelist file)
%  A symbolic system matrix
%  Can be e.g transformed into state-space presentation.
%  (This is not quite straightforward, but the data is there). 
%
% IN MATLAB, TYPE:
% >> H=flow_tf('filename.suffix');
% YOU CAN SEE THE TRANSFER FUNCTION(S) BY TYPING E.G:
% >> H.sym{1}
% 
% More info:
% [1] manual.pdf
% [2] Neitola, Rahkonen (2005):
%   "A Fully Automated Flowgraph Analysis Tool for Matlab", in proc. ECCTD2005
% [3] Ruston, Bordogna (1966): 
%   "Electric Networks: Functions, Filters, Analysis", Chap. 6
% [4] Vlach, Singal (1983): 
%   "Computer Methods for Circuit Analysis and Design", Chp.14
% 
% Authors:
%  Marko Neitola, Timo Rahkonen 
%   Electronics laboratory
%   Dept. of Electrical and Information Engineering
%   University of Oulu, Finland. 
% 
% See also FOPEN FREAD SYM SYMS TF ZPK STRUCT

if nargout == 0
    warning('missing output variable, result will be in variable ans')
    pause(5)
end
if nargin == 0
    error('enter filename as input')
end
if ischar(flwfilename)==0
    error(['enter filename as character array ' char(39) 'name.suffix' char(39) ' !'])
end

fid=fopen(flwfilename,'r');
try,a=fread(fid);catch,error('File not found. Perhaps a typo, missing/wrong suffix or wrong directory.'),end
fclose(fid);

%--------prepare file contents for parsing

mask=find(a==13);a(mask)=10;clear mask % inapproppriate ascii symbol
mask_sp2=findstr(char(a)',[char(32) char(32)]);a(mask_sp2)=[]; %remove [space space]
mask_sp1=findstr(char(a)',[char(32) char(10)]);a(mask_sp1)=[]; %remove  [space EOL]
mask_sp3=findstr(char(a)',[char(10) char(10)]);a(mask_sp3)=[]; %remove [EOL EOL]
mask_das=findstr(char(a)','* *');a(mask_das)=[]; %remove ['* *']

mask_com_start=findstr(char(a)',[char(32) '*']);   %this makes
for ind=length(mask_com_start):-1:1                %the
mtemp=findstr(char(a)',char(10));                  %commenting
meol =mtemp(find((mtemp-mask_com_start(ind))>0));  %after
mask_com_end=meol(1)-1;                            %a node-descripition line
a(mask_com_start(ind):mask_com_end)=[];            %possible
end
a=[10;a]; %add EOL to the beginning

%--------end prepare

%--------parse the flw-file for flows, post/pre commands and tf definitions
mask=find(a==10);
indds=1;initial_commands='';post_commands='';
fprintf('*Parsing file %s:\n',flwfilename)

for ind = 1:(length(mask))
    try
        LINE=((a(mask(ind):mask(ind+1))'));
        LINE(find(LINE==10))=[];
        LINE=char(LINE);
    catch
        LINE = char(a((mask(ind)+1):end)');
    end
    if isempty(LINE)==0
        mki = LINE(1);
        if mki ~= '*'
            fprintf('%s\n',LINE)
            if mki == '.'
                if LINE(2:4)=='pre'
                    initial_commands=[initial_commands, LINE(6:end),';', char(10) ];
                elseif LINE(2:5)=='post'
                    post_commands=[post_commands ,LINE(7:end),';', char(10)];
                elseif LINE(2:3)=='tf'
                    clear tmp
                    tmp=LINE(5:end);
                    v_out(indds)=sym(['n' tmp((find(tmp==' ')+1):(end))]); %end-1
                    v_in(indds)= sym(['n' tmp(1:(find(tmp==' ')-1))]);    %)-1))
                    indds=indds+1;
                end
            else % flows - symbolic variables named for node names and values 'nNAME'
                spc=find(LINE==' ');
                innod(ind) = sym(['n' (LINE(1:(spc(1)-1)))]); %n: in case nodename is a number...
                try
                    outnod(ind) = sym(['n' (LINE((spc(1)+1):(spc(2)-1)))]);
                catch
                    error_msg = ['missing a nodename or value: line ' int2str(ind) ' on file ' flwfilename];
                    error(error_msg)
                end
                value(ind) = str2sym((LINE((spc(2)+1):end)));
            end
        end
   end
end

%--------end parse
innod_s=innod;outnod_s=outnod;
% --------give input and output nodes a corresponding number--------
ym_ind=1;
nodelist = [];
for ind = 1:length(innod)
    temp1=char(outnod(ind));
    if temp1(1)=='n'
        if exist(char(outnod(ind)))~=1
            temp=sym(char(outnod(ind)));
            nodelist = [nodelist;temp];
            innod(find(innod==temp))=ym_ind; 
            outnod(find(outnod==temp))=ym_ind;
            ym_ind=ym_ind+1;
        end
    end
end
for ind = 1:length(innod)
    temp2=char(innod(ind));
    if temp2(1)=='n'
        if exist(char(innod(ind)))~=1
            temp=sym(char(innod(ind)));
            innod(find(innod==temp))=ym_ind; 
            outnod(find(outnod==temp))=ym_ind;
            ym_ind=ym_ind+1;
            nodelist = [nodelist;temp];
        end
    end
end

N_tf=indds-1;
innod=eval(innod);outnod=eval(outnod);

if max([length(innod) length(outnod)])<N_tf
    error('Number of transfer functions should be LESS than number of nodes')
end
% --------end numbering

%non-character-zeroes indicate an empty line
mask=find(innod==0);
innod(mask)=[];outnod(mask)=[];value(mask)=[];
innod_s(mask)=[];outnod_s(mask)=[];

%The most important 5 lines are here:
L=max([innod outnod]);          % transmittance [3] matrix size is LxL
ym=diag(sym('1')*ones(L,1));
for ind=1:length(innod)
    ym(innod(ind),outnod(ind)) = ym(innod(ind),outnod(ind))-value(ind);
end

%initials (run the ".pre" commands)
try
    eval(initial_commands),
catch
    error_msg=['Un-initialized variable found in file: ' flwfilename];
    error(error_msg)
end

try
    if N_tf>0
        Z=inv(ym); 
    end
 catch
     if not(isempty(findstr(char(ym(:)),'*')))
         error_msg=['Something needs to be un-commented (file: ' flwfilename ')'];
     else
         error_msg=['Un-initialized variable found in file: ' flwfilename];
     end
     error(error_msg)
end

tfmask = [];
% evaluate the symbolic transfer functions
if N_tf>0
for ind = 1:N_tf
    maskin=innod(find(innod_s==v_in(ind)));
    maskout=outnod(find(outnod_s==v_out(ind)));  %innod_s
    try
        h_symb{ind} = simplify(Z(maskin(1),maskout(1)));
        
        tfmask = [tfmask;maskin(1),maskout(1)];
    catch
        error('Error. False nodename OR missing a signal path')
    end
    htf.sym=h_symb;
end
else
    fprintf('\nMissing the .tf command line.'),htf=0;
    fprintf('\nNo transfer functions are calculated.\n')
end

% optional: calculate the tf object transfer functions
if length(post_commands)>0 && N_tf>0
	clear z s
    try
        eval(post_commands);
    catch
        if exist('lti')~=2
            warning('You do not have Control System Toolbox.')
        else
            warning('Cant evaluate your post command.')
        end
        return
    end
	for ind = 1:length(htf.sym)
        try
            eval(['char(htf.sym{' int2str(ind) '});'])
            try
                htf.tf{ind} = minreal(zpk(eval(eval(['char(htf.sym{' int2str(ind) '});'])))); % is tf object
            catch
                htf.tf{ind} = eval(eval(['char(htf.sym{' int2str(ind) '});'])); % is sym or vector or something
            end
        catch ME
            ME
            error('Error. Propably >1 or none variables in post tf evaluation (just one required in zpk-object).')
        end
	end
end

if N_tf>0
    if nargout == 1
        varargout{1}=htf;
    elseif nargout == 2
        varargout{1}=htf;
        varargout{2}=ym;
    elseif nargout == 3
        varargout{1}=htf;
        varargout{2}=ym;
        varargout{3}=tfmask;
    elseif nargout == 4
        varargout{1}=htf;
        varargout{2}=ym;
        varargout{3}=tfmask;
        varargout{4}=nodelist;
    end
else
    if nargout == 1
        varargout{1}=ym;
    elseif nargout == 2
        varargout{1}=ym;
        varargout{2}=tfmask;
    elseif nargout == 3
        varargout{1}=htf;
        varargout{1}=ym;
        varargout{2}=tfmask;
        varargout{3}=nodelist;
    end
end
