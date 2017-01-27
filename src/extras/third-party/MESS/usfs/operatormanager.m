function op =operatormanager(name,flagmulA,flagpremulA,flagpostmulA,flagmulE,flagpremulE,flagpostmulE,flagsize,flagpresize,flagpostsize,flagsolA,flagpresolA,flagpostsolA,flagsolE,flagpresolE,flagpostsolE,flagsolApE,flagpresolApE,flagpostsolApE,flagmulApE,flagpremulApE,flagpostmulApE,flaginit,flaginitres,flagpreinitres,flagpostinitres,flaggetritzvals)
%% function op = operatormanager(name,flagAadd,flagEadd,flagmulA,flagmulE,
%% flagsizeA,flagsizeE,flagsolA,flagsolE,flagsolApE,flagmulApE,flaginit,flaginitres)
%
%  Generate structure with function handles.
%
%  Return structure with functionhandles that are implemented in folder
%  name. The flags are optional and specifiy the function handles to
%  return.
%
%  Calling sequence:
%
%    op = operatormanager(name)
%    op = operatormanager(name,flagAadd,flagEadd,flagmulA,flagmulE,
%       flagsizeA,flagsizeE,flagsolA,flagsolE,flagsolApE,flagmulApE)
%
%  Input:
%
%    name    name of folder which contains implemented function handles
%    
%    flag*   numeric (optional) flag~=0 op contain * function handle 
%                               flag==0 op does not contain addition 
%                               function handle
%
%  Output:
%
%    op      struct contains function handles     
%

%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, see <http://www.gnu.org/licenses/>.
%
% Copyright (C) Jens Saak, Martin Koehler, Peter Benner and others 
%               2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
%

%% check inputarguments
na = nargin;
if(~(na==1 || na==27))
   error('MESS:control_data','Number of input arguments has to be 1 or 27.'); 
end

if(~ischar(name))
   error('MESS:control_data','input argument name has to be a string.');
end

if(na==27&&~(isnumeric(flagmulA)&&isnumeric(flagpremulA)&&isnumeric(flagpostmulA)...
   &&isnumeric(flagmulE)&&isnumeric(flagpremulE)&&isnumeric(flagpostmulE)&&...
   isnumeric(flagsize)&&isnumeric(flagpresize)&&isnumeric(flagpostsize)...
   &&isnumeric(flagsolA)&&isnumeric(flagpresolA)&&isnumeric(flagpostsolA)...
   &&isnumeric(flagsolE)&&isnumeric(flagpresolE)&&isnumeric(flagpostsolE)...
   &&isnumeric(flagsolApE)&&isnumeric(flagpresolApE)&&isnumeric(flagpostsolApE)...
   &&isnumeric(flagmulApE)&&isnumeric(flagpremulApE)&&isnumeric(flagpostmulApE)...
   &&isnumeric(flaginit)&&isnumeric(flaginitres)&&isnumeric(flagpreinitres)...
   &&isnumeric(flagpostinitres)&&isnumeric(flaggetritzvals)));
error('MESS:control_data','flag has to be numeric');
end

%% check path to function handles
pa = which('operatormanager');
pa = strrep(pa,'operatormanager.m','');
pa = strcat(pa,name);

%% check whether the directory exists
if(~(exist(pa,'dir')))
    error('MESS:control_data','theres no folder %s',path);
end

%% check function in path
funcs = {'mul_A','mul_A_pre','mul_A_post','mul_E','mul_E_pre','mul_E_post',...
    'size','size_pre','size_post','sol_A','sol_A_pre','sol_A_post',...
    'sol_E','sol_E_pre','sol_E_post','sol_ApE','sol_ApE_pre','sol_ApE_post'...
    ,'mul_ApE','mul_ApE_pre','mul_ApE_post','init','init_res','init_res_pre','init_res_post'};
for f = funcs
   if(~exist(strcat(pa,'/',f{1},'_',name,'.m'),'file'))
    error('MESS:check_data','file %s does not exist',strcat(f{1},'_',name,'.m'));
   end
end

%% create op struct
op.name = name;
if(na==1)
    for f = funcs
       eval(sprintf('op.%s = @%s;',f{1},strcat(f{1},'_',name))); 
    end
    % add optional ritz value function.
    if(exist(strcat(pa,'/','get_ritz_vals_',name,'.m'),'file'))
      eval(sprintf('op.%s = @%s;','get_ritz_vals',strcat('get_ritz_vals_',name)));
    end
end

%here it is necessary that order is not changed!
if(na==27)
    flags=[flagmulA,flagpremulA,flagpostmulA,flagmulE,flagpremulE,flagpostmulE,...
        flagsize,flagpresize,flagpostsize,flagsolA,flagpresolA,flagpostsolA,...
        flagsolE,flagpresolE,flagpostsolE,flagsolApE,flagpresolApE,flagpostsolApE,...
        flagmulApE,flagpremulApE,flagpostmulApE,flaginit,...
        flaginitres,flagpreinitres,flagpostinitres,flaggetritzvals];
    i=1;
    for f = funcs
       if(flags(i))
       eval(sprintf('op.%s = @%s;',f{1},strcat(f{1},'_',name)));
       end
       i = i+1;
    end
    if(flags(i))
       eval(sprintf('op.%s = @%s;','get_ritz_vals',strcat('get_ritz_vals_',name)));
    end
end

end
