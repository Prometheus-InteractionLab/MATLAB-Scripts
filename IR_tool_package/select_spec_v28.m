function varargout = select_spec_v28(varargin)
% SELECT_SPEC_V28 MATLAB code for select_spec_v28.fig
%      SELECT_SPEC_V28, by itself, creates a new SELECT_SPEC_V28 or raises the existing
%      singleton*.
%
%      H = SELECT_SPEC_V28 returns the handle to a new SELECT_SPEC_V28 or the handle to
%      the existing singleton*.
%
%      SELECT_SPEC_V28('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SELECT_SPEC_V28.M with the given input arguments.
%
%      SELECT_SPEC_V28('Property','Value',...) creates a new SELECT_SPEC_V28 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before select_spec_v28_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to select_spec_v28_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help select_spec_v28

% Last Modified by GUIDE v2.5 11-Dec-2014 18:04:44


% v17: 17/04/2013: Ecart type des mesures de pénétration vapeur

% v20: 16/04/2014 - LMM
%   - Ajout de la génération des profils radiaux de vitesse, de fraction de
%   mélange, et de température adiabatique de flamme

% v22: 01/09/2014 - LMM
%   - Modification du calcul de Zst

% v23, v24: 19/09/2014 - LMM
%   - Correction de bugs

% v25: 19/09/2014 - LMM
%   - Lifting en profondeur du code

% v26: 19/10/2014 - LMM
%   - Amélioration du calcul des champs de vitesse et de fraction de
%   mélange

% v27: 28/10/2014 - LMM
%   - Possibilité d'utiliser une origine "virtuelle" pour le jet, située en
%   aval de l'orifice réel

% v28: 11/12/2014 - LMM
%   - Possibilité d'avoir un angle vapeur variable

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @select_spec_v28_OpeningFcn, ...
    'gui_OutputFcn',  @select_spec_v28_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

try handles = varargin{end}; end
try compute_data_ambient(handles); end
try compute_data_fuel(handles); end
try compute_data_injector(handles); end
try compute_data_lol(handles); end
try compute_data_mixing(handles); end
try compute_data_phi(handles); end
try plot_velocities(handles); end


% #########################################################################################
% #################################### GUI FUNCTIONS ######################################
% #########################################################################################

function select_spec_v28_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to select_spec_v28 (see VARARGIN)

% Choose default command line output for select_spec_v28
handles.output = hObject;

% Initialization
handles.change = 0;

data = [];

data.injector   = [];
data.fuel       = [];
data.ambient    = [];
data.options0d  = [];
data.lol        = [];
data.mixing     = [];

% Injector
data.injector.casel = get(handles.casel,'Value');
data.injector.xjet_m = [];
data.injector.Rjet_m = [];

% Fuel
data.fuel.autorho = 0;

% Ambient
data.ambient.xo2sel = 1;
data.ambient.pambsel = 1;

% Options of 0d model
data.options0d.unif_profile = get(handles.gaussprof,'Value');
data.options0d.alpha = str2double(get(handles.alphaman,'String'));
data.options0d.a_Siebers = str2double(get(handles.sieberscoeff,'String'));

% Lift-off length
data.lol.autolol = 0;
data.lol.coeflol = str2double(get(handles.coeflolval,'String'));

% Vapor penetration
data.mixing.adjust_soi = 0;
data.mixing.plotapprox = 0;

setappdata(0,'ListDialogAppData__',data);

% Update handles structure
guidata(hObject, handles);
function varargout = select_spec_v28_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
function editspec_CreateFcn(hObject, ~, ~)
% hObject    handle to editspec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
try rmappdata(0,'ListDialogAppData__'); end
delete(hObject);



% #########################################################################################
% #################################### MANUAL INPUTS ######################################
% #########################################################################################

% ############# Injector
function diamval_Callback(hObject, eventdata, handles)
% hObject    handle to diamval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of diamval as text
%        str2double(get(hObject,'String')) returns contents of diamval as a double
data = getappdata(0,'ListDialogAppData__');
data.injector.diameter_m = 1e-6*str2num(get(handles.diamval,'String'));
setappdata(0,'ListDialogAppData__',data);
function caval_Callback(hObject, eventdata, handles)
% hObject    handle to caval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of caval as text
%        str2double(get(hObject,'String')) returns contents of caval as a double
data = getappdata(0,'ListDialogAppData__');
data.injector.ca = str2num(get(handles.caval,'String'));
setappdata(0,'ListDialogAppData__',data);
function cdval_Callback(hObject, eventdata, handles)
% hObject    handle to cdval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cdval as text
%        str2double(get(hObject,'String')) returns contents of cdval as a double
 data = getappdata(0,'ListDialogAppData__');
data.injector.cd = str2num(get(handles.cdval,'String'));
setappdata(0,'ListDialogAppData__',data);
function thetaval_Callback(hObject, eventdata, handles)
% hObject    handle to thetaval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thetaval as text
%        str2double(get(hObject,'String')) returns contents of thetaval as a double

data = getappdata(0,'ListDialogAppData__');
data.injector.theta = str2num(get(handles.thetaval,'String'));
setappdata(0,'ListDialogAppData__',data);
function rjetval_Callback(hObject, eventdata, handles)
% hObject    handle to rjetval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rjetval as text
%        str2double(get(hObject,'String')) returns contents of rjetval as a double

data = getappdata(0,'ListDialogAppData__');
data.injector.Rjet_m = str2num(get(handles.rjetval,'String'));
setappdata(0,'ListDialogAppData__',data);
function xjetval_Callback(hObject, eventdata, handles)
% hObject    handle to xjetval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xjetval as text
%        str2double(get(hObject,'String')) returns contents of xjetval as a double

data = getappdata(0,'ListDialogAppData__');
data.injector.xjet_m = str2num(get(handles.xjetval,'String'))/1000;
setappdata(0,'ListDialogAppData__',data);
function mfrval_Callback(hObject, eventdata, handles)
% hObject    handle to mfrval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mfrval as text
%        str2double(get(hObject,'String')) returns contents of mfrval as a double

data = getappdata(0,'ListDialogAppData__');
data.injector.mfr_kg_s = 1e-3*str2num(get(handles.mfrval,'String'));
setappdata(0,'ListDialogAppData__',data);
function momentumval_Callback(hObject, eventdata, handles)
% hObject    handle to momentumval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of momentumval as text
%        str2double(get(hObject,'String')) returns contents of momentumval as a double

data = getappdata(0,'ListDialogAppData__');
data.injector.momentum_N = str2num(get(handles.momentumval,'String'));
setappdata(0,'ListDialogAppData__',data);
function uipanel6_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel6
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% handles = store_data(handles);
% dsp_input_data(handles,handles.store);
data = getappdata(0,'ListDialogAppData__');
data.injector.casel = get(handles.casel,'Value');
if data.injector.casel
    set(handles.caval,'Enable','on');
    set(handles.cdval,'Enable','on');
    set(handles.mfrval,'Enable','off');
    set(handles.momentumval,'Enable','off');
else
    set(handles.caval,'Enable','off');
    set(handles.cdval,'Enable','off');
    set(handles.mfrval,'Enable','on');
    set(handles.momentumval,'Enable','on');
end
setappdata(0,'ListDialogAppData__',data);
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename,pathname] = uigetfile('*.xls','Select the ROI data');
data.injector.roi_file = fullfile(pathname,filename);
[num,txt,raw] = xlsread(data.injector.roi_file);
data.injector.roi_time_mus = num(:,1);
data.injector.roi_injvelnorm = num(:,2);
try data.injector.roi_injvel_m_s = num(:,3); end

% ############ Fuel
function tfuelval_Callback(hObject, eventdata, handles)
% hObject    handle to tfuelval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tfuelval as text
%        str2double(get(hObject,'String')) returns contents of tfuelval as a double
data = getappdata(0,'ListDialogAppData__');
data.fuel.T_K = str2num(get(handles.tfuelval,'String'));
setappdata(0,'ListDialogAppData__',data);
function rhofuelval_Callback(hObject, eventdata, handles)
% hObject    handle to rhofuelval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rhofuelval as text
%        str2double(get(hObject,'String')) returns contents of rhofuelval as a double

data = getappdata(0,'ListDialogAppData__');
data.fuel.rho_kg_m3 = str2num(get(handles.rhofuelval,'String'));
setappdata(0,'ListDialogAppData__',data);
function pfuelval_Callback(hObject, eventdata, handles)
% hObject    handle to pfuelval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pfuelval as text
%        str2double(get(hObject,'String')) returns contents of pfuelval as a double
data = getappdata(0,'ListDialogAppData__');
data.fuel.p_bar = str2num(get(handles.pfuelval,'String'));
setappdata(0,'ListDialogAppData__',data);
function stoichval_Callback(hObject, eventdata, handles)
% hObject    handle to stoichval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stoichval as text
%        str2double(get(hObject,'String')) returns contents of stoichval as a double
data = getappdata(0,'ListDialogAppData__');
data.fuel.pco_air = str2num(get(handles.stoichval,'String'));
setappdata(0,'ListDialogAppData__',data);
function autorhofuel_Callback(hObject, eventdata, handles)
% hObject    handle to autorhofuel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of autorhofuel
data = getappdata(0,'ListDialogAppData__');
data.fuel.autorho = get(handles.autorhofuel,'Value');
setappdata(0,'ListDialogAppData__',data);

% ############ Ambient
function tambval_Callback(hObject, eventdata, handles)
% hObject    handle to tambval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tambval as text
%        str2double(get(hObject,'String')) returns contents of tambval as a double
data = getappdata(0,'ListDialogAppData__');
data.ambient.T_K = str2num(get(handles.tambval,'String'));
setappdata(0,'ListDialogAppData__',data);
function rhoambval_Callback(hObject, eventdata, handles)
% hObject    handle to rhoambval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rhoambval as text
%        str2double(get(hObject,'String')) returns contents of rhoambval as a double

data = getappdata(0,'ListDialogAppData__');
data.ambient.rho_kg_m3 = str2num(get(handles.rhoambval,'String'));
setappdata(0,'ListDialogAppData__',data);
function pambval_Callback(hObject, eventdata, handles)
% hObject    handle to pambval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pambval as text
%        str2double(get(hObject,'String')) returns contents of pambval as a double
data = getappdata(0,'ListDialogAppData__');
data.ambient.p_bar = str2num(get(handles.pambval,'String'));
setappdata(0,'ListDialogAppData__',data);
function xo2val_Callback(hObject, eventdata, handles)
% hObject    handle to xo2val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xo2val as text
%        str2double(get(hObject,'String')) returns contents of xo2val as a double
data = getappdata(0,'ListDialogAppData__');
data.ambient.xo2 = str2num(get(handles.xo2val,'String'));
setappdata(0,'ListDialogAppData__',data);
function yo2val_Callback(hObject, eventdata, handles)
% hObject    handle to yo2val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yo2val as text
%        str2double(get(hObject,'String')) returns contents of yo2val as a double
data = getappdata(0,'ListDialogAppData__');
data.ambient.yo2 = str2num(get(handles.yo2val,'String'));
setappdata(0,'ListDialogAppData__',data);
function Mfg_Callback(hObject, eventdata, handles)
% hObject    handle to Mfg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Mfg as text
%        str2double(get(hObject,'String')) returns contents of Mfg as a double
data = getappdata(0,'ListDialogAppData__');
data.ambient.Mfg = str2num(get(handles.Mfg,'String'));
setappdata(0,'ListDialogAppData__',data);
function uipanel15_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel15 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
% --- Executes during object creation, after setting all properties.
data = getappdata(0,'ListDialogAppData__');
data.ambient.xo2sel = get(handles.xo2sel,'Value');
if data.ambient.xo2sel
    set(handles.xo2val,'Enable','on');
    set(handles.yo2val,'Enable','off');
else
    set(handles.xo2val,'Enable','off');
    set(handles.yo2val,'Enable','on');
end
setappdata(0,'ListDialogAppData__',data);
function uipanel17_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel17 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
data = getappdata(0,'ListDialogAppData__');
data.ambient.pambsel = get(handles.pambsel,'Value');
if data.ambient.pambsel
    set(handles.pambval,'Enable','on');
    set(handles.Mfg,'Enable','off');
else
    set(handles.pambval,'Enable','off');
    set(handles.Mfg,'Enable','on');
end
setappdata(0,'ListDialogAppData__',data);
function autostoich_Callback(hObject, eventdata, handles)
% hObject    handle to autostoich (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data = getappdata(0,'ListDialogAppData__');

% Calcul du PCO
prompt = {'xCxHyOzNt:','yCxHyOzNt:','zCxHyOzNt:','tCxHyOzNt:'};
dlg_title = 'Fuel Composition';
num_lines = 1;
def = {'','','',''};
answer = inputdlg(prompt,dlg_title,num_lines,def);

xCxHyOzNt = str2num(answer{1});
yCxHyOzNt = str2num(answer{2});
zCxHyOzNt = str2num(answer{3});
tCxHyOzNt = str2num(answer{4});

% Useful constants
pco_air = function_pco_air(xCxHyOzNt,yCxHyOzNt,zCxHyOzNt,tCxHyOzNt);
data.fuel.pco_air = pco_air;
setappdata(0,'ListDialogAppData__',data);
set(handles.stoichval,'String',num2str(data.fuel.pco_air,'%2.2f'));

% ############ Options 0D
function totallengthval_Callback(hObject, eventdata, handles)
% hObject    handle to totallengthval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of totallengthval as text
%        str2double(get(hObject,'String')) returns contents of totallengthval as a double
data = getappdata(0,'ListDialogAppData__');
data.options0d.total_length_m = 1e-3*str2num(get(handles.totallengthval,'String'));
setappdata(0,'ListDialogAppData__',data);
function dxval_Callback(hObject, eventdata, handles)
% hObject    handle to dxval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dxval as text
%        str2double(get(hObject,'String')) returns contents of dxval as a double
data = getappdata(0,'ListDialogAppData__');
data.options0d.dx_m = 1e-3*str2num(get(handles.dxval,'String'));
setappdata(0,'ListDialogAppData__',data);
function alphaman_Callback(hObject, eventdata, handles)
% hObject    handle to alphaman (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alphaman as text
%        str2double(get(hObject,'String')) returns contents of alphaman as a double
data = getappdata(0,'ListDialogAppData__');
data.options0d.alpha = str2num(get(handles.alphaman,'String'));
setappdata(0,'ListDialogAppData__',data);
function sieberscoeff_Callback(hObject, eventdata, handles)
% hObject    handle to sieberscoeff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sieberscoeff as text
%        str2double(get(hObject,'String')) returns contents of sieberscoeff as a double
data = getappdata(0,'ListDialogAppData__');
data.options0d.a_Siebers = str2num(get(handles.sieberscoeff,'String'));
setappdata(0,'ListDialogAppData__',data);
function uipanel7_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel7
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

data = getappdata(0,'ListDialogAppData__');
data.options0d.unif_profile = get(handles.unifprof,'Value');
if data.options0d.unif_profile
    set(handles.alphaman,'Enable','off');
    set(handles.sieberscoeff,'Enable','on');
    set(handles.alphaman,'String',num2str(1e+10,'%1.3f'));
else
    set(handles.alphaman,'Enable','on');
    set(handles.sieberscoeff,'Enable','off');
    try set(handles.alphaman,'String',num2str(data.options0d.alpha,'%1.3f')); end
end
setappdata(0,'ListDialogAppData__',data);
function xoffset_Callback(hObject, eventdata, handles)
% hObject    handle to xoffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xoffset as text
%        str2double(get(hObject,'String')) returns contents of xoffset as a double
data = getappdata(0,'ListDialogAppData__');
data.options0d.xoffset_m = 1e-3*str2num(get(handles.xoffset,'String'));
setappdata(0,'ListDialogAppData__',data);
function toffset_Callback(hObject, eventdata, handles)
% hObject    handle to toffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of toffset as text
%        str2double(get(hObject,'String')) returns contents of toffset as a double
data = getappdata(0,'ListDialogAppData__');
data.options0d.toffset_s = 1e-3*str2num(get(handles.toffset,'String'));
setappdata(0,'ListDialogAppData__',data);

% ########### Lift-off Length
function coeflolval_Callback(hObject, eventdata, handles)
% hObject    handle to coeflolval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of coeflolval as text
%        str2double(get(hObject,'String')) returns contents of coeflolval as a double
data = getappdata(0,'ListDialogAppData__');
data.lol.coeflol = str2num(get(handles.coeflolval,'String'));
setappdata(0,'ListDialogAppData__',data);
function autolol_Callback(hObject, eventdata, handles)
% hObject    handle to autolol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of autolol
data = getappdata(0,'ListDialogAppData__');
data.lol.autolol = get(handles.autolol,'Value');
setappdata(0,'ListDialogAppData__',data);
function estimlolval_Callback(hObject, eventdata, handles)
% hObject    handle to estimlolval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of estimlolval as text
%        str2double(get(hObject,'String')) returns contents of estimlolval as a double
data = getappdata(0,'ListDialogAppData__');
data.lol.lolxp = str2num(get(handles.estimlolval,'String'));
setappdata(0,'ListDialogAppData__',data);

% ########### Vapor penetration
function penvapbutton_Callback(hObject, eventdata, handles)
% hObject    handle to penvapbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = getappdata(0,'ListDialogAppData__');
[filename,pathname] = uigetfile('*.xls','Select the vapor penetration file');
data.penvap_xp.folder = fullfile(pathname,filename);
[num,txt,raw] = xlsread(data.penvap_xp.folder);
t_mus = num(:,1);
penvap_mm = num(:,2);
try std_penvap_mm = num(:,3); end
try
    data.penvap_xp.time_mus = t_mus;
    data.penvap_xp.penvap_mm = penvap_mm;
    data.penvap_xp.std_penvap_mm = std_penvap_mm;
end
setappdata(0,'ListDialogAppData__',data);
function adjsoi_Callback(hObject, eventdata, handles)
% hObject    handle to adjsoi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of adjsoi
data = getappdata(0,'ListDialogAppData__');
data.mixing.adjust_soi = get(handles.adjsoi,'Value');
setappdata(0,'ListDialogAppData__',data);
function plotapprox_Callback(hObject, eventdata, handles)
% hObject    handle to plotapprox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotapprox
data = getappdata(0,'ListDialogAppData__');
data.penvap_xp.plotapprox = get(handles.plotapprox,'Value');
setappdata(0,'ListDialogAppData__',data);
function t1fit_Callback(hObject, eventdata, handles)
% hObject    handle to t1fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t1fit as text
%        str2double(get(hObject,'String')) returns contents of t1fit as a double
data = getappdata(0,'ListDialogAppData__');
data.penvap_xp.t1fit = str2num(get(handles.t1fit,'String'));
[a,i1fit] = min(abs(data.penvap_xp.time_mus/1000-data.penvap_xp.t1fit));
data.penvap_xp.i1fit = i1fit;
setappdata(0,'ListDialogAppData__',data);
function t2fit_Callback(hObject, eventdata, handles)
% hObject    handle to t2fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t2fit as text
%        str2double(get(hObject,'String')) returns contents of t2fit as a double
data = getappdata(0,'ListDialogAppData__');
data.penvap_xp.t2fit = str2num(get(handles.t2fit,'String'));
[a,i2fit] = min(abs(data.penvap_xp.time_mus/1000-data.penvap_xp.t2fit));
data.penvap_xp.i2fit = i2fit;
setappdata(0,'ListDialogAppData__',data);
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data = getappdata(0,'ListDialogAppData__');

if data.mixing.adjust_soi
    [theta,soi_corr] = error_function(data);
else
    myfun = @(x) error_function_1param(data,x);
    theta = fminsearch(myfun,20,optimset('TolX',1,'TolFun',1));
    soi_corr = 0;
end

data.mixing.theta = theta;
data.mixing.soi_corr = soi_corr;

set(handles.anglemuscval,'String',num2str(data.mixing.theta,'%2.1f'));
set(handles.soimuscval,'String',num2str(data.mixing.soi_corr,'%2.1f'));

setappdata(0,'ListDialogAppData__',data);



% #########################################################################################
% #################################### COMPUTATIONS #######################################
% #########################################################################################

% Injector
function compute_data_injector(handles)

data = getappdata(0,'ListDialogAppData__');
if data.injector.casel %Ca et Cd connus
    try ca = data.injector.ca; end
    try cd = data.injector.cd; end
    try pamb = data.ambient.p_bar; end
    try pfuel = data.fuel.p_bar; end
    try rhofuel = data.fuel.rho_kg_m3; end
    try d0 = data.injector.diameter_m; end    
    
    try
        u0 = cd/ca*sqrt(2*(pfuel-pamb)*1e5/rhofuel);        
        data.injector.u0_m_s = u0;
        set(handles.u0val,'String',num2str(u0,'%4.1f'));
        setappdata(0,'ListDialogAppData__',data);
    catch
        set(handles.u0val,'String','!');
    end

    try
        deff = sqrt(ca) * d0;
        data.injector.eff_diameter_m = deff;
        set(handles.effdiamval,'String',num2str(deff*1e6,'%2.1f'));        
        setappdata(0,'ListDialogAppData__',data);
    catch
        set(handles.effdiamval,'String','!');
    end
  
    try
        cv = cd/ca;
        data.injector.cv = cv;
        set(handles.cvval,'String',num2str(cv,'%1.2f'));        
        setappdata(0,'ListDialogAppData__',data);
    catch
        set(handles.cvval,'String','!');
    end

    try
        mfr = cd*sqrt(2*rhofuel*(pfuel-pamb)*1e5)*pi*(d0^2)/4 ;
        data.injector.mfr_kg_s = mfr;
        set(handles.mfrval,'String',num2str(mfr*1000,'%2.2f'));
        setappdata(0,'ListDialogAppData__',data);
    catch
        set(handles.mfrval,'String','!');
    end
    
    try
        M = mfr * cd/ca*sqrt(2*(pfuel-pamb)*1e5/rhofuel);
        data.injector.momentum_N = M;
        set(handles.momentumval,'String',num2str(M,'%2.2f'));
        setappdata(0,'ListDialogAppData__',data);
    catch
        set(handles.momentumval,'String','!');
    end
else
        try mfr = data.injector.mfr_kg_s; end
        try M = data.injector.momentum_N; end
        try pamb = data.ambient.p_bar; end
        try rhofuel = data.fuel.rho_kg_m3; end        
        try pfuel = data.fuel.p_bar; end
        try d0 = data.injector.diameter_m; end    
    
        try
            deff = sqrt(4*mfr^2/pi/rhofuel/M);
            data.injector.eff_diameter_m = deff;
            set(handles.effdiamval,'String',num2str(deff*1e6,'%2.1f'));
            setappdata(0,'ListDialogAppData__',data);
        catch
            set(handles.effdiamval,'String','!');
        end
    
        try
            ca = deff^2/d0^2;
            data.injector.ca = ca;
            set(handles.caval,'String',num2str(ca,'%1.2f'));            
            setappdata(0,'ListDialogAppData__',data);
        catch
            set(handles.caval,'String','!');
        end
    
        try
            u0 = M/mfr;
            data.injector.u0_m_s = u0;
            set(handles.u0val,'String',num2str(u0,'%4.1f'));
            setappdata(0,'ListDialogAppData__',data);
        catch
            set(handles.u0val,'String','!');
        end
    
        try
            uth = sqrt(2*(pfuel-pamb)*1e5/rhofuel);
            cd = ca * u0/uth;
            data.injector.cd = cd;
            set(handles.cdval,'String',num2str(cd,'%1.2f'));
            setappdata(0,'ListDialogAppData__',data);
        catch
            set(handles.cdval,'String','!');
        end
end

% Calcul de l'angle de siebers
if data.options0d.unif_profile
    try        
        theta = data.injector.theta;
        a_Siebers = data.options0d.a_Siebers;
        thetacorr=2*atand(tand(theta/2)*a_Siebers);
        data.injector.thetacorr = thetacorr;
        set(handles.thetacorrval,'String',num2str(thetacorr,'%2.2f'));
        setappdata(0,'ListDialogAppData__',data);
    catch
        set(handles.thetacorrval,'String','!');
    end
end

% Fuel
function compute_data_fuel(handles)

data = getappdata(0,'ListDialogAppData__');

if data.fuel.autorho
    try 
        set(handles.rhofuelval,'Enable','off');
        [Rho_f, Mu_f] = C12Characteristics(data.ambient.p_bar/10, data.fuel.T_K); 
        data.fuel.rho_kg_m3 = Rho_f;
        set(handles.rhofuelval,'String',num2str(Rho_f,'%4.1f'));
        setappdata(0,'ListDialogAppData__',data);
    catch
        set(handles.rhofuelval,'String','!');
    end
else
    set(handles.rhofuelval,'Enable','on');
end

pco_air = data.fuel.pco_air;
Mfg = data.ambient.Mfg;
[pco_fg,Zst] = function_pco_fg(pco_air,data.ambient.xo2/100);
data.fuel.pco_fg = pco_fg;
data.fuel.Zst = Zst;

set(handles.stoichcorrval,'String',num2str(data.fuel.pco_fg,'%2.2f'));
set(handles.zstval,'String',num2str(data.fuel.Zst,'%2.4f'));   

setappdata(0,'ListDialogAppData__',data);
     
% Ambient
function compute_data_ambient(handles)

data = getappdata(0,'ListDialogAppData__');

tamb = data.ambient.T_K;
rhoamb = data.ambient.rho_kg_m3;
if data.ambient.pambsel
    pamb = data.ambient.p_bar;
    Mfg = 8.314/(pamb*1e5/rhoamb/tamb)*1000; 
    data.ambient.Mfg = Mfg;
    set(handles.Mfg,'String',num2str(data.ambient.Mfg,'%2.3f'));   
else
    Mfg = data.ambient.Mfg;
    pamb = rhoamb*8.314/Mfg*tamb/1e2;
    data.ambient.p_bar = pamb;
    set(handles.pambval,'String',num2str(data.ambient.p_bar,'%2.2f'));  
end
setappdata(0,'ListDialogAppData__',data);

Mfg = data.ambient.Mfg;
if data.ambient.xo2sel
    try
        yo2 = function_xo2_2_yo2(data.ambient.xo2,Mfg);
        data.ambient.yo2 = yo2;
        set(handles.yo2val,'String',num2str(data.ambient.yo2,'%2.2f'));
    catch
        set(handles.yo2val,'String','!');
    end
else
    try
        xo2 = function_yo2_2_xo2(data.ambient.yo2,Mfg);
        data.ambient.xo2 = xo2;
        set(handles.xo2val,'String',num2str(data.ambient.xo2,'%2.2f'));
    catch
        set(handles.yo2val,'String','!');
    end
end
setappdata(0,'ListDialogAppData__',data);

% Lift-off length
function compute_data_lol(handles)

data = getappdata(0,'ListDialogAppData__');

if data.lol.autolol
    set(handles.estimlolval,'Enable','off');
    try
        coef = str2num(get(handles.coeflolval,'String'));
        u0 = data.injector.u0_m_s;
        tamb = data.ambient.T_K;
        rhoamb = data.ambient.rho_kg_m3;
        d0 = data.injector.diameter_m*1e6;
        zst = data.fuel.Zst;
        lol = coef*tamb^-3.74*rhoamb^-0.85*d0^0.34*u0/zst;
        data.lol.lolestim = lol;
        set(handles.estimlolval,'String',num2str(lol,'%2.2f'));
    catch
        set(handles.estimlolval,'String','!');
    end
else
    set(handles.estimlolval,'Enable','on');
    lol = data.lol.lolxp;
    set(handles.estimlolval,'String',num2str(lol,'%2.2f'));
end
setappdata(0,'ListDialogAppData__',data);

% Mixing field
function compute_data_mixing(handles)

data = getappdata(0,'ListDialogAppData__');

diameter    = data.injector.eff_diameter_m;
xjet        = data.injector.xjet_m;
total_length= data.options0d.total_length_m;
dx          = data.options0d.dx_m;
alpha       = data.options0d.alpha;
u0          = data.injector.u0_m_s;
m0          = data.injector.mfr_kg_s;
M0          = data.injector.momentum_N;
rhoa        = data.ambient.rho_kg_m3;
rhof        = data.fuel.rho_kg_m3;
pco_fg      = data.fuel.pco_fg;
theta       = data.injector.theta;
xoffset     = data.options0d.xoffset_m;
toffset     = data.options0d.toffset_s;

Rjet = cumsum([xjet(1) diff(xjet)].*tand(theta/2));
diameter = data.injector.eff_diameter_m;% + 2*tand(theta/2)*xoffset*0;


xx_m = [dx/2:dx:total_length];        % coordinates of CV centers, relative to nozzle exit
ivardens = 0;
[velocity,mole_fraction,mass_fraction,FA,EntRate,rhoac,Tadiab,t_s,alpha_vect,r_m,R,area_m2] = function_mixing_fields(xx_m,theta,xjet,Rjet,diameter,m0,M0,rhof,rhoa,alpha,ivardens);

t0 = interp1([0 xx_m],t_s,xoffset);
t0 = toffset - t0;

phi.mean = FA.mean*pco_fg; % Steady jet equivalence ratio
phi.axial = FA.axial*pco_fg;
phi.field = FA.field *pco_fg;

data.mixing.x_m = xx_m;
data.mixing.r_m = r_m;
data.mixing.jet_radius = R;
data.mixing.section_area = area_m2;
data.mixing.t_s = t_s+t0;
data.mixing.alpha_vect = alpha_vect;
data.mixing.phi = phi;
data.mixing.mass_fraction = mass_fraction;
data.mixing.mole_fraction = mole_fraction;
data.mixing.FA = FA;
data.mixing.velocity = velocity;
data.mixing.rhoac = rhoac;
data.mixing.Tadiab = Tadiab;
data.mixing.EntRate=EntRate;
% data.mixing.AirMFR=AirMFR;

setappdata(0,'ListDialogAppData__',data);

% Equivalence ratio at lift-off length
function compute_data_phi(handles)

data = getappdata(0,'ListDialogAppData__');

if data.lol.autolol
    lol = data.lol.lolestim;
else
    lol = data.lol.lolxp;
end

philol_axial = interp1(data.mixing.x_m*1e3,data.mixing.phi.axial,lol);
data.mixing.philol_axial = philol_axial;
set(handles.phimusculusaxial,'String',num2str(philol_axial,'%2.2f'));


philol_mean = interp1(data.mixing.x_m*1e3,data.mixing.phi.mean,lol);
data.mixing.philol_mean = philol_mean;
set(handles.phimusculusmean,'String',num2str(philol_mean,'%2.2f'));
setappdata(0,'ListDialogAppData__',data);


% #########################################################################################
% ############################### PERSONAL FUNCTIONS ######################################
% #########################################################################################
function [cp_func,cv_func,h_func,r,mu_func,lambda_func] = Lecture_Prop_Gaz_Single(dirprop,gaz_file)

gaz_file_list = {gaz_file};

cp=[];
cv=[];
r=[];
mu=[];
lambda=[];

cp_func=[];
cv_func=[];

for ii=1:length(gaz_file_list)
    
    gaz_file = gaz_file_list{ii};
    
    %% Lecture du fichier
    DataBrutes = textread(fullfile(dirprop,gaz_file),'%s');
    DataTraitees = [];
    DataName=[];
    count = 0;
    for jj=1:length(DataBrutes)
        if isempty(str2num(DataBrutes{jj}))==0
            %         if length(str2num(DataBrutes{jj})>0)
            count = count +1;
            DataTraitees{count,1} = str2num(DataBrutes{jj});
            DataName=[];
        else
            DataName=[DataName ' ' DataBrutes{jj}];
            DataTraitees{count,2} = DataName;
        end
    end
    
    %% Janaf
    if length(DataTraitees)>7 % Janaf
        jj=1;
        while jj<length(DataTraitees)
            str=DataTraitees{jj,2};
            i_up = findstr(str,'UPPER');
            if isempty(i_up)==0
                i_up = jj;
                break
            end
            jj = jj+1;
        end
        
        while jj<length(DataTraitees)
            str=DataTraitees{jj,2};
            i_dw = findstr(str,'LOWER');
            if isempty(i_dw)==0
                i_dw = jj;
                break
            end
            jj = jj+1;
        end
        
        coef_up=[];
        for jj=i_up:(i_up+5)
            coef_up = [coef_up;DataTraitees{jj,1}];
        end
        
        coef_dw=[];
        for jj=i_dw:(i_dw+5)
            coef_dw = [coef_dw;DataTraitees{jj,1}];
        end
        
        %         T_inf = DataTraitees{3,1};
        %         T_sup = DataTraitees{4,1};
        T_mid = DataTraitees{5,1};
        
        r(ii) = 8.314 / DataTraitees{2,1} * 1e3;
        cp_func = @(x) r(ii)*dot([x'.^0 x'.^1 x'.^2 x'.^3 x'.^4]',[coef_dw(1:end-1)*(x<=T_mid)+coef_up(1:end-1)*(x>T_mid)]);
        cv_func = @(x) cp_func{ii}(x) - r(ii);
        h_func = @(x) r(ii)*dot([x'.^1 1/2*x'.^2 1/3*x'.^3 1/4*x'.^4 1/5*x'.^5 x'.^0]',[coef_dw*(x<=T_mid)+coef_up*(x>T_mid)]);
        
        mu_func = @(x)  1e-7*dot([x'.^0 x'.^1 x'.^2]',[DataTraitees{end-5,1};DataTraitees{end-4,1};DataTraitees{end-3,1}]);
        
        lambda_func = @(x)  dot([x'.^0 x'.^1 x'.^2]',[DataTraitees{end-2,1};DataTraitees{end-1,1};DataTraitees{end,1}]);
        
    else %% cte
        r = DataTraitees{5,1};
        cp_func = @(x) 0*x + DataTraitees{2,1};
        cv_func = @(x) cp_func{ii}(x) - r(ii);
        h_func = @(x) cp_func{ii}(x)*x;
        mu_func = @(x)  0*x + DataTraitees{3,1};
        lambda_func = @(x)  0*x + DataTraitees{4,1};
    end
end
function pco_air = function_pco_air(x,y,z,t)

% 09/01/2014 - L.M. MALBEC
% This function computes the pco (mair/mfuel at stoechiometry) knowing the
% composition of the fuel
%
% Inputs
%   - x: number of C moles in fuel
%   - y: number of H moles in fuel
%   - z: number of O moles in fuel
%   - t: number of N moles in fuel


% Notation:
% Rst: stoechiometric ratio: number of mole of air necessary to burn
% completely one mole of fuel
% Rst = x+y/4-z/2
% alpha: molar fraction of n2 to o2 in air (alpha = 79/21 = 3.7619)
%
% Computation:
% PCO = mair/mfuel = nair*Mair/(nfuel*Mfuel) = Rst*(Mo2+ alpha*Mn2)/Mfuel

M_C   = 12.0;       %% molar weight of carbon
M_H   = 1.008;      %% molar weight of hydrogen
M_O  = 16.0;        %% molar weight of oxygen
M_N   = 14.008;     %% molar weight of nitrogene
alpha = 79/21;

Rst = x+y/4-z/2;
pco_air = Rst*(M_O*2+alpha*M_N*2)/(x*M_C+y*M_H+z*M_O+t*M_N);
function [pco_fg,zst] = function_pco_fg(pco_air,xo2)

% 09/01/2014 - L.M. MALBEC
% This function computes the pco in fresh gasses (mgf/mfuel at stoechiometry) knowing the
% pco in air and the o2 mass fraction in fresh gasses
%
% Inputs
%   - pco_air
%   - xo2: mass fraction of O2 in fresh gasses

% Outputs:
%   pco_fg = mfuel/m_fg
%   zst = mfuel/(mfuel+m_fg);

% Notation:
% pco_air = mair/mfuel = mo2/xo2air/mfuel
% pco_gf = mgf/mfuel
% mgf = mo2/xo2
% => pco_fg = pco_air * xo2air/xo2;

xo2air = 0.2330;
pco_fg = pco_air * xo2air/xo2;
zst = 1/(1+pco_fg);
function yo2 = function_xo2_2_yo2(xo2,Mmixt)

% 16/09/2014 - L.M. MALBEC
% This function converts mass fraction (xo2) into  mole fraction (yo2) for O2

% Inputs
%   - xo2: mass fraction of o2 in the mixture
%   - Mmixt: Molar mass of the mixture

yo2 = xo2*Mmixt/32;
function xo2 = function_yo2_2_xo2(yo2,Mmixt)

% 16/09/2014 - L.M. MALBEC
% This function converts mole fraction (yo2) into mass fraction (xo2) for O2

% Inputs
%   - yo2: mole fraction of o2 in the mixture
%   - Mmixt: Molar mass of the mixture

xo2 = yo2*32/Mmixt;
function plot_velocities(handles)

data = getappdata(0,'ListDialogAppData__');

%% Axial velocity
delete(get(handles.axes2, 'Children'));
axes(handles.axes2);
hold on
plot(data.mixing.x_m*1e3,data.mixing.velocity.tip,'r','linewidth',2)
plot(data.mixing.x_m*1e3,data.mixing.velocity.axial,':r','linewidth',2)
% plot(data.mixing.x_m*1e3,data.mixing.velocity.tip,'--r','linewidth',2)
set(gca,'Box','on','FontName','Eras Light ITC','linewidth',2,'fontsize',10,'fontweight','bold','color',[1 1 1]*1)
grid on
xlabel('Axial dist. [mm]')
ylabel('Axial Velocity [m/s]')
hold off

% Velocity field
delete(get(handles.axes8, 'Children'));
axes(handles.axes8);
imagesc(data.mixing.x_m*1e3,data.mixing.r_m*1e3,data.mixing.velocity.field'),
set(gca,'Box','on','FontName','Eras Light ITC','linewidth',2,'fontsize',10,'fontweight','bold','color',[1 1 1]*1)
grid on
caxis([0 max(data.mixing.velocity.axial)/4])

%% Axial afr
delete(get(handles.axes3, 'Children'));
axes(handles.axes3);
hold on
plot(data.mixing.x_m*1e3,data.mixing.phi.mean,'r','linewidth',2)
plot(data.mixing.x_m*1e3,data.mixing.phi.axial,':r','linewidth',2)
set(gca,'Box','on','FontName','Eras Light ITC','linewidth',2,'fontsize',10,'fontweight','bold','color',[1 1 1]*1)
grid on
ylim = get(gca,'Ylim');
% xlim = get(gca,'Xlim');
ylim(2) = min(ylim(2),10);
set(gca,'Ylim',[0 10]);
set(gca,'Xlim',[0 round(max(data.mixing.x_m*1e3))],'Ylim',ylim)
xlabel('Axial dist. [mm]')
ylabel('AFR [-]')
hold off

% AFR field
delete(get(handles.axes9, 'Children'));
axes(handles.axes9);
imagesc(data.mixing.x_m*1e3,data.mixing.r_m*1e3,data.mixing.phi.field'),
set(gca,'ydir','reverse')
set(gca,'Box','on','FontName','Eras Light ITC','linewidth',2,'fontsize',10,'fontweight','bold','color',[1 1 1]*1)
grid on
caxis([0 5])


%% Axial fuel mass fraction
delete(get(handles.axes4, 'Children'));
axes(handles.axes4);
hold on
plot(data.mixing.x_m*1e3,data.mixing.mass_fraction.mean,'r','linewidth',2)
plot(data.mixing.x_m*1e3,data.mixing.mass_fraction.axial,':r','linewidth',2)
set(gca,'Box','on','FontName','Eras Light ITC','linewidth',2,'fontsize',10,'fontweight','bold','color',[1 1 1]*1)
grid on
% ylim = get(gca,'Ylim');
% xlim = get(gca,'Xlim');
% ylim(2) = min(ylim(2),10/data.fuel.pco_fg);
set(gca,'Ylim',[0 0.2])
set(gca,'Xlim',[0 round(max(data.mixing.x_m*1e3))])
xlabel('Axial dist. [mm]')
ylabel('Fuel mass fraction [-]')
% set(gca,'ydir','reverse')
hold off

% Fuel mass fraction field
delete(get(handles.axes11, 'Children'));
axes(handles.axes11);
imagesc(data.mixing.x_m*1e3,data.mixing.r_m*1e3,data.mixing.mass_fraction.field'),
set(gca,'ydir','reverse')
set(gca,'Box','on','FontName','Eras Light ITC','linewidth',2,'fontsize',10,'fontweight','bold','color',[1 1 1]*1)
grid on
caxis([0 0.2])

%% Adiabatic mixing temperature
delete(get(handles.axes6, 'Children'));
axes(handles.axes6);
hold on
plot(data.mixing.x_m*1e3,data.mixing.Tadiab.mean,'r','linewidth',2)
plot(data.mixing.x_m*1e3,data.mixing.Tadiab.axial,':r','linewidth',2)
set(gca,'Box','on','FontName','Eras Light ITC','linewidth',2,'fontsize',10,'fontweight','bold','color',[1 1 1]*1)
grid on
ylim = get(gca,'Xlim');
% xlim = get(gca,'Ylim');
% xlim(2) = min(xlim(2),1000);
set(gca,'Ylim',[300 1000]);
set(gca,'Xlim',[0 round(max(data.mixing.x_m*1e3))],'Xlim',ylim)
xlabel('Axial dist. [mm]')
ylabel('Adiab. mix. temp. [K]')
% set(gca,'ydir','reverse')
hold off

% Temperature field
delete(get(handles.axes10, 'Children'));
axes(handles.axes10);
imagesc(data.mixing.x_m*1e3,data.mixing.r_m*1e3,data.mixing.Tadiab.field'),
caxis([max(data.mixing.Tadiab.field(:))-200 max(data.mixing.Tadiab.field(:))])
set(gca,'ydir','reverse')
set(gca,'Box','on','FontName','Eras Light ITC','linewidth',2,'fontsize',10,'fontweight','bold','color',[1 1 1]*1)
grid on


%% Air Entrainment
delete(get(handles.axes7, 'Children'));
axes(handles.axes7);
hold on
% plot(data.mixing.x_m*1e3,data.mixing.AirMFR.raw*1e3,'r','linewidth',2)
plot(data.mixing.x_m*1e3,data.mixing.EntRate.raw,'r--','linewidth',2)
plot(data.mixing.x_m*1e3,data.mixing.EntRate.norm,'r:','linewidth',2)
set(gca,'Box','on','FontName','Eras Light ITC','linewidth',2,'fontsize',10,'fontweight','bold','color',[1 1 1]*1)
grid on
xlabel('Distance [mm]')
ylabel('Air MFR [g/s]')
hold off

%% Spray penetration
delete(get(handles.axes5, 'Children'));
axes(handles.axes5);
hold on
try
    plot(data.penvap_xp.time_mus/1000,data.penvap_xp.penvap_mm,'ok','linewidth',2);
    plot(data.penvap_xp.time_mus/1000,data.penvap_xp.penvap_mm+data.penvap_xp.std_penvap_mm,':k','linewidth',1);
    plot(data.penvap_xp.time_mus/1000,data.penvap_xp.penvap_mm-data.penvap_xp.std_penvap_mm,':k','linewidth',1);
end

try
     plot((data.penvap_xp.time_mus+data.mixing.soi_corr)/1000,data.penvap_xp.penvap_mm,'+k','linewidth',1);
end

try
    i1fit = data.penvap_xp.i1fit;
    i2fit = data.penvap_xp.i2fit;
    plot(data.penvap_xp.time_mus(i1fit:i2fit)/1000,data.penvap_xp.penvap_mm(i1fit:i2fit),'o','linewidth',2,'color',[1 1 1]*0.7);
end

try
    i1fit = data.penvap_xp.i1fit;
    i2fit = data.penvap_xp.i2fit;
    plot((data.penvap_xp.time_mus(i1fit:i2fit)-data.mixing.soi_corr)/1000,data.penvap_xp.penvap_mm(i1fit:i2fit),'+','linewidth',2,'color',[1 1 1]*0.7);
end
plot(data.mixing.t_s*1e3,[0 data.mixing.x_m]*1e3,'r','linewidth',2)
set(gca,'Box','on','FontName','Eras Light ITC','linewidth',2,'fontsize',10,'fontweight','bold','color',[1 1 1]*1)
grid on
xlim = get(gca,'Xlim');
ylim = get(gca,'Ylim');
set(gca,'Xlim',[0 100],'Ylim',[0 6]);
set(gca,'Xlim',xlim,'Ylim',ylim)
xlabel('Time [ms]')
ylabel('Distance [mm]')
ht=title(data.penvap_xp.folder);
set(ht,'fontSize',6,'fontname','eras light ITC');
hold off
function [theta,soi_corr] = error_function(data)

diameter    = data.injector.eff_diameter_m;
total_length= data.options0d.total_length_m;
dx          = data.options0d.dx_m;
alpha       = data.options0d.alpha;
m0          = data.injector.mfr_kg_s;
M0          = data.injector.momentum_N;
rhoa       = data.ambient.rho_kg_m3;
rhof       = data.fuel.rho_kg_m3;
theta       = data.injector.theta;

xx_m = (dx/2:dx:total_length);        % coordinates of CV centers, relative to nozzle exit
ivardens = 0;

% Computation of mean slope of experimental points
dpvxp = diff(data.penvap_xp.penvap_mm(i1fit:i2fit))./diff(data.penvap_xp.time_mus(i1fit:i2fit));
dpvxp_m = mean(dpvxp)*100;

% Computation of mean slope of 0d  points
time0d = interp1([0 x_m]*1e3,t*1e6,data.penvap_xp.penvap_mm(i1fit:i2fit));
time0d(isnan(time0d))=0;
isel = zeros(1,length(time0d));
for ii=1:length(time0d)
    [~,isel(ii)] = min(abs(time0d(ii)-t*1e6));
end
dpv0d = diff(x_m(isel-1)*1e3)./diff(t(isel)*1e6);
dpv0d_m = mean(dpv0d)*100;

% Theta is modified until the xp and 0d slopes are close enough
while abs(dpv0d_m-dpvxp_m)>0.1
    theta = theta + (dpv0d_m-dpvxp_m)*5
    
    [~,~,~,~,~,~,~,t_s,~,~,~,~] = function_mixing_fields(xx_m,theta,diameter,m0,M0,rhof,rhoa,alpha,ivardens);
    
    time0d = interp1([0 xx_m]*1e3,t_s*1e6,data.penvap_xp.penvap_mm(i1fit:i2fit));
    time0d(isnan(time0d))=0;
    isel = zeros(1,length(time0d));
    for ii=1:length(time0d)
        [~,isel(ii)] = min(abs(time0d(ii)-t*1e6));
    end
    dpv0d = diff(x_m(isel-1)*1e3)./diff(t(isel)*1e6);
    dpv0d_m = mean(dpv0d)*100;
end

% Computation of mean soi correction
soi_corr = mean(time0d) - mean(data.penvap_xp.time_mus(i1fit:i2fit));
function err = error_function_1param(data,theta2)


diameter    = data.injector.eff_diameter_m;
xjet        = data.injector.xjet_m;
theta       = data.injector.theta;
total_length= data.options0d.total_length_m;
dx          = data.options0d.dx_m;
alpha       = data.options0d.alpha;
u0          = data.injector.u0_m_s;
m0          = data.injector.mfr_kg_s;
M0          = data.injector.momentum_N;
rhoa       = data.ambient.rho_kg_m3;
rhof       = data.fuel.rho_kg_m3;
pco_fg      = data.fuel.pco_fg;
xoffset     = data.options0d.xoffset_m;
toffset     = data.options0d.toffset_s;

theta(end) = theta2;

% Virtual origine of the spray
Rjet = cumsum([xjet(1) diff(xjet)].*tand(theta/2));
diameter = diameter;% + 2*tand(theta/2)*xoffset;

xx_m = [dx/2:dx:total_length];        % coordinates of CV centers, relative to nozzle exit
ivardens = 0;
[velocity,mole_fraction,mass_fraction,FA,EntRate,rhoac,Tadiab,t_s,alpha_vect,r_m,R,area_m2] =  function_mixing_fields(xx_m,theta,xjet,Rjet,diameter,m0,M0,rhof,rhoa,alpha,ivardens);
xx_m = xx_m + xoffset;
t_s  = t_s  + toffset;

i1fit = data.penvap_xp.i1fit;
i2fit = data.penvap_xp.i2fit;
penvapxp = data.penvap_xp.penvap_mm(i1fit:i2fit);

penvap0d = interp1(t_s*1e6,[xoffset xx_m],data.penvap_xp.time_mus(i1fit:i2fit))*1e3;
penvap0d(isnan(penvap0d))=0;

theta
err = sum((penvap0d-penvapxp).^2)
function [Rho_f, Mu_f] = C12Characteristics(p, t)
%
% C12Characteristics(p, t)
%
% This function calculates the density and viscosity of the dodecane as a
% function of temperature and pressure. The data have been extracted from
% Caudwell D.R. [2003] and adjusted to the correct pressure with a linear
% interpolation.
%
% A two dimensional interpolation is used to calculate the characteristics
% of the fuel at specific conditions in the range 0.1 to 200 MPa for the
% pressure and 298.15 to 473.15 K for the temperature.
%
% "p" and "t" are respectively the pressure [MPa] and temperature [K] at
% which the density and viscosity have to be calculated.

% Pressure and temperature of the tests (ajusted)
P = [0.1 40 80 120 160 200];
T = [298.15 323.15 348.15 373.15 398.15 423.15 448.15 473.15];

% Experimental density measured by Caudwell (adjusted)
D(1,:) = [746.0 771.5 790.8 806.1 819.8 832.1];
D(2,:) = [727.4 756.8 777.3 794.5 808.9 821.3];
D(3,:) = [708.9 741.1 764.4 782.3 797.7 810.2];
D(4,:) = [689.8 726.7 750.9 770.9 787.7 800.7];
D(5,:) = [670.8 712.1 738.7 759.9 777.2 789.8];
D(6,:) = [651.1 698.1 727.3 749.4 768.2 780.1];
D(7,:) = [630.3 683.9 714.8 738.7 758.3 770.4];
D(8,:) = [608.9 670.2 704.3 729.0 748.6 761.2];
% Experimental viscosity measured by Caudwell (adjusted)
nu(1,:) = [1.344 2.139 3.154 4.497 6.128 7.882];
nu(2,:) = [0.911 1.422 2.012 2.762 3.627 4.642];
nu(3,:) = [0.659 1.006 1.428 1.886 2.452 3.124];
nu(4,:) = [0.503 0.777 1.069 1.406 1.809 2.269];
nu(5,:) = [0.401 0.623 0.855 1.123 1.395 1.728];
nu(6,:) = [0.324 0.512 0.706 0.915 1.158 1.454];
nu(7,:) = [0.264 0.429 0.588 0.763 0.966 1.251];
nu(8,:) = [0.218 0.367 0.510 0.657 0.813 1.047];

% % Density vs pressure is plotted for different temperatures
% figure(1), plot(P, D(1,:)), hold all;
% plot(P, D(2,:)), plot(P, D(3,:)), plot(P, D(4,:)), plot(P, D(5,:));
% plot(P, D(6,:)), plot(P, D(7,:)), plot(P, D(8,:));
% title('Density [kg/m^3]')
% legend('298.15 K', '323.15 K', '348.15 K', '373.15 K', '398.15 K', '423.15 K',...
%     '448.15 K', '473.15 K', 'location', 'southeast')
% % Viscosity vs pressure is plotted for different temperatures
% figure(2), plot(P, nu(1,:)), hold all;
% plot(P, nu(2,:)), plot(P, nu(3,:)), plot(P, nu(4,:)), plot(P, nu(5,:));
% plot(P, nu(6,:)), plot(P, nu(7,:)), plot(P, nu(8,:));
% title('Dynamic Viscosity [mPa.s]')
% legend('298.15 K', '323.15 K', '348.15 K', '373.15 K', '398.15 K', '423.15 K',...
%     '448.15 K', '473.15 K', 'location', 'northwest')

% Density and viscosity are interpolated to have an exact! result at
% specific conditions
Rho_f = interp2(P, T, D, p, t, 'spline');
Mu_f = interp2(P, T, nu, p, t, 'spline');
function dsp_input_data(handles, data)

% Injector
try set(handles.diamval,'String',num2str(1e6*data.injector.diameter_m)); end% [m] Injector orifice diameter
try 
    set(handles.casel,'Value',data.injector.casel); 
    set(handles.radiobutton2,'Value',1-data.injector.casel);
end
if data.injector.casel
    try set(handles.caval,'String',num2str(data.injector.ca,'%1.2f')); end  % [-] Injector area contraction coefficient
    try set(handles.cdval,'String',num2str(data.injector.cd,'%1.2f')); end  % [m] Injector discharge coefficient
    try set(handles.caval,'Enable','on'); end
    try set(handles.cdval,'Enable','on'); end
    try set(handles.mfrval,'Enable','off'); end
    try set(handles.momentumval,'Enable','off'); end
else
    try set(handles.mfrval,'String',num2str(1e3*data.injector.mfr_kg_s)); end             % [kg/s] steady state mass flow rate
    try set(handles.momentumval,'String',num2str(data.injector.momentum_N)); end        % [N] steady state momentum
    try set(handles.caval,'Enable','off'); end
    try set(handles.cdval,'Enable','off'); end
    try set(handles.mfrval,'Enable','on'); end
    try set(handles.momentumval,'Enable','on'); end
end
try set(handles.thetaval,'String',num2str(data.injector.theta)); end           % [°] Injector cone angle (vapor)
try set(handles.rjetval,'String',num2str(data.injector.Rjet_m)); end          % [m] rather than using single angle, specify jet radial width
try set(handles.xjetval,'String',num2str(data.injector.xjet_m*1000)); end          % [m] axial distances corresponding to jet radial width

% Ambient
try set(handles.tambval,'String',num2str(data.ambient.T_K)); end           % [K] ambient temperature (for liquid length calc only)
try set(handles.rhoambval,'String',num2str(data.ambient.rho_kg_m3)); end    % [kg/m^3] Ambient density

try 
    set(handles.pambsel,'Value',data.ambient.pambsel); 
    set(handles.radiobutton19,'Value',1-data.ambient.pambsel);
end
if data.ambient.pambsel
    try set(handles.pambval,'String',num2str(data.ambient.p_bar,'%3.1f')); end  % [-] Injector area contraction coefficient
    try set(handles.pambval,'Enable','on'); end
    try set(handles.Mfg,'Enable','off'); end
else
    try set(handles.Mfg,'String',num2str(data.ambient.Mfg,'%3.3f')); end             % [kg/s] steady state mass flow rate
    try set(handles.pambval,'Enable','off'); end
    try set(handles.Mfg,'Enable','on'); end
end

try 
    set(handles.xo2sel,'Value',data.ambient.xo2sel); 
    set(handles.radiobutton16,'Value',1-data.ambient.xo2sel);
end
if data.ambient.xo2sel
    try set(handles.xo2val,'String',num2str(data.ambient.xo2,'%3.1f')); end  % [-] Injector area contraction coefficient
    try set(handles.xo2val,'Enable','on'); end
    try set(handles.yo2val,'Enable','off'); end
else
    try set(handles.yo2val,'String',num2str(data.ambient.yo2,'%3.1f')); end             % [kg/s] steady state mass flow rate
    try set(handles.xo2val,'Enable','off'); end
    try set(handles.yo2val,'Enable','on'); end
end

% Fuel
try set(handles.pfuelval,'String',num2str(data.fuel.p_bar)); end
try set(handles.tfuelval,'String',num2str(data.fuel.T_K)); end
try set(handles.stoichval,'String',num2str(data.fuel.pco_air)); end    
if ~data.fuel.autorho
    set(handles.autorhofuel,'Value',data.fuel.autorho);
    set(handles.rhofuelval,'Enable','on');
    try set(handles.rhofuelval,'String',num2str(data.fuel.rho_kg_m3,'%3.1f')); end        % [kg/m^3] Fuel density
else
    set(handles.autorhofuel,'Value',data.fuel.autorho);
    set(handles.rhofuelval,'Enable','off');
end

% Lif-off
if data.lol.autolol
    set(handles.autolol,'Value',data.lol.autolol);
    set(handles.estimlolval,'Enable','off');
else
    try set(handles.estimlolval,'String',num2str(data.lol.lolxp)); end        % [kg/m^3] Fuel density
end

% Options 0d
try set(handles.totallengthval,'String',num2str(data.options0d.total_length_m*1e3)); end	% [meters] Domain size
try set(handles.dxval,'String',num2str(data.options0d.dx_m*1e3)); end          	% [meters] CV widths

if data.options0d.unif_profile % [-] set to 0 for uniform velocity and 1 for Abramovich (real jet)
    set(handles.gaussprof,'Value',0);
    set(handles.unifprof,'Value',1);
    set(handles.alphaman,'Enable','off');
    set(handles.sieberscoeff,'Enable','on');
else
    set(handles.gaussprof,'Value',1);
    set(handles.unifprof,'Value',0);
    set(handles.alphaman,'Enable','on');
    set(handles.sieberscoeff,'Enable','off');
end

try set(handles.alphaman,'String',num2str(data.options0d.alpha)); end

% Vapor penetration
try set(handles.adjsoi,'Value',data.mixing.adjust_soi); end
try set(handles.t1fit,'String',num2str(data.penvap_xp.t1fit,'%1.1f')); end
try set(handles.t2fit,'String',num2str(data.penvap_xp.t2fit,'%1.1f')); end
try set(handles.anglemuscval,'String',num2str(data.mixing.theta,'%2.1f')); end
try set(handles.soimuscval,'String',num2str(data.mixing.soi_corr,'%2.1f')); end

function [velocity,mole_fraction,mass_fraction,FA,EntRate,rhoac,Tadiab,t,alpha_vect,rad_dist,R,area_m2] = function_mixing_fields(xx_m,theta,xjet,Rjet,d0,m0,M0,rhof,rhoa,alpha,ivardens)

% This function computes the mixture fraction and velocity fields in the
% Diesel sprays, following the assumption of the 1d spray model.
% ATTENTION: This is an EXACT solution, not an approximation as proposed by
% Musculus in SAE 2009-01-XXXX (Entrainment Waves)
% If the ambient density is chosen constant, then rhoac = @(x) 0*x

u0 = M0/m0;
x0 = d0/2/tand(theta(1)/2);
area0_m2 = pi*(d0^2)/4;
% area_m2 = pi*((xx_m+x0)*tand(theta/2)).^2;

R = interp1([0 xjet],[0 Rjet],xx_m)+d0/2;
area_m2 = pi*R.^2;

rhoac = 0*xx_m;
alpha_vect = 0*xx_m+alpha;

[uc,yf,fa,xf,b,c,delta] = function_axial_velocity(rhoac,area_m2,m0,M0,rhof,rhoa,alpha_vect);


if ivardens
    [rhoac,Tadiab] = function_adiab_mixing(fa);
    rhoac = rhoac - rhoa;
    [uc,yf,fa,xf,b,c,delta] = function_axial_velocity(rhoac,area_m2,m0,M0,rhof,rhoa,alpha_vect);
end

%% Transition region: If uc>u0, alphaha is modified so that uc=u0
i0 = find(uc>u0);

for ii=i0
    alpha1 = 1;
    alpha2 = 1e4;
    [uc(ii),yf(ii),fa(ii),xf(ii),b(ii),c(ii),delta] = function_axial_velocity(rhoac(ii),area_m2(ii),m0,M0,rhof,rhoa,alpha1);
    err1 = u0^2+b(ii)*u0+c(ii);
    
    [uc(ii),yf(ii),fa(ii),xf(ii),b(ii),c(ii),delta] = function_axial_velocity(rhoac(ii),area_m2(ii),m0,M0,rhof,rhoa,alpha2);
    err2 = u0^2+b(ii)*u0+c(ii);
    
    err = err1;
    while abs(err)>0.1
        alpha_adp = (alpha1+alpha2)/2;
        [uc(ii),yf(ii),fa(ii),xf(ii),b(ii),c(ii),delta] = function_axial_velocity(rhoac(ii),area_m2(ii),m0,M0,rhof,rhoa,alpha_adp);
        err = u0^2+b(ii)*u0+c(ii);
        if err>0
            alpha2 = alpha_adp;
        else
            alpha1 = alpha_adp;
        end
    end
    alpha_vect(ii) = alpha_adp;
end

%% Vapor penetration versus time + Mixing fields
dr = diff(xx_m)*1e3;dr=dr(1);
rad_dist = [-30:dr:30]*1e-3;
I2 = function_stand_int(alpha_vect,2);
% R = R(2:end);

velocity.mean = 2*I2.*uc;    % [m/s] Steady-state axial velocity
velocity.axial = uc;
velocity.field=[];
for ii=1:length(velocity.axial)
    velocity.field(ii,:) = velocity.axial(ii)*(1-min(abs(rad_dist)/R(ii),1).^alpha_vect(ii)).^2;
end

density.axial = rhoa+rhoac;
density.mean = 2*I2.*density.axial;
density.field=[];
for ii=1:length(density.axial)
    density.field(ii,:) = rhoa+rhoac(ii)*(1-min(abs(rad_dist)/R(ii),1).^alpha_vect(ii)).^2;
end

mole_fraction.mean =  2*I2.*yf;  % [-] fuel volume fraction in center of CVs
mole_fraction.axial = yf;
mole_fraction.field=[];
for ii=1:length(mole_fraction.axial)
    mole_fraction.field(ii,:) = mole_fraction.axial(ii)*(1-min(abs(rad_dist)/R(ii),1).^alpha_vect(ii)).^2;
end

FA.mean = 2*I2.*fa;
FA.axial = fa;
FA.field = rhof*mole_fraction.field./(density.field.*(1-mole_fraction.field));

rhoac = [];
[rhoac.mean,Tadiab.mean] = function_adiab_mixing(FA.mean);
[rhoac.axial,Tadiab.axial] = function_adiab_mixing(FA.axial);
[rhoac.field,Tadiab.field] = function_adiab_mixing(FA.field);
if ivardens 
    rhoac.mean = rhoac.mean - rhoa;
    rhoac.axial = rhoac.axial - rhoa;
    rhoac.field = rhoac.field - rhoa;
else
    rhoac.mean = 0*rhoac.mean;
    rhoac.axial = 0*rhoac.axial;
    rhoac.field = 0*rhoac.field;
end

mass_fraction.mean = FA.mean./(FA.mean+1);
mass_fraction.axial = FA.axial./(FA.axial+1);
mass_fraction.field = FA.field./(FA.field+1);

I4 = function_stand_int(alpha_vect,4);
I6 = function_stand_int(alpha_vect,6);
% dMdx = [0 (rhoa+(rhof-rhoa)*mole_fraction.mean).*area_m2.*velocity.mean]; Formule approximé de Mark
dMdx = [0 2*area_m2.*velocity.axial.*((rhof-rhoa)*mole_fraction.axial.*I4+rhoa*I2)]; % Formule exacte
M = cumtrapz([0 xx_m],dMdx); % integrated momentum
t = M/M0;  % time to deliver integrated momentum at maximum nozzle Mdot rate
velocity.tip = velocity.axial .* ((rhof-rhoa).*mole_fraction.axial.*I6+rhoa.*I4)./((rhof-rhoa).*mole_fraction.axial.*I4+rhoa.*I2);

EntRate.raw = diff(dMdx)./diff([0 xx_m]); % [kg/s/m] Air Entrainment
EntRate.norm = diff(dMdx)./diff([0 xx_m])/(m0*1000); % [kg/s/m] Air Entrainment
function [uc,yf,FA,xf,b,c,delta] = function_axial_velocity(rhoac,area_m2,m0,M0,rhof,rhoa,alpha)

u0 = M0/m0;

I2 = function_stand_int(alpha,2);
I4 = function_stand_int(alpha,4);
I6 = function_stand_int(alpha,6);
I8 = function_stand_int(alpha,8);
  
% Exact solution
b = m0*((rhof-rhoa)*I6-rhoac.*I8)./(2*I4*rhof.*(rhoa*I4+rhoac.*I6).*area_m2);
c = -M0./(2*area_m2.*(rhoa*I4+rhoac.*I6));
  
% Mark's approximation
% b = m0*(rhof-rhoa)*I2./(I4*rhof.*rhoa.*area_m2);
% c = -M0./(2*area_m2.*rhoa.*I4);

delta = b.^2-4*c;
uc = (-b+sqrt(delta))/2;
yf =  m0./I4/rhof/2./uc./area_m2;
FA =  rhof*yf./((rhoa+rhoac).*(1-yf));
xf =  FA./(FA+1);
function In = function_stand_int(alpha,n)

% In = integrale(0,1,x*(1-x^alpha)^n dx)
% In= alpha^n*Produit(0,n-1,(n-i)/(2+i*alpha))*1/(n*alpha+2)
In = alpha.^n./(n*alpha+2);
for ii=0:n-1
    In = In*(n-ii)./(2+ii*alpha);
end
function [rhoa_axial,Tadiab_axial] = function_adiab_mixing(FA_axial)

data = getappdata(0,'ListDialogAppData__');

Tint = 473;
Tair =  data.ambient.T_K;
Tfuel = data.fuel.T_K;
Pcell = data.ambient.p_bar;

R = 8.314; % J/mol/K
Mfuel = 170; %g/mol
Mair = data.ambient.Mfg; % g/mol
r_fuel = R/Mfuel*1e3; % J/kg/K
r_air = R/Mair*1e3; % J/kg/K

[cpa_func,cv_func,hl_func,r,mu_func,lambda_func] = Lecture_Prop_Gaz_Single(fullfile(cd,'DataFuel\'),'air.data');
[cpv_func,cv_func,hl_func,r,mu_func,lambda_func] = Lecture_Prop_Gaz_Single(fullfile(cd,'DataFuel\'),'Fuel.data');
[cpl_func,cv_func,hl_func,r,mu_func,lambda_func] = Lecture_Prop_Gaz_Single(fullfile(cd,'DataFuel\'),'Fuel_liquid.data');
[hvap_func,cv_func,hl_func,r,mu_func,lambda_func] = Lecture_Prop_Gaz_Single(fullfile(cd,'DataFuel\'),'hvap.data');
Tmix_adiab = [Tint:0.1:Tair];

% Int(cpa(T)dT,Ta,T)
intcpa = (Tmix_adiab(2)-Tmix_adiab(1)).*(cumsum(cpa_func(Tmix_adiab))-cpa_func(Tmix_adiab(1)));
intcpa = abs(intcpa(end)-intcpa);

% Int(cpv(T)dT,Tint,T)
intcpv = (Tmix_adiab(2)-Tmix_adiab(1)).*(cumsum(cpv_func(Tmix_adiab))-cpv_func(Tmix_adiab(1)));

% Int(cpl(T)dT,Tfuel,Tint)
intcpl = (Tmix_adiab(2)-Tmix_adiab(1)).*(sum(cpl_func([Tfuel:0.1:Tint]))-cpl_func(Tfuel));

% hvap(Tint)
hvap = hvap_func(Tint);

FAtheo = intcpa./(intcpl+hvap+intcpv);

Tadiab_axial = interp1(FAtheo,Tmix_adiab,FA_axial);
Tadiab_axial(isnan(Tadiab_axial)) = Tint;
rhoa_axial = Pcell*1e5/r_air./Tadiab_axial;

% #########################################################################################
% #################################### BUTTONS ACTIONS ####################################
% #########################################################################################
function browse_Callback(hObject, eventdata, handles)
% hObject    handle to browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try rmappdata(0,'ListDialogAppData__'); end
[filename,pathname] = uigetfile('*.mat','Select the injector''s spec. to load');
set(handles.editspec,'String',filename);
load(fullfile(pathname,filename));
setappdata(0,'ListDialogAppData__',data);
try dsp_input_data(handles,data); end
function exportbutton_Callback(hObject, eventdata, handles)
data = getappdata(0,'ListDialogAppData__');

[filename,pathname] = uiputfile;


figure,grid on,hold on
set(gcf,'color','w','position',[411 46 964 932])


%% Axial velocity
subplot(5,2,3)
hold on
plot(data.mixing.x_m*1e3,data.mixing.velocity.mean,'r','linewidth',2)
plot(data.mixing.x_m*1e3,data.mixing.velocity.axial,':r','linewidth',2)
plot(data.mixing.x_m*1e3,data.mixing.velocity.tip,'--r','linewidth',2)
set(gca,'Box','on','FontName','Eras Light ITC','linewidth',2,'fontsize',10,'fontweight','bold','color',[1 1 1]*1,'Xlim',[0 55])
grid on
xlabel('Axial dist. [mm]')
ylabel('Axial Velocity [m/s]')
legend('mean','axial','tip')
hold off

% Velocity field
subplot(5,2,1)
imagesc(data.mixing.x_m*1e3,data.mixing.r_m*1e3,data.mixing.velocity.field'),
set(gca,'Box','on','FontName','Eras Light ITC','linewidth',2,'fontsize',10,'fontweight','bold','color',[1 1 1]*1,'Ylim',[-10 10])
grid on
caxis([0 max(data.mixing.velocity.axial)/4])
ylabel('Velocity Distribution')

%% Axial afr
subplot(5,2,4)
hold on
plot(data.mixing.x_m*1e3,data.mixing.phi.mean,'r','linewidth',2)
plot(data.mixing.x_m*1e3,data.mixing.phi.axial,':r','linewidth',2)
set(gca,'Box','on','FontName','Eras Light ITC','linewidth',2,'fontsize',10,'fontweight','bold','color',[1 1 1]*1)
grid on
ylim = get(gca,'Ylim');
xlim = get(gca,'Xlim');
%ylim(2) = min(ylim(2),10);
%set(gca,'Ylim',[0 10]);
%set(gca,'Xlim',[0 round(max(data.mixing.x_m*1e3))],'Ylim',ylim)
xlabel('Axial dist. [mm]')
ylabel('AFR [-]')
hold off

% AFR field
subplot(5,2,2)
imagesc(data.mixing.x_m*1e3,data.mixing.r_m*1e3,data.mixing.phi.field'),
set(gca,'ydir','reverse')
set(gca,'Box','on','FontName','Eras Light ITC','linewidth',2,'fontsize',10,'fontweight','bold','color',[1 1 1]*1)
grid on
caxis([0 5])


%% Axial fuel mass fraction
subplot(5,2,7)
hold on
plot(data.mixing.x_m*1e3,data.mixing.mass_fraction.mean,'r','linewidth',2)
plot(data.mixing.x_m*1e3,data.mixing.mass_fraction.axial,':r','linewidth',2)
set(gca,'Box','on','FontName','Eras Light ITC','linewidth',2,'fontsize',10,'fontweight','bold','color',[1 1 1]*1)
grid on
ylim = get(gca,'Ylim');
xlim = get(gca,'Xlim');
%ylim(2) = min(ylim(2),10/data.fuel.pco_fg);
set(gca,'Xlim',[0 round(max(data.mixing.x_m*1e3))],'Ylim',ylim)
xlabel('Axial dist. [mm]')
ylabel('Fuel mass fraction [-]')
% set(gca,'ydir','reverse')
hold off

% Fuel mass fraction field
subplot(5,2,5)
imagesc(data.mixing.x_m*1e3,data.mixing.r_m*1e3,data.mixing.mass_fraction.field'),
set(gca,'ydir','reverse')
set(gca,'Box','on','FontName','Eras Light ITC','linewidth',2,'fontsize',10,'fontweight','bold','color',[1 1 1]*1)
grid on
caxis([0 5]/data.fuel.pco_fg)

%% Adiabatic mixing temperature
subplot(5,2,8)
hold on
plot(data.mixing.x_m*1e3,data.mixing.Tadiab.mean,'r','linewidth',2)
plot(data.mixing.x_m*1e3,data.mixing.Tadiab.axial,':r','linewidth',2)
set(gca,'Box','on','FontName','Eras Light ITC','linewidth',2,'fontsize',10,'fontweight','bold','color',[1 1 1]*1)
grid on
ylim = get(gca,'Xlim');
xlim = get(gca,'Ylim');
xlim(2) = min(xlim(2),1000);
set(gca,'Ylim',[300 1000]);
set(gca,'Xlim',[0 round(max(data.mixing.x_m*1e3))],'Xlim',ylim)
xlabel('Axial dist. [mm]')
ylabel('Adiab. mix. temp. [K]')
% set(gca,'ydir','reverse')
hold off

% Temperature field
subplot(5,2,6)
imagesc(data.mixing.x_m*1e3,data.mixing.r_m*1e3,data.mixing.Tadiab.field'),
caxis([max(data.mixing.Tadiab.field(:))-200 max(data.mixing.Tadiab.field(:))])
set(gca,'ydir','reverse')
set(gca,'Box','on','FontName','Eras Light ITC','linewidth',2,'fontsize',10,'fontweight','bold','color',[1 1 1]*1)
grid on


%% Air Entrainment
subplot(5,2,10)
hold on
%plot(data.mixing.x_m*1e3,data.mixing.AirMFR.raw*1e3,'r','linewidth',2)
plot(data.mixing.x_m*1e3,data.mixing.EntRate.raw,'r--','linewidth',2)
plot(data.mixing.x_m*1e3,data.mixing.EntRate.norm,'r:','linewidth',2)
set(gca,'Box','on','FontName','Eras Light ITC','linewidth',2,'fontsize',10,'fontweight','bold','color',[1 1 1]*1)
grid on
xlabel('Distance [mm]')
ylabel('Air MFR [g/s]')
hold off

%% Spray penetration
subplot(5,2,9)
hold on
try
    plot(data.penvap_xp.time_mus/1000,data.penvap_xp.penvap_mm,'ok','linewidth',2);
    plot(data.penvap_xp.time_mus/1000,data.penvap_xp.penvap_mm+data.penvap_xp.std_penvap_mm,':k','linewidth',1);
    plot(data.penvap_xp.time_mus/1000,data.penvap_xp.penvap_mm-data.penvap_xp.std_penvap_mm,':k','linewidth',1);
end

try
    i1fit = data.penvap_xp.i1fit;
    i2fit = data.penvap_xp.i2fit;
    plot(data.penvap_xp.time_mus(i1fit:i2fit)/1000,data.penvap_xp.penvap_mm(i1fit:i2fit),'o','linewidth',2,'color',[1 1 1]*0.7);
end

try
    plot((data.penvap_xp.time_mus+data.mixing.soi_corr)/1000,data.penvap_xp.penvap_mm,'+k','linewidth',2);
end

try
    i1fit = data.penvap_xp.i1fit;
    i2fit = data.penvap_xp.i2fit;
    plot((data.penvap_xp.time_mus(i1fit:i2fit)-data.penvap_xp.soi_corr)/1000,data.penvap_xp.penvap_mm(i1fit:i2fit),'+','linewidth',2,'color',[1 1 1]*0.7);
end
plot(data.mixing.t_s*1e3,[0 data.mixing.x_m]*1e3,'r','linewidth',2)
set(gca,'Box','on','FontName','Eras Light ITC','linewidth',2,'fontsize',10,'fontweight','bold','color',[1 1 1]*1)
grid on
xlim = get(gca,'Xlim');
ylim = get(gca,'Ylim');
set(gca,'Xlim',[0 100],'Ylim',[0 6]);
set(gca,'Xlim',xlim,'Ylim',ylim)
xlabel('Time [ms]')
ylabel('Distance [mm]')
%ht=title(data.penvap_xp.folder);
%set(ht,'fontSize',6,'fontname','eras light ITC');
hold off

saveas(gcf,fullfile(pathname,[filename '.fig']),'fig');
function savebutton_Callback(hObject, eventdata, handles)
% hObject    handle to savebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data = getappdata(0,'ListDialogAppData__');
try data.mixing = rmfield(data.mixing, 't'); end
[filename,pathname] = uiputfile('*.mat','Save data as...');
save(fullfile(pathname,filename),'data');
set(handles.editspec,'String',filename);


% #########################################################################################
% ################################### CREATE OBJECTS ######################################
% #########################################################################################

function thetacorrval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thetacorrval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function phimusculusaxial_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phimusculusaxial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function phiexactaxial_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phiexactaxial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function phisiebersaxial_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phisiebersaxial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function phimusculusmean_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phimusculusmean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function phiexactmean_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phiexactmean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function phisiebersmean_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phisiebersmean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function anglemuscval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to anglemuscval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function export_plot(handles)

data = getappdata(0,'ListDialogAppData__');

figure,
set(gcf,'color','w','Position',[203 29 1331 956]);

% Axial velocity
subplot(3,2,1),grid on
hold on
if get(handles.musculusvel,'Value')
    plot(data.outputs.x_m*1e3,data.outputs.U_mean_musculus,'b','linewidth',2)
    plot(data.outputs.x_m*1e3,data.outputs.U_axial_musculus,':b','linewidth',2)
    plot(data.outputs.x_m*1e3,data.outputs.U_tip_musculus,'--b','linewidth',2)
end

if get(handles.exactvel,'Value');
    plot(data.outputs.x_m*1e3,data.outputs.U_mean_exact,'r','linewidth',2)
    plot(data.outputs.x_m*1e3,data.outputs.U_axial_exact,':r','linewidth',2)
    plot(data.outputs.x_m*1e3,data.outputs.U_tip_exact,'--r','linewidth',2)
end

if get(handles.siebersvel,'Value');
    plot(data.outputs.x_m*1e3,data.outputs.U_mean_siebers,'--','Color',[0 0.5 0],'linewidth',2)
    plot(data.outputs.x_m*1e3,data.outputs.U_axial_siebers,'--','Color',[0 0.5 0],'linewidth',2)
end

set(gca,'Box','on','FontName','Calibri','linewidth',2,'fontsize',10,'fontweight','bold','color',[1 1 1]*1)
grid on
xlabel('Distance [mm]')
ylabel('Axial Velovity [m/s]')
hold off

% Axial afr
subplot(3,2,3),grid on
hold on
if get(handles.musculusvel,'Value')
    plot(data.outputs.x_m*1e3,data.outputs.phi_mean_musculus,'b','linewidth',2)
    plot(data.outputs.x_m*1e3,data.outputs.phi_axial_musculus,':b','linewidth',2)
    %     plot(data.outputs.x_m*1e3,data.outputs.phi_tip_musculus,'--b','linewidth',2)
end

if get(handles.exactvel,'Value');
    plot(data.outputs.x_m*1e3,data.outputs.phi_mean_exact,'r','linewidth',2)
    plot(data.outputs.x_m*1e3,data.outputs.phi_axial_exact,':r','linewidth',2)
    %     plot(data.outputs.x_m*1e3,data.outputs.phi_tip_exact,'--r','linewidth',2)
end

if get(handles.siebersvel,'Value');
    plot(data.outputs.x_m*1e3,data.outputs.phi_mean_siebers,'--','Color',[0 0.5 0],'linewidth',2)
    plot(data.outputs.x_m*1e3,data.outputs.phi_axial_siebers,'--','Color',[0 0.5 0],'linewidth',2)
end

set(gca,'Box','on','FontName','Calibri','linewidth',2,'fontsize',10,'fontweight','bold','color',[1 1 1]*1)
grid on
xlim = get(gca,'Xlim');
ylim = get(gca,'Ylim');
ylim(2) = min(ylim(2),10);
set(gca,'Ylim',[0 10]);
set(gca,'Xlim',xlim,'Ylim',ylim)
xlabel('Distance [mm]')
ylabel('AFR [-]')
hold off


% Axial fa
subplot(3,2,2),grid on
hold on
if get(handles.musculusvel,'Value')
    plot(data.outputs.x_m*1e3,data.outputs.FA_mean_musculus,'b','linewidth',2)
    plot(data.outputs.x_m*1e3,data.outputs.FA_axial_musculus,':b','linewidth',2)
    %     plot(data.outputs.x_m*1e3,data.outputs.FA_tip_musculus,'--b','linewidth',2)
end

if get(handles.exactvel,'Value');
    plot(data.outputs.x_m*1e3,data.outputs.FA_mean_exact,'r','linewidth',2)
    plot(data.outputs.x_m*1e3,data.outputs.FA_axial_exact,':r','linewidth',2)
    %     plot(data.outputs.x_m*1e3,data.outputs.FA_tip_exact,'--r','linewidth',2)
end

if get(handles.siebersvel,'Value');
    plot(data.outputs.x_m*1e3,data.outputs.FA_mean_siebers,'--','Color',[0 0.5 0],'linewidth',2)
    plot(data.outputs.x_m*1e3,data.outputs.FA_axial_siebers,'--','Color',[0 0.5 0],'linewidth',2)
end

set(gca,'Box','on','FontName','Calibri','linewidth',2,'fontsize',10,'fontweight','bold','color',[1 1 1]*1)
grid on
xlim = get(gca,'Xlim');
ylim = get(gca,'Ylim');
ylim(2) = min(ylim(2),10/data.fuel.pco_aircorr);
set(gca,'Ylim',[0 10]/data.fuel.pco_aircorr)
set(gca,'Xlim',xlim,'Ylim',ylim)
xlabel('Distance [mm]')
ylabel('F/A [-]')
hold off

% Spray penetration
subplot(3,2,4),grid on
hold on
try
    plot(data.penvap_xp.time_mus/1000,data.penvap_xp.penvap_mm,'ok','linewidth',2);
end

try
    plot((data.penvap_xp.time_mus-data.penvap_xp.soiadj_mus_exact)/1000,data.penvap_xp.penvap_mm,'+k','linewidth',2);
end

try
    if data.penvap_xp.plotapprox
        plot(data.penvap_xp.t_musculus_s*1000+data.penvap_xp.soiadj_mus_musculus/1000,data.penvap_xp.x_m*1000,':b');
        plot(data.penvap_xp.t_exact_s*1000+data.penvap_xp.soiadj_mus_exact/1000,data.penvap_xp.x_m*1000,':r');
        plot(data.penvap_xp.t_siebers_s*1000+data.penvap_xp.soiadj_mus_siebers/1000,data.penvap_xp.x_m*1000,':g');
    end
end

if get(handles.musculusvel,'Value')
    plot(data.outputs.t_musculus*1e3,[0 data.outputs.x_m]*1e3,'b','linewidth',2)
end

if get(handles.exactvel,'Value');
    plot(data.outputs.t_exact*1e3,[0 data.outputs.x_m]*1e3,'r','linewidth',2)
end

if get(handles.siebersvel,'Value');
    plot(data.outputs.t_siebers*1e3,[0 data.outputs.x_m]*1e3,'Color',[0 0.5 0],'linewidth',2)
end

set(gca,'Box','on','FontName','Calibri','linewidth',2,'fontsize',10,'fontweight','bold','color',[1 1 1]*1)
grid on
xlim = get(gca,'Xlim');
ylim = get(gca,'Ylim');
set(gca,'Xlim',[0 6],'Ylim',[0 100]);
set(gca,'Xlim',xlim,'Ylim',ylim)
xlabel('Time [ms]')
ylabel('Distance [mm]')
hold off

% Adiabatic mixing temperature
subplot(3,2,5),grid on
hold on
if get(handles.musculusvel,'Value')
    plot(data.outputs.x_m*1e3,data.outputs.Tfuel_mean_musculus,'b','linewidth',2)
    plot(data.outputs.x_m*1e3,data.outputs.Tfuel_axial_musculus,':b','linewidth',2)
end

if get(handles.exactvel,'Value');
    plot(data.outputs.x_m*1e3,data.outputs.Tfuel_mean_exact,'r','linewidth',2)
    plot(data.outputs.x_m*1e3,data.outputs.Tfuel_axial_exact,':r','linewidth',2)
end

if get(handles.siebersvel,'Value');
    plot(data.outputs.x_m*1e3,data.outputs.Tfuel_mean_siebers,'Color',[0 0.5 0],'linewidth',2)
    plot(data.outputs.x_m*1e3,data.outputs.Tfuel_axial_siebers,':','Color',[0 0.5 0],'linewidth',2)
end

set(gca,'Box','on','FontName','Calibri','linewidth',2,'fontsize',10,'fontweight','bold','color',[1 1 1]*1)
grid on
xlim = get(gca,'Xlim');
ylim = get(gca,'Ylim');
ylim(2) = min(ylim(2),1000);
set(gca,'Ylim',[300 1000]);
set(gca,'Xlim',xlim,'Ylim',ylim);
xlabel('Distance [mm]')
ylabel('Adiab. mix. temp. [K]')
hold off

% Air Entrainment
subplot(3,2,6),grid on
hold on
if get(handles.musculusvel,'Value')
    plot(data.outputs.x_m*1e3,data.outputs.AirMFR_musculus*1e3,'b','linewidth',2)
    plot(data.outputs.x_m*1e3,data.outputs.EntRate_musculus,'b','linewidth',2)
end

if get(handles.exactvel,'Value');
    plot(data.outputs.x_m*1e3,data.outputs.AirMFR_exact*1e3,'r','linewidth',2)
    plot(data.outputs.x_m*1e3,data.outputs.EntRate_exact,'r','linewidth',2)
end

if get(handles.siebersvel,'Value');
    plot(data.outputs.x_m*1e3,data.outputs.AirMFR_siebers*1e3,'Color',[0 0.5 0],'linewidth',2)
    plot(data.outputs.x_m*1e3,data.outputs.EntRate_siebers,'Color',[0 0.5 0],'linewidth',2)
end

set(gca,'Box','on','FontName','Calibri','linewidth',2,'fontsize',10,'fontweight','bold','color',[1 1 1]*1)
grid on
xlabel('Distance [mm]')
ylabel('Air MFR [g/s]')
hold off
function angleexval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to angleexval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function anglesiebval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to anglesiebval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function soimuscval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to soimuscval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function soiexval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to soiexval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function cvval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cvval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function zstval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zstval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function soisiebval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to soisiebval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit59_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zstval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function effdiamval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to effdiamval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function stoichcorrval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stoichcorrval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function u0val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to u0val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function t1fit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t1fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function t2fit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t2fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function coeflolval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to coeflolval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function estimlolval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to estimlolval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function totallengthval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to totallengthval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function dxval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dxval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function alphaman_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alphaman (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function sieberscoeff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sieberscoeff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function tambval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tambval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function xo2val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xo2val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pambval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pambval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function rhoambval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rhoambval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Mfg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Mfg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function yo2val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yo2val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function rhofuelval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rhofuelval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function stoichval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stoichval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pfuelval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pfuelval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function tfuelval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tfuelval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function diamval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to diamval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function caval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to caval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function cdval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cdval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function thetaval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thetaval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function rjetval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rjetval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function xjetval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xjetval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function mfrval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mfrval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function momentumval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to momentumval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function xoffset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xoffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function toffset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to toffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
