function varargout = FEMIC_inverse5(varargin)
% FEMIC_INVERSE5 MATLAB code for FEMIC_inverse5.fig
%      FEMIC_INVERSE5, by itself, creates a new FEMIC_INVERSE5 or raises the existing
%      singleton
%
%      H = FEMIC_INVERSE5 returns the handle to a new FEMIC_INVERSE5 or the handle to
%      the existing singleton*.
%
%      FEMIC_INVERSE5('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FEMIC_INVERSE5.M with the given input arguments.
%
%      FEMIC_INVERSE5('Property','Value',...) creates a new FEMIC_INVERSE5 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FEMIC_inverse5_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FEMIC_inverse5_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FEMIC_inverse5

% Last Modified by GUIDE v2.5 12-Jul-2013 11:32:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FEMIC_inverse5_OpeningFcn, ...
                   'gui_OutputFcn',  @FEMIC_inverse5_OutputFcn, ...
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
% End initialization code - DO NOT EDIT

% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function text_dt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_dt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes just before FEMIC_inverse5 is made visible.
function FEMIC_inverse5_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FEMIC_inverse5 (see VARARGIN)

% Choose default command line output for FEMIC_inverse5
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FEMIC_inverse5 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FEMIC_inverse5_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function uitable_offset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitable_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
freq=[1530;8250;20370;33030;47970;63090;93090];
offsets=1.66*ones(7,1);
data=[freq offsets];
set(hObject, 'Data', data);
guidata(hObject, handles);


% --- Executes on button press in checkbox_plotlcurve.
function checkbox_plotlcurve_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_plotlcurve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_plotlcurve


% --- Executes on button press in pushbutton_startInversion.
function pushbutton_startInversion_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_startInversion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% code from original FEMIC:
% ***[p.model,p.mu,p.err,p.G]=FEMIC_inverse2D ...
%   (mdlparms.mdlECa, mdlparms.mdlD, dataparms.obs, dataparms.obsstd, ....
%   dataparms.f, dataparms.r, mdlparms.muh, mdlparms.muv, cnvparms.ecatol, ...
%   cnvparms.errtol, cnvparms.mxiter, dataparms.dtypes, mdlparms.pmin,...
%   mdlparms.pmax, mdlparms.bar, invType, ...
%   [cnvparms.plotInter cnvparms.plotL cnvparms.wprof]);');

% [sound_id, layer_thick, Eca, constraint_wgt]
initmodel_data=get(handles.uitable_initmodel, 'Data');
max_id=max(initmodel_data(:, 1));

%get priori information - all model data where constraint weight > 0
all_priori=initmodel_data(initmodel_data(:, 5)>0, :);
%[sound_id, layer number, conductivity, constraint weight]
priori=[all_priori(:, 1), all_priori(:, 2), all_priori(:, 4), all_priori(:, 5)];
 
 
%format [x, y, z, sound_id, freq, STD_error, HDM/VDM flag, Eca]
input_data=get(handles.uitable_datainput, 'Data');


%get all sound_ids
sound_id=unique(input_data(:, 4));

%validate - make sure same number of soundings in data and model
if length(sound_id)> max_id
    msgbox('There is a greater number of soundings in data input than initial model');
    return;
elseif length(sound_id) < max_id
    msgbox('There is a greater number of initial model soundings than in data input.');
    return;
end

oneD=get(handles.radiobutton_1D, 'Value');
twoD=get(handles.radiobutton_2D, 'Value');
threeD=get(handles.radiobutton_3D, 'Value');
if length(sound_id)>1 && oneD==1
    msgbox('There is more than one sounding in your data.  Select 2D or 3D inversion.');
    return;
elseif length(sound_id)==1
    if oneD==0
        msgbox('There is only one sounding in your data.  Select 1D inversion.');
        return;
    end        
end
 
 
meas=[];
STD_error=[];
coords=[];
for ii=1:size(sound_id, 1) 
    %set up model conductivity/layer thicknesses
    model_cond{ii}=initmodel_data((initmodel_data(:, 1)==ii), 4); %mdlECa
    layer_thick{ii}=initmodel_data((initmodel_data(:, 1)==ii), 3); %mdlD
    
    data_sound=input_data(input_data(:, 4)==sound_id(ii), :); %data for each sounding
    data_sound=sortrows(data_sound, 5); %sort by freq
    data_sound=sortrows(data_sound, 7); %sort by vdm/hdm flag
    
    add_meas=imag(data_sound(:, 8));
    meas=[meas, add_meas]; %obs 
    
    add_error=data_sound(:, 6);
    STD_error=[STD_error, add_error]; %obstd
    
    s_id=data_sound(1, 4);
    x=data_sound(1,1);
    y=data_sound(1,2);
    z=data_sound(1,3);
    add_coords=[s_id, x, y, z];
    coords=[coords; add_coords];
end
%dlmwrite('meass.dat',meas);
offset_data=get(handles.uitable_offset, 'Data'); 
offset_data=sortrows(offset_data, 1);

%get unique frequencies
datainput_freq=unique(input_data(:, 5));

freq=offset_data(:, 1); %f

%if length(datainput_freq)~=length(freq)
 %   msgbox('Check Tx-Rx offset and frequency data.  Number of frequencies in data input does not match number frequency/offset data.');
  %  return;
%end

offset=offset_data(:, 2); %r
 
horizreg=str2num(get(handles.edit_horizreg, 'String')); %muh

vertreg_cell= get(handles.listbox_vertreg, 'String'); %muv
s=size(vertreg_cell, 1);
for ii=1:s
   vertreg_cell{ii}=str2num(vertreg_cell{ii});
end
vertreg=cell2mat(vertreg_cell);
el=str2num(get(handles.edit_elevation, 'String'));
sx=str2num(get(handles.edit_sx, 'String'));
sy=str2num(get(handles.edit_sy, 'String'));
sz=str2num(get(handles.edit_sy, 'String'));
ecatol=str2num(get(handles.edit_ecatol, 'String')); %ecatol
errtol=str2num(get(handles.edit_errtol, 'String')); %errtol
maxiter=str2num(get(handles.edit_maxiter, 'String')); %mxiter

dt_str=get(handles.text_dt, 'String'); %dtypes
switch dt_str
    case 'VMD'
      dt=1;
    case 'HMD'
      dt=2;
    case 'Both VMD and HMD'  
        dt=3;
end
data=get(handles.uitable_datainput, 'Data');
dtt=data(:,7);
pmin=str2num(get(handles.edit_minparam, 'String')); %pmin
pmax=str2num(get(handles.edit_maxparam, 'String')); %pmax
bar=str2num(get(handles.edit_initbarrier, 'String')); %bar

%see which inversion type is selected
invType=get(handles.popupmenu_regcriteria, 'Value');

%plotParms
%plotInter=get(handles.checkbox_interimPlot, 'Value');
plotL=get(handles.checkbox_plotlcurve, 'Value');
plotSens=get(handles.checkbox_plotsensit, 'Value');
%plotNprof=str2num(get(handles.edit_profilenum, 'String'));

plotParms=[plotSens plotL]; % [plotInter plotL wprof]

%see which radio button - inversion dimension - is selected
%oneD, twoD and threeD defined above for validation check

statusUpdate = @updateStatusBar;

if oneD==1 
    statusUpdate('Starting 1D Inversion....');
    [p_final, mu_final, rms_error, G]=FEMIC_inverse1D(model_cond, layer_thick, ...
        meas, STD_error,coords,freq,offset,horizreg,vertreg, ...
        ecatol, errtol, maxiter, dt, pmin, pmax, bar, invType, priori, sz, ...
        plotParms, statusUpdate,el);
    results.p_final=p_final;
    results.mu_final=mu_final;
    results.rms_error=rms_error;
    results.G=G;
elseif twoD==1    
    statusUpdate('Starting 2D Inversion....');
    vall2=str2num(get(handles.parallel,'String'));
    plotdoi=get(handles.checkbox_DOI, 'Value');
    perc=str2num(get(handles.percent,'String'));
    [p_final, mu_final, rms_error, G, xi, ei]=FEMIC_inverse2D(model_cond, layer_thick,...
        meas, STD_error, coords, freq, offset, horizreg, vertreg, ...
        ecatol, errtol, maxiter, dtt, pmin, pmax, bar, invType, priori, sx, sz, ...
        statusUpdate, plotSens,vall2,plotdoi,perc,el); 
    results.p_final=p_final;
    results.mu_final=mu_final;
    results.rms_error=rms_error;
    results.xi=xi;
    results.ei=ei;
    results.G=G;
     
elseif threeD==1 
    statusUpdate('Starting 3D Inversion....');
    vall=str2num(get(handles.parallel,'String'));
    b=findResource();
    b=b.ClusterSize;
    plotdoi=get(handles.checkbox_DOI, 'Value');
    perc=str2num(get(handles.percent,'String'));
    if vall>b
        prompt('error: number of processor exceeded');end  
    warning('off');
    [p_final, mu_final, rms_error, G]=FEMIC_inverse3D(model_cond, layer_thick, ...
        meas,STD_error, freq, offset,coords, horizreg,vertreg,...
        ecatol, errtol, maxiter, dt, pmin, pmax, bar, invType, priori, sx, sy, sz, ...
        statusUpdate, plotSens,vall,plotdoi, perc,el);
    
    results.p_final=p_final;
    results.mu_final=mu_final;
    results.rms_error=rms_error;
    results.G=G;
end

save('result.mat', '-struct', 'results');

delete('results.xls');
%dlmwrite('model.dat', p_final);
%dlmwrite('rms_error.dat', rms_error);
%dlmwrite('estimted_data.dat', G);

%save modeel/inversion parameters
fid=fopen('modelInv_parms.dat', 'w');

fprintf(fid, '%4.2f\r\n', pmin);
fprintf(fid, '%4.2f\r\n', pmax);
fprintf(fid, '%4.2f\r\n', bar);
fprintf(fid, '%6.15f\r\n', horizreg);
for ii=1:size(vertreg, 1)
    fprintf(fid, '%6.15f\r\n', vertreg(ii));
end
fprintf(fid, '%4.2f\r\n', maxiter);
fprintf(fid, '%4.2f\r\n', errtol);
fprintf(fid, '%4.2f\r\n', ecatol);
fprintf(fid, '%4.2f\r\n', sx);
fprintf(fid, '%4.2f\r\n', sy);
fprintf(fid, '%4.2f\r\n', sz);
fprintf(fid, '%f\r\n', invType);
fprintf(fid, '%f\r\n', oneD);
fprintf(fid, '%f\r\n', twoD);
fprintf(fid, '%f\r\n', threeD);
fclose(fid);
%end save modeel/inversion parameters
guidata(hObject, handles);

function updateStatusBar(msg)
    h = findobj('Tag','text_status');
    set(h, 'String', msg);

        
% --- Executes on button press in pushbutton_Cancel.
function pushbutton_Cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1);


% --- Executes on button press in pushbutton_modelparamFile.
function pushbutton_modelparamFile_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_modelparamFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[f_modelparam,pname,fIndex] = uigetfile('*.*', 'Select model-inversion parameters input file:');
if fIndex~=0
   
    fid=fopen(strcat(pname, f_modelparam), 'r');
    f = fread(fid,Inf);
    s = char(f');  % read in the file into a single character string
    [modelparam, status] = str2num(s);
    if status==0
        msgbox('Your model-inversion parameters file selection contains strings.  Select a different file.', 'Incorrect file format');
        return;
    end
    fclose(fid);
    
    pmin = modelparam(1,1);
    pmax = modelparam(2,1);
    bar = modelparam(3,1);
    muh = modelparam(4,1);
    muv = modelparam(5:end-10, 1);
    muv={muv}; %convert to cell
    
    mxiters = modelparam(end-9,1);
    errtol = modelparam(end-8,1);
    ecatol = modelparam(end-7,1);
    sx = modelparam(end-6,1);
    sy = modelparam(end-5,1);
    sz = modelparam(end-4,1);
    inv_type=modelparam(end-3,1);
    oneD=modelparam(end-2, 1);
    twoD=modelparam(end-1, 1);
    threeD=modelparam(end, 1);
    
    
    set(handles.edit_minparam, 'String', pmin);
    set(handles.edit_maxparam, 'String', pmax);
    set(handles.edit_initbarrier, 'String', bar);
    set(handles.edit_horizreg, 'String', muh);
    set(handles.listbox_vertreg, 'String', muv);
    set(handles.edit_maxiter, 'String', mxiters);
    set(handles.edit_errtol, 'String', errtol);
    set(handles.edit_ecatol, 'String', ecatol);    
    
    set(handles.text_modelparamFile, 'String', f_modelparam);
    set(handles.popupmenu_regcriteria, 'Value', inv_type);
        
    set(handles.radiobutton_1D, 'Value', oneD);
    set(handles.radiobutton_2D, 'Value', twoD);
    set(handles.radiobutton_3D, 'Value', threeD);
       
    if oneD==1
        if inv_type==4
            set(handles.checkbox_plotlcurve, 'Visible', 'on');  
            set(handles.uipanel11, 'Visible', 'on');
        else
            set(handles.percent,'Visible','off')
            set(handles.checkbox_plotlcurve, 'Visible', 'off');  
            set(handles.uipanel11, 'Visible', 'off');
        end
        set(handles.checkbox_plotsensit, 'Visible', 'off');
        set(handles.text_smooth, 'String', 'Smoothness (z)');    
        set(handles.edit_sx, 'Visible', 'off');
        set(handles.edit_sy, 'Visible', 'off');
    elseif twoD==1
        set(handles.percent,'Visible','on')
        set(handles.uipanel11, 'Visible', 'on');
        set(handles.checkbox_plotlcurve, 'Visible', 'off');
        set(handles.checkbox_plotsensit, 'Visible', 'on');
        set(handles.text_smooth, 'String', 'Smoothness (x, z)');
        set(handles.edit_sx, 'Visible', 'on'); 
        set(handles.edit_sy, 'Visible', 'off'); 
    elseif threeD==1
        set(handles.percent,'Visible','on')
        set(handles.uipanel11, 'Visible', 'on');
        set(handles.checkbox_plotlcurve, 'Visible', 'off');
        set(handles.checkbox_plotsensit, 'Visible', 'on');
        set(handles.text_smooth, 'String', 'Smoothness (x, y, z)');
        set(handles.edit_sx, 'Visible', 'on'); 
        set(handles.edit_sy, 'Visible', 'on'); 
    end
    
end


% --- Executes on button press in pushbutton_viewData.
function pushbutton_viewData_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_viewData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h=figure('WindowStyle', 'modal', 'Name', 'View data input');

%create table viewing parameters
data=get(handles.uitable_datainput, 'Data');

s=size(data, 2);
isEditable=zeros(1,s);
isEditable = logical(isEditable);

t_data=uitable(h, 'Tag', 'uitable_datainput', 'Position', [30 80 500 300],...
    'Data', get(handles.uitable_datainput, 'Data'), ...
    'FontSize', 10,...
    'ColumnEditable', isEditable, ...
    'ColumnWidth', 'auto',...
    'ColumnName', {'x', 'y', 'z', 'Sounding id', 'Freq', 'STD Error', 'VMD[1],HMD[2]', 'Quadrature'});  

cancel_but=uicontrol(h, 'Style', 'pushbutton', 'String', 'Cancel',...
    'FontSize', 10, 'Position', [460 30 70 30],...
    'Callback', 'close');
%make editable?
%saveandclose_but=uicontrol(h, 'Style', 'pushbutton', 'String', 'Save and Close',...
%    'FontSize', 10, 'Position', [305 30 150 30],...
%    'Callback', '@updateuitable_data');

% --- Executes on button press in pushbutton_getDatafile.
function pushbutton_getDatafile_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_getDatafile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fname,pname,fIndex] = uigetfile('*.*', 'Select data input');%fname is filename, pname is directory
if fIndex~=0

    %fid=fopen(strcat(pname), 'r');
    %f = fread(fid,Inf);
    data = GEM2FEM(strcat(pname,fname)); %added by Neil
    %dlmwrite('data.dat',data)
    %s = char(f');  % read in the file into a single character string

    %[data, status] = str2num(s);
    %if status==0
    %    msgbox('Your data input file selection contains strings.  Select a different file.', 'Incorrect file format');
    %    return;
    %end
      
    set(handles.uitable_datainput, 'Data', data);
    set(handles.text_dataFile, 'String', fname);
     
    %populate the uitable_coords [sound_id x y z]
    %get all sound_ids
    sound_id=data(:,4); 
    s=unique(sound_id); 
    coords=[];
    for ii=1:length(s)
        ind=find(data(:, 4)==s(ii));
        x=data(ind(1), 1);
        y=data(ind(1), 2);
        z=data(ind(1), 3);
        add_coords=[s(ii) x y z];
        coords(ii, 1:4)=[add_coords];
    end
    set(handles.uitable_coords, 'Data', coords);
    
    
    %see if VDM/HDM or both
    flag=data(:, 7);   
    vdm=find(flag==1);
    hdm=find(flag==2);
    if ~isempty(vdm) && ~isempty(hdm)
        set(handles.text_dt, 'String', 'Both VMD and HMD');
    elseif ~isempty(hdm)
        set(handles.text_dt, 'String', 'HMD');
    elseif ~isempty(vdm)
        set(handles.text_dt, 'String', 'VMD');
    end
    
    %see if initial model has same number of soundings as data input
    data=get(handles.uitable_datainput, 'Data');
    sound_id=unique(data(:, 4));
    num_soundings=length(sound_id);
    
    init_model=get(handles.uitable_initmodel, 'Data');
    model_soundid=unique(init_model(:, 1));
    num_modelsoundings=length(model_soundid);
    
    %{
    if num_soundings~=num_modelsoundings
        str={'Number of soundings in data input'; 'does not match'; 'number of soundings in initial model.'; 'Review initial model.'};
        warndlg(str);
    end
    %}
    
end
% --- Executes on button press in pushbutton_viewOffset.
function pushbutton_viewOffset_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_viewOffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h=figure('WindowStyle', 'modal', 'Name', 'View Tx-Rx offset / Frequency data');

%create table viewing parameters
data=get(handles.uitable_offset, 'Data');
s=size(data, 2);
isEditable=zeros(1,s);
isEditable = logical(isEditable);

t_data=uitable(h, 'Tag', 'uitable_STDerror', 'Position', [30 80 500 300],...
    'Data', get(handles.uitable_offset, 'Data'), ...
    'FontSize', 10,...
    'ColumnEditable', isEditable,...
    'ColumnName', {'Frequency (Hz)', 'Tx-Rx offset (m)'});  

cancel_but=uicontrol(h, 'Style', 'pushbutton', 'String', 'Cancel',...
    'FontSize', 10, 'Position', [460 30 70 30],...
    'Callback', 'close');


% --- Executes on button press in pushbutton_getOffsetfile.
function pushbutton_getOffsetfile_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_getOffsetfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[f_offset,pname,fIndex] = uigetfile('*.*', 'Select Tx-Rx offsets and frequencies');
if fIndex~=0   
    fid=fopen(strcat(pname, f_offset), 'r');
    f = fread(fid,Inf);
    s = char(f');  % read in the file into a single character string
    [data, status] = str2num(s);
    if status==0
        msgbox('Your Tx-Rx offset and frequency file contains strings.  Select a different file.', 'Incorrect file format');
        return;
    end
    
    set(handles.uitable_offset, 'Data', data);
    set(handles.text_offsetFile, 'String', f_offset);
    
end



function edit_minparam_Callback(hObject, eventdata, handles)
% hObject    handle to edit_minparam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_minparam as text
%        str2double(get(hObject,'String')) returns contents of edit_minparam as a double
val=get(hObject, 'String');
[val, status]=str2num(val);
if ~status
    msgbox('Entry must be numeric', 'Incorrect data format');
    set(hObject, 'String', '0.0');
    return;
end

% --- Executes during object creation, after setting all properties.
function edit_minparam_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_minparam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_maxparam_Callback(hObject, eventdata, handles)
% hObject    handle to edit_maxparam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_maxparam as text
%        str2double(get(hObject,'String')) returns contents of edit_maxparam as a double
val=get(hObject, 'String');
[val, status]=str2num(val);
if ~status
    msgbox('Entry must be numeric', 'Incorrect data format');
    set(hObject, 'String', '1000.0');
    return;
end


% --- Executes during object creation, after setting all properties.
function edit_maxparam_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_maxparam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_initbarrier_Callback(hObject, eventdata, handles)
% hObject    handle to edit_initbarrier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_initbarrier as text
%        str2double(get(hObject,'String')) returns contents of edit_initbarrier as a double
val=get(hObject, 'String');
[val, status]=str2num(val);
if ~status
    msgbox('Entry must be numeric', 'Incorrect data format');
    set(hObject, 'String', '2.0');
    return;
end

% --- Executes during object creation, after setting all properties.
function edit_initbarrier_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_initbarrier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_horizreg_Callback(hObject, eventdata, handles)
% hObject    handle to edit_horizreg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_horizreg as text
%        str2double(get(hObject,'String')) returns contents of edit_horizreg as a double
val=get(hObject, 'String');
[val, status]=str2num(val);
if ~status
    msgbox('Entry must be numeric', 'Incorrect data format');
    set(hObject, 'String', '10e-6');
    return;
end

% --- Executes during object creation, after setting all properties.
function edit_horizreg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_horizreg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_vertreg.
function popupmenu_vertreg_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_vertreg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_vertreg contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_vertreg


% --- Executes during object creation, after setting all properties.
function popupmenu_vertreg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_vertreg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: text_dt controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_maxiter_Callback(hObject, eventdata, handles)
% hObject    handle to edit_maxiter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_maxiter as text
%        str2double(get(hObject,'String')) returns contents of edit_maxiter as a double
val=get(hObject, 'String');
[val, status]=str2num(val);
if ~status
    msgbox('Entry must be numeric', 'Incorrect data format');
    set(hObject, 'String', '10');
    return;
end

% --- Executes during object creation, after setting all properties.
function edit_maxiter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_maxiter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_errtol_Callback(hObject, eventdata, handles)
% hObject    handle to edit_errtol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_errtol as text
%        str2double(get(hObject,'String')) returns contents of edit_errtol as a double
val=get(hObject, 'String');
[val, status]=str2num(val);
if ~status
    msgbox('Entry must be numeric', 'Incorrect data format');
    set(hObject, 'String', '0.1');
    return;
end

% --- Executes during object creation, after setting all properties.
function edit_errtol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_errtol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ecatol_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ecatol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ecatol as text
%        str2double(get(hObject,'String')) returns contents of edit_ecatol as a double
val=get(hObject, 'String');
[val, status]=str2num(val);
if ~status
    msgbox('Entry must be numeric', 'Incorrect data format');
    set(hObject, 'String', '0.3');
    return;
end

% --- Executes during object creation, after setting all properties.
function edit_ecatol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ecatol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_plotsensit.
function checkbox_plotsensit_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_plotsensit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_plotsensit



% --- Executes during object creation, after setting all properties.
function uitable_initmodel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitable_initmodel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
initsound_id=[1;1;1;1;1];
initlayer_num=[1;2;3;4;5];
initmodelcond=log10(1/100)*ones(5,1);
initmodelthick=[1;2;3;4;5];
constraint_wgt=[0;0;0;0;0];

data=[initsound_id initlayer_num initmodelthick initmodelcond constraint_wgt];
set(hObject, 'Data', data);


% --- Executes on button press in pushbutton_getInitModelFile.
function pushbutton_getInitModelFile_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_getInitModelFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(handles.text_dataFile, 'String'), 'Default data')
    msgbox('Enter data input first before creating layered model', 'Warning');
    return;
end

[f_modelc,pname,fIndex] = uigetfile('*.*', 'Select initial model file');

if fIndex~=0
    fid=fopen(strcat(pname, f_modelc), 'r');
    f = fread(fid,Inf);
    s = char(f');  % read in the file into a single character string
    [model_parms, status] = str2num(s);
    if status==0
        msgbox('Your model parameter file selection contains strings.  Select a different file.', 'Incorrect file format');
        return;
    end
    fclose(fid);
    %load data into uitable_initmodel
    %[1] , [2], and [3]  are x, y and z coordinates
    %[3]=sounding_id
    %[4]=layer_num
    %[5]=layer_thickness
    %[6]=Eca (S/m)
    %[7]=constraint_wgt
    
    sound_id=model_parms(:, 4);
    layer_num=model_parms(:, 5);
    layer_thick=model_parms(:, 6);
    Eca=model_parms(:, 7);
    constraint_wgt=model_parms(:, 8);
    
    
    %make sure number of model soundings matches number of data input
    %soundings
    data=get(handles.uitable_datainput, 'Data');
    datasound_id=unique(data(:, 4));
    num_soundings=length(datasound_id);
    
    datamodel_id=unique(sound_id);
    num_model=length(datamodel_id);
    
    if num_soundings~=num_model
        str={'Number of soundings in data input'; 'does not match'; 'number of soundings in initial model.'; 'Select a different file.'};
        warndlg(str);
        return;
    end
    
    data=[sound_id layer_num layer_thick Eca constraint_wgt];
    set(handles.uitable_initmodel, 'Data', data);   
   
    show=model_parms(model_parms(:, 4)==1, 6:7);
    set(handles.uitable_initmodel_show, 'Data', show);
    set(handles.edit_soundingid, 'String', 1);
    
    %get coordinates
    coords=get(handles.uitable_coords, 'Data');
    x=coords(1, 2);
    y=coords(1,3);
    z=coords(1,4);
    str={['Coordinates'], ['(', num2str(x), ', ', num2str(y), ', ' num2str(z), ')']};
    set(handles.text_coords, 'String', str);
    
    %set currently selected row and column to zero
    set(handles.text_currentRow, 'String', 0);
    set(handles.text_currentColumn, 'String', 0);

end

% --- Executes on button press in pushbutton_resetModel.
function pushbutton_resetModel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_resetModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
initsound_id=[1;1;1;1;1];
initlayer_num=[1;2;3;4;5];
initmodelcond=log10(1/100)*ones(5,1);
initmodelthick=[1;2;3;4;5];
constraint_wgt=[0;0;0;0;0];

data=[initsound_id initlayer_num initmodelthick initmodelcond constraint_wgt];
set(handles.uitable_initmodel, 'Data', data);

data=[initmodelcond initmodelthick];
set(handles.uitable_initmodel_show, 'Data', data);
set(handles.edit_soundingid, 'String', 1);

%modify coordinates
str={['Coordinates'], ['(0,0)']};
set(handles.text_coords, 'String', str);

%set currently selected row and column to zero
set(handles.text_currentRow, 'String', 0);
set(handles.text_currentColumn, 'String', 0);

%default data [x, y, z, sound_id, freq, STD error, VDM/HDM flag, Eca]
x=zeros(7, 1);
y=zeros(7, 1);
z=zeros(7, 1);
sound_id = ones(7, 1);
freq=[1530;8250;20370;33030;47970;63090;93090];
STD_error=ones(7,1);
flag=ones(7,1);
Eca=[56.2430; 195.2661;  334.0948;  428.5091;  511.9754;  580.6218;  693.2474];

data=[x y z sound_id freq STD_error flag Eca];
set(handles.uitable_datainput, 'Data', data);
set(handles.text_dataFile,'String', 'Default data');

%reset offset and freq data
freq=[1530;8250;20370;33030;47970;63090;93090];
offsets=1.66*ones(7,1);
data=[freq offsets];
set(handles.uitable_offset, 'Data', data);
set(handles.text_offsetFile,'String', 'Default Parameters');

set(handles.uipanel_constraint, 'Visible', 'off');
set(handles.radiobutton_1D, 'Value', 1);
set(handles.edit_sx, 'Visible', 'off');
set(handles.edit_sy, 'Visible', 'off');
set(handles.percent,'Visible','off')
set(handles.text_smooth, 'String', 'Smoothness (z)');

set(handles.text_modelparamFile, 'String', 'Default data');
set(handles.edit_minparam, 'String', '0.0');
set(handles.edit_maxparam, 'String', '1000.0');
set(handles.edit_initbarrier, 'String', '2.0');
set(handles.edit_maxiter, 'String', '10');
set(handles.edit_errtol, 'String', '0.1');
set(handles.edit_ecatol, 'String', '0.3');
set(handles.edit_horizreg, 'String', '10e-6');

muv={'10';'15'};
set(handles.listbox_vertreg, 'String', muv);

set(handles.edit_sx, 'String', '1');
set(handles.edit_sx, 'Visible', 'off');

set(handles.edit_sy, 'String', '1');
set(handles.edit_sy, 'Visible', 'off');

set(handles.edit_sz, 'String', '1');
set(handles.edit_sz, 'Visible', 'on');

set(handles.checkbox_plotlcurve, 'Visible', 'off')
set(handles.checkbox_plotsensit, 'Visible', 'off');
set(handles.uipanel11, 'Visible', 'off');


msgbox('All values have been reset to default', 'Data reset');
guidata(hObject, handles);

%CREATE LAYERED MODEL
% --- Executes on button press in pushbutton_createLayModel.
function pushbutton_createLayModel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_createLayModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(handles.text_dataFile, 'String'), 'Default data')
    msgbox('Enter data input first before creating layered model', 'Warning');
    return;
end

set(handles.uipanel_constraint, 'Visible', 'off');

%get number of soundings from input data
data=get(handles.uitable_datainput, 'Data');
sound_id=unique(data(:, 4));
num_soundings=length(sound_id);

str= strcat('Data input has:  ', num2str(num_soundings), ' soundings.  Enter how many layers in each sounding:');
prompt = {str, 'Enter layer thicknesses separated by a space:', 'Enter layer conductivities log10(S/m) separated by a space'};
dlg_title = 'Create Layered Model';
num_lines = 1;
def = {'5', '2  3  4  4  2', '10  5  7  5  10'};

TextInfo.FontSize = 12;

layer_model = inputdlg(prompt,dlg_title,num_lines,def);

if ~isempty(layer_model)
    %number of layers
    [num_layers, status]=str2num(layer_model{1});
    if ~status
        msgbox('All entries must only contain numeric values', 'Incorrect data format');
        return;
    end
        
    %thicknesses:  layer_model{2}
    [thickness, status]=str2num(layer_model{2});
    num_thick=length(thickness);
    if ~status
        msgbox('All entries must only contain numeric values', 'Incorrect data format');
        return;
    end
    
    %model cond: layer_model{3}
    [model_cond, status]=str2num(layer_model{3});
    if ~status
        msgbox('All entries must only contain numeric values', 'Incorrect data format');
        return;
    end
    
    %validation of user entry
    if length(model_cond)~=length(thickness)
        msgbox('Number of model conductivities must equal number of layer thicknesses.');
        return;
    else
        if length(model_cond)~=num_layers
            num_layers=length(model_cond);
        end
    end
    
    %make row for each layers in each sounding
    sound_id=ones(num_layers, 1);

    constraint_wgt=zeros(num_layers,1);
    layer_num=[1:1:num_layers];
    data=[sound_id, layer_num', thickness', model_cond', constraint_wgt];
    
    for ii=1:num_soundings-1
        add_data=[sound_id+ii, layer_num', thickness', model_cond', constraint_wgt];
        data=[data; add_data];
    end
    
   
    set(handles.uitable_initmodel, 'Data', data)
    
    %what to display - sounding number - uitable
    set(handles.edit_soundingid, 'String', 1);
    
    data_show=data((data(:,1)==1), 3:4);
    set(handles.uitable_initmodel_show, 'Data', data_show);
    
   %modify coordinates
    data_coords=get(handles.uitable_coords, 'Data');
    [coords]=data_coords(data_coords(:, 1)==1, 1:3);
    x=coords(1);
    y=coords(2);
    %z=coords(3);
    str={['Coordinates'], ['(' num2str(x) ', ' num2str(y), ')']};
    %str={['Coordinates'], ['(0, 0)']};
    set(handles.text_coords, 'String', str);
    
    %set currently selected row and column to zero
    set(handles.text_currentRow, 'String', 0);
    set(handles.text_currentColumn, 'String', 0);

    guidata(hObject, handles);
end


function edit_soundingid_Callback(hObject, ~, handles)
% hObject    handle to edit_soundingid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_soundingid as text
%        str2double(get(hObject,'String')) returns contents of edit_soundingid as a double
sound_id=get(hObject, 'String');
[sound_id, status]=str2num(sound_id);
if status
    data=get(handles.uitable_initmodel, 'Data');
    min_id=1;
    max_id=max(data(:,1));
    
    if sound_id >= min_id && sound_id <= max_id
        data_show=data((data(:,1)==sound_id), 3:4);
        set(handles.uitable_initmodel_show, 'Data', data_show);
        
        %modify coordinates
        data_coords=get(handles.uitable_coords, 'Data');
        [coords]=data_coords(data_coords(:, 1)==sound_id, 2:4);
        x=coords(1);
        y=coords(2);
        z=coords(3);
        str={['Coordinates'], ['(', num2str(x), ', ', num2str(y), ')']};
        set(handles.text_coords, 'String', str);

    else
        msgbox('This value is out of range', 'Out of range');
        set(hObject, 'String', 1);
        data_show=data((data(:,1)==1), 2:3);
        set(handles.uitable_initmodel_show, 'Data', data_show);    
    end
    button_val=unique(data((data(:,1)==sound_id), 4));
else
    msgbox('A numeric value is required for this field', 'Incorrect data format');
    set(hObject, 'String', 1);
    data=get(handles.uitable_initmodel, 'Data');
    button_val=unique(data((data(:,1)==1), 4));
    
    data_show=data((data(:,1)==1), 2:3);
    set(handles.uitable_initmodel_show, 'Data', data_show);
end
guidata(hObject, handles);    

% --- Executes during object creation, after setting all properties.
function edit_soundingid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_soundingid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_goUp.
function pushbutton_goUp_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_goUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
sound_id=str2num(get(handles.edit_soundingid, 'String'));
sound_id=sound_id+1;

data=get(handles.uitable_initmodel, 'Data');

min_id=1;
max_id=max(data(:,1));

if sound_id >= min_id && sound_id <= max_id
    data_show=data((data(:,1)==sound_id), 3:4);
    set(handles.uitable_initmodel_show, 'Data', data_show);
    
    %modify coordinates
    data_coords=get(handles.uitable_coords, 'Data');
    [coords]=data_coords(data_coords(:, 1)==sound_id, 2:4);
    x=coords(1);
    y=coords(2);
    z=coords(3);
    str={['Coordinates'], ['(', num2str(x), ', ', num2str(y), ', ' num2str(z), ')']};
    set(handles.text_coords, 'String', str);
    
    set(handles.edit_soundingid, 'String', sound_id);
    
    %set currently selected row and column to zero
    set(handles.text_currentRow, 'String', 0);
    set(handles.text_currentColumn, 'String', 0);
    set(handles.uipanel_constraint, 'Visible', 'off');
    
    guidata(hObject, handles);
end


% --- Executes on button press in pushbutton_delete.
function pushbutton_delete_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data=get(handles.uitable_initmodel, 'Data');
sound_id=str2num(get(handles.edit_soundingid, 'String'));
ind_sounding=find(data(:, 1)==sound_id); %entire sounding rows
tot_layers=length(ind_sounding);

if tot_layers==2
    warndlg('Sounding must have at least two layers.');
    return;
end

prompt = {'Enter layer to delete:'};
dlg_title = 'Delete Layer';
num_lines = 1;
def = {'1'};
layer_num= inputdlg(prompt,dlg_title,num_lines,def);
 
if isempty(layer_num)
   return; %user hit 'cancel' 
end
[layer_num, status]=str2num(layer_num{1});

if ~status
    msgbox('Entry must be numeric.', 'Incorrrect data format');
    return
end

ind=(ind_sounding(1)-1)+layer_num; % delete row

%save data_before and data_after
if ind==1;
    data_before=[];
    
    data(ind+1:layer_num, 2)=data(ind+1:layer_num, 2)-1;%change layer_num
    data_after=data(ind+1:end, :);
elseif ind==tot_layers
    data_before=data(1:ind-1, :);
    data_after=data(ind+1:end, :);
else
    data_before=data(1:ind-1, :);
    data(ind+1:layer_num, 2)=data(ind+1:layer_num, 2)-1;%change layer_num
    data_after=data(ind+1:end, :);
end

data=[data_before; data_after];
set(handles.uitable_initmodel, 'Data', data);

%refresh screen
if ~isempty(data)
    data_show=data((data(:,1)==sound_id), 3:4);
else
    data_show=[];
end
set(handles.uitable_initmodel_show, 'Data', data_show);

guidata(hObject, handles);

% --- Executes on button press in pushbutton_goDown.
function pushbutton_goDown_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_goDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
sound_id=str2num(get(handles.edit_soundingid, 'String'));
sound_id=sound_id-1;

data=get(handles.uitable_initmodel, 'Data');

min_id=1;
max_id=max(data(:,1));

if sound_id >= min_id && sound_id <= max_id
    data_show=data((data(:,1)==sound_id), 3:4);
    set(handles.uitable_initmodel_show, 'Data', data_show);

    %modify coordinates
    data_coords=get(handles.uitable_coords, 'Data');
    [coords]=data_coords(data_coords(:, 1)==sound_id, 2:4);
    x=coords(1);
    y=coords(2);
    z=coords(3);
    str={['Coordinates'], ['(', num2str(x), ', ', num2str(y), ', ' num2str(z), ')']};
    set(handles.text_coords, 'String', str);
    
    set(handles.edit_soundingid, 'String', sound_id);
    
    %set currently selected row and column to zero
    set(handles.text_currentRow, 'String', 0);
    set(handles.text_currentColumn, 'String', 0);
end
set(handles.uipanel_constraint, 'Visible', 'off');
guidata(hObject, handles);

% --- Executes on button press in pushbutton_insert.
function pushbutton_insert_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_insert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prompt = {'Enter layer number position to insert:', 'Enter layer thickness:', ...
    'Enter layer conductivity (S/m):', 'Enter constraint weight (0 to 1):'};
dlg_title = 'Insert Layer';
num_lines = 1;
def = {'1', '', '', '1'};
insert_layer = inputdlg(prompt,dlg_title,num_lines,def);

data=get(handles.uitable_initmodel, 'Data');
sound_id=str2num(get(handles.edit_soundingid, 'String'));
sound_index_all=find(data(:, 1)==sound_id);
sound_index=sound_index_all(1);
  
data_sounding=data(data(:, 1)==sound_id, :);
num_layers=size(data_sounding,1);

if ~isempty(insert_layer)    
    [ins_num, status]=str2num(insert_layer{1});
    if ~status
       msgbox('All entries must be numeric.', 'Incorrect data format');
       return;
    elseif ins_num > num_layers + 1
        %check to make sure valid number was entered
        str={'Insertion layer number must be less than ', num2str(num_layers+1)};
        msgbox(str);
        return;
    end
    
    [thickness, status]=str2num(insert_layer{2});
    if ~status
       msgbox('All entries must be numeric.', 'Incorrect data format');
       return;
    end
        
    [model_cond, status]=str2num(insert_layer{3});
    if ~status
       msgbox('All entries must be numeric.', 'Incorrect data format');
       return;
    end
    
    if length(thickness)~=length(model_cond) || length(thickness)>1 || length(model_cond)>1
        msgbox('Must have one layer thickness and one model conductivity.');
        return;
    end
    
    [constraint_wgt, status]=str2num(insert_layer{4});
    if ~status
       msgbox('All entries must be numeric.', 'Incorrect data format');
       return;
    elseif constraint_wgt > 1
        msgbox('Constraint weight must be between 0 and 1');
        return;
    end
    
    
    data_add=[sound_id, ins_num, thickness, model_cond, constraint_wgt];
    row_selected=(sound_index-1)+ins_num;
    
    if sound_index==1 && row_selected==1 %add to beginning
        data_before=[];
    else %insert
        data_before=data(1:row_selected-1, :);   
    end
    
    if ins_num==size(data, 1) %add to end
        data_after=[];
    else %insert
        data_after=data(row_selected:end, :);   
        data_after(1:num_layers-ins_num+1, 2)=data_after(1:num_layers-ins_num+1, 2)+1; %increment layer_num
    end
    
    data=[data_before; data_add; data_after]; 
     
    set(handles.uitable_initmodel, 'Data', data);
    
    data_show=data(data(:,1)==sound_id, 3:4);
    set(handles.uitable_initmodel_show, 'Data', data_show);
    
    if constraint_wgt > 0
        set(handles.uipanel_constraint, 'Visible', 'on');
        set(handles.edit_constraintWgt, 'String', constraint_wgt);
        set(handles.slider_constraintWgt, 'Value', constraint_wgt);
    end
    
    set(handles.text_currentRow, 'String', row_selected); %save current row
        
end
guidata(hObject, handles);

% --- Executes on selection change in listbox_vertreg.
function listbox_vertreg_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_vertreg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_vertreg contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_vertreg
if strcmp(get(handles.figure1, 'SelectionType'), 'open') %user hit enter or double-clicked
    prompt = {'Enter vertical regularization coefficient:'};
    dlg_title = 'Edit vertical regularization';
    num_lines = 1;
    
    index_selected = get(hObject,'Value');
    list = get(hObject,'String');
    item_selected = list{index_selected};
    
    def = {item_selected };
    vert_reg = inputdlg(prompt,dlg_title,num_lines,def);

    if ~isempty(vert_reg)
        list(index_selected)=vert_reg;
        set(hObject, 'String', list);
    end
end
guidata(hObject, handles);  

% --- Executes during object creation, after setting all properties.
function listbox_vertreg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_vertreg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on selection change in popupmenu_regcriteria.
function popupmenu_regcriteria_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_regcriteria (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_regcriteria contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_regcriteria
val=get(hObject,'Value');
dim=get(handles.radiobutton_1D, 'Value');
if val==4 
    if dim==1
        set(handles.checkbox_plotlcurve, 'Visible', 'on')
        set(handles.percent,'Visible','off')
    else
        set(handles.checkbox_plotlcurve, 'Visible', 'off')
        set(handles.checkbox_DOI, 'Visible', 'off')
        set(handles.percent,'Visible','on')

    end
    set(handles.checkbox_plotsensit, 'Visible', 'off');
    set(handles.checkbox_DOI, 'Visible', 'off');
    set(handles.percent,'Visible','off')
    set(handles.uipanel11, 'Visible', 'on');
else
    set(handles.checkbox_plotlcurve, 'Visible', 'off')
    set(handles.checkbox_plotsensit, 'Visible', 'on');
    set(handles.checkbox_DOI, 'Visible', 'on');
    set(handles.uipanel11, 'Visible', 'on');
    set(handles.percent,'Visible','on')
end

    
% --- Executes during object creation, after setting all properties.
function popupmenu_regcriteria_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_regcriteria (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: text_dt controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when entered data in editable cell(s) in uitable_initmodel_show.
function uitable_initmodel_show_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable_initmodel_show (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
str=eventdata.EditData;
[val, status] = str2num(str);
if (status)
    %user can enter constraint weight
    prompt = {'Enter constraint weight of layer (value from 0 to 1):'};
    dlg_title = 'Constraint Weight';
    num_lines = 1;
    def = {'1'};
    val1 = inputdlg(prompt,dlg_title,num_lines,def);

    if isempty(val1)
        constraint_wgt=0;
    else
        [constraint_wgt, status]=str2num(val1{1});
       if constraint_wgt > 1    
            msgbox('Value must be between 0 and 1.  Defaulting to 1.', 'Out of range');
            set(handles.edit_constraintWgt, 'String', '1');
            set(handles.slider_constraintWgt, 'Value', 1);
       else
            set(handles.slider_constraintWgt, 'Value', constraint_wgt);
            set(handles.edit_constraintWgt, 'String', constraint_wgt);
       end
    end
              
    sound_id=str2num(get(handles.edit_soundingid, 'String'));
     
    data=get(handles.uitable_initmodel, 'Data');
    ind_edit=min(find(data(:, 1)==sound_id));
    
    row=(ind_edit-1) + eventdata.Indices(1);
    set(handles.text_currentRow, 'String', row); %save current row
    
    col=eventdata.Indices(2)+2; %two more columns in model_uitable 
    set(handles.text_currentColumn, 'String', col); %save current column
    
    data(row, col)=val;
    data(row, 5)=constraint_wgt;
    set(handles.uitable_initmodel, 'Data', data);
    
    %set slider and constraint_wgt
    set(handles.uipanel_constraint, 'Visible', 'on');
else
    msgbox('A numeric value is required for this field', 'Incorrect data format');
    %eventdata.EditData=eventdata.PreviousData;    
    data = get(handles.uitable_initmodel_show,'Data');
    data(eventdata.Indices(1),eventdata.Indices(2)) = eventdata.PreviousData; 
    set(handles.uitable_initmodel_show,'Data',data);    
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function uitable_initmodel_show_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitable_initmodel_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
initmodelthick=[2;2;2;2;2];
initmodelcond=log10(1/100)*ones(5,1);

data=[initmodelthick initmodelcond];
set(hObject, 'Data', data);
set(hObject, 'ColumnEditable', [true, true]);
guidata(hObject, handles);

% --- Executes on button press in pushbutton_goFirst.
function pushbutton_goFirst_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_goFirst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=get(handles.uitable_initmodel, 'Data');
data_show=data((data(:,1)==1), 3:4);
set(handles.uitable_initmodel_show, 'Data', data_show);
set(handles.edit_soundingid, 'String', 1);
button_val=unique(data((data(:,1)==1), 4));

%modify coordinates
data_coords=get(handles.uitable_coords, 'Data');
[coords]=data_coords(data_coords(:, 1)==1, 2:4);
x=coords(1);
y=coords(2);
z=coords(3);
str={['Coordinates'], ['(', num2str(x), ', ', num2str(y), ', ' num2str(z), ')']};
set(handles.text_coords, 'String', str);

%set currently selected row and column to zero
set(handles.text_currentRow, 'String', 0);
set(handles.text_currentColumn, 'String', 0);
guidata(hObject, handles);

% --- Executes on button press in pushbutton_goEnd.
function pushbutton_goEnd_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_goEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=get(handles.uitable_initmodel, 'Data');
max_id=max(data(:,1));
data_show=data((data(:,1)==max_id), 3:4);
set(handles.uitable_initmodel_show, 'Data', data_show);
set(handles.edit_soundingid, 'String', max_id);

%modify coordinates
data_coords=get(handles.uitable_coords, 'Data');
[coords]=data_coords(data_coords(:, 1)==max_id, 2:4);
x=coords(1);
y=coords(2);
z=coords(3);
str={['Coordinates'], ['(', num2str(x), ', ', num2str(y), ', ' num2str(z), ')']};
set(handles.text_coords, 'String', str);

%set currently selected row and column to zero
set(handles.text_currentRow, 'String', 0);
set(handles.text_currentColumn, 'String', 0);

set(handles.uipanel_constraint, 'Visible', 'off');
guidata(hObject, handles);


% --- Executes on button press in pushbutton_addvertreg.
function pushbutton_addvertreg_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_addvertreg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prompt = {'Enter vertical regularization coefficient:'};
dlg_title = 'Add vertical regularization';
num_lines = 1;

def = {''};
vert_reg = inputdlg(prompt,dlg_title,num_lines,def);
 
if ~isempty(vert_reg)
    vert_reg=cell2mat(vert_reg);
    [val, status]=str2num(vert_reg);
    if status 
        list = get(handles.listbox_vertreg,'String');
        cnt=length(list);
        list{cnt+1}=vert_reg;
        set(handles.listbox_vertreg, 'String', list);
    else
        msgbox('You must enter a numeric value.', 'Incorrect data format');
    end
end
guidata(hObject, handles);

% --- Executes on button press in pushbutton_insvert_reg.
function pushbutton_insvert_reg_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_insvert_reg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prompt = {'Enter vertical regularization coefficient:'};
dlg_title = 'Insert vertical regularization';
num_lines = 1;

def = {''};
vert_reg = inputdlg(prompt,dlg_title,num_lines,def);
 
if ~isempty(vert_reg)
    vert_reg=cell2mat(vert_reg);
    [val, status]=str2num(vert_reg);
    if status 
        index_selected = get(handles.listbox_vertreg,'Value');
        list = get(handles.listbox_vertreg,'String');
        cnt=length(list);
        
        add=vert_reg;
        if index_selected==1
            list=[add; list];
        else
            list_before=list(1:index_selected-1);
            list_after=list(index_selected: cnt);
            list=[list_before; add; list_after];
        end
        set(handles.listbox_vertreg, 'String', list);
    else
        msgbox('You must enter a numeric value.', 'Incorrect data format');
    end
end
guidata(hObject, handles);

% --- Executes on button press in pushbutton_delvertreg.
function pushbutton_delvertreg_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_delvertreg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
index_selected = get(handles.listbox_vertreg,'Value');
list = get(handles.listbox_vertreg,'String');
cnt=length(list);

if cnt==1
    msgbox('You need at least one vertical regularization parameter.', 'Error');
    return;
end

if index_selected==1
    list=list(2:cnt);
elseif index_selected==cnt
    %list=list(1:cnt-1);
    list(cnt)=[];
else
    list_before=list(1:index_selected-1);
    list_after=list(index_selected+1: cnt);
    list=[list_before; list_after];
end

set(handles.listbox_vertreg, 'String', list, 'Value', 1);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function pushbutton_delvertreg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_delvertreg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on slider movement.
function slider_constraintWgt_Callback(hObject, eventdata, handles)
% hObject    handle to slider_constraintWgt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

val=get(hObject,'Value');
round_val= (round(100.*val))/100;
set(handles.edit_constraintWgt, 'String', round_val);

%save new entry to data table
current_row=str2num(get(handles.text_currentRow, 'String'));
data=get(handles.uitable_initmodel, 'Data');
data(current_row, 5)=val;
set(handles.uitable_initmodel, 'Data', data);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider_constraintWgt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_constraintWgt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_constraintWgt_Callback(hObject, eventdata, handles)
% hObject    handle to edit_constraintWgt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_constraintWgt as text
%        str2double(get(hObject,'String')) returns contents of edit_constraintWgt as a double
val=get(hObject, 'String');
[val, status]=str2num(val);
if status
    if val > 1
        msgbox('Value must be between 0 and 1', 'Out of range');
        set(hObject, 'String', 0);
        set(handles.slider_constraintWgt, 'Value', 0);
        return;
    end
    set(handles.slider_constraintWgt, 'Value', val);
    
    %save new entry to data table
    current_row=str2num(get(handles.text_currentRow, 'String'));
    data=get(handles.uitable_initmodel, 'Data');
    data(current_row, 5)=val;
    set(handles.uitable_initmodel, 'Data', data);
    
else
    msgbox('A numeric value is required', 'Incorrect data format');
    set(hObject, 'String', 0);
    set(handles.slider_constraintWgt, 'Value', 0);
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_constraintWgt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_constraintWgt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected cell(s) is changed in uitable_initmodel_show.
function uitable_initmodel_show_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitable_initmodel_show (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)

if numel(eventdata.Indices)==0
    return;
end
r=eventdata.Indices(1);
c=eventdata.Indices(2);

sound_id=str2num(get(handles.edit_soundingid, 'String'));
all_data=get(handles.uitable_initmodel, 'Data');
data=all_data(all_data(:, 1)==sound_id, :);
 

ind_edit=min(find(all_data(:, 1)==sound_id));
row=(ind_edit-1) + r;
set(handles.text_currentRow, 'String', row); %save current row
set(handles.text_currentColumn, 'String', c); %save current column

constraint_wgt=data(r, 5);
if constraint_wgt > 0
    set(handles.uipanel_constraint, 'Visible', 'on');
    set(handles.slider_constraintWgt,'Value', constraint_wgt);
    round_val= (round(100.*constraint_wgt))/100;
    set(handles.edit_constraintWgt, 'String', round_val);
else
    set(handles.uipanel_constraint, 'Visible', 'off');
end 
guidata(hObject, handles);

% --- Executes on button press in pushbutton_saveInitModel.
function pushbutton_saveInitModel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_saveInitModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   
[fname,pname] = uiputfile({'*.dat';'*.txt';'*.*'}, 'Save initial inversion model');
%save modeling parameters
% [x, y, sound_id, layer_num, lyr_thick, Eca, constraint_wgt]
if fname==0
    return;
end
fid=fopen(strcat(pname, fname), 'w');

data=get(handles.uitable_initmodel, 'Data'); %[sound_id, layer_num thick, Eca, wgt]
model_parms=data;
coords=get(handles.uitable_coords, 'Data');
num_soundings=size(unique(data(:,1)),1);

%get coordinates
for ii=1:num_soundings
    ind_sounding=find(data(:,1)==ii);
    sound_coords=coords(coords(:,1)==ii, 2:4);
    
    for jj=1:size(ind_sounding, 1)
        model_parms(ind_sounding(jj), 6:8)=sound_coords;
    end
end
 
model_parms=[model_parms(:, 6:8), model_parms(:, 1:5)];
fprintf(fid, '%4.2f  %4.2f  %4.2f %d %d %6.4f  %6.4f  %6.4f\r\n',  model_parms');
fclose(fid);
msgbox('Initial model has been saved');


% --- Executes during object creation, after setting all properties.
function uitable_datainput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitable_datainput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

%default data [x, y, z, sound_id, freq, STD error, VDM/HDM flag, Eca]
x=zeros(7, 1);
y=zeros(7, 1);
z=zeros(7, 1);
sound_id = ones(7, 1);
freq=[1530;8250;20370;33030;47970;63090;93090];
STD_error=ones(7,1);
flag=ones(7,1);
Eca=[56.2430; 195.2661;  334.0948;  428.5091;  511.9754;  580.6218;  693.2474];

data=[x y z sound_id freq STD_error flag Eca];
set(hObject, 'Data', data);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function uitable_coords_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitable_coords (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

sound_id=1;
x=0;
y=0;
z=0;
data=[sound_id x y z];
set(hObject, 'Data', data);
guidata(hObject, handles);


function edit_sx_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sx as text
%        str2double(get(hObject,'String')) returns contents of edit_sx as a double
val=get(hObject, 'String');
[val, status]=str2num(val);
if ~status
    msgbox('Entry must be numeric', 'Incorrect data format');
    set(hObject, 'String', '1.0');
    return;
end

% --- Executes during object creation, after setting all properties.
function edit_sx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_sy_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sy as text
%        str2double(get(hObject,'String')) returns contents of edit_sy as a double
val=get(hObject, 'String');
[val, status]=str2num(val);
if ~status
    msgbox('Entry must be numeric', 'Incorrect data format');
    set(hObject, 'String', '1.0');
    return;
end

% --- Executes during object creation, after setting all properties.
function edit_sy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in uipanel_dimension.
function uipanel_dimension_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_dimension 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read onl
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
str=get(eventdata.NewValue, 'String');

if strcmp(str, '1D')
    set(handles.text_smooth, 'String', 'Smoothness (z)');
    set(handles.edit_sx, 'Visible', 'off');
    set(handles.edit_sy, 'Visible', 'off');
    set(handles.percent,'Visible','off')
    invType=get(handles.popupmenu_regcriteria, 'Value');
    if invType==4
        set(handles.uipanel11, 'Visible', 'on');
        set(handles.checkbox_plotlcurve, 'Visible', 'on');
            else
        set(handles.uipanel11, 'Visible', 'off');
        set(handles.checkbox_plotlcurve, 'Visible', 'off');
    end 
    set(handles.checkbox_plotsensit, 'Visible', 'off');
elseif strcmp(str,'2D') 
    set(handles.percent,'Visible','on')
    set(handles.text_smooth, 'String', 'Smoothness (x, z)');
    set(handles.edit_sx, 'Visible', 'on'); 
    set(handles.edit_sy, 'Visible', 'off'); 
    set(handles.checkbox_plotlcurve, 'Visible', 'off');
    set(handles.checkbox_plotsensit, 'Visible', 'on');
    set(handles.uipanel11, 'Visible', 'on');
elseif strcmp(str,'3D')
    set(handles.percent,'Visible','on')
    set(handles.text_smooth, 'String', 'Smoothness (x, y, z)');
    set(handles.edit_sx, 'Visible', 'on'); 
    set(handles.edit_sy, 'Visible', 'on'); 
    set(handles.uipanel11, 'Visible', 'on');
    set(handles.checkbox_plotlcurve, 'Visible', 'off');
    set(handles.checkbox_plotsensit, 'Visible', 'on');
end
        


% --- Executes during object creation, after setting all properties.
function uipanel_dimension_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel_dimension (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function edit_sz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in uibuttongroup1.
function uibuttongroup1_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup1 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_constraintwgt.
function pushbutton_constraintwgt_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_constraintwgt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prompt = {'Enter layer number to modify constraint weight:',...
    'Enter constraint weight (0 to 1):'};
dlg_title = 'Modify/Edit Constraint Weight';
num_lines = 1;
def = {'', ''};
wgt_layer = inputdlg(prompt,dlg_title,num_lines,def);

data=get(handles.uitable_initmodel, 'Data');
sound_id=str2num(get(handles.edit_soundingid, 'String'));

sound_index_all=find(data(:, 1)==sound_id);
sound_index=sound_index_all(1);

data_sounding=data(data(:, 1)==sound_id, :);
ind_sounding=find(data(:, 1)==sound_id); %entire sounding rows
num_layers=size(data_sounding,1);

if ~isempty(wgt_layer)    
    [layer_num, status]=str2num(wgt_layer{1});
    if ~status
       msgbox('All entries must be numeric.', 'Incorrect data format');
       return;
    elseif layer_num > num_layers
        %check to make sure valid number was entered
        str={'Layer number must be less than ', num2str(num_layers)};
        msgbox(str);
        return;
    end
    
    [constraint_wgt, status]=str2num(wgt_layer{2});
    if ~status
       msgbox('All entries must be numeric.', 'Incorrect data format');
       return;
    elseif constraint_wgt > 1
        msgbox('Constraint weight must be between 0 and 1');
        return;
    end
    
    %modify constraint wgt
    ind=(sound_index-1)+layer_num;
    data(ind, 5)=constraint_wgt;
    set(handles.uitable_initmodel, 'Data', data);
    set(handles.text_currentRow, 'String', ind); %save current row

    data_show=data(data(:,1)==sound_id, 3:4);
    set(handles.uitable_initmodel_show, 'Data', data_show);
    
    if constraint_wgt > 0
        set(handles.uipanel_constraint, 'Visible', 'on');
        set(handles.edit_constraintWgt, 'String', constraint_wgt);
        set(handles.slider_constraintWgt, 'Value', constraint_wgt);
    end

end
guidata(hObject, handles);



function edit_sz_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sz as text
%        str2double(get(hObject,'String')) returns contents of edit_sz as a double
val=get(hObject, 'String');
[val, status]=str2num(val);
if ~status
    msgbox('Entry must be numeric', 'Incorrect data format');
    set(hObject, 'String', '1.0');
    return;
end

function parallel_Callback(hObject, eventdata, handles)
% hObject    handle to parallel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of parallel as text
%        str2double(get(hObject,'String')) returns contents of parallel as a double
vall=str2num(get(handles.parallel,'String'));
        b=findResource();
        b=b.ClusterSize;
    if vall>b
        msgbox(['Error: maximum allowed number of parallel processors is ', num2str(b)]);return;end 


% --- Executes during object creation, after setting all properties.
function parallel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to parallel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_DOI.
function checkbox_DOI_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_DOI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_DOI


% --- Executes on button press in lmodel.
% --- Executes on button press in modeel.
function lmodel_Callback(hObject, eventdata, handles)
% hObject    handle to modeel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[f_model,pname,fIndex] = uigetfile('*.*', 'Select model');
if fIndex~=0   
    fid=fopen(strcat(pname, f_model), 'r');
    f = fread(fid,Inf);
    s = char(f');  % read in the file into a single character string
    [data, status] = str2num(s);
    if status==0
        msgbox('Your model contains strings.  Select a different file.', 'Incorrect file format');
        return;
    end
    
   set(handles.uitable_model, 'Data', data);
    set(handles.modeel, 'String', f_model);
    
end
function pushbutton64_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton64 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h=figure('WindowStyle', 'modal', 'Name', 'View conductivity model');

%create table viewing parameters
data=get(handles.uitable_model, 'Data');
s=size(data, 2);
isEditable=zeros(1,s);
isEditable = logical(isEditable);

t_data=uitable(h, 'Tag', 'uitable_model', 'Position', [30 80 500 300],...
    'Data', get(handles.uitable_model, 'Data'), ...
    'FontSize', 10,...
    'ColumnEditable', isEditable,...
    'ColumnName', {'X (m)', 'Y (m)', 'Depth (m)', 'Conductivity (S/m)'});  

cancel_but=uicontrol(h, 'Style', 'pushbutton', 'String', 'Cancel',...
    'FontSize', 10, 'Position', [240 30 70 30],...
    'Callback', 'close');


%CALIBRATION
% --- Executes on button press in pushbutton65
function pushbutton65_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton65 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
initmodel_data=get(handles.uitable_initmodel, 'Data');de=initmodel_data(:,3);
data1=get(handles.uitable_datainput, 'Data');
data2=get(handles.uitable_model, 'Data');
%get depth
rr=get(handles.uitable_offset, 'Data');r=rr(:,2);f=rr(:,1);%str2double(get(handles.coil,'String'));
%get q
fff=data1(:,7);q=fff;vall2=str2num(get(handles.parallel,'String'));
%[rawGem calGem filtGem]=FEMIC_cal(data2,data1,f,r,q,vall2);filtGem=filtGem';filtGem=filtGem(:);
[rawGem calGem filtGem]=FEMIC_cal(data2,data1,f,r,q(1),vall2);
filtGem=filtGem;filtGem=filtGem(:);
%replace original data with filtered data
%dlmwrite('Filtered_data.dat',filtGem);

data=get(handles.uitable_datainput, 'Data');data(:,8)=filtGem;
set(handles.uitable_datainput, 'Data', data);



function percent_Callback(hObject, eventdata, handles)
% hObject    handle to percent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of percent as text
%        str2double(get(hObject,'String')) returns contents of percent as a double
perc=str2num(get(handles.percent,'String'));


% --- Executes during object creation, after setting all properties.
function percent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to percent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_5_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
open('FEMIC_Manual.pdf')


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_6_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox({'                FEMIC version 1.0, 2012'; '   '; ' Developed by Rutgers Newark and USGS'; ''},'About FEMIC')



function edit_elevation_Callback(hObject, eventdata, handles)
% hObject    handle to edit_elevation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val=get(hObject, 'String');
[val, status]=str2num(val);
if ~status
    msgbox('Entry must be numeric', 'Incorrect data format');
    set(hObject, 'String', '-1.0');
    return;
end

% Hints: get(hObject,'String') returns contents of edit_elevation as text
%        str2double(get(hObject,'String')) returns contents of edit_elevation as a double


% --- Executes during object creation, after setting all properties.
function edit_elevation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_elevation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
