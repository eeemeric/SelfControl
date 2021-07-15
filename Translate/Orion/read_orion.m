function varargout = read_orion(varargin)
%READ_ORION reads behavioral data from the Orion ACCDB (*.mdb;*.accdb;)
%   TRIALBLOB=READ_ORION
%       If no argument is provided, this function opens a file dialog to
%       get the path to an ACCDB file and returns trialblobs.
%   VARARGOUT=READ_ORION(VARARGIN)
%       When you have multiple blob types to read, you can give the list of
%       them as VARARGIN. The possible types are 'event','trial','eye' and
%       so on (see the code for the complete list). If one of the argument
%       indicates a existing file, it is considered a path to the ACCDB.
%       The return variables in VARARGOUT will be assigned as the order of
%       the blob types in VARARGIN. For example,
%       [TRIAL,EYE]=READ_ORION('trial','eye');
%       [TRIAL,EYE]=READ_ORION(ACCDBPATH,'trial','eye');
%       [TRIAL,EYE]=READ_ORION('trial','eye',ACCDBPATH);
%
%   Apr 16, 2008        written by Jaewon Hwang (jaewon.hwang@gmail.com)
%   Jul 20, 2008        added two more options, 'ver' & 'len'
%   May 20, 2009        added ACCDB support (MS ACCESS 2007 format)
%   May 31, 2009        added an option, 'info' which returns the neural
%                       recording file name and the session note when
%                       available
%   May 30, 2013        added an option for reading ArmBlob, 'arm'
%   Nov 21, 2013        added an option for reading Params, 'param'

% parse input arguments
dbfield = [];
for m=1:nargin
    if ischar(varargin{m})
        switch varargin{m}
            case 'event'
                dbfield(end+1) = 1; %#ok<AGROW>
            case 'trial'
                dbfield(end+1) = 2; %#ok<AGROW>
            case 'eye'
                dbfield(end+1) = 3; %#ok<AGROW>
            case 'joy'
                dbfield(end+1) = 4; %#ok<AGROW>
            case 'arm'
                dbfield(end+1) = 5; %#ok<AGROW>
            case 'ver'
                dbfield(end+1) = 6; %#ok<AGROW>
            case 'len'
                dbfield(end+1) = 7; %#ok<AGROW>
            case 'info'
                dbfield(end+1) = 8; %#ok<AGROW>
            case 'param'
                dbfield(end+1) = 9; %#ok<AGROW>
            otherwise
                if exist(varargin{m},'file')
                    dbfile = varargin{m};
                else
                    error('read_orion:varargin','Unknown option or file');
                end
        end
    else
        error('read_orion:varargin','This function does not take any numeric option');
    end
end
% default query when there is no query
if isempty(dbfield), dbfield = 2; end

% get the dbfile path if it is not given
if ~exist('dbfile','var')
    [filename,path] = uigetfile( { '*.mdb;*.accdb;','Microsoft Access database (*.mdb;*.accdb)' } );
    cd(path);
    dbfile = [path,filename];
end
disp(dbfile)
% trial blob version correction
[~,n,e] = fileparts(dbfile);
filename = [n e];
switch lower(filename)
    case 'pb021709.mdb', trialblob_ver = 200;
    case 'pb030309.mdb', trialblob_ver = 401;
    otherwise, trialblob_ver = 0;
end

% connect the database (i.e., the mdb or accdb file saved by Orion)
conn = actxserver('adodb.connection');
rs  = actxserver('adodb.recordset');
% conn.Open(['Provider=Microsoft.Jet.OLEDB.4.0; Persist Security Info=False; Data Source=',dbfile]); % old connection string
conn.Open(['Provider=Microsoft.ACE.OLEDB.12.0; Persist Security Info=False; Data Source=',dbfile]);

% read the database table.
nfield = length(dbfield);
varargout = cell(1,nfield);
for m=1:nfield
    switch dbfield(m)
        case 1  % event
            query = 'SELECT EventName, Code FROM Events ORDER BY ID';
            rs.Open(query,conn,1,3,1);
            Events = rs.GetRows(-1)';
            rs.Close();
            query = 'SELECT TrialID, Stage, EventID, TimeOffset FROM TrialEvents ORDER BY ID';
            rs.Open(query,conn,1,3,1);
            TrialEvents = rs.GetRows(-1)';
            rs.Close();
            varargout{m}.codelist         = Events;
            varargout{m}.trialevents      = zeros(size(TrialEvents));
            varargout{m}.trialevents(:,1) = [TrialEvents{:,1}]';
            varargout{m}.trialevents(:,2) = [TrialEvents{:,2}]';
            varargout{m}.trialevents(:,3) = [Events{[TrialEvents{:,3}],2}]';
            varargout{m}.trialevents(:,4) = [TrialEvents{:,4}]';
        case 2  % trial
            query = 'SELECT TrialBlob FROM BlockTrials ORDER BY BlockID, ID';
            rs.Open(query,conn,1,3,1);
            BlockTrials = rs.GetRows(-1)';
            rs.Close();
            varargout{m} = parse_trialblob(BlockTrials,trialblob_ver);
        case 3  % eye
            query = 'SELECT EyeBlob FROM BlockTrials ORDER BY BlockID, ID';
            rs.Open(query,conn,1,3,1);
            BlockTrials = rs.GetRows(-1)';
            rs.Close();
            varargout{m} = parse_xyblob(BlockTrials);
        case {4,5}  % joy & arm
            query = 'SELECT JoyBlob FROM BlockTrials ORDER BY BlockID, ID';
            rs.Open(query,conn,1,3,1);
            BlockTrials = rs.GetRows(-1)';
            rs.Close();
            varargout{m} = parse_xyblob(BlockTrials);
        case 6  % ver
            if 0~=trialblob_ver
                varargout{m} = trialblob_ver;
            else
                query = 'SELECT TOP 1 TrialBlob FROM BlockTrials ORDER BY BlockID, ID';
                rs.Open(query,conn,1,3,1);
                BlockTrials = rs.GetRows(-1)';
                rs.Close();
                varargout{m} = hex2int(BlockTrials{1}(5:8));
            end
        case 7  % len
            query = 'SELECT TOP 1 TrialBlob FROM BlockTrials ORDER BY BlockID, ID';
            rs.Open(query,conn,1,3,1);
            BlockTrials = rs.GetRows(-1)';
            rs.Close();
            varargout{m} = hex2int(BlockTrials{1}(1:4));
        case 8  % info
            varargout{m}.filename = filename;
            varargout{m}.Type = 'Behavior';
            try
                query = 'SELECT TOP 1 PlexonFile, PlexonWindow, TDTServer, TDTTank, TDTBlock From SessionDetails ORDER BY SessionID';
                rs.Open(query,conn,1,3,1);
                if ~rs.EOF
                    SessionDetails = rs.GetRows(-1)';
                    if strcmp(SessionDetails{1},'')
                        varargout{m} = cell2struct(SessionDetails(3:5),{'TDTServer','TDTTank','TDTBlock'},2);
                        varargout{m}.Type = 'TDT';
                    else
                        varargout{m} = cell2struct(SessionDetails(1:2),{'PlexonFile','PlexonWindow'},2);
                        varargout{m}.Type = 'Plexon';
                    end
                end
                rs.Close();
            catch %#ok<CTCH>
                query = 'SELECT TOP 1 PlexonFile, PlexonWindow From SessionDetails ORDER BY SessionID';
                rs.Open(query,conn,1,3,1);
                if ~rs.EOF
                    SessionDetails = rs.GetRows(-1)';
                    varargout{m} = cell2struct(SessionDetails,{'PlexonFile','PlexonWindow'},2);
                    varargout{m}.Type = 'Plexon';
                end
                rs.Close();
            end
            query = 'SELECT TOP 1 Message From SessionMessages ORDER BY SessionID';
            rs.Open(query,conn,1,3,1);
            if ~rs.EOF, SessionMessages = rs.GetRows(-1)'; varargout{m}.Message = SessionMessages{1}; else varargout{m}.Message = ''; end
            rs.Close();
            if ~strcmp('Behavior',varargout{m}.Type)
                query = 'SELECT ID, ManifoldID From ManifoldChannels';
                rs.Open(query,conn,1,3,1);
                if ~rs.EOF, varargout{m}.Ch2Manifold = cell2mat(rs.GetRows(-1)'); end
                rs.Close();
                query = 'SELECT ChannelID, Depth, Unit0, Unit1, Unit2, Unit3 From SessionDepths ORDER BY SessionID';
                rs.Open(query,conn,1,3,1);
                if ~rs.EOF, varargout{m}.SessionDepths = cell2mat(rs.GetRows(-1)'); end
                rs.Close();
                query = 'SELECT ManifoldID, Number1, Number2, IsPolarCoords From SessionManifolds ORDER BY SessionID';
                rs.Open(query,conn,1,3,1);
                if ~rs.EOF
                    SessionManifolds = rs.GetRows(-1)';
                    varargout{m}.SessionManifolds = [double([SessionManifolds{:,1}]') cell2mat(SessionManifolds(:,2:3)) [SessionManifolds{:,4}]'];
                end
                rs.Close();
            end
        case 9  % param
            query = 'SELECT ID, ParamName FROM Params ORDER BY ID';
            rs.Open(query,conn,1,3,1);
            Params = rs.GetRows(-1)';
            rs.Close();
            query = 'SELECT ParameterID, FirstTrialID, Value FROM BlockValues ORDER BY ParameterID, FirstTrialID';
            rs.Open(query,conn,1,3,1);
            BlockValues = rs.GetRows(-1)';
            rs.Close();
            param = [true; 0<diff([BlockValues{:,1}]')];
            value = 0==strcmp(BlockValues(:,3),circshift(BlockValues(:,3),1));
            row = param | value;
            paramID = cond_find([Params{:,1}]',[BlockValues{row,1}]');
            varargout{m} = [BlockValues(row,1) Params([paramID{:}]',2) BlockValues(row,2:3)];
    end
end

% disconnect the DB connection
conn.Close();

% end of read_orion
end



function trialblob = parse_trialblob(raw_blob,blob_ver)
%PARSE_TRIALBLOB reads the user-defined blob structure.
trialblob = [];
raw_blob = cell2mat(raw_blob);

% return early if the blob is empty
if isempty(raw_blob), return; end

% get the length and version of the blob
trialblob.len = hex2int(raw_blob(1,1:4));
trialblob.ver = hex2int(raw_blob(1,5:8));

% manual correction of the blob version written in the data file
if 0~=blob_ver, trialblob.ver = blob_ver; end

% define your blob structure here.
% Branch a 'case' statement according to your blob version below and define
% a n-by-2 cell variable named 'format'. The first column of 'format'
% should have the variable types in the blob and the second column, the
% number of the variable.
% For example, if the version or your blob is 0 and it contains three
% booleans and one interger, you can define the 'format' like the following.
% case 0
%   format = { 'bool', 3; 'int', 1 };
%    or
%   format = { 'bool', 1; 'bool', 1; 'bool', 1; 'int', 1 };
%
% The data types you can use are
% 'double', 'float', 'unsigned long', 'long', 'unsigned int', 'int',
% 'unsigned short', 'short', 'unsigned char', 'char', and 'bool'.
switch trialblob.ver
    case 0	 % Example
        format = { 'unsigned int', 3; 'bool', 3; 'int', 1 };
    case 100 % Movie Task
        format = { 'unsigned int', 2; 'int', 4; 'unsigned long', 2; ...
            'bool', 8; 'char', 40; };
    case 200 % NMS Training
        format = { 'unsigned int', 2; 'int', 6; 'unsigned long', 7; ...
            'char', 80; 'bool', 16; 'int', 8; };
    case 201 % NMS Training
        format = { 'unsigned int', 2; 'int', 6; 'unsigned long', 7; ...
            'bool', 16; 'char', 80; 'int', 8; };
    case 300 % HT Traing
        format = { 'unsigned int', 2; 'int', 6; 'unsigned long', 7; ...
            'bool', 16; 'char', 100; };
    case {350,351} % HTA & HTB
        format = { 'unsigned int', 2; 'int', 4; 'unsigned long', 7; ...
            'bool', 16; 'char', 100; };
    case {400,401,402,403} % Cocktail Party or Face vs Non-face
        switch trialblob.len
            case 180    % 400
                format = { 'unsigned int', 2; 'int', 6; 'unsigned long', 8; ...
                    'bool', 16; 'char', 90; };
            case 188    % 401
                format = { 'unsigned int', 2; 'int', 6; 'unsigned long', 8; ...
                    'bool', 16; 'char', 90; 'int', 2; };
            case 196    % 401, 402, 403
                format = { 'unsigned int', 2; 'int', 6; 'unsigned long', 8; ...
                    'bool', 16; 'char', 100; 'int', 2; };
        end
    case 404 % Face vs Non-face
        format = { 'unsigned int', 2; 'int', 5; 'unsigned long', 5; ...
            'bool', 16; 'char', 100; };
    case 501 % Differential Reward
        format = { 'unsigned int', 2; 'int', 5; 'unsigned long', 7; ...
            'bool', 16; 'char', 60; };
    case 1000 % DR training
        format = { 'unsigned int', 2; 'unsigned long', 47; };
    case 1001 % DR training
        format = { 'unsigned int', 2; 'unsigned long', 19; };
    case 1002 % Delayed Reward
        switch trialblob.len
            case 108, format = { 'unsigned int', 2; 'unsigned long', 23; };
            case 112, format = { 'unsigned int', 2; 'unsigned long', 24; };
            case 120, format = { 'unsigned int', 2; 'unsigned long', 26; };
        end
    case 1003 % GB training
        format = { 'unsigned int', 2; 'unsigned long', 18; };
    case 1004 % Value Calibration
        format = { 'unsigned int', 2; 'unsigned long', 32; };
    case 1005 % Delayed Gratification
        format = { 'unsigned int', 2; 'unsigned long', 21; };
    case 1006 % Bundled Reward
        format = { 'unsigned int', 2; 'unsigned long', 20; };
    case 1007 % Self Control
        switch trialblob.len
            case 100, format = { 'unsigned int', 2; 'unsigned long', 21; };
            case 104, format = { 'unsigned int', 2; 'unsigned long', 22; };
            case 112, format = { 'unsigned int', 2; 'unsigned long', 24; };
            case 120, format = { 'unsigned int', 2; 'unsigned long', 25; 'int', 1; };
        end
    otherwise
        error('read_orion:parse_trialblob','The blob format of version %d is not defined',trialblob.ver);
end

% read the blob
count = 0;
pos = 8;
trialblob.data = zeros( size(raw_blob,1), sum([format{:,2}]) );
for m=1:size(format,1)
    switch lower(format{m,1})
        case 'double',         nbyte = 8;   signed = 2;
        case 'float',          nbyte = 4;   signed = 2;
        case 'unsigned long',  nbyte = 4;   signed = 0;
        case 'long',           nbyte = 4;   signed = 1;
        case 'unsigned',       nbyte = 4;   signed = 0;
        case 'unsigned int',   nbyte = 4;   signed = 0;
        case 'int',            nbyte = 4;   signed = 1;
        case 'unsigned short', nbyte = 2;   signed = 0;
        case 'short',          nbyte = 2;   signed = 1;
        case 'unsigned char',  nbyte = 1;   signed = 0;
        case 'char',           nbyte = 1;   signed = 1;
        case 'bool',           nbyte = 1;   signed = 1;
        otherwise, error('read_orion:parse_trialblob','Incorrect blob format');
    end
    pos = ceil(pos/nbyte)*nbyte;
    trialblob.data(:,count+1:count+format{m,2}) = hex2real(raw_blob(:,pos+1:pos+nbyte*format{m,2}),nbyte,signed);
    count = count + format{m,2};
    pos = pos + nbyte*format{m,2};
end

% end of parse_trialblob
end



function xyblob = parse_xyblob(raw_blob)
% This function reads XYblob structure such as eyeblob and joyblob. The
% return value 'xyblob' is a structure and contains 4 members.
% xyblob.ver : version of the xyblob
% xyblob.t_per_sample : sampling interval between data points in milliseconds
% xyblob.data : n-by-1 cell. n is the number of trials. Each cell contains
%               XY coordinates (in pixel) sampled in each trial.
% xyblob.mark : n-by-1 cell of eventcode markers.
%               Each cell consists of a 3-column matrix.
%               The 1st column is the event codes.
%               The 2nd column is extra information.
%               The numbers in the 3rd column point the first XY coordinate
%               in xyblob.data after the event occurred. You can use this
%               information to get eye positions during a particular
%               period. For example, let's say your eyeblob.mark looks like
%               the following.
%                  row 1 :    20     0     513
%                  row 2 :    94     1    1018
%               Then, you can get the eye trace between the event 20 and 94
%               like this.
%                  eyeblob.data(513:1017,:)
xyblob = [];

% return early if the blob is empty
if isempty(raw_blob{1}), return; end

% for the conversion from binary to number
mark = hex2int([190 186 173 222],4,1);

% get the version of the blob
xyblob.ver = hex2int(raw_blob{1}(5:8));
xyblob.t_per_samp = hex2double(raw_blob{1}(9:16));

for k=1:length(raw_blob)
    coordinate = reshape( hex2int(raw_blob{k}(17:end),4,1), 2, [] )';
    xyblob.data{k,1} = coordinate(mark~=coordinate(:,1),:);
    mark_list = find(mark==coordinate(:,1));
    mark_len = length(mark_list);
    xyblob.mark{k,1} = [mod(coordinate(mark_list,2),65536) floor(coordinate(mark_list,2)/65536) mark_list-(0:mark_len-1)'];
    %if 16384<=xyblob.mark{k,1}(1,1),   xyblob.mark{k,1}(1,1)   = xyblob.mark{k,1}(1,1) - 16384;   end
    %if 16384<=xyblob.mark{k,1}(end,1), xyblob.mark{k,1}(end,1) = 32767 - xyblob.mark{k,1}(end,1); end
end

% end of parse_xyblob
end



function n = hex2real(h,nbyte,type)
% h     : byte sequence you read from Orion blob
% nbyte : length of the data type
% type  : type of the number
%         0: unsigned fixed-point number (unsigned long, unsigned int, ...)
%         1: signed fixed-point number (long, int, short, char, bool)
%         2: single or double precision number (float or double)
% The default mode is 4-byte unsigned fixed-point number (e.g., unsigned
% long or unsigned int).
if nargin<2, nbyte=4; end
if nargin<3, type=0; end

if 2==type          % floating-point number
    if 4==nbyte     % single precision (float)
        n = hex2float(h);
    else            % double precision (double)
        n = hex2double(h);
    end
else                % fixed-point number
    n = hex2int(h,nbyte,type);
end

end



function n = hex2int(h,nbyte,signed)

if nargin<2, nbyte=4; end
if nargin<3, signed=0; end

h = double(h);
n = h(:,nbyte:nbyte:end);
for m = nbyte-1:-1:1
	n = n*256 + h(:,m:nbyte:end);
end

if 0~=signed
    ix = 128<=h(:,nbyte:nbyte:end);
    n(ix) = n(ix) - 2^(8*nbyte);
end

end



function n = hex2float(h)

h = double(h);
s = 128<=h(:,4:4:end);
e = mod(h(:,4:4:end),128)*2 + floor(h(:,3:4:end)/128);
f = h(:,1:4:end);
f = f/256 + h(:,2:4:end);
f = ( f/256 + mod(h(:,3:4:end),128) ) / 128;

n = (-1).^s .* 2.^(e-127) .* (1+f);

end



function n = hex2double(h)

h = double(h);
s = 128<=h(:,8:8:end);
e = mod(h(:,8:8:end),128)*16 + floor(h(:,7:8:end)/16);
f = h(:,1:8:end);
for m = 2:6
    f = f/256 + h(:,m:8:end);
end
f = ( f/256 + mod(h(:,7:8:end),16) ) / 16;

n = (-1).^s .* 2.^(e-1023) .* (1+f);

end
