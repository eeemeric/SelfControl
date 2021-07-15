READ_ORION reads behavioral data from the Orion ACCDB (*.mdb;*.accdb;)
  TRIALBLOB=READ_ORION
      If no argument is provided, this function opens a file dialog to
      get the path to an ACCDB file and returns trialblobs.
  VARARGOUT=READ_ORION(VARARGIN)
      When you have multiple blob types to read, you can give the list of
      them as VARARGIN. The possible types are 'event','trial','eye' and
      so on (see the code for the complete list). If one of the argument
      indicates a existing file, it is considered a path to the ACCDB.
      The return variables in VARARGOUT will be assigned as the order of
      the blob types in VARARGIN. For example,
      [TRIAL,EYE]=READ_ORION('trial','eye');
      [TRIAL,EYE]=READ_ORION(ACCDBPATH,'trial','eye');
      [TRIAL,EYE]=READ_ORION('trial','eye',ACCDBPATH);

  Apr 16, 2008        written by Jaewon Hwang (jaewon.hwang@gmail.com)
  Jul 20, 2008        added two more options, 'ver' & 'len'
  May 20, 2009        added ACCDB support (MS ACCESS 2007 format)
  May 31, 2009        added an option, 'info' which returns the neural
                      recording file name and the session note when
                      available
  May 30, 2013        added an option for reading ArmBlob, 'arm'
  Nov 21, 2013        added an option for reading Params, 'param'
