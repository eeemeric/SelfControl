# Extract spikes times, spike waveforms, and behavioral events from a TDT tank.

## TDTBIN2MAT  TDT tank data extraction.
data = TDTbin2mat(BLOCK_PATH), where BLOCK_PATH is a string, retrieve all data from specified block directory in struct format. This reads the binary tank data and requires no Windows-based software.

  - data.epocs      contains all epoc store data (onsets, offsets, values)
  - data.snips      contains all snippet store data (timestamps, channels, and raw data)
  - data.streams    contains all continuous data (sampling rate and raw data)
  - data.scalars    contains all scalar data (samples and timestamps)
  - data.info       contains additional information about the block

  'parameter', value pairs
     'T1'         scalar, retrieve data starting at T1 (default = 0 for
                      beginning of recording)
     'T2'         scalar, retrieve data ending at T2 (default = 0 for end
                      of recording)
     'SORTNAME'   string, specify sort ID to use when extracting snippets
     'TYPE'       array of scalars or cell array of strings, specifies
                      what type of data stores to retrieve from the tank
                    1: all (default)
                    2: epocs
                    3: snips
                    4: streams
                    5: scalars
                    TYPE can also be cell array of any combination of
                      'epocs', 'streams', 'scalars', 'snips', 'all'
                    examples:
                      data = TDTbin2mat(BLOCKPATH,'TYPE',[1 2]);
                          > returns only epocs and snips
                      data = TDTbin2mat(BLOCKPATH,'TYPE',{'epocs','snips'});
                          > returns only epocs and snips
     'RANGES'     array of valid time range column vectors
     'NODATA'     boolean, only return timestamps, channels, and sort 
                      codes for snippets, no waveform data (default = false)
     'STORE'      string, specify a single store to extract
                  cell of strings, specify cell arrow of stores to extract
     'CHANNEL'    integer or array, choose a single channel or array of
                      channels to extract from stream or snippet events
                      Default is 0, to extract all channels.
     'BITWISE'    string, specify an epoc store or scalar store that 
                      contains individual bits packed into a 32-bit 
                      integer. Onsets/offsets from individual bits will
                      be extracted.
     'HEADERS'    var, set to 1 to return only the headers for this
                      block, so that you can make repeated calls to read
                      data without having to parse the TSQ file every
                      time. Or, pass in the headers using this parameter.
                  example:
                      heads = TDTbin2mat(BLOCK_PATH, 'HEADERS', 1);
                      data = TDTbin2mat(BLOCK_PATH, 'HEADERS', heads, 'TYPE', {'snips'});
                      data = TDTbin2mat(BLOCK_PATH, 'HEADERS', heads, 'TYPE', {'streams'});
     'COMBINE'    cell, specify one or more data stores that were saved 
                      by the Strobed Data Storage gizmo in Synapse (or an
                      Async_Stream_Store macro in OpenEx). By default,
                      the data is stored in small chunks while the strobe
                      is high. This setting allows you to combine these
                      small chunks back into the full waveforms that were
                      recorded while the strobe was enabled.
                  example:
                      data = TDTbin2mat(BLOCK_PATH, 'COMBINE', {'StS1'});

