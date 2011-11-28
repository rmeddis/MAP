function x = MTprofile12_41hr04_Oct_2011
%created: 12_41hr04_Oct_2011

x.BFs = [250   500  1000  2000  4000  8000];

x.LongTone = [40.1      42.7      53.1      81.2      92.3      82.3];
x.ShortTone = [49.6      48.3      64.4        89      95.5      92.3];

x.Gaps = [0.01      0.03      0.05      0.07      0.09];
x.TMCFreq = [250   500  1000  2000  4000  8000];
x.TMC = [
60.6	58.5	77.7	NaN	NaN	NaN	 
62.8	56	85.2	NaN	NaN	NaN	 
60.8	60.2	81.6	NaN	NaN	NaN	 
75.2	67.5	85.4	NaN	NaN	NaN	 
83.9	77.4	93.2	NaN	NaN	NaN	 
];
x.TMC = x.TMC';

x.MaskerRatio = [0.5      0.7      0.9        1      1.1      1.3      1.6];
x.IFMCFreq = [250   500  1000  2000  4000  8000];
x.IFMCs = [
60.9	81.2	89.1	105	NaN	NaN	 
57.6	69	85.4	105	NaN	NaN	 
52.4	56.4	77.2	106	NaN	NaN	 
52.8	57.7	78.3	105	NaN	NaN	 
51.8	57.2	78.6	NaN	NaN	NaN	 
57.2	62.1	89.3	NaN	NaN	NaN	 
66.8	93.5	103	NaN	NaN	NaN	 
];
x.IFMCs = x.IFMCs';
