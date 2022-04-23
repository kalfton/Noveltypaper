function Hist_Smooth=Smooth_Histogram(Hist_Raw,Smoothing)
switch Smoothing
case 1
%   disp('Convolving with BOXCAR of BINSIZE=20ms')
   BinSize=20;
   Half_BW=(round(BinSize)/2);
   BinSize=Half_BW*2;
   Kernel=ones(1,(Half_BW*2)+1);
   Kernel=Kernel.*1000/sum(Kernel);
case 2
%   disp('Convolving with GAUSSIAN of SIGMA=5ms')
   Sigma=100;
   Kernel=[-5*Sigma:5*Sigma];
   BinSize=length(Kernel);
   Half_BW=(BinSize-1)/2;
   Kernel=[-BinSize/2:BinSize/2];   
   Factor=1000/(Sigma*sqrt(2*pi));
   Kernel=Factor*(exp(-(0.5*((Kernel./Sigma).^2))));
case 3
%   disp('Convolving with EPSP of GROWTH=1ms and DECAY=20ms')
   Growth=1;
   Decay=20;
   Half_BW=round(Decay*8);
   BinSize=(Half_BW*2)+1;
   Kernel=[0:Half_BW];
   Half_Kernel=(1-(exp(-(Kernel./Growth)))).*(exp(-(Kernel./Decay)));
   Half_Kernel=Half_Kernel./sum(Half_Kernel);
   Kernel(1:Half_BW)=0;
   Kernel(Half_BW+1:BinSize)=Half_Kernel;
   Kernel=Kernel.*1000;
otherwise
end
Hist_Raw=Hist_Raw';%Make Hist_Raw Rows(which are Trials) to columns
Kernel=Kernel'; %make Kernel to a column vector 
%Now Convolve Hist_raw columns with Column Kernel
Hist_Smooth=convn(Hist_Raw,Kernel,'same');
[mrows mcols]=size(Hist_Smooth);
Hist_Smooth=(Hist_Smooth)';
