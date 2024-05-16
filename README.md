 SAND2024-01061O The HiFiPDV-2 software is similar to other photonic Doppler velocimetry (PDV) 
 analysis codes in that it converts signals from intensity-time space to frequency-time space 
 using a short-time Fourier transform (STFT). It then analyzes peaks of the resulting power 
 spectral density information to extract the velocity-time history. The unique aspects of 
 HiFiPDV and HiFiPDV-2 are that this reduction process is repeated for an array of reasonable 
 STFT input variables. The resulting population distribution of output velocity-time histories 
 is evaluated to determine the most likely velocity history and its associated uncertainty. 
 This uncertainty is defined as the systematic uncertainty of the PDV signal as it is a 
 function of input values to the reduction process. The HiFiPDV-2 program brings new 
 functionality and efficiency to these calculations—most notably significant reductions in 
 memory usage and calculation times, the ability to evaluate frequency-shifted signals, and 
 the ability to filter out baseline signals. Developers recommend that users have at least 
 2 GB of RAM per CPU core in their system. Sandia National Laboratories is a multimission 
 laboratory managed and operated by National Technology & Engineering Solutions of Sandia, LLC,
 a wholly owned subsidiary of Honeywell International Inc., for the U.S. Department of Energy’s
 National Nuclear Security Administration under contract DE-NA0003525.
