# nuDSP
A library of digital signal processing functions written in C.

## Secondary Functions:
- `dsp_printInputSampleValues`: Simply prints the sample values of the opened file.
- `lerp`: Performs linear interpolation between two samples.
- `dsp_get_gain`: Returns the peak dB in an audio file.
- `sinc`: Applies the sinc function to an input.
- `convolve`: Performs convolution on two input arrays.
- `normalize_double`: The same as dsp_normalize, but takes in doubles instead of floats.
- `createIRfromFFT_BSF`: Creates an impulse response of an ideal brickwall band-stop filter in the frequency domain, performs an inverse Fast Fourier Transform (IFFT) to convert it back to the time domain, and then fills the amplitude and phase response arrays with the real and imaginary parts of the resulting signal, respectively.
- `chebyshev_subroutine`: Computes intermediate filter coefficients for a specified pole of a Chebyshev filter, given the filter type, cutoff frequency, ripple factor, total number of poles, and current pole number.
- `chebyshev_get_coeffs`: Calculates and normalizes the final coefficients for a Chebyshev filter of a specified type, using a given cutoff frequency, ripple factor, and number of poles, through multiple calls to `chebyshev_subroutine`.

## Main Functions:
- `dsp_reverse`: Reverses the audio file.
- `dsp_gainChange`: Creates a gain change. The change amount is passed in using the dBChange variable.
- `dsp_normalize`: Normalizes an audio file. The threshold is passed in using the normalizeThresholdInDB variable.
- `dsp_fadeIn`: Applies a fade in to an audio file. Linear, equal power, or s-shaped.
- `dsp_fadeOut`: Applies a fade out to an audio file. Linear, equal power, or s-shaped.
- `dsp_simpleSquarewave`: Generates a square wave.
- `dsp_simpleSinewave`: Generates a sine wave.
- `dsp_simpleTrianglewave`: Generates a triangle wave.
- `dsp_rampSinewave`: Generates a sine wave with a changing frequency.
- `dsp_additiveSquarewave`: Generates a square wave using additive synthesis.
- `dsp_additiveTrianglewave`: Generates a triangle wave using additive synthesis.
- `dsp_tremolo`: Applies a tremolo effect to an input file.
- `dsp_ampModulation`: Applies an amplitude modulation effect to an input file.
- `dsp_ringModulation`: Applies a ring modulation effect to an input file.
- `dsp_flanger`: Applies a sine flanger effect to an input file.
- `fir_lpf`: Applies an FIR lowpass filter to an audio file.
- `fir_bandstop`: Applies an FIR bandstop filter to an audio file.
- `iir_4stage_lpf`: Applies a 4-stage IIR lowpass filter to an audio file.
- `iir_bpf`: Applies an IIR bandpass filter to an audio file.
- `chebyshev_filter`: Applies a Chebyshev filter to an audio file.
- `dspa_pitchChange`: Applies a pitch change to an audio file.
