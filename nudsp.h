#include "audioLib.h"
#include "AnalysisDisplay.h"
#include "Spectrogram.h"
#include "time.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SUCCESS 0;
#define INVALID_OUTPUT_POINTER -1;
#define INVALID_NUM_SAMPLES -2;
#define INVALID_FREQ -3;
#define INVALID_DB -4;
#define INVALID_SAMPLE_RATE -5;
#define INVALID_DURATION -6;
#define INVALID_INPUT_POINTER -7;
#define INVALID_DEPTH -8;
#define INVALID_AMP -9;
#define INVALID_RIPPLE -10;
#define INVALID_POLES -11;
#define INVALID_FILTER -12;

#pragma mark FUNCTION_DECLARATIONS
//.................................................................................................................. dsp_printInputSampleValues
// FUNCTION:    dsp_printInputSampleValues
// DECRIPTION:  Simply prints the sample values of the opened file.
//
// PARAMETERS:  iAudioPtr        Pointer to the input audio buffer
//              iNumSamples      Number of sample frames in the input file
//
// RETURN:       0      SUCCESS
//
//..................................................................................................................
int dsp_printInputSampleValues(float* iAudioPtr, int iNumSamples);

//.................................................................................................................. dsp_reverse
// FUNCTION:    dsp_reverse
// DECRIPTION:  Reverses the audio file.
//
// PARAMETERS:  iAudioPtr        Pointer to the input audio buffer
//              iNumSamples      Number of sample frames in the input file
//              oAudioPtr        Pointer to an output buffer. The pointer is created by the calling class,
//                               but the buffer is allocated within this function. It must be freed by the
//                               caller.
//
//
// RETURN:       0      SUCCESS
//              -1      FAILED TO ALLOCATE MEMORY FOR OUTPUT
//              -2      INVALID INPUT PARAMETER
//
//..................................................................................................................
int dsp_reverse(float* iAudioPtr, int iNumSamples, float* oAudioPtr);

//.................................................................................................................. dsp_gainChange
// FUNCTION:    dsp_gainChange
// DECRIPTION:  Creates a gain change. The change amount is passed in using the dBChange variable.
//
// PARAMETERS:  iAudioPtr        Pointer to the input audio buffer
//              iNumSamples      Number of sample frames in the input file
//              oAudioPtr        Pointer to an output buffer. The pointer is created by the calling class,
//                               but the buffer is allocated within this function. It must be freed by the
//                               caller.
//              dBChange         A float that indicates the change in decibels. (-60 dB to +60 dB)
//
//
// RETURN:       0      SUCCESS
//              -1      FAILED TO ALLOCATE MEMORY FOR OUTPUT
//              -2      INVALID INPUT PARAMETER
//
//..................................................................................................................
int dsp_gainChange(float* iAudioPtr, int iNumSamples, float* oAudioPtr, float dBChange);

//.................................................................................................................. dsp_normalize
// FUNCTION:    dsp_normalize
// DECRIPTION:  Normalizes an audio file. The threshold is passed in using the normalizeThresholdInDB variable.
//
// PARAMETERS:  iAudioPtr                       Pointer to the input audio buffer
//              iNumSamples                     Number of sample frames in the input file
//              oAudioPtr                       Pointer to an output buffer. The pointer is created by the calling class,
//                                              but the buffer is allocated within this function. It must be freed by the
//                                               caller.
//              normalizeThresholdInDB          A float that indicates the threshold in decibels. (<=0dB)
//
//
// RETURN:       0      SUCCESS
//              -1      FAILED TO ALLOCATE MEMORY FOR OUTPUT
//              -2      INVALID INPUT PARAMETER
//
//..................................................................................................................
int dsp_normalize(float* iAudioPtr, int iNumSamples, float* oAudioPtr, float normalizeThresholdInDB);

//.................................................................................................................. dsp_fadeIn
// FUNCTION:    dsp_fadeIn
// DECRIPTION:  Applies a fade in to an audio file.
//
// PARAMETERS:  iAudioPtr                       Pointer to the input audio buffer
//              iNumSamples                     Number of sample frames in the input file
//              oAudioPtr                       Pointer to an output buffer. The pointer is created by the calling class,
//                                              but the buffer is allocated within this function. It must be freed by the
//                                              caller.
//              durationInMS                    The length of the fade in milliseconds.
//              sampleRate                      The sample rate of the input file.
//              fadeType                        The type of fade. Linear (1), equal power (2), and s-shaped (3).
//
//
// RETURN:       0      SUCCESS
//              -1      INVALID INPUT FILE
//              -2      CORRUPT OUTPUT FILE
//              -3      MISSING SAMPLE DATA
//              -4      INVALID FADE TYPE
//              -5      INVALID DURATION
//              -6      INVALID SAMPLE RATE
//..................................................................................................................
int dsp_fadeIn(float* iAudioPtr, int iNumSamples, float* oAudioPtr, int durationInMS, int sampleRate, short fadeType);

//.................................................................................................................. dsp_fadeOut
// FUNCTION:    dsp_fadeOut
// DECRIPTION:  Applies a fade out to an audio file.
//
// PARAMETERS:  iAudioPtr                       Pointer to the input audio buffer
//              iNumSamples                     Number of sample frames in the input file
//              oAudioPtr                       Pointer to an output buffer. The pointer is created by the calling class,
//                                              but the buffer is allocated within this function. It must be freed by the
//                                              caller.
//              durationInMS                    The length of the fade in milliseconds.
//              sampleRate                      The sample rate of the input file.
//              fadeType                        The type of fade. Linear (1), equal power (2), and s-shaped (3).
//
//
// RETURN:       0      SUCCESS
//              -1      INVALID INPUT FILE
//              -2      CORRUPT OUTPUT FILE
//              -3      MISSING SAMPLE DATA
//              -4      INVALID FADE TYPE
//              -5      INVALID DURATION
//              -6      INVALID SAMPLE RATE
//
//..................................................................................................................
int dsp_fadeOut(float* iAudioPtr, int iNumSamples, float* oAudioPtr, int durationInMS, int sampleRate, short fadeType);


//.................................................................................................................. dsp_simpleSquarewave
// FUNCTION:    dsp_simpleSquarewave
// DECRIPTION:  Generates a square wave.
//
// PARAMETERS:  oAudioPtr                       Pointer to an output buffer. The pointer is created by the calling class,
//                                              but the buffer is allocated within this function. It must be freed by the
//                                              caller.
//              iNumSamples                     Number of sample frames to output.
//              freq                            The frequency of the wave.
//              amp                             The amplitude of the wave.
//              sampleRate                      The sample rate of the output file.
//
// RETURN:       0      SUCCESS
//              -1      INVALID OUTPUT POINTER
//              -2      INVALID NUMBER OF SAMPLES
//              -3      INVALID FREQUENCY
//              -4      INVALID AMPLITUDE
//              -5      INVALID SAMPLE RATE
//
//..................................................................................................................
int dsp_simpleSquarewave(float* oAudioPtr, int iNumSamples, float freq, float amp, int sampleRate);

//.................................................................................................................. dsp_simpleSinewave
// FUNCTION:    dsp_simpleSinewave
// DECRIPTION:  Generates a sine wave.
//
// PARAMETERS:  oAudioPtr                       Pointer to an output buffer. The pointer is created by the calling class,
//                                              but the buffer is allocated within this function. It must be freed by the
//                                              caller.
//              iNumSamples                     Number of sample frames to output.
//              freq                            The frequency of the wave.
//              amp                             The amplitude of the wave.
//              sampleRate                      The sample rate of the output file.
//
// RETURN:       0      SUCCESS
//              -1      INVALID OUTPUT POINTER
//              -2      INVALID NUMBER OF SAMPLES
//              -3      INVALID FREQUENCY
//              -4      INVALID AMPLITUDE
//              -5      INVALID SAMPLE RATE
//
//..................................................................................................................
int dsp_simpleSinewave(float* oAudioPtr, int iNumSamples, float freq, float amp, int sampleRate);

//.................................................................................................................. dsp_simpleTrianglewave
// FUNCTION:    dsp_simpleTrianglewave
// DECRIPTION:  Generates a triangle wave.
//
// PARAMETERS:  oAudioPtr                       Pointer to an output buffer. The pointer is created by the calling class,
//                                              but the buffer is allocated within this function. It must be freed by the
//                                              caller.
//              iNumSamples                     Number of sample frames to output.
//              freq                            The frequency of the wave.
//              amp                             The amplitude of the wave.
//              sampleRate                      The sample rate of the output file.
//
// RETURN:       0      SUCCESS
//              -1      INVALID OUTPUT POINTER
//              -2      INVALID NUMBER OF SAMPLES
//              -3      INVALID FREQUENCY
//              -4      INVALID AMPLITUDE
//              -5      INVALID SAMPLE RATE
//
//..................................................................................................................
int dsp_simpleTrianglewave(float* oAudioPtr, int iNumSamples, float freq, float amp, int sampleRate);

//.................................................................................................................. dsp_rampSinewave
// FUNCTION:    dsp_rampSinewave
// DECRIPTION:  Generates a sine wave with a changing frequency.
//
// PARAMETERS:  oAudioPtr                       Pointer to an output buffer. The pointer is created by the calling class,
//                                              but the buffer is allocated within this function. It must be freed by the
//                                              caller.
//              nSamples                        Number of sample frames to output.
//              startingFreq                    The starting frequency of the wave.
//              endingFreq                      The ending frequency of the wave.
//              gain_dB                         The level of the wave in decibels.
//              sampleRate                      The sample rate of the output file.
//
// RETURN:       0      SUCCESS
//              -1      INVALID OUTPUT POINTER
//              -2      INVALID NUMBER OF SAMPLES
//              -3      INVALID FREQUENCY
//              -4      INVALID DECIBELS
//              -5      INVALID SAMPLE RATE
//
//..................................................................................................................
int dsp_rampSinewave(float* oAudioPtr, int nSamples, float startingFreq, float endingFreq, float gain_dB, int sampleRate);

//.................................................................................................................. dsp_additiveSquarewave
// FUNCTION:    dsp_additiveSquarewave
// DECRIPTION:  Generates a square wave using additive synthesis.
//
// PARAMETERS:  oAudioPtr                       Pointer to an output buffer. The pointer is created by the calling class,
//                                              but the buffer is allocated within this function. It must be freed by the
//                                              caller.
//              nSamples                        Number of sample frames to output.
//              freq                            The frequency of the wave.
//              gain_dB                         The level of the wave in decibels.
//              sampleRate                      The sample rate of the output file.
//
// RETURN:       0      SUCCESS
//              -1      INVALID OUTPUT POINTER
//              -2      INVALID NUMBER OF SAMPLES
//              -3      INVALID FREQUENCY
//              -4      INVALID DECIBELS
//              -5      INVALID SAMPLE RATE
//
//..................................................................................................................
int dsp_additiveSquarewave(float* oAudioPtr, int nSamples, float freq, float gain_dB, int sampleRate);

//.................................................................................................................. dsp_additiveTrianglewave
// FUNCTION:    dsp_additiveTrianglewave
// DECRIPTION:  Generates a triangle wave using additive synthesis.
//
// PARAMETERS:  oAudioPtr                       Pointer to an output buffer. The pointer is created by the calling class,
//                                              but the buffer is allocated within this function. It must be freed by the
//                                              caller.
//              nSamples                        Number of sample frames to output.
//              freq                            The frequency of the wave.
//              gain_dB                         The level of the wave in decibels.
//              sampleRate                      The sample rate of the output file.
//
// RETURN:       0      SUCCESS
//              -1      INVALID OUTPUT POINTER
//              -2      INVALID NUMBER OF SAMPLES
//              -3      INVALID FREQUENCY
//              -4      INVALID DECIBELS
//              -5      INVALID SAMPLE RATE
//
//..................................................................................................................
int dsp_additiveTrianglewave(float* oAudioPtr, int nSamples, float freq, float gain_dB, int sampleRate);

//.................................................................................................................. dsp_tremolo
// FUNCTION:    dsp_tremolo
// DECRIPTION:  Applies a tremolo effect to an input file.
//
// PARAMETERS:  iAudioPtr                       Pointer to the input audio buffer.
//              oAudioPtr                       Pointer to an output buffer. The pointer is created by the calling class,
//                                              but the buffer is allocated within this function. It must be freed by the
//                                              caller.
//              iNumSamples                     Number of samples in the input file.
//              lfoStartRate                    The starting frequency of the LFO.
//              lfoEndRate                      The ending frequency of the LFO.
//              lfoDepth                        The depth of the effect in ampltitude.
//              sampleRate                      The sample rate of the input/output file.
//
// RETURN:       0      SUCCESS
//              -1      INVALID OUTPUT POINTER
//              -2      INVALID NUMBER OF SAMPLES
//              -3      INVALID FREQUENCY
//              -5      INVALID SAMPLE RATE
//              -7      INVALID INPUT POINTER
//              -8      INVALID DEPTH
//
//..................................................................................................................
int dsp_tremolo(float* iAudioPtr, int iNumSamples, float* oAudioPtr, float lfoStartRate, float lfoEndRate, float lfoDepth, int sampleRate);

//.................................................................................................................. dsp_ampModulation
// FUNCTION:    dsp_ampModulation
// DECRIPTION:  Applies an amplitude modulation effect to an input file.
//
// PARAMETERS:  iAudioPtr                       Pointer to the input audio buffer.
//              oAudioPtr                       Pointer to an output buffer. The pointer is created by the calling class,
//                                              but the buffer is allocated within this function. It must be freed by the
//                                              caller.
//              iNumSamples                     Number of samples in the input file.
//              modFrequency                    The frequency of the modulator.
//              modAmplitude                    The amplitude of the modulator.
//              sampleRate                      The sample rate of the input/output file.
//
// RETURN:       0      SUCCESS
//              -1      INVALID OUTPUT POINTER
//              -2      INVALID NUMBER OF SAMPLES
//              -3      INVALID FREQUENCY
//              -5      INVALID SAMPLE RATE
//              -7      INVALID INPUT POINTER
//              -9      INVALID AMPLITUDE
//
//..................................................................................................................
int dsp_ampModulation(float* iAudioPtr, int iNumSamples, float* oAudioPtr, float modFrequency, float modAmplitude, int sampleRate);

//.................................................................................................................. dsp_ringModulation
// FUNCTION:    dsp_ringModulation
// DECRIPTION:  Applies a ring modulation effect to an input file.
//
// PARAMETERS:  iAudioPtr                       Pointer to the input audio buffer.
//              oAudioPtr                       Pointer to an output buffer. The pointer is created by the calling class,
//                                              but the buffer is allocated within this function. It must be freed by the
//                                              caller.
//              iNumSamples                     Number of samples in the input file.
//              modFrequency                    The frequency of the modulator.
//              modAmplitude                    The amplitude of the modulator.
//              sampleRate                      The sample rate of the input/output file.
//
// RETURN:       0      SUCCESS
//              -1      INVALID OUTPUT POINTER
//              -2      INVALID NUMBER OF SAMPLES
//              -3      INVALID FREQUENCY
//              -5      INVALID SAMPLE RATE
//              -7      INVALID INPUT POINTER
//              -9      INVALID AMPLITUDE
//
//..................................................................................................................
int dsp_ringModulation(float* iAudioPtr, int iNumSamples, float* oAudioPtr, float modFrequency, float modAmplitude, int sampleRate);

//.................................................................................................................. dsp_flanger
// FUNCTION:    dsp_flanger
// DECRIPTION:  Applies a sine flanger effect to an input file.
//
// PARAMETERS:  iAudioPtr                       Pointer to the input audio buffer.
//              oAudioPtr                       Pointer to an output buffer. The pointer is created by the calling class,
//                                              but the buffer is allocated within this function. It must be freed by the
//                                              caller.
//              iNumSamples                     Number of samples in the input file.
//              modFrequency                    The frequency of the sine LFO. (<= 20 Hz)
//              modAmplitude                    The amplitude of the LFO. (> 0, <= 1) Determines delay time. (> 0 ms, <= 20 ms)
//              sampleRate                      The sample rate of the input/output file. Must be 44.1 kHz.
//
// RETURN:       0      SUCCESS
//              -1      INVALID OUTPUT POINTER
//              -2      INVALID NUMBER OF SAMPLES
//              -3      INVALID FREQUENCY
//              -5      INVALID SAMPLE RATE
//              -7      INVALID INPUT POINTER
//              -9      INVALID AMPLITUDE
//
//..................................................................................................................
int dsp_flanger(float* iAudioPtr, int iNumSamples, float* oAudioPtr, float modFrequency, float modAmplitude, int sampleRate);

//.................................................................................................................. fir_lpf
// FUNCTION:    fir_lpf
// DECRIPTION:  Applies an FIR lowpass filter to an audio file.
//
// PARAMETERS:  iAudioPtr                       Pointer to the input audio buffer.
//              oAudioPtr                       Pointer to an output buffer. The pointer is created by the calling class,
//                                              but the buffer is allocated within this function. It must be freed by the
//                                              caller.
//              oNumSamples                     Number of samples to output.
//              cutoff                          The cutoff frequency.
//              firSize                         The size of the filter window.
//              sampleRate                      The sample rate of the input/output file. Must be 44.1 kHz.
//
// RETURN:       0      SUCCESS
//              -1      INVALID OUTPUT POINTER
//              -2      INVALID NUMBER OF SAMPLES
//              -3      INVALID FREQUENCY
//              -5      INVALID SAMPLE RATE
//              -7      INVALID INPUT POINTER
//              -9      INVALID AMPLITUDE
//
//..................................................................................................................
int fir_lpf(float* iAudioPtr, int iNumSamples, float* oAudioPtr, int oNumSamples, double cutoff, int firSize, int sampleRate);

//.................................................................................................................. fir_bandstop
// FUNCTION:    fir_bandstop
// DECRIPTION:  Applies an FIR bandstop filter to an audio file.
//
// PARAMETERS:  iAudioPtr                       Pointer to the input audio buffer.
//              oAudioPtr                       Pointer to an output buffer. The pointer is created by the calling class,
//                                              but the buffer is allocated within this function. It must be freed by the
//                                              caller.
//              oNumSamples                     Number of samples to output.
//              lowerEdgeInHz                   The lower cutoff frequency.
//              upperEdgeInHz                   The upper cutoff frequency.
//              firSize                         The size of the filter window.
//              sampleRate                      The sample rate of the input/output file. Must be 44.1 kHz.
//
// RETURN:       0      SUCCESS
//              -1      INVALID OUTPUT POINTER
//              -2      INVALID NUMBER OF SAMPLES
//              -3      INVALID FREQUENCY
//              -5      INVALID SAMPLE RATE
//              -7      INVALID INPUT POINTER
//              -9      INVALID AMPLITUDE
//
//..................................................................................................................
int fir_bandstop(float* iAudioPtr, int iNumSamples, float* oAudioPtr, int oNumSamples, float lowerEdgeInHz, float upperEdgeInHz, int firSize, int sampleRate);

//.................................................................................................................. iir_4stage_lpf
// FUNCTION:    iir_4stage_lpf
// DECRIPTION:  Applies a 4-stage IIR lowpass filter to an audio file.
//
// PARAMETERS:  iAudioPtr                       Pointer to the input audio buffer.
//              oAudioPtr                       Pointer to an output buffer. The pointer is created by the calling class,
//                                              but the buffer is allocated within this function. It must be freed by the
//                                              caller.
//              oNumSamples                     Number of samples to output.
//              cutoff                          The cutoff frequency.
//              sampleRate                      The sample rate of the input/output file. Must be 44.1 kHz.
//
// RETURN:       0      SUCCESS
//              -1      INVALID OUTPUT POINTER
//              -2      INVALID NUMBER OF SAMPLES
//              -3      INVALID FREQUENCY
//              -5      INVALID SAMPLE RATE
//              -7      INVALID INPUT POINTER
//              -9      INVALID AMPLITUDE
//
//..................................................................................................................
int iir_4stage_lpf(float* iAudioPtr, int iNumSamples, float* oAudioPtr, int oNumSamples, double cutoff, int sampleRate);

//.................................................................................................................. iir_bpf
// FUNCTION:    iir_bpf
// DECRIPTION:  Applies an IIR bandpass filter to an audio file.
//
// PARAMETERS:  iAudioPtr                       Pointer to the input audio buffer.
//              oAudioPtr                       Pointer to an output buffer. The pointer is created by the calling class,
//                                              but the buffer is allocated within this function. It must be freed by the
//                                              caller.
//              oNumSamples                     Number of samples to output.
//              center                          The center frequency.
//              bandwidth                       The bandwidth.
//              sampleRate                      The sample rate of the input/output file. Must be 44.1 kHz.
//
// RETURN:       0      SUCCESS
//              -1      INVALID OUTPUT POINTER
//              -2      INVALID NUMBER OF SAMPLES
//              -3      INVALID FREQUENCY
//              -5      INVALID SAMPLE RATE
//              -7      INVALID INPUT POINTER
//              -9      INVALID AMPLITUDE
//
//..................................................................................................................
int iir_bpf(float* iAudioPtr, int iNumSamples, float* oAudioPtr, int oNumSamples, double center, double bandwidth, int sampleRate);

//.................................................................................................................. chebyshev_filter
// FUNCTION:    chebyshev_filter
// DECRIPTION:  Applies a Chebyshev filter to an audio file.
//
// PARAMETERS:  iAudioPtr                       Pointer to the input audio buffer.
//              oAudioPtr                       Pointer to an output buffer. The pointer is created by the calling class,
//                                              but the buffer is allocated within this function. It must be freed by the
//                                              caller.
//              oNumSamples                     Number of samples to output.
//              filterType                      (0) - LPF | (1) - HPF
//              freq                            The cutoff freq.
//              ripple                          The ripple percentage.
//              numPoles                        The number of poles in the filter.
//              sampleRate                      The sample rate of the input/output file. Must be 44.1 kHz.
//
// RETURN:       0      SUCCESS
//              -1      INVALID OUTPUT POINTER
//              -2      INVALID NUMBER OF SAMPLES
//              -3      INVALID FREQUENCY
//              -5      INVALID SAMPLE RATE
//              -7      INVALID INPUT POINTER
//              -9      INVALID AMPLITUDE
//              -10     INVALID_RIPPLE
//              -11     INVALID POLES
//              -12     INVALID FILTER
//
//..................................................................................................................
int chebyshev_filter(float* iAudioPtr, int iNumSamples, float* oAudioPtr, short filterType, double freq, double ripple, int numPoles, int sampleRate);

//.................................................................................................................. dspa_pitchChange
// FUNCTION:    dspa_pitchChange
// DECRIPTION:  Applies a pitch change to an audio file.
//
// PARAMETERS:  iAudioPtr                       Pointer to the input audio buffer.
//              oAudioPtr                       Pointer to an output buffer. The pointer is created by the calling class,
//                                              but the buffer is allocated within this function. It must be freed by the
//                                              caller.
//              oNumSamples                     Number of samples to output.
//              filterType                      (0) - LPF | (1) - HPF
//              numerator                       The numerator of the pitch change percentage.
//              denominator                     The denominator of the pitch change percentage.
//              sampleRate                      The sample rate of the input/output file. Must be 44.1 kHz.
//
// RETURN:       0      SUCCESS
//              -1      INVALID OUTPUT POINTER
//              -2      INVALID NUMBER OF SAMPLES
//              -5      INVALID SAMPLE RATE
//              -7      INVALID INPUT POINTER
//
//..................................................................................................................
int dspa_pitchChange(float* iAudioPtr, int iNumSamples, float* oAudioPtr, int oNumSamples, int numerator, int denominator, int sampleRate);

#pragma mark FUNCTION_IMPLEMENTATIONS

// Linear Interpolation Helper
float lerp(float v0, float v1, float t) {
    return v0 + t * (v1 - v0);
}

int dsp_printInputSampleValues(float* iAudioPtr, int iNumSamples)
{
    int returnCode = -1; 
    float* i_ptr = iAudioPtr;
    
    for(int i = 0; i < iNumSamples; i++)
    {
        printf("a: %f\n", i_ptr[i]);
    }
    return returnCode;
}

int dsp_reverse(float* iAudioPtr, int iNumSamples, float* oAudioPtr)
{
    int returnCode = 0;
    float* input = iAudioPtr;
    float* output = oAudioPtr;

    for (int i = 0; i < iNumSamples; i++)
    {
        output[iNumSamples - 1 - i] = input[i];
    }
    printf("reverse\n");
    return returnCode;
}

int dsp_gainChange(float* iAudioPtr, int iNumSamples, float* oAudioPtr, float dBChange)
{
    int returnCode = 0;

    // check inputs
    if (iAudioPtr == NULL || oAudioPtr == NULL || iNumSamples <= 0) {
        printf("Error: Input file is corrupt or invalid!\n");
        returnCode = -2;
        return returnCode;
    }

    float* input = iAudioPtr;
    float* output = oAudioPtr;
    float ampChange = 0;

    // convert dB change to amplitude change
    ampChange = powf(10, (dBChange / 20));

    // apply amplitude change factor to all samples in file
    for (int i = 0; i < iNumSamples; i++)
    {
        output[i] = input[i] * ampChange;
    }
    printf("gain change\n");
    return returnCode;
}

int dsp_get_gain(float* iAudioPtr, int iNumSamples) {
    float peakAmp = 0;
    float peakdB = 0;

    // find peak amplitude sample in file
    for (int i = 0; i < iNumSamples; i++)
    {
        if (iAudioPtr[i] > peakAmp) {
            peakAmp = iAudioPtr[i];
        }
    }

    // convert peak amplitude to dB
    peakdB = 20 * log10f(fabs(peakAmp));
    return peakdB;
}

int dsp_normalize(float* iAudioPtr, int iNumSamples, float* oAudioPtr, float normalizeThresholdInDB)
{
    int returnCode = 0;

    // check inputs
    if (iAudioPtr == NULL || oAudioPtr == NULL || iNumSamples <= 0) {
        printf("Error: Input file is corrupt or invalid!\n");
        returnCode = INVALID_INPUT_POINTER;
        return returnCode;
    }

    float* input = iAudioPtr;
    float* output = oAudioPtr;
    float dBChange = 0;
    float peakAmp = 0;
    float peakdB = 0;

    // find peak amplitude sample in file
    for (int i = 0; i < iNumSamples; i++)
    {
        if (input[i] > peakAmp) {
            peakAmp = input[i];
        }
    }

    // convert peak amplitude to dB
    peakdB = 20 * log10f(fabs(peakAmp));

    // check range of threshold parameter
    if (normalizeThresholdInDB > 0) {
        printf("Error: Threshold must be less than or equal to 0 dB!\n");
        returnCode = -2;
        dBChange = 0;
    }

    // calculate dB change needed to normalize to threshold
    else {
        dBChange = normalizeThresholdInDB - peakdB;
    }

    // pass needed dB change to gain change function
    returnCode = dsp_gainChange(iAudioPtr, iNumSamples, oAudioPtr, dBChange);
    printf("normalize\n");
    return returnCode;
}

int dsp_fadeIn(float* iAudioPtr, int iNumSamples, float* oAudioPtr, int durationInMS, int sampleRate, short fadeType) {
    int returnCode = 0;

    // check inputs
    if (iAudioPtr == NULL) {
        printf("Error: Input file is corrupt or invalid!\n");
        returnCode = -1;
        return returnCode;
    }
    if (oAudioPtr == NULL) {
        printf("Error: Output file is corrupt or invalid!\n");
        returnCode = -2;
        return returnCode;
    }
    if (iNumSamples <= 0) {
        printf("Error: No sample data found!\n");
        returnCode = -3;
        return returnCode;
    }
    if (durationInMS <= 0) {
        printf("Error: Fade duration is too short!\n");
        returnCode = -5;
        return returnCode;
    }
    if (sampleRate < 44100 || sampleRate > 192000) {
        printf("Error: Invalid sample rate!\n");
        returnCode = -6;
        return returnCode;
    }

    float* input = iAudioPtr;
    float* output = oAudioPtr;

    int samples = durationInMS * (sampleRate / 1000);

    if (samples > iNumSamples) {
        printf("Error: Fade duration is too long!\n");
        returnCode = -5;
        return returnCode;
    }

    switch (fadeType) {
    case 1:
        for (int i = 0; i < iNumSamples; i++) {
            if (i < samples) {
                output[i] = input[i] * ((float)i / samples);
            }
            else {
                output[i] = input[i];
            }
        }
        break;
    case 2:
        output[0] = 0;
        for (int i = 1; i < iNumSamples; i++) {
            if (i < samples) {
                output[i] = input[i] * (log10f(i) / log10f(samples));
            }
            else {
                output[i] = input[i];
            }
        }
        break;
    case 3:
        for (int i = 0; i < iNumSamples; i++) {
            if (i < samples) {
                output[i] = input[i] * (-0.5 * cosf((M_PI / samples) * i) + 0.5);
            }
            else {
                output[i] = input[i];
            }
        }
        break;
    default:
        printf("Error: Invalid fade type!\n");
        returnCode = -4;
    }
    return returnCode;
}

int dsp_fadeOut(float* iAudioPtr, int iNumSamples, float* oAudioPtr, int durationInMS, int sampleRate, short fadeType) {
    int returnCode = 0;

    if (iAudioPtr == NULL) {
        printf("Error: Input file is corrupt or invalid!\n");
        returnCode = -1;
        return returnCode;
    }
    if (oAudioPtr == NULL) {
        printf("Error: Output file is corrupt or invalid!\n");
        returnCode = -2;
        return returnCode;
    }
    if (iNumSamples <= 0) {
        printf("Error: No sample data found!\n");
        returnCode = -3;
        return returnCode;
    }
    if (durationInMS <= 0) {
        printf("Error: Fade duration is too short!\n");
        returnCode = -5;
        return returnCode;
    }
    if (sampleRate < 44100 || sampleRate > 192000) {
        printf("Error: Invalid sample rate!\n");
        returnCode = -6;
        return returnCode;
    }

    float* input = iAudioPtr;
    float* output = oAudioPtr;

    int samples = durationInMS * (sampleRate / 1000);
    if (samples > iNumSamples) {
        printf("Error: Fade duration is too long!\n");
        returnCode = -5;
        return returnCode;
    }

    int x = 0;

    switch (fadeType) {
    case 1:
        for (int i = iNumSamples; i > -1; i--) {
            if (i > (iNumSamples - samples)) {
                output[i] = input[i] * ((float)x / samples);
            }
            else {
                output[i] = input[i];
            }
            x++;
        }
        break;
    case 2:
        for (int i = iNumSamples; i > -1; i--) {
            if (i > (iNumSamples - samples)) {
                output[i] = input[i] * (log10f(x) / log10f(samples));
            }
            else {
                output[i] = input[i];
            }
            x++;
        }
        break;
    case 3:
        for (int i = iNumSamples; i > -1; i--) {
            if (i > (iNumSamples - samples)) {
                output[i] = input[i] * (-0.5 * cosf((M_PI / samples) * x) + 0.5);
            }
            else {
                output[i] = input[i];
            }
            x++;
        }
        break;
    default:
        printf("Error: Invalid fade type!\n");
        returnCode = -4;
    }
    return returnCode;
}

int dsp_simpleSquarewave(float* oAudioPtr, int iNumSamples, float freq, float amp, int sampleRate) {
    int returnCode = 0;

    if (oAudioPtr == NULL) {
        printf("Error: Output file is corrupt or invalid!\n");
        returnCode = -1;
        return returnCode;
    }
    if (iNumSamples <= 0) {
        printf("Error: Number of samples is invalid!\n");
        returnCode = -2;
        return returnCode;
    }
    if (freq <= 0 || freq >= (sampleRate / 2)) {
        printf("Error: Frequency is invalid!\n");
        returnCode = -3;
        return returnCode;
    }
    if (amp > 1 || amp == 0 || amp < -1) {
        printf("Error: Amplitude is invalid!\n");
        returnCode = -4;
        return returnCode;
    }
    if (sampleRate < 44100 || sampleRate > 192000) {
        printf("Error: Sample rate is invalid!\n");
        returnCode = -5;
        return returnCode;
    }

    float* output = oAudioPtr;
    float sine = amp;

    for (int i = 0; i < iNumSamples; i++) {
        sine = sinf((float)(freq * (2 * M_PI) * i) / sampleRate);
        output[i] = amp * ((sine > 0) - (sine < 0));
    }
    return returnCode;
}

int dsp_simpleSinewave(float* oAudioPtr, int iNumSamples, float freq, float amp, int sampleRate) {
    int returnCode = 0;

    if (oAudioPtr == NULL) {
        printf("Error: Output file is corrupt or invalid!\n");
        returnCode = -1;
        return returnCode;
    }
    if (iNumSamples <= 0) {
        printf("Error: Number of samples is invalid!\n");
        returnCode = -2;
        return returnCode;
    }
    if (freq <= 0 || freq >= (sampleRate / 2)) {
        printf("Error: Frequency is invalid!\n");
        returnCode = -3;
        return returnCode;
    }
    if (amp > 1 || amp == 0 || amp < -1) {
        printf("Error: Amplitude is invalid!\n");
        returnCode = -4;
        return returnCode;
    }
    if (sampleRate < 44100 || sampleRate > 192000) {
        printf("Error: Sample rate is invalid!\n");
        returnCode = -5;
        return returnCode;
    }

    float* output = oAudioPtr;

    for (int i = 0; i < iNumSamples; i++) {
        output[i] = amp * sinf((float)(freq * (2 * M_PI) * i) / sampleRate);
    }
    return returnCode;
}

int dsp_simpleTrianglewave(float* oAudioPtr, int iNumSamples, float freq, float amp, int sampleRate) {
    int returnCode = 0;

    if (oAudioPtr == NULL) {
        printf("Error: Output file is corrupt or invalid!\n");
        returnCode = -1;
        return returnCode;
    }
    if (iNumSamples <= 0) {
        printf("Error: Number of samples is invalid!\n");
        returnCode = -2;
        return returnCode;
    }
    if (freq <= 0 || freq >= (sampleRate / 2)) {
        printf("Error: Frequency is invalid!\n");
        returnCode = -3;
        return returnCode;
    }
    if (amp > 1 || amp == 0 || amp < -1) {
        printf("Error: Amplitude is invalid!\n");
        returnCode = -4;
        return returnCode;
    }
    if (sampleRate < 44100 || sampleRate > 192000) {
        printf("Error: Sample rate is invalid!\n");
        returnCode = -5;
        return returnCode;
    }

    float* output = oAudioPtr;
    int period = sampleRate / freq;

    for (int i = 0; i < iNumSamples; i++) {
        output[i] = ((4 * amp) / period) * fabs((i % period) - (period / 2)) - amp;
    }
    return returnCode;
}

int dsp_rampSinewave(float* oAudioPtr, int nSamples, float startingFreq, float endingFreq, float gain_dB, int sampleRate) {
    int returnCode = 0;

    if (oAudioPtr == NULL) {
        printf("Error: Output file is corrupt or invalid!\n");
        returnCode = INVALID_OUTPUT_POINTER;
        return returnCode;
    }
    if (nSamples <= 0) {
        printf("Error: Number of samples is invalid!\n");
        returnCode = INVALID_NUM_SAMPLES;
        return returnCode;
    }
    if (startingFreq <= 0 || endingFreq <= 0 || startingFreq >= (sampleRate / 2) || endingFreq >= (sampleRate / 2)) {
        printf("Error: Frequency is invalid!\n");
        returnCode = INVALID_FREQ;
        return returnCode;
    }
    if (gain_dB > 0) {
        printf("Error: Gain must be less than or equal to 0 dB!\n");
        returnCode = INVALID_DB;
        return returnCode;
    }
    if (sampleRate < 44100 || sampleRate > 192000) {
        printf("Error: Sample rate is invalid!\n");
        returnCode = INVALID_SAMPLE_RATE;
        return returnCode;
    }

    double amp = powf(10, (gain_dB / 20));
    double TWO_PI = M_PI * 2;
    double freq_inc = (endingFreq - startingFreq) / nSamples;
    double current_freq = startingFreq;
    double phase_inc = (TWO_PI * current_freq) / sampleRate;
    double phase = 0.0;

    if (oAudioPtr != NULL) {
        float* output = oAudioPtr;
        for (int i = 0; i < nSamples; i++) {
            *output = amp * sin(phase);
            phase_inc = (TWO_PI * current_freq) / sampleRate;
            phase += phase_inc;
            current_freq += freq_inc;
            output++;
        }
    }

    return returnCode;
}

int dsp_additiveSquarewave(float* oAudioPtr, int nSamples, float freq, float gain_dB, int sampleRate) {
    int returnCode = 0;

    // check inputs
    if (oAudioPtr == NULL) {
        printf("Error: Output file is corrupt or invalid!\n");
        returnCode = INVALID_OUTPUT_POINTER;
        return returnCode;
    }
    if (nSamples <= 0) {
        printf("Error: Number of samples is invalid!\n");
        returnCode = INVALID_NUM_SAMPLES;
        return returnCode;
    }
    if (freq <= 0 || freq >= (sampleRate / 2)) {
        printf("Error: Frequency is invalid!\n");
        returnCode = INVALID_FREQ;
        return returnCode;
    }
    if (gain_dB > 0) {
        printf("Error: Gain must be less than or equal to 0 dB!\n");
        returnCode = INVALID_DB;
        return returnCode;
    }
    if (sampleRate < 44100 || sampleRate > 192000) {
        printf("Error: Sample rate is invalid!\n");
        returnCode = INVALID_SAMPLE_RATE;
        return returnCode;
    }

    // initialize helpful variables
    float* output = oAudioPtr;
    double TWO_PI = M_PI * 2;
    double nyquist = ((float)sampleRate / 2);
    double fun_amp = 0.5;

    // initialize harmonic number counter
    int harm_num = 1;

    // initialize harmonic frequency and amplitude variables
    double harm_freq = freq;
    double harm_amp = fun_amp;

    // write fundamental frequency sine wave to output
    for (int i = 0; i < nSamples; i++) {
        output[i] = fun_amp * sin((harm_freq * TWO_PI * i) / sampleRate);
    }

    // increment harmonic number by 2 (odds only starting at 1)
    harm_num += 2;

    // calculate frequency of next odd harmonic
    harm_freq = freq * harm_num;

    // calculate amplitude of next odd harmonic
    harm_amp = harm_amp / harm_num;

    // sum additional harmonics to output until nyquist is reached
    while (harm_freq < nyquist) {
        for (int i = 0; i < nSamples; i++) {
            output[i] += harm_amp * sin((harm_freq * TWO_PI * i) / sampleRate);
        }
        harm_num += 2;
        harm_freq = freq * harm_num;
        harm_amp = (fun_amp / harm_num);
    }

    // normalize output
    returnCode = dsp_normalize(output, nSamples, output, gain_dB);

    return returnCode;
}

int dsp_additiveTrianglewave(float* oAudioPtr, int nSamples, float freq, float gain_dB, int sampleRate) {
    int returnCode = 0;

    // check inputs
    if (oAudioPtr == NULL) {
        printf("Error: Output file is corrupt or invalid!\n");
        returnCode = INVALID_OUTPUT_POINTER;
        return returnCode;
    }
    if (nSamples <= 0) {
        printf("Error: Number of samples is invalid!\n");
        returnCode = INVALID_NUM_SAMPLES;
        return returnCode;
    }
    if (freq <= 0 || freq >= (sampleRate / 2)) {
        printf("Error: Frequency is invalid!\n");
        returnCode = INVALID_FREQ;
        return returnCode;
    }
    if (gain_dB > 0) {
        printf("Error: Gain must be less than or equal to 0 dB!\n");
        returnCode = INVALID_DB;
        return returnCode;
    }
    if (sampleRate < 44100 || sampleRate > 192000) {
        printf("Error: Sample rate is invalid!\n");
        returnCode = INVALID_SAMPLE_RATE;
        return returnCode;
    }

    // initialize helpful variables
    float* output = oAudioPtr;
    double TWO_PI = M_PI * 2;
    double nyquist = ((float)sampleRate / 2);
    double fun_amp = 0.5;

    // initialize harmonic number counter
    int harm_num = 1;

    // initialize harmonic frequency and amplitude variables
    double harm_freq = freq;
    double harm_amp = fun_amp;

    // initialize phase inverter
    int phase = -1;

    // write fundamental frequency sine wave to output
    for (int i = 0; i < nSamples; i++) {
        output[i] = fun_amp * sin((harm_freq * TWO_PI * i) / sampleRate);
    }

    // increment harmonic number by 2 (odds only starting at 1)
    harm_num += 2;

    // calculate frequency of next odd harmonic
    harm_freq = freq * harm_num;

    // calculate amplitude of next odd harmonic
    harm_amp = phase * (fun_amp / pow(harm_num, 2));

    // sum additional harmonics to output until nyquist is reached
    while (harm_freq < nyquist) {
        for (int i = 0; i < nSamples; i++) {
            output[i] += harm_amp * sin((harm_freq * TWO_PI * i) / sampleRate);
        }
        harm_num += 2;
        harm_freq = freq * harm_num;

        // invert phase of every other odd harmonic
        phase *= -1;

        harm_amp = phase * (fun_amp / pow(harm_num, 2));
    }

    // normalize output
    returnCode = dsp_normalize(output, nSamples, output, gain_dB);

    return returnCode;
}

int dsp_tremolo(float* iAudioPtr, int iNumSamples, float* oAudioPtr, float lfoStartRate, float lfoEndRate, float lfoDepth, int sampleRate) {
    int returnCode = 0;

    // check inputs
    if (iAudioPtr == NULL) {
        printf("Error: Input file is corrupt or invalid!\n");
        returnCode = INVALID_INPUT_POINTER;
        return returnCode;
    }
    if (oAudioPtr == NULL) {
        printf("Error: Output file is corrupt or invalid!\n");
        returnCode = INVALID_OUTPUT_POINTER;
        return returnCode;
    }
    if (iNumSamples <= 0) {
        printf("Error: Number of samples is invalid!\n");
        returnCode = INVALID_NUM_SAMPLES;
        return returnCode;
    }
    if (lfoStartRate < 0 || lfoStartRate > 20) {
        printf("Error: LFO start frequency is invalid!\n");
        returnCode = INVALID_FREQ;
        return returnCode;
    }
    if (lfoEndRate < 0 || lfoEndRate > 20) {
        printf("Error: LFO end frequency is invalid!\n");
        returnCode = INVALID_FREQ;
        return returnCode;
    }
    if (lfoDepth < 0 || lfoDepth > 100) {
        printf("Error: LFO depth is invalid! 0-100 range supported.\n");
        returnCode = INVALID_DEPTH;
        return returnCode;
    }
    if (sampleRate != 44100) {
        printf("Error: Sample rate is invalid! Only 44.1 kHz is supported by this function.\n");
        returnCode = INVALID_SAMPLE_RATE;
        return returnCode;
    }

    float* output = oAudioPtr;
    float* input = iAudioPtr;
    double depth = lfoDepth / 100;
    double TWO_PI = M_PI * 2;
    double freq_inc = (lfoEndRate - lfoStartRate) / iNumSamples;
    double current_freq = lfoStartRate;
    double phase_inc = (TWO_PI * current_freq) / sampleRate;
    double phase = 0.0;

    for (int i = 0; i < iNumSamples; i++) {
        output[i] = input[i] - (input[i] * (depth * (0.5 * sin(phase - 1.5) + 0.5)));
        phase_inc = (TWO_PI * current_freq) / sampleRate;
        phase += phase_inc;
        current_freq += freq_inc;
    }

    return returnCode;
}

int dsp_ampModulation(float* iAudioPtr, int iNumSamples, float* oAudioPtr, float modFrequency, float modAmplitude, int sampleRate) {
    int returnCode = 0;

    // check inputs
    if (iAudioPtr == NULL) {
        printf("Error: Input file is corrupt or invalid!\n");
        returnCode = INVALID_INPUT_POINTER;
        return returnCode;
    }
    if (oAudioPtr == NULL) {
        printf("Error: Output file is corrupt or invalid!\n");
        returnCode = INVALID_OUTPUT_POINTER;
        return returnCode;
    }
    if (iNumSamples <= 0) {
        printf("Error: Number of samples is invalid!\n");
        returnCode = INVALID_NUM_SAMPLES;
        return returnCode;
    }
    if (modFrequency <= 0 || modFrequency > (sampleRate / 2)) {
        printf("Error: Modulator frequency is invalid!\n");
        returnCode = INVALID_FREQ;
        return returnCode;
    }
    if (modAmplitude < 0 || modAmplitude > 1) {
        printf("Error: Modulator amplitude is invalid!\n");
        returnCode = INVALID_AMP;
        return returnCode;
    }
    if (sampleRate != 44100) {
        printf("Error: Sample rate is invalid! Only 44.1 kHz is supported by this function.\n");
        returnCode = INVALID_SAMPLE_RATE;
        return returnCode;
    }

    float* output = oAudioPtr;
    float* input = iAudioPtr;
    double TWO_PI = M_PI * 2;

    for (int i = 0; i < iNumSamples; i++) {
        output[i] = input[i] * (modAmplitude * -(0.5 * sin((modFrequency * TWO_PI * i) / sampleRate) + 0.5) + 1);
    }

    return returnCode;
}

int dsp_ringModulation(float* iAudioPtr, int iNumSamples, float* oAudioPtr, float modFrequency, float modAmplitude, int sampleRate) {
    int returnCode = 0;

    // check inputs
    if (iAudioPtr == NULL) {
        printf("Error: Input file is corrupt or invalid!\n");
        returnCode = INVALID_INPUT_POINTER;
        return returnCode;
    }
    if (oAudioPtr == NULL) {
        printf("Error: Output file is corrupt or invalid!\n");
        returnCode = INVALID_OUTPUT_POINTER;
        return returnCode;
    }
    if (iNumSamples <= 0) {
        printf("Error: Number of samples is invalid!\n");
        returnCode = INVALID_NUM_SAMPLES;
        return returnCode;
    }
    if (modFrequency <= 0 || modFrequency >= (sampleRate / 2)) {
        printf("Error: Modulator frequency is invalid!\n");
        returnCode = INVALID_FREQ;
        return returnCode;
    }
    if (modAmplitude < 0 || modAmplitude > 1) {
        printf("Error: Modulator amplitude is invalid!\n");
        returnCode = INVALID_AMP;
        return returnCode;
    }
    if (sampleRate != 44100) {
        printf("Error: Sample rate is invalid! Only 44.1 kHz is supported by this function.\n");
        returnCode = INVALID_SAMPLE_RATE;
        return returnCode;
    }

    float* output = oAudioPtr;
    float* input = iAudioPtr;
    double TWO_PI = M_PI * 2;

    for (int i = 0; i < iNumSamples; i++) {
        output[i] = input[i] * (modAmplitude * sin((modFrequency * TWO_PI * i) / sampleRate));
    }

    return returnCode;
}

int dsp_flanger(float* iAudioPtr, int iNumSamples, float* oAudioPtr, float modFrequency, float modAmplitude, int sampleRate) {
    int returnCode = 0;

    // check inputs
    if (iAudioPtr == NULL) {
        printf("Error: Input file is corrupt or invalid!\n");
        returnCode = INVALID_INPUT_POINTER;
        return returnCode;
    }
    if (oAudioPtr == NULL) {
        printf("Error: Output file is corrupt or invalid!\n");
        returnCode = INVALID_OUTPUT_POINTER;
        return returnCode;
    }
    if (iNumSamples <= 0) {
        printf("Error: Number of samples is invalid!\n");
        returnCode = INVALID_NUM_SAMPLES;
        return returnCode;
    }
    if (modFrequency <= 0 || modFrequency > 20) {
        printf("Error: Modulator frequency is invalid!\n");
        returnCode = INVALID_FREQ;
        return returnCode;
    }
    if (modAmplitude <= 0 || modAmplitude > 1) {
        printf("Error: Modulator amplitude is invalid!\n");
        returnCode = INVALID_AMP;
        return returnCode;
    }
    if (sampleRate != 44100) {
        printf("Error: Sample rate is invalid! Only 44.1 kHz is supported by this function.\n");
        returnCode = INVALID_SAMPLE_RATE;
        return returnCode;
    }

    float* output = oAudioPtr;
    float* input = iAudioPtr;

    // convert modAmplitude to a sample offset, used in the LFO
    double delayMS = modAmplitude * 20;
    double sampleOffset = 44.1 * delayMS;
    double halfSampleOffset = sampleOffset / 2;

    // initalize helpful PI related variables
    double TWO_PI = M_PI * 2;
    double HALF_PI = M_PI / 2;

    // initalize variable for the oscillating shift amount, as well as other variables used for linear interpolation
    double shift = 0;
    int lower = 0;
    int upper = 0;
    double between = 0;

    for (int i = 0; i < iNumSamples; i++) {
        if (i > sampleOffset) {

            // calculate current shift amount in number of samples
            shift = halfSampleOffset * sin((modFrequency * TWO_PI * i) / sampleRate) + halfSampleOffset;

            // calculate quantities needed for linear interpolation
            lower = floor(shift);
            upper = ceil(shift);
            between = shift - floor(shift);
        }
        else {
            shift = 0;
        }

        // sum the processed signal with the original input signal
        // lerp() is a helper function that calculates an interpolated amplitude value at the current shift amount
        output[i] = input[i] + lerp(input[i - lower], input[i - upper], between);
    }

    // normalize the output to avoid clipping
    returnCode = dsp_normalize(output, iNumSamples, output, 0);

    return returnCode;
}

double sinc(double x) {
    if (x == 0) {
        return 1;
    }
    else {
        return sin(x) / x;
    }
}

float* convolve(float* h, double* x, int lenH, int lenX, int lenY) {
    int i, j, h_start, x_start, x_end;
    float* y = (float*)calloc(lenY, sizeof(float));
    for (i = 0; i < lenY; i++)
    {
        h_start = fmin(i, lenH - 1);
        x_start = fmax(0, i - lenH + 1);
        x_end = fmin(i + 1, lenX);
        for (j = x_start; j < x_end; j++)
        {
            y[i] += h[h_start--] * x[j];
        }
    }
    return y;
}

int normalize_double(double* input, int samples, double* output, float threshold) {
    int error = 0;
    double max = 0;
    double absolute;
    for (int i = 0; i < samples; i++) {
        absolute = abs(input[i]);
        if (absolute > max) {
            max = absolute;
        }
    }
    double scalar = threshold / max;
    for (int i = 0; i < samples; i++) {
        output[i] = input[i] * scalar;
    }
    return error;
}

int fir_lpf(float* iAudioPtr, int iNumSamples, float* oAudioPtr, int oNumSamples, double cutoff, int firSize, int sampleRate) {
    int returnCode = 0;

    // check inputs
    if (iAudioPtr == NULL) {
        printf("Error: Input file is corrupt or invalid!\n");
        returnCode = INVALID_INPUT_POINTER;
        return returnCode;
    }
    if (oAudioPtr == NULL) {
        printf("Error: Output file is corrupt or invalid!\n");
        returnCode = INVALID_OUTPUT_POINTER;
        return returnCode;
    }
    if (iNumSamples <= 0) {
        printf("Error: Number of samples is invalid!\n");
        returnCode = INVALID_NUM_SAMPLES;
        return returnCode;
    }
    if (cutoff > sampleRate / 2) {
        printf("Error: Invalid cutoff frequency!\n");
        returnCode = INVALID_FREQ;
        return returnCode;
    }
    if (sampleRate != 44100) {
        printf("Error: Sample rate is invalid! Only 44.1 kHz is supported by this function.\n");
        returnCode = INVALID_SAMPLE_RATE;
        return returnCode;
    }

    float* input = iAudioPtr;
    float* output = oAudioPtr;

    double fc = cutoff / sampleRate;

    double* impulse_response = (double*)malloc(firSize * sizeof(double));

    returnCode = normalize_double(impulse_response, firSize, impulse_response, 1);

    for (int i = 0; i < firSize; i++) {
        impulse_response[i] = sinc(2 * M_PI * fc * (i - (firSize / 2))) * (0.42 - (0.5 * cos((2 * M_PI * i) / firSize)) + (0.08 * cos((2 * M_PI * i) / firSize)));
    }

    int output_buffer_size = iNumSamples + firSize - 1;

    float* convolved = convolve(input, impulse_response, iNumSamples, firSize, output_buffer_size);

    for (int i = 0; i < output_buffer_size; i++) {
        output[i] = convolved[i];
    }

    free(convolved);

    returnCode = dsp_normalize(output, iNumSamples, output, 0);
    return returnCode;
}

//..................................................................................................... createIRfromFFT_BSF
void createIRfromFFT_BSF(int firSize, int nyquist, float lower_cutoff, float upper_cutoff, std::vector<double>& a_resp, std::vector<double>& p_resp)
{
    int frameSize = firSize;
    int numBins = frameSize / 2;
    float binWidth = nyquist / numBins;

    std::vector<std::complex<float>> inputFrame(frameSize);
    std::vector<std::complex<float>> outputFrame(frameSize);
    int i;

    // fill complex frames with real and imag parts
    // we're creating spectral data in the shape of an 'ideal brickwall filter'
    float binEdge = 0;
    for (i = 0; i < numBins; i++)
    {
        binEdge = (i + 1) * binWidth; // the upper edge of this bin, in Hz
        inputFrame[i].imag(0.0);
        inputFrame[frameSize - 1 - i].imag(0.0);
        if (binEdge < lower_cutoff || binEdge > upper_cutoff)
        {
            inputFrame[i].real(1.0);
            inputFrame[frameSize - 1 - i].real(1.0); 
        }
        else
        {
            inputFrame[i].real(0.0);
            inputFrame[frameSize -1 - i].real(0.0);
        }
    }

    int order = log10(firSize) / log10(2);

    // perform the IFFT
    juce::dsp::FFT fft(order); // the (#) is fftOrder. It's the exponent that 2 is raised to. 2^7 = 128, 2^8 = 256, etc...
    fft.perform(inputFrame.data(), outputFrame.data(), true); // this is an IFFT. When the last param is false it's an FFT.

    // pack amplitude and phase into arrays
    for (i = 0; i < firSize; i++)
    {
        a_resp.push_back(outputFrame[i].real());
        p_resp.push_back(outputFrame[i].imag());
    }
}

//..................................................................................................... fir_bandstop
int fir_bandstop(float* iAudioPtr, int iNumSamples, float* oAudioPtr, int oNumSamples, float lowerEdgeInHz, float upperEdgeInHz, int firSize, int sampleRate)
{
    int returnCode = 0;

    // check inputs
    if (iAudioPtr == NULL) {
        printf("Error: Input file is corrupt or invalid!\n");
        returnCode = INVALID_INPUT_POINTER;
        return returnCode;
    }
    if (oAudioPtr == NULL) {
        printf("Error: Output file is corrupt or invalid!\n");
        returnCode = INVALID_OUTPUT_POINTER;
        return returnCode;
    }
    if (iNumSamples <= 0) {
        printf("Error: Number of samples is invalid!\n");
        returnCode = INVALID_NUM_SAMPLES;
        return returnCode;
    }
    if (lowerEdgeInHz > (sampleRate / 2) || upperEdgeInHz > (sampleRate / 2)) {
        printf("Error: Invalid cutoff frequency!\n");
        returnCode = INVALID_FREQ;
        return returnCode;
    }
    if (lowerEdgeInHz > upperEdgeInHz) {
        printf("Error: Invalid cutoff frequency!\n");
        returnCode = INVALID_FREQ;
        return returnCode;
    }
    if (lowerEdgeInHz <= 0 || upperEdgeInHz <= 0) {
        printf("Error: Invalid cutoff frequency!\n");
        returnCode = INVALID_FREQ;
        return returnCode;
    }
    if (sampleRate != 44100) {
        printf("Error: Sample rate is invalid! Only 44.1 kHz is supported by this function.\n");
        returnCode = INVALID_SAMPLE_RATE;
        return returnCode;
    }

    int i;

    // create IR from FFT for band stop
    std::vector<double> a_resp;
    std::vector<double> p_resp;
    createIRfromFFT_BSF(firSize, sampleRate / 2, lowerEdgeInHz, upperEdgeInHz, a_resp, p_resp);

    // mirror IR around nyquist
    float* fftIR = (float*)malloc(firSize * sizeof(float));

    for (i = 0; i < firSize; i++)
    {
        if (i < firSize / 2)
            fftIR[i] = a_resp.at(i + (firSize / 2));
        else
            fftIR[i] = a_resp.at(i - (firSize / 2));
    }

    // create window
    double* Blackman = (double*)malloc(firSize * sizeof(double));
    double f = 0;
    for (i = 0; i < firSize; i++)
    {
        f = (2 * M_PI * i) / (firSize);
        Blackman[i] = 0.42 - 0.5 * cos(f) + 0.08 * cos(2 * f);
    }

    // window the IR
    double* WIR = (double*)malloc(firSize * sizeof(double));
    for (i = 0; i < firSize; i++)
    {
        WIR[i] = fftIR[i] * Blackman[i];
    }

    // convolve WIR with input
    float* convolved = convolve(iAudioPtr, WIR, iNumSamples, firSize, oNumSamples);

    // write result to output
    for (int i = 0; i < oNumSamples; i++) {
        oAudioPtr[i] = convolved[i];
    }

    // normalize output
    returnCode = dsp_normalize(oAudioPtr, iNumSamples, oAudioPtr, 0);

    free(fftIR);
    free(Blackman);
    free(WIR);
    free(convolved);

    return returnCode;
}

int iir_4stage_lpf(float* iAudioPtr, int iNumSamples, float* oAudioPtr, int oNumSamples, double cutoff, int sampleRate) {
    int returnCode = 0;

    // check inputs
    if (iAudioPtr == NULL) {
        printf("Error: Input file is corrupt or invalid!\n");
        returnCode = INVALID_INPUT_POINTER;
        return returnCode;
    }
    if (oAudioPtr == NULL) {
        printf("Error: Output file is corrupt or invalid!\n");
        returnCode = INVALID_OUTPUT_POINTER;
        return returnCode;
    }
    if (iNumSamples <= 0) {
        printf("Error: Number of samples is invalid!\n");
        returnCode = INVALID_NUM_SAMPLES;
        return returnCode;
    }
    if (cutoff > sampleRate / 2) {
        printf("Error: Invalid cutoff frequency!\n");
        returnCode = INVALID_FREQ;
        return returnCode;
    }
    if (sampleRate != 44100) {
        printf("Error: Sample rate is invalid! Only 44.1 kHz is supported by this function.\n");
        returnCode = INVALID_SAMPLE_RATE;
        return returnCode;
    }

    double fc = cutoff / sampleRate;
    double x = exp(-14.445 * fc);
    double a = pow((1 - x), 4);
    double b1 = 4 * x;
    double b2 = -6 * pow(x, 2);
    double b3 = 4 * pow(x, 3);
    double b4 = -1 * pow(x, 4);

    for (int i = 0; i < 4; i++) {
        oAudioPtr[i] = 0.0;
    }

    for (int i = 4; i < iNumSamples; i++) {
        oAudioPtr[i] = (a * iAudioPtr[i]) + (b1 * oAudioPtr[i - 1]) + (b2 * oAudioPtr[i - 2]) + (b3 * oAudioPtr[i - 3]) + (b4 * oAudioPtr[i - 4]);
    }
    return returnCode;
}

int iir_bpf(float* iAudioPtr, int iNumSamples, float* oAudioPtr, int oNumSamples, double center, double bandwidth, int sampleRate) {
    int returnCode = 0;

    // check inputs
    if (iAudioPtr == NULL) {
        printf("Error: Input file is corrupt or invalid!\n");
        returnCode = INVALID_INPUT_POINTER;
        return returnCode;
    }
    if (oAudioPtr == NULL) {
        printf("Error: Output file is corrupt or invalid!\n");
        returnCode = INVALID_OUTPUT_POINTER;
        return returnCode;
    }
    if (iNumSamples <= 0) {
        printf("Error: Number of samples is invalid!\n");
        returnCode = INVALID_NUM_SAMPLES;
        return returnCode;
    }
    if (center > sampleRate / 2 || bandwidth > sampleRate / 2) {
        printf("Error: Invalid cutoff frequency!\n");
        returnCode = INVALID_FREQ;
        return returnCode;
    }
    if (sampleRate != 44100) {
        printf("Error: Sample rate is invalid! Only 44.1 kHz is supported by this function.\n");
        returnCode = INVALID_SAMPLE_RATE;
        return returnCode;
    }

    double fc = center / sampleRate;
    double bw = bandwidth / sampleRate;

    double r = 1 - (3 * bw);
    double k = (1 - (2 * r * cos(2 * M_PI * fc)) + pow(r, 2)) / (2 - (2 * cos(2 * M_PI * fc)));

    double a0 = 1 - k;
    double a1 = -2 * (k - r) * cos(2 * M_PI * fc);
    double a2 = pow(r, 2) - k;

    double b1 = 2 * r * cos(2 * M_PI * fc);
    double b2 = -1 * pow(r, 2);

    returnCode = dsp_normalize(iAudioPtr, iNumSamples, iAudioPtr, 0);

    for (int i = 0; i < 2; i++) {
        oAudioPtr[i] = 0.0;
    }

    for (int i = 2; i < iNumSamples; i++) {
        oAudioPtr[i] = (a0 * iAudioPtr[i]) + (a1 * iAudioPtr[i - 1]) + (a2 * iAudioPtr[i - 2]) + (b1 * oAudioPtr[i - 1]) + (b2 * oAudioPtr[i - 2]);
    }

    returnCode = dsp_normalize(oAudioPtr, oNumSamples, oAudioPtr, 0);

    return returnCode;
}

int chebyshev_subroutine(short filterType, double cf, double ripple, int numPoles, short currentPole, double* coeffArray) {
    int returnCode = 0;

    double RP = -cos((M_PI / (numPoles * 2)) + (currentPole - 1) * (M_PI / numPoles));
    double IP = sin((M_PI / (numPoles * 2)) + (currentPole - 1) * (M_PI / numPoles));
    double ES, VX, KX, K;

    if (ripple != 0) {
        ES = sqrt(pow(100.0 / (100 - ripple), 2) - 1);
        VX = (1.0 / numPoles) * log((1.0 / ES) + sqrt(pow(1.0 / ES, 2) + 1));
        KX = (1.0 / numPoles) * log((1.0 / ES) + sqrt(pow(1.0 / ES, 2) - 1));
        KX = (exp(KX) + exp(-KX)) / 2;
        RP = RP * ((exp(VX) - exp(-VX)) / 2) / KX;
        IP = IP * ((exp(VX) + exp(-VX)) / 2) / KX;
    }

    double T = 2 * tan(0.5);
    double W = 2 * M_PI * cf;
    double M = pow(RP, 2) + pow(IP, 2);
    double D = 4 - 4 * RP * T + M * pow(T, 2);
    double x0 = pow(T, 2) / D;
    double x1 = 2 * pow(T, 2) / D;
    double x2 = pow(T, 2) / D;
    double y1 = (8 - 2 * M * pow(T, 2)) / D;
    double y2 = (-4 - 4 * RP * T - M * pow(T, 2)) / D;

    if (filterType == 0) {
        K = sin(0.5 - (W / 2)) / sin(0.5 + (W / 2));
    }
    else if (filterType == 1) {
        K = -cos((W / 2) + 0.5) / cos((W / 2) - 0.5);
    }
    else {

    }

    D = 1 + (y1 * K) - y2 * pow(K, 2);
    coeffArray[0] = (x0 - x1 * K + x2 * pow(K, 2)) / D;
    coeffArray[1] = (-2 * x0 * K + x1 + x1 * pow(K, 2) - 2 * x2 * K) / D;
    coeffArray[2] = (x0 * pow(K, 2) - x1 * K + x2) / D;
    coeffArray[3] = (2 * K + y1 + y1 * pow(K, 2) - 2 * y2 * K) / D;
    coeffArray[4] = (-pow(K, 2) - y1 * K + y2) / D;

    if (filterType == 1) {
        coeffArray[1] = -1 * coeffArray[1];
        coeffArray[3] = -1 * coeffArray[3];
    }

    return returnCode;
}

int chebyshev_get_coeffs(short filterType, double cf, double ripple, int numPoles, double* a_coeffs, double* b_coeffs) {
    int returnCode = 0;

    double TA[23];
    double TB[23];

    for (int i = 0; i < 23; i++) {
        a_coeffs[i] = 0.0;
        b_coeffs[i] = 0.0;
    }

    a_coeffs[2] = 1.0;
    b_coeffs[2] = 1.0;

    for (int p = 1; p < numPoles / 2; p++) {
        double temp[5];

        chebyshev_subroutine(filterType, cf, ripple, numPoles, p, temp);

        for (int i = 0; i < 23; i++) {
            TA[i] = a_coeffs[i];
            TB[i] = b_coeffs[i];
        }

        for (int i = 2; i < 23; i++) {
            a_coeffs[i] = (temp[0] * TA[i]) + (temp[1] * TA[i - 1]) + (temp[2] * TA[i - 2]);
            b_coeffs[i] = TB[i] - (temp[3] * TB[i - 1]) - (temp[4] * TB[i - 2]);
        }
    }

    b_coeffs[2] = 0.0;

    for (int i = 0; i < 21; i++) {
        a_coeffs[i] = a_coeffs[i + 2];
        b_coeffs[i] = -b_coeffs[i + 2];
    }

    double SA = 0;
    double SB = 0;

    for (int i = 0; i < 21; i++) {
        if (filterType == 0) {
            SA = SA + a_coeffs[i];
            SB = SB + b_coeffs[i];
        }
        else if (filterType == 1) {
            SA = SA + (a_coeffs[i] * pow(-1, i));
            SB = SB + (b_coeffs[i] * pow(-1, i));
        }
        else {

        }
    }

    double gain = SA / (1 - SB);

    for (int i = 0; i < 21; i++) {
        a_coeffs[i] = a_coeffs[i] / gain;
    }

    return returnCode;
}

int chebyshev_filter(float* iAudioPtr, int iNumSamples, float* oAudioPtr, short filterType, double freq, double ripple, int numPoles, int sampleRate) {
    int returnCode = 0;

    // check inputs
    if (iAudioPtr == NULL) {
        printf("Error: Input file is corrupt or invalid!\n");
        returnCode = INVALID_INPUT_POINTER;
        return returnCode;
    }
    if (oAudioPtr == NULL) {
        printf("Error: Output file is corrupt or invalid!\n");
        returnCode = INVALID_OUTPUT_POINTER;
        return returnCode;
    }
    if (iNumSamples <= 0) {
        printf("Error: Number of samples is invalid!\n");
        returnCode = INVALID_NUM_SAMPLES;
        return returnCode;
    }
    if (sampleRate != 44100) {
        printf("Error: Sample rate is invalid! Only 44.1 kHz is supported by this function.\n");
        returnCode = INVALID_SAMPLE_RATE;
        return returnCode;
    }
    if (ripple < 0 || ripple > 29) {
        printf("Error: Invalid ripple percentage!\n");
        returnCode = INVALID_RIPPLE;
        return returnCode;
    }
    if (numPoles < 0 || numPoles > 20 || numPoles % 2 != 0) {
        printf("Error: Invalid number of poles!\n");
        returnCode = INVALID_POLES;
        return returnCode;
    }
    if (filterType < 0 || filterType > 1) {
        printf("Error: Invalid filter type!\n");
        returnCode = INVALID_FILTER;
        return returnCode;
    }

    double cf = freq / sampleRate;

    if (cf < 0 || cf > 0.5) {
        printf("Error: Invalid cutoff frequency!\n");
        returnCode = INVALID_FREQ;
        return returnCode;
    }

    double a_coeffs[23];
    double b_coeffs[23];
    double a;
    double b;

    chebyshev_get_coeffs(filterType, cf, ripple, numPoles, a_coeffs, b_coeffs);

    for (int i = 0; i < 22; i++) {
        oAudioPtr[i] = 0;
    }

    for (int i = 22; i < iNumSamples; i++) {
        a = 0;
        b = 0;

        a = a + (a_coeffs[0] * iAudioPtr[i]);
        for (int j = 1; j < 23; j++) {
            a = a + (a_coeffs[j] * iAudioPtr[i - j]);
            b = b + (b_coeffs[j] * oAudioPtr[i - j]);
        }
        oAudioPtr[i] = a + b;
    }

    returnCode = dsp_normalize(oAudioPtr, iNumSamples, oAudioPtr, dsp_get_gain(iAudioPtr, iNumSamples));

    return returnCode;
}

int dspa_pitchChange(float* iAudioPtr, int iNumSamples, float* oAudioPtr, int oNumSamples, int numerator, int denominator, int sampleRate) {
    int returnCode = 0;

    // check inputs
    if (iAudioPtr == NULL) {
        printf("Error: Input file is corrupt or invalid!\n");
        returnCode = INVALID_INPUT_POINTER;
        return returnCode;
    }
    if (oAudioPtr == NULL) {
        printf("Error: Output file is corrupt or invalid!\n");
        returnCode = INVALID_OUTPUT_POINTER;
        return returnCode;
    }
    if (iNumSamples <= 0) {
        printf("Error: Number of samples is invalid!\n");
        returnCode = INVALID_NUM_SAMPLES;
        return returnCode;
    }
    if (sampleRate != 44100) {
        printf("Error: Sample rate is invalid! Only 44.1 kHz is supported by this function.\n");
        returnCode = INVALID_SAMPLE_RATE;
        return returnCode;
    }
    if (((double)numerator / (double)denominator) < 0.5 || ((double)numerator / (double)denominator) > 2) {
        printf("Error: Invalid duration ratio!\n");
        returnCode = INVALID_DURATION;
        return returnCode;
    }

    double nyquist = sampleRate / 2;
    double lpf_cutoff = nyquist;

    if (numerator > denominator) {
        lpf_cutoff = nyquist / (double)numerator;
    }
    else {
        lpf_cutoff = nyquist / (double)denominator;
    }

    int num_stuff = numerator * iNumSamples;
    float* stuffer_buffer = (float*)malloc(num_stuff * sizeof(float));

    int i = 0;
    for (int j = 0; j < num_stuff; j++) {
        if (j % numerator == 0) {
            stuffer_buffer[j] = iAudioPtr[i];
            i++;
        }
        else {
            stuffer_buffer[j] = 0;
        }
    }

    float* filter_buffer = (float*)malloc(num_stuff * sizeof(float));

    returnCode = fir_lpf(stuffer_buffer, num_stuff, filter_buffer, num_stuff, lpf_cutoff, 64, sampleRate);

    i = 0;
    for (int j = 0; j < num_stuff; j++) {
        if (j % denominator == 0) {
            oAudioPtr[i] = filter_buffer[j];
            i++;
        }
    }

    return returnCode;
}