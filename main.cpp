#include <cstddef>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "ps2000.h"
#include "ps3000aApi.h"
#include <math.h>
#include <rfftw.h>
#include <ftw.h>
#include <sndfile.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

// defaults
#define DEFAULT_METHOD 2
#define SAMPLES_PER_SECOND 1000000
#define FRAME_SIZE_IN_SEC 1
#define FREQPARAL 0
#define FREQPARAH 2
#define MAX_SHIFT_IN_NSEC 10000
#define WND_BRD_IN_SEC 0.1
#define SIGNAL_GEN_VOLTAGE_IN_mV 2000
#define MNUMBERPERFREQ 10

#define PS3000_BUFFER_SIZE 100000
static int16_t ps3000_driver_buffer[2][PS3000_BUFFER_SIZE];

typedef struct
{
	double c0;
	double c1;
	double c2;
	double c3;
} cspline_para_t;

static int g_method = DEFAULT_METHOD;
static int g_voltage_range = PS2000_2V;
static long g_totalSamples = 0;
static int16_t g_overflow = 0;
static SNDFILE* g_wavfile = NULL;
static unsigned int g_signal_generator_voltage_in_mV = SIGNAL_GEN_VOLTAGE_IN_mV;
static long g_samples_per_second = SAMPLES_PER_SECOND;
static long g_frame_size_in_seconds = FRAME_SIZE_IN_SEC;
static long g_total_sample_nbr = 0;
static int g_maxn = MNUMBERPERFREQ;
static double g_filter_para_low = FREQPARAL;
static double g_filter_para_high = FREQPARAH;
static unsigned int g_max_shift_in_nsec = MAX_SHIFT_IN_NSEC;
static int g_freqs_default[] =
{
	1000, 1252, 1568, 1964, 2460, 3080, 3857, 4831, 6050, 7576, 9488,
	11882, 14880, 18634, 23336, 29224, 36598, 45833, 57397
};

static void convolution(double* s1, double* s2, double* rs, int n, int cn)
{
	int m = n - 2 * cn;

	for (int k = -cn; k <= cn; k++)
	{
		rs[k + cn] = 0;

		for (int i = 0; i < m; i++)
		{
			rs[k + cn] += s1[cn + i] * s2[cn + i + k];
		}
	}
}

static void cubic_spline(double* s, int n, cspline_para_t* p)
{
	typedef struct
	{
		double* c;
		double* g;
		double* diag;
		double* offdiag;
	} cspline_state_t;

	gsl_spline *spline_ptr = gsl_spline_alloc(gsl_interp_cspline, n);

	double* x_array = new double [n];

	for (int i = 0; i < n; i++) x_array[i] = i;

	gsl_spline_init(spline_ptr, x_array, s, n);

	const cspline_state_t* state = (const cspline_state_t*)(spline_ptr->interp->state);

	for (int i = 0; i < n - 1; i++)
	{
		p[i].c0 = s[i];
		double y_hi = s[i + 1];
		double dy = y_hi - p[i].c0;
		double c_ip1 = state->c[i + 1];
		p[i].c2 = state->c[i];
		p[i].c1 = dy - (c_ip1 + 2.0 * p[i].c2) / 3.0;
		p[i].c3 = (c_ip1 - p[i].c2) / 3.0;
	}

	gsl_spline_free(spline_ptr);

	delete [] x_array;
}

static double find_argmax(double* s, int n)
{
	cspline_para_t* p = new cspline_para_t [n];

	cubic_spline(s, n, p);

	double maxy = s[0];
	double maxx = 0;

	for (int i = 0; i < n - 1; i++)
	{
		double vs = p[i].c2 * p[i].c2 - 3 * p[i].c1 * p[i].c3;

		if (vs >= 0)
		{
			for (int j = 0; j < 2; j++)
			{
				double x = (-p[i].c2 - (2 * j - 1) * sqrt(vs)) / (3 * p[i].c3);

				if ((x >= 0) && (x <= 1))
				{
					double y = p[i].c0 + p[i].c1 * x + p[i].c2 * x * x + p[i].c3 * x * x * x;

					if (y > maxy)
					{
						maxy = y;
						maxx = i + x;
					}
				}
			}
		}

		double x = - p[i].c1 / (2 * p[i].c2);

		if ((x >= 0) && (x <= 1))
		{
			double y = p[i].c0 + p[i].c1 * x + p[i].c2 * x * x;

			if (y > maxy)
			{
				maxy = y;
				maxx = i + x;
			}
		}

		if (s[i] > maxy)
		{
			maxy = s[i];
			maxx = i;
		}

		if (s[i + 1] > maxy)
		{
			maxy = s[i + 1];
			maxx = i + 1;
		}
	}

	delete [] p;

	return maxx;
}

static void normalize(double* s, int n)
{
	double min = s[0];
	double max = s[0];

	for (int i = 0; i < n; i++)
	{
		if (s[i] < min) min = s[i];

		if (s[i] > max) max = s[i];
	}

	double c0, c1;

	c0 = -(max + min) / (max - min);
	c1 = 2.0 / (max - min);

	for (int i = 0; i < n; i++)
	{
		s[i] = c0 + c1 * s[i];
	}
}

void write_debug_wav(const char* filename, double* s, int n, int sr)
{
	SF_INFO sfinfo;

	sfinfo.channels = 1;
	sfinfo.format = SF_FORMAT_WAV | SF_FORMAT_DOUBLE;
	sfinfo.samplerate = sr;
	SNDFILE* sndfp = sf_open(filename, SFM_WRITE, &sfinfo);

	if (sndfp == NULL)
	{
		fprintf(stderr, "error: could not create file '%s'\n", filename);
		return;
	}

	sf_write_double(sndfp, s, n);

	sf_close(sndfp);
}

static void close_wav(SNDFILE** sndf)
{
	if ((sndf != NULL) && (*sndf != NULL)) sf_close(*sndf);

	*sndf = 0;
}

static SNDFILE* open_wav(const char* filename, int sr)
{
	SF_INFO sfinfo;

	if (filename == 0)
	{
		fprintf(stderr, "error: no output file given\n");
		return NULL;
	}

	sfinfo.channels = 2;
	sfinfo.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;
	sfinfo.samplerate = sr;
	SNDFILE* sndf = sf_open(filename, SFM_WRITE, &sfinfo);

	return sndf;
}

void write_wav(const char* filename, int sr, int n, double* s1, double* s2)
{
	SNDFILE* sndf = open_wav(filename, sr);

	if (sndf == NULL)
	{
		fprintf(stderr, "error: could not create file '%s'\n", filename);
		return;
	}

	short* buffer = new short [2 * n];

	for (int i = 0; i < n; i++)
	{
		buffer[i * 2] = (short)(s1[i] * 0x7fff);
		buffer[i * 2 + 1] = (short)(s2[i] * 0x7fff);
	}

	sf_write_short(sndf, buffer, 2 * n);

	delete [] buffer;

	close_wav(&sndf);
}

void PREF4 PS2000_StreamingReadyCallback(int16_t** overviewBuffers,
		int16_t overflow, uint32_t triggeredAt, int16_t triggered,
		int16_t auto_stop, uint32_t nValues)
{
	g_totalSamples += nValues;
	g_overflow = overflow;

	short* buffer = new short [2 * nValues];

	for (uint32_t i = 0; i < nValues; i++)
	{
		buffer[2 * i + 0] = overviewBuffers[0][i];
		buffer[2 * i + 1] = overviewBuffers[2][i];
	}

	sf_write_short(g_wavfile, buffer, 2 * nValues);

	delete [] buffer;
}

void PREF4 PS3000_StreamingReadyCallback(int16_t handle, int32_t noOfSamples,
		uint32_t startIndex, int16_t overflow,
		uint32_t triggerAt, int16_t triggered,
		int16_t autoStop, void* pParameter)
{
	g_totalSamples += noOfSamples;
	g_overflow = overflow;

	int16_t* buffer = new int16_t [2 * noOfSamples];

	for (int32_t i = 0; i < noOfSamples; i++)
	{
		// Channel B seems to be shifted compared to channel A.
		// This becomes clear when both channels of the wave file are filled with
		// data from the same channel. The calculated shift is then almost zero.
		// This shows that the algorithm works.
		buffer[2 * i + 0] = ps3000_driver_buffer[0][startIndex + i];
		buffer[2 * i + 1] = ps3000_driver_buffer[1][startIndex + i];
	}

	sf_write_short(g_wavfile, buffer, 2 * noOfSamples);

	delete [] buffer;
}

static void read_wav(char* name, SF_INFO* sfinfo, double** s1, double** s2)
{
	SNDFILE* sndf = sf_open(name, SFM_READ, sfinfo);

	if (!sfinfo->channels)
	{
		fprintf(stderr, "error: could not open file '%s'\n", name);
		exit(EXIT_FAILURE);
	}

	if (sfinfo->channels != 2)
	{
		fprintf(stderr, "error: wav-file must have two channels\n");
		exit(EXIT_FAILURE);
	}

	if (sfinfo->format != (SF_FORMAT_WAV | SF_FORMAT_PCM_16))
	{
		fprintf(stderr, "error: wav-file format is invalid\n");
		exit(EXIT_FAILURE);
	}

	if ((s1 != NULL) && (s2 != NULL))
	{
		short* buffer = new short [2 * sfinfo->frames];
		sf_count_t rb = sf_readf_short(sndf, buffer, sfinfo->frames);

		sf_close(sndf);

		if (rb != sfinfo->frames)
		{
			fprintf(stderr, "error: could not read input file\n");
			exit(EXIT_FAILURE);
		}

		*s1 = new double [sfinfo->frames];
		*s2 = new double [sfinfo->frames];

		for (int i = 0; i < sfinfo->frames; i++)
		{
			(*s1)[i] = buffer[i * 2];
			(*s2)[i] = buffer[i * 2 + 1];
			(*s1)[i] /= 0x7fff;
			(*s2)[i] /= 0x7fff;
		}

		delete [] buffer;
	}
}

static void find_fft_max(double* buf, int n, double sr, double* freq, double* phase)
{
	fftw_complex* in = new fftw_complex [n];
	fftw_complex* spec = new fftw_complex [n];
	rfftw_plan p;

	for (int i = 0; i < n; i++)
	{
		in[i].re = buf[i];
		in[i].im = 0;
		spec[i].re = 0;
		spec[i].im = 0;
	}

	p = fftw_create_plan(n, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_one(p, in, spec);
	fftw_destroy_plan(p);

	double maxa = -1;
	int maxi = -1;

	for (int i = 1; i < n / 2; i++)
	{
		double a = spec[i].re * spec[i].re + spec[i].im * spec[i].im;

		if (maxa <= a)
		{
			maxa = a;
			maxi = i;
		}
	}

	double freqres = sr / n;
	*freq = maxi * freqres;

	double a = spec[maxi].re;
	double b = spec[maxi].im;
	double r = sqrt(a * a + b * b);
	*phase = 0;

	if (b >= 0)
	{
		*phase = acos(a / r);
	}
	else
	{
		*phase = 2 * M_PI - acos(a / r);
	}

	delete [] in;
	delete [] spec;
}

static void filter_buf(double* buf, int n, int fl, int fh)
{
	fftw_complex* in = new fftw_complex [n];
	fftw_complex* spec = new fftw_complex [n];
	rfftw_plan p;

	for (int i = 0; i < n; i++)
	{
		in[i].re = buf[i];
		in[i].im = 0;
		spec[i].re = 0;
		spec[i].im = 0;
	}

	p = fftw_create_plan(n, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_one(p, in, spec);
	fftw_destroy_plan(p);

	spec[0].re = 0;
	spec[0].im = 0;
	spec[n / 2].re = 0;
	spec[n / 2].im = 0;

	for (int i = 1; i < n / 2; i++)
	{
		if ((i > fh) || (i < fl))
		{
			spec[i].re = 0;
			spec[i].im = 0;
			spec[n - i].re = 0;
			spec[n - i].im = 0;
		}
	}

	p = fftw_create_plan(n, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_one(p, spec, in);
	fftw_destroy_plan(p);

	for (int i = 0; i < n; i++) buf[i] = in[i].re / n;

	delete [] in;
	delete [] spec;
}

static double wnd_func(int i, int n, int m)
{
	double f;

	if ((m == 0) || (n < m) || (i >= n))
	{
		fprintf(stderr, "error: invalid parameters in function %s\n", __func__);
		exit(EXIT_FAILURE);
	}

	if (i <= m)
	{
		f = 0.5 * (1.0 - cos((2 * M_PI * i) / (2 * m)));
	}
	else if (i >= n - m)
	{
		f = 0.5 * (1.0 - cos((2 * M_PI * (i - (n - m) + m)) / (2 * m)));
	}
	else
	{
		f = 1.0;
	}

	return f;
}

static void full_filter(double* data, int n, int fl, int fh, int m)
{
	// apply window
	for (int i = 0; i < n; i++)
	{
		data[i] = wnd_func(i, n, m) * data[i];
	}

	// process filter
	filter_buf(data, n, fl, fh);
}

static void full_filter_wav(char* src, char* dest, int cfreq_low, int cfreq_high)
{
	double* s1;
	double* s2;
	SF_INFO sfinfo = {0};
	read_wav(src, &sfinfo, &s1, &s2);

	int fl = (int)((1.0 * sfinfo.frames * cfreq_low) / sfinfo.samplerate);
	int fh = (int)((1.0 * sfinfo.frames * cfreq_high) / sfinfo.samplerate);
	int m = (int)(1.0 * WND_BRD_IN_SEC * sfinfo.samplerate);

	if ((fl > fh) || (fh - fl < 10))
	{
		fprintf(stderr, "error: invalid filter frequencies\n");
		exit(EXIT_FAILURE);
	}

	full_filter(s1, sfinfo.frames, fl, fh, m);
	full_filter(s2, sfinfo.frames, fl, fh, m);
	write_wav(dest, sfinfo.samplerate, sfinfo.frames -  2 * m, &(s1[m]), &(s2[m]));
	delete [] s1;
	delete [] s2;
}

static void calc_shift(char* name, double* shift, int freq, int n)
{
	if (g_method == 1)
	{
		// FFT filter
		char filter_name [0x200];
		snprintf(filter_name, sizeof(filter_name), "./debug/%i_%i_filter.wav", freq, n);
		full_filter_wav(name, filter_name, freq * g_filter_para_low, freq * g_filter_para_high);

		double* s1;
		double* s2;
		SF_INFO sfinfo = {0};
		read_wav(filter_name, &sfinfo, &s1, &s2);

		// calculate convolution
		int cn = (int)(g_max_shift_in_nsec * 1e-9 * sfinfo.samplerate);
		int cnbs = 2 * cn + 1;
		double* conv = new double [cnbs];

		convolution(s1, s2, conv, sfinfo.frames, cn);
		normalize(conv, cnbs);

		char conv_name [0x200];
		snprintf(conv_name, sizeof(conv_name), "./debug/%i_%i_conv.wav", freq, n);
		write_debug_wav(conv_name, conv, cnbs, sfinfo.samplerate);

		*shift = 1e9 * (find_argmax(conv, cnbs) - cn) / sfinfo.samplerate;

		delete [] conv;
		delete [] s1;
		delete [] s2;
	}
	else if (g_method == 2)
	{
		double freq1, phase1, freq2, phase2;
		double* s1;
		double* s2;
		SF_INFO sfinfo = {0};
		read_wav(name, &sfinfo, &s1, &s2);

		// apply window
		int m = (int)(1.0 * WND_BRD_IN_SEC * sfinfo.samplerate);

		for (int i = 0; i < n; i++)
		{
			s1[i] = wnd_func(i, sfinfo.frames, m) * s1[i];
			s2[i] = wnd_func(i, sfinfo.frames, m) * s2[i];
		}

		// calculate FFT and find phases for the maximum amplitude
		find_fft_max(s1, sfinfo.frames, sfinfo.samplerate, &freq1, &phase1);
		find_fft_max(s2, sfinfo.frames, sfinfo.samplerate, &freq2, &phase2);

		if (freq1 != freq2)
		{
			fprintf(stderr, "error: channel A has frequency %f and channel A frequency: %f\n", freq1, freq2);
		}

		double dphase = phase1 - phase2;

		if (dphase > M_PI) dphase -= 2 * M_PI;

		if (dphase < -M_PI) dphase += 2 * M_PI;

		*shift = dphase / (2 * M_PI * freq1) * 1e9;

		delete [] s1;
		delete [] s2;
	}
	else
	{
		fprintf(stderr, "error: method not supported\n");
		exit(EXIT_FAILURE);
	}
}

static void perform_measurement_ps3000(int16_t handle, int freq, char* name)
{
	if (!handle)
	{
		fprintf(stderr, "error: invalid handle\n");
		exit(EXIT_FAILURE);
	}

	// we use only the both first channels
	if (ps3000aSetChannel(handle, PS3000A_CHANNEL_A, true, PS3000A_DC, (PS3000A_RANGE)g_voltage_range, 0) != PICO_OK)
	{
		fprintf(stderr, "error: could not initialize channel A\n");
		exit(EXIT_FAILURE);
	}

	if (ps3000aSetChannel(handle, PS3000A_CHANNEL_B, true, PS3000A_DC, (PS3000A_RANGE)g_voltage_range, 0) != PICO_OK)
	{
		fprintf(stderr, "error: could not initialize channel B\n");
		exit(EXIT_FAILURE);
	}

	// disable trigger
	if (ps3000aSetSimpleTrigger(handle, false, PS3000A_CHANNEL_A, 0, PS3000A_ABOVE, 0, 0) != PICO_OK)
	{
		fprintf(stderr, "error: could not disable trigger\n");
		exit(EXIT_FAILURE);
	}

	// setup signal generator
	if (ps3000aSetSigGenBuiltInV2(handle, 0, g_signal_generator_voltage_in_mV * 1000, PS3000A_SINE,
								  (double) freq, (double) freq,	0, 0,
								  (PS3000A_SWEEP_TYPE)0, (PS3000A_EXTRA_OPERATIONS)0, 0, 0,
								  (PS3000A_SIGGEN_TRIG_TYPE)0, (PS3000A_SIGGEN_TRIG_SOURCE)0, 0) != PICO_OK)
	{
		fprintf(stderr, "error: could not initialize signal generator\n");
		exit(EXIT_FAILURE);
	}

	if (ps3000aSetDataBuffer(handle, PS3000A_CHANNEL_A, ps3000_driver_buffer[0],
							 PS3000_BUFFER_SIZE, 0, PS3000A_RATIO_MODE_NONE) != PICO_OK)
	{
		fprintf(stderr, "error: could not set data buffer for channel A\n");
		exit(EXIT_FAILURE);
	}

	if (ps3000aSetDataBuffer(handle, PS3000A_CHANNEL_B, ps3000_driver_buffer[1],
							 PS3000_BUFFER_SIZE, 0, PS3000A_RATIO_MODE_NONE) != PICO_OK)
	{
		fprintf(stderr, "error: could not set data buffer for channel B\n");
		exit(EXIT_FAILURE);
	}

	uint32_t sample_time_in_ns_in = 1000000000 / g_samples_per_second;
	uint32_t sample_time_in_ns =  sample_time_in_ns_in;

	if (ps3000aRunStreaming(handle, &sample_time_in_ns, PS3000A_NS, 0,
							2 * g_total_sample_nbr, 0, 1,
							PS3000A_RATIO_MODE_NONE, PS3000_BUFFER_SIZE) != PICO_OK)
	{
		fprintf(stderr, "error: could not start streaming\n");
		exit(EXIT_FAILURE);
	}

	g_total_sample_nbr = (long)(1000000000.0 / sample_time_in_ns * g_frame_size_in_seconds);

	g_wavfile = open_wav(name, (int)(1000000000.0 / sample_time_in_ns));

	if (g_wavfile == NULL)
	{
		fprintf(stderr, "error: could not create wavfile '%s'\n", name);
		exit(EXIT_FAILURE);
	}

	g_totalSamples = 0;

	while (g_totalSamples < g_total_sample_nbr)
	{
		ps3000aGetStreamingLatestValues(handle, PS3000_StreamingReadyCallback, NULL);
	}

	ps3000aStop(handle);

	ps3000aCloseUnit(handle);

	close_wav(&g_wavfile);

	if (g_overflow)
	{
		fprintf(stderr, "error: there was an overflow\n");
		exit(EXIT_FAILURE);
	}
}

static void perform_measurement_ps2000(int16_t handle, int freq, char* name)
{
	if (!handle)
	{
		fprintf(stderr, "error: invalid handle\n");
		exit(EXIT_FAILURE);
	}

	// disable ETS
	if (ps2000_set_ets(handle, PS2000_ETS_OFF, 0, 0))
	{
		fprintf(stderr, "error: could not disable ETS mode\n");
		exit(EXIT_FAILURE);
	}

	if (!ps2000_set_sig_gen_built_in(handle,
									 0, g_signal_generator_voltage_in_mV * 1000,
									 PS2000_SINE, freq, freq,
									 0, 0, PS2000_UPDOWN, 0))
	{
		printf("could not initialize sinus sweep\n");
		exit(EXIT_FAILURE);
	}

	// setup channels
	if (!ps2000_set_channel(handle, 0, true, true, (PS2000_RANGE)g_voltage_range))
	{
		fprintf(stderr, "error: setup of channel 1 failed\n");
		exit(EXIT_FAILURE);
	}

	if (!ps2000_set_channel(handle, 1, true, true, (PS2000_RANGE)g_voltage_range))
	{
		fprintf(stderr, "error: setup of channel 2 failed\n");
		exit(EXIT_FAILURE);
	}

	// disable trigger
	if (!ps2000_set_trigger(handle, PS2000_NONE, 0, 0, 0, 0))
	{
		fprintf(stderr, "error: could not disable trigger\n");
		exit(EXIT_FAILURE);
	}

	g_wavfile = open_wav(name, g_samples_per_second);

	if (g_wavfile == NULL)
	{
		fprintf(stderr, "error: could not create wavfile '%s'\n", name);
		exit(EXIT_FAILURE);
	}

	// setup sampling rate
	g_totalSamples = 0;

	if (!ps2000_run_streaming_ns(handle, 1000000000 / g_samples_per_second,
								 PS2000_NS, g_samples_per_second, false,
								 1, g_samples_per_second))
	{
		fprintf(stderr, "error: ps2000_run_streaming_ns failed\n");
		exit(EXIT_FAILURE);
	}

	while (g_totalSamples < g_total_sample_nbr)
	{
		ps2000_get_streaming_last_values(handle, PS2000_StreamingReadyCallback);
	}

	ps2000_stop(handle);

	ps2000_close_unit(handle);

	close_wav(&g_wavfile);

	if (g_overflow)
	{
		fprintf(stderr, "error: there was an overflow\n");
		exit(EXIT_FAILURE);
	}
}

static void perform_measurement(int freq, char* name)
{
	int16_t handle;

	PICO_STATUS status = ps3000aOpenUnit(&handle, NULL);

	if (status != PICO_NOT_FOUND)
	{
		if (status == PICO_USB3_0_DEVICE_NON_USB3_0_PORT)
		{
			fprintf(stderr, "error: PicoScope 3000 is not connected with an USB3 port\n");
			exit(EXIT_FAILURE);
		}

		if (status != PICO_OK)
		{
			fprintf(stderr, "error: could not open PicoScope 3000 because of error 0x%08lx\n", (uint64_t)status);
			exit(EXIT_FAILURE);
		}

		perform_measurement_ps3000(handle, freq, name);
		return;
	}

	handle = ps2000_open_unit();

	if (handle)
	{
		perform_measurement_ps2000(handle, freq, name);
		return;
	}

	fprintf(stderr, "error: no device\n");
	exit(EXIT_FAILURE);
}

static int remove_cb(const char *fpath, const struct stat *sb, int typeFlag, struct FTW *ftwbuf)
{
	if (ftwbuf->level) remove(fpath);

	return 0;
}

static void clean_up_debug(void)
{
	nftw("./debug", remove_cb, 10, FTW_DEPTH);
}

static void process(int* freqs, int freq_number, int maxn)
{
	char wavename [0x200];
	double shift;

	for (int i = 0; i < freq_number; i++)
	{
		for (int j = 0; j < maxn; j++)
		{
			snprintf(wavename, sizeof(wavename), "./debug/%i_%i_signal.wav",  freqs[i], (int)j);
			perform_measurement(freqs[i], wavename);
			calc_shift(wavename, &shift, freqs[i], j);
			FILE* pf = fopen("./debug/result.txt", "a");

			if (!((i == 0) && (j == 0))) fprintf(pf, ",");

			fprintf(pf, "{%i,%f}", freqs[i], shift);
			printf("%i, freq=%i, shift=%f ns\n", j, freqs[i], shift);
			fclose(pf);
		}
	}

	FILE* pf = fopen("./debug/result.txt", "a");
	fprintf(pf, "};\n");
	fclose(pf);
}

static void print_help()
{
	printf("available options:\n");
	printf("s - samplerate (default: %i)\n", SAMPLES_PER_SECOND);
	printf("d - total measurement time in seconds (default: %i)\n", FRAME_SIZE_IN_SEC);
	printf("f - signal generator frequency (default is an array or frequencies)\n");
	printf("n - number of measurements per frequency (default: %i)\n", MNUMBERPERFREQ);
	printf("v - scope voltage range (2V=%i, 20V=%i, default: %i)\n", PS2000_2V, PS2000_20V, PS2000_2V);
	printf("m - maximum expected shift in nanoseconds (default: %i)\n", MAX_SHIFT_IN_NSEC);
	printf("g - signal generator voltage in mV (default: %i)\n", SIGNAL_GEN_VOLTAGE_IN_mV);
	printf("l - filter parameter. lower border = p*frequency (default: %1.2f)\n", (float)FREQPARAL);
	printf("u - filter parameter. upper border = p*frequency (default: %1.2f)\n", (float)FREQPARAH);
	printf("M - phase estimation method (default: %i)\n", DEFAULT_METHOD);
	printf("C - string for documentation\n");
}

int main(int argc, char *argv[])
{
	int opt;
	int* freqs = g_freqs_default;
	int freq_number = sizeof(g_freqs_default) / sizeof(g_freqs_default[0]);
	unsigned int maxn = g_maxn;

	while ((opt = getopt(argc, argv, "s:d:f:n:v:p:m:g:l:u:M:C:h")) != -1)
	{
		switch (opt)
		{
			case 'h':
				print_help();
				exit(EXIT_SUCCESS);

			case 's':
				if (sscanf(optarg, "%li", &g_samples_per_second) != 1)
				{
					fprintf(stderr, "error: invalid samplerate %s\n", optarg);
					exit(EXIT_FAILURE);
				}

				break;

			case 'M':
				if (sscanf(optarg, "%i", &g_method) != 1)
				{
					fprintf(stderr, "error: invalid method parameter %s\n", optarg);
					exit(EXIT_FAILURE);
				}

				break;

			case 'd':
				if (sscanf(optarg, "%li", &g_frame_size_in_seconds) != 1)
				{
					fprintf(stderr, "error: invalid frame size %s\n", optarg);
					exit(EXIT_FAILURE);
				}

				break;

			case 'f':
				int frequency;

				if ((sscanf(optarg, "%i", &frequency) != 1) || (frequency < 0) || (frequency > 100000))
				{
					fprintf(stderr, "error: invalid frequency %s\n", optarg);
					exit(EXIT_FAILURE);
				}
				else
				{
					freqs[0] = frequency;
					freq_number = 1;
				}

				break;

			case 'n':
				if (sscanf(optarg, "%ui", &maxn) != 1)
				{
					fprintf(stderr, "error: invalid number of measurements per frequency %s\n", optarg);
					exit(EXIT_FAILURE);
				}

				break;

			case 'm':
				if (sscanf(optarg, "%ui", &g_max_shift_in_nsec) != 1)
				{
					fprintf(stderr, "error: invalid value for maximum expected shift in nanoseconds %s\n", optarg);
					exit(EXIT_FAILURE);
				}

				break;

			case 'g':
				if (sscanf(optarg, "%ui", &g_signal_generator_voltage_in_mV) != 1)
				{
					fprintf(stderr, "error: invalid value for signal generator voltage %s\n", optarg);
					exit(EXIT_FAILURE);
				}

				break;

			case 'l':
				if ((sscanf(optarg, "%lf", &g_filter_para_low) != 1) || (g_filter_para_low < 0) || (g_filter_para_low >= 1))
				{
					fprintf(stderr, "error: invalid filter parameter %s\n", optarg);
					exit(EXIT_FAILURE);
				}

				break;

			case 'u':
				if ((sscanf(optarg, "%lf", &g_filter_para_high) != 1) || (g_filter_para_high <= 1))
				{
					fprintf(stderr, "error: invalid filter parameter %s\n", optarg);
					exit(EXIT_FAILURE);
				}

				break;

			case 'v':
				if (sscanf(optarg, "%i", &g_voltage_range) != 1)
				{
					fprintf(stderr, "error: invalid voltage range %s\n", optarg);
					exit(EXIT_FAILURE);
				}

				break;

			case 'C':
				break;

			default:
				fprintf(stderr, "error: invalid arguments\n");
				print_help();
				exit(EXIT_FAILURE);
		}
	}

	g_total_sample_nbr = g_samples_per_second * g_frame_size_in_seconds;

	if ((g_total_sample_nbr > 200000000) || (g_total_sample_nbr == 0))
	{
		fprintf(stderr, "error: invalid total sample number (%li)\n", g_total_sample_nbr);
		exit(EXIT_FAILURE);
	}

	if (g_samples_per_second < 1000000)
	{
		fprintf(stderr, "error: sample rate %li is too low\n", g_samples_per_second);
		exit(EXIT_FAILURE);
	}

	struct stat st = {0};

	clean_up_debug();

	if (stat("./debug", &st) == -1) mkdir("./debug", 0700);

	FILE* pf = fopen("./debug/result.txt", "a");
	fprintf(pf, "(* ");

	for (int i = 0; i < argc; i++)
	{
		fprintf(pf, "%s ", argv[i]);
	}

	fprintf(pf, "*)\n");
	fprintf(pf, "\ndata = {");
	fclose(pf);

	process(freqs, freq_number, maxn);

	return EXIT_SUCCESS;
}
