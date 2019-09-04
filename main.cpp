/*
    Program for controlling a PicoScope 2204A to measure the phase shift 
    between the two inputs.
	
    Copyright (C) 2019  Steffen Kühn / steffen.kuehn@quantino-theory.de

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <cstddef>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <termios.h>
#include <sys/ioctl.h>
#include <unistd.h>
#include <sndfile.h>
#include <time.h>
#include "ps2000.h"
#include <sndfile.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <rfftw.h>

#define FRAME_SIZE_IN_SEC 60.0
#define SAMPLES_PER_SECOND 1000000
#define CHANNELS 2
#define MEASUREMENT_TIME_IN_MSEC 65000
#define MAX_SHIFT_IN_SEC 1e-5

//#define MUSIC

#ifdef MUSIC
#define VOLTAGE_RANGE PS2000_20V
const int maxn = 1;
const int freqs[] = {5000};
#else
#define VOLTAGE_RANGE PS2000_2V
const int maxn = 5;
const int freqs[] =
{
	1000, 1252, 1568, 1964, 2460, 3080, 3857, 4831, 6050, 7576, 9488,
	11882, 14880, 18634, 23336, 29224, 36598, 45833, 57397
};
#endif

typedef struct
{
	double c0;
	double c1;
	double c2;
	double c3;
} cspline_para_t;

static uint32_t g_totalSamples = 0;
static int16_t g_overflow = 0;
static SNDFILE* g_wavfile = NULL;

static double average(double* b, int n)
{
	double v = 0;

	for (int i = 0; i < n; i++)
	{
		v += b[i];
	}

	return v / n;
}

static double standarddev(double* b, int n)
{
	double a = average(b, n);
	double v = 0;

	for (int i = 0; i < n; i++)
	{
		v += b[i] * b[i];
	}

	return sqrt(v / n - a * a);
}

static void convolution(double* s1, double* s2, int j, double* rs, int n, int cn)
{
	for (int k = -cn; k <= cn; k++)
	{
		rs[k + cn] = 0;

		for (int i = 0; i < n; i++)
		{
			rs[k + cn] += s1[j + i] * s2[j + i + k];
		}

		rs[k + cn] /= n;
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
		printf("error: could not create file '%s'\n", filename);
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
		printf("error: no output file given\n");
		return NULL;
	}

	sfinfo.channels = CHANNELS;
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
		printf("error: could not create file '%s'\n", filename);
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

void PREF4 StreamingReadyCallback(int16_t** overviewBuffers,
								  int16_t overflow, uint32_t triggeredAt, int16_t triggered,
								  int16_t auto_stop, uint32_t nValues)
{
	g_totalSamples += nValues;
	g_overflow = overflow;

	short* buffer = new short [2 * nValues];

	for (uint32_t i = 0; i < nValues; i++)
	{
		buffer[2 * i] = overviewBuffers[0][i];
		buffer[2 * i + 1] = overviewBuffers[2][i];
	}

	sf_write_short(g_wavfile, buffer, 2 * nValues);

	delete [] buffer;
}

static uint32_t get_mclock_in_ms()
{
	struct timespec ts;
	uint32_t tick = 0;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	tick  = ts.tv_nsec / 1000000;
	tick += ts.tv_sec * 1000;
	return tick;
}

static void scope_measurement(int freq, int nbr, uint32_t time_in_ms, char* name)
{
	// open picoscope 2204A
	int16_t handle = ps2000_open_unit();

	if (!handle)
	{
		printf("error: unable to open device\n");
		exit(1);
	}

	// disable ETS
	if (ps2000_set_ets(handle, PS2000_ETS_OFF, 0, 0))
	{
		printf("error: could not disable ETS mode\n");
		exit(1);
	}

	if (!ps2000_set_sig_gen_built_in(handle,
									 0, 1800000, // 1.8 volt
									 PS2000_SINE, freq, freq,
									 0, 0, PS2000_UPDOWN, 0))
	{
		printf("could not initialize sinus sweep\n");
		exit(1);
	}

	// setup channels
	if (!ps2000_set_channel(handle, 0, true, true, VOLTAGE_RANGE))
	{
		printf("error: setup of channel 1 failed\n");
		exit(1);
	}

	if (!ps2000_set_channel(handle, 1, true, true, VOLTAGE_RANGE))
	{
		printf("error: setup of channel 2 failed\n");
		exit(1);
	}

	// disable trigger
	if (!ps2000_set_trigger(handle, PS2000_NONE, 0, 0, 0, 0))
	{
		printf("error: could not disable trigger\n");
		exit(1);
	}

	// setup sampling rate SAMPLES_PER_SECOND
	if (!ps2000_run_streaming_ns(handle, 1000000000 / SAMPLES_PER_SECOND, PS2000_NS, SAMPLES_PER_SECOND, false, 1, SAMPLES_PER_SECOND))
	{
		printf("error: ps2000_run_streaming_ns failed\n");
		exit(1);
	}

	g_wavfile = open_wav(name, SAMPLES_PER_SECOND);

	if (g_wavfile == NULL)
	{
		printf("could not create wavfile '%s'\n", name);
		exit(1);
	}

	uint32_t start_time = get_mclock_in_ms();

	for (;;)
	{
		ps2000_get_streaming_last_values(handle, StreamingReadyCallback);

		if (get_mclock_in_ms() - start_time > time_in_ms) break;
	}

	ps2000_stop(handle);

	ps2000_close_unit(handle);

	close_wav(&g_wavfile);

	if (g_overflow)
	{
		printf("error: there was an overflow\n");
		exit(1);
	}
}

static void read_wav(char* name, SF_INFO* sfinfo, double** s1, double** s2)
{
	SNDFILE* sndf = sf_open(name, SFM_READ, sfinfo);

	if (!sfinfo->channels)
	{
		printf("error: could not open file '%s'\n", name);
		exit(1);
	}

	if (sfinfo->channels != 2)
	{
		printf("error: wav-file must have two channels\n");
		exit(1);
	}

	if (sfinfo->format != (SF_FORMAT_WAV | SF_FORMAT_PCM_16))
	{
		printf("error: wav-file format is invalid\n");
		exit(1);
	}

	short* buffer = new short [2 * sfinfo->frames];
	sf_count_t rb = sf_readf_short(sndf, buffer, sfinfo->frames);

	sf_close(sndf);

	if (rb != sfinfo->frames)
	{
		printf("error: could not read input file\n");
		exit(1);
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

static void calc_shift(char* name, double* avg, double* std, int freq)
{
	double* s1;
	double* s2;
	SF_INFO sfinfo = {0};
	read_wav(name, &sfinfo, &s1, &s2);

	int frame_size = (int)(FRAME_SIZE_IN_SEC * sfinfo.samplerate);
	int step_size = (int)(FRAME_SIZE_IN_SEC * 0.5 * sfinfo.samplerate);

	int cn = (int)(MAX_SHIFT_IN_SEC * sfinfo.samplerate);
	int cnbs = 2 * cn + 1;
	double* conv = new double [cnbs];

	int n = 0;

	for (int i = cn; i + frame_size < sfinfo.frames - cn; i += step_size) n++;

	double* mpos = new double [n];

	for (int i = 0; i < n; i++)
	{
		convolution(s1, s2, i * step_size + cn, conv, frame_size, cn);
		normalize(conv, cnbs);
		char dbgn [0x200];
		snprintf(dbgn, sizeof(dbgn), "./debug/%i_%i_conv.wav", i, freq);
		write_debug_wav(dbgn, conv, cnbs, sfinfo.samplerate);
		mpos[i] = 1e9 * (find_argmax(conv, cnbs) - cn) / sfinfo.samplerate;
	}

	*avg = average(mpos, n);
	*std = standarddev(mpos, n);

	delete [] conv;
	delete [] s1;
	delete [] s2;
	delete [] mpos;
}

static void filter_buf(double* buf, int n, int co)
{
	fftw_complex in[n];
	fftw_complex spec[n];
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

	for (int i = 1; i <= n / 2; i++)
	{
		if (i > co)
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
}

void filter(double* data, int n, int sr, int cofreq)
{
	double* input = new double [n];
	memcpy(input, data, n * sizeof(double));
	memset(data, 0, n * sizeof(double));
	const int bs = 8192;
	double wnd [bs];
	int co = (int)((1.0 * bs * cofreq) / sr);

	for (int i = 0; i < n - bs; i += bs / 2)
	{
		for (int j = 0; j < bs; j++)
		{
			double f;

			if (j < bs / 2)
			{
				f = (2.0 * j) / bs;
			}
			else
			{
				f = (2 - (2.0 * j) / bs);
			}

			wnd[j] = f * input[i + j];
		}

		filter_buf(wnd, bs, co);

		for (int j = 0; j < bs; j++)
		{
			data[i + j] += wnd[j];
		}
	}

	delete [] input;
}

void filter_wav(char* src, char* dest, int cofreq)
{
	double* s1;
	double* s2;
	SF_INFO sfinfo = {0};
	read_wav(src, &sfinfo, &s1, &s2);
	filter(s1, sfinfo.frames, sfinfo.samplerate, cofreq);
	filter(s2, sfinfo.frames, sfinfo.samplerate, cofreq);
	write_wav(dest, sfinfo.samplerate, sfinfo.frames, s1, s2);
	delete [] s1;
	delete [] s2;
}

int main(void)
{
	char measurements [0x200];
	char filter [0x200];
	double avg;
	double std;

	FILE* pf = fopen("result.txt", "w");
	fclose(pf);

	int n = 0;

	for (int j = 0; j < maxn; j++)
	{
		for (size_t i = 0; i < sizeof(freqs) / sizeof(freqs[0]); i++)
		{
			snprintf(measurements, sizeof(measurements), "measurements_%i_%i_%i.wav", (int)j, (int)i, freqs[i]);
			snprintf(filter, sizeof(filter), "filter_%i_%i_%i.wav", (int)j, (int)i, freqs[i]);
			scope_measurement(freqs[i], n, MEASUREMENT_TIME_IN_MSEC, measurements);
			filter_wav(measurements, filter, 2 * freqs[i]);
			calc_shift(filter, &avg, &std, freqs[i]);
			FILE* pf = fopen("result.txt", "a");
			fprintf(pf, "{%i,%f,%f},\n", freqs[i], avg, std);
			printf("%i: freq=%i, shift=%f ns, std=%f ns\n", n + 1, freqs[i], avg, std);
			fclose(pf);
			n++;
		}
	}
}
