extern "C" {
#include <libavcodec/avcodec.h>
#include <libavcodec/avfft.h>
#include <libavformat/avformat.h>
#include <libavutil/time.h>
#include <libswresample/swresample.h>
}

#include <stdio.h>
#include <complex>

using std::complex;

typedef struct AudioParams {
	int freq;
	int channels;
	int64_t channel_layout;
	enum AVSampleFormat fmt;
} AudioParams;

AVCodec *codec;
AVCodecContext *codec_ctx;
AVFormatContext *format_ctx;
AudioParams src, tgt;
AVStream *audio_st;
AVPacket pkt, pkt1;
AVFrame frame1;
uint8_t buf[288000], *buf1;
unsigned int buf_size, buf1_size;
#define frame_size 2048
#define rdft_bits 9
#define hop_size 256
#define ar_order 30
#define warp_order 200
#define warp_c 0.63
#define max_frames 100000
#define num_chan 51
#define num_sub 5
const double sub_weigh[num_sub+1] = {0, 1, 0.5, 0.333333, 0.25, 0.2};
SwrContext *swr_ctx;
double audio_clock;
const double pi = 3.14159265358979324;
int frame_cnt = 1;

double data[frame_size], data1[frame_size];
FFTSample data2[frame_size];
RDFTContext *rdft_ctx;
complex<double> c[2*ar_order+1], r[2*ar_order], r2[2*ar_order];
double dp[max_frames][ar_order*num_sub+num_chan];
int dpt[max_frames][ar_order*num_sub+num_chan];
typedef struct Result {
	double freq, amp;
} Result;
Result res[max_frames][ar_order*num_sub];
int reslen[max_frames], track[max_frames];
double ts[max_frames];
double b[ar_order+1];
double d[ar_order+1][warp_order];

int cmp(const void* a, const void* b) {
	return ((Result*)a)->amp > ((Result*)b)->amp ? -1 : 1;
}

void analyze_root(int order) {
	int i, j;
	for (i = 0; i < order; i++) {
		int flag = 0;
		if (isnan(r[i].real())) flag = 1;
		for (j = i+1; !flag && j < order; j++)
			if (std::abs(r[i] - r[j]) < 1e-8) flag = 1;
		if (flag) r2[i] = std::exp(complex<double>(0, 2 * pi * i/ order));
		else r2[i] = r[i];
	}
	while (1) {
		for (i = 0; i < order; i++) {
			complex<double> x = 0;
			complex<double> y = 1;
			for (j = order; j >= 0; j--)
				x = x * r2[i] + c[j];
			for (j = i+1; j < order; j++) y *= r2[i] - r2[j];
			for (j = 0; j < i; j++) y *= r2[i] - r[j];
			r[i] = r2[i] - x / y;
		}
		//double max = 0;
		int flag = 0;
		for (i = 0; i < order; i++) {
			if (std::abs(r2[i] - r[i]) > 1e-8) flag = 1;
			//if (cabs(r2[i] - r[i]) > max) max = cabs(r2[i] - r[i]);
		}
		//printf("%le\n", max);
		if (!flag) break;
		for (i = 0; i < order; i++) r2[i] = r[i];
	}
}

double loudness(double amp, double freq) {
	if (amp < 1e-5) amp = 1;
	amp = log(amp);
	double relfreq = fabs(log(freq / 2000.0));
	amp = amp * (1 + 0.017 * pow(relfreq, 2.5)) - 0.0745 * pow(relfreq, 3.5);
	if (amp < 0) amp = 0;
	return amp;
}

double a[ar_order];
void analyze() {
	//printf("%lf: ", audio_clock);
	int i, j, k;
	for (j = 0; j < frame_size; j++)
		data1[j] = data[j];
	for (j = 0; j < frame_size; j++)
		data1[j] = data1[j] * exp(2.0 * j / frame_size - 2.0) * (0.5 - 0.5 * cos(2 * pi * j / (frame_size - 1)));
	/*
	double x[warp_order];
	if (a[0] < 0 || a[0] >= 0) {
		for (j = 0; j < warp_order; j++) {
			x[j] = 0;
			for (k = 0; k < ar_order; k++) x[j] += a[k] * d[k][j];
		}
		for (j = frame_size-1; j >= 0; j--)
			for (k = 0; k < warp_order && k <= j; k++)
				data1[j] -= data1[j - k] * x[k] * 0.8;
	}
	*/

	ts[frame_cnt] = audio_clock;
	double s, g1[warp_order], g[ar_order+1], f[ar_order][ar_order+1];
	s = 0;
	for (j = 0; j < frame_size; j++)
		s += data1[j];
	s /= frame_size;
	for (j = 0; j < frame_size; j++)
		data1[j] -= s;
	for (i = 0; i < warp_order; i++) {
		g1[i] = 0;
		for (j = 0; j < frame_size-i; j++)
			g1[i] += data1[j] * data1[j+i];
	}
	for (i = 0; i <= ar_order; i++) {
		g[i] = 0;
		for (j = 0; j < warp_order; j++)
			g[i] += d[i][j] * g1[j];
	}
	for (i = 0; i < ar_order; i++) {
		for (j = 0; j < ar_order; j++)
			f[i][j] = g[abs(i-j)];
		f[i][ar_order] = g[i+1];
	}
	/*
	for (i = 0; i < ar_order; i++)
		for (j = 0; j <= ar_order; j++)
			f[i][j] = 0;
	for (i = ar_order-1; i < frame_size-1; i++) {
		int hi = i + 1;
		for (j = 0; j < ar_order; j++) {
			int hj = i - j;
			f[j][ar_order] += data1[hi] * data1[hj];
			for (k = j; k < ar_order; k++)
				f[j][k] += data1[hj] * data1[i-k];
		}
	}
	for (i = 0; i < ar_order; i++)
		for (j = i; j <= ar_order; j++) {
			f[i][j] /= frame_size - ar_order;
			if (j < ar_order) f[j][i] = f[i][j];
		}
	*/
	for (i = 0; i < ar_order; i++) {
		int maxi = i;
		for (j = i+1; j < ar_order; j++)
			if (fabs(f[j][i]) > fabs(f[maxi][i])) maxi = j;
		for (j = i; j <= ar_order; j++) {
			double tmp = f[i][j]; f[i][j] = f[maxi][j], f[maxi][j] = tmp;
		}
		for (j = i+1; j < ar_order; j++)
			for (k = ar_order; k >= i; k--)
				f[j][k] -= f[j][i] / f[i][i] * f[i][k];
	}
	for (i = ar_order-1; i >= 0; i--) {
		a[i] = f[i][ar_order];
		for (j = i+1; j < ar_order; j++)
			a[i] -= f[i][j] * a[j];
		a[i] /= f[i][i];
	}
	b[0] = 1;
	for (i = 0; i < ar_order; i++) {
		b[0] += a[i] * a[i];
		b[i+1] = -2 * a[i];
	}
	for (i = 0; i < ar_order; i++)
		for (j = i+1; j < ar_order; j++)
			b[j-i] += 2 * a[i] * a[j];
	for (i = 1; i <= ar_order; i++) c[ar_order-i] = i * b[i];
	c[ar_order] = 0;
	for (i = 1; i <= ar_order; i++) c[ar_order+i] = -i * b[i];
	for (i = 0; i <= 2*ar_order; i++) c[i] /= c[2*ar_order];
	analyze_root(2*ar_order);

	double gain = g[0];
	for (i = 0; i < ar_order; i++) gain -= a[i] * g[i+1];
	gain = sqrt(gain);

	/*if (fabs(audio_clock - 34.832) < 0.001) {
		printf("[");
		for (i = 0; i < frame_size; i++)
			printf("%lf%c", data1[i], ",]"[i == frame_size-1]);
		printf("\n%lf\n[", gain);
		for (i = 0; i < ar_order; i++)
			printf("%lf%c", a[i], ",]"[i == ar_order-1]);
		puts("");
		exit(0);
	}*/

	int l1 = 0;
	for (i = 0; i < 2*ar_order; i++) {
		if (fabs(std::abs(r[i])-1) > 1e-6) continue;
		double th = std::arg(r[i]), check = 0, val = 0;
		if (th < 0) continue;
		double freq = (th - 2*atan(warp_c*sin(th)/(1+warp_c*cos(th)))) * 22050 / (2 * pi);
		if (!(freq <= 2000 && freq >= 110)) continue;
		for (j = 0; j <= ar_order; j++)
			check += -j * j * b[j] * cos(j * th);
		if (check < 0) continue;
		for (j = 1; j <= num_sub; j++)
			if (freq / j >= 110)
				res[frame_cnt][l1++].freq = freq / j;
	}
	for (i = 0; i < l1; i++) {
		double freq = res[frame_cnt][i].freq, val, sum = 0;
		for (j = 1; j <= num_sub; j++) {
			if (freq * j > 2000) continue;
			double val = 0, th = freq * j / 22050 * 2 * pi;
			th += 2*atan(warp_c*sin(th)/(1-warp_c*cos(th)));
			for (k = 0; k <= ar_order; k++)
				val += b[k] * cos(k * th);
			val = 1.0 / val - 100;
			if (val > 0) sum += val * sub_weigh[j];
		}
		sum = gain * sqrt(sum) / frame_size;
		//if (sum > 0) printf("%lf,%lf ", freq, sum);
		res[frame_cnt][i].amp = sum;
	}
	j = 0;
	for (i = 0; i < l1; i++)
		if (res[frame_cnt][i].amp > 0) {
			res[frame_cnt][j].freq = res[frame_cnt][i].freq;
			res[frame_cnt][j].amp = res[frame_cnt][i].amp;
			j++;
		}
	l1 = j;
	//puts("");

	reslen[frame_cnt] = l1;
	int l2 = reslen[frame_cnt-1];
	double *f1 = dp[frame_cnt], *f2 = dp[frame_cnt-1];
	int *t1 = dpt[frame_cnt];
	Result *r1 = res[frame_cnt], *r2 = res[frame_cnt-1];
	double freq1[ar_order*num_sub+num_chan], amp1[ar_order*num_sub+num_chan];
	double freq2[ar_order*num_sub+num_chan], amp2[ar_order*num_sub+num_chan];
	for (i = 0; i < l1 + num_chan; i++)
		if (i < l1) freq1[i] = r1[i].freq, amp1[i] = loudness(r1[i].amp, freq1[i]);
		else freq1[i] = exp((i - l1) / 12.0 * log(2) + log(110)), amp1[i] = 0;
	for (j = 0; j < l2 + num_chan; j++)
		if (j < l2) freq2[j] = r2[j].freq, amp2[j] = loudness(r2[j].amp, freq2[j]);
		else freq2[j] = exp((j - l2) / 12.0 * log(2) + log(110)), amp2[j] = 0;
	for (i = 0; i < l1 + num_chan; i++) {
		f1[i] = -1e30;
		for (j = 0; j < l2 + num_chan; j++) {
			double eval = f2[j] + amp1[i] - pow(fabs(log(freq1[i] / freq2[j])), 0.8) * (amp1[i] * amp2[j] + 300) * 0.1;
			if (eval >= f1[i]) {
				f1[i] = eval;
				t1[i] = j;
			}
		}
	}

	frame_cnt++;
}

void prepare() {
	d[0][0] = 1;
	d[1][0] = -warp_c;
	double tmp = 1 - warp_c*warp_c;
	int i, j, k;
	for (i = 1; i < warp_order; i++)
		d[1][i] = tmp, tmp *= warp_c;
	for (i = 2; i <= ar_order; i++)
		for (j = 0; j < warp_order; j++)
			for (k = 0; k <= j; k++)
				d[i][j] += d[i-1][k] * d[1][j-k];
}

void finalize() {
	int p = frame_cnt-1, q = 0, i;
	double max = -1e30;
	for (i = 0; i <= reslen[p]; i++)
		if (dp[p][i] > max) {
			max = dp[p][i];
			q = i;
		}
	while (p > 0) {
		track[p] = q;
		q = dpt[p][q], p--;
	}
	int cnt = 0;
	for (i = 1; i < frame_cnt; i++)
		if (track[i] < reslen[i]) cnt++;
	double *ori = (double *)malloc(sizeof(double) * cnt);
	int l = 0;
	for (i = 1; i < frame_cnt; i++)
		if (track[i] < reslen[i])
			ori[l++] = log(res[i][track[i]].freq / 220) / log(2) * 12;
	printf("%d\n", cnt);
	for (i = 0; i < cnt; i++)
		printf("%.2lf%c", ori[i], " \n"[i == cnt-1]);
	/*
	for (i = 1; i < frame_cnt; i++)
		if (track[i] >= reslen[i])
			printf("%lf %.1lf 0 0\n", ts[i], exp((track[i] - reslen[i]) / 12.0 * log(2) + log(110)));
		else
			printf("%lf %.1lf %.0lf %.0lf\n", ts[i], res[i][track[i]].freq, res[i][track[i]].amp, loudness(res[i][track[i]].amp, res[i][track[i]].freq) * 8.7);
			*/
}

int main(int argc, char **argv) {
	
	av_register_all();
	avformat_network_init();

	format_ctx = avformat_alloc_context();
	if (avformat_open_input(&format_ctx, argv[1], av_find_input_format(argv[1]), NULL) < 0)
		return -1;
	if (avformat_find_stream_info(format_ctx, NULL) < 0)
		return -1;

	int audio_stream = -1, i;
	for (i = 0; i < format_ctx->nb_streams; i++)
		if (format_ctx->streams[i]->codec->codec_type == AVMEDIA_TYPE_AUDIO && audio_stream < 0)
			audio_stream = i;
	if (audio_stream == -1)
		return -1;

	audio_st = format_ctx->streams[audio_stream];
	codec_ctx = format_ctx->streams[audio_stream]->codec;
	if ((codec = avcodec_find_decoder(codec_ctx->codec_id)) == NULL)
		return -1;
	codec_ctx->codec_id = codec->id;
	if (avcodec_open2(codec_ctx, codec, NULL) < 0)
		return -1;
	src.fmt = AV_SAMPLE_FMT_S16;
	src.freq = 22050;
	src.channel_layout = codec_ctx->channel_layout;
	//src.channels = av_get_channel_layout_nb_channels(codec_ctx->channel_layout);
	src.channels = codec_ctx->channels;
	if (!src.channel_layout)
		src.channel_layout = av_get_default_channel_layout(src.channels);
	tgt = src;

	prepare();

	int x = 0;

	/*
	for (; x < frame_size; x++)
		if (x > 200) data[x] = sin(0.2 * (x - 200)) * 10000 + sin(0.4 * (x - 200)) * 5000;
		else data[x] = 0;
	analyze();
	return 0;
	*/

	while (1) {
		if (pkt.data)
			av_free_packet(&pkt);
		int ret = av_read_frame(format_ctx, &pkt);
		if (ret < 0)
			if (ret == AVERROR_EOF || url_feof(format_ctx->pb))
				break;
			else if (format_ctx->pb && format_ctx->pb->error)
				return -1;
			else {
				av_usleep(0.01 * 1000000);
				continue;
			}

		pkt1 = pkt;
		int got_frame = 0;
		while (pkt1.size > 0) {
			int len1 = avcodec_decode_audio4(codec_ctx, &frame1, &got_frame, &pkt1);
			if (len1 < 0) break;
			pkt1.data += len1;
			pkt1.size -= len1;

			if (!got_frame) continue;

			int data_size = av_samples_get_buffer_size(NULL, av_frame_get_channels(&frame1), frame1.nb_samples, (AVSampleFormat)frame1.format, 1);

			int64_t dec_channel_layout = frame1.channel_layout;
			if (!dec_channel_layout || av_frame_get_channels(&frame1) != av_get_channel_layout_nb_channels(frame1.channel_layout))
				dec_channel_layout = av_get_default_channel_layout(av_frame_get_channels(&frame1));

			if (frame1.format != src.fmt || dec_channel_layout != src.channel_layout || frame1.sample_rate != src.freq) {
				src.channel_layout = dec_channel_layout;
				src.channels = av_frame_get_channels(&frame1);
				src.freq = frame1.sample_rate;
				src.fmt = (AVSampleFormat)frame1.format;
				if (swr_ctx) swr_free(&swr_ctx);
				swr_ctx = swr_alloc_set_opts(NULL, tgt.channel_layout, tgt.fmt, tgt.freq,
						src.channel_layout, src.fmt, src.freq, 0, NULL);
				if (swr_ctx == NULL || swr_init(swr_ctx) < 0)
					return -1;
			}

			if (swr_ctx) {
				const uint8_t **in = (const uint8_t **)frame1.extended_data;
				uint8_t **out = &buf1;
				int out_count = (int64_t)frame1.nb_samples * tgt.freq / frame1.sample_rate + 256;
				int out_size = av_samples_get_buffer_size(NULL, tgt.channels, out_count, tgt.fmt, 0);
				
				av_fast_malloc(&buf1, &buf1_size, out_size);
				int len2 = swr_convert(swr_ctx, out, out_count, in, frame1.nb_samples);
				if (len2 < 0)
					break;
				data_size = len2 * tgt.channels * av_get_bytes_per_sample(tgt.fmt);
				memcpy(buf, buf1, data_size);
			} else
				memcpy(buf, frame1.data[0], data_size);

			buf_size = data_size;

			int i, j;
			int len = data_size / tgt.channels / av_get_bytes_per_sample(tgt.fmt);
			for (i = 0; i < len; i++) {
				data[x] = 0;
				for (j = 0; j < tgt.channels; j++)
					data[x] += ((short*)buf)[i * tgt.channels + j];
				x++;
				if (x >= frame_size) {
					analyze();
					for (j = 0; j < frame_size - hop_size; j++)
						data[j] = data[j + hop_size];
					x = frame_size - hop_size;
					audio_clock = av_q2d(audio_st->time_base) * pkt.pts + (double)(i - frame_size + hop_size) / tgt.freq;
				}
			}
		}
	}

	finalize();

	return 0;
}
