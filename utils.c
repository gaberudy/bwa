/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Contact: Heng Li <lh3@sanger.ac.uk> */
#define FSYNC_ON_FLUSH

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <errno.h>
#ifdef FSYNC_ON_FLUSH
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#endif
#include <sys/resource.h>
#include <sys/time.h>
#include "utils.h"

#include "ksort.h"
#define pair64_lt(a, b) ((a).x < (b).x || ((a).x == (b).x && (a).y < (b).y))
KSORT_INIT(128, pair64_t, pair64_lt)
KSORT_INIT(64, uint64_t, ks_lt_generic)

#include "kseq.h"

/// Expand out this macro so it's easy to debug and inspect
// KSEQ_INIT2(, gzFile, err_gzread)
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
typedef struct __kstream_t
{
	unsigned char *buf;
	int begin, end, is_eof;
	gzFile f;
} kstream_t;

static inline kstream_t *ks_init(gzFile f)
{
	kstream_t *ks = (kstream_t *)calloc(1, sizeof(kstream_t));
	ks->f = f;
	ks->buf = (unsigned char *)malloc(16384);
	return ks;
}
static inline void ks_destroy(kstream_t *ks)
{
	if (ks)
	{
		free(ks->buf);
		free(ks);
	}
}
static inline int ks_getc(kstream_t *ks)
{
	if (ks->is_eof && ks->begin >= ks->end)
		return -1;
	if (ks->begin >= ks->end)
	{
		ks->begin = 0;
		ks->end = err_gzread(ks->f, ks->buf, 16384);
		if (ks->end == 0)
		{
			ks->is_eof = 1;
			return -1;
		}
	}
	return (int)ks->buf[ks->begin++];
}

/// delimiter is one of:
/// 0: isspace(c)
/// 1: isspace(c) && c != ' '
/// 2: c == '\n'
/// > 2: the delimiter as a char itself
static int ks_getuntil2(kstream_t *ks, int delimiter, kstring_t *str, int *dret, int append)
{
	int gotany = 0;
	if (dret)
		*dret = 0;
	str->l = append ? str->l : 0;
	for (;;)
	{
		int i;
		if (ks->begin >= ks->end)
		{
			if (!ks->is_eof)
			{
				ks->begin = 0;
				ks->end = err_gzread(ks->f, ks->buf, 16384);
				if (ks->end == 0)
				{
					ks->is_eof = 1;
					break;
				}
			}
			else
				break;
		}
		if (delimiter == 2 /* line separator: "\n" (Unix) or "\r\n" (Windows)*/)
		{
			for (i = ks->begin; i < ks->end; ++i)
				if (ks->buf[i] == '\n')
					break;
		}
		else if (delimiter > 2)
		{
			for (i = ks->begin; i < ks->end; ++i)
				if (ks->buf[i] == delimiter)
					break;
		}
		else if (delimiter == 0 /* isspace(): \t, \n, \v, \f, \r*/)
		{
			for (i = ks->begin; i < ks->end; ++i)
				if (isspace(ks->buf[i]))
					break;
		}
		else if (delimiter == 1 /* isspace() && !' '*/)
		{
			for (i = ks->begin; i < ks->end; ++i)
				if (isspace(ks->buf[i]) && ks->buf[i] != ' ')
					break;
		}
		else
			i = 0; /* never come to here! */
		if (str->m - str->l < (size_t)(i - ks->begin + 1))
		{
			str->m = str->l + (i - ks->begin) + 1;
			/// Round up to next power of 2
			/// See https://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
			(--(str->m), (str->m) |= (str->m) >> 1, (str->m) |= (str->m) >> 2, (str->m) |= (str->m) >> 4, (str->m) |= (str->m) >> 8, (str->m) |= (str->m) >> 16, ++(str->m));
			str->s = (char *)realloc(str->s, str->m);
		}
		gotany = 1;
		memcpy(str->s + str->l, ks->buf + ks->begin, i - ks->begin);
		str->l = str->l + (i - ks->begin);
		ks->begin = i + 1;
		if (i < ks->end)
		{
			if (dret)
				*dret = ks->buf[i];
			break;
		}
	}
	if (!gotany && ((ks)->is_eof && (ks)->begin >= (ks)->end))
		return -1;
	if (str->s == 0)
	{
		str->m = 1;
		str->s = (char *)calloc(1, 1);
	}
	else if (delimiter == 2 /* line separator: "\n" (Unix) or "\r\n" (Windows)*/ && str->l > 1 && str->s[str->l - 1] == '\r')
		--str->l;
	str->s[str->l] = '\0';
	return str->l;
}

static inline int ks_getuntil(kstream_t *ks, int delimiter, kstring_t *str, int *dret)
{
	return ks_getuntil2(ks, delimiter, str, dret, 0);
}

typedef struct
{
	kstring_t name, comment, seq, qual;
	int last_char;
	kstream_t *f;
} kseq_t;

kseq_t *kseq_init(gzFile fd)
{
	kseq_t *s = (kseq_t *)calloc(1, sizeof(kseq_t));
	s->f = ks_init(fd);
	return s;
}
void kseq_destroy(kseq_t *ks)
{
	if (!ks)
		return;
	free(ks->name.s);
	free(ks->comment.s);
	free(ks->seq.s);
	free(ks->qual.s);
	ks_destroy(ks->f);
	free(ks);
}

/// Reads a FASTA entry of a FASTQ seq + qual pair, re-using the "seq"
/// data-structure that may have been used in the previous entry
///
/// Return value:
///   >=0  length of the sequence (normal) -1   end-of-file -2   truncated
///   quality string
int kseq_read(kseq_t *seq)
{
	int c;
	kstream_t *ks = seq->f;

	/// then jump to the next header line (read until we have the first char of the line)
	if (seq->last_char == 0)
	{
		/// Scan for '>' or '@', but sets -1 for EOF
		while ((c = ks_getc(ks)) != -1 && c != '>' && c != '@')
			;
		if (c == -1)
			return -1; /// EOF
		seq->last_char = c;
	}
	/// else: the first header char has been read in the previous call

	/// Reset our strings to length = 0, but note we are re-using their memory
	/// buffors already allocated and expanded to largest previous seq
	seq->comment.l = seq->seq.l = seq->qual.l = 0;

	/// Get until whitepace
	if (ks_getuntil(ks, 0, &seq->name, &c) < 0)
		return -1; /// This is the EOF normal exit point when we hit last sequence

	/// If not EOL, then get comment until newline
	/// 2 means line separator: "\n" (Unix) or "\r\n" (Windows)
	if (c != '\n')
		ks_getuntil(ks, 2, &seq->comment, 0);

	/// If our seq string is non-initialize, pre-alloc 256b as estimate of what
	/// is needed, but ks_getuntil2 will re-alloc the kstring seq->seq as-needed
	if (seq->seq.s == 0)
	{
		seq->seq.m = 256;
		seq->seq.s = (char *)malloc(seq->seq.m);
	}
	// While the beginning of line char is not a >, +, or @ and not EOL
	while ((c = ks_getc(ks)) != -1 && c != '>' && c != '+' && c != '@')
	{
		/// skip blank lines
		if (c == '\n')
			continue;
		/// Save our current scan char
		/// this is safe: we always have enough space for 1 char
		seq->seq.s[seq->seq.l++] = c;

		/// read the rest of the line (2: EOL delim)
		ks_getuntil2(ks, 2, &seq->seq, 0, 1);
	}
	/// the first header char has been read
	if (c == '>' || c == '@')
		seq->last_char = c;

	/// We need to ensure that subsequent seq->seq.s[seq->seq.l] access is not
	/// out of bounds (generally won't take this branch)
	if (seq->seq.l + 1 >= seq->seq.m)
	{
		/// Clone of code in ks_getuntil2 that decides a new capacity for the
		/// string and then realloc to that new capacity.
		seq->seq.m = seq->seq.l + 2;
		(--(seq->seq.m), (seq->seq.m) |= (seq->seq.m) >> 1, (seq->seq.m) |= (seq->seq.m) >> 2, (seq->seq.m) |= (seq->seq.m) >> 4, (seq->seq.m) |= (seq->seq.m) >> 8, (seq->seq.m) |= (seq->seq.m) >> 16, ++(seq->seq.m));
		seq->seq.s = (char *)realloc(seq->seq.s, seq->seq.m);
	}

	/// null terminated string
	seq->seq.s[seq->seq.l] = 0;

	/// FASTA has no quality line, so we return here.
	if (c != '+')
		return seq->seq.l;

	/// allocate memory for qual in case insufficient (qual is 1:1 with seq)
	if (seq->qual.m < seq->seq.m)
	{
		seq->qual.m = seq->seq.m;
		seq->qual.s = (char *)realloc(seq->qual.s, seq->qual.m);
	}

	/// skip the rest of '+' line (we don't care about its contents)
	while ((c = ks_getc(ks)) != -1 && c != '\n')
		;

	/// error: no quality string  (EOL hit prematurely)
	if (c == -1)
		return -2;

	/// Read lines into qual bufer (break if qual >= seq length)
	while (ks_getuntil2(ks, 2, &seq->qual, 0, 1) >= 0 && seq->qual.l < seq->seq.l)
		;

	/// 0 is our sentinal that we have not come to the next header line and
	/// next call into this function should scan can until we hit it.
	seq->last_char = 0;
	/// Strict expectation of length of seq and qual matching
	if (seq->seq.l != seq->qual.l)
		return -2;
	/// On success, return seq length
	return seq->seq.l;
}
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// KSEQ_INIT2(, gzFile, err_gzread)

/********************
 * System utilities *
 ********************/

FILE *err_xopen_core(const char *func, const char *fn, const char *mode)
{
	FILE *fp = 0;
	if (strcmp(fn, "-") == 0)
		return (strstr(mode, "r")) ? stdin : stdout;
	if ((fp = fopen(fn, mode)) == 0)
	{
		err_fatal(func, "fail to open file '%s' : %s", fn, strerror(errno));
	}
	return fp;
}

FILE *err_xreopen_core(const char *func, const char *fn, const char *mode, FILE *fp)
{
	if (freopen(fn, mode, fp) == 0)
	{
		err_fatal(func, "fail to open file '%s' : %s", fn, strerror(errno));
	}
	return fp;
}

gzFile err_xzopen_core(const char *func, const char *fn, const char *mode)
{
	gzFile fp;
	if (strcmp(fn, "-") == 0)
	{
		fp = gzdopen(fileno((strstr(mode, "r")) ? stdin : stdout), mode);
		/* According to zlib.h, this is the only reason gzdopen can fail */
		if (!fp)
			err_fatal(func, "Out of memory");
		return fp;
	}
	if ((fp = gzopen(fn, mode)) == 0)
	{
		err_fatal(func, "fail to open file '%s' : %s", fn, errno ? strerror(errno) : "Out of memory");
	}
	return fp;
}

void err_fatal(const char *header, const char *fmt, ...)
{
	va_list args;
	va_start(args, fmt);
	fprintf(stderr, "[%s] ", header);
	vfprintf(stderr, fmt, args);
	fprintf(stderr, "\n");
	va_end(args);
	exit(EXIT_FAILURE);
}

void err_fatal_core(const char *header, const char *fmt, ...)
{
	va_list args;
	va_start(args, fmt);
	fprintf(stderr, "[%s] ", header);
	vfprintf(stderr, fmt, args);
	fprintf(stderr, " Abort!\n");
	va_end(args);
	abort();
}

void _err_fatal_simple(const char *func, const char *msg)
{
	fprintf(stderr, "[%s] %s\n", func, msg);
	exit(EXIT_FAILURE);
}

void _err_fatal_simple_core(const char *func, const char *msg)
{
	fprintf(stderr, "[%s] %s Abort!\n", func, msg);
	abort();
}

size_t err_fwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream)
{
	size_t ret = fwrite(ptr, size, nmemb, stream);
	if (ret != nmemb)
		_err_fatal_simple("fwrite", strerror(errno));
	return ret;
}

size_t err_fread_noeof(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
	size_t ret = fread(ptr, size, nmemb, stream);
	if (ret != nmemb)
	{
		_err_fatal_simple("fread", ferror(stream) ? strerror(errno) : "Unexpected end of file");
	}
	return ret;
}

int err_gzread(gzFile file, void *ptr, unsigned int len)
{
	int ret = gzread(file, ptr, len);

	if (ret < 0)
	{
		int errnum = 0;
		const char *msg = gzerror(file, &errnum);
		_err_fatal_simple("gzread", Z_ERRNO == errnum ? strerror(errno) : msg);
	}

	return ret;
}

int err_fseek(FILE *stream, long offset, int whence)
{
	int ret = fseek(stream, offset, whence);
	if (0 != ret)
	{
		_err_fatal_simple("fseek", strerror(errno));
	}
	return ret;
}

long err_ftell(FILE *stream)
{
	long ret = ftell(stream);
	if (-1 == ret)
	{
		_err_fatal_simple("ftell", strerror(errno));
	}
	return ret;
}

int err_printf(const char *format, ...)
{
	va_list arg;
	int done;
	va_start(arg, format);
	done = vfprintf(stdout, format, arg);
	int saveErrno = errno;
	va_end(arg);
	if (done < 0)
		_err_fatal_simple("vfprintf(stdout)", strerror(saveErrno));
	return done;
}

int err_fprintf(FILE *stream, const char *format, ...)
{
	va_list arg;
	int done;
	va_start(arg, format);
	done = vfprintf(stream, format, arg);
	int saveErrno = errno;
	va_end(arg);
	if (done < 0)
		_err_fatal_simple("vfprintf", strerror(saveErrno));
	return done;
}

int err_fputc(int c, FILE *stream)
{
	int ret = putc(c, stream);
	if (EOF == ret)
	{
		_err_fatal_simple("fputc", strerror(errno));
	}

	return ret;
}

int err_fputs(const char *s, FILE *stream)
{
	int ret = fputs(s, stream);
	if (EOF == ret)
	{
		_err_fatal_simple("fputs", strerror(errno));
	}

	return ret;
}

int err_puts(const char *s)
{
	int ret = puts(s);
	if (EOF == ret)
	{
		_err_fatal_simple("puts", strerror(errno));
	}

	return ret;
}

int err_fflush(FILE *stream)
{
	int ret = fflush(stream);
	if (ret != 0)
		_err_fatal_simple("fflush", strerror(errno));

#ifdef FSYNC_ON_FLUSH
	/* Calling fflush() ensures that all the data has made it to the
	   kernel buffers, but this may not be sufficient for remote filesystems
	   (e.g. NFS, lustre) as an error may still occur while the kernel
	   is copying the buffered data to the file server.  To be sure of
	   catching these errors, we need to call fsync() on the file
	   descriptor, but only if it is a regular file.  */
	{
		struct stat sbuf;
		if (0 != fstat(fileno(stream), &sbuf))
			_err_fatal_simple("fstat", strerror(errno));

		if (S_ISREG(sbuf.st_mode))
		{
			if (0 != fsync(fileno(stream)))
				_err_fatal_simple("fsync", strerror(errno));
		}
	}
#endif
	return ret;
}

int err_fclose(FILE *stream)
{
	int ret = fclose(stream);
	if (ret != 0)
		_err_fatal_simple("fclose", strerror(errno));
	return ret;
}

int err_gzclose(gzFile file)
{
	int ret = gzclose(file);
	if (Z_OK != ret)
	{
		_err_fatal_simple("gzclose", Z_ERRNO == ret ? strerror(errno) : zError(ret));
	}

	return ret;
}

/*********
 * Timer *
 *********/

double cputime()
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

double realtime()
{
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}
