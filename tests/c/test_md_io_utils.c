/* Unit tests for md_libs/md_io_utils.c (built and run by tests/test_c_units.py). */
#include "md_io_utils.h"

static int failures = 0;

static void check_complex(const char *text, int expect_ok, double re, double im)
{
	COMPLEX z = 12345.0;
	int ret = parse_complex(text, &z);
	if (expect_ok && (ret != 0 || fabs(creal(z) - re) > 1e-15 || fabs(cimag(z) - im) > 1e-15)) {
		printf("FAIL parse_complex(\"%s\"): ret=%d value=%g%+gi, expected %g%+gi\n", text, ret, creal(z), cimag(z), re, im);
		failures++;
	}
	if (!expect_ok && ret == 0) {
		printf("FAIL parse_complex(\"%s\") accepted (value %g%+gi)\n", text, creal(z), cimag(z));
		failures++;
	}
}

static void check_file_readers(void)
{
	const char *name = "unit_test_input.txt";
	FILE *fp = fopen(name, "w");
	double tab[4];
	int n;
	COMPLEX z;
	char s[SIZE_STR_BUFFER];
	fprintf(fp, "# comment line with N = 99\nN = 12 ; NS = 3\nname = abc  # trailing comment\n"
	            "nu = 1.5 - i0.25\narr = 1 2.5\n -3e-1 4\n");
	fclose(fp);
	fp = fopen(name, "r");
	if (lire_int(fp, "N", &n) || n != 12)   { printf("FAIL lire_int N -> %d\n", n); failures++; }
	if (lire_int(fp, "NS", &n) || n != 3)   { printf("FAIL lire_int NS -> %d\n", n); failures++; }
	if (lire_string(fp, "name", s) || strcmp(s, "abc")) { printf("FAIL lire_string -> %s\n", s); failures++; }
	if (lire_complex(fp, "nu", &z) || creal(z) != 1.5 || cimag(z) != -0.25) { printf("FAIL lire_complex\n"); failures++; }
	if (lire_int(fp, "missing", &n) == 0)   { printf("FAIL lire_int accepted a missing key\n"); failures++; }
	fclose(fp);
	if (lire_tab(name, "arr", tab, 4) || tab[0] != 1 || tab[1] != 2.5 || tab[2] != -0.3 || tab[3] != 4) {
		printf("FAIL lire_tab\n"); failures++;
	}
	remove(name);
}

int main(void)
{
	check_complex("1.3 + i7.1", 1, 1.3, 7.1);
	check_complex("1.3+i7.1", 1, 1.3, 7.1);
	check_complex("  1.3 - i7.1", 1, 1.3, -7.1);
	check_complex("1.3-i7.1", 1, 1.3, -7.1);
	check_complex("1.3 + i-7.1", 1, 1.3, -7.1);
	check_complex("1.3 - i-7.1", 1, 1.3, 7.1);
	check_complex("1.3 + i 7.1", 1, 1.3, 7.1);
	check_complex("1.3 + 7.1i", 1, 1.3, 7.1);
	check_complex("1.3-7.1*i", 1, 1.3, -7.1);
	check_complex("1.3 - 7.1 * I", 1, 1.3, -7.1);
	check_complex("1.5", 1, 1.5, 0.0);
	check_complex("-2e-3 + i1e2", 1, -2e-3, 100.0);
	check_complex("1.5 ; next = 2", 1, 1.5, 0.0);
	check_complex("1.3 + 7.1", 0, 0, 0);
	check_complex("abc", 0, 0, 0);
	check_complex("", 0, 0, 0);
	check_complex("1.3 + i", 0, 0, 0);
	check_file_readers();
	if (failures == 0) printf("all md_io_utils unit tests passed\n");
	return failures != 0;
}
