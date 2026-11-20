/* Per-call-site counts of halo exchanges (diagnostic for halo_timing, enabled by halo_timing_sites
   in namelist.config). Each record is keyed by the return addresses of the first frames above the
   exchange wrapper; the dump lists them as (object file, offset) for addr2line, with the exchange
   group, the number of calls, fields and seconds. */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <execinfo.h>
#include <stdio.h>
#include <string.h>

#define NFR 6
#define NSITE 4096
typedef struct { void *pc[NFR]; int cat, end; long calls, fields; double sec; } site_t;
static site_t sites[NSITE];
static int nsite = 0, overflow = 0;

void fesom_halo_site(int cat, int end, int nfields, double sec)
{
    void *bt[NFR + 2];
    int n = backtrace(bt, NFR + 2), i, k;
    void *pc[NFR] = {0};
    for (i = 2; i < n && i - 2 < NFR; i++) pc[i - 2] = bt[i];   /* skip this function and the wrapper */
    for (k = 0; k < nsite; k++)
        if (sites[k].cat == cat && sites[k].end == end && !memcmp(sites[k].pc, pc, sizeof pc)) break;
    if (k == nsite) {
        if (nsite == NSITE) { overflow++; return; }
        memcpy(sites[k].pc, pc, sizeof pc); sites[k].cat = cat; sites[k].end = end; nsite++;
    }
    sites[k].calls++; sites[k].fields += nfields; sites[k].sec += sec;
}

void fesom_halo_site_dump(int rank)
{
    char fn[64];
    int k, j;
    FILE *f;
    snprintf(fn, sizeof fn, "halo_sites_%05d.txt", rank);
    if (!(f = fopen(fn, "w"))) return;
    fprintf(f, "# cat end calls fields seconds | frames as object+offset (overflow %d)\n", overflow);
    for (k = 0; k < nsite; k++) {
        fprintf(f, "%d %d %ld %ld %.6f |", sites[k].cat, sites[k].end, sites[k].calls, sites[k].fields, sites[k].sec);
        for (j = 0; j < NFR && sites[k].pc[j]; j++) {
            Dl_info di;
            if (dladdr(sites[k].pc[j], &di) && di.dli_fname)
                fprintf(f, " %s+0x%lx", di.dli_fname,
                        (unsigned long)((char *)sites[k].pc[j] - (char *)di.dli_fbase) - 1);
            else
                fprintf(f, " ?+%p", sites[k].pc[j]);
        }
        fprintf(f, "\n");
    }
    fclose(f);
}
