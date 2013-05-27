#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "meds.h"

int main(int argc, char **argv)
{
    if (argc < 4) {
        printf("usage: meds-get-cutout meds_file outfile iobj [icutout]\n"
               "  extract a mosaic of all cutouts for the object of index\n"
               "  iobj.  If icutout is sent, only extract the cutout \n"
               "  at that index.\n"
               "  iobj are zero offset, equal to NUMBER-1 in the \n"
               "  catalog.  icutout are zero offset with the coadd cutout\n"
               "  at icutout==0\n");
        exit(1);
    }

    const char *meds_file=argv[1];
    const char *outfile=argv[2];
    int index=atoi(argv[3]);

    struct meds *meds=meds_open(meds_file);
    if (!meds) {
        fprintf(stderr,"error reading meds, exiting\n");
        exit(1);
    }

    struct meds_cutout *image = NULL;
    struct meds_icutout *seg    = NULL;
    if (argc > 4) {
        int icutout=atoi(argv[4]);
        image = meds_get_cutout(meds, index, icutout);
        seg   = meds_get_seg_cutout(meds, index, icutout);
    } else {
        image = meds_get_mosaic(meds, index);
        seg   = meds_get_seg_mosaic(meds, index);
    }

    if (image) {
        int status=0;
        int clobber=1;
        meds_cutout_write_fits(image, outfile, clobber, &status);
    }
    
    if (seg) {
        int status=0;
        int clobber=1;
	char soutfile[strlen(outfile)+1];
	strcpy(soutfile+1,outfile);
	soutfile[0]='s';
        meds_icutout_write_fits(seg, soutfile, clobber, &status);
    }

    image=meds_cutout_free(image);
    seg=meds_icutout_free(seg);
    meds=meds_free(meds);
}
