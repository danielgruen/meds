#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "meds.h"

int main(int argc, char **argv)
{
    if (argc < 5) {
        printf("usage: meds-get-collec meds_file outfile iobj_fits iobs_last\n"
               "  extract a mosaic of all cutouts for a collection of objects of index\n"
               "  iobj_fits-iobj_last.\n"
               "  iobj are zero offset, equal to NUMBER-1 in the \n"
               "  catalog.\n");
        exit(1);
    }
    
    const char *meds_file=argv[1];
    const char *outfile=argv[2];
    int index_first=atoi(argv[3]);
    int index_last=atoi(argv[4]);
    int n=index_last-index_first+1;
    if(n<1) {
              printf("usage: meds-get-collec meds_file outfile iobj_fits iobs_last\n"
               "  extract a mosaic of all cutouts for a collection of objects of index\n"
               "  iobj_fits-iobj_last.\n"
               "  iobj are zero offset, equal to NUMBER-1 in the \n"
               "  catalog.\n");
              exit(1);
    }
    
    struct meds *meds=meds_open(meds_file,READONLY);
    if (!meds) {
        fprintf(stderr,"error reading meds, exiting\n");
        exit(1);
    }

    struct meds_cutout *image[n];
    struct meds_icutout *seg[n];
    struct meds_icutout *crmask[n];

    int in=0;
    for(int index=index_first; index<=index_last; index++)
    {
      image[in] = meds_get_mosaic(meds, index);
      seg[in]   = meds_get_seg_mosaic(meds, index);
      crmask[in] = meds_get_crmask_mosaic(meds, index);
      in++;
    }

    int status=0;
    int clobber=1;
    meds_cutouts_write_fits(image, n, outfile, clobber, &status);	

    char soutfile[strlen(outfile)+1];
    strcpy(soutfile+1,outfile);
    soutfile[0]='s';
    meds_icutouts_write_fits(seg, n, soutfile, clobber, &status);
        
    strcpy(soutfile+1,outfile);
    soutfile[0]='m';
    meds_icutouts_write_fits(crmask, n, soutfile, clobber, &status);
    
    for(in=0; in<n; in++)
    {
      image[in]=meds_cutout_free(image[in]);
      seg[in]=meds_icutout_free(seg[in]);
      crmask[in]=meds_icutout_free(crmask[in]);
    }

    meds=meds_free(meds);
}

