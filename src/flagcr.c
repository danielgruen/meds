#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "meds.h"

int max(int a, int b)
{
 if(a>=b) return a;
 return b; 
}

int min(int a, int b)
{
 if(a>=b) return b;
 return a; 
}

double median(double **data, double **weight, int n)
// return median of n data points, discarding ones with weight <= 0
{
    double sorted[n];
    const int N = n;
    int i=0, ii=0;

    while(i<N) {
	if((*(weight[i]))>0) {
	        sorted[ii] = (*(data[i]));
		++ii;
	}
	else {
		--n;
	}
	++i;
    }

    for (int i = n - 1; i > 0; --i) {
        for (int j = 0; j < i; ++j) {
            if (sorted[j] > sorted[j+1]) {
                double d = sorted[j];
                sorted[j] = sorted[j+1];
                sorted[j+1] = d;
            }
        }
    }

    // Middle or average of middle values in the sorted array.
    double dm = 0.0;
    if ((n % 2) == 0) {
        dm = (sorted[n/2] + sorted[(n/2) - 1])/2.0;
    } else {
        dm = sorted[n/2];
    }
    return dm;
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        printf("usage: flagcr meds_file\n");
        exit(1);
    }

    const char *meds_file=argv[1];

    printf("opening meds file: %s\n", meds_file);
    struct meds *meds=meds_open(meds_file,READWRITE);

    if (!meds) {
        fprintf(stderr,"error reading meds, exiting\n");
        exit(1);
    }

    //struct meds_cutout *mosaic=NULL;
    //struct meds_icutout *smosaic=NULL;

    long nobj=meds_get_size(meds);
    printf("running through %ld objects\n", nobj);

    long nrow[1]; long ncol[1];
    double rowcen=0,colcen=0;

    for (long iobj=0; iobj<nobj; iobj++) {
	printf("\nprocessing object %ld\n",iobj);
        //mosaic=meds_get_mosaic(meds, iobj);
        //smosaic=meds_get_seg_mosaic(meds, iobj);

	long ncutout=meds_get_ncutout(meds, iobj);

	double * vcutouts[ncutout-1];
	int    * vsegcutouts[ncutout-1];
	double * vwtcutouts[ncutout-1];
	char   update[ncutout-1];
	int    vicutout[ncutout-1];
	int ngood = 0;
	
	for(long icutout=1; icutout<ncutout; ++icutout) {
		//printf("getting cutout %ld %ld\n",iobj,icutout);
		meds_get_cutout_cen(meds, iobj, icutout, &rowcen, &colcen);
	        //printf("center pixel [%lf,%lf]\n", rowcen, colcen);
		double * vcutout = meds_get_cutoutp(meds, iobj, icutout, nrow, ncol);
		double dr = rowcen-(*nrow)/2.;
		double dc = colcen-(*ncol)/2.;
		if(dr>1 || dr<0 || dc>1 || dc<0) {
			fprintf(stderr,"skipping off-centered object\n");
			free(vcutout);
			continue;
		}
		//printf("sounds good, saving %ld\n",icutout); 
		vcutouts[ngood]    = vcutout;
		//printf("got cutout\n"); 
		vsegcutouts[ngood] = meds_get_seg_cutoutp(meds, iobj, icutout, nrow, ncol);
		//printf("got seg\n"); 
		vwtcutouts[ngood] = meds_get_weight_cutoutp(meds, iobj, icutout, nrow, ncol);
		//printf("got weight\n"); 
		update[ngood] = 0;
		vicutout[ngood] = icutout;
		ngood++;

		// TODO: make sure all are of the same size
		// TODO: make sure WCS is not causing us trouble
	}

	int npix = (*nrow)*(*ncol);
	int mask[] = { -100000, -100000, -100000 };
	for(int ii=0; ii<npix; ++ii)
	{
		//if(ii%100==0) printf("\npixel %d: ",ii);
		double mu; mu = median(vcutouts,vwtcutouts,ngood);
		double amu; amu = fabs(mu);
		//if(ii%100==0) printf("median %lf\n(val,out): ",mu);
		for (int icutout=0; icutout<ngood; icutout++)
		{
		  double wt=*(vwtcutouts[icutout]);
		  if(wt>0) // it's a good pixel
		  {
		    double delta = fmax(fabs(*(vcutouts[icutout])-mu)-0.3*amu,0)*sqrt(wt);
		    //if(ii%100==0) printf("(%lf,%lf) ",*(vcutouts[icutout]),delta);
		    if(delta>5) { // 5sigma outlier, mask small box
		      update[icutout]=1;
		      int col=ii%(*nrow);
		      int row=(ii-col)/(*nrow);
		      for (int irow=max(row-1,0); irow<=min(row+1,(*nrow)-1); irow++)
		      {
			memcpy(vsegcutouts[icutout]-ii+irow*(*nrow)+max(col-1,0),mask,sizeof(int)*min(3,(*ncol)-col+1));
		      }
		    }
		  }
		  
		  vcutouts[icutout]++;
		  vwtcutouts[icutout]++;
		  vsegcutouts[icutout]++;
		}
	}
	

	
	printf("\nneed to update? ");
	for(int icutout=0; icutout<ngood; icutout++) 
	{
	 printf("%d ",(int)(update[icutout]));
	 if(update[icutout])
	   meds_update_seg_cutout(meds, iobj, vicutout[icutout], vsegcutouts[icutout]-npix);
	 
	 vcutouts[icutout]   -= npix;
	 vwtcutouts[icutout] -= npix;
	 vsegcutouts[icutout]-= npix;
	 
	 free(vcutouts[icutout]);
	 free(vsegcutouts[icutout]);
	 free(vwtcutouts[icutout]);
	}
	
    }

    meds=meds_free(meds);
}
