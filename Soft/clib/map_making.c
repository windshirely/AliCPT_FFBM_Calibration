#include<stdio.h>
#include<stdlib.h>

void rectang_mapmaking(double *az, double *el, double *s, double *w, double *m, double *mw, double xmin, double ymin, int nx, int ny, double xstep, double ystep, int len_t)
{
	int xind, yind;
	int i;
//	printf("len_tch_sum=%d\n", len_t);
	for(i=0;i<len_t;i++)
	{
//	    printf("i=%d\n", i);
//        printf("az=%f\n", *(az+i));
//        printf("el=%f\n", *(el+i));
//        printf("xmin=%f\n", xmin);
//        printf("ymin=%f\n", ymin);
		xind = (int)((*(az+i)-xmin)/xstep);
		yind = ny-(int)((*(el+i)-ymin)/ystep);
		if(xind<0 || yind<0){
			printf("xind=%d\n", xind);
		    printf("yind=%d\n", yind);
		}
        if(yind*nx+xind>nx*ny){
            printf("idx out of range!!!\n");
        }
//		printf("ind2=%d\n", yind*nx+xind);
		*(m+yind*nx+xind) += *(s+i);
		*(mw+yind*nx+xind) += 1.;
//		getchar();
//		*(m+xind*ny+yind) += *(s+i);
//		*(mw+xind*ny+yind) += 1.;
//		printf("hah\n");
	}
}

