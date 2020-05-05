#include<stdio.h>
#include<stdlib.h>

void beam_convolve(long long *bm_idx_array, long long *src_pix_idx, double *Ba, double *src, double *scale, double *da, int len_t, int len_bm, int len_src)
{
	int i, j, k;
	for(i=0;i<len_t;i++){
	    for(j=0;j<len_src;j++){
//	        printf("src%lld\n", *(src_pix_idx+j));
//	        getchar();
	        for(k=0;k<len_bm;k++){
	            // find the pix which are equal
//	            printf("%d,%d,%d\n", i,j,k);
//	            printf("bms%lld\n", *(bm_idx_array+k*len_t+i));
//	            getchar();
//	            if(*(src_pix_idx+j)==*(bm_idx_array+i*len_bm+k)){
	            if(*(src_pix_idx+j)==*(bm_idx_array+k*len_t+i)){
	                *(da+i) += *(src+j)*(*(Ba+k))*(*(scale+i));
//	                printf("%d,%d,%d\n", i,j,k);
//	                printf("%f\n", *(da+i));
//	                printf("bme%lld\n", *(bm_idx_array+k*len_t+i));
//	                getchar();
//	                if((*(da+i))>1000) {
//	                    printf("%d,%d,%d\n", i,j,k);
//	                    printf("bme%d\n", *(bm_idx_array+k*len_t+i));
//	                    printf("%f\n", *(da+i));
//	                }
	            }
	        }
	    }
	}
}