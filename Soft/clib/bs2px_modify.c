#include<stdio.h>
#include<math.h>
#define PI acos(-1.)

void bs2px_with_dk(float *bs_az, float *bs_el, float *bs_dk, \
                 float *pxa_r, float *pxa_theta, float *pxb_r, float *pxb_theta, \
                 float *pxa_az, float *pxa_el, float *pxb_az, float *pxb_el, float *px_dk,\
                 int lent, int lenfp)
{
	#pragma omp parallel for num_threads(10)
	for(int i=0;i<lenfp;i++)
	{
		for(int j=0;j<lent;j++)
		{
			//printf("%d\n", j);
			float alpha_a = *(pxa_theta+i)-(PI/2-*(bs_dk+j));
			float alpha_b = *(pxb_theta+i)-(PI/2-*(bs_dk+j));
			// calculate the elevation for each pixel
			float sin_pxa_el = sin(*(bs_el+j))*cos(*(pxa_r+i))+cos(*(bs_el+j))*sin(*(pxa_r+i))*cos(alpha_a);
			float sin_pxb_el = sin(*(bs_el+j))*cos(*(pxb_r+i))+cos(*(bs_el+j))*sin(*(pxb_r+i))*cos(alpha_b);
//			if(i==0){
//			    printf("%f\n", alpha_a);
////			    printf("%f\n", sin_pxa_el);
//			}
			float el_pxa = asin(sin_pxa_el);
			float el_pxb = asin(sin_pxb_el);
			*(pxa_el+i*lent+j) = el_pxa;
			*(pxb_el+i*lent+j) = el_pxb;
			// calculate the azimuth for each pixel
			float sin_beta_a = sin(alpha_a)*sin(*(pxa_r+i))/cos(el_pxa);
			float cos_beta_a = (cos(*(pxa_r+i))-sin(*(bs_el+j))*sin(el_pxa))/(cos(*(bs_el+j))*cos(el_pxa));
			*(pxa_az+i*lent+j) = *(bs_az+j)+atan2(sin_beta_a, cos_beta_a);
			float sin_beta_b = sin(alpha_b)*sin(*(pxb_r+i))/cos(el_pxb);
			float cos_beta_b = (cos(*(pxb_r+i))-sin(*(bs_el+j))*sin(el_pxb))/(cos(*(bs_el+j))*cos(el_pxb));
			*(pxb_az+i*lent+j) = *(bs_az+j)+atan2(sin_beta_b, cos_beta_b);
			// calculate the polrization angle
			float sin_BPZ_a = sin(alpha_a)/cos(el_pxa)*cos(*(bs_el+j));
			float cos_BPZ_a = (sin(el_pxa)-sin(*(bs_el+j))*cos(*(pxa_r+i)))/(cos(*(bs_el+j))*sin(*(pxa_r+i)));
			*(px_dk+i*lent+j) = -atan2(sin_BPZ_a, cos_BPZ_a)+PI;
		}
	}
}