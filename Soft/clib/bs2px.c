#include<stdio.h>
#include<math.h>
#define PI acos(-1.)

void bs2px_altaz(double *bs_az, double *bs_el, double *bs_dk, \
                 double *pxa_r, double *pxa_theta, double *pxb_r, double *pxb_theta, \
                 double *pxa_az, double *pxa_el, double *pxb_az, double *pxb_el, \
                 int lent, int lenfp)
{
	#pragma omp parallel for num_threads(10)
	for(int i=0;i<lenfp;i++)
	{
		for(int j=0;j<lent;j++)
		{
			//printf("%d\n", j);
			double alpha_a = *(pxa_theta+i)-(PI/2-*(bs_dk+j));
			double alpha_b = *(pxb_theta+i)-(PI/2-*(bs_dk+j));
			// calculate the elevation for each pixel
			double sin_pxa_el = sin(*(bs_el+j))*cos(*(pxa_r+i))+cos(*(bs_el+j))*sin(*(pxa_r+i))*cos(alpha_a);
			double sin_pxb_el = sin(*(bs_el+j))*cos(*(pxb_r+i))+cos(*(bs_el+j))*sin(*(pxb_r+i))*cos(alpha_b);
//			if(i==0){
//			    printf("%f\n", alpha_a);
////			    printf("%f\n", sin_pxa_el);
//			}
			double el_pxa = asin(sin_pxa_el);
			double el_pxb = asin(sin_pxb_el);
			*(pxa_el+i*lent+j) = el_pxa;
			*(pxb_el+i*lent+j) = el_pxb;
			// calculate the azimuth for each pixel
			double sin_beta_a = sin(alpha_a)*sin(*(pxa_r+i))/cos(el_pxa);
			double cos_beta_a = (cos(*(pxa_r+i))-sin(*(bs_el+j))*sin(el_pxa))/(cos(*(bs_el+j))*cos(el_pxa));
			*(pxa_az+i*lent+j) = *(bs_az+j)-atan2(sin_beta_a, cos_beta_a);
			double sin_beta_b = sin(alpha_b)*sin(*(pxb_r+i))/cos(el_pxb);
			double cos_beta_b = (cos(*(pxb_r+i))-sin(*(bs_el+j))*sin(el_pxb))/(cos(*(bs_el+j))*cos(el_pxb));
			*(pxb_az+i*lent+j) = *(bs_az+j)-atan2(sin_beta_b, cos_beta_b);
		}
	}
}


void bs2px_with_dk(double *bs_az, double *bs_el, double *bs_dk, \
                 double *pxa_r, double *pxa_theta, double *pxb_r, double *pxb_theta, \
                 double *pxa_az, double *pxa_el, double *pxb_az, double *pxb_el, double *px_dk,\
                 int lent, int lenfp)
{
	#pragma omp parallel for num_threads(10)
	for(int i=0;i<lenfp;i++)
	{
		for(int j=0;j<lent;j++)
		{
			//printf("%d\n", j);
			double alpha_a = *(pxa_theta+i)-(PI/2-*(bs_dk+j));
			double alpha_b = *(pxb_theta+i)-(PI/2-*(bs_dk+j));
			// calculate the elevation for each pixel
			double sin_pxa_el = sin(*(bs_el+j))*cos(*(pxa_r+i))+cos(*(bs_el+j))*sin(*(pxa_r+i))*cos(alpha_a);
			double sin_pxb_el = sin(*(bs_el+j))*cos(*(pxb_r+i))+cos(*(bs_el+j))*sin(*(pxb_r+i))*cos(alpha_b);
//			if(i==0){
//			    printf("%f\n", alpha_a);
////			    printf("%f\n", sin_pxa_el);
//			}
			double el_pxa = asin(sin_pxa_el);
			double el_pxb = asin(sin_pxb_el);
			*(pxa_el+i*lent+j) = el_pxa;
			*(pxb_el+i*lent+j) = el_pxb;
			// calculate the azimuth for each pixel
			double sin_beta_a = sin(alpha_a)*sin(*(pxa_r+i))/cos(el_pxa);
			double cos_beta_a = (cos(*(pxa_r+i))-sin(*(bs_el+j))*sin(el_pxa))/(cos(*(bs_el+j))*cos(el_pxa));
			*(pxa_az+i*lent+j) = *(bs_az+j)-atan2(sin_beta_a, cos_beta_a);
			double sin_beta_b = sin(alpha_b)*sin(*(pxb_r+i))/cos(el_pxb);
			double cos_beta_b = (cos(*(pxb_r+i))-sin(*(bs_el+j))*sin(el_pxb))/(cos(*(bs_el+j))*cos(el_pxb));
			*(pxb_az+i*lent+j) = *(bs_az+j)-atan2(sin_beta_b, cos_beta_b);
			// calculate the polrization angle
			double sin_BPZ_a = sin(alpha_a)/cos(el_pxa)*cos(*(bs_el+j));
			double cos_BPZ_a = (sin(el_pxa)-sin(*(bs_el+j))*cos(*(pxa_r+i)))/(cos(*(bs_el+j))*sin(*(pxa_r+i)));
			*(px_dk+i*lent+j) = atan2(sin_BPZ_a, cos_BPZ_a);
		}
	}
}

void bs2px_bm(double *bs_az, double *bs_el, double *bs_dk, \
                 double *r, double *theta, \
                 double *px_az, double *px_el, \
                 int lent, int lenfp)
{
	#pragma omp parallel for num_threads(10)
	for(int i=0;i<lenfp;i++)
	{
		for(int j=0;j<lent;j++)
		{
			//printf("%d\n", j);
//			double alpha = *(theta+i)-(PI/2-*(bs_dk+j));
			double alpha = *(bs_dk+j)+*(theta+i)-PI/2;
			// calculate the elevation for each pixel
			double sin_px_el = sin(*(bs_el+j))*cos(*(r+i))+cos(*(bs_el+j))*sin(*(r+i))*cos(alpha);
//			if(i==0){
//			    printf("%f\n", alpha_a);
////			    printf("%f\n", sin_pxa_el);
//			}
			double el_px = asin(sin_px_el);
			*(px_el+i*lent+j) = el_px;
			// calculate the azimuth for each pixel
			double sin_beta = sin(alpha)*sin(*(r+i))/cos(el_px);
			double cos_beta = (cos(*(r+i))-sin(*(bs_el+j))*sin(el_px))/(cos(*(bs_el+j))*cos(el_px));
			*(px_az+i*lent+j) = *(bs_az+j)-atan2(sin_beta, cos_beta);
		}
	}
}


void bs2px(double *RA, double *DEC, double *DK, double *R, double *THETA, double *CHI, double *RA_CH, double *DEC_CH, double *PSI_CH, int lent, int lenfp)
{
	//int i, j;
	//double alpha;
	//double sin_dec_ch, decch;
	//double sin_beta, cos_beta, beta;
	//double sin_BPZ, cos_BPZ, BPZ;
	#pragma omp parallel for num_threads(8)
	for(int i=0;i<lenfp;i++)
	{
		for(int j=0;j<lent;j++)
		{
			//printf("%d\n", j);
			double alpha = *(THETA+i)-(PI/2-*(DK+j));
			// calculate the dec for channel
			double sin_dec_ch = sin(*(DEC+j))*cos(*(R+i))+cos(*(DEC+j))*sin(*(R+i))*cos(alpha);
			double decch = asin(sin_dec_ch);
			//double decch = asin(sin(*(DEC+j))*cos(*(R+i))+cos(*(DEC+j))*sin(*(R+i))*cos(alpha));
			*(DEC_CH+i*lent+j) = decch;
			// calculate the ra
			double sin_beta = sin(alpha)*sin(*(R+i))/cos(decch);
			double cos_beta = (cos(*(R+i))-sin(*(DEC+j))*sin(decch))/(cos(*(DEC+j))*cos(decch));
			//beta = atan2(sin_beta, cos_beta);
			*(RA_CH+i*lent+j) = *(RA+j)+atan2(sin_beta, cos_beta);
			//*(RA_CH+i*lent+j) = *(RA+j)+atan2(sin(alpha)*sin(*(R+i))/cos(decch), cos(*(R+i))-sin(*(DEC+j))*sin(decch))/(cos(*(DEC+j))*cos(decch));
			// calculate the polarization angle
			double sin_BPZ = sin(alpha)/cos(decch)*cos(*(DEC+j));
			double cos_BPZ = (sin(decch)-sin(*(DEC+j))*cos(*(R+i)))/(cos(*(DEC+j))*sin(*(R+i)));
			//double BPZ = atan2(sin_BPZ, cos_BPZ);
			*(PSI_CH+i*lent+j) = *(CHI+i)-atan2(sin_BPZ, cos_BPZ)+PI;
			//*(PSI_CH+i*lent+j) = *(CHI+i)-atan2(sin(alpha)/cos(decch)*cos(*(DEC+j)), sin(decch)-sin(*(DEC+j))*cos(*(R+i)))/(cos(*(DEC+j))*sin(*(R+i)))+PI;
		}
	}
}
//
//void src2bmcoord(double *pxaz, double *pxel, double *pxdk, double *bm_az, double *bm_el, double src_az, double src_el, int lent)
//{
//    #pragma omp parallel for num_threads(8)
//    for(int i=0;i<lent;i++)
//    {
//        // in triangle SZP
//        double szp = saz-paz;
//        double zp = PI/2-pel;
//        double sz = PI/2-sel;
//        double cos_sp = cos(sz)*cos(zp)+sin(sz)*sin(zp)*cos(szp);
//        double sp = acos(cos_sp);
//        double sin_spz = sin(szp)/sin(sp)*sin(sz);
//        double cos_spz = (cos(sz)-cos(zp)*cos(sp))/(sin(zp)*sin(sp));
//        double spz = atan2(sin_spz, cos_spz);
//        double zpf = theta-dk;
//        double spf = zpf-spz;
//        double sin_sf = sin(spf)*sin(sp);
//        double sf = asin(sin_sf);
//        double cos_fp = cos(sp)/cos(sf);
//        double sin_fp = (cos(sf)-cos(sp)*cos(fp))/(sin(sp)*cos(spf));
//        double fp = atan2(sin_fp, cos_fp);
//    }
//}
