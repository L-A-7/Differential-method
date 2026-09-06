#include<stdio.h>
#include<math.h>

int main()
{
	int i1,i2,i3,j1,j2,j3,N;
	double ax,ay,az,bx,by,bz,cos_theta;
	N=10;

	for(i1=0;i1<=N;i1++){
		for(i2=0;i2<=N;i2++){
			for(i3=0;i3<=N;i3++){
				for(j1=0;j1<=N;j1++){
					for(j2=0;j2<=N;j2++){
						for(j3=0;j3<=N;j3++){
	
							ax=sqrt(i1);
							ay=i2;
							az=i3;
							bx=j1;
							by=sqrt(j2);
							bz=j3;
	
							cos_theta= (ax*bx+ay*by+az*bz)/(sqrt(ax*ax+ay*ay+az*az)*sqrt(bx*bx+by*by+bz*bz));
							if (cos_theta ==0.5){
								printf("%d %d %d %d %d %d\n",i1,i2,i3,j1,j2,j3);
							}
						}
					}
				}
			}
		}
	}






	return 0;
}
