#include<stdio.h>
#include<math.h>
#define grid 129
#define Re 400.0
#define emax .00001
int main()
{
FILE *u,*v,*s,*w;

int i,j,k=0;
float S[200][200],temp_S[200][200],W[200][200],temp_W[200][200],U[1][200],V[1][200];
float X[200],Y[200],dx=(1/128),dy=(1/128);
float N=0,D=0,e=0,er=1;
printf("Obseravtions For Reynold's no.:400\n");

for(i=0;i<grid;i++)
{
for(j=0;j<grid;j++)
{
S[i][j]=0;
temp_S[i][j]=0;
if(i==0||i==128)
{
W[i][j]=0;
temp_W[i][j]=0;
}
else
{
if(j==128)
{
W[i][j]=-2/dy;
temp_W[i][j]=-2/dy;
}
else
{
W[i][j]=0;
temp_W[i][j]=0;
}
}
}
}
//do
//
while(er>emax)
{
k=k+1;
printf("k=%d\n",k);
for(i=1;i<(grid-1);i++)
{
for(j=1;j<(grid-1);j++)
{
temp_S[i][j]=(temp_S[i+1][j]+temp_S[i-1][j]+temp_S[i][j+1]+temp_S[i][j-1]+(dx*dx*W[i][j]))*0.25;
}
}

for(j=1;j<(grid-1);j++)
{
W[0][j]=(-2*temp_S[1][j])/(dx*dx);
temp_W[0][j]=W[0][j];
W[128][j]=(-2*temp_S[127][j])/(dx*dx);
temp_W[128][j]=W[128][j];
}

for(i=1;i<(grid-1);i++)
{
W[i][0]=(-2*temp_S[i][1])/(dy*dy);
temp_W[i][0]=W[i][0];
W[i][128]=((-2*temp_S[i][127])-(2*dy))/(dy*dy);
temp_W[i][128]=W[i][128];
}
for(i=1;i<(grid-1);i++)
{
for(j=1;j<(grid-1);j++)
{
temp_W[i][j]=(temp_W[i+1][j]+temp_W[i-1][j]+W[i][j+1]+temp_W[i][j-1]-(0.25*Re*(W[i+1][j]-temp_W[i-1][j])*(S[i][j+1]-S[i][j-1]))+(0.25*Re*(W[i][j+1]-temp_W[i][j-1])*(S[i+1][j]-S[i-1][j])))*0.25;
}
}

for(i=0;i<grid;i++)
{
for(j=0;j<grid;j++)
{
N=N+fabs(temp_W[i][j]-W[i][j]);
D=D+fabs(temp_W[i][j]);
}
}
e=e+(N/D);
//printf("%f",e);
for(i=0;i<grid;i++)
{
for(j=0;j<grid;j++)
{
S[i][j]=temp_S[i][j];
W[i][j]=temp_W[i][j];
}
}
er=e;
e=0;
N=0;
D=0;
}
//while(er>emax);
printf("\n\n");
printf("Error:%f\n",er);
printf("Total no. of Iteration till convergence:%d",k);

X[0]=0;
Y[0]=0;
for(j=0;j<grid;j++)
{
Y[j+1]=Y[j]+dy;
}
for(i=0;i<grid;i++)
{
X[i+1]=X[i]+dx;
}

s=fopen("Stream Function.dat","w");
fprintf(s,"X\tY\tStream function\n",grid,grid);
for(int p=0;p<grid;p++)
{
for(int q=0;q<grid;q++)
{
fprintf(s,"%f\t%f\t%f\n",X[p],Y[q],S[p][q]);
}
}
fclose(s);


w=fopen("Vorticity.dat","w");
fprintf(w,"X\tY\tVorticity\n",grid,grid);
for(int p=0;p<grid;p++)
{
for(int q=0;q<grid;q++)
{
fprintf(w,"%f\t%f\t%f\n",X[p],Y[q],W[p][q]);
}
}
fclose(w);


for(j=0;j<grid;j++)
{
U[0][j]=(S[64][j+1]-S[64][j-1])/(2*dy);
V[0][j]=(S[65][j]-S[63][j])/(2*dx);
}

u=fopen("U-Velocity.dat","w");
fprintf(u,"X\tY\tU-velocity",grid,grid);
for(int q=0;q<grid;q++)
{
fprintf(u,"l/2\t%f\t%f\n",Y[q],U[0][q]);
}
fclose(u);

v=fopen("V-Velocity.dat","w");
fprintf(v,"X\tY\tU-velocity",grid,grid);
for(int q=0;q<grid;q++)
{
fprintf(v,"l/2\t%f\t%f\n",Y[q],V[0][q]);
}
fclose(v);

}
