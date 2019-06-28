#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#define N
#define L


float setbox(float l, int n, float *posicion);
float aleatorio(void);
float gaussiana(float mu, float sigma);
float velocidad(float *v, int n, float T);
int fuerzas(int n, float l, float *posicion, float *fuerza, float rc, float *fuerza_vieja);
int siguiente_paso(int n, float h, float *fuerza, float *velocidad, float rij_2, float dr);
int cc(int n, float l, float *posicion);
float velocidad_cuadrada(float *velocidad, int i);
float distancia_cuadrada(float *posicion, int i, int j, float *dr);
float energia_cinetica(int n, float *velocidad);
float potencial(int n, float l, float *posicion, float *dr);
float coeficiente_verlet(float l, int n, float *posicion);
float lennardjones(float r,float *potencial);
float fuerza_ln(float *fuerza, float *potencial, float rij_2, float *dr,int i, int j);
int velocity_verlet(float *posicion, float *velocidad, float *fuerza, float *fuerza_vieja,int n);


int main(int argc, char* argv[])
{

  int n;
  float l;
  float *v;
  float T;
  float h;
  float *posicion, *fuerza, rc, *fuerza_vieja;
 
  v=(float*)malloc(3*n*sizeof(float)); 
  fuerza=(float*)malloc(n*sizeof(float));
  fuerza_vieja=(float*)malloc(n*sizeof(float));
  posicion=(float*)malloc(n*sizeof(float));
 
h=0;
rc=0;
  

  if (argc==3)
   {
    sscanf(argv[1], "%d", &n);
    sscanf(argv[2], "%f", &l);
   }
  else printf("\nWrong number of parameters....\n");

setbox(l, n, posicion);
velocidad(v, n, T);

for(h=0; h<100; h=h+0.001)
{

	fuerzas(n, l, posicion, fuerza, rc, fuerza_vieja);
	velocity_verlet(posicion, v,fuerza,fuerza_vieja, n);
}
	//savelamps
 return 0;

}


//creo la cajita
float setbox(float l, int n, float *posicion) //buscar un nombre mas entendible jaja

{   

    
    int m;
    int i, j, k;
    float d;
   

 
    m=cbrt(n);
    d=l/m;
    i=0;
   
    for (i=0;i<m;i++)
        {
        for (j=0;j<m;j++)
          {
            for (k=0;k<m;k++)
                 {
                    *(posicion+3*i+0+(3*m)*j+(3*m*m)*k)=d*(i+0.5); 
                    *(posicion+3*i+1+(3*m)*j+(3*m*m)*k+1)=d*(j+0.5);
                    *(posicion+3*i+2+(3*m)*j+(3*m*m)*k+2)=d*(k+0.5);
                   
                    
                	 }
             }
         }
 
return l;

}


float aleatorio(void)

{

   float a;
   a=((float)(rand()))/((float)(RAND_MAX));

return a;

}


float gaussiana(float mu, float sigma)

{

int i,n;
float z;
n=10;
z=0;

for (i=0; i<n; i++)
   {
     z+=aleatorio();
    }

z=sqrt(12*n)*(z/n -0.5);

return z*sigma+mu;

}


float velocidad(float *v, int n, float T) //funcion para generar las velocidades de las particulas

{
    float sigma=sqrt(T);
    int i,j,k;
    float *vcm;
   // float var[3];


    v=(float*)malloc(3*n*sizeof(float));
    vcm=(float*)malloc(3*n*sizeof(float));
   
//    var=(float*)malloc(3*n*sizeof(float));
 
    for(i=0; i<3*n; i++)
       {
          *(v+i)=gaussiana(0,sigma);
        
      }

    for(i=0; i<n; i++)
       {
          for(j=0; j<3; j=j+3)
             {
               *(vcm+j)=*(v+3*i+j)/n;
             }
       }


for (i=0;i<3*n;i++)

   {  
     for(j=0;k<3;k++)
      {
         *(v+3*i+j)-=*(vcm+k);
       }

    }

return 0.0;

}





int fuerzas(int n, float l, float *posicion, float *fuerza, float rc, float *fuerza_vieja) 

{
    int i, j, k;
    float rij_1, rij_2;
    float *dr,*potencial;
    
	potencial=(float*)malloc(sizeof(float));
     dr=(float*)malloc(3*sizeof(float));
	
	for(i=0;i<3;i++){
		*(dr+i)=0.0;
		
	}
  

  rij_2=0;
	rij_1=0;  

    for(i=0; i<3*n; i++)
         {
           (*(fuerza_vieja+i))=(*(fuerza+i));
           *(fuerza+i)=0;
         }

       for(i=0; i<n-1; i++)
         {
           for(j=i+1; j<n; j++) 
            {
               rij_2=0;

               for(k=0; k<3; k++)
                  {
                    *(dr+k)=(*(posicion +i*3+k))-(*(posicion +j*3+k));
 
                    if (*(dr+k)>l/2)
                        {
                          *(dr+k)-=l;
    
                        }
               
                      else if(*(dr+k)<-l/2)
                           {
                               *(dr+k)+=l;
                            }

                 
                    }
                   rij_2=distancia_cuadrada(posicion, i, j, dr);
     		rij_1=sqrt(rij_2);
             

              if (rc>rij_1)
                    {
                     fuerza_ln(fuerza, potencial, rij_2, dr, i, j);
                     }
              }

               

        }
    

return 0;

}

int siguiente_paso(int n, float h, float *fuerza, float *v, float rij_2, float dr)

{  
    int i,j;
    
    for(i=0; i<3*n; i++)

       { for (j=0; j<3*n; j++)
  
        {

        fuerza_ln(fuerza, potencial, rij_2, dr, i, j);

         }
       }
return 0;

}

int cc(int n, float l, float *posicion) //para las condiciones de contorno

{
    int i;
    
    for(i=0; i<3*n; i++)

      
       {
        (*(posicion+i))=(*(posicion+i))-l*floor(*(posicion+i)/l);//mayor valor

       }

return 0;

}


//completo la energia
float velocidad_cuadrada(float *v, int i)
{
   return (v[3*i+0]*v[3*i+0])+(v[3*i+1]*v[3*i+1])+(v[3*i+2]*v[3*i+2]);
}

float distancia_cuadrada(float *posicion, int i, int j, float *dr)
{
    int k;
    float r_cuadrado=0.0;
	for(k=0;k<3;k++){
		*(dr+k)=*(posicion+3*i+k)-*(posicion+3*j+k);
		r_cuadrado+=(*(dr+k))*(*(dr+k));
	}
return r_cuadrado;

}



float energia_cinetica(int n, float *v)

{
    int i;
    float e_cinetica;
	e_cinetica=0.0;
    
    for(i=0; i<n; i++)
      {   
        e_cinetica+=velocidad_cuadrada(v,i)/2.0;

        }

return e_cinetica;

}


float coeficiente_verlet(float l, int n, float *posicion) 
{
    int i, m;
    float d, lambda_x, lambda_y, lambda_z, lambda, c;
    float PI;
    PI=3.14;
    d=l/n; 
    m=cbrt(n);
    c=(2*PI)/d; 
    lambda_x=0;
    lambda_y=0;
    lambda_z=0;
       
    for (i=0; i<=n-2;i++) 
        {
          lambda_x+=cos(c*(*(posicion+3*i-m)));
          lambda_y+=cos(c*(*(posicion+3*i+1-m)));
          lambda_z+=cos(c*(*(posicion+3*i+2-m)));
        }

    lambda=(lambda_x+lambda_y+lambda_z)/(3*n); 

return lambda;

}
/*
float Bolztmann(float m, float *velocidad, float T)
{
    int i;
    float h;
    h=0;

    for (i=0;i<3*m;i++) 
          {
                h+=H(*(velocidad+i),T);
         }

return h;
}

float H(float veloc, float T)
 {
    float b, b_a;
    float PI;
    PI=3.14;
    b=(4.0*PI)*pow((2.0*PI*T),(-3.0/2.0))*(*(velocidad))*(*(velocidad))*exp(-(veloc*veloc)/(2.0*T));
    b_a=b*log(b)*0.05;

return b_a;

}

*/



float fuerza_ln(float *fuerza, float *potencial, float rij_2, float *dr,int i, int j)
{

float b,e;
int k;
 float rsext=pow(1/rij_2,3);
 float rdoc=pow(rsext,2);
  
  e=1;

 *potencial=4*e*(rdoc-rsext);
for(k=0;k<3;k++)
    {
	b=(*(potencial)+8*e*rdoc)*(6*(*(dr+k))/rij_2);
	*(fuerza+3*i+k)+=b;
	*(fuerza+3*j+k)-=b;	
    }

return 0;
}


int velocity_verlet(float *posicion, float *v, float *fuerza, float *fuerza_vieja,int n)
{

//float v_auxiliar;
int i;
float h=0.01;




for(i=0; i<3*n; i++)
   {
    *(posicion+i)+=*(v+i)*h+0.5*h*h*(*(fuerza+i));
	//(v_auxiliar)=*(velocidad+i)*h+0.5*h*(*(fuerza+i));
    *(v+i)+=0.5*h*((*(fuerza+i))+(*(fuerza_vieja+i)));
   }
    
return 0;

}



