#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <omp.h>
#include <openacc.h>

#define END " \033[1;32m  =)\033[0m "
#define WAR " \033[1;33m [!]\033[0m "
#define ERR " \033[1;31m [ERROR]\033[0m "
#define OK " \033[1;35m [OK]\033[0m "

#define PI 3.1415926535
#define _g_ 9.81
#define epsilon 1.0E-14
#define _gamma_ 1.4

#define TOL4 1.0E-4
#define TOL8 1.0E-8
#define TOL14 1.0E-14
#define MIN(x,y) (x < y ? x : y) 
#define MAX(x,y) (x > y ? x : y)
#define ABS(x) (x < 0 ? -x : x)

#define WENO 1
#define OPT_WEIGHTS 0 

//Solvers
#define LINEAR_TRANSPORT 0 
#define EULER 1
#define SW 0
#define HLLE 1
#define HLLC 0
#define ROE 0

//Debug code 
#define DEBUG_MESH 0 //0: no debug; 1: screen info; 
//Debug cases
#define DEBUG_FRONTOGENESIS 0  //0: no debug; 1: debug;

#define NTHREADS 27

////////////////////////////////////////////////////
//////////////  S T R U C T U R E S  ///////////////
////////////////////////////////////////////////////

typedef struct t_node_ t_node;
typedef struct t_cell_ t_cell; 
typedef struct t_wall_ t_wall;
typedef struct t_mesh_ t_mesh;
typedef struct t_sim_ t_sim;


struct t_node_{
	int id;
	double x,y; 

};


struct t_cell_{
	int id;
	int l,m;
	double *U; 
	double *U_aux;
	double pres,u_int; 
	double dx,dy;
	double xc,yc;
	int n1,n2,n3,n4; //ID's of the nodes of the cell
	int w1_id,w2_id,w3_id,w4_id; //ID's of the walls of the cell
	t_wall *w1, *w2, *w3, *w4;
	int st_sizeX, st_sizeY;		//stencil size
	int stX[9], stY[9];		//id's of the cells in the X and Y stencil.
					
};


struct t_wall_{
	int id;
	int stencil;
	double *UL, *UR; //array of extrapolated values at the left and right side of the wall
	double *fR_star,*fL_star;
	int cellR_id, cellL_id; //id of the right and left cell
	t_cell *cellR, *cellL; //pointers to the right and left cells of the wall
	double nx, ny;

};



struct t_mesh_{
	int xcells, ycells; //values given by the user
      double dx, dy;
	int ncells;
	int nwalls;
	int nnodes;
      int bc[4]; 
	int flux_bc_flag,cell_bc_flag; 
      int periodicX,periodicY;
	t_cell *cell;
	t_wall *wall;
	t_node *node;
	double lambda_max;
	double mass;
	t_sim *sim;
};

struct t_sim_{
	double dt,t,CFL; 
	double tf, tVolc; 
	int rk_steps;
	int order; 
	int nvar; //number of variables

};


////////////////////////////////////////////////////
//////  F U N C T I O N   D E F I N I T I O N //////
////////////////////////////////////////////////////


int create_mesh(t_mesh *mesh,t_sim *sim);
int update_initial(t_mesh *mesh);
int update_cell_boundaries(t_mesh *mesh);
int update_dt(t_mesh *mesh,t_sim *sim);

int compute_wallX(t_cell *cell);
int compute_wallY(t_cell *cell);
int half_update(t_cell *cell); 

double weno3L(double *phi);
double weno3R(double *phi);

void compute_solver(t_wall *wall);
void compute_transport(t_wall *wall);

void update_cell(t_mesh *mesh, t_sim *sim);

void update_cellK1(t_mesh *mesh, t_sim *sim);
void update_cellK2(t_mesh *mesh, t_sim *sim);
void update_cellK3(t_mesh *mesh, t_sim *sim);

void compute_euler_HLLE(t_wall *wall,double *lambda_max);

void mass_calculation(t_mesh *mesh, t_sim *sim);

////////////////////////////////////////////////////
//////  P R E - P R O C.   F U N C T I O N S ///////
////////////////////////////////////////////////////

int create_mesh(t_mesh *mesh, t_sim *sim){
	int i,l,m,k,aux,p;
	int xcells,ycells; 
	t_cell *cell;
	t_wall *wall;
	t_node *node;

	int semiSt; 
	
	mesh->sim=sim;

	//Cells
	xcells=mesh->xcells;
	ycells=mesh->ycells;
	mesh->ncells=xcells*ycells;
	mesh->cell=(t_cell*)malloc(mesh->ncells*sizeof(t_cell));
	cell=mesh->cell;
      
	for(m=0;m<ycells;m++){
		for(l=0;l<xcells;l++){
			
			k=xcells*m+l;
			cell[k].id=999;
			cell[k].l=l;
			cell[k].m=m;
                  
			cell[k].w1_id=999;
			cell[k].w2_id=999;
			cell[k].w3_id=999;
			cell[k].w4_id=999;

			cell[k].dx=mesh->dx;
			cell[k].dy=mesh->dy;
                  
                  cell[k].xc=999;
			cell[k].yc=999;
                  
			cell[k].n1=999;
			cell[k].n2=999;
			cell[k].n3=999;
			cell[k].n4=999;
		}
	}
      
	//Walls
	mesh->nwalls=2*mesh->ncells+xcells+mesh->ycells;
	mesh->wall=(t_wall*)malloc(mesh->nwalls*sizeof(t_wall));
      
#ifdef _OPENACC
    acc_attach((void**)&mesh->sim);
#pragma acc update device(mesh->ncells,mesh->nwalls)
#pragma acc enter data create(mesh->cell[:mesh->ncells],mesh->wall[:mesh->nwalls])
#endif

	wall=mesh->wall;
	for(k=0;k<mesh->nwalls;k++){
		wall[k].id=k;
		wall[k].fR_star=(double*)malloc(sim->nvar*sizeof(double));
		wall[k].fL_star=(double*)malloc(sim->nvar*sizeof(double));
		wall[k].UR=(double*)malloc(sim->nvar*sizeof(double));
		wall[k].UL=(double*)malloc(sim->nvar*sizeof(double));
#pragma acc enter data create(wall[k].fR_star[:sim->nvar], wall[k].fL_star[:sim->nvar])
	}

	//Walls amd nodes of the cells
	for(m=0;m<ycells-1;m++){
		for(l=0;l<xcells-1;l++){	
			k=xcells*m+l;
			cell[k].id=k;
			cell[k].l=l;
			cell[k].m=m;	

			cell[k].w1_id=2*(k)+m;
			cell[k].w2_id=cell[k].w1_id+3;
			cell[k].w3_id=cell[k].w1_id+xcells*2+1;
			cell[k].w4_id=cell[k].w1_id+1;

			cell[k].w1=&(mesh->wall[cell[k].w1_id]);
			cell[k].w2=&(mesh->wall[cell[k].w2_id]);
			cell[k].w3=&(mesh->wall[cell[k].w3_id]);
			cell[k].w4=&(mesh->wall[cell[k].w4_id]);
                  
#ifdef _OPENACC
#pragma acc update device(cell[k:1])
                  acc_attach((void**)&cell[k].w1);
                  acc_attach((void**)&cell[k].w2);
                  acc_attach((void**)&cell[k].w3);
                  acc_attach((void**)&cell[k].w4);
#endif

			cell[k].dx=mesh->dx;
			cell[k].dy=mesh->dy;
                  
                  cell[k].xc=(l+0.5)*cell[k].dx;
			cell[k].yc=(m+0.5)*cell[k].dy;
                  
                 	cell[k].n1=k+m;
                 	cell[k].n2=cell[k].n1+1;
                 	cell[k].n3=cell[k].n2+xcells+1;
                	cell[k].n4=cell[k].n2+xcells;
		}
	}
		for(m=0;m<mesh->ycells-1;m++){
		
		l=xcells-1;
            k=xcells*m+l;
            
            cell[k].id=k;
            cell[k].l=l;
            cell[k].m=m;

            cell[k].w1_id=2*(k)+m;
            cell[k].w2_id=cell[k].w1_id+3-1;
            cell[k].w3_id=cell[k].w1_id+xcells*2+1;
            cell[k].w4_id=cell[k].w1_id+1;
	
	      cell[k].w1=&(mesh->wall[cell[k].w1_id]);
		cell[k].w2=&(mesh->wall[cell[k].w2_id]);
		cell[k].w3=&(mesh->wall[cell[k].w3_id]);
		cell[k].w4=&(mesh->wall[cell[k].w4_id]);
            
#ifdef _OPENACC
#pragma acc update device(cell[k:1])
                  acc_attach((void**)&cell[k].w1);
                  acc_attach((void**)&cell[k].w2);
                  acc_attach((void**)&cell[k].w3);
                  acc_attach((void**)&cell[k].w4);
#endif

            cell[k].dx=mesh->dx;
            cell[k].dy=mesh->dy;
            
            cell[k].xc=(l+0.5)*cell[k].dx;
            cell[k].yc=(m+0.5)*cell[k].dy;
            
            cell[k].n1=k+m;
            cell[k].n2=cell[k].n1+1;
            cell[k].n3=cell[k].n2+xcells+1;
            cell[k].n4=cell[k].n2+xcells;

	}
      
	aux=2*(xcells*(ycells-1)+xcells-1)+(ycells-1)+3-1;

	for(l=0;l<mesh->xcells-1;l++){
		
		m=ycells-1;
		k=xcells*m+l;
            
		cell[k].id=k;
		cell[k].l=l;
		cell[k].m=m;

		cell[k].w1_id=2*(k)+m;
		cell[k].w2_id=cell[k].w1_id+3;
		cell[k].w3_id=aux+l+1;
		cell[k].w4_id=cell[k].w1_id+1;
	
	      cell[k].w1=&(mesh->wall[cell[k].w1_id]);
		cell[k].w2=&(mesh->wall[cell[k].w2_id]);
		cell[k].w3=&(mesh->wall[cell[k].w3_id]);
		cell[k].w4=&(mesh->wall[cell[k].w4_id]);
            
#ifdef _OPENACC
#pragma acc update device(cell[k:1])
            acc_attach((void**)&cell[k].w1);
            acc_attach((void**)&cell[k].w2);
            acc_attach((void**)&cell[k].w3);
            acc_attach((void**)&cell[k].w4);
#endif

		cell[k].dx=mesh->dx;
		cell[k].dy=mesh->dy;

		cell[k].xc=(l+0.5)*cell[k].dx;
		cell[k].yc=(m+0.5)*cell[k].dy;

		cell[k].n1=k+m;
		cell[k].n2=cell[k].n1+1;
		cell[k].n3=cell[k].n2+xcells+1;
		cell[k].n4=cell[k].n2+xcells;
	}

	m=ycells-1;
	l=xcells-1;
	k=xcells*m+l;

	cell[k].id=k;
	cell[k].l=l;
	cell[k].m=m;

	cell[k].w1_id=2*(k)+m;
	cell[k].w2_id=cell[k].w1_id+3-1;
	cell[k].w3_id=cell[k-1].w3_id+1;
	cell[k].w4_id=cell[k].w1_id+1;
	
      cell[k].w1=&(mesh->wall[cell[k].w1_id]);
	cell[k].w2=&(mesh->wall[cell[k].w2_id]);
	cell[k].w3=&(mesh->wall[cell[k].w3_id]);
	cell[k].w4=&(mesh->wall[cell[k].w4_id]);
#ifdef _OPENACC
#pragma acc update device(cell[k:1])
            acc_attach((void**)&cell[k].w1);
            acc_attach((void**)&cell[k].w2);
            acc_attach((void**)&cell[k].w3);
            acc_attach((void**)&cell[k].w4);
#endif 

	cell[k].dx=mesh->dx;
	cell[k].dy=mesh->dy;

	cell[k].xc=(l+0.5)*cell[k].dx;
	cell[k].yc=(m+0.5)*cell[k].dy;

	cell[k].n1=k+m;
	cell[k].n2=cell[k].n1+1;
	cell[k].n3=cell[k].n2+xcells+1;
	cell[k].n4=cell[k].n2+xcells;
	
  	for(k=0;k<mesh->ncells;k++){
		mesh->wall[cell[k].w1_id].nx=0.0;
		mesh->wall[cell[k].w4_id].nx=1.0;
		mesh->wall[cell[k].w1_id].ny=1.0;
		mesh->wall[cell[k].w4_id].ny=0.0;
		mesh->wall[cell[k].w3_id].nx=0.0;
		mesh->wall[cell[k].w2_id].nx=1.0;
		mesh->wall[cell[k].w3_id].ny=1.0;
		mesh->wall[cell[k].w2_id].ny=0.0;
	}
		
	//Nodes
	mesh->nnodes=(mesh->xcells+1)*(mesh->ycells+1);
	mesh->node=(t_node*)malloc(mesh->nnodes*sizeof(t_node));
	
	node=mesh->node;
	//Loop for nodes
	for(m=0;m<ycells+1;m++){
		for(l=0;l<xcells+1;l++){
                	k=(xcells+1)*m+l;
                 	node[k].id=k;
                 	node[k].x=l*mesh->dx;
                 	node[k].y=m*mesh->dy;
		}
	}
      
	if(mesh->bc[0]==99 && mesh->bc[1]==99 && mesh->bc[2]==99 && mesh->bc[3]==99){ 
		mesh->cell_bc_flag=1; 
 	}else{
		mesh->cell_bc_flag=0;
	}
	if(mesh->bc[0]!=1 || mesh->bc[1]!=1 || mesh->bc[2]!=1 || mesh->bc[3]!=1){
		mesh->flux_bc_flag=1; 
 	}else{
		mesh->flux_bc_flag=0;
	}

	if(mesh->bc[1]==1 && mesh->bc[3]==1){
		mesh->periodicX=1;
	}else{
		mesh->periodicX=0;
	}
	if(mesh->bc[0]==1 && mesh->bc[2]==1){
		mesh->periodicY=1;
	}else{
		mesh->periodicY=0;

	}

	//Set cell stencils
	//Initially all the cells has a stencil of size order
	for(k=0;k<mesh->ncells;k++){
		cell[k].st_sizeX=sim->order;
		cell[k].st_sizeY=sim->order;
	}

	//But there are special cases at boundary cells	
	semiSt=(sim->order-1)*0.5;
	for(m=0;m<ycells;m++){
		for(l=0;l<xcells;l++){
			k=xcells*m+l;
			if(mesh->periodicX==0){
				//x setencils
				if(l<semiSt){
					cell[k].st_sizeX=MIN(cell[k].st_sizeX,2*l+1);
				}
				if(xcells-(l+1)<semiSt){
					cell[k].st_sizeX=MIN(cell[k].st_sizeX,2*(xcells-(l+1))+1);
				}
			}
			if(mesh->periodicY==0){
				//y stencils
				if(m<semiSt){
					cell[k].st_sizeY=MIN(cell[k].st_sizeY,2*m+1);
				}
				if(ycells-(m+1)<semiSt){
					cell[k].st_sizeY=MIN(cell[k].st_sizeY,2*(ycells-(m+1))+1);
				}
			}
		}
	}
	
	//Initialization of cell stencils
	for(k=0;k<mesh->ncells;k++){
		for(p=0;p<9;p++){
			cell[k].stX[p]=-1;
			cell[k].stY[p]=-1;
		}
	}
	
	for(m=0;m<ycells;m++){
		for(l=0;l<xcells;l++){
			
			k=xcells*m+l;
			//x setencils
			for(p=0;p<cell[k].st_sizeX;p++){
				if(mesh->periodicX==0){
					cell[k].stX[p]=l-((cell[k].st_sizeX-1)/2)+p;
				}else{
					cell[k].stX[p]=l-((cell[k].st_sizeX-1)/2)+p;
					if(cell[k].stX[p]<0){
						cell[k].stX[p]+=xcells*(1);
					}
					if(cell[k].stX[p]>xcells-1){
						cell[k].stX[p]+=xcells*(-1);
					}
				}
                        cell[k].stX[p]+=xcells*m;
			}
                 
			
			//y stencils 
			for(p=0;p<cell[k].st_sizeY;p++){
				if(mesh->periodicY==0){
					cell[k].stY[p]=m-((cell[k].st_sizeY-1)/2)+p;
				}else{
					cell[k].stY[p]=m-((cell[k].st_sizeY-1)/2)+p;
					if(cell[k].stY[p]<0){
						cell[k].stY[p]+=ycells*(1);
					}
					if(cell[k].stY[p]>ycells-1){
						cell[k].stY[p]+=ycells*(-1);
					}
				}
                        cell[k].stY[p]=xcells*cell[k].stY[p]+l;
			}
		}
	}
	
	//Assigment of wall's neighbour cells 
	for(m=0;m<ycells;m++){
		for(l=0;l<xcells;l++){
			
			k=xcells*m+l;
			
			cell[k].w1->cellR_id=cell[k].id;
			cell[k].w4->cellR_id=cell[k].id;
			cell[k].w2->cellL_id=cell[k].id;
			cell[k].w3->cellL_id=cell[k].id;

			cell[k].w1->cellR=&(cell[k]);
			cell[k].w4->cellR=&(cell[k]);
			cell[k].w2->cellL=&(cell[k]);
			cell[k].w3->cellL=&(cell[k]);

			if(m==0){
				cell[k].w1->cellL_id=cell[xcells*(ycells-1)+l].id;
				cell[k].w1->cellL=&(cell[xcells*(ycells-1)+l]);
			}
			
			if(m==ycells-1){
				cell[k].w3->cellR_id=cell[l].id;
				cell[k].w3->cellR=&(cell[l]);
			}

			if(l==0){		
				cell[k].w4->cellL_id=cell[(m+1)*xcells-1].id;		
				cell[k].w4->cellL=&(cell[(m+1)*xcells-1]);		
			}
			if(l==xcells-1){	
				cell[k].w2->cellR_id=cell[m*xcells].id;		
				cell[k].w2->cellR=&(cell[m*xcells]);		
			}
                           
#ifdef _OPENACC
#pragma acc update device(cell[k].w1->cellR_id)
#pragma acc update device(cell[k].w4->cellR_id)
#pragma acc update device(cell[k].w2->cellL_id)
#pragma acc update device(cell[k].w3->cellL_id)
            acc_attach((void**)&cell[k].w1->cellR);
            acc_attach((void**)&cell[k].w4->cellR);
            acc_attach((void**)&cell[k].w2->cellL);
            acc_attach((void**)&cell[k].w3->cellL);
#endif

		}
	}

#if EULER
	for(k=0;k<mesh->ncells;k++){		
		mesh->cell[k].U=(double*)malloc(sim->nvar*sizeof(double));
		mesh->cell[k].U_aux=(double*)malloc(sim->nvar*sizeof(double));
#pragma acc enter data create(mesh->cell[k].U[:sim->nvar])
#pragma acc enter data create(mesh->cell[k].U_aux[:sim->nvar])
	}
	for(k=0;k<mesh->nwalls;k++){
		mesh->wall[k].UR=(double*)malloc(sim->nvar*sizeof(double));
		mesh->wall[k].UL=(double*)malloc(sim->nvar*sizeof(double));
	}
#endif
	
	return 1;
}

int update_initial(t_mesh *mesh){
	int l,m,k,n,RP;
      double d,aux1,xaux,yaux,r,umax,ut,w,wp,r_d;
      double p,u,v,rho,phi;
	t_cell *cell;

	cell=mesh->cell;

#if EULER
      for(k=0;k<mesh->ncells;k++){
            u=0.0;
            v=0.0;
            r=sqrt((cell[k].yc-0.5)*(cell[k].yc-0.5)+(cell[k].xc-0.5)*(cell[k].xc-0.5));
            if (r<0.2499) {
                  p=3.0;
                  rho=1.0  ;  
                  phi=1.0;
            }else{
                  p=1.0;
                  rho=5.0;
                  phi=2.0;
            }
		
		cell[k].U[0]=rho;
            cell[k].U[1]=u*cell[k].U[0];
		cell[k].U[2]=v*cell[k].U[0];
		cell[k].U[3]=p/(_gamma_-1.0)+0.5*rho*(u*u + v*v);
            cell[k].U[4]=phi;

	}
#endif		
	return 1;
}




////////////////////////////////////////////////////
//// C A L C U L A T I O N  F U N C T I O N S //////
////////////////////////////////////////////////////


int compute_fluxes(t_mesh *mesh, t_sim *sim){	
	double phi3[3],phi5[5],phi7[7]; 
	double lambdaMax;
	int order;
	int n,i,j,k;
	int st[9]; 
	t_wall *wall;

	mesh->lambda_max=0.0;
	lambdaMax=mesh->lambda_max;
#pragma omp parallel for default(none) private(wall,phi3,phi5,phi7,order,i,j,k,st) shared(sim,mesh) reduction(max:lambdaMax)
	for(n=0;n<mesh->nwalls;n++){
		wall=&(mesh->wall[n]);

		//RIGHT RECONSTRUCTION
		if(wall->nx<TOL4){
			order=wall->cellR->st_sizeY;
			for(j=0;j<9;j++){
				st[j]=wall->cellR->stY[j];
			}
		}else{
			order=wall->cellR->st_sizeX;
			for(j=0;j<9;j++){
				st[j]=wall->cellR->stX[j];
			}
		}
		if(order==1){
			for(k=0;k<sim->nvar;k++){
				wall->UR[k]=wall->cellR->U[k];
			}
		}else {

			for(k=0;k<sim->nvar;k++){
				for(i=0;i<order;i++){
					phi3[i]=mesh->cell[st[i]].U[k];
				}
				wall->UR[k]=weno3R(phi3);
			}
		}

		//LEFT RECONSTRUCTION
		if(wall->nx<TOL4){
			order=wall->cellL->st_sizeY;
			for(j=0;j<9;j++){
				st[j]=wall->cellL->stY[j];
			}
		}else{
			order=wall->cellL->st_sizeX;
			for(j=0;j<9;j++){
				st[j]=wall->cellL->stX[j];
			}
		}

		if(order==1){
			for(k=0;k<sim->nvar;k++){
				wall->UL[k]=wall->cellL->U[k];
			}
		}else{
			for(k=0;k<sim->nvar;k++){
				for(i=0;i<order;i++){
					phi3[i]=mesh->cell[st[i]].U[k];
				}
				wall->UL[k]=weno3L(phi3);
			}
		}
  	
	#if LINEAR_TRANSPORT==0
		#if EULER
			#if HLLE
				compute_euler_HLLE(wall,&lambdaMax);
			#endif
			#if ROE
				compute_euler_Roe(wall,&lambdaMax);
			#endif
		#endif
	#endif
		compute_transport(wall);
	}
	mesh->lambda_max=lambdaMax;

	return 1;
}


double weno3R(double *phi){

	double b0, b1;		
	double a0, a1, a_sum;	
	double g0, g1;		
	double w0, w1;		
	double UR;
		
	g0=2.0/3.0;
	g1=1.0/3.0;
    
#if OPT_WEIGHTS == 0  

	b0=(phi[1]-phi[0])*(phi[1]-phi[0]);
	b1=(phi[2]-phi[1])*(phi[2]-phi[1]);

	a0=g0/((b0+epsilon)*(b0+epsilon));
	a1=g1/((b1+epsilon)*(b1+epsilon));
	
	a_sum=a0+a1;

	w0=a0/a_sum;
	w1=a1/a_sum;
    
#else  
    
    w0=g0;
	w1=g1;
    
#endif  
	
	UR=w0*(0.5*phi[1]+0.5*phi[0])+w1*(-0.5*phi[2]+1.5*phi[1]);


	return UR;
}


double weno3L(double *phi){

	double b0, b1;		
	double a0, a1, a_sum;	
	double g0, g1;		
	double w0, w1;		
	double UL;
		
	g0=1.0/3.0;
	g1=2.0/3.0;
    
#if OPT_WEIGHTS == 0   

	b0=(phi[1]-phi[0])*(phi[1]-phi[0]);
	b1=(phi[2]-phi[1])*(phi[2]-phi[1]);

	a0=g0/((b0+epsilon)*(b0+epsilon));
	a1=g1/((b1+epsilon)*(b1+epsilon));
	
	a_sum=a0+a1;

	w0=a0/a_sum;
	w1=a1/a_sum;
    
#else  
    
    w0=g0;
	w1=g1;
    
#endif  
	
	UL=w0*(-0.5*phi[0]+1.5*phi[1])+w1*(0.5*phi[1]+0.5*phi[2]);
      
	return UL;
}




void compute_euler_HLLE(t_wall *wall,double *lambda_max){

	int m;
	double WR[4], WL[4]; /**<Auxiliar array of variables rotated for the 1D problem**/
	double uL, uR, vL, vR, pL, pR, HL, HR, cL, cR;
	double raizrhoR, raizrhoL, sumRaizRho;
	double u_hat, v_hat, H_hat, c_hat;
	double S1, S2, diffS, maxS;
	double FR[4], FL[4];
	double F_star[4];
	
	WR[0]=wall->UR[0];
	WL[0]=wall->UL[0];

	WR[1]=wall->UR[1]*wall->nx+wall->UR[2]*wall->ny;
	WL[1]=wall->UL[1]*wall->nx+wall->UL[2]*wall->ny;

	WR[2]=-wall->UR[1]*wall->ny+wall->UR[2]*wall->nx;
	WL[2]=-wall->UL[1]*wall->ny+wall->UL[2]*wall->nx;
      
	WR[3]=wall->UR[3];
	WL[3]=wall->UL[3];

	uL=WL[1]/WL[0];
	uR=WR[1]/WR[0];
	vL=WL[2]/WL[0];
	vR=WR[2]/WR[0];

	pL=(_gamma_-1.0)*(WL[3]-0.5*WL[0]*(uL*uL+vL*vL));
	pR=(_gamma_-1.0)*(WR[3]-0.5*WR[0]*(uR*uR+vR*vR));
	
	HL=(WL[3]+pL)/WL[0];
	HR=(WR[3]+pR)/WR[0];
	
	cL=sqrt(_gamma_*pL/WL[0]);
	cR=sqrt(_gamma_*pR/WR[0]);
	
	raizrhoL=sqrt(WL[0]);
	raizrhoR=sqrt(WR[0]);
	sumRaizRho=raizrhoR+raizrhoL;

	u_hat=(uR*raizrhoR+uL*raizrhoL)/sumRaizRho;
	v_hat=(vR*raizrhoR+vL*raizrhoL)/sumRaizRho;
	H_hat=(HR*raizrhoR+HL*raizrhoL)/sumRaizRho;

	c_hat=sqrt((_gamma_-1)*(H_hat-0.5*(u_hat*u_hat+v_hat*v_hat)));

	FR[0]=WR[1];
	FL[0]=WL[1];

	FR[1]=WR[1]*uR+pR;
	FL[1]=WL[1]*uL+pL;

	FR[2]=WR[1]*vR;
	FL[2]=WL[1]*vL;

	FR[3]=uR*(WR[3]+pR);
	FL[3]=uL*(WL[3]+pL);

	
	S1=MIN(uL-cL,u_hat-c_hat);
	S2=MAX(uR+cR, u_hat+c_hat);

	maxS=MAX(ABS(S1),ABS(S2));
	diffS=S2-S1;

	for(m=0;m<4;m++){
		if(S1>=0){
			F_star[m]=FL[m];
		}else if(S2<=0){
			F_star[m]=FR[m];
		}else{
                  F_star[m]=(S2*FL[m]-S1*FR[m]+S1*S2*(WR[m]-WL[m]))/(diffS);
            }
	}
	
	wall->fR_star[0]=F_star[0]; 
	wall->fR_star[1]=F_star[1]*wall->nx - F_star[2]*wall->ny;
	wall->fR_star[2]=F_star[1]*wall->ny + F_star[2]*wall->nx;
	wall->fR_star[3]=F_star[3]; 
	
      for(m=0;m<4;m++){
            wall->fL_star[m]=wall->fR_star[m];
      }
      
      
	*lambda_max=MAX(*lambda_max,maxS);
		
}

void compute_transport(t_wall *wall){
	if(wall->fR_star[0]<TOL14){
		wall->fR_star[4]=wall->fR_star[0]*wall->UR[4]/wall->UR[0];
		wall->fL_star[4]=wall->fR_star[4];
	}else{
		wall->fR_star[4]=wall->fL_star[0]*wall->UL[4]/wall->UL[0];
		wall->fL_star[4]=wall->fR_star[4];
	}

}

int update_dt(t_mesh *mesh,t_sim *sim){
	double dl;
	dl=MIN(mesh->dx,mesh->dy);
	sim->dt=sim->CFL*dl/mesh->lambda_max;
	if(sim->dt+sim->t>sim->tf){
		sim->dt=sim->tf-sim->t+TOL8;
	}
	return 1;
}

void update_cell(t_mesh *mesh, t_sim *sim){	
	int i,k;
	t_cell *cell;
	cell=mesh->cell;
	for(i=0;i<mesh->ncells;i++){
		for(k=0;k<sim->nvar;k++){
			cell->U[k]-=sim->dt*((cell->w2->fL_star[k]-cell->w4->fR_star[k])/cell->dx + (cell->w3->fL_star[k]-cell->w1->fR_star[k])/cell->dy);
		}           
		cell++;
	}
}

void update_cellK1(t_mesh *mesh, t_sim *sim){
	
	int i,k;
	t_cell *cell;
#pragma acc parallel loop present(mesh,sim) private(cell,k)
#pragma omp parallel for default(none) private(k,cell) shared(sim,mesh)
	for(i=0;i<mesh->ncells;i++){
		cell=&(mesh->cell[i]);
		for(k=0;k<sim->nvar;k++){
			cell->U_aux[k]=cell->U[k];
			cell->U[k]-=sim->dt*((cell->w2->fL_star[k]-cell->w4->fR_star[k])/cell->dx + (cell->w3->fL_star[k]-cell->w1->fR_star[k])/cell->dy);
		}
	}
}

void update_cellK2(t_mesh *mesh, t_sim *sim){
	
	int i,k;
	t_cell *cell;
#pragma acc parallel loop present(mesh,sim) private(cell,k)
#pragma omp parallel for default(none) private(k,cell) shared(sim,mesh)
	for(i=0;i<mesh->ncells;i++){
            cell=&(mesh->cell[i]);
		for(k=0;k<sim->nvar;k++){
			cell->U[k]=0.75*cell->U_aux[k]+0.25*cell->U[k]-0.25*sim->dt*((cell->w2->fL_star[k]-cell->w4->fR_star[k])/cell->dx + (cell->w3->fL_star[k]-cell->w1->fR_star[k])/cell->dy);

		}
	}
}

void update_cellK3(t_mesh *mesh, t_sim *sim){
	
	int i,k;
	t_cell *cell;
#pragma acc parallel loop present(mesh,sim) private(cell,k)
#pragma omp parallel for default(none) private(k,cell) shared(sim,mesh)    
	for(i=0;i<mesh->ncells;i++){
            cell=&(mesh->cell[i]);
		for(k=0;k<sim->nvar;k++){
			cell->U[k]=(1.0/3.0)*cell->U_aux[k]+(2.0/3.0)*cell->U[k]-(2.0/3.0)*sim->dt*((cell->w2->fL_star[k]-cell->w4->fR_star[k])/cell->dx + (cell->w3->fL_star[k]-cell->w1->fR_star[k])/cell->dy);

		}
	}
}

void mass_calculation(t_mesh *mesh, t_sim *sim){
	int i;
	double massAux;
	double area;
	area=mesh->dx*mesh->dy;
	massAux=0.0;
#pragma omp parallel for default(none) shared(area,mesh) reduction(+:massAux)
	for(i=0;i<mesh->ncells;i++){
		massAux+=mesh->cell[i].U[0]*area;
	}
	mesh->mass=massAux;

}


////////////////////////////////////////////////////
//////////////////// M A I N ///////////////////////
////////////////////////////////////////////////////

int main(int argc, char * argv[]){
	
	int i, j, k, p;
	t_mesh *mesh;
	t_sim *sim;
	char vtkfile[1024];
	double tf,t;
	int nIt; //counter for iterations
	double timeac; //accumulated time before dump
	t_cell *cell;
      

	//Mesh allocation
	mesh=(t_mesh*)malloc(sizeof(t_mesh));
	//Simulation allocation
	sim=(t_sim*)malloc(sizeof(t_sim));


	////////////////////////////////////////////////////
	////////////// P R E - P R O C E S S ///////////////
	////////////////////////////////////////////////////
	//
	// S I M U L A T I O N   P A R A M E T E R S
	sim->tf=2.0;
	sim->tVolc=0.2;
	sim->order=3;
	sim->CFL=0.4;

#if EULER
	sim->nvar=5;
#endif
#ifdef _OPENMP 
	omp_set_num_threads(NTHREADS);
     	printf("The number of threads is set to %d.\n",NTHREADS);
#endif
#pragma omp parallel

	// M E S H   P A R A M E T E R S
	mesh->xcells= 100;
	mesh->ycells= 100;
	mesh->dx= 1.0/mesh->xcells;
	mesh->dy= mesh->dx;//1.0/mesh->ycells;;
	mesh->bc[0]=1;
	mesh->bc[1]=1;
	mesh->bc[2]=1;
	mesh->bc[3]=1;
      
	timeac=0.0;
		
	// S I M U L A T I O N  S E T - U P
	if(sim->order==1){
		sim->rk_steps=1;
	}else{
		sim->rk_steps=3;
	}
	
      
#pragma acc enter data copyin(mesh[:1]) create(sim[:1])

	// S E T  F U N C T I O N S
	create_mesh(mesh,sim);
	update_initial(mesh);
	
 
	////////////////////////////////////////////////////
	////////////// C A L C U L A T I O N ///////////////
	////////////////////////////////////////////////////
	tf=sim->tf;
      sim->t=0.0;
	t=0.0;
	nIt=0;
	while(t<tf){
		
		cell=mesh->cell;
		for(k=1;k<=sim->rk_steps;k++){
			if(k==1){
				compute_fluxes(mesh,sim);
				update_dt(mesh,sim);
				if(sim->order==1){
					update_cell(mesh,sim);
				}else{
					update_cellK1(mesh,sim);
				}

			}else if(k==2){
				compute_fluxes(mesh,sim);
				update_cellK2(mesh,sim);
			}else{ //k=3
				compute_fluxes(mesh,sim);
				update_cellK3(mesh,sim);
			}
		}
            
		////////////////////////////////////////////////////
		///////////// P O S T - P R O C E S S //////////////
		////////////////////////////////////////////////////
		
            timeac=timeac+sim->dt;
		if(timeac>sim->tVolc){			
                  printf("\n");
                  printf("T= %lf, dt= %lf\n",t+sim->dt,sim->dt);
		  mass_calculation(mesh,sim);
		  printf("Total mass: %14.14e\n",mesh->mass);
	        printf("\n");
              nIt++;		
	        timeac=0.0;
		}

		
		t+=sim->dt;
		sim->t=t;
            
	}
	
      printf(" \n");
      printf("Final time is T= %14.14e \n \n",sim->t);
      
      

	return 1;

}

