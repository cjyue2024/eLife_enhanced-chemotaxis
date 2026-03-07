#include "MersenneTwister.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <thread>
#include <omp.h>

using namespace std;


const int N = 100;              // 细胞数目 

// Swim related parameters
const double size = 160.0;    // simulation size chanelwidth*400 um^2
const double v0 = 20.0;       // swimming speed um/s
const double pi = 3.14159;

const double D0 = 0.062;     // rotation diffusion coefficient. (rad^2/s)0.062

// Chemotaxis related parameters.   SPEC model, Tu, Plos conputational biology,2010.
const double kR = 0.005; 
const double a0 = 0.3333;
const double kTumble = 5.0; // switch rate of Run to Tumble 

const double kL = 0.05;  // gradient of attractant conc. uM/um     0.1uM/um
const double L0 = 25.0;  // attractant conc. at x = 0 um.          15.0 uM

void simulation(double cwdh, double rc, int Ns) {
			
	MTRand myrandx,myrandy,myrandtheta,myrandRunOrTumble,myrand6,myrand8;
    	
	double lig, a, kRun; // Linear gradient field with L(x) = kL*x/size + L0.  // receptor activity  // switch rate of Run to Tumble 1.3636
    
	double x[N],y[N],theta[N],vx[N],vy[N],m[N],Run[N],B0[N];
	double tdetention = 0.0; // 侧壁的运动时间 ****************/
	double dt=0.05;       // time step
	
	int trackindex[N];
	int tracki = N;
	
	char s[200];		
	//cout<<"Simulating chanelwidth = "<<cwdh<<"um / L0 = "<<L0<<"/ kL = "<<kL<<"/ N = "<<N<<"/ "<<Ns<<endl; 
	sprintf(s, "temp_chanelwidth%04.0f_D%1.3f_r%1.2f_L0%.0f_kL%1.3f_N%04d_%03d.txt",cwdh,D0,rc,L0,kL, N,Ns);
	ofstream outfile1(s);
//*************************细菌状态初始化**************************************************************************
	for(int i=0;i<N;i++) {
        trackindex[i] = i+1;
		// Direction bias for single cell. (rad/s)
		if(rc==0.0){
			B0[i] = 0.0;
		}else{
			B0[i] = -v0/rc;  
		}
		
		// swim related.
	  	x[i] = size*myrandx.rand(); //  Random the x position [0,size].
		y[i] = cwdh*myrandy.rand(); //  Random the x position [0,cwdh].
		theta[i] = pi*(2.0*myrandtheta.randExc()-1.0); // Random the velocity direction [-pi,pi).
		
		// chemotaxis related.
		lig = kL*x[i]+L0; // ============================================================ Initial ligand conc. uM
		m[i] = 1.0 - ((log(1.0/a0 -1.0))/4.6 - (log((1+lig/1.7)/(1+lig/12.0))))/1.7; //  Methylation level at steady state a0.
		kRun = (pow(7.86*a0,10.3)/(pow(3.1,10.3)+pow(7.86*a0,10.3)))/0.11; // =========== Swith rate of CCW (Run) to CW (Tumble).
		
        kRun = kRun/2.8;  /*****************/ //run length = 2.0 s   14 PRL
        
		Run[i] = myrandx.randInt(1);       // =========================================== Swim state [Run(1) or tumble(0)] at steady state a0.
		if(Run[i]==1){
			if(kRun*dt>=myrandRunOrTumble.randDblExc()){Run[i] = 0;}
		}
		else{
			if(kTumble*dt>=myrandRunOrTumble.randDblExc()){Run[i] = 1;}
		}
	}

//*****************************循环更新位置与速度信息***************************************************************** 
	double t=0.0;         // current time
	while(t<150.0) {
		for(int i=0;i<N;i++) {
	    	// 确定此时第i个细菌的运动状态是run还是tumble （Run = 1 or 0）;
	    	
	          lig = kL*x[i]+L0;
			  a = 1.0/(1.0+exp(4.6*(1.7*(1.0-m[i])+log((1+lig/1.7)/(1+lig/12.0)))));  //yp = 7.86*a;//cwbias=pow(7.86*a,10.3)/(pow(7.86*a,10.3)+pow(3.1,10.3));
		      kRun = (pow(7.86*a,10.3)/(pow(7.86*a,10.3)+pow(3.1,10.3)))/0.11; // kRun = cwbias/0.11;
		      
		      kRun = kRun/2.8;   /*****************/ //run length = 2.0 s   14 PRL

			 double XXX,YYY,aa,bb,cc,aa1,bb1,cc1;
		     MTRand randXXX, randYYY;
	
			  if(Run[i]==1){
			  	if(kRun*dt>=myrandRunOrTumble.randDblExc()){ // Run to Tumble change angle follow tumble angle distribution.
					   Run[i] = 0;
				  }
				  else {
					   Run[i] = 1;
					   if(y[i]<cwdh and y[i]>0.0){
					   	 theta[i] += B0[i]*dt+sqrt(2.0*D0*dt)*myrand6.randNorm(0.0,1.0);
				  	   }
		      	  }
			  }
			  else{
			  	if(kTumble*dt>=myrandRunOrTumble.randDblExc()){  // tumble to run  change direction.
				  	Run[i] = 1;
					if(y[i]<cwdh and y[i]>0.0){
						if(y[i]>3.0 and y[i]<cwdh-3.0){bb = 1.35;}
						else{bb = 0.25;} 
						//解析式计算tumble角度分布，从均匀随机分布转换而来 P = 1.0/bb/(1.0-exp(-pi/bb))*exp(-x/bb)
						theta[i] += (2.0*randXXX.randInt(1)-1.0)*bb*log((exp(-pi/bb)-1.0)*randYYY.rand()+1.0);
					}
				 }
		      }
	         m[i] += kR*((1.0-a)-(1.0-a0)*a/a0)*dt;
	        		
			tdetention = 2.0;
			bb = 0.25;
			if(dt/tdetention>=myrand8.randDblExc()){
				if(y[i]==0){
					//Run[i] = 1;	
					if(cos(theta[i])>0.0){
						theta[i] += -bb*log((exp(-pi/bb)-1.0)*randYYY.rand()+1.0);				
					}else{
						theta[i] += bb*log((exp(-pi/bb)-1.0)*randYYY.rand()+1.0);
					}
				}
			
				if(y[i]==cwdh){
					//Run[i] = 1;				
					if(cos(theta[i])>0.0){
						theta[i] += bb*log((exp(-pi/bb)-1.0)*randYYY.rand()+1.0);
					}else{
						theta[i] += -bb*log((exp(-pi/bb)-1.0)*randYYY.rand()+1.0);
					}
				}
			}
			
		// 通过以上的所有影响确定当前细菌的位置与速度方向状态 	  	
			if(theta[i]<-pi){theta[i]+=2.0*pi;}
			if(theta[i]>=pi){theta[i]-=2.0*pi;}
			x[i] += Run[i]*v0*cos(theta[i])*dt;
			y[i] += Run[i]*v0*sin(theta[i])*dt;
			
		//	x-direction 周期性边界条件 	
		  	if(x[i]>size) {
		  		tracki++;
		  		trackindex[i] = tracki;
				x[i]-=size;
				y[i] = cwdh*myrandy.rand();
				lig = kL*x[i]+L0;
				m[i] = 1.0 - ((log(1.0/a0 -1.0))/4.6 - (log((1+lig/1.7)/(1+lig/12.0))))/1.7;
			}
			 
		    if(x[i]<0.0) {
		    	tracki++;
		  		trackindex[i] = tracki;
			  	x[i]+=size;
			  	y[i] = cwdh*myrandy.rand();
			  	lig = kL*x[i]+L0;
			   	m[i] = 1.0 - ((log(1.0/a0 -1.0))/4.6 - (log((1+lig/1.7)/(1+lig/12.0))))/1.7;
			}
		
			// y-direction solid wall aliengment. 
			if(y[i]>cwdh) {
				y[i]=cwdh;
				if(cos(theta[i])>0.0){theta[i]=0.0;}
				else{theta[i]=-pi;}	
			}
				
			if(y[i]<0.0) {
				y[i]=0.0;
				if(cos(theta[i])>0.0){theta[i]=0.0;}
				else{theta[i]=-pi;}
			}
			outfile1<<t<<' '<<x[i]<<' '<<y[i]<<' '<<trackindex[i]<<endl;	
		}
		 
	    t+=dt;
    }
	outfile1.close();	
} 


int main(){
	double channelwidth[17] = {4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,20.0,25.0,35.0,50.0};     // chanel width in um.
	double r_circle[8] = {4.0,6.0,8.0,10.0,12.0,14.0,16.0,0.0};

//	double channelwidth[8] = {4.0,6.0,8.0,10.0,12.0,15.0,25.0,50.0};     // chanel width in um.
//	double r_circle[1] = {10.0};
	int kindex = 1;
	for(int i=0;i<50;i=i+1){ // number of simulation for each parameters.
		for	(int j=0;j<17;j=j+1){ // number of r_circle.
			cout << "8 threads finished: "<<"index = "<<kindex<<"/50"<< std::endl;

//			std::thread t01(simulation, channelwidth[0],r_circle[j],i);
//			std::thread t02(simulation, channelwidth[1],r_circle[j],i);
//			std::thread t03(simulation, channelwidth[2],r_circle[j],i);
//			std::thread t04(simulation, channelwidth[3],r_circle[j],i);
//			std::thread t05(simulation, channelwidth[4],r_circle[j],i);
//			std::thread t06(simulation, channelwidth[5],r_circle[j],i);
//			std::thread t07(simulation, channelwidth[6],r_circle[j],i);
//			std::thread t08(simulation, channelwidth[7],r_circle[j],i);


			std::thread t01(simulation, channelwidth[j],r_circle[0],i);
			std::thread t02(simulation, channelwidth[j],r_circle[1],i);
			std::thread t03(simulation, channelwidth[j],r_circle[2],i);
			std::thread t04(simulation, channelwidth[j],r_circle[3],i);
			std::thread t05(simulation, channelwidth[j],r_circle[4],i);
			std::thread t06(simulation, channelwidth[j],r_circle[5],i);
			std::thread t07(simulation, channelwidth[j],r_circle[6],i);
			std::thread t08(simulation, channelwidth[j],r_circle[7],i);
	
			// Wait for the threads to finish
		    t01.join();
		    t02.join();
		    t03.join();
		    t04.join();
		    t05.join();
		    t06.join();
		    t07.join();
		    t08.join();
	
		    kindex++;
	   }
    }
    
    return 0;
}

/*
// 在表面驻留的时间 
if(dt/6.0>=myrand8.randDblExc()){
	tracki++;
	trackindex[i] = tracki;
	x[i] = size*myrandx.rand(); 
	y[i] = cwdh*myrandy.rand(); 
	lig = kL*x[i]+L0;
   	m[i] = 1.0 - ((log(1.0/a0 -1.0))/4.6 - (log((1+lig/1.7)/(1+lig/12.0))))/1.7;	
}
*/

//蒙特卡罗随机tumble角度分布 P = aa*exp(-x/bb)+cc
/*		
aa = 1.05; 
bb = 0.71;
cc = 0.08;

int kkk = 0;
while(kkk==0) {	
	YYY = (2.0*randYYY.rand()-1.0)*pi;
	XXX = randXXX.rand()*(aa+cc);
	if(XXX<=aa*exp(-fabs(YYY)/bb)+cc){
	   		kkk=1;
			theta[i] += YYY;
	}
}	
*/
