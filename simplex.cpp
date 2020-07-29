
/*
SIMPLEX METJOD for Linear Programming
augmented with an algorithm to find the binding constraints
and thus reducing the number of constraints and subsequently
the computing cost
*/

/*
To compile run the command on the unix/linux shell
 
g++ simplex.cpp -o simplex

*/

/*

To run the executable type on the unix/linux shell:

./simplex

*/

#include <stdio.h>
#include <math.h>
#include <unordered_set>
#include <bits/stdc++.h> 

#include <iostream>
#include <vector>
//using namespace std; 


#define CMAX  25  //max. number of variables in economic function
#define VMAX  25  //max. number of constraints

void Data();
void Pivot();
void Simplex();
void Optimize();
void Formula();
void Results();
bool sameSign(double, double);
void popout(std::vector<int> &);




  int NC, NV, NOPTIMAL,P1,P2,XERR;
  double TS[CMAX][VMAX];
  int count_prod = 0;
  int count_div = 0 ;

  void Data() {
    double R1,R2;
    char R;
    int I,J;
    printf("\n ---------------------------------------------------------\n");
    printf("\n Linear programming optimized with binding constraints\n");
    printf("\n ---------------------------------------------------------\n");
    printf(" Maximize (Y/N) ? "); scanf("%c", &R);
    printf("\n Number of variables of Objective function ? "); scanf("%d", &NV);
    printf("\n Number of constrints ? "); scanf("%d", &NC);

    if (R == 'Y' || R=='y')
      R1 = 1.0;
    else
      R1 = -1.0;

    printf("\n Input coefficients of Objective function:\n");

    for (J = 1; J<=NV; J++) {
      printf("       #%d ? ", J); 
      scanf("%lf", &R2);
      TS[1][J+1] = R2 * R1; // c-vector
      count_prod++;
    }

    printf("-------Right hand side ? "); scanf("%lf", &R2);

    TS[1][1] = R2 * R1;
    count_prod++;

    for (I = 1; I<=NC; I++) {
      printf("\n CONSTRAINT #%d:\n", I);
      for (J = 1; J<=NV; J++) {
        printf("       #%d ? ", J); scanf("%lf", &R2);
        TS[I + 1][J + 1] = -R2;
      }
      printf("------Right hand side ? "); scanf("%lf", &TS[I+1][1]);
    }

    
    printf("\n\n RESULTS:\n\n");
    for(J=1; J<=NV; J++)  TS[0][J+1] = J;
    for(I=NV+1; I<=NV+NC; I++)  TS[I-NV+1][0] = I;
    
  }

  void Simplex() {

e10: Pivot();
     Formula();
     Optimize();
     if (NOPTIMAL == 1) goto e10;
  }

  void Pivot() {

    double RAP,V,XMAX;
    int I,J;

    XMAX = 0.0;
    for(J=2; J<=NV+1; J++) {
	if (TS[1][J] > 0.0 && TS[1][J] > XMAX) {
        XMAX = TS[1][J];
        P2 = J;
      }
    }
    RAP = 999999.0;
    for (I=2; I<=NC+1; I++) {
      if (TS[I][P2] >= 0.0) goto e10;
      V = fabs(TS[I][1] / TS[I][P2]);
      count_div++;
      if (V < RAP) {
        RAP = V;
        P1 = I;
      }
e10:;}
    V = TS[0][P2]; TS[0][P2] = TS[P1][0]; TS[P1][0] = V;
  }

  void Formula() {;
    //Labels: e60,e70,e100,e110;
    int I,J;

    for (I=1; I<=NC+1; I++) {
      if (I == P1) goto e70;
      for (J=1; J<=NV+1; J++) {
        if (J == P2) goto e60;
        TS[I][J] -= TS[P1][J] * TS[I][P2] / TS[P1][P2];
        count_div++;
e60:;}
e70:;}
    TS[P1][P2] = 1.0 / TS[P1][P2];
    count_div++;
    for (J=1; J<=NV+1; J++) {
      if (J == P2) goto e100;
      TS[P1][J] *= fabs(TS[P1][P2]);
      count_prod++;

e100:;}
    for (I=1; I<=NC+1; I++) {
      if (I == P1) goto e110;
      TS[I][P2] *= TS[P1][P2];
      count_prod++;

e110:;}
  }   

  void Optimize() {
    int I,J;
    for (I=2; I<=NC+1; I++)
      if (TS[I][1] < 0.0)  XERR = 1;
    NOPTIMAL = 0;
    if (XERR == 1)  return;
    for (J=2; J<=NV+1; J++)
      if (TS[1][J] > 0.0)  NOPTIMAL = 1;
  }

  void Results() {
    //Labels: e30,e70,e100;
    int I,J;
    
    if (XERR == 0) goto e30;
    printf(" NO SOLUTION.\n"); goto e100;
e30:for (I=1; I<=NV; I++)
    for (J=2; J<=NC+1; J++) {
      if (TS[J][0] != 1.0*I) goto e70;
      printf("       VARIABLE #%d: %f\n", I, TS[J][1]);
e70:  ;}
    printf("\n       ECONOMIC FUNCTION: %f\n", TS[1][1]);
e100:printf("\n");
  }


bool sameSign(double num1, double num2)
{
    return num1 > 0 && num2 > 0 || num1 < 0 && num2 < 0;
}

 
void popout(std::vector<int> &v)
{
 
    std::unordered_set<int> s(v.begin(), v.end());
    v.assign(s.begin(), s.end());

}



void Bindings()
{

    double A[NC][NV];
    double c[NV];
    double b[NC];    
    double lambda[NC];
    double zlambda[NC];
    double zbx[NC];
    double bx[NC];
    double sum,r_kj,r_ij,minm;
    int sgn_1, sgn_2;
    int mindx;
    int binding[NC] = {0};
    std::vector<int> bindconstraints;
        
    std::cout << "RHS vector: \n" ;

    for (int I = 0;I < NC ; I++){
      b[I] = TS[I+2][1]; 
      std::cout << b[I] << "\n";    
    }

    std::cout << "Objective function coefficients\n";
    for (int J= 0 ;J < NV; J++){
      c[J] = TS[1][J+2];
      std::cout << c[J] << "\n";
    }


    for (int I = 0;I < NC; I++)
      for (int J= 0;J <  NV; J++){
        A[I][J] =  (-1)*TS[I+2][J+2];
        count_prod++; }



    std::cout << "Constraints matrix\n";
    for (int I = 0;I < NC; I++)
    {
        printf("\n ");
      for (int J= 0; J< NV ; J++) {
      printf("%lf ", A[I][J]); 
    }   
  }

  std::cout << "\n\n\n";
  std::cout << "The LAMBDA vector: \n";

  // Compute LAMBDA vector
  for (int I= 0; I < NC; I++){
    sum = 0.0;{
    for (int J=0; J < NV; J++)
      sum = sum + A[I][J];}
    lambda[I] = b[I]/sum;
    count_div++;
    std::cout << lambda[I] << "\n";
  } 

  // find lambda_min and k
  minm = lambda[0];

  for (int I= 0; I < NC; I++){
    if (lambda[I] < minm){
      minm = lambda[I];
      mindx = I;
    } }

  //std::cout << "Min. LAMBDA and its index\n";
  std::cout << "First binding constraint: " << " " << mindx+1 << "\n";
  
  int booltemp;  

  booltemp = mindx;

  bindconstraints.push_back(mindx);

  //std::cout << *bnd_constraint << std::endl;
    
  // compute z(lambda)
  std::cout << "Check z(lambda)vector \n";  

  for(int I = 0; I < NC; I++){
    sum = 0;
    for (int J = 0; J < NV; J++){
      sum +=  lambda[I]*c[J];
      count_prod++;

      }
    zlambda[I] = sum;
    std::cout << zlambda[I] << " ";
    }

  std::cout << std::endl;

  int count_binds = 1;

  // main loop

  for(int J = 0; J < NV; J++){ 
    for(int I = 0; I < NC; I++){
      bx[I] = b[I]/A[I][J];
      count_div++;
      zbx[I] = c[J]*bx[I];
      count_prod++;     
    } 
  
  r_kj = (zbx[mindx] - zlambda[mindx]) / (bx[mindx] - lambda[mindx]);
  count_div++;

  
   for (int I = 0; I < NC; I++){

    if (I != booltemp & I != mindx)
    {
      r_ij = (zbx[I] -zlambda[I])/(bx[I] - lambda[I]);
      count_div++;

    if(r_kj != 0 & r_ij !=0)
      if (sameSign(r_kj,r_ij) == 0) {
        booltemp = I;
        count_binds++;
        bindconstraints.push_back(I);

      }
    }
   }

  } //main loop

  
 
    popout(bindconstraints);

    sort(bindconstraints.begin(), bindconstraints.end()); 

    std::cout << "The binding constraints are:";
    for(int i = 0; i < bindconstraints.size(); i++)
      std::cout << bindconstraints[i]+1<< " ";

    int newsize = bindconstraints.size();

    
    printf("\n%s\n","Constraints RHS after deleting all non-binding constraints" );
     for (int I = 0;I < newsize ; I++){
      TS[I+2][1] = b[bindconstraints[I]]; ; 
      std::cout << TS[I+2][1] << "\n";    
    }


     std::cout << std::endl;

  }


int main()  {

  Data();
  Bindings();


  Simplex();
  Results();

 

  std::cout << "The total number of DIVISION operations: " << count_div << "\n";
  std::cout << "The total number of MULTIPLICATION operations: " << count_prod << "\n";
}

