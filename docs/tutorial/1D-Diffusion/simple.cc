#include <stdio.h>
#include <iostream>

int max_iterations      = 1;  
int simulation_finished = 0;

int      iter = 0;
double  *nu, *xn, *xc, *unew, *uold, *dt, *gradient;   
int     numCells, numNodes, numFaces;

void CellCentroid();
void LeftBoundaryCondition();
void RightBoundaryCondition();
void TimeStep();
void Gradient();
void Advance();
void CollapseCondition();
void InitialConditions();
void Update();

int main()
{
  
  numNodes = 11;
  // Generate the node Coordinates ...

  xn = new double[numNodes];
  double dx = 1.0/(numNodes-1);
  for(int  i = 0; i < numNodes; i++) {
    xn[i] = (double) i*dx;
  }
   
  // Calculate the cell centroids ...
  CellCentroid();

  // At each cell-centroid, provide initial values ...
  InitialConditions();
 
  while(!simulation_finished ) {
    iter++;
    LeftBoundaryCondition();
    RightBoundaryCondition();
    TimeStep();
    Advance();
    CollapseCondition();
    Update();
  }

  cout << " FINAL SOLUTION " << endl;

  for( int icell = 0; icell < numCells; icell++)
       cout << unew[icell] << endl;
 
}

//****************************************************************

void CellCentroid()
{
  numCells = numNodes-1;

  xc = new double[numCells];
  
  for( int icell = 0; icell < numCells; icell++) {
    xc[icell] = 0.5*(xn[icell]+xn[icell+1]);    
  }
   
}
//****************************************************************
void Advance()
{
  int     icell, lface, rface;
  double  dudx;

  Gradient();
  
  for( icell = 0; icell < numCells; icell++) {
    rface       =  icell+1;
    lface       =  icell-1;
    dudx        = (gradient[rface]-gradient[lface])/
      (xn[rface]-xn[lface]);
    unew[icell] =  uold[icell] + nu[icell]*dt[icell]*dudx;    
  }

}

//****************************************************************

void LeftBoundaryCondition()
{
  gradient[0] = -1;
  
}

//**************************************************************

void RightBoundaryCondition()
{
  gradient[numNodes-1] =   uold[numCells-1]/
                           (xc[numCells-1]-xn[numNodes-1]); 
}
//***************************************************************

void CollapseCondition()
{
  double diffval, maxError = 0.0;

  for( int icell = 0; icell < numCells; icell++){
    diffval = unew[icell]-uold[icell];
    if( fabs(diffval) > maxError)
      maxError = diffval;
  }

  if( maxError < 1.0E-06 || iter >= max_iterations) 
      simulation_finished = 1;
  
  cout << "ITER = " << iter << " MAX ERROR " << maxError << endl;
  
}
//***************************************************************

void Gradient()
{
 
  int  iface, rCell, lCell;
  
  for( iface = 1; iface < numFaces-1; iface++) {
    rCell           = iface;
    lCell           = iface-1;
    gradient[iface] = (uold[rCell] - uold[lCell])/
                      (xc[rCell]-xc[lCell]);
  }
                                        
}

//****************************************************************

double f(double x) {
  return 0.0;
}

//****************************************************************

void InitialConditions()
{
  uold     = new double[numCells];
  unew     = new double[numCells];
  gradient = new double[numNodes];
  dt       = new double[numCells];
  nu       = new double[numCells];

  numFaces = numNodes;
  
  for(int icell=0; icell < numCells; icell++){
    uold[icell] =   f(xc[icell]);
    nu[icell]   =   1.0;
  }
  
}

//****************************************************************

void TimeStep()
{

  int    icell;
  double dx, dtmin;
  
  for(icell = 0; icell < numCells; icell++){
    dx        = xc[icell+1]-xc[icell-1];
    dt[icell] = dx*dx*nu[icell]/4.0;
  }

  dtmin = 1.0E+10;
  for( icell = 0; icell < numCells; icell++)
    if( dt[icell] < dtmin ) dtmin = dt[icell];

  for( icell = 0; icell < numCells; icell++)
    dt[icell] = dtmin;

}
//****************************************************************

void Update()
{

  for( int icell = 0; icell < numCells; icell++)
       uold[icell] = unew[icell];
  
}
//****************************************************************
