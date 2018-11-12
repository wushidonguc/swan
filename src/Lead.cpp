#include "Lead.h"

Lead::~Lead(){}

//cx_mat Lead::GetSurfaceGreensFunction(const complex<double> energy, const double eta, const double phase) const
//{
//  cx_mat es = mpMaterial->GetPrincipalLayerHamiltonian(phase);
//  cx_mat ee = mpMaterial->GetPrincipalLayerHamiltonian(phase);
//  cx_mat a = mpMaterial->GetPrincipalLayerInteractionX(phase);
//  cx_mat b = mpMaterial->GetPrincipalLayerInteractionX(phase).t();
//
//  int nRows = es.n_rows;
//  int nCols = es.n_cols;
// 
//  cx_mat es_temp(nRows, nCols);
//  cx_mat ee_temp(nRows, nCols);
//  cx_mat a_temp(nRows, nCols);
//  cx_mat b_temp(nRows, nCols);
//  
//  cx_mat E(nRows, nCols);
//  cx_mat g00(nRows, nCols);
//
//  E = (energy + eta * I) * eye(nRows, nCols);
//
//  while ((abs(a)).max() > mCriteria or (abs(b)).max() > mCriteria)
//  {
//    a_temp = a, b_temp = b, es_temp = es, ee_temp = ee;
//    a = a_temp * (E - ee_temp).i() * a_temp;
//    b = b_temp * (E - ee_temp).i() * b_temp;
//    es = es_temp + a_temp * (E - ee_temp).i() * b_temp;
//    ee = ee_temp + a_temp * (E - ee_temp).i() * b_temp + b_temp * (E - ee_temp).i() * a_temp;
//  }
//
//  g00 = (E - es).i();
//
//  return g00;
//}


//// Testing
//#include "Graphene.h"
//int main()
//{
//  // Input parameters
//  Parameters parameters;
//  parameters.ParseInputFile();
//  parameters.Print();  
//
//  Graphene* pGraphene = new Graphene();
//  pGraphene->Init();
//  
//  Lead graphene_lead(parameters, pGraphene);
//  cout << graphene_lead.GetSurfaceGreensFunction(1.0, 0.01, 0.) << endl;
//
//  return 0;
//}
