#include "headers/Transport.h"
#include "headers/Source.h"
#include "headers/Utility.h"
#include "headers/Output.h"
#include "headers/OuterSource.h"
#include "headers/Dispersion.h"
#include "headers/Dose.h"
#include "TApplication.h"
#include "TMath.h"
#include <time.h>

void StandaloneApplication(int argc, char** argv) {
  clock_t tStart = clock();
  // ==>> here the ROOT macro is called

  OuterSource outer_source(false);

  Dispersion dispersion(outer_source);
  dispersion.Go();

  Dose dose(dispersion);
  dose.Go();

  cout << "\nElapsed Time: " << (double)(clock() - tStart) / CLOCKS_PER_SEC << endl;

}
  // This is the standard "main" of C++ starting
  // a ROOT application
int main(int argc, char** argv) {

  TApplication app("ROOT Application", &argc, argv);
  StandaloneApplication(app.Argc(), app.Argv());
  app.Run();

  return 0;
}